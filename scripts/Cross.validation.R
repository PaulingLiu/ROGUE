library(SingleCellExperiment);library(tidyverse);library(reticulate);library(tidyverse)
library(ggplot2);library(scmap);library(fmsb);library(ggsci);library(scibet)
library(Seurat);library(M3Drop);library(ROCR);library(cluster);library(parallel);library(ROGUE)

expr <- readr::read_csv("/data1/pauling/04_SEmodel/07_NC_revision/03.data/01.clustering.ari/GSM3618014_gene_count.csv.gz")
gene.name <- expr$gene_id
expr <- as.matrix(expr[,-1])
rownames(expr) <- gene.name
expr <- t(expr)

cda <- readr::read_csv("/data1/pauling/04_SEmodel/07_NC_revision/03.data/01.clustering.ari/01.10x.5cl/sc_10x_5cl.metadata.csv.gz")
cda <- cda %>% 
  dplyr::select(sample, cell_line) %>% 
  dplyr::rename(label = cell_line)

expr <- expr[cda$sample,]

m3d_fun <- function(expr){ 
  expr <- as.matrix(expr)
  expr <- t(expr)
  norm <- M3DropConvertData(expr, is.counts=TRUE)
  DEgenes <- M3Drop::M3DropFeatureSelection(norm, suppress.plot = T, mt_threshold = 2)
  DEgenes <- DEgenes %>% 
    dplyr::arrange(p.value)
  return(DEgenes)
}
HVG_fun <- function(expr){  
  expr <- t(expr)
  hvg.res <- BrenneckeGetVariableGenes(expr, suppress.plot = T, fdr = 2)
  return(hvg.res)
}
Gini_fun <- function(expr){ 
  
  calcul.gini = function(x, unbiased = TRUE, na.rm = FALSE){
    if (!is.numeric(x)){
      warning("'x' is not numeric; returning NA")
      return(NA)
    }
    if (!na.rm && any(na.ind = is.na(x)))
      stop("'x' contain NAs")
    if (na.rm)
      x = x[!na.ind]
    n = length(x)
    mu = mean(x)
    N = if (unbiased) n * (n - 1) else n * n
    ox = x[order(x)]
    dsum = drop(crossprod(2 * 1:n - n - 1,  ox))
    dsum / (mu * N)
  }
  
  expr <- t(expr)
  ExprM.RawCounts <- expr
  
  
  minCellNum = 0
  minGeneNum = 0
  expressed_cutoff = 1
  gini.bi = 0
  log2.expr.cutoffl = 0
  log2.expr.cutoffh = 30
  Gini.pvalue_cutoff = 0.0001
  Norm.Gini.cutoff = 1
  span = 0.9
  outlier_remove = 0.75
  GeneList = 1
  Gamma = 0.9
  diff.cutoff = 1
  lr.p_value_cutoff = 1e-5
  CountsForNormalized = 100000
  
  ExpressedinCell_per_gene=apply(ExprM.RawCounts,1,function(x) length(x[x > expressed_cutoff ]))
  nonMir = grep("MIR|Mir", rownames(ExprM.RawCounts), invert = T)  # because Mir gene is usually not accurate 
  Genelist = intersect(rownames(ExprM.RawCounts)[nonMir],rownames(ExprM.RawCounts)[ExpressedinCell_per_gene >= minCellNum])
  ExpressedGene_per_cell=apply(ExprM.RawCounts[Genelist,],2,function(x) length(x[x>0]))
  ExprM.RawCounts.filter = ExprM.RawCounts[Genelist,ExpressedGene_per_cell >= 0]
  
  if(gini.bi==0){
    gini = apply(as.data.frame(ExprM.RawCounts.filter), 1, function(x){calcul.gini(as.numeric(x)) } )    #theoretically, gini have very low chance to have a 1 value
    GiniIndex = as.data.frame(cbind(1:dim(ExprM.RawCounts.filter)[1], gini))
  } else {
    GiniIndex1 <- as.data.frame(apply(ExprM.RawCounts.filter, 1, function(x){calcul.gini(as.numeric(x)) } ) )
    GiniIndex2 <- as.data.frame(apply(ExprM.RawCounts.filter+0.00001, 1, function(x){calcul.gini(as.numeric(1/x)) } ) ) #bi directional
    GiniIndex  <- cbind(GiniIndex1, GiniIndex2)
    colnames(GiniIndex)=c("gini1","gini2")
    GiniIndex$gini2_sign = 0 - GiniIndex$gini2;
    GiniIndex$gini = apply(GiniIndex, 1, max)
    GiniIndex <- na.omit(GiniIndex)
    GiniIndex$gini_sign = GiniIndex$gini
    for(genei in 1:dim(GiniIndex)[1])
    {
      GiniIndex[genei, 5] = ifelse(  GiniIndex[genei, 1] > GiniIndex[genei,2], "up-regulation", "down-regulation") 
    }
  }
  
  Maxs          = apply(ExprM.RawCounts.filter,1,max)
  Means         = apply(ExprM.RawCounts.filter,1,mean)
  log2.Maxs     = log2(Maxs+0.1)
  ExprM.Stat1   = as.data.frame(cbind(Maxs,GiniIndex$gini,log2.Maxs))
  colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs")
  ExprM.Stat1 = ExprM.Stat1[ExprM.Stat1$log2.Maxs>log2.expr.cutoffl & ExprM.Stat1$log2.Maxs<=log2.expr.cutoffh ,]  # is this necessary?
  log2.Maxs = ExprM.Stat1$log2.Maxs
  Gini      = ExprM.Stat1$Gini
  Maxs      = ExprM.Stat1$Maxs
  
  # .3 fitting in max-gini space 
  Gini.loess.fit        = loess(Gini~log2.Maxs, span=span, degree=1)
  Normlized.Gini.Score  = Gini.loess.fit$residuals   #residuals = Gini - Gini.fitted
  Gini.fitted           = Gini.loess.fit$fitted    
  ExprM.Stat1           = as.data.frame(cbind(ExprM.Stat1[,c("Maxs","Gini", "log2.Maxs")], Normlized.Gini.Score, Gini.fitted))
  colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs", "Norm.Gini", "Gini.fitted")
  
  
  ### remove 25% of first round outlier genes, do second round loess
  Gini.loess.fit.residual = residuals(Gini.loess.fit)                               
  thresh.outlier = quantile(Gini.loess.fit.residual[Gini.loess.fit.residual>0], outlier_remove) 
  id.genes.loess.fit = which(Gini.loess.fit.residual < thresh.outlier)               
  id.outliers.loess.fit = which(Gini.loess.fit.residual >= thresh.outlier)          
  log2.Maxs.genes = log2.Maxs[id.genes.loess.fit]                                   
  log2.Maxs.outliers = log2.Maxs[id.outliers.loess.fit]                            
  Gini.loess.fit.2 = loess(Gini[id.genes.loess.fit]~log2.Maxs[id.genes.loess.fit], span=span, degree = 1)
  Gini.loess.fit.2.predict = predict(Gini.loess.fit.2)  
  
  Gini.loess.fit.2.x.y = cbind(log2.Maxs.genes,Gini.loess.fit.2.predict)
  Gini.loess.fit.2.x.y.uniq = as.data.frame(unique(Gini.loess.fit.2.x.y))
  Gini.loess.fit.2.x.y.uniq = Gini.loess.fit.2.x.y.uniq[order(Gini.loess.fit.2.x.y.uniq[,1]),]
  log2.Maxs.genes.sorted = log2.Maxs.genes[order(log2.Maxs.genes)]                   
  Gini.loess.fit.2.predict.sorted = Gini.loess.fit.2.predict[order(log2.Maxs.genes)] 
  #using Gini.loess.fit.2 as model, predict gini value for those outlier which are not used for build model.
  #for each max in outliers set, find the id of max value which is most close in fitted data set
  loc.outliers = apply(matrix(log2.Maxs.outliers),1,function(x){
    if(x<max(log2.Maxs.genes.sorted)){
      return(which(log2.Maxs.genes.sorted>=x)[1])
    }else{
      return(which.max(log2.Maxs.genes.sorted))
    }})                
  #check the results
  outlier_max_in_fit <- cbind(log2.Maxs.outliers, loc.outliers, log2.Maxs.genes.sorted[loc.outliers])
  
  #based on Gini.loess.fit.2, predict outliers which was not used for fitting
  Gini.outliers.predict = apply(cbind(seq(length(log2.Maxs.outliers)),log2.Maxs.outliers),1,function(x){
    id = x[1]
    value = x[2]
    if(value == log2.Maxs.genes.sorted[loc.outliers[id]]){
      return(as.numeric(Gini.loess.fit.2.x.y.uniq[which(Gini.loess.fit.2.x.y.uniq$log2.Maxs.genes>=value)[1],2]))
    }else{
      if(loc.outliers[id]>1){
        return(Gini.loess.fit.2.predict.sorted[loc.outliers[id]-1]+(Gini.loess.fit.2.predict.sorted[loc.outliers[id]]-Gini.loess.fit.2.predict.sorted[loc.outliers[id]-1])*(value-log2.Maxs.genes.sorted[loc.outliers[id]-1])/(log2.Maxs.genes.sorted[loc.outliers[id]]-log2.Maxs.genes.sorted[loc.outliers[id]-1]))
      }else{
        return(Gini.loess.fit.2.predict.sorted[2]-(Gini.loess.fit.2.predict.sorted[2]-Gini.loess.fit.2.predict.sorted[1])*(log2.Maxs.genes.sorted[2]-value)/(log2.Maxs.genes.sorted[2]-log2.Maxs.genes.sorted[1]))
      }
    }
  })
  
  #plot outliers predict results
  outliers.precit.x.y.uniq = as.data.frame(unique(cbind(log2.Maxs.outliers, Gini.outliers.predict)))
  #plot(outliers.precit.x.y.uniq)
  #plot whole fit2 
  colnames(outliers.precit.x.y.uniq) = colnames(Gini.loess.fit.2.x.y.uniq)
  Gini.loess.fit.2.full.x.y.uniq = rbind(Gini.loess.fit.2.x.y.uniq, outliers.precit.x.y.uniq)
  #plot(Gini.loess.fit.2.full.x.y.uniq)
  
  #calcualte Normlized.Gini.Score2
  Normlized.Gini.Score2                        = rep(0,length(Gini.loess.fit.residual))               
  Normlized.Gini.Score2[id.genes.loess.fit]    = residuals(Gini.loess.fit.2)                         
  Normlized.Gini.Score2[id.outliers.loess.fit] = Gini[id.outliers.loess.fit] - Gini.outliers.predict 
  
  Gini.fitted2           = Gini - Normlized.Gini.Score2         
  ExprM.Stat1            = as.data.frame(cbind(ExprM.Stat1[,c("Maxs","Gini", "log2.Maxs", "Gini.fitted", "Norm.Gini" )], Gini.fitted2, Normlized.Gini.Score2))
  colnames(ExprM.Stat1)  = c("Maxs","Gini","log2.Maxs", "Gini.fitted","Norm.Gini",  "Gini.fitted2", "Norm.Gini2")
  Gini.pvalue            = pnorm(-abs(scale(ExprM.Stat1$Norm.Gini2, center=TRUE,scale=TRUE)))
  ExprM.Stat2            = cbind(ExprM.Stat1, Gini.pvalue)  #first time use ExprM.Stat2
  colnames(ExprM.Stat2)  = c("Maxs","Gini","log2.Maxs", "Gini.fitted","Norm.Gini",  "Gini.fitted2", "Norm.Gini2", "p.value")
  
  ExprM.Stat2 %>%
    tibble::rownames_to_column(var = 'Gene') %>%
    dplyr::arrange(p.value)
}
sct_fun <- function(expr){
  expr <- t(expr)
  colnames(expr) <- paste0("Cell",1:ncol(expr))
  sce <- CreateSeuratObject(counts = expr)
  sce <- Seurat::SCTransform(sce, verbose = FALSE, do.center = F, variable.features.n = nrow(sce))
  sce@assays$SCT@meta.features %>% 
    tibble::rownames_to_column(var = "Gene") %>%
    dplyr::arrange(desc(residual_variance)) -> sct.res
  
  return(sct.res)
}
FanoFactor_fun <- function(expr){
  #calculate Fano factor, choose top 1000 Fano genes
  Fano <- apply(expr,2,function(x) var(x)/mean(x))
  Fano <- Fano %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Gene")
  
  colnames(Fano) <- c("Gene","fano")
  Fano %>% dplyr::filter(!is.na(fano)) %>% dplyr::arrange(desc(fano))
}
raceid_fun <- function(x, mthr = -1){
  uvar  <- function(x,fit){
    err <- coef(summary(fit))[, "Std. Error"]
    2**(coef(fit)[1] + err[1] + log2(x)*(coef(fit)[2] + err[2]) + (coef(fit)[3] + err[3]) * log2(x)**2)
  }
  
  x <- t(x)
  m <- apply(x, 1, mean)
  v <- apply(x, 1, var)
  ml <- log2(m)
  vl <- log2(v)
  f <- ml > -Inf & vl > -Inf
  ml <- ml[f]
  vl <- vl[f]
  mm <- -8
  repeat {
    fit <- lm(vl ~ ml + I(ml^2))
    if (coef(fit)[3] >= 0 | mm >= mthr) {
      break
    }
    mm <- mm + 0.5
    f <- ml > mm
    ml <- ml[f]
    vl <- vl[f]
  }
  vln <- log2(v) - log2(sapply(m, FUN = uvar, fit = fit))
  tibble(
    Gene = names(vln),
    vln = vln
  ) %>%
    dplyr::filter(!is.na(vln)) %>%
    dplyr::arrange(desc(vln))
}

rf.fun <- function(expr, label, index_id){
  use_python('/home/pauling/anaconda3/bin/python')
  rf <- import('sklearn.ensemble')
  RF <- rf$RandomForestClassifier
  
  X_label <- label$label[index_id]
  Y_label <- label$label[-index_id]
  clf <- RF(n_estimators = 200L)
  
  res1 <- SE_fun(t(expr[index_id,]))
  res2 <- m3d_fun(expr[index_id,])
  res3 <- Gini_fun(expr[index_id,])
  res4 <- HVG_fun(expr[index_id,])
  res5 <- sct_fun(expr[index_id,])
  res6 <- FanoFactor_fun(expr[index_id,])
  res7 <- raceid_fun(expr[index_id,])
  
  expr <- log2(expr+1)
  
  tibble(res = list(res1,res2,res3,res4,res5,res6,res7),
         method = c('SE','M3Drop','Gini','HVG','sct','fano','raceid')) %>%
    dplyr::mutate(
      auc = purrr::map(
        .x = res,
        .f = function(.x){
          tmp <- c()
          n <- 0
          for (num in c(30,50,100,200,300,500,1000,5000)) {
            n <- n + 1
            X_train <- as.matrix(expr[index_id, .x$Gene[1:num]])
            X_test <- as.matrix(expr[-index_id, .x$Gene[1:num]])
            clf$fit(X_train, X_label)
            prd <- clf$predict(X_test)
            CohenKappa <- tibble(
              ori = as.character(Y_label),
              prd = prd)
            
            CohenKappa %>%
              dplyr::filter(ori == prd) %>%
              nrow(.) -> accuracy
            
            tmp[n] <- accuracy/nrow(CohenKappa)
          }
          
          tibble(
            num = c(30,50,100,200,300,500,1000,5000),
            acc = tmp
          ) -> tmp
          
          return(tmp)
        }
      )
    ) %>% 
    dplyr::select(-res) -> fin_res
  
  print("done")
  return(fin_res)
}

acc <- list()
sample.index <- list()

for (i in 1:50) {
  cda %>%
    dplyr::mutate(ID = 1:nrow(.)) %>%
    dplyr::sample_frac(0.7) %>%
    dplyr::pull(ID) -> sample.index[[i]]
}
 

for (i in 1:50) {
  acc[[i]] <- rf.fun(expr, cda, sample.index[[i]])
  print(i)
}

Reduce(rbind, acc) %>%
  dplyr::mutate(
    auc = purrr::map2(
      .x = method,
      .y = auc,
      .f = function(.x, .y){
        .y %>%
          dplyr::mutate(method = .x)
      }
    )
  ) -> res

res %>% readr::write_rds("/home/pauling/projects/04_SEmodel/07_NC_revision/01.res/01.featureSelection/01.10x.classification/02.5.cell.line.rds.gz", compress = "gz")