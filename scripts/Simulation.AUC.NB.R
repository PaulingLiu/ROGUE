## Simulation --- AUC

#library
library(SingleCellExperiment);library(tidyverse);library(reticulate);library(tidyverse)
library(ggplot2);library(scmap);library(fmsb);library(ggsci);library(scibet);library(ggpubr)
library(Seurat);library(M3Drop);library(ROCR);library(cluster);library(parallel)

#Functions
simul_da <- function(gene_means, r = 2, n_gene = 10000, n_cell = 2000, ZINB = F){
  sda <- matrix(data=NA, nrow=n_cell, ncol=n_gene, byrow=FALSE, dimnames=NULL)
  gene_means <- gene_means[gene_means > 0]
  u <- median(gene_means)
  for (i in 1:n_gene) {
    p <- gene_means[i]/(gene_means[i]+r)
    tmp <- rnbinom(n=n_cell, prob = 1-p, size = r)
    if(isTRUE(ZINB)){
      x <- -(mean(tmp)/u - 1.5)
      p <- 1/(1+exp(-x))
      n <- ceiling(n_cell*p)
      tmp[sample(n_cell,n)] <- 0
    }
    sda[,i] <- tmp
  }
  
  colnames(sda) <- paste("Gene", 1:ncol(sda), sep = '')
  rownames(sda) <- paste("Cell", 1:nrow(sda), sep = '')
  sda <- as.data.frame(sda)
  sda <- lapply(sda, as.numeric) %>% do.call("data.frame", .)
  return(sda)
}
simul_diff <- function(gene_means, fc, r = 2, n_gene = 10000, n_cell = 2000, n_diff = 200, sub = 0.5, ZINB = F){    
  
  #gene_means <- exp(rnorm(n_gene, 0, sd = 2))
  
  sda1 <- simul_da(gene_means = gene_means, r = r, n_gene = n_gene, n_cell = n_cell, ZINB = ZINB)
  diff1 <- gene_means[1:n_diff]
  #fc <- exp(rnorm(n_diff, mean = 0,sd = 2))
  
  tmp <- tibble(
    mean.expr1 = diff1,
    mean.expr2 = diff1*fc,
    fc = fc,
    Gene = paste("Gene", 1:n_diff, sep = "")
  )
  
  u <- median(gene_means)
  
  simul_expr <- function(.x, .l){
    p <- .x/(.x+r)
    tmp <- rnbinom(.l, prob = 1-p, size = r)
    if(isTRUE(ZINB)){
      x <- -(mean(tmp)/u - 1.5)
      p <- 1/(1+exp(-x))
      n <- ceiling(.l*p)
      tmp[sample(.l, n)] <- 0
    }
    
    return(tmp)
  }
  for (i in 1:nrow(tmp)) {
    expr1 <- simul_expr(tmp[i,]$mean.expr1, .l = ceiling(n_cell*sub))
    expr2 <- simul_expr(tmp[i,]$mean.expr2, .l = n_cell - ceiling(n_cell*sub))
    sda1[,i] <- c(expr1, expr2)
  }
  sda1 <- as.data.frame(sda1)
  sda <- list(sda1, tmp)
  return(sda)
}
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
  hvg.res <- hvg.res %>% dplyr::arrange(p.value)
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
  
  sct.res <- sct.res %>% 
    dplyr::mutate(
      p.value = 1-pnorm(
        sct.res$residual_variance,
        mean = mean(sct.res$residual_variance),
        sd = sd(sct.res$residual_variance)))
  
  return(sct.res)
}
FanoFactor_fun <- function(expr){
  #calculate Fano factor
  Fano <- apply(expr,2,function(x) var(x)/mean(x))
  Fano <- Fano %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Gene")
  
  colnames(Fano) <- c("Gene","fano")
  Fano <- Fano %>% dplyr::filter(!is.na(fano)) %>% dplyr::arrange(desc(fano))
  
  Fano <- Fano %>% 
    dplyr::mutate(
      p.value = 1-pnorm(
        Fano$fano,
        mean = mean(Fano$fano),
        sd = sd(Fano$fano)))
  
  return(Fano)
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
    dplyr::arrange(desc(vln)) -> race.res
  
  race.res <- race.res %>% 
    dplyr::mutate(
      p.value = 1-pnorm(
        race.res$vln,
        mean = mean(race.res$vln),
        sd = sd(race.res$vln)))
  return(race.res)
}

cal_auc <- function(.x, gene){
  .x <- .x %>% dplyr::mutate(diff = ifelse(Gene %in% gene, 0, 1))
  pred <- prediction(.x$p.value, .x$diff)
  perf <- performance(pred,'auc')
  auc <- perf@y.values[[1]]
  return(auc)
}

out.path <- "/home/pauling/projects/04_SEmodel/05_data_phase2/01.gene.selection.benchmark/02.sim.auc.data"

sub <- 0.5  # sub <- 0.2 / 0.1 / 0.01
AUC <- list()
count <- 0

tibble(rep = 1:200) %>%
  dplyr::mutate(
    mean = purrr::map(
      .x = rep,
      .f = function(.x){
        exp(rnorm(10000, 0, sd = 2))
      }
    )
  ) %>%
  dplyr::mutate(
    fc = purrr::map(
      .x = rep,
      .f = function(.x){
        exp(rnorm(200, mean = 0, sd = 2))
      }
    )
  ) -> sda.pre

for (r in c(5, 10, 15, 20)) {
  for (i in 1:50) {
    count <- count + 1
    sda <- simul_diff(sda.pre$mean[[count]], sda.pre$fc[[count]], sub = sub, r = r, ZINB = F)  # ZINB = T

    res1 <- SE_fun(t(sda[[1]]), span = 0.1)
    res2 <- m3d_fun(sda[[1]])
    res3 <- HVG_fun(sda[[1]])
    res4 <- Gini_fun(sda[[1]])
    res5 <- sct_fun(sda[[1]])
    res6 <- FanoFactor_fun(sda[[1]])
    res7 <- raceid_fun(sda[[1]])
    
    diff.gene <- sda[[2]] %>% 
      dplyr::filter(fc <= 1.5 | fc >= 1.5) %>% 
      dplyr::pull(Gene)
    
    auc1 <- cal_auc(res1, diff.gene)
    auc2 <- cal_auc(res2, diff.gene)
    auc3 <- cal_auc(res3, diff.gene)
    auc4 <- cal_auc(res4, diff.gene)
    auc5 <- cal_auc(res5, diff.gene)
    auc6 <- cal_auc(res6, diff.gene)
    auc7 <- cal_auc(res7, diff.gene)
    
    AUC[[count]] <- c(auc1, auc2, auc3, auc4, auc5, auc6, auc7, r, sub)
    print(paste0("r = ",r,", i = ",i))
  }
}

pda <- Reduce(rbind,AUC) %>% as.matrix() %>% as.tibble()
colnames(pda) <- c("SE","M3Drop","HVG","Gini","SCT","Fano","RaceID","r","sub")

pda %>% readr::write_rds(file.path(out.path, paste("02.diff.gene.sim.auc", sub, ".rds.gz", sep = "")),compress = "gz")

auc.plot <- function(pda, r.value){
  pda %>%
    as.tibble() %>%
    dplyr::filter(r == r.value) %>% 
    dplyr::select(-r,-sub) %>%
    tidyr::gather(key = "method", value = "AUC") %>%
    ggplot(aes(factor(method, levels = c("SE","RaceID","SCT","HVG","Gini","Fano","M3Drop")), AUC)) +
    geom_boxplot(aes(colour = method), outlier.shape = NA, lwd = 0.6) +
    geom_jitter(aes(colour = method), width = 0.04, size = 0.8) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 15, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1)) +
    labs(
      x = " ",
      y = "AUC"
    ) +
    scale_colour_d3() -> p_boxplot
  

  return(p_boxplot)
}

p1 <- auc.plot(pda, r.value = 5)
p2 <- auc.plot(pda, r.value = 10)
p3 <- auc.plot(pda, r.value = 15)
p4 <- auc.plot(pda, r.value = 20)
mp <- ggarrange(p1,p2,p3,p4, ncol = 4)

ggsave(filename = "auc.sub.0.1.pdf",
       plot = mp, 
       path = "/home/pauling/projects/04_SEmodel/07_NC_revision/02.figures/03.simulation.auc/01.nb", 
       width = 16, height = 4)
