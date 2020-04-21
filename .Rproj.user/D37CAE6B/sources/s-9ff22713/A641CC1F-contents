#' Calculate the expression entropy of each gene
#' @name Entropy
#' @usage Entropy(expr, r = 1)
#' @param expr The expression matrix. Rows should be genes and columns should be cells.
#' @param r A small fixed value to avoid log(0) of mean gene expression levels. The default value of r is set to 1, but can also be set to other values such as 0.1 and 0.01.
#'
#' @return A tibble object with three columns 'Gene', 'mean.expr' and 'entropy'.
#' @export
#'
#' @examples
Entropy <- function(expr, r = 1){
  tmp <- log(expr+1)
  entropy <- Matrix::rowMeans(tmp)
  mean.expr <- log(Matrix::rowMeans(expr)+r)

  ent_res <- tibble(
    Gene = rownames(expr),
    mean.expr = mean.expr,
    entropy = entropy
  )

  return(ent_res)
}


#' Fit the relationship between expression entropy and mean gene expression
#' @name entropy_fit
#' @description Fit the relationship between expression entropy and mean gene expression using loess regression.
#' @usage entropy_fit(.x, span = 0.5, mt.method = c("fdr","BH"))
#' @param .x A tibble object returned from the Entropy function.
#' @param span The parameter α which controls the degree of smoothing.
#' @param mt.method The multiple testing method used in p.adjust.
#' @return A tibble object with six columns.
#' @export
#'
#' @examples ent.res <- Entropy(expr)
#' ent.res <- entropy_fit(ent.res, span = 0.3, mt.method = "fdr")
entropy_fit <- function(.x, span = 0.5, mt.method = "fdr"){
  .x <- .x %>% dplyr::filter(is.finite(mean.expr)) %>% dplyr::filter(entropy > 0)
  fit <- loess(entropy~mean.expr, data = .x, span=span)
  prd <- predict(fit, .x$mean.expr)
  .x %>%
    dplyr::mutate(fit = prd) %>%
    dplyr::mutate(ds = fit - entropy) %>%
    dplyr::mutate(pv = 1-pnorm(.$ds, mean = mean(.$ds), sd = sd(.$ds))) %>%
    dplyr::filter(pv > 0.1) -> tmp

  fit <- loess(entropy~mean.expr, data = tmp, span=span)
  prd <- predict(fit, .x$mean.expr)
  .x %>%
    dplyr::mutate(fit = prd) %>%
    dplyr::mutate(ds = fit - entropy) %>%
    dplyr::filter(is.finite(ds)) %>%
    dplyr::mutate(pv = 1-pnorm(.$ds, mean = mean(.$ds), sd = sd(.$ds))) %>%
    dplyr::filter(pv > 0.1) -> tmp

  fit <- loess(entropy~mean.expr, data = tmp, span=span)
  prd <- predict(fit, .x$mean.expr)

  .x %>%
    dplyr::mutate(fit = prd) %>%
    dplyr::mutate(ds = fit - entropy) %>%
    dplyr::filter(is.finite(ds)) -> .x

  .x <- .x %>% dplyr::mutate(p.value = 1-pnorm(.x$ds, mean = mean(.x$ds), sd = sd(.x$ds)))
  p.adj <- p.adjust(.x$p.value, method = mt.method)
  .x <- .x %>% dplyr::mutate(p.adj = p.adj) %>% dplyr::arrange(desc(ds))
}


#' Identify highly informative genes using S-E model
#' @description Use S-E curve to identify highly informative genes.
#' @param expr The expression matrix. Rows should be genes and columns should be cells.
#' @param span The parameter α which controls the degree of smoothing.
#' @param r A small fixed value to avoid log(0) of mean gene expression levels. The default value of r is set to 1, but can also be set to other values such as 0.1 and 0.01.
#' @param mt.method The multiple testing method used in p.adjust.
#' @param if.adj Whether to apply multiple testing method to adjust p.value.
#' @return A tibble object with seven columns:
#' @return * Gene, the gene name.
#' @return * mean.expr, the mean expression levels of genes.
#' @return * entropy, the expected expression entropy from a given mean gene expression.
#' @return * fit, the mean expression levels of genes.
#' @return * ds, the entropy reduction against the null expectation.
#' @return * p.value, the significance of ds.
#' @return * p.adj, adjusted P value.
#'
#' @export
#'
#' @examples ent.res <- SE_fun(expr, span = 0.1, r = 1, mt.method = "fdr")
SE_fun <- function(expr, span = 0.5, r = 1, mt.method = "fdr", if.adj = T){
  ent_res <- ROGUE::Entropy(expr, r = r)
  ent_res <- ROGUE::entropy_fit(ent_res, span = span, mt.method = mt.method)
  if(!isTRUE(if.adj)){
    ent_res <- ent_res %>% dplyr::mutate(p.adj = p.value)
  }
  return(ent_res)
}


#' S-E plot
#' @description Draws a point plot of the relationship between S and E.
#' @param .x A tibble object returned from the SE_fun or entropy_fit function.
#' @param point_size Point size for geom_point.
#' @param geom_line Logical, whether to show the expected expression entropy.
#' @param p.adj Logical, whether to highlight significantly varied genes.
#' @param cutoff The threshold (adjusted P value) for identifying significantly varied genes.
#'
#' @return A ggplot object
#' @export
#'
#' @examples ent.res <- SE_fun(expr, span = 0.1, r = 1, mt.method = "fdr")
#' SEplot(ent.res)
SEplot <- function(.x, point_size = 1, geom_line = T, p.adj = T, cutoff = 0.05){
  if(isFALSE(p.adj)){
    if(geom_line){
      .x %>%
        ggplot(aes(mean.expr, entropy)) +
        geom_point(colour = '#1E90FF', size = point_size) +
        geom_line(aes(mean.expr, fit), lwd = 0.7) +
        theme_bw() +
        theme(
          axis.title = element_text(size = 15,color="black"),
          axis.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 0),
          legend.text = element_text(size = 0),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black")
        ) +
        labs(
          x = "log(mean expression)",
          y = "expression entropy"
        ) -> p
    }
    else{
      .x %>%
        ggplot(aes(mean.expr, entropy)) +
        geom_point(colour = '#1E90FF', size = point_size) +
        #geom_line(aes(mean.expr, fit), lwd = 0.7) +
        theme_bw() +
        theme(
          axis.title = element_text(size = 15,color="black"),
          axis.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 0),
          legend.text = element_text(size = 0),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black")
        ) +
        labs(
          x = "log(mean expression)",
          y = "expression entropy"
        ) -> p
    }
  }
  if(isTRUE(p.adj)){
    .x <- .x %>% dplyr::mutate(sig = ifelse(p.adj <= cutoff, 1, 0))

    if(geom_line){
      .x %>%
        ggplot(aes(mean.expr, entropy)) +
        geom_point(aes(colour = factor(sig)), size = point_size) +
        geom_line(aes(mean.expr, fit), lwd = 0.7) +
        scale_color_manual(values = c("#1E90FF", "red")) +
        theme_bw() +
        theme(
          legend.position = "none",
          axis.title = element_text(size = 15,color="black"),
          axis.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 0),
          legend.text = element_text(size = 0),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black")
        ) +
        labs(
          x = "log(mean expression)",
          y = "expression entropy"
        ) -> p
    }

    else{
      .x %>%
        ggplot(aes(mean.expr, entropy)) +
        geom_point(aes(colour = factor(sig)), size = point_size) +
        #geom_line(aes(mean.expr, fit), lwd = 0.7) +
        scale_color_manual(values = c("#1E90FF", "red")) +
        theme_bw() +
        theme(
          legend.position = "none",
          axis.title = element_text(size = 15,color="black"),
          axis.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 0),
          legend.text = element_text(size = 0),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black")
        ) +
        labs(
          x = "log(mean expression)",
          y = "expression entropy"
        ) -> p
    }
  }
  return(p)
}


#' Filtering out low-abundance genes and low-quality cells
#'
#' @param expr The expression matrix. Rows should be genes and columns should be cells.
#' @param min.cells Include genes detected in at least this many cells.
#' @param min.genes Include cells where at least this many genes are detected.
#'
#' @return A filtered gene expression matrix.
#' @export
#'
#' @examples expr <- matrix(data = rbinom(n = 100, size = 20, prob = 0.5 ), nrow = 10)
#' expr
#' filtered.expr <- matr.filter(expr, min.cells = 3, min.genes = 3)
#' filtered.expr
matr.filter <- function(expr, min.cells = 10, min.genes = 10){
  gene_count <- colSums(expr > 0, na.rm = T)
  cell_count <- rowSums(expr > 0, na.rm = T)

  lq1 <- cell_count < min.cells
  lq2 <- gene_count < min.genes

  return(expr[!lq1, !lq2])
}


#' ROGUE calculation
#' @description  Using ROGUE to assess the purity of single cell population.
#' @param .x A tibble object returned from the SE_fun or entropy_fit function.
#' @param platform The platform ("UMI" or "full-length") used for generating the tested dataset.
#' @param cutoff The threshold (adjusted P value) for identifying significant ds. The default threshold is 0.05.
#' @param k The scaling factor for calculating ROGUE. The default value of K is set to 45 and 500 for droplet-based ("UMI") and "full-length" based datasets, respectively. When specifying a custom k value, the "platform" argument is redundant.
#' @param features Use these features to calculate ROGUE.
#' @details By taking advantage of the wide applicability of S-E model to scRNA-seq data, we introduce the statistic ROGUE to measure the purity of a cell population as:
#' @details \deqn{ROGUE=1-∑sig.ds/(∑sig.ds+K)}
#' @details where K is an important parameter that constrains the ROGUE value between 0 and 1. A cell population with no significant ds for all genes will receive a ROGUE value of 1, while a population with maximum summarization of significant ds is supposed to yield a purity score of ~0.
#' @return A value of ROGUE.
#' @export
#'
#' @examples ent.res <- SE_fun(expr, span = 0.1, r = 1, mt.method = "fdr")
#' CalculateRogue(ent.res, platform = "UMI")
#' CalculateRogue(ent.res, k = 30)
CalculateRogue <- function(.x, platform = NULL, cutoff = 0.05, k = NULL, features = NULL){
  if(is.null(k)){
    if(is.null(platform)){
      warning("Please provide a \"platform\" argument or specify a k value")
    }else if(platform == "UMI"){
      k = 45
    }else if(platform == "full-length"){
      k = 500
    }else if(!is.null(platform) & !(platform %in% c("UMI","full-length"))){
      warning("Please provide valid \"platform\" argument")
    }
  }else if(!is.null(k)){
    k <- k
  }

  if(!is.null(features)){
    .x <- .x %>% dplyr::filter(Gene %in% features)
    sig_value <- sum(abs(.x$ds))
    Rogue <- 1-sig_value/(sig_value+k)
    return(Rogue)
  }else{
    sig_value <- abs(.x$ds[.x$p.adj < cutoff & .x$p.value < cutoff])
    sig_value <- sum(sig_value)
    Rogue <- 1-sig_value/(sig_value+k)
    return(Rogue)
  }
}


#' Remove outlier cells when calculating ROGUE
#' @usage ent.toli(ent, expr, n = 2, span = 0.5, r = 1, mt.method = c("fdr","BH"))
#' @param ent A tibble object returned from the SE_fun or entropy_fit function.
#' @param expr The expression matrix. Rows should be genes and columns should be cells.
#' @param n Remove this many outlier cells.
#' @param span The parameter α which controls the degree of smoothing.
#' @param r A small fixed value to avoid log(0) of mean gene expression levels. The default value of r is set to 1, but can also be set to other values such as 0.1 and 0.01.
#' @param mt.method The multiple testing method used in p.adjust.
#'
#' @return A tibble object with seven columns as 'ent' object.
#' @export
#'
#' @examples ent.toli(ent.res, expr, n = 2, mt.method = "fdr")
ent.toli <- function(ent, expr, n = 2, span = 0.5, r = 1, mt.method = "fdr"){
  sig.gene <- ent %>% dplyr::filter(p.adj < 0.05) %>% dplyr::pull(Gene)
  ng <- length(sig.gene)
  expr <- expr[sig.gene,]

  mean.v <- c()
  entr.v <- c()
  for (i in 1:ng) {
    .x <- as.numeric(expr[i,])
    .x <- base::sort(.x, decreasing = T)
    .x <- .x[-c(1:n)]
    mean.v[i] <- log(mean(.x)+r)
    entr.v[i] <- mean(log(.x+1))
  }

  mean.cut <- min(ent$mean.expr)

  ent$mean.expr[1:ng] <- mean.v
  ent$entropy[1:ng] <- entr.v

  ent <- ent %>% dplyr::select(-p.adj) %>% dplyr::filter(mean.expr > mean.cut)
  ent <- entropy_fit(ent, span = span, mt.method = "fdr")
  return(ent)
}


#' Calculate the ROGUE value of each putative cluster for each sample.
#' @usage rogue(expr, labels, samples, platform = NULL, k= NULL, min.cell.n = 10, remove.outlier.n = 2, span = 0.5, r = 1, mt.method = c("fdr","BH"))
#' @param expr The expression matrix. Rows should be genes and columns should be cells.
#' @param labels A vector of cell cluster lables for all cells corresponding to 'expr' argument.
#' @param samples A vector of samples (e.g. patients) to which each cell belongs, corresponding to 'expr' argument.
#' @param min.cell.n Only clusters containing at least this many cells will receive ROGUE values.
#' @param remove.outlier.n Remove this many outlier cells when calculating ROGUE.
#' @param span The parameter α which controls the degree of smoothing.
#' @param r A small fixed value to avoid log(0) of mean gene expression levels. The default value of r is set to 1, but can also be set to other values such as 0.1 and 0.01.
#' @param filter Logical, whether to filter out low-abundance genes and low-quality cells.
#' @param min.cells if parameter filter is "TRUE", include genes detected in at least this many cells.
#' @param min.genes if parameter filter is "TRUE", Include cells where at least this many genes are detected.
#' @param mt.method The multiple testing method used in p.adjust.
#'
#' @return A dataframe where rows represent samples, cols represent clusters, and values represent corresponding ROGUEs.
#' @export
#'
#' @examples
rogue <- function(expr, labels, samples, platform = NULL, k = NULL, min.cell.n = 10, remove.outlier.n = 2, span = 0.5, r = 1, filter = F, min.cells = 10, min.genes = 10, mt.method = "fdr"){
  clusters <- unique(labels)
  meta <- tibble(CellID = 1:ncol(expr), ct = labels, sample = samples)
  sample.rogue <- function(meta, cluster){
    tmp <- meta %>% dplyr::filter(ct == cluster)
    s <- unique(samples)
    rogue <- c()
    for (i in 1:length(s)) {
      index1 <- tmp %>% dplyr::filter(sample == s[i]) %>% dplyr::pull(CellID)
      if(length(index1) >= min.cell.n){
        tmp.matr <- expr[,index1]
        if(isTRUE(filter)){
          print("Filtering out low-abundance genes and low-quality cells")
          tmp.matr <- matr.filter(tmp.matr, min.cells = min.cells, min.genes = min.genes)
        }else{
          tmp.matr <- tmp.matr
        }
        tmp.res <- SE_fun(tmp.matr, span = span, r = r)
        tmp.res <- ent.toli(tmp.res, tmp.matr, span = span, r = r, n = remove.outlier.n)
        rogue[i] <- CalculateRogue(tmp.res, platform = platform, k = k)
      }
      else{
        rogue[i] <- NA
      }
    }
    return(rogue)
  }

  res <- list()

  for (i in 1:length(clusters)) {
    res[[i]] <- sample.rogue(meta, clusters[i])
  }

  res.tibble <- Reduce(rbind, res) %>% as.matrix() %>% t() %>% as.data.frame()
  colnames(res.tibble) <- clusters
  rownames(res.tibble) <- unique(samples)
  return(res.tibble)
}


#' Visualize ROGUE values on a boxplot
#' @description Draws a boxplot of the ROGUE values for each cluster in different samples.
#' @param res.rogue A dataframe returned from the 'rogue' function.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples res.rogue <- rogue(expr, labels, samples)
#' rogue.boxplot(res.rogue)
rogue.boxplot <- function(res.rogue){
  res.rogue %>%
    tidyr::gather(key = clusters, value = ROGUE) %>%
    ggplot(aes(clusters, ROGUE)) +
    geom_boxplot(color = "#FF3E96", outlier.shape = NA) +
    geom_point(color = "#FF3E96", size = 1.5) +
    theme_bw() +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 13, colour = "black")) +
    labs(
      x = "Clusters",
      y = "ROGUE"
    ) -> p
  return(p)
}


#' Calculate the value of the reference factor K
#' @description Determine the value of the reference factor K by specifying a highly heterogeneous dataset
#' @param expr The expression matrix. Rows should be genes and columns should be cells.
#' @param span The parameter α which controls the degree of smoothing.To improve the performance.
#' @param r A small fixed value to avoid log(0) of mean gene expression levels. The default value of r is set to 1, but can also be set to other values such as 0.1 and 0.01.
#' @param mt.method The multiple testing method used in p.adjust.
#' @param if.adj Whether to apply multiple testing method to adjust p.value.
#' @return The K value
#'
#' @export
#'
#' @examples k <- DetermineK(expr, span = 0.5, r = 1, mt.method = "fdr")
DetermineK <- function(expr, span = 0.5, r = 1, mt.method = "fdr", if.adj = T){
  ent_res <- ROGUE::Entropy(expr, r = r)
  ent_res <- ROGUE::entropy_fit(ent_res, span = span, mt.method = mt.method)
  if(!isTRUE(if.adj)){
    ent_res <- ent_res %>% dplyr::mutate(p.adj = p.value)
  }
  k <- ent_res %>% dplyr::filter(p.adj < 0.05) %>% dplyr::pull(ds) %>% sum()
  k <- k/2
  return(k)
}
