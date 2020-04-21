library(tidyverse)
library(ROGUE)

out.path <- "/home/pauling/projects/04_SEmodel/01_data/04_rogue/09_simulation_rogue/01.vary.gene/01.NB"

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
simul_diff <- function(r = 2, n_gene = 10000, n_cell = 2000, n_diff = 200, sub = 0.5, ZINB = F){   
  
  gene_means <- exp(rnorm(n_gene, 0, sd = 2))
  
  sda1 <- simul_da(gene_means = gene_means, r = r, n_gene = n_gene, n_cell = n_cell, ZINB = ZINB)
  diff1 <- gene_means[1:n_diff]
  fc <- exp(rnorm(n_diff, mean = 0,sd = 2))
  
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


sub <- 1/3   ##### Different values
prop <- c(0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5)
tibble(prop = rep(prop, 20)) %>%
  dplyr::mutate(
    rogue = purrr::map_dbl(
      .x = prop,
      .f = function(.x){
        print(.x)
        n_diff <- .x*20000
        sda <- simul_diff(r = 6, n_gene = 20000, n_cell = 2000, n_diff = n_diff, sub = sub)
        ent_res <- SE_fun(t(sda[[1]]), span = 0.1)
        CalculateRogue(ent_res, platform = "UMI")
      }
    )
  ) -> res

res %>% readr::write_rds(file.path(out.path,paste0(sub,".rds.gz")), compress = "gz")