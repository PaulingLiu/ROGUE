library(tidyverse)
library(ROGUE)
library(ggplot2)

simul_diff <- function(r = 2, n_gene = 10000, n_cell = 2000, n_diff = 200, sub = 0.5, ZINB = F, sd = 2){   
  
  gene_means <- exp(rnorm(n_gene, 0, sd = 2))
  
  sda1 <- simul_da(gene_means = gene_means, r = r, n_gene = n_gene, n_cell = n_cell, ZINB = ZINB)
  diff1 <- gene_means[1:n_diff]
  fc <- exp(rnorm(n_diff, mean = 0,sd = sd))
  
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

sda <- simul_diff(n_diff = 100, n_cell = 2000, sd = 1)
sim.expr <- sda[[1]]
diff2 <- exp(rnorm(500, mean = 0, sd = 2))

tmp <- tibble(
  diff2 = diff2,
  Gene = sample(colnames(sda[[1]])[-c(1:50)], 500)
)

expr3 <- sda[[1]][1:1000,]
r <- 2

simul_expr <- function(.x, .l){
  p <- .x/(.x+r)
  tmp <- rnbinom(.l, prob = 1-p, size = r)
  return(tmp)
}

for (i in 1:nrow(tmp)) {
  a <- mean(sim.expr[,tmp[i,]$Gene])
  expr2 <- simul_expr(a*tmp[i,]$diff2, .l = 1000)
  expr3[,tmp[i,]$Gene] <- expr2
}

sim.expr <- rbind(sim.expr, expr3)
rownames(sim.expr) <- paste("cell", 1:nrow(sim.expr), sep = "")
mat_T <- CreateSeuratObject(counts = t(sim.expr))
mat_T <- SCTransform(mat_T, verbose = FALSE)
mat_T <- RunPCA(mat_T, verbose = FALSE)
mat_T <- RunTSNE(mat_T, dims = 1:45, verbose = F)
#mat_T <- RunUMAP(mat_T, dims = 1:20, verbose = F)
mat_T <- FindNeighbors(mat_T, dims = 1:30, verbose = FALSE)
mat_T <- FindClusters(mat_T, verbose = FALSE, resolution = 0.2)
DimPlot(mat_T, label = TRUE, reduction = "tsne") + NoLegend()

mat_T@reductions$tsne@cell.embeddings %>%
  as.tibble() %>%
  dplyr::mutate(label = mat_T$SCT_snn_res.1.2) %>%
  ggplot(aes(tSNE_1, tSNE_2)) +
  geom_point(aes(colour = label)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#FFBBFF", "#5CACEE", "#FF3E96", "#FF6A6A", "#D15FEE")) -> p1

ggsave("01.5clusters.pdf", plot = p1, path = "/home/pauling/projects/04_SEmodel/04_figures/05_ROGUE/07.silhouette/01.sim1", width = 5.8, height = 4.7, units = "in")

pca.data <- mat_T@reductions$pca@cell.embeddings
dd <- dist(pca.data)

si1 <- summary(silhouette(as.numeric(mat_T$SCT_snn_res.0.001), dd))$avg.width
si2 <- summary(silhouette(as.numeric(mat_T$SCT_snn_res.0.1), dd))$avg.width
si3 <- summary(silhouette(as.numeric(mat_T$SCT_snn_res.1), dd))$avg.width
si4 <- summary(silhouette(as.numeric(mat_T$SCT_snn_res.1.2), dd))$avg.width

clusters <- c('1_cluster',"2_clusters","3_clusters","4_clusters","5_clusters")

sim.expr <- t(sim.expr)

tibble(
  resolution = c(0.001, 0.1, 1, 1.2)
) %>%
  dplyr::mutate(
    rogue = purrr::map(
      .x = resolution,
      .f = function(.x){
        print(.x)
        mat_T <- FindClusters(mat_T, verbose = FALSE, resolution = .x)
        se_cluster <- as.numeric(mat_T$seurat_clusters)
        tmp.res <- rogue(sim.expr, labels = se_cluster, samples = rep("a", length(se_cluster)),platform = "UMI")
        unlist(tmp.res)
      }
    )
  ) -> res

ggplot(aes(1:5, c(NA,si1,si2,si3,si4)), data = NULL) +
  geom_line(colour = '#1E90FF') +
  geom_point(colour = '#1E90FF', size = 2) +
  scale_x_continuous(
    breaks = c(1:5),
    label = clusters
  ) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 0),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 40, hjust = 1)
  ) +
  labs(
    x = " ",
    y = "Silhouette"
  ) -> p2

ggsave("02.silh.pdf", plot = p2, path = "/home/pauling/projects/04_SEmodel/04_figures/05_ROGUE/07.silhouette/01.sim1", width = 5.5, height = 3.5, units = "in")

ggplot(aes(1:5, c(0.4465,mean(res$rogue[[1]]),mean(res$rogue[[2]]),mean(res$rogue[[3]]),mean(res$rogue[[4]]))), data = NULL) +
  geom_line(colour = '#FF3E96') +
  geom_point(colour = '#FF3E96', size = 2) +
  scale_x_continuous(
    breaks = c(1:5),
    label = paste0(1:5,"_clusters")
  ) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 0),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 40, hjust = 1)
  ) +
  labs(
    x = " ",
    y = "ROGUE"
  ) -> p2

ggsave("02.rogue.pdf", plot = p2, path = "/home/pauling/projects/04_SEmodel/04_figures/05_ROGUE/07.silhouette/01.sim1", width = 5.5, height = 3.5, units = "in")

