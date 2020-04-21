
library(SingleCellExperiment);library(tidyverse);library(reticulate);library(tidyverse)
library(ggplot2);library(scmap);library(fmsb);library(ggsci);library(ROGUE)
library(Seurat);library(M3Drop);library(ROCR);library(cluster);library(parallel);library(scran)

# Load data
lung.sce <- readr::read_rds("/data1/pauling/04_SEmodel/01_data/08_lung_Bcells/com.lung.liver.Bcell.rds.gz")
meta <- readr::read_rds("/home/pauling/projects/04_SEmodel/01_data/08_lung_Bcells/meta.lung.liver.Bcell.rds.gz")
expr <- lung.sce@assays$RNA@counts

# OUT path
fig.path <- "/home/pauling/projects/04_SEmodel/07_NC_revision/02.figures/07.Bcell"
out.path <- "/home/pauling/projects/04_SEmodel/07_NC_revision/03.data/07.Bcell"

# Select variable genes using S-E
expr <- as.matrix(expr)
ent.res <- SE_fun(expr)

# Perform clustering with Seurat
lung.se.sce <- RunPCA(lung.sce, verbose = FALSE, features = ent.res$Gene[1:3000])
pca <- lung.se.sce@reductions$pca@cell.embeddings
use_python("/home/heyao/data/tools/basic/anaconda3/bin/python")
anndata = import("anndata",convert=FALSE)
sc = import("scanpy.api",convert=FALSE)
np = import("numpy",convert=FALSE)  
bbknn = import("bbknn", convert=FALSE)
adata.com = anndata$AnnData(X=pca, obs=meta$Donor)
sc$tl$pca(adata.com)

adata.com$obsm$X_pca = pca
bbknn$bbknn(adata.com,batch_key=0)

sc$tl$umap(adata.com, random_state = 3316L)
sc$tl$leiden(adata.com, resolution = 0.8)
umap.com = py_to_r(adata.com$obsm$X_umap)

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata.com$obs$leiden)), 
  UMAP1 = umap.com[,1], 
  UMAP2 = umap.com[,2]) %>%
  dplyr::mutate(Donor = meta$Donor) %>%
  dplyr::mutate(tissue = meta$tissue) %>%
  dplyr::mutate(project = meta$project) -> bb.res

umap.plot <- bb.res %>%
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(colour = leiden.res.1.bbknn), size = 2) +
  theme_bw() +
  theme_void() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("#FF83FA", "#66CDAA", "#AB82FF", "#00C5CD", "#00B2EE", "#FF6A6A", "#D02090"))

tis.dist.plot <- bb.res %>%
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(colour = tissue), size = 2) +
  theme_bw() +
  theme_void() +
  theme(legend.position = "top") +
  scale_colour_manual(values =  c("#00EEEE", "#00B2EE", "#FF69B4", "#FF3E96", "#48D1CC", "#FF8247", "#90EE90"))


# Confusion heatmap
conf.matr <- table(bb.res$leiden.res.1.bbknn, meta$leiden.res.1.bbknn)
conf.matr <- as.matrix(conf.matr)
conf.matr <- conf.matr/rowSums(conf.matr)
conf.matr %>%
  as.data.frame() %>%
  dplyr::rename(SE = Var1, SCT = Var2) %>%
  ggplot(aes(SCT, SE)) +
  geom_tile(aes(fill = Freq)) +
  scale_fill_viridis_c() +
  theme_minimal() -> conf.plot


#ROGUE calculation
bb.res <- bb.res %>% dplyr::mutate(ID = 1:nrow(.))
clusters <- unique(bb.res$leiden.res.1.bbknn)
patient.rogue <- function(info, cluster){
  tmp <- bb.res %>% dplyr::filter(leiden.res.1.bbknn == cluster)
  patients <- unique(bb.res$Donor)
  rogue <- c()
  for (i in 1:length(patients)) {
    print(i)
    index1 <- tmp %>% dplyr::filter(Donor == patients[i]) %>% dplyr::pull(ID)
    if(length(index1) >= 20){
      tmp.matr <- matr[,index1]
      #tmp.matr <- matr.filter(tmp.matr)
      tmp.res <- SE_fun(tmp.matr)
      rogue[i] <- CalculateRogue(tmp.res)
    }
    else{
      rogue[i] <- NA
    }
  }
  return(rogue)
}

res <- list()

for (i in 1:length(clusters)) {
  res[[i]] <- patient.rogue(bb.res, clusters[i])
}

res.tibble <- Reduce(rbind, res) %>% as.matrix() %>% t() %>% as.tibble()
colnames(res.tibble) <- clusters
res.tibble %>%
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
  ) -> rogue.plot


# Marker heatmap
Idents(lung.sce) <- bb.res$leiden.res.1.bbknn
b.markers <- FindAllMarkers(lung.sce, only.pos = TRUE, min.pct = 0.25, test.use = "t")
top10 <- b.markers %>% 
  dplyr::mutate(cluster = as.character(cluster)) %>%
  group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% 
  dplyr::arrange(cluster)

DoHeatmap(lung.sce, top10$gene)
marker.matr <- function(tmp.sce, gene){
  
  matr <- as.matrix(tmp.sce@assays$SCT@scale.data)
  gene <- intersect(gene, rownames(matr))
  matr <- matr[gene,]
  
  tibble(cluster = as.numeric(unique(Idents(tmp.sce)))-1) %>%
    dplyr::mutate(
      expr = purrr::map(
        .x = cluster,
        .f = function(.x){
          tmp.matr <- matr[,Idents(tmp.sce) == .x]
          rowMeans(tmp.matr)
        }
      )
    ) -> tmp.res
  
  matr <- Reduce(rbind, tmp.res$expr)
  matr <- as.data.frame(matr)
  rownames(matr) <- paste("C", tmp.res$cluster, sep = "")
  return(matr)
}
pda <- marker.matr(lung.sce, top10$gene)
pda %>%
  tibble::rownames_to_column(var = "Cluster") %>%
  dplyr::mutate_if(is.numeric, funs((. - mean(.))/sd(.))) %>%
  tidyr::gather(key = "Gene", value = "expr", -Cluster) %>%
  ggplot(aes(Cluster, factor(Gene, levels = rev(unique(top10$gene))))) +
  geom_tile(aes(fill = expr)) +
  #scale_fill_distiller(palette = "Spectral") +
  theme(strip.text.x = element_blank(), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 9), 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13), 
        axis.text.y = element_text(color = "black"), 
        axis.text.x = element_text(color = "black"), 
        panel.background = element_rect(colour = "black", fill = "white"), 
        panel.grid = element_line(colour = "grey", linetype = "dashed"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2)) + 
  #facet_grid(. ~ group, scales = "free", space = "free") + 
  scale_fill_distiller(palette = "RdBu") + 
  labs(x = "", y = "") -> marker.heatmap


# R O/E calculation
meta.sub <- bb.res %>% dplyr::filter(project == "liver")

cluster.table <- table(meta.sub$leiden.res.1.bbknn, meta.sub$tissue)
cluster.table %>%
  as.data.frame() %>%
  dplyr::mutate(
    p.value = purrr::pmap_dbl(
      list(
        .x = Var1,
        .y = Var2,
        .z = Freq
      ),
      .f = function(.x, .y, .z){
        a <- .z
        b <- sum(cluster.table[,.y]) - a
        c <- sum(cluster.table[.x,]) - a
        d <- sum(cluster.table) - a - b - c
        
        #o <- fisher.test(matrix(c(a, b,c, d), ncol = 2), alternative = "greater")
        #o$estimate
        o <- chisq.test(matrix(c(a, b, c, d), ncol = 2))
        oe <- o$observed/o$expected
        oe[1,1]
      }
    )
  ) -> enrich.res

#adj.p.value <- p.adjust(enrich.res$p.value, method = "BH")
#enrich.res <- enrich.res %>% dplyr::mutate(adj.p.value = adj.p.value)
enrich.res %>%
  dplyr::rename(`Ro/e` = p.value) %>%
  #dplyr::mutate(p.value = ifelse(p.value < 1, -1/p.value, p.value)) %>%
  #dplyr::mutate(`-log10(adj.P-value)` = -log10(adj.p.value)) %>%
  ggplot(aes(Var2, Var1, fill = `Ro/e`)) +
  geom_tile(colour = "white", lwd = 0.8) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black" ,angle = 45, hjust = 1)) +
  scale_fill_distiller(palette = "Spectral") -> roe.plot


# Ro/e barplot
tissues <- unique(bb.res$tissue)
a <- bb.res %>% dplyr::filter(leiden.res.1.bbknn == 2) %>% dplyr::filter(project == "liver")  # [== 2] # [== "Lung"] #
a$tissue <- as.factor(a$tissue)
sda <- table(a$Donor, a$tissue) %>% as.data.frame()
sda %>%
  as.data.frame() %>%
  tidyr::spread(key = "Var2", value = "Freq") %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "Var1") -> sda

sda2 <- t(sda) %>% as.data.frame() %>% dplyr::mutate_all(funs((./sum(.))))
rownames(sda2) <- colnames(sda)

tibble(Tissue = tissues[1:5]) %>%
  dplyr::mutate(mean = purrr::map_dbl(.x = Tissue, function(.x){mean(unlist(sda2[.x,]))})) %>%
  dplyr::mutate(sem = purrr::map_dbl(.x = Tissue, function(.x){sd(unlist(sda2[.x,]))/sqrt(5)})) %>%
  dplyr::mutate(ymax = mean+sem) %>%
  dplyr::mutate(ymin = mean-sem) %>%
  ggplot(aes(factor(Tissue, levels = c("Tumor","Lymphnode","Blood", "Normal", "Ascites")), mean)) +
  geom_col(aes(fill = Tissue), colour = "black") +
  scale_fill_manual(values = c("#009ACD", "#EE3A8C", "#EE82EE", "#20B2AA", "#FF7F50")) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 10)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black")) +
  labs(
    x = "",
    y = "Fraction"
  ) +
  ylim(0,0.9) -> liver.c2

rownames(umap.com) <- rownames(lung.sce@reductions$umap@cell.embeddings)
colnames(umap.com) <- colnames(lung.sce@reductions$umap@cell.embeddings)
lung.sce@reductions$umap@cell.embeddings <- umap.com
FeaturePlot(lung.sce, features = c("MKI67","STMN1"), pt.size = 0.8, cols = c("#EEE9BF", "red"))



# Save
bb.res %>% readr::write_rds(file.path(out.path,"01.clustering.res.rds.gz"), compress = "gz")
ggsave(plot = umap.plot, filename = "01.umap.plot.pdf", path = fig.path, width = 6, height = 4.5)
ggsave(plot = tis.dist.plot, filename = "02.tissue.distribution.plot.pdf", path = fig.path, width = 6, height = 5)
ggsave(plot = rogue.plot, filename = "03.rogue.boxplot.pdf", path = fig.path, width = 6, height = 4.5)
ggsave(plot = conf.plot, filename = "04.confution.heatmap.pdf", path = fig.path, width = 6, height = 4.5)
ggsave(plot = marker.heatmap, filename = "05.marker.heatmap.pdf", path = fig.path, width = 8, height = 8)
ggsave(plot = roe.plot, filename = "06.roe.heatmap.pdf", path = fig.path, width = 4, height = 5)
ggsave(plot = liver.c4, filename = "07.liver.c4.pdf", path = fig.path, width = 5, height = 3.7)
ggsave(plot = liver.c2, filename = "08.liver.c2.pdf", path = fig.path, width = 5, height = 3.7)
