library(ROGUE)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)

fig.path <- "/home/pauling/projects/04_SEmodel/07_NC_revision/02.figures/10.brain"
my.co <- c("#FF8C00", "#D02090", "#FFE7BA", "#00B2EE","#00C5CD",  "#B0E2FF", "#FF83FA","#9FB6CD", "#AB82FF", "#FFC1C1","#00868B","#FF4040", "#3CB371", "#CD9B9B")

load("/data1/et/muris/facs/FACS_all.Robj")
tiss_FACS@meta.data %>%
  tibble::rownames_to_column(var = "barcode") %>%
  dplyr::filter(tissue == "Brain_Non-Myeloid") %>%
  dplyr::rename(label = cell_ontology_class) %>%
  dplyr::select(barcode, tissue, label, mouse.id) -> cda

expr <- tiss_FACS@data
expr <- expr[,cda$barcode]
expr <- exp(expr)-1
expr <- as.matrix(expr)

row.expr <- tiss_FACS@raw.data
row.expr <- row.expr[,cda$barcode]
row.expr <- as.matrix(row.expr)

rogue.brain <- rogue(
  expr,
  labels = cda$label,
  samples = cda$mouse.id,
  min.cell.n = 15,
  platform = "full-length",
  remove.outlier.n = 0
)

brain.clusters <- colnames(rogue.brain)
brain.clusters <- brain.clusters[c(7,5,6,1,3,2,4)]  ## Order

rogue.brain %>%
  tidyr::gather(key = clusters, value = ROGUE) %>% 
  ggplot(aes(factor(clusters, levels = brain.clusters), ROGUE)) + 
  geom_boxplot(color = "#FF3E96", outlier.shape = NA) + 
  geom_point(color = "#FF3E96", size = 1.5) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12,colour = "black"), 
        axis.title = element_text(size = 13,colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 40, hjust = 1)) + 
  labs(x = "Clusters", y = "ROGUE") -> p

ggsave(plot = p, filename = "01.overall.rogue.pdf", path = fig.path, width = 3.97, height = 5.1)

cda.oligo <- cda %>% dplyr::filter(label == "oligodendrocyte")
expr.oligo <- expr[,cda$barcode]

#################################### Clustering ################################################

# Feature selection

genes <- rownames(expr.oligo)
genes <- stringr::str_replace_all(genes, "_", "-")
rownames(expr.oligo) <- genes

ent.res <- SE_fun(expr.oligo)    ## Feature selection with S-E

cda <- as.data.frame(cda)
rownames(cda) <- cda$barcode

sce <- CreateSeuratObject(counts = expr.oligo, meta.data = cda)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1000000)
sce <- ScaleData(sce, features = genes)
sce <- RunPCA(sce, features = ent.res$Gene[1:1000])
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, label = T)


rogue.oligo <- rogue(
  expr.oligo,
  labels = Idents(sce),
  samples = sce$mouse.id,
  min.cell.n = 15,
  platform = "full-length",
  remove.outlier.n = 0
)


rogue.oligo %>% 
  tidyr::gather(key = clusters, value = ROGUE) %>% 
  dplyr::mutate(clusters = stringr::str_replace_all(clusters,"Cluster","Cluster ")) %>%
  ggplot(aes(clusters, ROGUE)) + 
  geom_boxplot(color = "#FF3E96",outlier.shape = NA) + 
  geom_point(color = "#FF3E96", size = 1.5) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12,colour = "black"), 
        axis.title = element_text(size = 13,colour = "black"),
        axis.text.x = element_text(size = 12,colour = "black", angle = 45, hjust = 1)) + 
  labs(x = "", y = "ROGUE") -> p

ggsave(plot = p, filename = "05.reclustering.rogue.pdf", path = fig.path, width = 5.67, height = 5.1)

################################

sce %>% readr::write_rds("/home/pauling/projects/04_SEmodel/07_NC_revision/03.data/09.brain/01.sce.rds.gz", compress = "gz")
new.cluster.ids <- paste0("Cluster", 1:length(levels(sce)))
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, label = F, pt.size = 1) + scale_colour_manual(values = my.co) + theme_void() + NoLegend()

markers <- FindAllMarkers(sce, logfc.threshold = 0.25, min.pct = 0.25,test.use = "t")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(sce, features = top10$gene, size = 4) + NoLegend() 
markers %>% readr::write_rds("/home/pauling/projects/04_SEmodel/07_NC_revision/03.data/09.brain/02.markers.rds.gz", compress = "gz")

markers <- markers %>% dplyr::distinct(gene, .keep_all = T)
entr <- bitr(markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb='org.Mm.eg.db')
markers <- markers %>% dplyr::inner_join(entr, by = c("gene" = "SYMBOL"))


kegg.fun <- function(.x){
  cluster1 <- markers %>% dplyr::filter(cluster == .x)
  
  enrich.pathway <- enrichKEGG(
    cluster1$ENTREZID,
    organism = "mmu",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )
  
  enrich.pathway@result %>% dplyr::arrange(p.adjust) -> res
  return(res)
}

res <- kegg.fun("Cluster5")

res %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(factor(Description, levels = rev(res$Description)), -log10(p.adjust))) +
  geom_col(aes(fill = -log10(p.adjust))) +
  theme_bw() +
  scale_fill_distiller(palette = "Spectral") +
  coord_flip() +
  labs(
    x = " ",
    y = "-log10(p.value)"
  ) +
  theme(
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )


expr <- sce@assays$RNA@counts
expr <- t(expr)
expr <- 1000000*expr/rowSums(expr)
expr <- as.matrix(expr)
expr <- as.data.frame(expr)
expr <- expr %>% dplyr::mutate(label = Idents(sce))
etest.gene <- scibet::SelectGene(expr, k = 50, r = F)
Marker_heatmap(expr = expr, gene = etest.gene)
