library(SingleCellExperiment);library(tidyverse);library(reticulate);library(tidyverse)
library(ggplot2);library(scmap);library(fmsb);library(ggsci);library(ROGUE);library(loomR)
library(Seurat);library(M3Drop);library(ROCR);library(cluster);library(parallel);library(scran)

# Load data
datapath <- "/home/pauling/projects/04_SEmodel/07_NC_revision/03.data/06.fibroblast"
cell.info <- readr::read_csv(file.path(datapath,"Thienpont_Tumors_52k_v4_R_fixed.cellInfo.txt"))
gene <- readr::read_rds("/home/pauling/projects/02_data/09_Gene/coding_gene.rds.gz")

dloom <- connect(file.path(datapath,"fibroblasts"), mode = "r")
a <- cell.info %>% dplyr::filter(stringr::str_detect(ClusterName,"fibroblasts"))
a <- as.data.frame(a)
rownames(a) <- a$CellID


sce <- as.Seurat(dloom)
expr <- sce@assays$RNA@counts

overlap.gene <- intersect(gene$gene_name, rownames(expr))
filt.expr <- expr[overlap.gene,]
a <- a[colnames(filt.expr),]

#############
pauling.theme <- function(poi = "right", size = 12, title.size = 13){
  theme(
    legend.position = poi,
    axis.title = element_text(size = title.size,color="black"),
    axis.text = element_text(size = size,color="black"),
    legend.title = element_text(size = size),
    legend.text = element_text(size = size),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))
}
#############

# OUT path
fig.path <- "/home/pauling/projects/04_SEmodel/07_NC_revision/02.figures/06.fibroblasts"
out.path <- "/home/pauling/projects/04_SEmodel/07_NC_revision/03.data/06.fibroblast"

# Clustering (S-E)
filt.expr <- filt.expr[rownames(sce),]
filt.expr <- as.matrix(filt.expr)
ent.res <- SE_fun(filt.expr)

sce.se <- RunPCA(sce, verbose = FALSE, features = ent.res$Gene[1:3000])
sce.se <- RunUMAP(sce.se, dims = 1:30, verbose = FALSE)
pca <- sce.se@reductions$pca@cell.embeddings

adata.com = anndata$AnnData(X=pca, obs=a$PatientNumber)
sc$tl$pca(adata.com)

adata.com$obsm$X_pca = pca
bbknn$bbknn(adata.com,batch_key=0)

sc$tl$umap(adata.com)
sc$tl$leiden(adata.com, resolution = 1.1)
umap.com = py_to_r(adata.com$obsm$X_umap)

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata.com$obs$leiden)), 
  UMAP1 = umap.com[,1], 
  UMAP2 = umap.com[,2]) -> bb.res2

sce.se <- loadumap(sce.se, umap.com, bb.res2)

new.cluster.ids <- paste0("BF_C",c(2,3,1,4,10,9,0,5,6,7,8))
names(new.cluster.ids) <- levels(sce.se)
sce.se <- RenameIdents(sce.se, new.cluster.ids)
DimPlot(sce.se, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# save
readr::write_rds(sce, path = file.path(out.path,"01.sce.11.clusters.rds.gz"), compress = "gz")

# Clustering (SCT)
sce <- CreateSeuratObject(counts = filt.expr)
sce <- SCTransform(sce, verbose = FALSE, do.center = F)
sce <- RunPCA(sce, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:30, verbose = FALSE)
pca <- sce@reductions$pca@cell.embeddings

use_python("/home/heyao/data/tools/basic/anaconda3/bin/python")
anndata = import("anndata",convert=FALSE)
sc = import("scanpy.api",convert=FALSE)
np = import("numpy",convert=FALSE)  
bbknn = import("bbknn", convert=FALSE)
adata.com = anndata$AnnData(X=pca, obs=a$PatientNumber)
sc$tl$pca(adata.com)

adata.com$obsm$X_pca = pca
bbknn$bbknn(adata.com,batch_key=0)

sc$tl$umap(adata.com)
sc$tl$leiden(adata.com, resolution = 1.1)
umap.com = py_to_r(adata.com$obsm$X_umap)

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata.com$obs$leiden)), 
  UMAP1 = umap.com[,1], 
  UMAP2 = umap.com[,2]) -> bb.res

loadumap <- function(.sce, umap, bb.res){
  rownames(umap) <- rownames(.sce@reductions$umap@cell.embeddings)
  colnames(umap) <- colnames(.sce@reductions$umap@cell.embeddings)
  .sce@reductions$umap@cell.embeddings <- umap
  Idents(.sce) <- bb.res$leiden.res.1.bbknn
  return(.sce)
}

sce <- loadumap(sce, umap.com, bb.res)

## Confustion matrix
conf.matr <- table(bb.res$leiden.res.1.bbknn, bb.res2$leiden.res.1.bbknn)
conf.matr <- as.matrix(conf.matr)
conf.matr <- conf.matr/rowSums(conf.matr)
conf.matr %>%
  as.data.frame() %>%
  dplyr::rename(SE = Var1, SCT = Var2) %>%
  ggplot(aes(factor(SCT, levels = c(2,5,10,1,0,3,6,4,7,8,9)), SE)) +
  geom_tile(aes(fill = Freq)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    x = "SCT",
    y = "SE"
  ) +
  pauling.theme() -> conf.plot

## UMAP plot
my.cols <- c("#9F79EE", "#00CDCD", "#008B8B", "#FF6EB4", "#B452CD", "#7EC0EE", "#D02090", "#F5DEB3", "#8B864E", "#3CB371","#9BCD9B")
umap.plot1 <- DimPlot(sce.se, label = T, pt.size = 2.3, cols = my.cols) + theme_void() + NoLegend()

ori.sce <- sce.se
Idents(ori.sce) <- a$ClusterID
umap.plot2 <- DimPlot(ori.sce, label = T, pt.size = 2.3, cols = my.cols) + theme_void() + NoLegend()


## Cell type-specific genes
fib.markers <- FindAllMarkers(sce.se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "t")

top10 <- fib.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
cell.type.genes.plot <- DoHeatmap(sce.se, features = top10$gene) + NoLegend()

## Marker gene expression
mef2c.plot <- FeaturePlot(sce.se, features = "MEF2C", cols = c("#F0E68C","red"), pt.size = 1.5) + theme_void()
myh11.plot <- FeaturePlot(sce.se, features = "MYH11", cols = c("#F0E68C","red"), pt.size = 1.5) + theme_void()


## ROGUE calculation
rogue.res <- rogue(filt.expr, 
                   labels = bb.res2$leiden.res.1.bbknn, 
                   samples = a$PatientNumber, 
                   platform = "UMI", 
                   remove.outlier.n = 5, 
                   filter = T,
                   min.cell.n = 20)

av.rogue <- c()
for (i in 1:11) {
  tmp.r <- rogue.res[,i]
  tmp.r <- tmp.r[!is.na(tmp.r)]
  av.rogue[i] <- mean(tmp.r)
}

## ROGUE plot
rogue.plot <- rogue.boxplot(rogue.res) + ylim(0.6,0.97)


## pathway analysis
gs <- readr::read_tsv("/data1/pauling/02_data/18_Enrichment/h.all.v6.2.symbols.gmt", col_names = F)
gs <- gs[,-2]
gs <- gs %>% tidyr::gather(key = "pathway", value = "sets", -X1)
gs <- gs %>% dplyr::filter(sets %in% rownames(sce.se))
table(gs$sets) %>% 
  as.data.frame() %>% 
  as.tibble() %>% 
  dplyr::arrange(desc(Freq)) %>% 
  dplyr::filter(Freq == 1) %>% 
  dplyr::pull(Var1) -> unique.gene
gs <- gs %>% dplyr::filter(sets %in% unique.gene)

terms <- unique(gs$X1)
geneSets <- list()
for (i in 1:length(terms)) {
  genesets <- gs %>% dplyr::filter(X1 == terms[i]) %>% dplyr::pull(sets)
  geneSets[[terms[i]]] <- genesets
}

y <- as.matrix(sce.se@assays$SCT@data)
gsva_es <- gsva(y, geneSets, mx.diff=1)

sce.gsva <- CreateSeuratObject(counts = gsva_es)
Idents(sce.gsva) <- Idents(sce.se)
fib.markers <- FindAllMarkers(sce.gsva, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.1, test.use = "t")
top10 <- fib.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% dplyr::mutate(gene = stringr::str_replace_all(gene, "-", "_"))
sce.gsva@assays$RNA@scale.data <- gsva_es

tibble(cluster = levels(sce.se)) %>%
  dplyr::mutate(
    expr = purrr::map(
      .x = cluster,
      .f = function(.x){
        tmp.matr <- gsva_es[,Idents(sce.se) == .x]
        rowMeans(tmp.matr)
      }
    )
  ) -> tmp.res

tmp.res <- Reduce(rbind, tmp.res$expr)
tmp.res <- as.data.frame(tmp.res)
rownames(tmp.res) <- levels(sce.se)
tmp.res <- tmp.res[,top10$gene]
tmp.res %>%
  tibble::rownames_to_column(var = "Cluster") %>%
  dplyr::mutate_if(is.numeric, funs((. - mean(.))/sd(.))) %>%
  tidyr::gather(key = "Gene", value = "expr", -Cluster) %>%
  ggplot(aes(factor(Cluster, levels = unique(top10$cluster)), factor(Gene, levels = rev(unique(top10$gene))))) +
  geom_tile(aes(fill = expr), color = "white", lwd = 0.8) +
  #scale_fill_distiller(palette = "Spectral") +
  theme(strip.text.x = element_blank(), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 9), 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13), 
        axis.text.y = element_text(color = "black"), 
        axis.text.x = element_text(color = "black",angle = 45, hjust = 1), 
        panel.background = element_rect(colour = "black", fill = "white"), 
        panel.grid = element_line(colour = "grey", linetype = "dashed"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2)) + 
  #facet_grid(. ~ group, scales = "free", space = "free") + 
  scale_fill_distiller(palette = "Spectral") + 
  labs(x = "", y = "") -> pathway.plot


#############
ggsave(plot = conf.plot, filename = "01.confution.heatmap.pdf", path = fig.path, width = 6, height = 4.5)
ggsave(plot = rogue.plot, filename = "02.se.rogue.box.pdf", path = fig.path, width = 5, height = 4)
ggsave(plot = umap.plot1, filename = "03.re.clustering.se.umap.pdf", path = fig.path, width = 6, height = 4.5)
ggsave(plot = umap.plot2, filename = "04.ori.label.umap.pdf", path = fig.path, width = 6, height = 4.5)
ggsave(plot = mef2c.plot, filename = "05.mef2c.gene.pdf", path = fig.path, width = 5, height = 4)
ggsave(plot = myh11.plot, filename = "05.myh11.gene.pdf", path = fig.path, width = 5, height = 4)
ggsave(plot = cell.type.genes.plot, filename = "06.cell.type.gene.heatmap.pdf", path = fig.path, width = 8, height = 11)
ggsave(plot = pathway.plot, filename = "07.pathway.plot.pdf", path = fig.path, width = 8, height = 8)
