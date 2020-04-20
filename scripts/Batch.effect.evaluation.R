library(tidyverse)
library(ROGUE)
library(ggplot2)

# Load data
matr1 <- readr::read_rds("/home/pauling/projects/04_SEmodel/01_data/04_rogue/10_batch/matr1_con.rds.gz")
matr2 <- readr::read_rds("/home/pauling/projects/04_SEmodel/01_data/04_rogue/10_batch/matr2_sti.rds.gz")
mda <- rbind(matr1, matr2)

info <- readr::read_tsv("/home/pauling/projects/04_SEmodel/01_data/04_rogue/10_batch/nbt.4096-S3.txt")
info <- info %>% dplyr::filter(Cell_Name %in% rownames(mda))
info1 <- info %>% dplyr::filter(Stim_Condition == "ctrl")
info2 <- info %>% dplyr::filter(Stim_Condition == "stim")

matr1 <- t(matr1[info1$Cell_Name,])
matr2 <- t(matr2[info2$Cell_Name,])
mda <- t(mda[info$Cell_Name,])


###############################################################################
###                                                                         ###
###                         Batch ~ interferon-beta                         ###
###                                                                         ###  
###############################################################################


# ROGUE calculation
rogue.res1 <- rogue(
  matr1,
  labels = info1$Cluster_ID,
  samples = info1$Patient,
  min.cell.n = 20,
  platform = "UMI",
  remove.outlier.n = 0
)

rogue.res2 <- rogue(
  matr2,
  labels = info2$Cluster_ID,
  samples = info2$Patient,
  min.cell.n = 20,
  platform = "UMI",
  remove.outlier.n = 0
)

rogue.res.comb <- rogue(
  mda,
  labels = info$Cluster_ID,
  samples = info$Patient,
  min.cell.n = 20,
  platform = "UMI",
  remove.outlier.n = 0
)

# Group
rogue.res1 <- rogue.res1 %>% dplyr::mutate(Condition = "Batch 1 (ctrl)")
rogue.res2 <- rogue.res2 %>% dplyr::mutate(Condition = "Batch 2 (stim)")
rogue.res.comb <- rogue.res.comb %>% dplyr::mutate(Condition = "Aggregated")

rogue.res <- rogue.res1 %>% dplyr::bind_rows(rogue.res2, rogue.res.comb)


# Boxplot
fig.path <- "/home/pauling/projects/04_SEmodel/07_NC_revision/02.figures/09.batch/01.pbmc/01.stimulation"
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

rogue.res$Condition <- factor(rogue.res$Condition, levels = c("Batch 1 (ctrl)","Batch 2 (stim)","Aggregated"))

#rogue.res %>% readr::write_rds("/home/pauling/projects/04_SEmodel/07_NC_revision/03.data/08.batch/rogue.rds.gz", compress = "gz")

cell.types <- c("CD14 Mono","CD16 Mono","CD4 Memory T","CD4 Naive T","CD8 T", "DC")
for (i in 1:6) {
  rogue.res %>%
    tidyr::gather(key = "CellTypes", value = "ROGUE", -Condition) %>%
    dplyr::filter(CellTypes == cell.types[i]) %>%
    ggplot(aes(Condition, ROGUE, colour = Condition)) +
    geom_boxplot(aes(fill = Condition), alpha = 0.5) +
    geom_point(aes(colour = Condition)) +
    theme_classic() +
    pauling.theme(poi = "NULL") +
    labs(x = "") +
    scale_colour_manual(values = c("#00688B", "#008B8B", "#CD6839")) +
    scale_fill_manual(values = c("#00688B", "#008B8B", "#CD6839")) +
    labs(
      title = cell.types[i]
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 30, hjust = 1)) -> p
  
  ggsave(plot = p, filename = paste0(i,".",cell.types[i],".pdf"), path = fig.path, width = 3.3, height = 4.39)
}


# t.test

pv <- list()
for (i in 1:6) {
  t1 <- t.test(rogue.res1[,cell.types[i]], rogue.res.comb[,cell.types[i]], alternative = "greater")
  t2 <- t.test(rogue.res2[,cell.types[i]], rogue.res.comb[,cell.types[i]], alternative = "greater")
  
  p1 <- t1$p.value
  p2 <- t2$p.value

  pv[[i]] <- c(p1,p2)
}

pv <- Reduce(rbind, pv)
pv <- as.matrix(pv)

colnames(pv) <- c("P.value-1","P.value-2")
rownames(pv) <- use.cell.types


###############################################################################
###                                                                         ###
###                           Batch ~ Patient                               ###
###                                                                         ###  
###############################################################################
fig.path2 <- "/home/pauling/projects/04_SEmodel/07_NC_revision/02.figures/09.batch/01.pbmc/02.individual"

cell.types <- c("CD14 Mono","CD16 Mono","CD4 Memory T","CD4 Naive T","CD8 T", "DC")
patient.number <- length(unique(info1$Patient))

complete.rogue <- list()
for(i in 1:length(cell.types)){
  rep.rogue <- c()
  print(cell.types[i])
  for (n in 1:20) {
    info1 %>%
      dplyr::filter(Cluster_ID == cell.types[i]) %>%
      dplyr::group_by(Patient) %>%
      dplyr::sample_frac(1/patient.number) %>%
      dplyr::pull(Cell_Name) -> cell.names
    
    tmp.matr <- matr1[,cell.names]
    tmp.ent <- SE_fun(tmp.matr)
    rep.rogue[n] <- CalculateRogue(tmp.ent, platform = "UMI")
    print(rep.rogue[n])
  }
  complete.rogue[[i]] <- tibble(rogue = rep.rogue, cell.type = cell.types[i])
}

complete.rogue <- Reduce(rbind, complete.rogue)
complete.rogue <- complete.rogue %>% dplyr::mutate(Condition = "Aggregated")

rogue.ind <- rogue.res1 %>% 
  dplyr::mutate(Condition = "Individual") %>% 
  tidyr::gather(key = "cell.type", value = "rogue", -Condition)

rogue.mg <- rbind(rogue.ind, complete.rogue)
rogue.mg$Condition <- factor(rogue.mg$Condition, levels = c("Individual","Aggregated"))

for (i in 1:6) {
  rogue.mg %>%
    dplyr::filter(cell.type == cell.types[i]) %>%
    dplyr::rename(ROGUE = rogue) %>%
    ggplot(aes(Condition, ROGUE, colour = Condition)) +
    geom_boxplot(aes(fill = Condition), alpha = 0.5) +
    geom_point(aes(colour = Condition)) +
    theme_classic() +
    pauling.theme(poi = "NULL") +
    labs(x = "") +
    scale_colour_manual(values = c("#00688B", "#CD6839")) +
    scale_fill_manual(values = c("#00688B", "#CD6839")) +
    labs(
      title = cell.types[i]
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 30, hjust = 1)) -> p
  
  ggsave(plot = p, filename = paste0(i,".",cell.types[i],".pdf"), path = fig.path2, width = 3.3, height = 4.39)
}

