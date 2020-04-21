library(ggplot2)
library(tidyverse)


## For each dataset
fig.path <- "/home/pauling/projects/04_SEmodel/07_NC_revision/02.figures/01.Classification/01.10x"
res <- readr::read_rds("/home/pauling/projects/04_SEmodel/07_NC_revision/01.res/01.featureSelection/01.10x.classification/02.5.cell.line.rds.gz")
res <- Reduce(rbind, res$auc)

cols <- c(
  "Matisse" = "#1F77B4", "Flamenco" = "#FF7F0E",
  "ForestGreen" = "#2CA02C", "Punch" = "#D62728",
  "Wisteria" = "#9467BD", "SpicyMix" = "#8C564B",
  "Orchid" = "#E377C2", "Gray" = "#7F7F7F",
  "KeyLimePie" = "#BCBD22", "Java" = "#17BECF"
)

p <- res %>%
  dplyr::mutate(method = ifelse(method == "sct","SCT",method)) %>%
  dplyr::mutate(method = ifelse(method == "fano","Fano",method)) %>%
  dplyr::mutate(method = ifelse(method == "raceid","RaceID",method)) %>%
  ggplot(aes(factor(num), acc)) +
  geom_boxplot(aes(colour = factor(method, levels = c("SE","SCT","Fano","Gini","RaceID","M3Drop","HVG"))), outlier.shape = NA) +
  theme_bw() +
  theme(
    legend.position = 'top',
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 13),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_colour_manual(values = unname(cols)[c(7,6,1,2,5,4,3)]) +
  labs(
    y = "Classification accuracy",
    x = "Number of genes"
  ) +
  guides(col = guide_legend(nrow = 1))


ggsave(plot = p, filename = "02.5.cell.line.pdf", path = fig.path, width = 7.71, height = 4.3)