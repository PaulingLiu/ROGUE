library(ROGUE)
library(tidyverse)
library(ggplot2)


## Plot
pauling.theme <- function(poi = "right", size = 12, title.size = 15){
  theme(
    legend.position = poi,
    axis.title = element_text(size = title.size,color="black"),
    axis.text = element_text(size = size,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))
}


# Tabula Muris
path_matr <- '/data1/pauling/04_SEmodel/01_data/01_positive_datasets/02_10X_data/02_pbmc_facs/matr_facs.rds.gz'
expr <- readr::read_rds(path_matr)
expr <- t(expr)
expr <- as.matrix(expr)
colnames(expr) <- paste0("Cell",1:ncol(expr))

ent.res <- SE_fun(expr, if.adj = F)
y <- CalculateRogue(ent.res, platform = "UMI")
x <- length(ent.res$p.adj[ent.res$p.adj < 0.05])

tibble(number = x, ROGUE = y) -> sig.ds

number.of.genes <- c(50,100,200,500,1000,2000)
tibble(number = number.of.genes) %>%
  dplyr::mutate(
    ROGUE = purrr::map_dbl(
      .x = number,
      .f = function(.x){
        ent.res$p.adj[1:.x] <- 0 
        ent.res$p.adj[(.x+1):nrow(ent.res)] <- 1
        CalculateRogue(ent.res, platform = "UMI")
      }
    )
  ) -> res

res %>%
  dplyr::bind_rows(sig.ds) %>%
  ggplot(aes(number, 1-ROGUE)) +
  geom_point(colour = "#1E90FF", size = 2.5) +
  geom_line(colour = "#1E90FF", lwd = 0.8) +
  theme_bw() +
  pauling.theme() +
  labs(
    x = "Number of genes"
  ) +
  geom_point(aes(x, 1-y), colour = "red", size = 2)

# CRC T cells
matr <- readr::read_rds("/data1/pauling/02_data/02_CRC_Tcell/count_with_label.rds.gz")
matr <- matr %>% dplyr::select(-label)
matr <- as.matrix(matr)
matr <- t(matr)
colnames(matr) <- paste0("Cell",1:ncol(matr))

ent.res <- SE_fun(matr, if.adj = F)
y <- CalculateRogue(ent.res, platform = "full-length")
x <- length(ent.res$p.adj[ent.res$p.adj < 0.05])

tibble(number = x, ROGUE = y) -> sig.ds

number.of.genes <- c(50,100,200,500,1000,2000)
tibble(number = number.of.genes) %>%
  dplyr::mutate(
    ROGUE = purrr::map_dbl(
      .x = number,
      .f = function(.x){
        ent.res$p.adj[1:.x] <- 0 
        ent.res$p.adj[(.x+1):nrow(ent.res)] <- 1
        CalculateRogue(ent.res, platform = "full-length")
      }
    )
  ) -> res

res %>%
  dplyr::bind_rows(sig.ds) %>%
  ggplot(aes(number, 1-ROGUE)) +
  geom_point(colour = "#1E90FF", size = 2.5) +
  geom_line(colour = "#1E90FF", lwd = 0.8) +
  theme_bw() +
  pauling.theme() +
  labs(
    x = "Number of genes"
  ) +
  geom_point(aes(x, 1-y), colour = "red", size = 2)
  
