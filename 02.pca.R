knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE, 
  message = FALSE,
  eval = FALSE
)

set.seed(123)
library(here)
library(SummarizedExperiment)
library(BioNERO)
library(patchwork)
source(here("utils.R"))

###load data####
load(here::here("result_files", "log_positive.rda"))
metadata_positive <- read.csv(here("result_files", "metadata_positive.csv"), row.names = 1, header = TRUE)
metadata_positive[, c("Strain", "Ploidy")] <- lapply(metadata_positive[, c("Strain", "Ploidy")], as.factor)

load(here::here("result_files", "log_negative.rda"))
metadata_negative <- read.csv(here("result_files", "metadata_negative.csv"), row.names = 1, header = TRUE)
metadata_negative[, c("Strain", "Ploidy")] <- lapply(metadata_negative[, c("Strain", "Ploidy")], as.factor)

load(here::here("result_files", "log_gcms.rda"))
metadata_gcms <- read.csv(here("result_files", "metadata_gcms.csv"), row.names = 1, header = TRUE)
metadata_gcms[, c("Strain", "Ploidy")] <- lapply(metadata_gcms[, c("Strain", "Ploidy")], as.factor)


######PCA######
pca_pos <- plot_PCA(
  log_positive, metadata_positive[, "Ploidy", drop = FALSE]
)$data |>
  inner_join(
    metadata_positive |> rownames_to_column("Row.names")
  ) |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Strain, shape = Ploidy)) +
  theme_minimal() +
  ##ggsci::scale_color_jco() +
  scale_color_manual(values = c("s0013"="#20854E",
                                "s9316"="#7876B1", "s9242"="#E18727", "s9346"="#EE4C97")) +
  labs(
    x = "PC1 (48.3%)", y = "PC2 (15.9%)",
    title = "PCA of samples, ESI+"
  ) +
  stat_ellipse(aes(color= Strain), level = 0.95, show.legend = FALSE) 


pca_neg <- plot_PCA(
  log_negative, metadata_negative[, "Ploidy", drop = FALSE]
)$data |>
  inner_join(
    metadata_negative |> rownames_to_column("Row.names")
  ) |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Strain, shape = Ploidy)) +
  theme_minimal() +
  scale_color_manual(values = c("s0013"="#20854E",
                                "s9316"="#7876B1", "s9242"="#E18727", "s9346"="#EE4C97")) +
  labs(
    x = "PC1 (29.5%)", y = "PC2 (18.9%)",
    title = "PCA of samples, ESI-"
  ) +
  stat_ellipse(aes(color= Strain), level = 0.95, show.legend = FALSE) 

pca_gcms <- plot_PCA(
  log_gcms, metadata_gcms[, "Ploidy", drop = FALSE]
)$data |>
  inner_join(
    metadata_gcms |> rownames_to_column("Row.names")
  ) |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Strain, shape = Ploidy)) +
  theme_minimal() +
  scale_color_manual(values = c("s0013"="#20854E",
                                "s9316"="#7876B1", "s9242"="#E18727", "s9346"="#EE4C97")) +
  labs(
    x = "PC1 (17.3%)", y = "PC2 (11.3%)",
    title = "PCA of samples, GC-MS"
  ) +
  stat_ellipse(aes(color= Strain), level = 0.95, show.legend = FALSE) 

pca_gcms

samplepca_plots <- wrap_plots(
  wrap_plots(pca_pos, pca_neg, pca_gcms), 
    plot_layout(guides = "collect")
  
)
samplepca_plots
