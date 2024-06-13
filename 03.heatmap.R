knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE, 
  message = FALSE,
  eval = FALSE
)

set.seed(123)
library(here)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(BioNERO)
library(patchwork)
library(cowplot)
library(ggplot2)

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


###### Sample correlation (heatmap)######

coldata_colors <- list(
  Ploidy = setNames(
    c("#0073C2", "#EFC000"), # include colors here
    c("Diploid", "Tetraploid") # include levels of the variable here
  ),
  Strain = setNames(
    c("#20854E", "#E18727","#7876B1","#EE4C97"),
    c("s0013", "s9242", "s9316", "s9346")
  )
)

samplecor_pos <- plot_cor_heatmap(
  log_positive, metadata_positive, 
  row_title = "ESI+", show_colnames = FALSE,
  coldata_colors = coldata_colors,
  show_rownames = FALSE,
  cor_method = "pearson"
)

samplecor_neg <- plot_cor_heatmap(
  log_negative, metadata_negative,
  row_title = "ESI-", show_colnames = FALSE,
  coldata_colors = coldata_colors,
  show_rownames = FALSE,
  cor_method = "pearson"
)

samplecor_gcms <- plot_cor_heatmap(
  log_gcms, metadata_gcms, 
  row_title = "GC-MS", show_colnames = FALSE, 
  show_rownames = FALSE,
  cor_method = "pearson",
  coldata_colors = coldata_colors
  
)

samplecor_plots <- wrap_plots(
  ggplotify::as.ggplot(samplecor_pos), 
  ggplotify::as.ggplot(samplecor_neg),
  ggplotify::as.ggplot(samplecor_gcms),
  #plot_layout(guides = "collect"),
  ncol = 3
)



samplecor_plots
