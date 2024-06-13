knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE, 
  message = FALSE,
  eval = FALSE
)

set.seed(123)
library(here)
library(tidyverse)
library(SummarizedExperiment)
library(QFeatures)
library(limma)
library(BioNERO)
library(patchwork)

source(here("utils.R"))

log_negative <- read.csv(here("result_files","log_negative.csv"), header=T, row.names = 1, sep = ',')
metadata_negative <- read.csv(here("result_files","metadata_negative.csv"), row.names = 1, sep = ',', header=T)

log_positive <- read.csv(here("result_files","log_positive.csv"), header=T, row.names = 1, sep = ',')
metadata_positive <- read.csv(here("result_files","metadata_positive.csv"), row.names = 1, sep = ',', header=T)

log_gcms <-read.csv(here("result_files", "log_gcms.csv"), header=T, row.names = 1, sep = ',')
metadata_gcms <- read.csv(here("result_files","metadata_gcms.csv"), row.names = 1, sep = ',', header=T)


### pair-wise

###1. Diploid vs Tetraploid
###2. s9242 diploid VS s9242 tetraploid
###3. s9346 diploid VS s9346 tetraploid
###4. s9316 diploid VS s9316 tetraploid
###5. s0013 diploid VS s0013 tetraploid

## Reshape metadata and create design matrices
metadata_da <- metadata_positive |>
  mutate(
    Ploidy = factor(Ploidy, levels = c("Diploid", "Tetraploid")),
    PS = as.factor(str_c(Ploidy, Strain, sep = "_"))
  )

###check colnames
rownames(metadata_da) == colnames(log_negative)
rownames(metadata_da) == colnames(log_positive)
rownames(metadata_da) == colnames(log_gcms)


## Create design matrices for each contrast
d1 <- metadata_da |> mutate(Ploidy = fct_relevel(Ploidy, "Diploid"))
d1 <- model.matrix(~d1$Ploidy) # coef = 2

d2 <- metadata_da |> mutate(PS = fct_relevel(PS, "Diploid_s9242"))
d2 <- model.matrix(~d2$PS) # coef = 6

d3 <- metadata_da |> mutate(PS = fct_relevel(PS, "Diploid_s9346"))
d3 <- model.matrix(~d3$PS) # coef = 8

d4 <- metadata_da |> mutate(PS = fct_relevel(PS, "Diploid_s9316"))
d4 <- model.matrix(~d4$PS) # coef = 7

d5 <- metadata_da |> mutate(PS = fct_relevel(PS, "Diploid_s0013"))
d5 <- model.matrix(~d5$PS) # coef = 5


#### Perform differential abundance analyses for each contrast

C1_negative = lmFit(log_negative, d1) |> eBayes() |> topTable(coef = 2,  sort.by = "P", number = Inf)

length(which(C1_negative$adj.P.Val < 0.05))
length(which(C1_negative$logFC <1 & C1_negative$logFC > -1))

write.csv(C1_negative, file=here("DAFs", "FC_data","C1_negative.csv"))

C1_positive = lmFit(log_positive, d1) |> eBayes() |> topTable(coef = 2, sort.by = "P",  number = Inf)
write.csv(C1_positive, file=here("DAFs", "FC_data", "C1_positive.csv"))

s9242_negative = lmFit(log_negative, d2) |> eBayes() |> topTable(coef = 6, sort.by = "P", number = Inf)
write.csv(s9242_negative, file=here("DAFs", "FC_data", "s9242_negative.csv"))

s9242_positive = lmFit(log_positive, d2) |> eBayes() |> topTable(coef = 6, sort.by = "P", number = Inf)
write.csv(s9242_positive, file=here("DAFs", "FC_data", "s9242_positive.csv"))

s9346_negative = lmFit(log_negative, d3) |> eBayes() |> topTable(coef = 8, sort.by = "P", number = Inf)
write.csv(s9346_negative, file=here("DAFs", "FC_data","s9346_negative.csv"))

s9346_positive = lmFit(log_positive, d3) |> eBayes() |> topTable(coef = 8, sort.by = "P", number = Inf)
write.csv(s9346_positive, file=here("DAFs", "FC_data", "s9346_positive.csv"))

s9316_negative = lmFit(log_negative, d4) |> eBayes() |> topTable(coef = 7, sort.by = "P", number = Inf)
write.csv(s9316_negative, file=here("DAFs", "FC_data", "s9316_negative.csv"))

s9316_positive = lmFit(log_positive, d4) |> eBayes() |> topTable(coef = 7, sort.by = "P", number = Inf)
write.csv(s9316_positive, file=here("DAFs", "FC_data", "s9316_positive.csv"))

s0013_negative = lmFit(log_negative, d5) |> eBayes() |> topTable(coef = 5, sort.by = "P", number = Inf)
write.csv(s0013_negative, file=here("DAFs", "FC_data", "s0013_negative.csv"))

s0013_positive = lmFit(log_positive, d5) |> eBayes() |> topTable(coef = 5, sort.by = "P", number = Inf)
write.csv(s0013_positive, file=here("DAFs", "FC_data", "s0013_positive.csv"))

C1_gcms = lmFit(log_gcms, d1) |> eBayes() |> topTable(coef = 2,  sort.by = "P", number = Inf)
write.csv(C1_gcms, file=here("DAFs", "FC_data","C1_gcms.csv"))

s9242_gcms = lmFit(log_gcms, d2) |> eBayes() |> topTable(coef = 6, sort.by = "P", number = Inf)
write.csv(s9242_gcms, file=here("DAFs", "FC_data", "s9242_gcms.csv"))

s9346_gcms = lmFit(log_gcms, d3) |> eBayes() |> topTable(coef = 8, sort.by = "P", number = Inf)
write.csv(s9346_gcms, file=here("DAFs", "FC_data", "s9346_gcms.csv"))

s9316_gcms = lmFit(log_gcms, d4) |> eBayes() |> topTable(coef = 7, sort.by = "P", number = Inf)
write.csv(s9316_gcms, file=here("DAFs", "FC_data", "s9316_gcms.csv"))

s0013_gcms = lmFit(log_gcms, d5) |> eBayes() |> topTable(coef = 5, sort.by = "P", number = Inf)
write.csv(s0013_gcms, file=here("DAFs", "FC_data", "s0013_gcms.csv"))

