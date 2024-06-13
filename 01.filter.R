knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE, 
  message = FALSE,
  eval = FALSE
)

set.seed(123)
library(here)
library(SummarizedExperiment)
library(tidyverse)



# Load quantitative data
neg_quant <- read_csv(
  here("quantitative_data", "neg_all_compounds_normalized_dw.csv"),
  skip = 2, show_col_types = FALSE
) |>
  janitor::clean_names()

pos_quant <- read_csv(
  here("quantitative_data", "pos_all_compounds_normalized_dw.csv"),
  skip = 2, show_col_types = FALSE
) |>
  janitor::clean_names()

gcms_quant <- read_csv(
  here("quantitative_data", "GCMS_all_compounds_normalized_dw.csv"),
  #skip =1, 
  show_col_types = FALSE
) |>
  janitor::clean_names()


# Load sample metadata
neg_metadata <- read.csv(
  here("quantitative_data", "neg_metadata.csv")
)

pos_metadata <- read.csv(
  here("quantitative_data", "pos_metadata.csv")
)

gcms_metadata <- read.csv(
  here("quantitative_data", "GCMS_metadata.csv")
)

# assay (quantitative matrices)
## Negative
neg_assay_normalized <- as.matrix(neg_quant[, 16:84])
colnames(neg_assay_normalized) <- str_replace_all(
  colnames(neg_assay_normalized), c(
    "_121121.*" = "",
    "x" = "s"
  )
)
colnames(neg_assay_normalized)[17] <- "pooled_sample05"


## Positive
pos_assay_normalized <- as.matrix(pos_quant[, 16:84])
colnames(pos_assay_normalized) <- str_replace_all(
  colnames(pos_assay_normalized), c(
    "_161121.*" = "",
    "x" = "s"
  )
)


###GC_MS
gcms_assay_normalized <- as.matrix(gcms_quant[, 2:70])
colnames(gcms_assay_normalized) <- str_replace_all(
  colnames(gcms_assay_normalized), c(
    #"_161121.*" = "",
    "x" = "s",
    "mp0242_quinten" = "pooled"
  )
)


# colData (sample metadata)
neg_coldata <- data.frame(
  Sample = colnames(neg_assay_normalized)
) |>
  mutate(
    Strain = as.factor(str_replace_all(Sample, "_.*", "")),
    Ploidy = as.factor(case_when(
      str_detect(Sample, "tet") ~ "Tetraploid",
      str_detect(Sample, "dip") ~ "Diploid",
      TRUE ~ "Mixed"
    ))
  ) |>
  tibble::column_to_rownames("Sample")


pos_coldata <- data.frame(
  Sample = colnames(pos_assay_normalized)
) |>
  mutate(
    Strain = as.factor(str_replace_all(Sample, "_.*", "")),
    Ploidy = as.factor(case_when(
      str_detect(Sample, "tet") ~ "Tetraploid",
      str_detect(Sample, "dip") ~ "Diploid",
      TRUE ~ "Mixed"
    ))
  ) |>
  tibble::column_to_rownames("Sample")


gcms_coldata <- data.frame(
  Sample = colnames(gcms_assay_normalized)
) |>
  mutate(
    Strain = as.factor(str_replace_all(Sample, "_.*", "")),
    Ploidy = as.factor(case_when(
      str_detect(Sample, "tet") ~ "Tetraploid",
      str_detect(Sample, "dip") ~ "Diploid",
      TRUE ~ "Mixed"
    ))
  ) |>
  tibble::column_to_rownames("Sample")



# rowData (compound metadata)
neg_rowdata <- data.frame(
  row.names = neg_quant$compound,
  m_z = neg_quant$m_z,
  rt = neg_quant$retention_time_min,
  chromatographic_peak_width = neg_quant$chromatographic_peak_width_min
)

pos_rowdata <- data.frame(
  row.names = pos_quant$compound,
  m_z = pos_quant$m_z,
  rt = pos_quant$retention_time_min,
  chromatographic_peak_width = pos_quant$chromatographic_peak_width_min
)


gcms_rowdata <- data.frame(
  row.names = gcms_quant$sample,
  compound = gcms_quant$sample
)


# Create SummarizedExperiment objects
## Negative
se_neg <- SummarizedExperiment(
  assays = list(
    normalized = neg_assay_normalized
  ),
  rowData = neg_rowdata,
  colData = neg_coldata
)

## Positive
se_pos <- SummarizedExperiment(
  assays = list(
    normalized = pos_assay_normalized
  ),
  rowData = pos_rowdata,
  colData = pos_coldata
)

###GCMS
gcms <- SummarizedExperiment(
  assays = list(
    normalized = gcms_assay_normalized
  ),
  colData = gcms_coldata,
  rowData = gcms_rowdata
)


##########filter####################################################

# get a vector of the sample_names
rel_names<-colnames(assay(se_pos))[-which(str_detect(colnames(assay(se_pos)), "pool"))]
basenames <- unique(trimws(rel_names, whitespace = "[_r1/_r2/_r3/_r4/_r5/_r6/_r7/_r8]"))


### filter out lowly abundant features, given that we have 8 replicates the feature has to be present in at least 8 samples of the ploidy-strain combination to be retained
max_vect<-function(v1,v2){
  for (i in 1:length(v1)){
    if (v2[i] > v1[i]){
      v1[i] <- v2[i]
    }
  }
  return(v1)
}

remove_low_abundance <- function(qmatrix, basenames){
  non_zero <- rep(0,nrow(qmatrix))
  #summing al features per strain_ploidy combination and keeping the maximum
  for (name in basenames){
    name_df <-qmatrix[,which(str_detect(colnames(qmatrix), name))]
    nz_vect <-rowSums(name_df>0)
    non_zero <-max_vect(non_zero,nz_vect)
  }
  keep <- which(non_zero == 8)
  return(keep)
}

keep_pos <- remove_low_abundance(assay(se_pos, "normalized"), basenames)
keep_neg <- remove_low_abundance(assay(se_neg,"normalized"), basenames)
keep_gcms <- remove_low_abundance(assay(gcms,"normalized"), basenames)

se_neg1 <- se_neg[keep_neg, ]
se_pos1 <- se_pos[keep_pos, ]
gcms1 <- gcms[keep_gcms, ]

### 2 filter out features with a large within-group relative standard deviation (=coefficient of variation, CV) RSD=SD/mean (20% cut-off) 
RSD <- function(x){sd(x) / mean(x)}

filter2 <- function(qmatrix, basenames){
  max_rsd <- rep(0,nrow(qmatrix))
  for (name in basenames){
    name_df <- qmatrix[,which(str_detect(colnames(qmatrix), name))]
    print(name_df)
    rsd_vect <- apply(name_df,1,FUN=RSD)
    rsd_vect <- rsd_vect %>% replace(is.na(.), 0)
    max_rsd <-max_vect(max_rsd,rsd_vect)
  }
  plot(rank(max_rsd,ties.method = "first"),max_rsd )
  to_keep_pos_rsd <- which(rank(max_rsd,ties.method = "first") < length(rsd_vect)*0.8)
  #to_keep_pos_rsd <- which(max_rsd < 1)
  return(to_keep_pos_rsd)
}

keep_pos2 <- filter2(assay(se_pos1, "normalized"), basenames)
keep_neg2 <- filter2(assay(se_neg1, "normalized"), basenames)
keep_gcms2 <- filter2(assay(gcms1, "normalized"), basenames)

se_pos2 <- se_pos1[keep_pos2, ]
se_neg2 <- se_neg1[keep_neg2, ]
gcms2 <- gcms1[keep_gcms2, ]

### check the IQR and remove the lowest 20% (variables that are near constant trghout the conditions)

filter3 <- function(qmatrix){
  iqr_range <- apply(qmatrix, 1, FUN=IQR)
  plot(rank(iqr_range,ties.method = "first"),iqr_range)
  to_keep_pos_iqr <- which(rank(iqr_range,ties.method = "first") > length(iqr_range)*0.2)
  return(to_keep_pos_iqr)
}

keep_pos3 <- filter3(assay(se_pos2, "normalized"))
keep_neg3 <- filter3(assay(se_neg2, "normalized"))
keep_gcms3 <- filter3(assay(gcms2, "normalized"))

se_pos <- se_pos2[keep_pos3, ]
se_neg <- se_neg2[keep_neg3, ]
gcms3 <- gcms2[keep_gcms3, ] 

save(
  se_neg, compress = "xz",
  file = here("result_files", "se_neg.rda")
)

save(
  se_pos, compress = "xz",
  file = here("result_files", "se_pos.rda")
)

save(
  gcms3, compress = "xz",
  file = here("result_files", "gcms3.rda")
)


#######reload data for normalization (can start from here after filtering)#################

load(here::here("result_files", "se_pos.rda"))
load(here::here("result_files", "se_neg.rda"))
load(here::here("result_files", "gcms3.rda"))


# Get sample metadata
metadata_positive <- as.data.frame(colData(se_pos)) |> filter(Strain != "pooled")
metadata_negative <- as.data.frame(colData(se_neg)) |> filter(Strain != "pooled")
metadata_gcms <- as.data.frame(colData(gcms3)) |> filter(Strain != "pooled")

# Log transform data
log_positive <- log2(assay(se_pos, "normalized") + 1)
log_positive <- log_positive[, !str_detect(colnames(log_positive), "pool")]

log_negative <- log2(assay(se_neg, "normalized") + 1)
log_negative <- log_negative[, !str_detect(colnames(log_negative), "pool")]
log_negative <- log_negative[, colnames(log_positive)]

log_gcms <- log2(assay(gcms3, "normalized") + 1)
log_gcms <- log_gcms[, !str_detect(colnames(log_gcms), "pool")]
log_gcms <- log_gcms[, colnames(log_positive)]

write.csv(x= log_negative, file = here("result_files", "log_negative.csv"))
write.csv(x = log_positive, file = here("result_files", "log_positive.csv"))
write.csv(x=log_gcms, file = here("result_files", "log_gcms.csv"))

write.csv(x=metadata_positive, file = here("result_files", "metadata_positive.csv"))
write.csv(x=metadata_negative, file = here("result_files", "metadata_negative.csv"))
write.csv(x=metadata_gcms, file = here("result_files", "metadata_gcms.csv"))

save(
  log_positive, file = here("result_files", "log_positive.rda")
)

save(
  log_negative, file = here("result_files", "log_negative.rda")
)

save(
  log_gcms, file = here("result_files", "log_gcms.rda")
)


