knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE, 
  message = FALSE,
  eval = FALSE
)

set.seed(123)
library(here)
library(ggplot2)
library(vegan)
library(stringr)
library(MetaboAnalystR)
library(rstatix)
library(gplots)
library(extrafont)
library(patchwork)


##########################
### data preprocessing ###
##########################

### loading the data
neg_quant <- read.csv(here("quantitative_data","ESI_neg_raw_data.csv"))[,-1]
pos_quant <- read.csv(here("quantitative_data","ESI_pos_raw_data.csv"))[,-1]
gc_data <- read.csv(here("quantitative_data","GC_raw_data.csv"))[,-1]

### transform into a matrix
neg_mat <- as.matrix(neg_quant)
pos_mat <- as.matrix(pos_quant)
gc_mat <- as.matrix(gc_data)

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

basenames <- unique(str_sub(colnames(gc_mat), end=-4))

keep_pos <- remove_low_abundance(pos_mat, basenames)
keep_neg <- remove_low_abundance(neg_mat, basenames)
keep_gc <- remove_low_abundance(gc_mat, basenames)

neg1 <- neg_mat[keep_neg, ]
pos1 <- pos_mat[keep_pos, ]

gc1 <- gc_mat[keep_gc, ]

### 2 filter out features with a large within-group relative standard deviation (=coefficient ov variation, CV) RSD=SD/mean metaboanalyst uses a 20% cut-off 
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

keep_pos2 <- filter2(pos1, basenames)
keep_neg2 <- filter2(neg1, basenames)
keep_gc2 <- filter2(gc1, basenames)

pos2 <- pos1[keep_pos2, ]
neg2 <- neg1[keep_neg2, ]
gc2 <- gc1[keep_gc2, ]

### check the IQR and remove the lowest 20% (variables that are near constant trghout the conditions)

filter3 <- function(qmatrix){
  iqr_range <- apply(qmatrix, 1, FUN=IQR)
  plot(rank(iqr_range,ties.method = "first"),iqr_range)
  to_keep_pos_iqr <- which(rank(iqr_range,ties.method = "first") > length(iqr_range)*0.2)
  return(to_keep_pos_iqr)
}

keep_pos3 <- filter3(pos2)
keep_neg3 <- filter3(neg2)
keep_gc3 <- filter3(gc2)

neg3 <- neg2[keep_neg3, ]
pos3 <- pos2[keep_pos3, ]
gc3 <- gc2[keep_gc3, ]

### normalize all values
neg_mat_norm <- log2(neg3+1) 
pos_mat_norm <- log2(pos3+1) 
gc_mat_norm <- log2(gc3+1) 

write.csv(neg_mat_norm, here("result_files","neg_mat_norm.csv"), row.names=FALSE)
write.csv(pos_mat_norm, here("result_files", "pos_mat_norm.csv"), row.names=FALSE)
write.csv(gc_mat_norm, here("result_files", "gc_mat_norm.csv"), row.names=FALSE)



neg_mat_norm <- read.csv(here("result_files","neg_mat_norm.csv"), header=T, sep = ',')

### extract metadata
## neg
ploidy <- str_sub(colnames(neg_mat_norm), 1,3)
strain <- str_sub(colnames(neg_mat_norm), start=5, end=-3)
neg_md <- data.frame(ploidy,strain)

write.csv(neg_md, here("result_files", "neg_md.csv"), row.names=FALSE )

##pos
ploidy <- str_sub(colnames(pos_mat_norm), 1,3)
strain <- str_sub(colnames(pos_mat_norm), start=5, end=-4)
pos_md <- data.frame(ploidy,strain)
write.csv(pos_md, here("result_files", "pos_md.csv"), row.names=FALSE)

##gc-ms
ploidy <- str_sub(colnames(gc_mat_norm), 1,3)
strain <- str_sub(colnames(gc_mat_norm), start=5, end=-4)
gc_md <- data.frame(ploidy,strain)
write.csv(gc_md, here("result_files", "gc_md.csv"), row.names=FALSE)

#===============================================================================
################
### var part ###
################

vp_neg <- varpart(Y=t(neg_mat_norm), ~strain, ~ploidy, data=neg_md)
vp_pos <- varpart(Y=t(pos_mat_norm), ~strain, ~ploidy, data=pos_md)
vp_gc <- varpart(Y=t(gc_mat_norm), ~strain, ~ploidy, data=gc_md)
vpa_neg <- plot(vp_neg, cutoff = -Inf, cex = 0.7, bg = c("hotpink", "skyblue"),
     Xnames=c("strain","ploidy"))
vpa_pos <- plot(vp_pos,cutoff = -Inf, cex = 0.7, bg = c("hotpink", "skyblue"),
     Xnames=c("strain","ploidy"))
vpa_gcms <- plot(vp_gc,cutoff = -Inf, cex = 0.7, bg = c("hotpink", "skyblue"),
     Xnames=c("strain","ploidy"))


#===============================================================================
#################
### permanova ###
#################
pe_neg <- adonis2(t(neg_mat_norm)~strain*ploidy, data=neg_md, method="bray", permutations = 999) #i tried wit 10000 permutations, only the p values change
pe_pos <- adonis2(t(pos_mat_norm)~strain*ploidy, data=pos_md, method="bray", permutations = 999)
pe_gcms <- adonis2(t(gc_mat_norm)~strain*ploidy, data=gc_md, method="bray", permutations = 999)

#===============================================================================
#################################
### two way anova per feature ###
#################################

neg_df <- data.frame(neg_mat_norm)
Sample <- rownames(neg_df)
neg_df2 <- cbind(Sample, neg_df)
Sample <- colnames(neg_df)
neg_md <- cbind(Sample, neg_md)
write.csv(neg_df2, here("result_files", "neg_df.csv"), row.names=FALSE)
write.csv(neg_md,here("result_files", "neg_md.csv"), row.names=FALSE)

pos_df <- data.frame(pos_mat_norm)
Sample <- rownames(pos_df)
pos_df2 <- cbind(Sample, pos_df)
Sample <- colnames(pos_df)
pos_md <- cbind(Sample, pos_md)
write.csv(pos_df2, here("result_files", "pos_df.csv"), row.names=FALSE)
write.csv(pos_md, here("result_files", "pos_md.csv"), row.names=FALSE)

gc_df <- data.frame(gc_mat_norm)
Sample <- rownames(gc_df)
gc_df2 <- cbind(Sample, gc_df)
Sample <- colnames(gc_df)
gc_md <- cbind(Sample, gc_md)
write.csv(gc_df2,here("result_files", "gc_df.csv"), row.names=FALSE)
write.csv(gc_md, here("result_files", "gc_md.csv"), row.names=FALSE)

### create a mSetObj object for statistical data analysis
mSet_neg <-InitDataObjects("conc", "mf", FALSE)
mSet_neg<-SetDesignType(mSet_neg, "multi")
mSet_neg<-Read.TextDataTs(mSet_neg, here("result_files", "neg_df.csv"), "colmf")
mSet_neg <- ReadMetaData(mSetObj = mSet_neg, here("result_files", "neg_md.csv"))
mSet_neg<-SanityCheckData(mSet_neg)
mSet_neg<-ReplaceMin(mSet_neg)
mSet_neg<-SanityCheckMeta(mSet_neg, 1)
mSet_neg<-SetDataTypeOfMeta(mSet_neg)
mSet_neg<-SanityCheckData(mSet_neg)
#mSet_neg<-FilterVariable(mSet_neg, "none", -1, "F", "NULL", F)
mSet_neg<-PreparePrenormData(mSet_neg)
#mSet_neg<-Normalization(mSet_neg, "NULL", "NULL", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet_neg<-Normalization(mSet_neg, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

mSet_pos <-InitDataObjects("conc", "mf", FALSE)
mSet_pos<-SetDesignType(mSet_pos, "multi")
mSet_pos<-Read.TextDataTs(mSet_pos, here("result_files", "pos_df.csv"), "colmf")
mSet_pos <- ReadMetaData(mSetObj = mSet_pos, here("result_files", "pos_md.csv"))
mSet_pos<-SanityCheckData(mSet_pos)
mSet_pos<-ReplaceMin(mSet_pos)
mSet_pos<-SanityCheckMeta(mSet_pos, 1)
mSet_pos<-SetDataTypeOfMeta(mSet_pos)
mSet_pos<-SanityCheckData(mSet_pos)
#mSet_pos<-FilterVariable(mSet_pos, "none", -1, "F", "NULL", F)
mSet_pos<-PreparePrenormData(mSet_pos)
#mSet_pos<-Normalization(mSet_pos, "NULL", "NULL", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet_pos<-Normalization(mSet_pos, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

mSet_gc <-InitDataObjects("conc", "mf", FALSE)
mSet_gc<-SetDesignType(mSet_gc, "multi")
mSet_gc<-Read.TextDataTs(mSet_gc, here("result_files", "gc_df.csv"), "colmf")
mSet_gc <- ReadMetaData(mSetObj = mSet_gc, here("result_files", "gc_md.csv"))
mSet_gc<-SanityCheckData(mSet_gc)
mSet_gc<-ReplaceMin(mSet_gc)
mSet_gc<-SanityCheckMeta(mSet_gc, 1)
mSet_gc<-SetDataTypeOfMeta(mSet_gc)
mSet_gc<-SanityCheckData(mSet_gc)
#mSet_gc<-FilterVariable(mSet_gc, "none", -1, "F", "NULL", F)
mSet_gc<-PreparePrenormData(mSet_gc)
#mSet_gc<-Normalization(mSet_gc, "NULL", "NULL", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet_gc<-Normalization(mSet_gc, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

#### from source version4 https://github.com/xia-lab/MetaboAnalystR/blob/master/R/multifac_pca_anova2.R
## ANOVA based on anova_test rstatix
aov.2way <- function(x){
  res <- suppressMessages(anova_test(x ~ aov.facA * aov.facB, 
                                     data = data.frame(x = x, aov.facA = aov.facA, aov.facB = aov.facB)))
  res <- c(res$F, res$p)
  res <- unlist(res)
  return(res)
}

p.cor="fdr"
thresh=0.05
mSet <- mSet_gc
#remove(mSet)
aov.facB <<-mSet$dataSet$meta.info$strain
aov.facA <<-mSet$dataSet$meta.info$ploidy
aov.mat <- t(apply(as.matrix(mSet$dataSet$norm), 2, aov.2way))
#aov.mat <- t(apply(t(mSet$dataSet$norm), 2, aov.2way))
aov.mat2 <- cbind (aov.mat, p.adjust(aov.mat[,4], p.cor),
                   p.adjust(aov.mat[,5], p.cor),
                   p.adjust(aov.mat[,6], p.cor))
sig.facA <-(aov.mat2[,7] <= thresh)
sig.facB <-(aov.mat2[,8] <= thresh)
sig.intr <-(aov.mat2[,9] <= thresh)

# make table for display/download  
all.match <- cbind(sig.facA, sig.facB, sig.intr)
colnames(all.match) <- c("Ploidy", "Strain", "Interaction")
colnames(aov.mat2) <- c(paste("Ploidy", "(F.val)", sep = ""), 
                        paste("Strain", "(F.val)", sep = ""),
                        paste("Interaction", "(F.val)", sep = ""),
                        paste("Ploidy", "(raw.p)", sep = "") ,
                        paste("Strain", "(raw.p)", sep = ""),
                        paste("Interaction", "(raw.p)", sep = ""),
                        paste("Ploidy", "(adj.p)", sep = "") ,
                        paste("Strain", "(adj.p)", sep = ""),
                        paste("Interaction", "(adj.p)", sep = ""))

vennC_neg <-venn(all.match)
vennC_pos <-venn(all.match)
vennC_gc <-venn(all.match)
#===============================================================================
####################
### venn diagram ###
####################

library(eulerr)
fit_neg <- venn(c("strain" = 1953, "ploidy" = 199, "interaction" = 37,
                   "strain&ploidy" = 810, "strain&interaction" = 756, "ploidy&interaction" = 17,
                   "strain&ploidy&interaction" = 813))

fit_pos <- venn(c("strain" = 1300, "ploidy" = 163, "interaction" = 2,
                "strain&ploidy" = 481, "strain&interaction" = 47, "ploidy&interaction" = 2,
                "strain&ploidy&interaction" = 37))

fit_gc <- venn(c("strain" = 22, "ploidy" = 6, "interaction" = 0,
                   "strain&ploidy" = 12, "strain&interaction" = 1, "ploidy&interaction" = 1,
                   "strain&ploidy&interaction" = 2))
col <- c("orchid3","skyblue","goldenrod1")

venn_neg <- plot(fit_neg, quantities = TRUE, fills = list(fill = col, alpha = 0.7), labels = list (col = "black", font =4), edges = FALSE, main = list(label = "LC-MS ESI-", cex = 1, family = "Helvetica"))
venn_pos <- plot(fit_pos, quantities = TRUE, fills = list(fill = col, alpha = 0.7), labels = list (col = "black", font =4), edges = FALSE, main = list(label = "LC-MS ESI+", cex = 1, family = "Helvetica"))
venn_gcms <- plot(fit_gc, quantities = TRUE, fills = list(fill = col, alpha = 0.7),  labels = list (col = "black", font =4), edges = FALSE, main = list(label = "GC-MS", cex = 1, family = "Helvetica"))

samplevenn_plots <- wrap_plots(
  wrap_plots(venn_pos, venn_neg, venn_gcms) +
    plot_layout(guides = "collect")
  
)
samplevenn_plots
