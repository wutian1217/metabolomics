.libPaths(c(.libPaths(), "D:/Program Files/R/R-4.2.3/library"))
set.seed(123)
library(here)
library(ggplot2)
library(stringr)
library(limma)
library(vegan)
library(extrafont)
library(patchwork)

##########################
### data preprocessing ###
##########################

### loading the data
#neg_quant <- read.csv(here("quantitative_data","ESI_neg_raw_data.csv"))
#pos_quant <- read.csv(here("quantitative_data","ESI_pos_raw_data.csv"))
#gc_data <- read.csv(here("quantitative_data","GC_raw_data.csv"))

neg_quant <- read.csv(here("quantitative_data","neg_select.csv"))
pos_quant <- read.csv(here("quantitative_data","pos_select.csv"))
gc_data <- read.csv(here("quantitative_data","GC_select.csv"))

### per strain
nq_0013 <- neg_quant[,grep("13", colnames(neg_quant))]
nq_9316 <- neg_quant[,grep("9316", colnames(neg_quant))]
nq_9242 <- neg_quant[,grep("9242", colnames(neg_quant))]
nq_9346 <- neg_quant[,grep("9346", colnames(neg_quant))]
pq_0013 <- pos_quant[,grep("13", colnames(pos_quant))]
pq_9316 <- pos_quant[,grep("9316", colnames(pos_quant))]
pq_9242 <- pos_quant[,grep("9242", colnames(pos_quant))]
pq_9346 <- pos_quant[,grep("9346", colnames(pos_quant))]
gq_0013 <- gc_data[,grep("13", colnames(gc_data))]
gq_9316 <- gc_data[,grep("9316", colnames(gc_data))]
gq_9242 <- gc_data[,grep("9242", colnames(gc_data))]
gq_9346 <- gc_data[,grep("9346", colnames(gc_data))]

### remove features with zeros
nq_0013 <- nq_0013[apply(nq_0013!=0, 1, all),]
nq_9316 <- nq_9316[apply(nq_9316!=0, 1, all),]
nq_9242 <- nq_9242[apply(nq_9242!=0, 1, all),]
nq_9346 <- nq_9346[apply(nq_9346!=0, 1, all),]
pq_0013 <- pq_0013[apply(pq_0013!=0, 1, all),]
pq_9316 <- pq_9316[apply(pq_9316!=0, 1, all),]
pq_9242 <- pq_9242[apply(pq_9242!=0, 1, all),]
pq_9346 <- pq_9346[apply(pq_9346!=0, 1, all),]
gq_0013 <- gq_0013[apply(gq_0013!=0, 1, all),]
gq_9316 <- gq_9316[apply(gq_9316!=0, 1, all),]
gq_9242 <- gq_9242[apply(gq_9242!=0, 1, all),]
gq_9346 <- gq_9346[apply(gq_9346!=0, 1, all),]

### log2 transform
nq_0013 <- log2(nq_0013)
nq_9316 <- log2(nq_9316)
nq_9242 <- log2(nq_9242)
nq_9346 <- log2(nq_9346)
pq_0013 <- log2(pq_0013)
pq_9316 <- log2(pq_9316)
pq_9242 <- log2(pq_9242)
pq_9346 <- log2(pq_9346)
gq_0013 <- log2(gq_0013)
gq_9316 <- log2(gq_9316)
gq_9242 <- log2(gq_9242)
gq_9346 <- log2(gq_9346)

#===============================================================================
#######################
### calculate logFC ###
#######################
setwd("F:/Spirodela_project_ERC/Pictures/Metabolomics/new_analysis/workdir")
dn=matrix(c(rep(1,16),rep(0,8),rep(1,8)), ncol=2, byrow=FALSE)
dp=matrix(c(rep(1,16),rep(1,8),rep(0,8)), ncol=2, byrow=FALSE)

s0013_negative = lmFit(as.matrix(nq_0013), dn) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s0013_negative, file=here("cell_norm","s0013_ESI-.csv"))
#s0013_negative_i = lmFit(as.matrix(nq_0013), dp) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
#write.csv(s0013_negative_i, file="s0013_negative_inv.csv")
s9316_negative = lmFit(as.matrix(nq_9316), dp) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9316_negative, file=here("cell_norm", "s9316_ESI-.csv"))
s9242_negative = lmFit(as.matrix(nq_9242), dp) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9242_negative, file=here("cell_norm", "s9242_ESI-.csv"))
s9346_negative = lmFit(as.matrix(nq_9346), dp) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9346_negative, file=here("cell_norm", "s9346_ESI-.csv"))

s0013_positive = lmFit(as.matrix(pq_0013), dn) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s0013_positive, file=here("cell_norm", "s0013_ESI+.csv"))
s9316_positive = lmFit(as.matrix(pq_9316), dn) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9316_positive, file=here("cell_norm", "s9316_ESI+.csv"))
s9242_positive = lmFit(as.matrix(pq_9242), dp) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9242_positive, file=here("cell_norm", "s9242_ESI+.csv"))
s9346_positive = lmFit(as.matrix(pq_9346), dp) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9346_positive, file=here("cell_norm", "s9346_ESI+.csv"))

s0013_gas = lmFit(as.matrix(gq_0013), dn) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s0013_gas, file=here("cell_norm", "s0013_GC-MS.csv"))
s9316_gas = lmFit(as.matrix(gq_9316), dn) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9316_gas, file=here("cell_norm", "s9316_GC-MS.csv"))
s9242_gas = lmFit(as.matrix(gq_9242), dn) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9242_gas, file=here("cell_norm", "s9242_GC-MS.csv"))
s9346_gas = lmFit(as.matrix(gq_9346), dn) |> eBayes() |> topTable(coef = 2, sort.by = "P", number = Inf)
write.csv(s9346_gas, file=here("cell_norm", "s9346_GC-MS.csv"))

#===============================================================================
#############################
### correct and visualise ###
#############################

sv <- c("s0013","s9242","s9316","s9346")
tfv <- c(0.635,0.678,0.673,0.646)
colors <- c("#20854E","#E18727","#7876B1","#EE4C97")
df_vect <- c("s0013_ESI+","s0013_ESI-","s0013_GC-MS",
             "s9346_ESI+","s9346_ESI-","s9346_GC-MS",
             "s9316_ESI+","s9316_ESI-", "s9316_GC-MS",
             "s9242_ESI+","s9242_ESI-", "s9242_GC-MS")

for(i in df_vect){
  strain <- substring(i, 1,5)
  tf <- tfv[which(sv==strain)]
  kleur <- colors[which(sv==strain)]
  inf <- here("cell_norm", paste0(i, ".csv"))
  outf <- here("cell_norm", paste0(i, ".png"))
  met_df <- read.csv(inf)
  met_df$FC <- 2**met_df$logFC
  met_df$FCpc <- met_df$FC/tf
  st = paste0("mean dm: ", as.character(round(mean(met_df$FC), digits=3)), ", median dm: ", as.character(round(median(met_df$FC), digits=3)), 
  "\nmean cell: ", as.character(round(mean(met_df$FCpc), digits=3)), ", median cell: ", as.character(round(median(met_df$FCpc), digits=3)))
  
  fig <- ggplot(met_df)+
    geom_histogram(aes(x=FC, y=..density..), alpha=.5, fill="grey") +
    geom_density(aes(x=FC), alpha=.5, fill="grey") +
    geom_histogram(aes(x=FCpc, y=..density..), alpha=.5, fill=kleur) +
    geom_density(aes(x=FCpc), alpha=.5, fill=kleur, color=kleur) +
    xlim(c(0, 5))+
    theme_classic(18) +
    labs(title=i, subtitle=st) +
    labs(x = "Fold change") +
    geom_vline(xintercept = median(met_df$FC), color="grey", linetype = "dashed", size = 1.5)+
    geom_vline(xintercept = median(met_df$FCpc), color=kleur, linetype = "dashed", size = 1.5)+
    theme(plot.title = element_text(hjust=0.5, size = 22), 
          plot.subtitle = element_text(size=17),
          axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          axis.line.x = element_line(size = 0.7),
          axis.line.y = element_line(size = 0.7),
          )
  ggsave(outf, fig)
}
