###PLS_DA analysis

set.seed(123)
library(here)
library(ropls)
library(ggplot2)
library(mixOmics)
library(ggforce)
library(ggplot2)
library(ggsci)
library(tidyverse)


aspect.ratio <- 4/3

#######positive data############

data_pos <- read.csv(here("result_files","log_positive.csv"), header=T, row.names = 1, sep = ',')
data1_pos <- t(data_pos)
group_pos <- read.csv(here("result_files","metadata_positive.csv"), row.names = 1, sep = ',', header=T)

######data_strain
s0013_data_pos <- data1_pos[grepl("s0013", rownames(data1_pos)),]
s9242_data_pos <- data1_pos[grepl("s9242", rownames(data1_pos)),]
s9316_data_pos <- data1_pos[grepl("s9316", rownames(data1_pos)),]
s9346_data_pos <- data1_pos[grepl("s9346", rownames(data1_pos)),]

group_s0013_pos <- group_pos[grepl("s0013", rownames(group_pos)),]
group_s9242_pos <- group_pos[grepl("s9242", rownames(group_pos)),]
group_s9316_pos <- group_pos[grepl("s9316", rownames(group_pos)),]
group_s9346_pos <- group_pos[grepl("s9346", rownames(group_pos)),]


####s0013
pos_plsda <- opls(x = s0013_data_pos, y = group_s0013_pos[, 'Ploidy'], orthoI = 0)
pos_plsda
getSummaryDF(pos_plsda)
data2<- as.data.frame(pos_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s0013_pos$Strain
data2$ploidy=group_s0013_pos$Ploidy
x_lab <- pos_plsda@modelDF[1, "R2X"] * 100
y_lab <- pos_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s0013_pos <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
    labs(title = "s0013, ESI+", color = '') +
    theme(aspect.ratio=4/3)
    
f_s0013_pos

data_vip <- pos_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s0013_pos_vip.csv"))


####s9242
pos_plsda <- opls(x = s9242_data_pos, y = group_s9242_pos[, 'Ploidy'], orthoI = 0)
pos_plsda
getSummaryDF(pos_plsda)
data2<- as.data.frame(pos_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9242_pos$Strain
data2$ploidy=group_s9242_pos$Ploidy
x_lab <- pos_plsda@modelDF[1, "R2X"] * 100
y_lab <- pos_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9242_pos <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
    labs(title = "s9242, ESI+", color = '') +
    theme(aspect.ratio=4/3)

f_s9242_pos

data_vip <- pos_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s9242_pos_vip.csv"))



####s9316
pos_plsda <- opls(x = s9316_data_pos, y = group_s9316_pos[, 'Ploidy'], predI =2, orthoI = 0)
pos_plsda
data2<- as.data.frame(pos_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9316_pos$Strain
data2$ploidy=group_s9316_pos$Ploidy
x_lab <- pos_plsda@modelDF[1, "R2X"] * 100
y_lab <- pos_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9316_pos <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
    labs(title = "s9316, ESI+", color = '') +
    theme(aspect.ratio=4/3)

f_s9316_pos

data_vip <- pos_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s9316_pos_vip.csv"))



####s9346
pos_plsda <- opls(x = s9346_data_pos, y = group_s9346_pos[, 'Ploidy'], orthoI = 0)
pos_plsda
data2<- as.data.frame(pos_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9346_pos$Strain
data2$ploidy=group_s9346_pos$Ploidy
x_lab <- pos_plsda@modelDF[1, "R2X"] * 100
y_lab <- pos_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9346_pos <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
        labs(title = "s9346, ESI+", color = '') +
        theme(aspect.ratio=4/3)

f_s9346_pos

data_vip <- pos_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s9346_pos_vip.csv"))


####pos all data
pos_plsda <- opls(x = data1_pos, y = group_pos[, 'Ploidy'], orthoI = 0)
pos_plsda
data2<- as.data.frame(pos_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_pos$Strain
data2$ploidy=group_pos$Ploidy
x_lab <- pos_plsda@modelDF[1, "R2X"] * 100
y_lab <- pos_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_pos <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
    labs(title = "ESI+", color = '') +
    theme(aspect.ratio=4/3)

f_pos

data_vip <- pos_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","pos_vip.csv"))


################neg data###########################

data <- read.csv(here("result_files","log_negative.csv"), header=T, row.names = 1, sep = ',')
group <- read.csv(here("result_files","metadata_negative.csv"), row.names = 1, sep = ',', header=T)
data1 <- t(data)

######data_strain
s0013_data <- data1[grepl("s0013", rownames(data1)),]
s9242_data <- data1[grepl("s9242", rownames(data1)),]
s9316_data <- data1[grepl("s9316", rownames(data1)),]
s9346_data <- data1[grepl("s9346", rownames(data1)),]

group_s0013 <- group[grepl("s0013", rownames(group)),]
group_s9242 <- group[grepl("s9242", rownames(group)),]
group_s9316 <- group[grepl("s9316", rownames(group)),]
group_s9346 <- group[grepl("s9346", rownames(group)),]


####s0013
neg_plsda <- opls(x = s0013_data, y = group_s0013[, 'Ploidy'], orthoI = 0)
neg_plsda

data2<- as.data.frame(neg_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s0013$Strain
data2$ploidy=group_s0013$Ploidy
x_lab <- neg_plsda@modelDF[1, "R2X"] * 100
y_lab <- neg_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s0013_neg <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
        labs(title = "s0013, ESI-", color = '') +
        theme(aspect.ratio=4/3)

f_s0013_neg

data_vip <- neg_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s0013_neg_vip.csv"))


####s9242
neg_plsda <- opls(x = s9242_data, y = group_s9242[, 'Ploidy'], orthoI = 0)
neg_plsda

data2<- as.data.frame(neg_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9242$Strain
data2$ploidy=group_s9242$Ploidy
x_lab <- neg_plsda@modelDF[1, "R2X"] * 100
y_lab <- neg_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9242_neg <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
        labs(title = "s9242, ESI-", color = '') +
        theme(aspect.ratio=4/3)

f_s9242_neg

data_vip <- neg_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s9242_neg_vip.csv"))

####s9316
neg_plsda <- opls(x = s9316_data, y = group_s9316[, 'Ploidy'], predI =2, orthoI = 0)
neg_plsda

data2<- as.data.frame(neg_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9316$Strain
data2$ploidy=group_s9316$Ploidy
x_lab <- neg_plsda@modelDF[1, "R2X"] * 100
y_lab <- neg_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9316_neg <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
        labs(title = "s9316, ESI-", color = '') +
        theme(aspect.ratio=4/3)

f_s9316_neg

data_vip <- neg_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s9316_neg_vip.csv"))


####s9346
neg_plsda <- opls(x = s9346_data, y = group_s9346[, 'Ploidy'], orthoI = 0)
neg_plsda

data2<- as.data.frame(neg_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9346$Strain
data2$ploidy=group_s9346$Ploidy
x_lab <- neg_plsda@modelDF[1, "R2X"] * 100
y_lab <- neg_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9346_neg <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank()) +
        labs(title = "s9346, ESI-", color = '') +
         theme(aspect.ratio=4/3)

f_s9346_neg

data_vip <- neg_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s9346_neg_vip.csv"))


####neg all data
neg_plsda <- opls(x = data1, y = group[, 'Ploidy'], orthoI = 0)
neg_plsda
data2<- as.data.frame(neg_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group$Strain
data2$ploidy=group$Ploidy
x_lab <- pos_plsda@modelDF[1, "R2X"] * 100
y_lab <- pos_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_neg <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
    labs(title = "negative, ESI-", color = '') +
    theme(aspect.ratio=4/3)

f_neg

data_vip <- neg_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","neg_vip.csv"))



##############GC-MS data###############

data_gcms <- read.csv(here("result_files","log_gcms.csv"), header=T, row.names = 1, sep = ',')
data1_gcms <- t(data_gcms)
group_gcms <- read.csv(here("result_files", "metadata_gcms.csv"), row.names = 1, sep = ',', header=T)

######data_strain
s0013_data <- data1_gcms[grepl("s0013", rownames(data1_gcms)),]
s9242_data <- data1_gcms[grepl("s9242", rownames(data1_gcms)),]
s9316_data <- data1_gcms[grepl("s9316", rownames(data1_gcms)),]
s9346_data <- data1_gcms[grepl("s9346", rownames(data1_gcms)),]

group_s0013 <- group_gcms[grepl("s0013", rownames(group_gcms)),]
group_s9242 <- group_gcms[grepl("s9242", rownames(group_gcms)),]
group_s9316 <- group_gcms[grepl("s9316", rownames(group_gcms)),]
group_s9346 <- group_gcms[grepl("s9346", rownames(group_gcms)),]


s0013_plsda <- opls(x = s0013_data, y = group_s0013[, 'Ploidy'], orthoI = 0)
s0013_plsda

data2<- as.data.frame(s0013_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s0013$Strain
data2$ploidy=group_s0013$Ploidy
x_lab <- s0013_plsda@modelDF[1, "R2X"] * 100
y_lab <- s0013_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s0013_gcms <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
    labs(title = "s0013, GCMS", color = '') +
    theme(aspect.ratio=4/3)

f_s0013_gcms

data_vip <- s0013_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s0013_gcms_vip.csv"))

####s9242
s9242_plsda <- opls(x = s9242_data, y = group_s9242[, 'Ploidy'], orthoI = 0)
s9242_plsda

data2<- as.data.frame(s9242_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9242$Strain
data2$ploidy=group_s9242$Ploidy
x_lab <- s9242_plsda@modelDF[1, "R2X"] * 100
y_lab <- s9242_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9242_gcms <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
    labs(title = "s9242, gcms", color = '') +
    theme(aspect.ratio=4/3)

f_s9242_gcms

data_vip <- s9242_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA", "s9242_gcms_vip.csv"))

####s9316
s9316_plsda <- opls(x = s9316_data, y = group_s9316[, 'Ploidy'], predI =2, orthoI = 0)
s9316_plsda

data2<- as.data.frame(s9316_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9316$Strain
data2$ploidy=group_s9316$Ploidy
x_lab <- s9316_plsda@modelDF[1, "R2X"] * 100
y_lab <- s9316_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9316_gcms <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank(),
          legend.position="none") +
    labs(title = "s9316, gcms", color = '') +
    theme(aspect.ratio=4/3)

f_s9316_gcms

data_vip <- s9316_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA", "s9316_gcms_vip.csv"))


####s9346
s9346_plsda <- opls(x = s9346_data, y = group_s9346[, 'Ploidy'], orthoI = 0)
s9346_plsda

data2<- as.data.frame(s9346_plsda@scoreMN)
data2$sample = rownames(data2)
data2$strain=group_s9346$Strain
data2$ploidy=group_s9346$Ploidy
x_lab <- s9346_plsda@modelDF[1, "R2X"] * 100
y_lab <- s9346_plsda@modelDF[2, "R2X"] * 100
col=c("#0073C2", "#EFC000")

f_s9346_gcms <- ggplot(data2,aes(x=p1,y=p2, color=ploidy)) +
    geom_point(size=1.8)+
    theme_bw() +
    coord_fixed() +
    theme(panel.grid = element_blank()) +
    guides(color=guide_legend(title="ploidy")) +
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("P2 (",y_lab,"%)")) +
    stat_ellipse(data=data2,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=ploidy),
                 alpha=0.2,
                 show.legend = F)+
    scale_color_manual(values = col) +
    scale_fill_manual(values = c("#0073C2", "#EFC000"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank()) +
    labs(title = "s9346, gcms", color = '') +
    theme(aspect.ratio=4/3)

f_s9346_gcms

data_vip <- s9346_plsda@vipVn
data_vip
VIP <- data_vip[data_vip >1]
write.csv(x = VIP, file = here("PLS_DA","s9346_gcms_vip.csv"))



####plots

plsda_plots <- wrap_plots(f_s0013_pos, f_s9242_pos, f_s9316_pos, f_s9346_pos, f_s0013_neg, f_s9242_neg, f_s9316_neg, f_s9346_neg,
                          f_s0013_gcms,  f_s9316_gcms, f_s9346_gcms,
                        ncol = 4, nrow = 3) +
  plot_layout(guides = "collect") 

plsda_plots






