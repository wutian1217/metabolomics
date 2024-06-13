
set.seed(123)
library(here)
library(UpSetR)
library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggplotify)


pos_up_list <- read.csv(here("DAFs", "pos_up.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)
pos_down_list <- read.csv(here("DAFs", "pos_down.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)
pos_diff_list <- read.csv(here("DAFs", "pos_diff.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)

neg_up_list <- read.csv(here("DAFs", "neg_up.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)
neg_down_list <- read.csv(here("DAFs", "neg_down.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)
neg_diff_list <- read.csv(here("DAFs", "neg_diff.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)

gcms_up_list <- read.csv(here("DAFs", "gcms_up.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)
gcms_down_list <- read.csv(here("DAFs", "gcms_down.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)
gcms_diff_list <- read.csv(here("DAFs", "gcms_diff.list"), header = T, sep = ',', na.strings = "") %>% lapply(na.omit)



###plot
f_pos_up <- upset(fromList(pos_up_list), nsets=4,
      sets=c("s0013","s9242","s9316","s9346"),
      #queries = list(list(query = intersects, params = list('s0013', 's9242', 's9316', 's9346'), color = "darkorange", active =T)),
      keep.order = TRUE,
      #order.by = "freq",
      mb.ratio = c(0.6,0.4),
      #set_size.show = T,
      mainbar.y.label = "Intersection size",
      sets.x.label = "up-regulated DAFs, ESI+",
      text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
      set_size.angles = 0,
      #set_size.numbers_size = 8,
      #scale.sets = "identity",
      sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_pos_up

f_pos_down <- upset(fromList(pos_down_list), nsets=length(pos_down_list),
                  sets=c("s0013","s9242","s9316","s9346"),
                  keep.order = TRUE,
                  #order.by = "freq",
                  mb.ratio = c(0.6,0.4),
                  #set_size.show = T,
                  mainbar.y.label = "Intersection size",
                  sets.x.label = "down-regulated DAFs, ESI+",
                  text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
                  set_size.angles = 0,
                  sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_pos_down

f_pos_diff <- upset(fromList(pos_diff_list), nsets=length(pos_diff_list),
                    sets=c("s0013","s9242","s9316","s9346"),
                    keep.order = TRUE,
                    #order.by = "freq",
                    mb.ratio = c(0.6,0.4),
                    #set_size.show = T,
                    mainbar.y.label = "Intersection size",
                    sets.x.label = "total DAFs, ESI+",
                    text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
                    set_size.angles = 0,
                    sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_pos_diff

f_neg_up <- upset(fromList(neg_up_list), nsets=length(neg_up_list),
                  sets=c("s0013","s9242","s9316","s9346"),
                  queries = list(list(query = intersects, params = list('s0013', 's9242', 's9316', 's9346'), color = "darkorange", active =T)),
                  keep.order = TRUE,
                  #order.by = "freq",
                  mb.ratio = c(0.6,0.4),
                  #set_size.show = T,
                  mainbar.y.label = "Intersection size",
                  sets.x.label = "up-regulated DAFs, ESI-",
                  text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
                  set_size.angles = 0,
                  sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_neg_up

f_neg_down <- upset(fromList(neg_down_list), nsets=length(neg_down_list),
                  sets=c("s0013","s9242","s9316","s9346"),
                  queries = list(list(query = intersects, params = list('s0013', 's9242', 's9316', 's9346'), color = "darkorange", active =T)),
                  keep.order = TRUE,
                  #order.by = "freq",
                  mb.ratio = c(0.6,0.4),
                  #set_size.show = T,
                  mainbar.y.label = "Intersection size",
                  sets.x.label = "down-regulated DAFs, ESI-",
                  text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
                  set_size.angles = 0,
                  #set_size.numbers_size = 8,
                  #scale.sets = "identity",
                  sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_neg_down

f_neg_diff <- upset(fromList(neg_diff_list), nsets=length(neg_diff_list),
                    sets=c("s0013","s9242","s9316","s9346"),
                    queries = list(list(query = intersects, params = list('s0013', 's9242', 's9316', 's9346'), color = "darkorange", active =T)),
                    keep.order = TRUE,
                    #order.by = "freq",
                    mb.ratio = c(0.6,0.4),
                    #set_size.show = T,
                    mainbar.y.label = "Intersection size",
                    sets.x.label = "total DAFs, ESI-",
                    text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
                    set_size.angles = 0,
                    #set_size.numbers_size = 8,
                    #scale.sets = "identity",
                    sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_neg_diff

f_gcms_up <- upset(fromList(gcms_up_list), nsets=length(gcms_up_list),
                  sets=c("s0013","s9242","s9316","s9346"),
                  keep.order = TRUE,
                  #order.by = "freq",
                  mb.ratio = c(0.6,0.4),
                  #set_size.show = T,
                  mainbar.y.label = "Intersection size",
                  sets.x.label = "up-regulated DAFs, GC-MS",
                  text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
                  set_size.angles = 0,
                  #set_size.numbers_size = 8,
                  #scale.sets = "identity",
                  sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_gcms_up

f_gcms_down <- upset(fromList(gcms_down_list), nsets=length(gcms_down_list),
                   sets=c("s0013","s9242","s9316","s9346"),
                   keep.order = TRUE,
                   #order.by = "freq",
                   mb.ratio = c(0.6,0.4),
                   #set_size.show = T,
                   mainbar.y.label = "Intersection size",
                   sets.x.label = "down-regulated DAFs, GC-MS",
                   text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
                   set_size.angles = 0,
                   #set_size.numbers_size = 8,
                   #scale.sets = "identity",
                   sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_gcms_down

f_gcms_diff <- upset(fromList(gcms_diff_list), nsets=length(gcms_diff_list),
                     sets=c("s0013","s9242","s9316","s9346"),
                     keep.order = TRUE,
                     #order.by = "freq",
                     mb.ratio = c(0.6,0.4),
                     #set_size.show = T,
                     mainbar.y.label = "Intersection size",
                     sets.x.label = "total DAFs, GC-MS",
                     text.scale = c(1, 1, 1, 1, 1.3, 1.3), 
                     set_size.angles = 0,
                     #set_size.numbers_size = 8,
                     #scale.sets = "identity",
                     sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97")
)
f_gcms_diff


upset_plots <- wrap_plots(
  wrap_plots(
    as.ggplot(f_pos_up), as.ggplot(f_neg_up),as.ggplot(f_gcms_up),
    nrow = 1
  ),
  wrap_plots(
    as.ggplot(f_pos_down), as.ggplot(f_neg_down),as.ggplot(f_gcms_down), 
    nrow = 1
  ),
  wrap_plots(
    as.ggplot(f_pos_diff), as.ggplot(f_neg_diff),as.ggplot(f_gcms_diff), 
    nrow = 1
  ),
  nrow = 3
)



upset_plots
