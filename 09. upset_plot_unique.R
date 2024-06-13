set.seed(123)
library(here)
library(UpSetR)
library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggplotify)


########add unique Features#######
##neg_up
neg_uniq_up <- read.csv(here("DAFs", "unique_DAFs", "neg_up_adduniq.list"), header = T, na.strings = "") %>% lapply(na.omit) 
neg_uniq_up <- fromList(neg_uniq_up)

neg_uniq_up$unique[neg_uniq_up$unique == '1'] <- 'Y' 
neg_uniq_up$unique[neg_uniq_up$unique == '0'] <- 'N'

f_neg_up <- upset(neg_uniq_up, 
                  sets=c("s0013","s9242","s9316","s9346"),
                  query.legend = "top",
                  queries = list(list(query = elements, params = list("unique", "Y"), color = "#FF6633", active =T, query.name = "unique DAFs")),
                                 #list(query = intersects, params = list('s0013', 's9242', 's9316', 's9346'), active =T, query.name = "shared DAFs")),
                  keep.order = TRUE,
                  mb.ratio = c(0.65,0.35),
                  mainbar.y.label = "number of DAFs",
                  sets.x.label = "up-regulated DAFs, ESI-",
                  text.scale = c(2, 1.3, 1.3, 1.5, 2, 2), 
                  point.size = 3,
                  set_size.show = TRUE,
                  set_size.angles = 0,
                  show.numbers="yes",
                  sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97"),
                  #attribute.plots=list(gridrows=60,plots=list(list(plot=histogram, x="unique")))
)
f_neg_up

###neg_down
neg_uniq_down <- read.csv(here("DAFs", "unique_DAFs", "neg_down_adduniq.list"), header = T, na.strings = "") %>% lapply(na.omit) 
neg_uniq_down <- fromList(neg_uniq_down)

neg_uniq_down$unique[neg_uniq_down$unique == '1'] <- 'Y' 
neg_uniq_down$unique[neg_uniq_down$unique == '0'] <- 'N'

f_neg_down <- upset(neg_uniq_down, 
                  sets=c("s0013","s9242","s9316","s9346"),
                  query.legend = "top",
                  queries = list(list(query = elements, params = list("unique", "Y"), color = "#FF6633", active =T, query.name = "unique DAFs")),
                                 #list(query = intersects, params = list('s0013', 's9242', 's9316', 's9346'), active =T, query.name = "shared DAFs")),
                  keep.order = TRUE,
                  mb.ratio = c(0.65,0.35),
                  mainbar.y.label = "number of DAFs",
                  sets.x.label = "down-regulated DAFs, ESI-",
                  text.scale = c(2, 1.3, 1.3, 1.5, 2, 2), 
                  point.size = 3,
                  set_size.show = TRUE,
                  set_size.angles = 0,
                  sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97"),
                  #attribute.plots=list(gridrows=60,plots=list(list(plot=histogram, x="unique")))
)
f_neg_down

##neg_diff
neg_uniq_diff <- read.csv(here("DAFs", "unique_DAFs", "neg_diff_adduniq.list"), header = T, na.strings = "") %>% lapply(na.omit) 
neg_uniq_diff <- fromList(neg_uniq_diff)

neg_uniq_diff$unique[neg_uniq_diff$unique == '1'] <- 'Y' 
neg_uniq_diff$unique[neg_uniq_diff$unique == '0'] <- 'N'

f_neg_diff <- upset(neg_uniq_diff, 
                    sets=c("s0013","s9242","s9316","s9346"),
                    query.legend = "top",
                    queries = list(list(query = elements, params = list("unique", "Y"), color = "#FF6633", active =T, query.name = "unique DAFs")),
                                   #list(query = intersects, params = list('s0013', 's9242', 's9316', 's9346'), active =T, query.name = "shared DAFs")),
                    keep.order = TRUE,
                    mb.ratio = c(0.65,0.35),
                    mainbar.y.label = "number of DAFs",
                    sets.x.label = "total DAFs, ESI-",
                    text.scale = c(2, 1.3, 1.3, 1.5, 2, 2), 
                    point.size = 3,
                    set_size.show = TRUE,
                    set_size.angles = 0,
                    sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97"),
                    #attribute.plots=list(gridrows=60,plots=list(list(plot=histogram, x="unique")))
)
f_neg_diff


####pos_up
pos_uniq_up <- read.csv(here("DAFs", "unique_DAFs", "pos_up_adduniq.list"), header = T, na.strings = "") %>% lapply(na.omit) 
pos_uniq_up <- fromList(pos_uniq_up)

pos_uniq_up$unique[pos_uniq_up$unique == '1'] <- 'Y' 
pos_uniq_up$unique[pos_uniq_up$unique == '0'] <- 'N'

f_pos_up <- upset(pos_uniq_up, 
                  sets=c("s0013","s9242","s9316","s9346"),
                  query.legend = "top",
                  queries = list(list(query = elements, params = list("unique", "Y"), color = "#FF6633", active =T, query.name = "ploidy-specfic DAFs")),
                                 #list(query = intersects, params = list('s0013', 's9242', 's9316', 's9346'), active =T, query.name = "shared DAFs")),
                  keep.order = TRUE,
                  mb.ratio = c(0.65,0.35),
                  mainbar.y.label = "number of DAFs",
                  sets.x.label = "up-regulated DAFs, ESI+",
                  text.scale = c(2, 1.3, 1.3, 1.5, 2, 2), 
                  point.size = 3,
                  set_size.show = TRUE,
                  set_size.angles = 0,
                  sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97"),
                  #attribute.plots=list(gridrows=60,plots=list(list(plot=histogram, x="unique")))
)
f_pos_up

####pos_down
pos_uniq_down <- read.csv(here("DAFs", "unique_DAFs", "pos_down_adduniq.list"), header = T, na.strings = "") %>% lapply(na.omit) 
pos_uniq_down <- fromList(pos_uniq_down)

pos_uniq_down$unique[pos_uniq_down$unique == '1'] <- 'Y' 
pos_uniq_down$unique[pos_uniq_down$unique == '0'] <- 'N'

f_pos_down <- upset(pos_uniq_down, 
                  sets=c("s0013","s9242","s9316","s9346"),
                  query.legend = "top",
                  queries = list(list(query = elements, params = list("unique", "Y"), color = "#FF6633", active =T, query.name = "unique DAFs")),
                  keep.order = TRUE,
                  mb.ratio = c(0.65,0.35),
                  mainbar.y.label = "number of DAFs",
                  sets.x.label = "down-regulated DAFs, ESI+",
                  text.scale = c(2, 1.3, 1.3, 1.5, 2, 2), 
                  point.size = 3,
                  set_size.show = TRUE,
                  sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97"),
                  #attribute.plots=list(gridrows=60,plots=list(list(plot=histogram, x="unique")))
)
f_pos_down

####pos_diff
pos_uniq_diff <- read.csv(here("DAFs", "unique_DAFs", "pos_diff_adduniq.list"), header = T, na.strings = "") %>% lapply(na.omit) 
pos_uniq_diff <- fromList(pos_uniq_diff)

pos_uniq_diff$unique[pos_uniq_diff$unique == '1'] <- 'Y' 
pos_uniq_diff$unique[pos_uniq_diff$unique == '0'] <- 'N'

f_pos_diff <- upset(pos_uniq_diff, 
                    sets=c("s0013","s9242","s9316","s9346"),
                    query.legend = "top",
                    queries = list(list(query = elements, params = list("unique", "Y"), color = "#FF6633", active =T, query.name = "unique DAFs")),
                    keep.order = TRUE,
                    mb.ratio = c(0.65,0.35),
                    mainbar.y.label = "number of DAFs",
                    sets.x.label = "total DAFs, ESI+",
                    text.scale = c(2, 1.3, 1.3, 1.5, 2, 2), 
                    point.size = 3,
                    set_size.show = TRUE,
                    set_size.angles = 0,
                    sets.bar.color = c("#20854E", "#E18727", "#7876B1", "#EE4C97"),
                    #attribute.plots=list(gridrows=60,plots=list(list(plot=histogram, x="unique")))
)
f_pos_diff

