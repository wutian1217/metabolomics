set.seed(123)
library(here)
library(UpSetR)
library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggplotify)


#### Perform unique features selection for each column

####ESI-
log_negative <- read.csv(here("result_files","log_negative.csv"), header=T, row.names = 1, sep = ',')

rel_names<-colnames(log_negative)
basenames <- unique(trimws(rel_names, whitespace = "[_r1/_r2/_r3/_r4/_r5/_r6_r7/_r8]"))

group_patterns <- list(c(basenames[1], basenames[5]), 
                       c(basenames[6],basenames[3]), 
                       c(basenames[4],basenames[7]), 
                       c(basenames[8],basenames[2]))
all_unique_features <- data.frame(row_name = character(0), 
                                  zero_columns = character(0), 
                                  pattern = character(0), stringsAsFactors = FALSE)

row_names <- rownames(log_negative)


for (pattern_group in group_patterns){
  
  group1_pattern <- pattern_group[1]
  group2_pattern <- pattern_group[2]
  
  for (row_name in row_names) {
    group1_columns <- grep(group1_pattern, colnames(log_negative))
    group2_columns <- grep(group2_pattern, colnames(log_negative))
    
    group1_data <- log_negative[row_name, group1_columns]
    group2_data <- log_negative[row_name, group2_columns]
    
    if (all(group1_data == 0) & any(group2_data != 0)){
      print(paste("Condition statisfied for row:", row_name, "using pattern:", group1_pattern))
      zero_columns <- colnames(group1_data)[group1_data == 0]
      all_unique_features <- rbind(all_unique_features, data.frame(row_name = row_name, zero_columns = paste(zero_columns, collapse = ","), pattern = group1_pattern))
      
    } else if (any(group1_data != 0) & all(group2_data == 0)){
      print(paste("Condition statisfied for row:", row_name, "using pattern:", group2_pattern))
      zero_columns <- colnames(group2_data)[group2_data == 0]
      all_unique_features <- rbind(all_unique_features, data.frame(row_name = row_name, zero_columns = paste(zero_columns, collapse = ","), pattern = group2_pattern))
    } 
  }
}
print(all_unique_features)   

write.csv(all_unique_features, file=here("unique_DAFs","neg_any_unique_feature.csv"))

all_unique_features <- unique(all_unique_features)
write.csv(all_unique_features, file=here("unique_DAFs","neg_all_unique_feature.csv"))

####ESI+

log_positive <- read.csv(here("result_files","log_positive.csv"), header=T, row.names = 1, sep = ',')

rel_names<-colnames(log_positive)
basenames <- unique(trimws(rel_names, whitespace = "[_r1/_r2/_r3/_r4/_r5/_r6_r7/_r8]"))

group_patterns <- list(c(basenames[1], basenames[5]), c(basenames[6],basenames[3]), c(basenames[4],basenames[7]), c(basenames[8],basenames[2]))
all_unique_features <- data.frame(row_name = character(0), zero_columns = character(0), pattern = character(0), stringsAsFactors = FALSE)

row_names <- rownames(log_positive)

for (pattern_group in group_patterns){
  
  group1_pattern <- pattern_group[1]
  group2_pattern <- pattern_group[2]
  
  for (row_name in row_names) {
    group1_columns <- grep(group1_pattern, colnames(log_positive))
    group2_columns <- grep(group2_pattern, colnames(log_positive))
    
    group1_data <- log_positive[row_name, group1_columns]
    group2_data <- log_positive[row_name, group2_columns]
    
    if (all(group1_data == 0) & any(group2_data != 0)){
      print(paste("Condition statisfied for row:", row_name, "using pattern:", group1_pattern))
      zero_columns <- colnames(group1_data)[group1_data == 0]
      all_unique_features <- rbind(all_unique_features, data.frame(row_name = row_name, zero_columns = paste(zero_columns, collapse = ","), pattern = group1_pattern))
      
    } else if (any(group1_data != 0) & all(group2_data == 0)){
      print(paste("Condition statisfied for row:", row_name, "using pattern:", group2_pattern))
      zero_columns <- colnames(group2_data)[group2_data == 0]
      all_unique_features <- rbind(all_unique_features, data.frame(row_name = row_name, zero_columns = paste(zero_columns, collapse = ","), pattern = group2_pattern))
    } 
  }
}
print(all_unique_features)   
write.csv(all_unique_features, file=here("unique_DAFs", "pos_any_unique_feature.csv"))  

all_unique_features <- unique(all_unique_features)
write.csv(all_unique_features, file=here("unique_DAFs", "pos_all_unique_feature.csv"))
