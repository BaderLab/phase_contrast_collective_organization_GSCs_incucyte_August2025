
### PACKAGES TO LOAD
###------------------------------------------------------------------#### 
library(ggplot2)  
library(colorspace)
library(tidyr)
library(dplyr)
library(ggthemes)
library(ggpubr)
library(ggrepel)
library(effectsize)
library(ggthemes)
library(scales)
library(forcats)


### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 

### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
### directory paths
home <- path.expand('~')
path_to_repo <- paste0(home,'/Documents/professional/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/')
dataset_path <- paste0(path_to_repo, 'incucyte_validation/cp_output_validationMarch2025/cp_complete_Results/') 
script_path <- paste0(path_to_repo,'/scripts_final/')
out <- paste0(path_to_repo, 'out/')
save_dir_path <- paste0(out, 'script07_output_files/')
dir.create(save_dir_path)

#figures_filepath <- paste0(out, 'figures_final/')
#dir.create(figures_filepath)

#tables_path <-  paste0(out,'tables/')
#dir.create(tables_path)

### LOAD DATA
###------------------------------------------------------------------####
phase <- read.csv(paste0(out, 'script06_output_files/final_dataset_afterQC.csv'), row.names = 1)

source(paste0(script_path,"scripts_Aug2025_to_push/code/source_scripts/sourceData_General_variables_and_functions.R")) ### load general variables and functions

#SINCE THERE SEEMS TO BE A DEPENDENCE ON THE AREA ON THE FEATURE SIGNALS, DECIDED TO GROUP IMAGES WITH SIMILAR AVAREAGE AREAS. THIS WAS DETERMINED AFTER TRYING OUT DIFFERENT WAYS OF SPLITTING THE IMAGES.
#THERE WERE 9 BINS , EACH WITH AN INCREMENT OF 10. FOR EXAMPLE THE FIRST BIN WAS ANY IMAGE SET BELONGING TO THE SAME TIME-POINT WITH A MEDIAN AREA OCCUPANCY FALLING WITHIN 0-10%.
#FIRST IMAGES FROM THE SAME SAMPLE WAS SUBSETTED. NEXT THE IMAGES WERE GROUPED TOGETHER BY TIME-POINTS. NEXT THE MEDIAN AREA WAS CALCULATED PER TIME-POINT WAS CALCULATED. FROM THIS, THE IMAGES FALLING WITHIN EACH AREA BIN IN EACH TIME-POINT WAS SELECTED AS A WHOLE  GROUPS WITH SIMILAR CONFLUENCY TOGETHER. 

### reduce samples to smaller set of features
phase <- phase[shorter]
phase <- phase[!phase$Sample %in% c('g752rp', 'g477', 'g626'),]


phase$Sample[phase$Sample=='g411'] = 'G411'
phase$Sample[phase$Sample=='g683'] = 'G683'
phase$Sample[phase$Sample=='g637'] = 'G637'
phase$Sample[phase$Sample=='g620'] = 'G620'

phase %>% group_by(Sample) %>% count(.)

meta <- phase[1:5]

### plot the time dependent variation to complement S1B

granularity_names <- colnames(phase)[grep('Granularity', colnames(phase))]
txt_names <- setdiff(colnames(phase), c(colnames(meta),granularity_names))

phase_orig <- read.csv2('/Users/mystique27m/Documents/professional/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/out/script01_output_files/final_dataset_afterQC.csv', sep=',')
phase_orig <- phase_orig[shorter]

phase <- rbind(phase_orig, phase)

phase <- phase %>% mutate(Area=as.numeric(Area))

df1 <- phase[!(phase$Sample %in% c('G566','G583')),] %>% as.data.frame(.)

### samples imaged at 12hrly intervals so this was processed separately.
df2 <- phase[(phase$Sample %in% c('G566','G583')),] %>% as.data.frame(.)

### Prepare cellprofiler data for normalization 
df_list_tp_subsetted <- list()

for (confluency in 1:9){
  ### subset data from phase dataset and then perform z-normalization
  tp_AreaMax <- data.frame()
  for (c in 1:length(levels(factor(df1$Sample)))){
    S <- levels(factor(df1$Sample))[c]
    sample <- df1[df1$Sample==S,]
    tp <- sample %>% group_by(TimePt) %>% summarise(Avg=median(Area))
    t <- subset(tp, Avg > confluency*0.1 & Avg < (confluency+1)*0.1)
    s <- subset(sample, TimePt %in% t$TimePt)
    tp_AreaMax <- rbind(tp_AreaMax, s)
  }
  
  for (c in 1:length(levels(factor(df2$Sample)))){
    S <- levels(factor(df2$Sample))[c]
    sample <- df2[df2$Sample==S,]
    tp <- sample %>% group_by(TimePt) %>% summarise(Avg=median(Area))
    t <- subset(tp, Avg > confluency*0.1 & Avg < (confluency+1)*0.1)
    s <- subset(sample, TimePt %in% t$TimePt)
    tp_AreaMax <- rbind(tp_AreaMax, s)
  }
  
  df_list_tp_subsetted[[confluency]] <- tp_AreaMax
  
}


### plot the distribution of the confluencies

sample_names <- unique(unlist(lapply(df_list_tp_subsetted, function(x) x$Sample)))
colors_to_map <- c(cols_vivid, col_bold)
names(colors_to_map) <- sample_names


#### save nromalized data and proceed with actual biological analysis
saveRDS(df_list_tp_subsetted, paste0(save_dir_path, '/df_list_tp_subsetted.rds'))


remove(list=ls())

