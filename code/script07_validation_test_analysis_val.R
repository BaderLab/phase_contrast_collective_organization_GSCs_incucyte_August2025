
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
path_to_repo <- '/Users/mystique27m/Documents/professional/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/'
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

source(paste0(script_path,"code/source_scripts/sourceData_General_variables_and_functions.R")) ### load general variables and functions

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


meta <- phase[1:5]

### plot the time dependent variation to complement S1B

granularity_names <- colnames(phase)[grep('Granularity', colnames(phase))]
txt_names <- setdiff(colnames(phase), c(colnames(meta),granularity_names))

#grn_df <- phase[granularity_names]
#txt_df <- phase[txt_names]

### z-normalize 
#grn_df <- as.data.frame(apply(grn_df, 2, standardise))
#grn_df <- cbind(meta, grn_df)

#grn_melt <- reshape2::melt(grn_df, value.name = 'grn_score', id.vars=colnames(meta))
#grn_melt$TimePt <- factor(grn_melt$TimePt, levels=timePt_levels)

#txt_df <- as.data.frame((apply(txt_df, 2, standardise)))
#txt_df <- cbind(meta, txt_df)

#txt_melt <- reshape2::melt(txt_df, value.name = 'txt_score', id.vars=colnames(meta))
#txt_melt$TimePt <- factor(txt_melt$TimePt, levels = timePt_levels)

#library(ggpubr)

#dir.create(paste0(figures_filepath, '/Sup_Fig_S1/'))



# Plot using ggline with mean lines and SE bars only

#grn_plots <- list()
#for (sample in 1:length(unique(grn_melt$Sample))){
#  grn_sample <- grn_melt[grn_melt$Sample == unique(grn_melt$Sample)[sample], ]

#  grn_plots[[sample]] <- ggline(grn_sample, 
#                     x = "TimePt", 
#                     y = "grn_score", 
#                     group = "variable",  # Group by variable to get separate lines for each
#                     add = "mean_se",  # Add mean with SE bars
#                     #facet.by = "variable",  # Facet by variable
#                     ncol = 1,  # Single column layout
#                     color = "variable",  # Color by variable for clarity
#                     plot_type = "l",  # 'l' stands for line plot without points
#                     size=.3
#  ) +
#    theme_minimal() +
#    theme(
#      legend.position = "right",  # Remove the legend
#      legend.key.size = unit(1.0, 'cm'),
#      legend.text = element_text(size=10),
#      legend.key.spacing.y = unit(1, 'pt'),
#      legend.title = element_blank(),
#      axis.title = element_blank(),
#      axis.text.x = element_blank()
#    ) +
#    labs(
#      title = unique(grn_melt$Sample)[sample],
#      #x = "Time Point",
#      #y = "Z-normalized Granularity Score ± SE"
#    )



#pdf(file = paste0(figures_filepath,'Sup_Fig_S1/S1_B_support_grn.pdf'), width=15, height=3.5)
#ggarrange(plotlist = grn_plots, 
#          ncol= 7,
#          nrow=1,
#          common.legend = TRUE,
#          legend = 'right')
#dev.off()

# Plot using ggline with mean lines and SE bars only

#txt_plots <- list()
#for (sample in 1:length(unique(txt_melt$Sample))){
#  txt_sample <- txt_melt[txt_melt$Sample == unique(txt_melt$Sample)[sample],]

#  txt_plots[[sample]] <- ggline(txt_sample, 
#                     x = "TimePt", 
#                     y = "txt_score", 
#                     group = "variable",  # Group by variable to get separate lines for each
#                     add = "mean_se",  # Add mean with SE bars
#                     #facet.by = "Sample",  # Facet by variable
#                     ncol = 1,  # Single column layout
#                     color = "variable",  # Color by variable for clarity
#                     plot_type = "l",  # 'l' stands for line plot without points
#                     size=.3
#  ) +
#    theme_minimal() +
#    theme(
#      legend.position = "right",  # Remove the legend
#      legend.text = element_text(size=10),
#      legend.key.size = unit(1.0, 'cm'),
#      legend.key.spacing.y = unit(1, 'pt'),
#      legend.title = element_blank(),
#      axis.title = element_blank(),
#      axis.text.x = element_blank()
#    ) +
#  labs(
#    title = unique(txt_melt$Sample)[sample],
#    #x = "Time Point",
#    #y = "Z-normalized Granularity Score ± SE"
#  )
  
  

#}  

#pdf(file = paste0(figures_filepath,'Sup_Fig_S1/S1_B_support_txt.pdf'), width=15, height=3.5)
#ggarrange(plotlist = txt_plots, 
#          ncol=7,
#          nrow=1,
#          common.legend = TRUE,
#          legend = 'right')
#dev.off()


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

#glist <- list()

#for (i in 1:length(df_list_tp_subsetted)){  
#      glist[[i]] <- df_list_tp_subsetted[[i]] %>% 
#        ggplot(aes(x=Sample, y=Area, fill=Sample)) + 
#        geom_boxplot() + 
#        geom_point(size=0.1) +
#        theme_minimal() + 
#        scale_fill_manual(values=colors_to_map) +
#        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=21),
#              axis.text.y = element_text(size=21),
#              axis.title.x = element_blank(),
#              axis.title.y = element_text(size=24), 
#              legend.position='none') + 
#        labs(x='Sample', y='Area') +
#        ylim(0,1)
#}


###------------------------------------------------------------------####
### FIGURE S2A
###------------------------------------------------------------------####
#dir.create(paste0(figures_filepath, '/Sup_Fig_S2_repeat/'))

#pdf(file = paste0(figures_filepath,'Sup_Fig_S2_repeat/S2_A.pdf'), width=42, height=6)
#ggarrange(plotlist=glist, ncol=9, nrow=1)
#dev.off()

###------------------------------------------------------------------####
### SAVE INTERMEDIATE DATA
###------------------------------------------------------------------####
#dir.create(paste0(out, 'script01_1_image_data_by_confluencies/'))
#new_path <- paste0(out, 'script01_1_image_data_by_confluencies/')
#new_path <- paste0(new_dir, 'C')
#for (l in 1:length(df_list_tp_subsetted)){
#  write.csv(df_list_tp_subsetted[[l]], paste0(new_path, l, '.csv'))
#}


#### save nromalized data and proceed with actual biological analysis
saveRDS(df_list_tp_subsetted, paste0(save_dir_path, '/df_list_tp_subsetted.rds'))


remove(list=ls())

