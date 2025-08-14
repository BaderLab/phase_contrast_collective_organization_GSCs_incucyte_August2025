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


### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
### directory paths
path_to_repo <- '/Users/mystique27m/Documents/professional/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/'
dataset_path <- paste0(path_to_repo, 'incucyte_validation/cp_output_validationMarch2025/cp_complete_Results/') 
script_path <- paste0(path_to_repo,'/scripts_final/')
out <- paste0(path_to_repo, 'out/')
save_dir_path <- paste0(out, 'script06_output_files/')
dir.create(save_dir_path)

### file paths
#well_filepath <- paste0(path_to_repo,'datasets_fianl/wells_to_choose.csv') 


### source data 
###------------------------------------------------------------------#### 
source(paste0(script_path,"code/source_scripts/sourceData_General_variables_and_functions.R"))

### Original data
###------------------------------------------------------------------#### 
# read cellprofiler results
d1 <- read.csv(paste0(dataset_path,"/all_mediaOnly_fullDatasetImage.csv"))
d2 <- read.csv(paste0(dataset_path,"/G620_mediaOnly_fullDatasetImage.csv"))
d3 <- read.csv(paste0(dataset_path,"/G637_G683_mediaOnly_fullDatasetImage.csv"))

# read drug list
#wells <- read.csv(well_filepath)

setdiff(colnames(d3), colnames(d1))

### drop processing status from d2 and d3
d2 <- d2[!colnames(d2) %in% "ProcessingStatus"]
d3 <- d3[!colnames(d3) %in% "ProcessingStatus"]

#--------------------------------------------------------------------------------#
### --------------------------------- CODE START ------------------------------###
#--------------------------------------------------------------------------------#

imageData <- rbind(d1,d2,d3)


#g752 <- imageData[grepl("g752", imageData$FileName_Raw), ]

#imageData <- imageData[!grepl("g752", imageData$FileName_Raw), ]

#filenames <- unique(g752$FileName_Raw)
#filenames[1]


imageData$Metadata_Sample <- sub("_.*", "", imageData$FileName_Raw)  # Everything before first "_"
imageData$Metadata_Well <- sub(".*_([A-Z][0-9]+)_.*", "\\1", imageData$FileName_Raw)  # Extract well
imageData$Metadata_TimePt <- sub(".*_(T_[0-9]+)\\.png$", "\\1", imageData$FileName_Raw)  # Extract timepoint

imageData_verify <- imageData[c('FileName_Raw', 'Metadata_Sample', 'Metadata_Well', 'Metadata_TimePt')]

imageData <- subset(imageData, subset=Metadata_Sample != 'g304')
# remove objects once combined
rm(d1,d2,d3)


###RENAME FEATURE NAMES TO MAKE IT TIDIER AND NEATER
### Make column names shorter
colN <- colnames(imageData)
colN <- gsub('AreaOccupied_', 'Area', colN)
colN <- gsub('Threshold', 'Thr', colN)
colN <- gsub('ExecutionTime', 'ExecT', colN)
colN <- gsub('ImageQuality', 'ImgQly', colN)
colN <- gsub('LocalFocusScore', 'LFS', colN)
colN <- gsub('Intensity', 'Int', colN)
colN <- gsub('Texture_', '', colN)
colN <- gsub('AngularSecondMoment', 'Angular2ndMoment', colN)
colN <- gsub('MaskCropRaw', 'MskRaw', colN)
colN <- gsub('CropRaw', 'cRaw', colN)
colN <- gsub('Metadata_', '', colN)
colN <-gsub('_3_', '_', colN)
colN <- gsub('_256', '', colN)
colN <- gsub('Granularity_cRaw', 'Granularity_3_cRaw', colN)
colN <- gsub('Granularity_MskRaw', 'Granularity_3_MskRaw', colN)

colN[1:3] <- c('Area', 'Perimeter', 'TotalArea')
names(imageData) <- colN

###PREPARE DATA FOR ANALYSIS
###SEPARATE IMAGES ANALYSED FOR QUALITY CONTROL FEATURES (WITH AND WITHOUT MASKING)
final_df <-prepare_data_for_Analysis(imageData)

### masked data
df_mask <- final_df[["maskRaw"]]

#df_mask <- drug_added_to_data(df_mask, drug.list)
df_mask <- cbind(as.data.frame(rep('media', nrow(df_mask))), df_mask)
colnames(df_mask) <- c('Drug', colnames(final_df$maskRaw))

### raw data
df_raw <- final_df$raw
#df_raw <- drug_added_to_data(df_raw, wells)

df_mask <- df_mask %>% group_by(Sample) %>% 
  arrange(Sample, group_by=T) %>% 
  mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>%
  as.data.frame(.)


### QC from masked Raw images
qc_maskRaw <- final_df$QC_maskRaw
qc_maskRaw <- qc_maskRaw %>% group_by(Sample) %>% 
  arrange(Sample, group_by=T) %>% 
  mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>%
  as.data.frame(.)

### we remove images that are of low quality using the power log spectrum as described in materials & methods
### we then plot the quality of each image by patient sample and remove the low quality images with low pplog values 


###REMOVE LOW QUALITY IMAGES USING RAW UNMASKED IMAGES
### QC from Raw images
qc_Raw <- final_df$QC_raw
qc_Raw <- qc_Raw %>% group_by(Sample) %>% 
  arrange(Sample, group_by=T) %>% 
  mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>%
  as.data.frame(.)

n=15
g1 <- 
  ggplot(qc_Raw, aes(x=Sample, y=ImgQly_PowerLogLogSlope_cRaw, fill=Sample)) + 
  scale_fill_manual(values=c(col_bold, cols_vivid)) +
  #geom_boxplot(size=.3) +
  geom_violin() +
  geom_jitter(size=0.1, alpha=0.1) +
  scale_y_continuous(breaks=scales::breaks_width(0.1)) +
  theme_linedraw() + 
  theme(axis.text.x=element_text(size=9, face='bold', angle=30), 
        axis.text.y=element_text(size=6, angle=30),
        strip.text = element_text(size=12, colour='#80BA5A', face='bold'),
        legend.position = "none")

g1
qc_Raw_cleaned <- data.frame()


for (i in 1:15){
  levels <- levels(factor(qc_Raw$Sample))
  sample <- qc_Raw[qc_Raw$Sample==levels[i],]
  cleaned <- subset(sample, subset=!is.na(sample$ImgQly_PowerLogLogSlope_cRaw))
  cleaned <- subset(cleaned, subset=cleaned$ImgQly_PowerLogLogSlope_cRaw > quantile(cleaned$ImgQly_PowerLogLogSlope_cRaw)[2], na.rm=T)
  qc_Raw_cleaned <- rbind(qc_Raw_cleaned, cleaned)
}

g2 <- 
  ggplot(qc_Raw_cleaned, aes(x=Sample, y=ImgQly_PowerLogLogSlope_cRaw, fill=Sample)) + 
  scale_fill_manual(values=c(col_bold, cols_vivid)) +
  geom_violin() +
  geom_jitter(size=0.1, alpha=0.1) +
  ylim(-1.2, 2.15) +  
  theme_linedraw() + 
  theme(axis.text.x=element_text(size=9, face='bold', angle=30), 
        axis.text.y=element_blank(), 
        strip.text = element_text(size=12, colour='#80BA5A', face='bold'),
        legend.position = "none")

g2
### Figure was included in manuscript but described in paper

#if(dir.exists(paste0(figures_filepath,'Sup_Fig_S4'))==F){
#  dir.create(paste0(figures_filepath,'Sup_Fig_S4'))}

#pdf(file = paste0(figures_filepath,'Sup_Fig_S4/S4_A.pdf'), width=9, height=4 )
#ggarrange(g1,g2, nrow=1)
#dev.off()
### figure end



### OUTPUT DATA FILE
###------------------------------------------------------------------###
            # S1_outPutDataFile - SAVING CLEANED DATA #
###------------------------------------------------------------------###
final_qc_Mask_cleaned <- subset(qc_maskRaw, subset=rownames(qc_maskRaw) %in% rownames(qc_Raw_cleaned))
colN <- colnames(final_qc_Mask_cleaned) 
colN <- gsub('_MskRaw', '', colN) ### 
colnames(final_qc_Mask_cleaned) <- colN

write.csv(final_qc_Mask_cleaned, paste0(save_dir_path, "/final_qc_Mask_cleaned.csv"))


final_dataset_afterQC <- subset(df_mask, subset=rownames(df_mask) %in% rownames(qc_Raw_cleaned)) ### subset images with ilastik masks using whole images that passed QC from above
colN <- colnames(final_dataset_afterQC) 
colN <- gsub('_MskRaw', '', colN) ### 
colnames(final_dataset_afterQC) <- colN
qc_failed <- subset(df_mask, subset= rownames(df_mask) %in% setdiff(rownames(df_mask), rownames(qc_Raw_cleaned)))[2:4] %>% .[c(1,3,2)]

### Save QC failed image dataset
write.csv(qc_failed, paste0(save_dir_path, "/mediaOnly_failed_QC.csv"))
### output datafile save end


A=1066104 ### size of full image in pixel.sq
n=15 # sample number

### load imaging datasets
phase <- final_dataset_afterQC
phase <- phase[phase.c.nameOrder] ### reorder columns

### Save QC passed image dataset
s <- levels(factor(phase$Sample))
phase <- phase %>% mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>% arrange(TimePt)
phase$Area <- phase$Area/A #convert absolute area to fraction of total image area
write.csv(phase, paste0(save_dir_path, "/final_dataset_afterQC.csv"))


remove(list=ls())

