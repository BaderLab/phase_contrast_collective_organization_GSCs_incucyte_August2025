#SPLIT CELLPROFILER VARIABLES INTO GROUPS
#AREASHAPE
#TEXTURE
#GRANULARITY
areaShape <- colnames(scObj)[grep('AreaShape', colnames(scObj))]
areaShape2 <- c('AreaShape_Center_X', 'AreaShape_Center_Y') ### Re-arrange Area, Center_X and CEnter_Y - these parameters will not be used in analysis
areaShape <- c(areaShape2, setdiff(areaShape, areaShape2))
remove(areaShape2)

txt <- colnames(scObj)[grep('Texture', colnames(scObj))]
grn <- colnames(scObj)[grep('Granularity', colnames(scObj))]
metaData <- c('ImageNumber', 'ObjectNumber', 'Metadata_Sample','Metadata_Well', 'Metadata_time')

scObjF <- scObj[c(metaData, areaShape, txt, grn)]
cnames <- colnames(scObjF)
cnames <- gsub('Metadata_', '', cnames)
colnames(scObjF) <- cnames
#scObjF <- scObjF[c(1:5,7,8,6,9:ncol(scObjF))]
counts_scObjF <- dplyr::count(scObjF, Sample, Well, time)

#CHECK WHICH COLUMNS HAVE NA AND THE PERCENTAGE OF NA DATA POINTS PER COLUMN WITH NA
### compute NA by sample 

na_final <- as.data.frame(matrix(ncol=ncol(scObjF[8:ncol(scObjF)]), nrow=0))
colnames(na_final) <- colnames(scObjF[8:ncol(scObjF)])

for (s in 1:10){
  label <- levels(factor(scObjF$Sample))[s]
  df <- subset(scObjF, Sample==label)
  na <- ((colSums(is.na(df[8:ncol(df)]))/nrow(df[8:ncol(df)]))*100) %>% t(.) %>% as.data.frame(.)
  na_final <- rbind(na_final, na)
}
pheatmap::pheatmap(na_final,cluster_rows = F)


scObjF2 <- scObjF[,setdiff(colnames(scObjF),c("AreaShape_NormalizedMoment_0_0", "AreaShape_NormalizedMoment_0_1", "AreaShape_NormalizedMoment_1_0"))]
scObjF2 <- na.omit(scObjF2)

na_final <- as.data.frame(matrix(ncol=ncol(scObjF2[8:ncol(scObjF2)]), nrow=0))
for (s in 1:10){
  label <- levels(factor(scObjF2$Sample))[s]
  df <- subset(scObjF2, Sample==label)
  na <- ((colSums(is.na(df[8:ncol(df)]))/nrow(df[8:ncol(df)]))*100) %>% t(.) %>% as.data.frame(.)
  na_final <- rbind(na_final, na)
}

is.na(scObjF2) %>% sum(.)

#Remove features that are all NAs
scObjF_shortened <- scObjF2


colours_cor <- divergingx_hcl(21, 'Tropic')
colours_cor <- sequential_hcl(51, 'Rocket', rev=F)
z <-colnames(scObjF_shortened)[grep('Zernike', colnames(scObjF_shortened))]
scObjF_shortened_z <- scObjF_shortened[z] #%>% cor(.)

scObjF_shortened_as <- scObjF_shortened[colnames(scObjF_shortened) %in% areaShape] 
scObjF_shortened_as <- scObjF_shortened_as[!(colnames(scObjF_shortened_as) %in% z)] %>% .[3:ncol(.)] #%>% cor(.)
scObjF_shortened_txt <- scObjF_shortened[txt]
scObjF_shortened_grn <- scObjF_shortened[grn] %>% apply(.,2,FUN=function(x)(x-min(x))/(max(x)-min(x))) #%>% cbind(scObjF_shortened[1:5], .)
scObjF_shortened_txt_grn <- #scObjF_shortened_txt[grep('3_00', colnames(scObjF_shortened_txt))] %>% 
  apply(scObjF_shortened_txt,2,FUN=function(x)(x-min(x))/(max(x)-min(x))) %>% 
  cbind(scObjF_shortened[1:5], .) %>% cbind(., scObjF_shortened_grn)

#scObjF_shortened_as <- as.matrix(scObjF_shortened_as)
#scObjF_shortened_txt <- as.matrix(scObjF_shortened_txt)
#scObjF_shortened_grn <- as.matrix(scObjF_shortened_grn)
scObjF_shortened_as[is.na(scObjF_shortened_as)] <- 0