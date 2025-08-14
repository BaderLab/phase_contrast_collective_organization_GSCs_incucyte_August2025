###=========================== SET-UP 01 ============================###
#####################################################################

### PACKAGES REQUIRED
pkgs_reg <- c('SummarizedExperiment',
'GSVA',
'ggpubr',
'ggrepel',
'ggthemes',
'scales',
'tidyr',
'PCAtools',
'datawizard',
'forcats',
'effectsize',
'dplyr',
'ggplot2',
'colorspace',
'factoextra',
'NbClust',
'dendextend')

print("checking for packages if installed")
for (pkg in 1:length(pkgs_reg)){
  if(pkgs_reg[pkg] %in% rownames(installed.packages())==F){
    print(paste0(pkgs_reg[pkg], "needs to be installed"))
  }else{
      print('checked')
    }
}


###-------------------- VARIABLES & FUNCTIONS --------------------###
# COLOUR VECTORS
cols_vivid <- c('#E58606',
                         '#5D69B1',
                         '#52BCA3',
                         '#99C945',
                         '#CC61B0',
                         '#24796C',
                         '#DAA51B',
                         '#2F8AC4',
                         '#764E9F',
                         '#ED645A',
                         '#CC3A8E',
                         '#A5AA99',
                         '#A5B9D6',
                         '#F71B16',
                         '#AC9725')            
                         
col_bold <- c('#7F3C8D',
                       '#11A579',
                       '#3969AC',
                       '#F2B701',
                       '#E73F74',
                       '#80BA5A',
                       '#E68310',
                       '#008695',
                       '#CF1C90',
                       '#f97b72',
                       '#4b4b8f',
                       '#A5AA99')
                       
grad_colours <- divergingx_hcl(15, 'TealRose')


### Textural features reordered
phase.c.nameOrder<-c(
  'Drug',
  'Sample',
  'TimePt',
  'Well',
  'Area',
  'Granularity_1',
  'Granularity_2',
  'Granularity_3',
  'Granularity_4',
  'Granularity_5',
  'Granularity_6',
  'Granularity_7',
  'Granularity_8',
  'Granularity_9',
  'Granularity_10',
  'Granularity_11',
  'Granularity_12',
  'Granularity_13',
  'Granularity_14',
  'Granularity_15',
  'Granularity_16',
  'DifferenceEntropy_00',
  'DifferenceEntropy_01',
  'DifferenceEntropy_02',
  'DifferenceEntropy_03',
  'SumEntropy_00',
  'SumEntropy_01',
  'SumEntropy_02',
  'SumEntropy_03',
  'Entropy_00',
  'Entropy_01',
  'Entropy_02',
  'Entropy_03',
  'InfoMeas1_00',
  'InfoMeas1_01',
  'InfoMeas1_02',
  'InfoMeas1_03',
  'InfoMeas2_00',
  'InfoMeas2_01',
  'InfoMeas2_02',
  'InfoMeas2_03',
  'InverseDifferenceMoment_00',
  'InverseDifferenceMoment_01',
  'InverseDifferenceMoment_02',
  'InverseDifferenceMoment_03',
  'Angular2ndMoment_00',
  'Angular2ndMoment_01',
  'Angular2ndMoment_02',
  'Angular2ndMoment_03',
  'DifferenceVariance_00',
  'DifferenceVariance_01',
  'DifferenceVariance_02',
  'DifferenceVariance_03',
  'SumVariance_00',
  'SumVariance_01',
  'SumVariance_02',
  'SumVariance_03',
  'Variance_00',
  'Variance_01',
  'Variance_02',
  'Variance_03',
  'Correlation_00',
  'Correlation_01',
  'Correlation_02',
  'Correlation_03',
  'Contrast_00',
  'Contrast_01',
  'Contrast_02',
  'Contrast_03',
  'SumAverage_00',
  'SumAverage_01',
  'SumAverage_02',
  'SumAverage_03')

shorter<-c(
  'Drug',
  'Sample',
  'TimePt',
  'Well',
  'Area',
  'Granularity_1',
  'Granularity_2',
  'Granularity_3',
  'Granularity_4',
  'Granularity_5',
  'Granularity_6',
  'Granularity_7',
  'Granularity_8',
  'Granularity_9',
  'Granularity_10',
  'Granularity_11',
  'Granularity_12',
  'Granularity_13',
  'Granularity_14',
  'Granularity_15',
  'Granularity_16',
  'Entropy_00',
  'SumEntropy_00',
  'DifferenceEntropy_00',
  'Angular2ndMoment_00',
  'InverseDifferenceMoment_00',
  'SumAverage_00',
  'Variance_00',
  'SumVariance_00',
  'DifferenceVariance_00',
  'Correlation_00',
  'Contrast_00',
  'InfoMeas1_00',
  'InfoMeas2_00'
  )




###TIME VECTOR
timePt_levels <- c('T_0',
                   'T_4',
                   'T_8',
                   'T_12',
                   'T_16',
                   'T_20',
                   'T_24',
                   'T_28',
                   'T_32',
                   'T_36',
                   'T_40',
                   'T_44',
                   'T_48',
                   'T_52',
                   'T_56',
                   'T_60',
                   'T_64',
                   'T_68',
                   'T_72',
                   'T_76',
                   'T_80',
                   'T_84',
                   'T_88',
                   'T_92',
                   'T_96',
                   'T_100',
                   'T_104',
                   'T_108',
                   'T_112',
                   'T_116',
                   'T_120',
                   'T_124',
                   'T_128',
                   'T_132',
                   'T_136',
                   'T_140',
                   'T_144',
                   'T_148',
                   'T_152',
                   'T_156',
                   'T_160',
                   'T_164',
                   'T_168',
                   'T_172',
                   'T_176',
                   'T_180',
                   'T_184',
                   'T_188',
                   'T_192',
                   'T_196',
                   'T_200',
                   'T_204',
                   'T_208',
                   'T_212',
                   'T_216',
                   'T_220',
                   'T_224',
                   'T_228',
                   'T_232',
                   'T_236',
                   'T_240',
                   'T_244',
                   'T_248',
                   'T_252',
                   'T_256',
                   'T_260',
                   'T_264',
                   'T_268',
                   'T_272',
                   'T_276',
                   'T_280',
                   'T_284',
                   'T_288',
                   'T_292',
                   'T_296',
                   'T_300',
                   'T_304',
                   'T_308',
                   'T_312',
                   'T_316',
                   'T_320',
                   'T_324',
                   'T_328',
                   'T_332',
                   'T_336',
                   'T_340',
                   'T_344',
                   'T_348',
                   "T_352", 
                   "T_356", 
                   "T_360", 
                   "T_364", 
                   "T_368", 
                   "T_372", 
                   "T_376", 
                   "T_380", 
                   "T_384", 
                   "T_388", 
                   "T_392", 
                   "T_396",
                   "T_400")


###DEVELOPMENTAL-PROGENITOR/STEM CELLS
### Radial glial

rg_up <- c(
  "nowakowski_RG-div2_upreg",
  "nowakowski_MGE-RG1_upreg",       
  "nowakowski_MGE-RG2_upreg",      
  "nowakowski_RG-early_upreg",     
  "nowakowski_oRG_upreg",           
  "nowakowski_tRG_upreg",           
  "nowakowski_vRG_upreg",           
  "nowakowski_RG-div1_upreg")    

rg_down <- c(
  "nowakowski_RG-div2_downreg",     
  "nowakowski_MGE-RG1_downreg",     
  "nowakowski_MGE-RG2_downreg",     
  "nowakowski_RG-early_downreg",    
  "nowakowski_oRG_downreg",        
  "nowakowski_tRG_downreg",         
  "nowakowski_vRG_downreg",         
  "nowakowski_RG-div1_downreg")

gsc <- c(
"RNA.GSC.c1",                     
"RNA.GSC.c2",                    
"glioma.stem.cell",
"RNA_GSC_line_classifier",
"DEVELOPMENTAL_GSVA",
"INJURYRESPONSE_GSVA"
)

### NSCs
nsc_up <- c(
  "mizrak_aNSC1",                     
  "Neftel_NPC1", 
  "Neftel_NPC2", 
  "Zhong_NPCs_upreg",           
  "TCGA_PN")

nsc_down <- "Zhong_NPCs_downreg"


### Intermediate progenitor cells
ipcs_up <- c(
  "nowakowski_MGE-IPC1_upreg",      
  "nowakowski_MGE-IPC2_upreg",      
  "nowakowski_MGE-IPC3_upreg",      
  "nowakowski_IPC-nEN1_upreg",      
  "nowakowski_IPC-nEN2_upreg",      
  "nowakowski_IPC-div1_upreg",      
  "nowakowski_IPC-div2_upreg",     
  "nowakowski_IPC-nEN3_upreg")      

ipcs_down <- c(
  "nowakowski_MGE-IPC1_downreg",   
  "nowakowski_MGE-IPC2_downreg",    
  "nowakowski_MGE-IPC3_downreg", 
  "nowakowski_IPC-nEN1_downreg",   
  "nowakowski_IPC-nEN2_downreg",    
  "nowakowski_IPC-div1_downreg",    
  "nowakowski_IPC-div2_downreg",    
  "nowakowski_IPC-nEN3_downreg")  

### OPCs
opc_up <- c(
  "cahoy_OPC",                
  "Neftel_OPC",  
  "Zhong_OPC_upreg",            
  "nowakowski_OPC_upreg")           

opc_down <- c(         
  "Zhong_OPC_downreg",    
  "nowakowski_OPC_downreg")



### DEVELOPMENTAL DIFFERENTIATED
### astrocytes
astro_up <- c(
  "cahoy_astro_young",        
  "cahoy_astro_mature",       
  "cahoy_astro_in_vivo",      
  "cahoy_astrocyte",
  "mizrak_Astro1",                   
  "mizrak_Astro2",                    
  "mizrak_Astro3",                   
  "mizrak_Astro4",                    
  "mizrak_Astro5",  
  "mizrak_lateral_astro_M_upreg",     
  "mizrak_septal_astrocytes_M_upreg", 
  "mizrak_M_astro_lateral_upreg",     
  "mizrak_F_astro_septal_upreg",      
  "mizrak_lateral_astro_F_upreg",    
  "mizrak_septal_astrocytes_F_upreg", 
  "mizrak_M_astro_septal_upreg",      
  "mizrak_F_astro_lateral_upreg", 
  "Neftel_AC",   
  "Zhong_Astrocytes_upreg",     
  "nowakowski_Astrocyte_upreg",     
  "TCGA_CL")  

astro_down <- c(
  "Zhong_Astrocytes_downreg",   
  "nowakowski_Astrocyte_downreg")

### Neurons
neurons <- c(
  "cahoy_neuron",             
  "mizrak_Neuron1",                  
  "mizrak_Neuron2",                   
  "mizrak_Neuron3")         

interneurons_up <- c(
  "Zhong_Interneurons_upreg",   ### Interneurons
  "nowakowski_IN-STR_upreg",        
  "nowakowski_IN-CTX-CGE1_upreg",   
  "nowakowski_IN-CTX-CGE2_upreg",   
  "nowakowski_IN-CTX-MGE1_upreg",  
  "nowakowski_IN-CTX-MGE2_upreg",   
  "nowakowski_nIN1_upreg",          
  "nowakowski_nIN2_upreg",          
  "nowakowski_nIN3_upreg",          
  "nowakowski_nIN4_upreg",         
  "nowakowski_nIN5_upreg")          

interneurons_down <- c(
  "Zhong_Interneurons_downreg", 
  "nowakowski_IN-STR_downreg",     
  "nowakowski_IN-CTX-CGE1_downreg", 
  "nowakowski_IN-CTX-CGE2_downreg", 
  "nowakowski_IN-CTX-MGE1_downreg", 
  "nowakowski_IN-CTX-MGE2_downreg", 
  "nowakowski_nIN1_downreg",       
  "nowakowski_nIN2_downreg",        
  "nowakowski_nIN3_downreg",        
  "nowakowski_nIN4_downreg",        
  "nowakowski_nIN5_downreg")

excit.neurons_up <- c(
  "nowakowski_EN-PFC1_upreg",       ### Excitatory neurons
  "nowakowski_nEN-early2_upreg",    
  "nowakowski_nEN-late_upreg",     
  "nowakowski_EN-V1-1_upreg",       
  "nowakowski_EN-V1-2_upreg",       
  "nowakowski_EN-PFC2_upreg",       
  "nowakowski_nEN-early1_upreg",    
  "nowakowski_EN-PFC3_upreg",      
  "nowakowski_EN-V1-3_upreg")       

excit.neurons_down <- c(
  "nowakowski_EN-PFC1_downreg",     
  "nowakowski_nEN-early2_downreg",  
  "nowakowski_nEN-late_downreg",    
  "nowakowski_EN-V1-1_downreg",     
  "nowakowski_EN-V1-2_downreg",    
  "nowakowski_EN-PFC2_downreg",     
  "nowakowski_nEN-early1_downreg",  
  "nowakowski_EN-PFC3_downreg",     
  "nowakowski_EN-V1-3_downreg")     

### astroglia, mesenchymal
mes <- c(
  "cahoy_cultured_astroglia",
  "Neftel_MES2" ,
  "Neftel_MES1" ,
  "a1.astro" ,
  "a2.astro",
  "TCGA_MES")

### Oligodendrocytes
oligodendrocytes <- c(
  "cahoy_oligodendrocyte",    
  "cahoy_OL_myel",           
  "cahoy_MOG_pos")

#### Microglia
microglia_up <-  c(
  "mizrak_MicrogliaA",                
  "mizrak_MicrogliaB_1",             
  "mizrak_MicrogliaB_2", 
  "Zhong_Microglia_upreg",      
  "nowakowski_Microglia_upreg")     

microglia_down <-  c(   
  "Zhong_Microglia_downreg",  
  "nowakowski_Microglia_downreg")  

### Cycling
cycle <- c(
  "Neftel_G2/M", 
  "Neftel_G1/S",
  "nowakowski_MGE-div_upreg",       
  "nowakowski_MGE-div_downreg",
  "tirosh_s_phase",
  "tirosh_g2m") 

### Endothelial
endo <-  c(
  "mizrak_Endothelial1",              
  "mizrak_Endothelial2",             
  "mizrak_Endothelial3",
  "nowakowski_Endothelial_upreg",
  "mizrak_Pericyte",                  
  "mizrak_Fibroblast",                
  "mizrak_vSMC",                      
  "mizrak_Ependymal")             

### Others
other_up <- c(
  "mizrak_NB",                        
  "nowakowski_U1_upreg",            
  "nowakowski_U2_upreg",            
  "nowakowski_U3_upreg",            
  "nowakowski_U4_upreg",            
  "nowakowski_Mural_upreg",         
  "nowakowski_Glyc_upreg",         
  "nowakowski_Choroid_upreg",
  "Roulois_Interferon_Response")       

other_down <- c(
  "nowakowski_Endothelial_downreg",
  "nowakowski_U1_downreg",          
  "nowakowski_U2_downreg",          
  "nowakowski_U3_downreg",          
  "nowakowski_U4_downreg",          
  "nowakowski_Mural_downreg",       
  "nowakowski_Glyc_downreg",
  "nowakowski_Choroid_downreg")    


cell_types_up <- list(
  rg_up ,
  nsc_up ,
  gsc,
  ipcs_up ,
  opc_up,
  astro_up,  
  neurons,        
  interneurons_up, 
  excit.neurons_up,     
  mes,
  oligodendrocytes,
  microglia_up ,
  cycle,
  endo,
  other_up)

cell_types_down <- list(
  rg_down, 
  nsc_down,
  ipcs_down, 
  opc_down ,
  astro_down, 
  interneurons_down, 
  excit.neurons_down,     
  microglia_down, 
  other_down)


cell_type_names_up <- c(
  'rg_up' ,
  'nsc_up' , 'gsc', 'ipcs_up' ,
  'opc_up',
  'astro_up',  
  'neurons', 'interneurons_up', 'excit.neurons_up',     
  'mes',
  'oligodendrocytes',
  'microglia_up' ,
  'cycle',
  'endo',
  'other_up')

color_cellCLass_up <- c('#620069', # RG
                                 '#C700C5','#C700C5', '#C700C5', #GSC, #NSC, IPCS
                                 '#C700C5', #OPCs
                                 '#EC5D00', #Astro
                                 '#00938A','#1BAB85','#5EC276', #NEURONS, interneurons, excitatory neurons
                                 '#DDD7C6', # mes
                                 '#0093D7',   # oligodendrocytes
                                 '#FFD262', # Microglia 
                                 '#AC2500', # cell cycle 
                                 '#6E6E90', # endo
                                 '#6E6E90') # others
cell_type_names_down <- c(
  'rg_down', 
  'nsc_down',
  'ipcs_down', 
  'opc_down' ,
  'astro_down', 
  'interneurons_down', 
  'excit.neurons_down',     
  'microglia_down', 
  'other_down')

   
                                 
color_cellCLass_down <- c('#4C0017', '#69052C', '#862743', '#C25E76', '#3D516E', '#00866E', '#1E843A', '#5100B1', '#D01500')


###-------------------- QC FUNCTIONS ----------------------###

### Preparing datasets ###
prepare_data_for_Analysis <- function(data_image_CSV){
  QC_maskRaw = data_image_CSV %>% 
    select(contains('ImgQly')) %>% 
    select(contains('MskRaw'))
  
  QC_raw = data_image_CSV %>% 
    select(contains('ImgQly')) %>% 
    select(contains('cRaw'))
  
  total = data_image_CSV %>% select(!contains(c(colnames(QC_raw), colnames(QC_maskRaw)))) %>% select(contains(c('cRaw', 'MskRaw')))
  meta = data_image_CSV %>% select(!contains(c(colnames(total), colnames(QC_raw), colnames(QC_maskRaw))))
  maskRaw = total %>% select(matches('MskRaw'))
  raw = total %>% select(!contains(colnames(maskRaw)))
  
  raw = cbind(meta[c("Sample","TimePt", "Well", "Area", "Perimeter")], raw)
  maskRaw = cbind(meta[c("Sample","TimePt", "Well", "Area", "Perimeter")],maskRaw)
  QC_raw = cbind(meta[c("Sample","TimePt", "Well")],QC_raw)
  QC_maskRaw = cbind(meta[c("Sample","TimePt", "Well")],QC_maskRaw)
  
  data_list <- list(meta=meta, QC_raw=QC_raw, QC_maskRaw=QC_maskRaw, raw=raw, maskRaw=maskRaw)
}

### Adding drug name to dataset ###
drug_added_to_data <- function(dataset, druglist){
  
  drugs <- data.frame()
  samples <- data.frame()
  Well <- levels(factor(as.factor(dataset$Well)))
  for (w in 1:length(Well)){
    W <- Well[w]
    subset_data <- subset(dataset, subset=Well==W)
    n <- length(subset_data$Well)
    dw <- subset(druglist, subset=Well==W)
    d <- as.data.frame(lapply(dw[3], rep, n))
    drugs <- rbind(drugs,d)
    samples <- rbind(samples,subset_data)
  }
  names(drugs) <- 'Drug'
  dataset <- cbind(drugs, samples)
}

sc_features <-
c(
"ImageNumber",                                   
"ObjectNumber",                                  
"Sample",                                        
"Well",                                          
"time" ,                                        
"AreaShape_Center_X",                            
"AreaShape_Center_Y" ,                           
"AreaShape_Area"      ,                          
"AreaShape_BoundingBoxArea",                     
"AreaShape_BoundingBoxMaximum_X",               
"AreaShape_BoundingBoxMaximum_Y" ,               
"AreaShape_BoundingBoxMinimum_X"  ,              
"AreaShape_BoundingBoxMinimum_Y"   ,             
"AreaShape_CentralMoment_0_0",                   
"AreaShape_CentralMoment_0_1" ,                 
"AreaShape_CentralMoment_0_2"  ,                 
"AreaShape_CentralMoment_0_3"   ,                
"AreaShape_CentralMoment_1_0"    ,               
"AreaShape_CentralMoment_1_1"     ,              
"AreaShape_CentralMoment_1_2"      ,            
"AreaShape_CentralMoment_1_3"       ,            
"AreaShape_CentralMoment_2_0"        ,           
"AreaShape_CentralMoment_2_1"         ,          
"AreaShape_CentralMoment_2_2"          ,         
"AreaShape_CentralMoment_2_3"  ,
"AreaShape_Compactness"         ,                
"AreaShape_ConvexArea"           ,               
"AreaShape_Eccentricity"          ,              
"AreaShape_EquivalentDiameter"     ,             
"AreaShape_EulerNumber",                        
"AreaShape_Extent"      ,                        
"AreaShape_FormFactor"   ,                       
"AreaShape_HuMoment_0"    ,                      
"AreaShape_HuMoment_1"     ,                     
"AreaShape_HuMoment_2"      ,                   
#"AreaShape_HuMoment_3"                          
#"AreaShape_HuMoment_4"                          
#"AreaShape_HuMoment_5"                          
"AreaShape_HuMoment_6"       ,                   
"AreaShape_InertiaTensorEigenvalues_0",         
"AreaShape_InertiaTensorEigenvalues_1" ,         
"AreaShape_InertiaTensor_0_0",                   
"AreaShape_InertiaTensor_0_1" ,                  
"AreaShape_InertiaTensor_1_0"  ,                 
"AreaShape_InertiaTensor_1_1"   ,               
"AreaShape_MajorAxisLength"      ,               
"AreaShape_MaxFeretDiameter",                    
"AreaShape_MaximumRadius"   ,                    
"AreaShape_MeanRadius"       ,                   
"AreaShape_MedianRadius"      ,                 
"AreaShape_MinFeretDiameter"   ,                 
"AreaShape_MinorAxisLength"     ,                
"AreaShape_NormalizedMoment_0_2" ,               
"AreaShape_NormalizedMoment_0_3"  ,              
"AreaShape_NormalizedMoment_1_1"   ,            
#"AreaShape_NormalizedMoment_1_2"                
#"AreaShape_NormalizedMoment_1_3"                
"AreaShape_NormalizedMoment_2_0"  ,              
#"AreaShape_NormalizedMoment_2_1"                
#"AreaShape_NormalizedMoment_2_2"               
#"AreaShape_NormalizedMoment_2_3"                
#"AreaShape_NormalizedMoment_3_0"                
#"AreaShape_NormalizedMoment_3_1"                
#"AreaShape_NormalizedMoment_3_2"                
#"AreaShape_NormalizedMoment_3_3"  
"AreaShape_Orientation",                         
"AreaShape_Perimeter"   ,                        
"AreaShape_Solidity"     ,                       
"AreaShape_SpatialMoment_0_0",                   
"AreaShape_SpatialMoment_0_1" ,                 
#"AreaShape_SpatialMoment_0_2"                   
#"AreaShape_SpatialMoment_0_3"                   
#"AreaShape_SpatialMoment_1_0"                   
#"AreaShape_SpatialMoment_1_1"                   
#"AreaShape_SpatialMoment_1_2"                  
#"AreaShape_SpatialMoment_1_3"                   
#"AreaShape_SpatialMoment_2_0"                   
#"AreaShape_SpatialMoment_2_1"                   
"AreaShape_SpatialMoment_2_2"  ,                 
#"AreaShape_SpatialMoment_2_3"                  

"AreaShape_Zernike_0_0",                         
"AreaShape_Zernike_1_1" ,                        
"AreaShape_Zernike_2_0"  ,                       
"AreaShape_Zernike_2_2"   ,                      
"AreaShape_Zernike_3_1"    ,                    
"AreaShape_Zernike_3_3"     ,                    
"AreaShape_Zernike_4_0"      ,                   
#"AreaShape_Zernike_4_2"                         
#"AreaShape_Zernike_4_4"                         
"AreaShape_Zernike_5_1"       ,                 
"AreaShape_Zernike_5_3"        ,                 
"AreaShape_Zernike_5_5"         ,                
#"AreaShape_Zernike_6_0"                         
#"AreaShape_Zernike_6_2"                         
#"AreaShape_Zernike_6_4"                        
#"AreaShape_Zernike_6_6"                         
"AreaShape_Zernike_7_1",                         
"AreaShape_Zernike_7_3" ,                        
"AreaShape_Zernike_7_5"  ,                       
"AreaShape_Zernike_7_7"   ,                     
#"AreaShape_Zernike_8_0"                         
#"AreaShape_Zernike_8_2"                         
#"AreaShape_Zernike_8_4"                         
#"AreaShape_Zernike_8_6"                         
#"AreaShape_Zernike_8_8"                        
"AreaShape_Zernike_9_1"    ,                     
"AreaShape_Zernike_9_3"     ,                    
"AreaShape_Zernike_9_5"      ,                   
"AreaShape_Zernike_9_7"       ,                  
"AreaShape_Zernike_9_9"        ,                
"Texture_AngularSecondMoment_orig_3_00_256",     
"Texture_AngularSecondMoment_orig_3_01_256" ,    
"Texture_AngularSecondMoment_orig_3_02_256"  ,   
"Texture_AngularSecondMoment_orig_3_03_256"   ,  
"Texture_Contrast_orig_3_00_256"               ,
"Texture_Contrast_orig_3_01_256"                ,
"Texture_Contrast_orig_3_02_256"                ,
"Texture_Contrast_orig_3_03_256"                ,
"Texture_Correlation_orig_3_00_256"             ,
"Texture_Correlation_orig_3_01_256"            ,
"Texture_Correlation_orig_3_02_256"             ,
"Texture_Correlation_orig_3_03_256"             ,
"Texture_DifferenceEntropy_orig_3_00_256"       ,
"Texture_DifferenceVariance_orig_3_00_256"      ,
"Texture_DifferenceVariance_orig_3_01_256"      ,
"Texture_DifferenceVariance_orig_3_02_256"      ,
"Texture_DifferenceVariance_orig_3_03_256"     ,
"Texture_Entropy_orig_3_00_256"                 ,
"Texture_InfoMeas1_orig_3_00_256"              ,
"Texture_InfoMeas2_orig_3_00_256"               ,
"Texture_InverseDifferenceMoment_orig_3_00_256" ,
"Texture_InverseDifferenceMoment_orig_3_01_256" ,
"Texture_InverseDifferenceMoment_orig_3_02_256",
"Texture_InverseDifferenceMoment_orig_3_03_256" ,
"Texture_SumAverage_orig_3_00_256"              ,
"Texture_SumEntropy_orig_3_00_256"              ,
"Texture_SumVariance_orig_3_00_256"            ,
"Texture_SumVariance_orig_3_01_256"             ,
"Texture_SumVariance_orig_3_02_256"             ,
"Texture_SumVariance_orig_3_03_256"             ,
"Texture_Variance_orig_3_00_256"                ,
"Granularity_1_orig"                           ,
"Granularity_2_orig"                            ,
"Granularity_3_orig"                            ,
"Granularity_4_orig"                            ,
"Granularity_5_orig"                            ,
"Granularity_6_orig"                           ,
"Granularity_7_orig"                            ,
"Granularity_8_orig"  ,
"Granularity_9_orig",
"Granularity_10_orig",                           
"Granularity_11_orig" ,                          
"Granularity_12_orig"  ,                        
"Granularity_13_orig"   ,                        
"Granularity_14_orig"    ,                       
"Granularity_15_orig"     ,                      
"Granularity_16_orig")


