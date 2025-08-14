### import gene signature scores from bulk transcriptome datasets of GSCs
bulk_Dev_I <- read.csv('/data/SourceData3_Bulk_IR_DEV.csv', header=T) ### THIS SOURCE DATA CAN BE DOWNLOADED FROM THE RICHARDS PAPER FROM NATURE CANCER 2021

DI_grad_bar <- data.frame(colorbar=c(rep(3,15)), colorlabel= grad_colours)
bulk_Dev_I <- bulk_Dev_I[c(1,5,6)]
bulk_Dev_I$SAMPLE_ID <- gsub('_L', '', bulk_Dev_I$SAMPLE_ID) 
bulk_Dev_I <- subset(bulk_Dev_I, subset=SAMPLE_ID %in% levels(factor(phase$Sample))) ### reduce to matched sample set

### G800 does not have IR/Dev scores so create a new row with 0s
G800 <- data.frame(Sample='G800', Developmental_AUC=0, InjuryResponse_AUC=0, DI_grad=0)

### substract Dev-IR scores and standardize them to get a single score
bulk_Dev_I$DI_grad <- standardize(bulk_Dev_I$DEVELOPMENTAL_GSVA)-standardize(bulk_Dev_I$INJURYRESPONSE_GSVA) 
bulk_Dev_I <- arrange(bulk_Dev_I, SAMPLE_ID)
DI_grad_bar_withD_IR <- cbind(bulk_Dev_I, DI_grad_bar)

