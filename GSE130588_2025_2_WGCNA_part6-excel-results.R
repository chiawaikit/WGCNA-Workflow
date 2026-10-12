setwd("~/M5542_Bulk_RNA_Seq/GSE130588_WGNCA_05302025")
load("runi_GSE130588_WGNCA_qc_ardigen_analysis_parameters.RData")
setwd("~/M5542_Bulk_RNA_Seq/GSE130588_WGNCA_05302025/scripts/")
load("experimentName.recut.RData")
experimentDir<-paste0(analysisDir,experimentName.recut)
setwd(experimentDir)
load("final_results.DS3.RData")
geneInfo<-geneList.DS3.final.simple
colnames(geneInfo)[1]<-"moduleColor"
moduleCor<-moduleCor.DS3
rm(geneList.DS3.final.simple,moduleCor.DS3)

# you will need to make changes throughout this script

###############################################################################################################
# EXPORT MODULES TO EXCEL
options(java.parameters = "-Xmx16000m")
library(xlsx)
wb <- createWorkbook()

title.stem<-paste0(file.base,"_WGCNA_Final-Report_BulkRNASeq_2025-05-30")   ####CHANGE

##############################################
# CREATE SUMMARY PAGE

# Include at least:
# 1. GSE number / ArrayExpress number
# 2. Study title
# 3. PubMed ID (if referenced in a publication)
# 4. Type of sample
# 5. How sample was isolated
# 6. Timepoints in the study if there are more than one
# 7. What is being compared (e.g., Disease to healthy; disease to disease)
# 8. Whether or not there is in vitro stimulation of the cells
# 9. Numbers of samples: divided into healthy control and disease type if applicable
# 10. Ethnicity of patients
# 11. Clinical data (e.g., cohort, SLEDAI, drugs, etc.)
# 12. Platform used (e.g., Affy hgU133plus2)
# 13. ID of outliers (ie GSMxxxx)
# 14. Data normalization steps
# 15. Analysis files location
# 16. Version control (e.g., date of analysis and by whom)
# 17. Deepsplit used (i.e., 1, 2, 3, or 4)
# 18. Soft threshold used (e.g., 20)
# 19. Count of starting probes
# 20. Count of retained, analyzed probes
# 21. Correlation function type (i.e., Pearson, Bicor [biweight midcorrelation], or Spearman)
# 22. Signed or unsigned network
# 23. Were duplicates removed post-modularization, and if so, using what collapseRow method
# 24. Block size used in generating SFT and making TOMs

study <- data.frame(nrow=24, ncol=2)
study[1,] <- c("Study Numbers","NA")
study[2,] <- c("Study Title", "GSE130588_2025")
study[3,] <- c("Publication Status","Not Published")
study[4,] <- c("Cell Type","PBMC_MDDC")
study[5,] <- c("Cell Isolation","")
study[6,] <- c("Timepoints","Day_4")
study[7,] <- c("Disease Comparison","GSE130588 with disease activity score, treatment vs placebo and time of treatment")
study[8,] <- c("Cell Stimulation","")
study[9,] <- c("Cohorts","Treatment, None")
study[10,] <- c("Ethnicity","")
study[11,] <- c("Clinical Data","")
study[12,] <- c("Platform","Microarray")
study[13,] <- c("Outliers","None")
study[14,] <- c("Data Normalization","Raw counts data were log2 tranformed before being used as input into WGCNA")
study[15,] <- c("Files Location",".../")
study[16,] <- c("Version control","May 30th, 2025 by Wai-Kit Chia") 
study[17,] <- c("Deepsplit","4")
study[18,] <- c("Soft Threshold","10")
study[19,] <- c("Starting Features","15000+")
study[20,] <- c("Analyzed Features","33103")
study[21,] <- c("Correlation Function","Pearson")
study[22,] <- c("Network Type","Signed")
study[23,] <- c("Duplicates Removed?","NO")
study[24,] <- c("Block Size","10000")

ws <- createSheet(wb=wb, sheetName="WGCNA Summary")
addDataFrame(x=study, sheet=ws, row.names=FALSE,col.names=FALSE)
######################################################################
# NOW ADD Eigengene values from the desired deepsplit of choice
load(paste0(experimentDir,"/recut.DS-DS1-DS4.modules.and.genes.RData"))

me.table <- t(MEs.DS3)

ws <- createSheet(wb=wb, sheetName="Module Eigengenes")
addDataFrame(x=me.table, sheet=ws)

######################################################################
# ADD A TABLE OF CORRELATIONS FOR MODULE EIGENGENES AND CLINICAL TRAITS

ws <- createSheet(wb=wb, sheetName="Correlations")
moduleCor <- as.data.frame(moduleCor)
addDataFrame(x=moduleCor, sheet=ws, row.names=TRUE,col.names=TRUE)

######################################################################
# NOW ADD INDIVIDUAL MODULE SHEETS IN WORKBOOK

color_index <- levels(factor(geneInfo$moduleColor))
length(color_index)

library(dplyr)

geneInfo.list <- split(geneInfo, f = geneInfo$moduleColor) # Convert to a list object with elements based on module colors
data.list <- lapply(geneInfo.list,function(x){ # Remove all within module duplicates by taking the ones with highest kWithin
  as.data.frame(x) %>%
    group_by(symbol) %>% top_n(1,kWithin)})

library (plyr)

df <- ldply (data.list, data.frame) # Turn list back into a dataframe
colnames(df)[3]<-"symbol"
df <- df[,-1]
rownames(df) <- df$symbol #Now it looks just like geneInfo
#df <- df[!is.na(df$entrezgene),] # Remove all unannotated probes
df<-df[!is.na(df$symbol),]
#df<-df[-which(df$hgnc_symbol==""),]

for (i in 1:length(color_index) ) {
  single.mods <- df[ df$moduleColor==color_index[i], ]
  single.mods <- single.mods[ !is.na(single.mods$moduleColor), ]
  single.mods <- single.mods[ order(single.mods$MD_1_Tx.p), ] ######################MODIFY the primary clinical trait
  nrow(single.mods)
  ws <- createSheet(wb=wb, sheetName=paste(toupper(color_index[i])))
  addDataFrame(x=single.mods, sheet=ws, row.names=FALSE)
}

######################################################################
# SAVE THE WORKBOOK
setwd(experimentDir)
saveWorkbook(wb, paste0(title.stem,".xlsx"))
write.table(geneInfo, paste0(title.stem,".txt"))

write.table(geneInfo[,c(1,3)],quote = FALSE,sep = "\t" ,col.names = NA , file = "GSE130588_Module_Genes_DS_3_2025-05-30.txt")
write.table(df[,c(1,3)],quote = FALSE,sep = "\t" ,col.names = NA , file = "GSE130588_WGCNA_Module_Genes_NO_DUPS_DS_3_2025-05-30.txt")
