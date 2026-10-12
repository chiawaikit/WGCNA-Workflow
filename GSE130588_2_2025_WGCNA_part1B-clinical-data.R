setwd("~/M5542_Bulk_RNA_Seq/GSE130588_WGNCA_05302025")
load("runi_GSE130588_WGNCA_qc_ardigen_analysis_parameters.RData")

# Much of this script will need to be modified each time, no changes past the "end changes" section

# ENCODE CLINICAL DATA AS BINARY OR CONTINUOUS NUMERIC VALUES

### LOAD AND CONFIGURE LIBRARIES
library(Biobase)
library(WGCNA)
library(naniar)
library(dplyr)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

setwd(analysisDir)
# Encode binary designations and ensure numeric continuous values
load("part1-final-objects.RData")

# Load original pData from eset
pdata <- pData(eset.wgcna)
datTraits <- pdata
colnames(datTraits)
table(datTraits$MD_2_time)
# colnames(datTraits) <- c("sample_id_assay","experimental_perturbant_group")
#datTraits<-datTraits[,-c(1:3)]

### Factorize or make numeric all of the metadata desired
## For uninvolved skin, several categories were make NA to avoid sample pairing redundancy

datTraits$Ls_Dup_v_Ls_plac <- datTraits$MD_1_Tx
datTraits$Ls_Dup_v_Ls_plac <-gsub(1, 1,datTraits$Ls_Dup_v_Ls_plac)
datTraits$Ls_Dup_v_Ls_plac <-gsub("Ls_Dup_Wk0",2,datTraits$Ls_Dup_v_Ls_plac)
datTraits$Ls_Dup_v_Ls_plac <-as.numeric(datTraits$Ls_Dup_v_Ls_plac)
table(datTraits$Ls_Dup_v_Ls_plac)

# datTraits$Week0_v_Week4 <- datTraits$MD_2_time
# datTraits$Week0_v_Week4 <-gsub(1,1,datTraits$Week0_v_Week4)
# datTraits$Week0_v_Week4 <-gsub(2,2,datTraits$Week0_v_Week4)
# datTraits$Week0_v_Week4 <-gsub(3,NA,datTraits$Week0_v_Week4)
# datTraits$Week0_v_Week4 <-as.numeric(datTraits$Week0_v_Week4)
# table(datTraits$Week0_v_Week4)
# 
# datTraits$Week0_v_Week16 <- datTraits$MD_2_time
# datTraits$Week0_v_Week16 <-gsub(1,1,datTraits$Week0_v_Week16)
# datTraits$Week0_v_Week16 <-gsub(2,NA,datTraits$Week0_v_Week16)
# datTraits$Week0_v_Week16 <-gsub(3,2,datTraits$Week0_v_Week16)
# datTraits$Week0_v_Week16 <-as.numeric(datTraits$Week0_v_Week16)
# table(datTraits$Week0_v_Week16)




# treatment
# table(datTraits$tissue_treatment_time_activity_combo)
# # datTraits <- datTraits %>% mutate(treatment = replace(treatment, cohort == "skin_myositis_uninvolved", NA))
# 
# datTraits$tissue_treatment_time_activity_combo <-gsub("Healthy_cont_Wk0",1,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Ls_Dup_Wk0",2,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Ls_Dup_Wk4",3,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Ls_Dup_Wk16",4,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Ls_plac_Wk0",5,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Ls_plac_Wk4",6,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Ls_plac_Wk16",7,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Nl_Dup_Wk0",8,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Nl_Dup_Wk16",9,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Nl_plac_Wk0",10,datTraits$tissue_treatment_time_activity_combo)
# 
# datTraits$tissue_treatment_time_activity_combo<-gsub("Nl_plac_Wk16",11,datTraits$tissue_treatment_time_activity_combo)

# datTraits$experimental_perturbant_group<-gsub("Roca_12.5",11,datTraits$experimental_perturbant_group)
# datTraits$experimental_perturbant_group<-gsub("Roca_6.25",12,datTraits$experimental_perturbant_group)
# datTraits$experimental_perturbant_group<-gsub("Upa_12.5",13,datTraits$experimental_perturbant_group)
# datTraits$experimental_perturbant_group<-gsub("Upa_6.25",14,datTraits$experimental_perturbant_group)
# datTraits$experimental_perturbant_group<-gsub("None",15,datTraits$experimental_perturbant_group)
# datTraits$experimental_perturbant_group<-gsub("DMSO_6.25",16,datTraits$experimental_perturbant_group)
# datTraits$tissue_treatment_time_activity_combo<-as.numeric(datTraits$tissue_treatment_time_activity_combo)
# table(datTraits$tissue_treatment_time_activity_combo)



# datTraits<-datTraits[,-which(names(datTraits) %in% c("sample_id_assay","subject_id","sample_visit",
#                                                     "sample_collect_sequence","sample_collect_date",
#                                                     "sample_paired_tissues_available",
#                                                     "sample_paired_involved_tissue_available",
#                                                     "sample_disease",
#                                                     "subject_autoantibodies","subject_active_cancer",
#                                                     "subject_cancer_within_5yr","sample_accession",
#                                                     "subject_dmard_tx","sample_anatomical_entity","ID"))]
colnames(datTraits)

#datTraits <- datTraits[,c(1,2)]
datTraits <- datTraits[,c(13:29)]
identical(rownames(datTraits),rownames(datExpr)) # Just checking
#datTraits <- datTraits[ match(rownames(datTraits),rownames(datExpr)), ]



#### end changes ####
# There should be no more necessary changes past this point


datTraits.numeric <- data.matrix(datTraits)
datTraits <- data.frame(datTraits.numeric)
rm(datTraits.numeric)

# Plot dendrogram with clinical attributes
library(flashClust)
sampleTree <- flashClust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datTraits, signed = TRUE);

# Plot the sample dendrogram and the colors underneath

#par(mar = c(4,8,4,6)) # bottom, left, top, right margins
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = colnames(datTraits),marAll = c(1, 10, 1, 1),
                    cex.colorLabels = 0.8, cex.dendroLabels = 0.2, cex.rowText = 0.8,
                    main = paste(active.base,"WGCNA RNAseq log2TPM counts",sep=" "))

setwd(analysisDir)
# Overwrite the nascent frame created in part 1
save(datExpr, datTraits, eset.wgcna, badGenes, file="part1-final-objects.RData")

