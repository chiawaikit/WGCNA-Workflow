# renv::install("flashClust")
# renv::install("bioc::DEFormats")
# devtools::install_github("kassambara/ggpubr")
# renv::purge("Hmisc")
library(DESeq2)
library(edgeR)
library(DEFormats)
library(ggpubr)
library(dplyr)
library(plyr)
library(limma)
library(reshape)
library(xlsx)
library(plotly)
library(future)
library(RColorBrewer)
library(tibble)
library(ggforce)
library("stringr")
library(wacolors)
library(ggplot2)
library(NLP)
library(ggsignif)
require(RColorBrewer)
library(ggeasy)
library(DT)
library(arrow)
library(viridis)
library(WGCNA)
library(ProfilerAPI2)
library(tidyverse)
library(tidyselect)
library(tidyr)
library(arrow)

options(ggrepel.max.overlaps = Inf)



##### Authentication
# file.remove("~/.genedata-profiler/credentials")
# file.remove("~/.genedata-profiler/config")
api <- ProfilerAPI2::profiler_api()
# api$clear()
# api$authenticate()
# api$save_profile("default")
# api <- ProfilerAPI2::profiler_api(profile = "default")

# Set working directory
#### project parameters ####
# make sure all directories end in forwardslashes
globalDir<-"~/M5542_Bulk_RNA_Seq/GSE130588_WGNCA_05302025/"
projectDir<-paste0(globalDir,"")
scriptDir<-paste0(projectDir,"scripts/")
esetDir<-paste0(projectDir,"esets/")
file.base<-"runi_GSE130588_WGNCA_qc_ardigen" # for filenames etc.
title.base<-"GSE130588_WGNCA_2025" # for figure titles etc.
setwd(projectDir)
save(globalDir,projectDir,scriptDir,esetDir,file.base,title.base,file=paste0(file.base,"_project_parameters.RData"))
# getwd()
### Convert raw bulk RNAseq counts from XOP processed data into a DESeq object

### Load and prepare the expression data
setwd(esetDir)

# getwd()
GSE_dataset <- "GSE130588"
##### Dataset prep
# Four basic steps:
# 1. load and process count matrix. remove and record outliers
# 2. construct metadata. remove identified outliers from the metadata
# 3. test that the dataset can be loaded with the dataset mananger and that it works in the app
# 4. upload the prepared data files to Profiler for storage
# getwd()
# setwd("~/Th2_MLR_02122025")
##### Parameters
# Use the file IDs found in the source study folder within Genedata Profiler when gexp are stored
# assay.type <- "rnaseq-merck-uq" # "microarray", "rnaseq-merck-uq", "rnaseq-merck-deseq", "rnaseq-ardigen"
# If data are RNA-seq, supply the final normalized gexp view ID from Genedata Profiler

# Load dataset
# download_count_data <- api$data_lake$download_file("file-8afeea03-cb34-467e-b8d2-ffff0bc864ad", path = "~/M5542_Bulk_RNA_Seq/GSE130588_WGNCA_05212025/", overwrite = FALSE, show_progress = TRUE)
# m5542.count <- read_parquet("GSE130588_count_matrix.parquet")
# max(m5542.count$expression_value)
m5542.exprs <- read_parquet("GSE130588_expression_matrix.parquet")
m5542.exprs <- as.data.frame(m5542.exprs)
rownames(m5542.exprs) <- m5542.exprs[,1]
m5542.exprs <- m5542.exprs[,-1]
max(m5542.exprs)

# Load metadata
# download_metadata <- api$data_lake$download_file("file-7f88eaa8-109f-4560-8470-609a64d5ce01", path = "~/M5542_Bulk_RNA_Seq/GSE130588_WGNCA_05212025/", overwrite = FALSE, show_progress = TRUE)
m5542.meta <- read_parquet("GSE130588_metadata.parquet")
m5542.meta <- as.data.frame(m5542.meta)
colnames(m5542.meta)
rownames(m5542.meta) <- m5542.meta[,2]
all.equal(rownames(m5542.meta), colnames(m5542.exprs))
class(m5542.exprs)
all.equal(rownames(m5542.meta), colnames(m5542.exprs))


#Filter by study_id metadata if needed 
# gexp_view <- dplyr::tbl(api$conn, gexp_view_name)
# todays_date <- Sys.Date()
# output_filename_base <- "runi_AD_MLR_2025_gexpkallistobulkrnaseq"

#Choose appropriate control or other metadata category and replace NAs if necessary
# AD_meta <- xopcode::rnaseq_view_to_metadata(gexp_view)
# AD_tpm <- xopcode::rnaseq_view_to_matrix(gexp_view, data_type = "tpm")

### IMPORTANT: Transformations
# Per Horvath (WGCNA author), DESEQ2 normalized values should be VST'ed for use in WGCNA
# Similarly, log2 transformed normalized count data could also be used
# Therefore, VST or log2 transform your data!!!
# temp.tpm <- log2(m5542.exprs+1)
# temp.meta <- m5542.meta
temp.tpm <- m5542.exprs
temp.meta <- m5542.meta
max(temp.tpm)
# AD_count <- xopcode::rnaseq_view_to_matrix(gexp_view, data_type = "count")

## Get the column names of different fields
colnames(temp.tpm)

## Check row matches
all.equal(rownames(temp.meta), colnames(temp.tpm))
if(any(tolower(rownames(temp.meta)) != tolower(colnames(temp.tpm))))
  stop("colname mismatch")

# ## Data Formatting for expression data
# temp.tpm <- as.matrix(temp.tpm)
# class(temp.tpm)
# max(temp.tpm)
# 
# ## Data Formatting for metadata
# # temp.meta <- temp.meta[,c(2,4:ncol(temp.meta))]
# temp.meta <- temp.meta %>% distinct(sample_id,.keep_all = TRUE)
# temp.meta <- as.data.frame(temp.meta)
# temp.meta <- temp.meta[order(rownames(temp.meta)),]
# # rownames(temp.meta) <- temp.meta$sample_id_assay
# rownames(temp.meta) <- temp.meta$sample_id_assay
# temp.meta <- temp.meta[order(rownames(temp.meta)),]
# all.equal(colnames(temp.tpm), rownames(temp.meta))

# Concentration category
# temp.meta$treatment_time_activity_combo <- paste(temp.meta$treatment,temp.meta$time,sep = "_")
# temp.meta$treatment_time_activity_combo <- gsub("none_Week0","Normal_healthy",temp.meta$treatment_time_activity_combo)
# 
temp.meta$tissue_treatment_activity_combo <- paste(temp.meta$tissue,temp.meta$treatment,sep = "_")
temp.meta$tissue_treatment_activity_combo <- gsub("Normal_none","Normal_healthy",temp.meta$tissue_treatment_activity_combo)


# temp.meta$tissue_treatment_time_activity_combo <- paste(temp.meta$tissue_treatment_activity_combo,temp.meta$time,sep = "_")
# temp.meta$tissue_treatment_time_activity_combo <- gsub("Dupilumab","Dup",temp.meta$tissue_treatment_time_activity_combo)
# temp.meta$tissue_treatment_time_activity_combo <- gsub("Week","Wk",temp.meta$tissue_treatment_time_activity_combo)
# temp.meta$tissue_treatment_time_activity_combo <- gsub("placebo","plac",temp.meta$tissue_treatment_time_activity_combo)
# temp.meta$tissue_treatment_time_activity_combo <- gsub("Normal_healthy","Healthy_cont",temp.meta$tissue_treatment_time_activity_combo)

# Create numbering for treatment group: Dupilumab: 1; Placebo: 2
temp.meta$MD_1_Tx <- paste(temp.meta$treatment)
temp.meta$MD_1_Tx  <- gsub("Dupilumab", 1,temp.meta$MD_1_Tx)
temp.meta$MD_1_Tx  <- gsub("placebo",2,temp.meta$MD_1_Tx)

# Create numbering for time
temp.meta$MD_2_time <- paste(temp.meta$time)
temp.meta$MD_2_time <- gsub("Week0", 1, temp.meta$MD_2_time)
temp.meta$MD_2_time <- gsub("Week4", 2, temp.meta$MD_2_time)
temp.meta$MD_2_time <- gsub("Week16", 3, temp.meta$MD_2_time)

# Create numbering for race
temp.meta$MD_3_Race <- paste(temp.meta$subject_race)
temp.meta$MD_3_Race <- gsub("ASIAN", 1, temp.meta$MD_3_Race)
temp.meta$MD_3_Race <- gsub("BLACK OR AFRICAN AMERICAN", 2, temp.meta$MD_3_Race)
temp.meta$MD_3_Race <- gsub("WHITE", 3, temp.meta$MD_3_Race)
temp.meta$MD_3_Race <- gsub("OTHER", 4, temp.meta$MD_3_Race)

# Create numbering for sex
temp.meta$MD_4_Sex <- paste(temp.meta$subject_sex)
temp.meta$MD_4_Sex <- gsub("FEMALE", 1, temp.meta$MD_4_Sex)
temp.meta$MD_4_Sex <- gsub("MALE", 2, temp.meta$MD_4_Sex)



#Create subgroup of ls_Dupilumab and ls_placebo for comparison
temp.meta.subgroup <- temp.meta %>%
  filter(temp.meta$tissue_treatment_activity_combo %in% c("Ls_Dupilumab", "Ls_placebo"))
temp.tpm.subgroup <- temp.tpm[, colnames(temp.tpm) %in% temp.meta.subgroup$sample_id]
all.equal(rownames(temp.meta.subgroup), colnames(temp.tpm.subgroup))
# table(temp.meta.subgroup$treatment)

temp.meta.subgroup$MD_5_EASI <- round(temp.meta.subgroup$disease_activity)
temp.meta.subgroup$MD_6_SCORAD <- round(temp.meta.subgroup$sample_scorad)
all.equal(rownames(temp.meta.subgroup), colnames(temp.tpm.subgroup))

# Initialize the columns
temp.meta.subgroup$placebo_EASI <- NA
temp.meta.subgroup$placebo_SCORAD <- NA
temp.meta.subgroup$dupilumab_EASI <- NA
temp.meta.subgroup$dupilumab_SCORAD <- NA
# temp.meta.subgroup$placebo_EASI_time <- NA
# temp.meta.subgroup$placebo_SCORAD_time <- NA
# temp.meta.subgroup$dupilumab_EASI_time <- NA
# temp.meta.subgroup$dupilumab_SCORAD_time <- NA


# Loop through each row
for (i in 1:nrow(temp.meta.subgroup)) {
  if (temp.meta.subgroup$treatment[i] == "placebo") {
    temp.meta.subgroup$placebo_EASI[i] <- temp.meta.subgroup$MD_5_EASI[i]
    temp.meta.subgroup$placebo_SCORAD[i] <- temp.meta.subgroup$MD_6_SCORAD[i]
  } else if (temp.meta.subgroup$treatment[i] == "Dupilumab") {
    temp.meta.subgroup$dupilumab_EASI[i] <- temp.meta.subgroup$MD_5_EASI[i]
    temp.meta.subgroup$dupilumab_SCORAD[i] <- temp.meta.subgroup$MD_6_SCORAD[i]
  }
    else {
    # Assign NA for all other treatments
    temp.meta.subgroup$placebo_EASI[i] <- NA
    temp.meta.subgroup$placebo_SCORAD[i] <- NA
    temp.meta.subgroup$dupilumab_EASI[i] <- NA
    temp.meta.subgroup$dupilumab_SCORAD[i] <- NA
  }
}

# Check for empty fields
any(temp.meta.subgroup$placebo_EASI == "")
any(temp.meta.subgroup$placebo_SCORAD == "")
any(temp.meta.subgroup$dupilumab_EASI == "")
any(temp.meta.subgroup$dupilumab_SCORAD == "")


# Initialize the EASI_bin and SCORAD_bin columns
temp.meta.subgroup$EASI_bin <- NA  
temp.meta.subgroup$SCORAD_bin <- NA  

for (i in 1:nrow(temp.meta.subgroup)) {
  # Determine EASI_bin
  if (temp.meta.subgroup$disease_activity[i] < 1.0) {
    temp.meta.subgroup$EASI_bin[i] <- 1
  } else if (temp.meta.subgroup$disease_activity[i] <= 7.0 && 
             temp.meta.subgroup$disease_activity[i] >= 1.1) {
    temp.meta.subgroup$EASI_bin[i] <- 2
  } else if (temp.meta.subgroup$disease_activity[i] <= 21.0 && 
             temp.meta.subgroup$disease_activity[i] >= 7.1) {
    temp.meta.subgroup$EASI_bin[i] <- 3
  } else if (temp.meta.subgroup$disease_activity[i] <= 50.0 && 
             temp.meta.subgroup$disease_activity[i] >= 21.1) {
    temp.meta.subgroup$EASI_bin[i] <- 4
  } else {
    temp.meta.subgroup$EASI_bin[i] <- 5
  }
  
  # Determine SCORAD_bin
  if (temp.meta.subgroup$disease_activity[i] < 25.0) {
    temp.meta.subgroup$SCORAD_bin[i] <- 1
  } else if (temp.meta.subgroup$disease_activity[i] <= 50.0) {
    temp.meta.subgroup$SCORAD_bin[i] <- 2
  } else {
    temp.meta.subgroup$SCORAD_bin[i] <- 3
  }
}

temp.meta.subgroup$placebo_EASI_bin <- NA  # Initialize the column
temp.meta.subgroup$placebo_EASI_bin[temp.meta.subgroup$treatment == "placebo"] <- temp.meta.subgroup$EASI_bin[temp.meta.subgroup$treatment == "placebo"]

temp.meta.subgroup$placebo_SCORAD_bin <- NA  # Initialize the column
temp.meta.subgroup$placebo_SCORAD_bin[temp.meta.subgroup$treatment == "placebo"] <- temp.meta.subgroup$SCORAD_bin[temp.meta.subgroup$treatment == "placebo"]

temp.meta.subgroup$Dupilumab_EASI_bin <- NA  # Initialize the column
temp.meta.subgroup$Dupilumab_EASI_bin[temp.meta.subgroup$treatment == "Dupilumab"] <- temp.meta.subgroup$EASI_bin[temp.meta.subgroup$treatment == "Dupilumab"]

temp.meta.subgroup$Dupilumab_SCORAD_bin <- NA  # Initialize the column
temp.meta.subgroup$Dupilumab_SCORAD_bin[temp.meta.subgroup$treatment == "Dupilumab"] <- temp.meta.subgroup$SCORAD_bin[temp.meta.subgroup$treatment == "Dupilumab"]

temp.meta.subgroup$EASI_bin <- as.numeric(temp.meta.subgroup$EASI_bin)
temp.meta.subgroup$SCORAD_bin <- as.numeric(temp.meta.subgroup$SCORAD_bin)

# Check whether there is any NA or empty column
any(temp.meta.subgroup$EASI_bin == "" | is.na(temp.meta.subgroup$EASI_bin))
any(temp.meta.subgroup$SCORAD_bin == "" | is.na(temp.meta.subgroup$SCORAD_bin))


# temp.meta.subgroup$treatment_EASI_bin_combo <- paste(temp.meta.subgroup$treatment,temp.meta.subgroup$EASI_bin,sep = "_")
# temp.meta.subgroup$treatment_EASI_bin_time_combo <- paste(temp.meta.subgroup$treatment_EASI_bin_combo,temp.meta.subgroup$time,sep = "_")
# temp.meta.subgroup$treatment_SCORAD_bin_combo <- paste(temp.meta.subgroup$treatment,temp.meta.subgroup$SCORAD_bin,sep = "_")
# temp.meta.subgroup$treatment_SCORAD_bin_time_combo <- paste(temp.meta.subgroup$treatment_SCORAD_bin_combo,temp.meta.subgroup$time,sep = "_")

#### eset fields ####
setwd(esetDir)
# raw.counts.file <- temp.tpm
setwd(projectDir)
eset.name<-title.base
eset.lab<-"EMD Serono"
eset.contact<-"RU N&I: Discovery Biology"
eset.title<-"Atopic Dermatitis patients treated with Dupilumab or Placebo and Healthy control"
eset.abstract<-"See Chia-Chi Sun, Wai-Kit Chia, and Nick Geraci"
eset.url<-"emdserono.com"
eset.other<-list(notes="None")
eset.annotation<-paste(title.base," Public RNA microarray data generated in GSE130588, analyzed by EMD Serono",sep=" ")

pData <- temp.meta.subgroup

temp.tpm.filtered <- temp.tpm.subgroup %>%  as.data.frame() %>% dplyr::filter(rowSums(.) >= 1) # Only keep rows that have total counts above the cutoff
temp.tpm.filtered <- as.matrix(temp.tpm.filtered)

genes <- as.vector(rownames(temp.tpm.subgroup))
fData = data.frame(symbol=genes)
rownames(fData) <- fData$symbol

fData.meta <- data.frame(labelDescription = colnames(fData), row.names = colnames(fData) )
fData.eset <- new("AnnotatedDataFrame", data = data.frame(fData), varMetadata = fData.meta)

pData.meta <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData) )
pData.eset <- new("AnnotatedDataFrame", data = data.frame(pData), varMetadata = pData.meta)

# GENERATE ESET
expData <- new("MIAME", name = paste(eset.name,"RAW COUNTS",sep=" "),
               lab = eset.lab, 
               contact = eset.contact,
               title = eset.title,
               abstract = eset.abstract,
               url = eset.url,
               other = eset.other);

eset.GSE130588 <- ExpressionSet(assayData=temp.tpm.filtered,
                             phenoData = pData.eset,
                             # featureData = gene.annotations,
                             featureData = fData.eset,
                             experimentData = expData,
                             annotation = eset.annotation)

save(eset.GSE130588,file=paste0(esetDir,file.base,"_eset_UQnorm_log2TPM.RData"))

# GENERATE DENDROGRAM
library(flashClust)
# dev.off()
all.equal(rownames(temp.meta.subgroup),colnames(eset.GSE130588))
clus.data <- eset.GSE130588
colnames(clus.data) <- temp.meta.subgroup$sample_id
sampleTree <- flashClust(dist(t(exprs(clus.data))), method = "average");
sizeGrWindow(12,9)
pdf(file = "plot-sample-clustering.pdf", width = 12, height = 9);
par(cex = 0.9);
par(mar = c(2,6,4,2)) # bottom, left, top, right margins
plot(sampleTree, main = paste(title.base,"flashClust, average distance",sep=" "), sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 1.5,cex = 0.8)
dev.off()
