setwd("~/GSE130588-WGCNA-Workflow")
load("runi_GSE130588_WGNCA_qc_ardigen_project_parameters.RData")
active.base<-paste(title.base, "WGCNA", sep="_") # for figure titles etc. Could be Active Female, Male, Inactive...
analysisDir<- "~/GSE130588-WGCNA-Workflow/scripts/" # location of the current script, make sure it ends in a forwardslash
save(globalDir,projectDir,scriptDir,esetDir,analysisDir,file.base,title.base,active.base,file=paste0(file.base,"_analysis_parameters.RData"))
setwd(analysisDir)
# There should be no more necessary changes past this point

################################################################################################
### WGCNA DATA INPUT GENERATION
# The variance stabilizing transformed eset is loaded, counts are extracted to a datTraits
# frame, and clincal data to a datExpr frame. Note the counts frame must be transposed for
# input into WGCNA, where column names are probes keyed to the fData frame in the eset, and
# row names are sample IDs keyed to those in the eset pData frame and the newly generated
# datExpr frame. 
#
# datExpr requires additional encoding before correlation with downstream 
# module eigengenes, the encoding being completed in a separate script which can be performed
# while datTraits is used for generation of a transformed counts adjacency matrix and encoded 
# as a topological overlay matrix (TOM).
#
# CRITICAL SAMPLE FILTERING NOTE!
# This experiment only includes ACTIVE FEMALES
#
################################################################################################

### LOAD AND CONFIGURE LIBRARIES
library(Biobase)
library(limma)
library(WGCNA)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

################################################################################################
### Create the datExpr0 frame for WGCNA
setwd(esetDir)
load(paste0(file.base,"_eset_UQnorm_log2TPM.RData"))
eset <- eset.GSE130588
rm(eset.GSE130588)
fData <- fData(eset)

datExpr0 <- exprs(eset)
nrow(datExpr0)
ncol(datExpr0)

# IMPORTANT! TRANSPOSE THE ORIGINAL EXPRESSION MATRIX
# WGCNA requires columns be probes and rows be arrays
datExpr0 <- t(datExpr0)

# Check for gene entries with NAs and very low counts
gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

# Display and remove these bad entries
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    badGenes <- colnames(datExpr0)[!gsg$goodGenes];
  printFlush(paste("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
} else {
  badGenes<-NULL
}

# We must remove probe columns containing only one value. Failure to do so breaks downstream TOM creation,
# and it must be performed now before module colors are assigned.
datExpr1 <- datExpr0[,apply(datExpr0,2,function(x) any(c(FALSE,x[-length(x)]!=x[-1])))]
ncol(datExpr0) - ncol(datExpr1)
datExpr0 <- datExpr1

## Remove bad probes from the eset eliminated by good genes function. 07-19-17 RDR: 312 genes are removed
datExpr1 <- t(datExpr0)
eset <- eset[rownames(eset)%in%rownames(datExpr1),]
nrow(eset)==nrow(datExpr1)
rm(datExpr1,gsg)

# Save the modified eset which will accompany this WGCNA run
setwd(analysisDir)
eset.wgcna <- eset
rm(eset)
save(eset.wgcna,file=paste0(file.base,"_eset_VST_active_WGCNA.RData"))
eset <- eset.wgcna

#######################################################################################
# Plot PCA AFTER WGCNA FILTERING

genevars <- apply(eset, 1, var)
genevar_cutoff <- sort(genevars, decreasing=T)[1000]
keep <- genevars >= genevar_cutoff
pt <- ggtitle("PCA by Top 1000 Most Variable Genes")
library(PCAtools)
p <- PCAtools::pca(exprs(eset)[keep,], metadata = pData(eset))
pdf(file = paste0("pre_wgcna_pca.pdf"))
PCAtools::biplot(pcaobj = p,
                 colby = "sample_id",
                 # colkey = c("dermatomyositis" = 'red', "healthy" = 'blue'),
                 shape = "tissue_treatment_time_activity_combo",
                 shapekey = c("Healthy_cont_Wk0" = 1,
                              "Ls_Dup_Wk0" = 2,
                              "Ls_Dup_Wk4" = 3,
                              "Ls_Dup_Wk16" = 4,
                              "Ls_plac_Wk0" = 5,
                              "Ls_plac_Wk4" = 6,
                              "Ls_plac_Wk16" = 7,
                              "Nl_Dup_Wk0" = 8,
                              "Nl_Dup_Wk16" = 9,
                              "Nl_plac_Wk0" = 10,
                              "Nl_plac_Wk16" = 11),
                 legendPosition = "right",
                 lab = NULL,
                 pointSize = 2,
                 hline = 0, vline = 0) + pt
dev.off()

# As before, RGL 3D plots can only be run on a desktop version of RStudio
# plotPCA(eset, groups=as.numeric(eset$cohort), groupnames=levels(eset$cohort),
#         pcs=c(1,2,3), plot3d=TRUE, col=c("red","blue"), pch=c(16,17,17))

#######################################################################################
# Draw a quick dendrogram to look for outliers in the filtered set
# 07-18-17 RDR: The dendrogram is also identical
library(flashClust)
dev.off()
sampleTree <- flashClust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.9);
par(mar = c(2,6,4,2)) # bottom, left, top, right
plot(sampleTree, main = paste(active.base,"WGCNA flashclust post-var filter",sep=" "), sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 1.5)
dev.off()

datExpr <- datExpr0
nGenes <- ncol(datExpr) # number of columns (probes)
nSamples <- nrow(datExpr) # number of rows (arrays)
rm(sampleTree, datExpr0)

###########################################################################################
# SAVE DATEXPR FRAME
# DATTRAITS encoding can occur during TOM generation

###########################################################################################
## ENCODE CLINICAL TRAIT DATA (Applies to both filtered and unfiltered expression sets)

# Create a patients data frame analogous to expression data that will hold binary designations
# of trait membership, or a numerical measurement of a continuous clinical variable
datTraits <- pData(eset)
identical(rownames(datTraits),colnames(exprs(eset)))
colnames(datTraits)

# COMPLETE REMAINING CLINICAL DATA ENCODING IN A SUBSEQUENT SCRIPT DURING TOM GENERATION


###########################################################################################
### Clean up the dataset and display final dendrogram with corresponding clinical factors
collectGarbage();

###########################################################################################
### SAVE OBJECTS 
setwd(analysisDir)
save(datExpr, datTraits, eset.wgcna, badGenes, file="part1-final-objects.RData")

