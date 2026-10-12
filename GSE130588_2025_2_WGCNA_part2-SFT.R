setwd("~/M5542_Bulk_RNA_Seq/GSE130588_WGNCA_05302025")
load("runi_GSE130588_WGNCA_qc_ardigen_analysis_parameters.RData")

# There should be no more necessary changes past this point

### WGCNA STP CALCULATION

# This script uses the datExpr to establish a soft threshold
# USUALLY THIS MUST BE RUN FROM THE UNIX SHELL 
# i.e. $sudo Rscript part2-SFT.R

################################################################################################
### LOAD AND CONFIGURE WGCNA LIBRARY

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(5)

################################################################################################

# Load files
setwd(analysisDir)
load("part1-final-objects.RData")

# Choose a set of soft-thresholding powers
# First generate a sequence of powers to apply
powers <- c(c(1:10), seq(from = 12, to=30, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector = powers, verbose = 5, networkType="signed", corFnc = "bicor",
                         moreNetworkConcepts = TRUE, blockSize=10000) #bicor is the only option

# Plot the results:
sizeGrWindow(12, 6)
par(mfrow = c(1,2));
cex1 <- 0.8; # font size
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab=paste(active.base,"STP Selection",sep=" "),ylab=paste(active.base,"VST Scale Free Topology Model Fit, signed",sep=" "),type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

###SELECT THE THRESHOLD NUMBER AT WHICH THE SFT MODEL FIT PLATEAUS
###NOTE: If graphs are unusual, consider re-filtering, OR if sample numbers are low, use the SFT values on the WGCNA FAQ page table

save(sft, file="part2-output-SFT.RData")

sft.fits <- sft$fitIndices
write.table(sft.fits,file="part2-output-STP.txt", quote=FALSE, sep = "\t",row.names=FALSE, col.names=TRUE)
