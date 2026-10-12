setwd("~/M5542_Bulk_RNA_Seq/GSE130588_WGNCA_05302025")
load("runi_GSE130588_WGNCA_qc_ardigen_analysis_parameters.RData")

saveTOMFileBase = paste0(file.base,"_LOG2TPM_2025_05_30") #### change


# You may need to change recut parameters ~ line 70

### PART 4 - Recutting dendrogram trees

###############################################################################################################
### CONCEPTUAL OVERVIEW
# We first recut the blocks of clustered weighted, co-expression gene using all four levels
# of deep split (DS) and generally inspect the per-probe module assignment along with color
# bars of correlation to all clinical traits accompanying the dataset. It is these four DS
# recuts we scrutinize using subsequent part5 figures and calculations.

################################################################################################
### LOAD SUPPORTING AMPEL WGCNA SCRIPTS
# Load all scripts in shared scripts WGCNA directory
# scriptPaths <- list.files(pattern="[.]R$", path="/mnt/disks/rserver/projects/CopyOfshared_scripts/wgcna", full.names=TRUE)
# scriptPaths
# source(scriptPaths)
# rm(scriptPaths)

################################################################################################
### LOAD LIBRARIES AND CONFIGURE SESSION
library(WGCNA)
library(Hmisc) # contains capitalize function
library(Biobase)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

################################################################################################
### LOAD STUDY DATA
setwd(analysisDir)
load("part1-final-objects.RData")
eset <- eset.wgcna # Helps confirm we're loading the eset resulting from WGCNA variance filtering
rm(eset.wgcna)

################################################################################################
### LOAD TOM SETTINGS & MODULE AND COLOR ASSIGNMENTS
load(paste0(saveTOMFileBase,"_TOM_settings.RData"))
load(paste0(saveTOMFileBase,"_NET.RData"))

################################################################################################

projectName <- active.base
save(projectName,file="projectName.RData")
pData <- pData(eset)
fData <- fData(eset)
minModuleSize.original <- as.numeric(TOM.settings["minModuleSize"])
projectDir <- getwd() # Save original working directory to return to after directory traversal
TOMfiles <- net$TOMFiles

################################################################################################
## CALL THE SCRIPT TO ADD COLOR NAMES TO WGNCA MODULES GENERATED DURING TOM CREATION
## Note we'll be regenerating these anyways during the recut, so this is only
## really used for dianostics
source(paste0(scriptDir,"wgcna.net.colornames.R"))

###############################################################################################################
## RECUT THE TOM NETWORK BLOCKS USING FOUR LEVELS OF DEEP SPLIT
## This repeats blockwise module detection from the pre-calculated TOM
## Patience, these recuts take 2-3 minutes each, even on the cloud.

## SET POTENTIAL NEW RECUT PARAMETERS. SEE WGCNA DOCUMENTATION FOR ADDITIONAL INFORMATION
# Here we have the opportunity to experiment specficially usage of the PAM stage, adjustment of
# module cut height, minimum module size, and of course the four levels of deep split (DS1-4).
# Note PAM stage is conducted as a second pass to assign modules not first assigned by the general recut,
# and historically it's always been used.

recut.corType = "pearson" # pearson is safest. bicor is other common option but usually results in many unassigned probes
recut.pamStage <- TRUE # Usually true
recut.pamRespectsDendro <- TRUE # Always true if using PAM stage
recut.detectCutHeight <- 0.99 # Usually 1 for most experiments, adjust carefully!!
recut.minModuleSize <- 25 # Experiment with this carefully! Raise as needed for sets with large numbers of probes/genes
recut.mergeCutHeight=0.2 # Usually 0.2. Module merging can occur downstream

experimentName.recut <- paste("recut.","STP",TOM.settings["power"],".pam",recut.pamStage,".minMod-", recut.minModuleSize,sep="")
experimentName.recut
save(experimentName.recut, file="experimentName.recut.RData")

###############################################################################################################

###############################
### CALL THE RECUT BLOCKS SCRIPT. THIS MAY TAKE UP TO ~10 MINUTES EVEN ON THE CLOUD
source(paste0(scriptDir,"wgcna.blocks.recut.v2.R"))
###############################

# SAVE THE RECUT OBJECTS. IMPORTANT!!! ALL OUTPUT AND FIGURES RELATED TO A RECUT OF THE
# MASTER TOM ARE STORED IN A PROJECT SUBDIRECTORY!!!

setStorageDir <- dget(paste0(scriptDir,"set.storage.directory.R"))
outputStorage <- setStorageDir( experimentName.recut )

save( moduleColors, moduleColors.DS1, moduleColors.DS2, moduleColors.DS3, moduleColors.DS4,
      moduleLabels, moduleLabels.DS1, moduleLabels.DS2, moduleLabels.DS3, moduleLabels.DS4,
      MEs, MEs.labels, MEs.DS1, MEs.DS1.labels, MEs.DS2, MEs.DS2.labels, MEs.DS3, MEs.DS3.labels, MEs.DS4, MEs.DS4.labels,
      colorLevels, colorLevels.DS1, colorLevels.DS2, colorLevels.DS3, colorLevels.DS4,
      file="recut.DS-DS1-DS4.modules.and.genes.RData")

setwd(analysisDir)

###############################################################################################################
### PLOT THE BLOCK DENDROGRAMS OF RECUT PROBES

outputStorage <- setStorageDir( c(experimentName.recut,"plots") )

###      NOTE: This will fail if any of the traits have a single value only!!!!!!

pngWidth=1500 # width of output PNGs
pngHeight=1000 # height of output PNGs
marAll = c(1, 7, 3, 1) # bottom, left, top and right margins
source(paste0(scriptDir,"wgcna.blocks.dendro.with.traitbars.R"))

setwd(analysisDir)

