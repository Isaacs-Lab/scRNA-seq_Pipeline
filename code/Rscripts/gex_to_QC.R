#----------------------------------
# Multiome pipeline to QC
### Author: Gary Schweickart
### Created: 9/17/24
### Updated: 2/11/25
#----------------------------------

# Load required packages
library(GenomeInfoDb)
library(Matrix)
library(Seurat)
library(patchwork)
library(tidyverse)
library(dplyr)
library(readr)
library(magrittr)
library(readxl)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020) 
library(TFBSTools)


# Setup
args=commandArgs(trailingOnly = TRUE)
options(future.globals.maxSize= 1073741824000)
if (length(args)==0) {
  stop("All arguments need to be applied.n", call.=FALSE)
}

if (is.null(args[1]) | is.null(args[2]) | is.null(args[3]) | 
    is.null(args[4]) | is.null(args[5])) {
  stop("Missing at least one argument (besides seed).", call.=FALSE)
}

input.path = args[1]
out.base = args[2]
sample.csv = read.csv(args[3])
scrub = as.logical(args[4])
#seed = as.numeric(args[6])
project = args[5]

sample.names = sample.csv$sample_names
print(args)

#if (is.na(seed)) {
#  seed = 6  
#}
#set.seed(seed)
set.seed(6)

dir.create(paste0(out.base,"/objects"))
dir.create(paste0(out.base,"/objects/lists"))
dir.create(paste0(out.base,"/objects/individual"))

dir.create(paste0(out.base, "/QC"))
dir.create(paste0(out.base, "/QC/merged"))



## Load multiome functions
source("code/Rscripts/multiomeFunctions.R")

## Object creation
# Run create_sobj function
sobj.list <- lapply(sample.names, create_sobj, input.path, out.base=paste0(out.base,"/objects"), annotations=annotations, project = project)
names(sobj.list) <- sample.names
saveRDS(sobj.list, file = paste0(out.base, "/objects/lists/sobj.list.raw.rds"))
print("All saved.")



# cell cycle scoring

sobj.list <- lapply(sobj.list, cellCycleScoring)
names(sobj.list) <- sample.names
saveRDS(sobj.list, paste0(out.base, "/objects/lists/sobj.list.ccs.rds"))

# run qc function
sobj.list <- lapply(sobj.list, plotQC, frag.base.path = input.path, frag.full.path = NULL, out =paste0(out.base, "/QC"), 
                    sample=NULL, sample.column = "sample_id", find.markers= TRUE, 
                    run.macs2 = FALSE, macs2.path = NULL, group.by.macs2 = "seurat_clusters")

# save list
saveRDS(sobj.list, paste0(out.base, "/objects/lists/sobj.list.metrics.rds"))


# merge
dir.create(paste0(out.base, "/objects/merged"))
sobj.merged <- merge(sobj.list[[1]], sobj.list[2:length(sobj.list)])
#sobj.merged <- JoinLayers(sobj.merged)
saveRDS(sobj.merged, paste0(out.base,"/objects/merged/sobj.merged.rds"))

# merged QC

sobj.merged <- plotQC(sobj = sobj.merged, out = paste0(out.base, "/QC/merged"), merged = TRUE, sample = NULL, sample.column = NULL)

# add scrublet

if (scrub) {
  # Get doublet scores
  sobj.merged@meta.data$doublet_score <- "0"
  for (sample in sample.names) {
    pscrub <- read_table(paste0(out.base,"/scrub/", sample, "_doublet_scores.txt"),col_names = F)
    
    sobj.merged@meta.data[sobj.merged@meta.data$sample_id == sample, "doublet_score"] <- pscrub$X1
  }
  sobj.merged@meta.data$Predicted_Doublets <- ifelse(sobj.merged$doublet_score > 0.25, "Doublet","Singlet" )
  table(sobj.merged@meta.data$Predicted_Doublets)
  
}


saveRDS(sobj.merged, paste0(out.base, "/objects/merged/sobj.merged.metrics.rds"))

### END gex_to_qc
