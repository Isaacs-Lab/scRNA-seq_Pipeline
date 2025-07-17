#----------------------------------
# QC and Merge
### Author: Gary Schweickart
### Created: 7/15/25
### 
#----------------------------------

# Load required packages
library(GenomeInfoDb)
library(Matrix)
library(Seurat)
library(Signac)
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
    is.null(args[4])) {
  stop("Missing at least one argument.", call.=FALSE)
}

input.path = args[1]
out.base = args[2]
qc.csv = read.csv(args[3], row.names = 1)
integrate = as.logical(args[4])

print(args)

# Read in data
sobj.merged <- readRDS(paste0(out.base, "/objects/merged/sobj.merged.metrics.rds"))

cells.init <- dim(sobj.merged@meta.data)[1]

# assign qc cutoffs
sobj.merged <- subset(sobj.merged,
                      subset =
                        nFeature_RNA < qc.csv["nFeature_RNA", "score"] &
                        nCount_RNA < qc.csv["nCount_RNA", "score"] &
                        percent.mt < qc.csv["percent.mt", "score"] &
                        Predicted_Doublets == "Singlet")

cells.post_qc <- dim(sobj.merged@meta.data)[1]
# export qc'd merged object
saveRDS(sobj.merged,paste0(out.base, "/objects/merged/sobj.merged.qc.rds"))
write.csv(data.frame("Timepoint" = c("Initial Cell Count", "Post-QC Cell Count"), "Cells" = c(cells.init, cells.post_qc)), 
          file = paste0(out.base, "/QC/cell_count.csv"))



# Integrate using Harmony
if (integrate) {
# Perform log-normalization and feature selection, as well as SCT normalization on global object
DefaultAssay(sobj.merged) <- "RNA"
sobj.merged <- sobj.merged %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures=3000) %>%
  ScaleData() %>%
  SCTransform(vars.to.regress = NULL)
saveRDS(sobj.merged, file = paste0(out.base, "/objects/merged/sobj.mergedSCT.RDS"))



# Calculate PCs using variable features determined by SCTransform (3000 by default)
sobj.merged <- RunPCA(sobj.merged, assay = "SCT", npcs = 50)
saveRDS(sobj.merged, file = paste0(out.base, "/objects/sobj.mergedPCA.RDS"))


# Integrate with Harmony

sobj.harmony <- RunHarmony(sobj.merged,
                           group.by.vars = c("sample_id"),
                           reduction = "pca", assay.use = "SCT", reduction.save = "rna_harmony")



#sobj.harmony <- readRDS(paste0(out.base, "/objects/sobj.harmony.RDS"))

# Now find clusters using the seurat object normalized with harmony

DefaultAssay(sobj.harmony) <- "SCT"

sobj.harmony <- sobj.harmony  %>%
  RunUMAP(reduction = "rna_harmony", assay = "SCT", dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_') %>%
  FindNeighbors(reduction = "rna_harmony") %>%
  FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), verbose = TRUE)
}
