################################################################################
## GEX Functions
## Desc: Functions used to create the sc/snRNA-seq objects.
## Last Updated: 7/11/25
################################################################################

# Create seurat objects function
create_sobj <- function(name, input.path, out.base, annotations, input.format = "cellranger", project = "Seurat_Project") {
  
  # create paths
  if (!endsWith(out.base, "/")) {
    out.base <- paste0(out.base, "/")
  }
  # create paths
  if (!endsWith(input.path, "/")) {
    input.path <- paste0(input.path, "/")
  }
  
  # define paths
  
  if (input.format == "cellranger") {
    dir <- paste0(input.path, name, "/outs")
  } else {
    dir <- paste0(input.path, name)
  }
  
  output.path <- paste0(out.base, "individual/", name)
  mat.path <- paste(dir, "filtered_feature_bc_matrix", sep = "/")
  
  # read data
  rna.data <- Read10X(mat.path)

  
  # Create Seurat object
  sobj <-   CreateSeuratObject(counts = rna.data, project = project, min.cells = 3, min.features = 200)
  sobj@meta.data$sample_id <- name
  
  # Save RDS of seurat object
  output <- paste(output.path, ".rds", sep = "")
  saveRDS(object= sobj, file = output)
  
  message(paste0(name, " RDS saved."))
  return(sobj)
} # end create_sobj

#########################################################################################################
# cellCycleScoring

cellCycleScoring <- function(sobj) {
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  DefaultAssay(sobj) <- "RNA"
  sobj <- NormalizeData(sobj)
  sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
  sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
  return(sobj)
} # end cellCycleScoring

#########################################################################################################
# sobj_qc

plotQC <-function(sobj, out ="./", 
                  sample=NULL, sample.column = NULL, 
                  find.markers= TRUE,merged = FALSE,
                  input.format = "cellranger") {
  
  # need to have column "sample_id"
  # Takes in list of samples and generate gene expression quality control plots and metrics for each sample.
  if (!merged) {
    if (is.null(sample) & !is.null(sample.column)) {
      sample <- sobj@meta.data[,sample.column][1]
    } else if (is.null(sample) & is.null(sample.column)) {
      stop("Need sample or sample.column arguments.")
    }
    
    if (run.macs2 & is.null(macs2.path)) {
      stop("When running MACS2, the 'macs2.path' argument is required.")
    }
    
    # create paths
    if (!endsWith(out, "/")) {
      out <- paste0(out, "/")
    }
    
    out.path <- paste0(out, sample)
    dir.create(out.path)
    
    message(paste0("Running ", sample, "."))
    
    # QC for RNA
    
    DefaultAssay(sobj) <- "RNA"
    
    # store mitochondrial percentage in object meta data
    sobj <- PercentageFeatureSet(sobj, pattern = "MT-", col.name = "percent.mt")
    
    
   
    # Create Violin plot
    metrics <- c("nFeature_RNA","nCount_RNA","percent.mt")
    VlnPlot(object = sobj,
            features = metrics,
            pt.size = 0.1,
            ncol = 3,
    )
    ggsave(paste0(out.path,"/violinPlot.pdf"), units = "in",height = 4, width = 8)
    
    ##############################################################################
    # pre-processing and dimensional reduction on the RNA assay
    # using standard approach for RNA-seq data.
    
    # RNA analysis
    DefaultAssay(sobj) <- "RNA"
    sobj <- sobj %>%
      SCTransform(vst.flavor = 'v2') %>% 
      RunPCA(reduction.name = "pca", reduction.key = "PCArna_") %>%
      RunUMAP(reduction="pca", reduction.name = "umap.rna", reduction.key = "UMAPrna_", dims = 1:30) %>%
      FindNeighbors(dims = 1:30)
    
    # test cluster resolution
    sobj <- FindClusters(sobj, resolution=0.5)
    DimPlot(sobj, label = TRUE)
    ggsave((paste0(out.path,"/UMAP_RNA.pdf")))
    
    # Find all markers
    
    if (find.markers) {
      DefaultAssay(sobj) <- "RNA"
      markers.rna <- FindAllMarkers(sobj, only.pos = TRUE)
      write.csv(markers.rna,paste0(out.path, "/markers.rna.csv"))
    }
    
    # Find variable features
    DefaultAssay(sobj) <- "RNA"
    sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
    # plot variable features 
    LabelPoints(VariableFeaturePlot(sobj), points = head(VariableFeatures(sobj), 10), repel = TRUE)
    ggsave((paste0(out.path,"/variableFeaturesRNA.pdf")),height = 8.5, width = 11)
    # Feature plot with QC metrics
    
    sobj_copy <- sobj
    sobj_copy$nucleosome_signal[is.infinite(sobj_copy$nucleosome_signal)] <- 0
    FeaturePlot(sobj_copy, 
                features = metrics,
                pt.size = 0.4,
                label = TRUE
    )
    ggsave((paste0(out.path,"/FeaturePlot.pdf")), height = 8.5, width = 11)
    rm(sobj_copy)
    
    #top10 genes
    if(find.markers) {
      DefaultAssay(sobj) <- "SCT"
      markers.sobj <- markers.rna
      markers.sobj %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1) %>%
        slice_head(n = 10) %>%
        ungroup() -> top10
      DoHeatmap(sobj, features = top10$gene) + NoLegend()
      ggsave((paste0(out.path,"/topHeatmap.pdf")),width = 28,height = 14)
      
      
      markers.rna %>%
        group_by(cluster) %>%
        slice_head(n = 5) %>%
        ungroup() -> top
      rna.list <- lapply(levels(sobj), function(cluster) {
        cluster.markers <- markers.rna[which(markers.rna$cluster==cluster),"gene"]
      })
      names(rna.list) <- paste0("cluster", 0:(length(unique(sobj$seurat_clusters)) - 1))
      saveRDS(rna.list,paste0(out.path, "/rna.list.rds"))
    }
    
   
    # plot both umaps
    umap.rna.plot <- DimPlot(sobj, reduction="umap.rna", label=T) + NoLegend() 
    umap.atac.plot <- DimPlot(sobj, reduction="umap.atac", label=T) + NoLegend() 
    ggsave(paste0(out.path, "/umapRNA_ATAC.pdf"), 
           plot = umap.rna.plot + umap.atac.plot,
           units = "in",
           height = 8.5,
           width = 11)
  } else {
    # create paths
    if (!endsWith(out, "/")) {
      out <- paste0(out, "/")
    }
    
    out.path <- out
    dir.create(out.path)
    # Perform log-normalization and feature selection, as well as SCT normalization on global object
    DefaultAssay(sobj) <- "RNA"
    sobj <- sobj %>%
      NormalizeData() %>%
      FindVariableFeatures(selection.method = "vst", nfeatures=3000) %>%
      ScaleData() %>%
      SCTransform(vars.to.regress = NULL)
    
    # Calculate PCs using variable features determined by SCTransform (3000 by default)
    sobj <- RunPCA(sobj, assay = "SCT", npcs = 50)

    sobj <- sobj  %>%
      RunUMAP(assay = "SCT", dims = 1:50) %>%
      FindNeighbors() %>%
      FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), verbose = TRUE)
    # QC for snATAC-seq
    # store mitochondrial percentage in object meta data
    DefaultAssay(sobj) <- "RNA"
    sobj <- PercentageFeatureSet(sobj, pattern = "MT-", col.name = "percent.mt")
    
    
    # Feature plot with QC metrics
    metrics <- c("nFeature_RNA","nCount_RNA","percent.mt", "CC.Difference", "S.Score", "G2M.Score")
    
    FeaturePlot(sobj, 
                features = metrics,
                pt.size = 0.4,
                label = TRUE
    )
    ggsave(paste0(out.path,"/FeaturePlot.pdf"),
           units = "in",
           height = 8.5,
           width = 11)
    DimPlot(sobj, group.by ="SCT_snn_res.0.6",label = T) + ggtitle("res0.6 umap")
    ggsave(paste0(out.path,"/umap_merge.pdf"),
           units = "in",
           height = 8.5,
           width = 11)
    DimPlot(sobj, group.by ="SCT_snn_res.0.6",label = T) + 
      ggtitle("res0.6 UMAP and Feature plot of QC metrics") + 
      theme(legend.position = "none")+ 
      FeaturePlot(sobj, 
                  features = metrics,
                  pt.size = 0.4,
                  label = F)
    ggsave(paste0(out.path,"/metrics_features_merge.pdf"), 
           units = "in",
           height = 12,
           width = 20)
    
    VlnPlot(
      object = sobj,
      features = c("nFeature_RNA","nCount_RNA","percent.mt", "CC.Difference", "S.Score", "G2M.Score"),
      pt.size = 0.1,
      ncol = 3,
    )
    
    ggsave(paste0(out.path,"/violin_plot_merge.pdf"),
           units = "in",
           height = 10,
           width = 12
    )
    
    DefaultAssay(sobj) <- "RNA"
    sobj@meta.data$sample_trunc <- sub("_MW_.*", "", sobj@meta.data$sample)
    ggsave(paste0(out.path,"/violin_plot0.6_merge.pdf"), 
           plot = VlnPlot(
             object = sobj,
             features = c("TSS.enrichment","nucleosome_signal","blacklist_ratio",
                          "nFeature_ATAC","nCount_ATAC","percent.mt", "CC.Difference", "S.Score", "G2M.Score"),
             group.by ="sample_trunc",
             pt.size = 0.1,
             ncol = 3
           ),
           units = "in",
           height = 10,
           width = 12
    )
    sobj@meta.data$sample_trunc <- NULL
    
  }
  
  return(sobj)
}




