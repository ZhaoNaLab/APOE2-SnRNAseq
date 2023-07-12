# Load the necessary libraries
library("harmony")
library("Seurat")
library("tidyverse")

# Set the working directory
setwd("path/to/project")

# Load the Seurat object
SeuratObj <- readRDS("Step2.QC/SeuratObj.QC.rds")

# Perform various operations on the Seurat object

# Find variable features
SeuratObj <- FindVariableFeatures(SeuratObj, 
                                  selection.method = "vst", 
                                  nfeatures = 3000,
                                  verbose = F)

# Scale the data
SeuratObj <- ScaleData(SeuratObj, 
                       vars.to.regress = c("nFeature_RNA", "percent.mt", "decontX_contamination"), 
                       verbose = FALSE)

# Run PCA
SeuratObj <- RunPCA(SeuratObj, 
                    npcs = 50, 
                    pc.genes = SeuratObj@var.genes,
                    verbose = TRUE)

# Run Harmony for batch effect correction
SeuratObj <- SeuratObj %>% RunHarmony("batch", 
                                      plot_convergence = TRUE,
                                      max.iter.harmony = 50)

# Identify clusters using UMAP, find neighbors, and classify clusters
SeuratObj <- SeuratObj %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 1.0) %>% 
  identity()


##### Cluster cleaning
# visualize the clusters
DimPlot(SeuratObj, label = T) & NoLegend()

## Remove clusters expressing marker genes of different cell type
DefaultAssay(SeuratObj) <- "RNA"
SeuratObj <- NormalizeData(SeuratObj)

# marker gene expression
DotPlot(SeuratObj, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                                      "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                                      "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                                      "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                                      "PDGFRB"), assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 


## Define the major clusters based on markers
SeuratObj$MajorCluster <- NA
SeuratObj$MajorCluster[SeuratObj$seurat_clusters %in% 
                               c(3, 6, 9, 10, 11, 12, 13, 16, 18, 19, 
                                 20, 22, 26, 28, 35, 37, 38, 40, 44, 46, 47, 48, 53)] <- "Ex"

SeuratObj$MajorCluster[SeuratObj$seurat_clusters %in% c(0, 1, 2, 7, 23, 41)] <- "Olig"
SeuratObj$MajorCluster[SeuratObj$seurat_clusters %in% c(4, 21, 49)] <- "Ast"
SeuratObj$MajorCluster[SeuratObj$seurat_clusters %in% c(14, 15, 17, 27, 30, 32, 34, 39)] <- "In"
SeuratObj$MajorCluster[SeuratObj$seurat_clusters %in% c(5)] <- "Mic"
SeuratObj$MajorCluster[SeuratObj$seurat_clusters %in% c(8, 31, 36)] <- "OPC"
SeuratObj$MajorCluster[SeuratObj$seurat_clusters %in% c(33, 45)] <- "Vascular"
SeuratObj$MajorCluster[SeuratObj$seurat_clusters %in% c(24, 25, 29, 42, 43, 50, 51, 52)] <- "Mix"

# remove mixed cell population. 
SeuratObj <- subset(SeuratObj, MajorCluster != "Mix")

## The above clustering and cleaning processes will be done several 
## times until no more mixed cell population can be identified

# Save the Seurat object
saveRDS(SeuratObj, "Step3.Clustering/SeuratObj.rds")

