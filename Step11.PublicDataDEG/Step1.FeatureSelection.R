library("Seurat")
library("biomaRt")

setwd("path/to/project")
ProteinEncodingGenes <- readRDS("Step8.DEG/ProteinEncodingGenes.rds")


# load the data
Seurat.external.final <- readRDS("Step10.PublicDataIntegration/Seurat.external.final.rds")

# Calculate the percentage of cells expressing the genes
DefaultAssay(Seurat.external.final) <- "RNA"

MajorCellType <- unique(Seurat.external.final$MajorCluster)
FeatureSelect <- list()

for (i in 1:length(MajorCellType)){
  SeuratObj <- subset(Seurat.external.final, MajorCluster == MajorCellType[i])
  ExpMat <- GetAssayData(SeuratObj, assay = "RNA", slot = "counts")      
  PercMat <- as.matrix(rowMeans(ExpMat > 0))*100
  colnames(PercMat) <- "Percentage"

  ## include genes expressed by expressed by at least 10% of the cells
  FeaturesIncluded <- rownames(subset(as.data.frame(PercMat), Percentage >= 10))
  
  ## Use only protein-encoding genes
  FeaturesIncluded <- FeaturesIncluded[FeaturesIncluded %in% ProteinEncodingGenes$hgnc_symbol]
  FeatureSelect[[i]] <- FeaturesIncluded
}

names(FeatureSelect) <- MajorCellType

# save the data
saveRDS(FeatureSelect, "Step10.PublicDataIntegration/FeatureSelect.rds")









