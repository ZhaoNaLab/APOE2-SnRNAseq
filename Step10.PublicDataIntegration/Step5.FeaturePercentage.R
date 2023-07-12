library('RColorBrewer')
library('MAST')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

setwd(".../Step10.PublicDataIntegration")

# load the seurat object
Seurat.external.final <- readRDS("Seurat.external.final.rds")
DefaultAssay(Seurat.external.final) <- 'RNA'

# Calcualte the percentage of cells expression the features
ExpMat <- GetAssayData(Seurat.external.final, assay = "RNA", slot = "counts")      
PercMatOverall <- as.matrix(rowMeans(ExpMat > 0))*100
colnames(PercMatOverall) <- "Percentage.Overall"

# Calcualte the percentage of Major clusters expression the features
DefaultAssay(Seurat.external.final) <- "RNA"
MajorCellType <- unique(Seurat.external.final$MajorCluster)
FeaturePerc <- list()

for (i in 1:length(MajorCellType)){
  SeuratObj <- subset(Seurat.external.final, MajorCluster == MajorCellType[i])
  ExpMat <- GetAssayData(SeuratObj, assay = "RNA", slot = "counts")      
  PercMat <- as.matrix(rowMeans(ExpMat > 0))*100
  colnames(PercMat) <- paste0(MajorCellType[i], ".Percentage")
  FeaturePerc[[i]] <- PercMat
}

FeaturePerc.Comb <- do.call(cbind, FeaturePerc)

# all.equal(rownames(PercMatOverall), rownames(FeaturePerc.Comb))
FeatureAllinOne <- cbind(PercMatOverall, FeaturePerc.Comb)

## save the data
saveRDS(FeatureAllinOne, "Feature.Cell.Percent.rds")



