library('RColorBrewer')
library('MAST')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

setwd("path/to/project")

############ Convert the Seurat object to SCA obejct for Marker calculation #########
# load the seurat object
SeuratObj <- readRDS("Step3.Clustering/SeuratObj.rds")
DefaultAssay(SeuratObj) <- 'RNA'

# Calcualte the percentage of cells expression the features
ExpMat <- GetAssayData(SeuratObj, assay = "RNA", slot = "counts")      
PercMatOverall <- as.matrix(rowMeans(ExpMat > 0))*100
colnames(PercMatOverall) <- "Percentage.Overall"

# Calcualte the percentage of Major clusters expression the features
DefaultAssay(SeuratObj) <- "RNA"
MajorCellType <- unique(SeuratObj$MajorCluster)
FeaturePerc <- list()

for (i in 1:length(MajorCellType)){
  SeuratObj <- subset(SeuratObj, MajorCluster == MajorCellType[i])
  ExpMat <- GetAssayData(SeuratObj, assay = "RNA", slot = "counts")      
  PercMat <- as.matrix(rowMeans(ExpMat > 0))*100
  colnames(PercMat) <- paste0(MajorCellType[i], ".Percentage")
  FeaturePerc[[i]] <- PercMat
}

FeaturePerc.Comb <- do.call(cbind, FeaturePerc)

# all.equal(rownames(PercMatOverall), rownames(FeaturePerc.Comb))
FeatureAllinOne <- cbind(PercMatOverall, FeaturePerc.Comb)

## save the data
saveRDS(FeatureAllinOne, "Step3.Clustering/Feature.Cell.Percent.rds")
