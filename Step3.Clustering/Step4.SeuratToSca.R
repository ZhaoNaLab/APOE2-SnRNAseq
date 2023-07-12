library('RColorBrewer')
library('MAST')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

setwd("path/to/project")

# load the seurat object
SeuratObj <- readRDS("Step3.Clustering/SeuratObj.rds")

# loop the transformation 
MajorClusters <- unique(SeuratObj$MajorCluster)
scaObjs <- list()

for (i in 1:length(MajorClusters)){
  SeuratObj <- subset(SeuratObj, MajorCluster == MajorClusters[i])
  
  # filter the seurat object
  DefaultAssay(SeuratObj) <- 'RNA'
  
  ################################# construct singlecellassay object ####################################
  latent.vars <- SeuratObj@meta.data
  latent.vars$wellKey <- rownames(x = latent.vars)
  
  # make sure the expression matrix and the metadata have the same order
  latent.vars <- latent.vars[match(colnames(SeuratObj), latent.vars$wellKey), ]
  
  # prepare the fdat
  fdat <- data.frame(rownames(x = SeuratObj))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  
  # construct the SingleCellAssay object
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(SeuratObj[['RNA']]@counts),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  
  # double check the factor levels
  scaObjs[[i]] <- sca

}

names(scaObjs) <- MajorClusters
saveRDS(scaObjs, "Step3.Clustering/scaObjs.prefilter.rds")



