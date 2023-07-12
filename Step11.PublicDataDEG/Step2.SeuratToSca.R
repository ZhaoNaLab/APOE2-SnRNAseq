library('RColorBrewer')
library('MAST')
library('tidyverse')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

setwd("path/to/project")

# load the data
Seurat.external.final <- readRDS("Step10.PublicDataIntegration/Seurat.external.final.rds")

## Update cngeneson and zAge
Seurat.external.final$cngeneson <- scale(Seurat.external.final$nFeature_RNA)

Meta.external.ID <- Seurat.external.final@meta.data %>%
  dplyr::select(orig.ident, Age) %>%
  remove_rownames() %>%
  distinct()
  
Meta.external.ID$zAge <- scale(Meta.external.ID$Age)

# Add the zAge back to Meta.external.Cell
Meta.external.Cell <- Seurat.external.final@meta.data
Meta.external.Cell$zAge <- NULL
Meta.external.Cell <- Meta.external.Cell %>% 
  rownames_to_column("Cell") %>%
  left_join(Meta.external.ID %>%
              dplyr::select(orig.ident, zAge), 
            by = "orig.ident") %>%
  column_to_rownames("Cell")


# update the metadata
Seurat.external.final <- AddMetaData(Seurat.external.final, Meta.external.Cell)

# set the levels for factors: Pathological Dx
Seurat.external.final$Dx <- factor(Seurat.external.final$Dx, levels = c("Ctrl", "AD"))

# set the levels for factors: sex
Seurat.external.final$Sex <- factor(Seurat.external.final$Sex, levels = c("M", "F"))

# set the levels for factors: apoe
Seurat.external.final$Apoe[Seurat.external.final$Apoe == 23] <- "E2"
Seurat.external.final$Apoe[Seurat.external.final$Apoe == 33] <- "E3"
Seurat.external.final$Apoe[Seurat.external.final$Apoe == 34] <- "E4"
Seurat.external.final$Apoe[Seurat.external.final$Apoe == 44] <- "E4"
Seurat.external.final$Apoe <- factor(Seurat.external.final$Apoe, levels = c("E2", "E3", "E4"))

# load the gene filters
GeneSelect <- readRDS("Step11.PublicDataDEG/FeatureSelect.rds")
  
# loop the transformation 
MajorClusters <- unique(Seurat.external.final$MajorCluster)
scaObjs <- list()

for (i in 1:length(MajorClusters)){
  SeuratObj <- subset(Seurat.external.final, MajorCluster == MajorClusters[i])
  
  # filter the seurat object
  DefaultAssay(SeuratObj) <- 'RNA'
  SeuratObj <- subset(SeuratObj, features = GeneSelect[[MajorClusters[i]]])
  
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
saveRDS(scaObjs, "Step11.PublicDataDEG/scaObjs.rds")



