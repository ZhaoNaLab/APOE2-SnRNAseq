library('RColorBrewer')
library('MAST')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

setwd(".../Step10.PublicDataIntegration")

# load the seurat object
SeuratObj <- readRDS("Seurat.external.final.rds")

# set the levels for factors: Pathological Dx
SeuratObj$Dx <- factor(SeuratObj$Dx, levels = c("Ctrl", "AD"))

# set the levels for factors: Sex
SeuratObj$Sex <- factor(SeuratObj$Sex, levels = c("male", "female"))

# set the levels for factors: apoe
SeuratObj@meta.data$apoe <- NULL
SeuratObj@meta.data$apoe[SeuratObj@meta.data$Apoe %in% c(23)] <- "E2"
SeuratObj@meta.data$apoe[SeuratObj@meta.data$Apoe %in% c(33)] <- "E3"
SeuratObj@meta.data$apoe[SeuratObj@meta.data$Apoe %in% c(34, 44)] <- "E4"
SeuratObj@meta.data$apoe <- factor(SeuratObj@meta.data$apoe, levels = c("E2", "E3", "E4"))

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

saveRDS(sca, "ScaObjsInOne.rds")



