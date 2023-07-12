library('RColorBrewer')
library('MAST')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

# load the seurat object
ProteinEncodingGenes <- readRDS(".../Step4.DEG/ProteinEncodingGenes.rds")

InputDat <- c("Ex.Harmony.rds", "In.Harmony.rds", "Ast.Harmony.rds",
              "Mic.Harmony.rds", "Olig.Harmony.rds", "OPC.Harmony.rds")

MajorCluster <- c("Ex", "In", "Ast", "Mic", "Olig", "OPC")

## loop the creation of SCA files
Sca.subcluster.list <- list()

for (i in 1:length(MajorCluster)){
  Seurat.obj <- readRDS(paste0(".../Step6.SubclusterAnalysis/", MajorCluster[i], "/", InputDat[i]))
  DefaultAssay(Seurat.obj) <- "RNA"
  
  Seurat.obj@meta.data$SubCluster <- paste0(Seurat.obj@meta.data$MajorCluster, ".", Seurat.obj@meta.data$seurat_clusters)
  SubClusters <- unique(Seurat.obj@meta.data$SubCluster)

  for(j in 1:length(SubClusters)){
    Sub.Seurat.obj <- subset(Seurat.obj, SubCluster == SubClusters[j])
    ExpMat <- GetAssayData(Sub.Seurat.obj, assay = "RNA", slot = "counts")      
    PercMat <- as.matrix(rowMeans(ExpMat > 0))*100
    colnames(PercMat) <- "Percentage"
    
    ## include genes expressed by expressed by at least 15% of the cells
    FeaturesIncluded <- rownames(subset(as.data.frame(PercMat), Percentage >= 10))
    
    ## Use only protein-encoding genes
    FeaturesIncluded <- FeaturesIncluded[FeaturesIncluded %in% ProteinEncodingGenes$hgnc_symbol]
    # filter the seurat object
    DefaultAssay(Sub.Seurat.obj) <- 'RNA'
    
    SeuratObj.final <- subset(Sub.Seurat.obj, features = FeaturesIncluded)
    ################################# construct singlecellassay object ####################################
    latent.vars <- SeuratObj.final@meta.data
    latent.vars$wellKey <- rownames(x = latent.vars)
    
    # make sure the expression matrix and the metadata have the same order
    latent.vars <- latent.vars[match(colnames(SeuratObj.final), latent.vars$wellKey), ]
    
    # prepare the fdat
    fdat <- data.frame(rownames(x = SeuratObj.final))
    colnames(x = fdat)[1] <- "primerid"
    rownames(x = fdat) <- fdat[, 1]
    
    # construct the SingleCellAssay object
    sca <- MAST::FromMatrix(
      exprsArray = as.matrix(SeuratObj.final[['RNA']]@counts),
      check_sanity = FALSE,
      cData = latent.vars,
      fData = fdat)
    Sca.subcluster.list[[SubClusters[j]]] <- sca
  }
  
  setwd("...q/Step7.SubclusterDEG")
  saveRDS(Sca.subcluster.list, "sca.SubClusters.Objs.rds")
}


Sca.subcluster.list <- readRDS(".../Step7.SubclusterDEG/sca.SubClusters.Objs.rds")
