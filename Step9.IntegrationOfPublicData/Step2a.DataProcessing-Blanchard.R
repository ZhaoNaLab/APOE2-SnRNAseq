library("batchtools")


## set the working directory
setwd(".../Step9.IntegrationOfPublicData/Blanchard")
## use the bath tools for RNA contamination detection and doublet removal
reg <- makeExperimentRegistry("SnRNAseq.Registry",
                              packages = c("tidyverse", "Seurat", "data.table",
                                           "singleCellTK", "purrr", "DoubletFinder",
                                           "celda", "ggpubr", "biomaRt"), seed = 20221005)

reg <- loadRegistry("SnRNAseq.Registry", writeable = TRUE)

# load the combined SCE objects
setwd(".../Step9.IntegrationOfPublicData/Blanchard")
sce.all <- readRDS("sce.all.rds")

## Define instance
# define instance (data preparation)
subsample <- function(data, job, i) {
  sampleID <- unique(colData(data)$individualID)[i]
  sce <- subset(data, , individualID == sampleID)
  sce <- decontX(sce)
  Seurat.obj <- as.Seurat(sce, counts = "counts", data = NULL)
  
  # add another assay to the Seurat object
  Seurat.obj[["decontX"]] <- CreateAssayObject(counts = assay(sce, "decontXcounts"))
  
  # rename the assays
  Seurat.obj <- RenameAssays(object = Seurat.obj, originalexp = 'RNA')
  list(data = data, Seurat.obj = Seurat.obj, sampleID = sampleID)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "snRNAseqDat", 
           data = sce.all, 
           fun = subsample) 


# Define tasks
snRNAseqClean <- function(data, job, instance) { # instance is inherited from the last step
  data = instance$Seurat.obj
  sampleID = instance$sampleID
  
  ####### remove cells with features <200 and features expressed by less than 5 cells #######
  SelectedCells <- WhichCells(data, expression = nFeature_decontX > 200)
  SelectedFeatures <- rownames(data)[Matrix::rowSums(data) > 5]
  data <- subset(data, features = SelectedFeatures, cells = SelectedCells)
  
  # remove cell with a mitochondria percentage >= 10%
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data[["percent.rb"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
  SelectedMito <- WhichCells(data, expression = percent.mt < 5)
  data <- subset(data, cells = SelectedMito)
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data)
  data <- ScaleData(data, vars.to.regress = c("nFeature_decontX", "percent.mt"), verbose = FALSE)
  data <- RunPCA(data, verbose = F, npcs = 20)
  data <- RunUMAP(data, dims = 1:20, verbose = F)
  
  # run parameter optimization with paramSweep
  sweep.res <- paramSweep_v3(data) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # define the expected number of doublet cellscells.
  nExp <- round(ncol(data) * 0.04)  # expect 4% doublets
  data <- doubletFinder_v3(data, pN = 0.25, pK = mpK, nExp = nExp, PCs = 1:15)
  colnames(data@meta.data)[grepl("DF.classification", colnames(data@meta.data))] <- "DFClass"
  colnames(data@meta.data)[grepl("pANN", colnames(data@meta.data))] <- "pANN"
  
  # Save the individual Seurat file
  SaveDir <- ".../Step9.IntegrationOfPublicData/Blanchard/SeuratObj.Ind"
  saveRDS(data, paste0(SaveDir, "/Seurat.obj.", sampleID, ".rds"))
}


addAlgorithm(name = "snRNAseqClean", fun = snRNAseqClean)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design <- CJ(i = 1:32)

# final problem design
pdes <- list(snRNAseqDat = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(snRNAseqClean = data.table())

addExperiments(pdes, ades) # pdes set parameter to the data part; ades set parameters to the instance parts

# summarize the experiment
summarizeExperiments()

# Submitting and Collecting Results
options(batchtools.progress = FALSE)
submitJobs()
# getStatus()
# getJobTable()
# getErrorMessages()
# findExpired() 
