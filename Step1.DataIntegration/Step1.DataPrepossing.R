# Import the batchtools package
library("batchtools")

# Create a list of directory paths for data import
# There are 56 folders, each containing a dataset, named GEX_1R, GEX_2R, ..., GEX_56R
DirList <- list()
for (i in 1:56){
  DirList[[i]] <- paste0(".../GEX_", i, "R/")
}
DirList <- unlist(DirList)

# Change the working directory to snRNAseq/Step1.DataIntegration
setwd(".../Step1.DataIntegration")

# Initialize an experiment registry that allows all analysis steps to be logged and reproduced
# Load necessary packages and set the seed for reproducibility
reg <- makeExperimentRegistry("SnRNAseq.Registry",
                              packages = c("tidyverse", "Seurat", "data.table",
                                           "singleCellTK", "purrr", "DoubletFinder",
                                           "celda", "ggpubr", "biomaRt"), 
                              seed = 20221005)

# Load the registry, which allows to store data, definitions, and results in an organized way
reg <- loadRegistry("SnRNAseq.Registry", writeable = TRUE)

# The function subsample processes each sample from the list of directories
# It reads data, performs quality control, and prepares the data for downstream analyses
subsample <- function(data, job, i) {
  sampleID <- i
  sce <- importCellRanger(sampleDirs = data[i])
  sce.raw <- importCellRanger(sampleDirs = data[i], dataType = "raw")
  sce <- decontX(sce, background = sce.raw)
  rownames(sce) <- rowData(sce)$feature_name
  Seurat.obj <- as.Seurat(sce, counts = "counts", data = NULL)
  
  # add another assay to the Seurat object
  Seurat.obj[["decontX"]] <- CreateAssayObject(counts = assay(sce, "decontXcounts"))
  
  # rename the assays
  Seurat.obj <- RenameAssays(object = Seurat.obj, originalexp = 'RNA')
  Seurat.obj@meta.data <- Seurat.obj@meta.data[, c(1, 10, 9, 8, 7)]
  list(data = data, Seurat.obj = Seurat.obj, sampleID = sampleID)  # obtain different parameters that can be parsed to the next step
}

# Add problem definition for the data processing step
addProblem(name = "snRNAseqDat", 
           data = DirList, 
           fun = subsample) 

# Function snRNAseqClean performs further preprocessing and cleaning on each processed sample
snRNAseqClean <- function(data, job, instance) {
  data = instance$Seurat.obj
  sampleID = instance$sampleID
  
  ####### remove cells with features <200 and features expressed by less than 5 cells #######
  SelectedCells <- WhichCells(data, expression = nFeature_decontX > 200)
  SelectedFeatures <- rownames(data)[Matrix::rowSums(data) > 5]
  data <- subset(data, features = SelectedFeatures, cells = SelectedCells)
  
  # remove cell with a mitochondria percentage >= 10%
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data[["percent.rb"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
  SelectedMito <- WhichCells(data, expression = percent.mt < 10)
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
  colnames(data@meta.data)[grepl("pANN", colnames(data@meta.data))] <- "pANN"
  
  # Save the individual Seurat file
  SaveDir <- ".../Step1.DataIntegration/SeuratObj.Ind"
  saveRDS(data, paste0(SaveDir, "/Seurat.obj.S", sampleID, ".rds"))
}

# Add the cleaning algorithm to the registry
addAlgorithm(name = "snRNAseqClean", fun = snRNAseqClean)

# Define parameter settings for the data processing step
prob.design <- CJ(i = 1:56)

# Final problem design
pdes <- list(snRNAseqDat = prob.design)

# Algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(snRNAseqClean = data.table())

# Add experiments with defined problem and algorithm designs to the registry
addExperiments(pdes, ades)

# Summarize the experiment
summarizeExperiments()

# Submit jobs without showing the progress
options(batchtools.progress = FALSE)
submitJobs() 

# You can uncomment the following lines to check the status and results of the jobs
# getStatus()
# getJobTable()
# getErrorMessages()
# findExpired() 
