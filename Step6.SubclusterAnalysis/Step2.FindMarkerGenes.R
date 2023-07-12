library("batchtools")

## use the bath tools for RNA contamination detection and doublet removal
setwd(".../Step6.SubclusterAnalysis")
reg <- makeExperimentRegistry("Reclustering.Marker.Registry",
                              conf.file = ".../batchtools.conf.R", 
                              packages = c("tidyverse", "Seurat", "data.table", "harmony",
                                           "MAST", "lme4", "biomaRt", "RColorBrewer"), 
                              seed = 20221005)

reg <- loadRegistry("Reclustering.Marker.Registry", writeable = TRUE)


InputDat <- c("Ex.Harmony.rds", "In.Harmony.rds", "Ast.Harmony.rds", "Mic.Harmony.rds", "Olig.Harmony.rds", "OPC.Harmony.rds")

## Define instance
# define instance (data preparation)
subsample <- function(data, job, i) {
  MajorCluster <- c("Ex", "In", "Ast", "Mic", "Olig", "OPC")
  MajorCluster <- MajorCluster[i]
  
  # read the data
  data <- readRDS(paste0(".../Step6.SubclusterAnalysis/", MajorCluster, "/", data[i]))
  DefaultAssay(data) <- "RNA"
  list(data = data, MajorCluster = MajorCluster)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "MarkderFinder", 
           data = InputDat,
           fun = subsample) 


# Define tasks
FindMarkers <- function(data, job, instance) { # instance is inherited from the last step
  data = instance$data
  MajorCluster = instance$MajorCluster
  
  ## step1: Calculate the percentage of cells expressing the genes
  # Calcualte the percentage of cells expression the features
  ExpMat <- GetAssayData(data, assay = "RNA", slot = "counts")      
  PercMatOverall <- as.matrix(rowMeans(ExpMat > 0))*100
  colnames(PercMatOverall) <- "Percentage.Overall"
  
  # Calcualte the percentage of Major clusters expression the features
  DefaultAssay(data) <- "RNA"
  CellType <- unique(as.numeric(as.character(data$seurat_clusters)))
  FeaturePerc <- list()
  
  for (j in 1:length(CellType)){
    SeuratObj <- subset(data, seurat_clusters == CellType[j])
    ExpMat <- GetAssayData(SeuratObj, assay = "RNA", slot = "counts")      
    PercMat <- as.matrix(rowMeans(ExpMat > 0))*100
    colnames(PercMat) <- paste0("Cluster.", CellType[j], ".Percentage")
    FeaturePerc[[j]] <- PercMat
  }
  
  FeaturePerc.Comb <- do.call(cbind, FeaturePerc)
  
  # all.equal(rownames(PercMatOverall), rownames(FeaturePerc.Comb))
  FeatureAllinOne <- cbind(PercMatOverall, FeaturePerc.Comb)
  
  
  ########### construct singlecellassay object #################
  latent.vars <- data@meta.data
  latent.vars$wellKey <- rownames(x = latent.vars)
  
  # make sure the expression matrix and the metadata have the same order
  latent.vars <- latent.vars[match(colnames(data), latent.vars$wellKey), ]
  
  # prepare the fdat
  fdat <- data.frame(rownames(x = data))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  
  # construct the SingleCellAssay object
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(data[['RNA']]@counts),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  
  ####### FindMarkers #######
  for (k in 1:length(CellType)){
    tryCatch({
      Cluster <- CellType[k]
      Pos <- grep(paste0(Cluster), colnames(FeatureAllinOne))
      
      ## Select features expressed by at least 25% of the target cells
      FeatureSelected <- rownames(FeatureAllinOne)[FeatureAllinOne[, Pos] >= 25]
      sub.sca <- sca[FeatureSelected, ]
      
      ## Make a new varible CulsterOfInterest for Marker calculation
      colData(sub.sca)$CulsterOfInterest <- "Ref"
      colData(sub.sca)$CulsterOfInterest[colData(sub.sca)$seurat_clusters == Cluster] <- "Target"
      colData(sub.sca)$CulsterOfInterest <- factor(colData(sub.sca)$CulsterOfInterest,
                                                   level = c("Ref", "Target"))
      
      # MAST model for marker identification 
      zlmCond <- zlm(~ CulsterOfInterest + cngeneson + percent.mt + decontX_contamination,
                     parallel = 4,
                     sca = sub.sca)
      
      summaryCond <- summary(zlmCond, doLRT= c('CulsterOfInterestTarget'))
      # print out the data.table
      summaryDt <- summaryCond$datatable
      Markers <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                       summaryDt[component=='logFC', .(primerid, coef, contrast, z)], 
                       by=c('primerid', 'contrast')) #logFC coefficients
      
      Markers$contrast <- paste0(MajorCluster, ".", Cluster)
      
      # add expression level information
      Markers <- Markers %>% 
        left_join(FeatureAllinOne %>% as.data.frame() %>%
                    rownames_to_column("primerid"), by = "primerid")
      
      # save the result
      setwd("/research/labs/moleneurosci/bug/m198507/Projects/snRNAseq/Step12.SubclusterAnalysis")
      write.csv(Markers, paste0(MajorCluster, "/CellMarkers.", MajorCluster, ".", Cluster, ".csv"))
    }, error = function(x){})
  }
}

addAlgorithm(name = "FindMarkers", fun = FindMarkers)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design <- CJ(i = 1:6)

# final problem design
pdes <- list(MarkderFinder = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(FindMarkers = CJ())

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

