library("batchtools")

"""

We compared how Theta and the number of PCs can affect the final clustering results.
The optimal parameters were manually selected based on how well the subclusters can 
be seperated 


"""

setwd(".../Step12.SubclusterAnalysis")
reg <- makeExperimentRegistry("Reclustering.Registry",
                              conf.file = ".../Step6.SubclusterAnalysi/batchtools.conf.R",
                              packages = c("tidyverse", "Seurat", "data.table", "harmony"),
                              seed = 20221005)

reg <- loadRegistry("Reclustering.Registry", writeable = TRUE)

Dir <- ".../Step3.Clustering/SeuratObj.rds"

## Define instance
# define instance (data preparation)
subsample <- function(data, job, i) {
  MajorCluster <- c("Ex", "In", "Ast", "Olig", "OPC", "Mic", "Vascu")
  Cluster <- MajorCluster[i]
  
  # load the data
  data <- readRDS(data)
  data <- subset(data, MajorCluster == Cluster)
  
  data <- DietSeurat(data, assays = c("RNA", "decontX"), counts = TRUE)
  DefaultAssay(data) <- "RNA"
  list(data = data, Cluster = Cluster)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "ReclusterMajorCluster", 
           data = Dir, #Seurat.ZH.final
           fun = subsample) 

# Define tasks
HarmoneyRecluster <- function(data, job, instance, Theta, NPCs) { # instance is inherited from the last step
  data = instance$data
  Cluster = instance$Cluster
  
  ####### Harmony reclustering #######
  options(repr.plot.height = 2.5, repr.plot.width = 6)
  data <- FindVariableFeatures(data, 
                               selection.method = "vst", 
                               nfeatures = 3000,
                               verbose = F)
  data <- ScaleData(data, 
                    vars.to.regress = c("nFeature_RNA", "percent.mt", "decontX_contamination"), 
                    verbose = FALSE)
  
  data <- RunPCA(data, 
                 npcs = NPCs, 
                 pc.genes = data@var.genes,
                 verbose = TRUE)
  
  
  ## integration based on individual
  data <- data %>% RunHarmony("orig.ident", 
                              theta = Theta,
                              dims.use = 1:NPCs,
                              plot_convergence = TRUE,
                              max.iter.harmony = 50,
                              kmeans_init_nstart = 20, 
                              kmeans_init_iter_max = 250,
                              project.dim = F)
  
  
  ## Identify clusters
  data <- data %>% 
    RunUMAP(reduction = "harmony", dims = 1:NPCs, return.model = T) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:NPCs) %>% 
    FindClusters(resolution = 0.25) %>% 
    identity()
  
  
  # Save the individual Seurat file
  dir.create(paste0(".../Step6.SubclusterAnalysi/", Cluster),
                        showWarnings = FALSE)
  SaveDir <- paste0(".../Step6.SubclusterAnalysi/", Cluster)
  
  saveRDS(data, paste0(SaveDir, "/", Cluster, ".Harmony.", Theta, ".", NPCs, ".rds"))
}

addAlgorithm(name = "HarmoneyRecluster", fun = HarmoneyRecluster)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design <- CJ(i = 1:7)

# final problem design
pdes <- list(ReclusterMajorCluster = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(HarmoneyRecluster = CJ(Theta = c(1.0, 1.5, 2.0),
                                    NPCs = c(20, 30, 40, 50)))

addExperiments(pdes, ades) # pdes set parameter to the data part; ades set parameters to the instance parts

# summarize the experiment
summarizeExperiments()

# Submitting and Collecting Results
options(batchtools.progress = FALSE)
submitJobs() 
getStatus()
# getJobTable()
# getErrorMessages()
# findExpired() 

