library("Seurat")
library("tidyverse")
library("data.table")
library("batchtools")


setwd(".../Step10.PublicDataIntegration")

# Intergrated datasets at the individual level
Seurat.external.ind <- readRDS("Seurat.external.rds")

## Recluster with a higher resolution
Seurat.external.ind <- Seurat.external.ind %>%
  FindClusters(resolution = 1.0) %>%
  identity()

# visualize the clusters
DimPlot(Seurat.external.ind, label = T) & NoLegend()

## Remove clusters expressing marker genes of different cell type
DefaultAssay(Seurat.external.ind) <- "RNA"
Seurat.external.ind <- NormalizeData(Seurat.external.ind)

# marker gene expression
DotPlot(Seurat.external.ind, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1",
                                          "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                                          "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN",
                                          "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1",
                                          "PDGFRB"), assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")

## Define the major clusters based on markers
Seurat.external.ind$MajorCluster <- NA
Seurat.external.ind$MajorCluster[Seurat.external.ind$seurat_clusters %in%
                                   c(2, 5, 8, 10, 11, 17, 18, 20, 22, 25)] <- "Ex"
Seurat.external.ind$MajorCluster[Seurat.external.ind$seurat_clusters %in%
                                   c(0, 3, 15, 23, 24, 27, 28)] <- "Olig"
Seurat.external.ind$MajorCluster[Seurat.external.ind$seurat_clusters %in%
                                   c(1, 16, 26)] <- "Ast"
Seurat.external.ind$MajorCluster[Seurat.external.ind$seurat_clusters %in%
                                   c(6, 7, 9, 14, 19)] <- "In"
Seurat.external.ind$MajorCluster[Seurat.external.ind$seurat_clusters %in%
                                   c(4)] <- "Mic"
Seurat.external.ind$MajorCluster[Seurat.external.ind$seurat_clusters %in%
                                   c(7, 21)] <- "OPC"
Seurat.external.ind$MajorCluster[Seurat.external.ind$seurat_clusters %in%
                                   c(13)] <- "Vascular"
Seurat.external.ind$MajorCluster[Seurat.external.ind$seurat_clusters %in%
                                   c(12)] <- "Unknown"

# remove mixed cell population. With a notic of perct.mt and contamination on c
# remove mixed cell population and unknown cell population (it is not our focus in this study)
Seurat.external.ind <- subset(Seurat.external.ind, MajorCluster != "Unknown")
table(Seurat.external.ind$MajorCluster)

## Remove individuals with cell number of cell < 400, including S12, S26, S41, S80, S98 
Excluded.orig.ident <- c("S26", "S80", "S41", "S78", "S44", "S75")
Included.Orig.ident <- setdiff(unique(Seurat.external.ind$orig.ident), Excluded.orig.ident)

# examine the operational results
print(Included.Orig.ident)
intersect(Included.Orig.ident, Excluded.orig.ident)

Seurat.external.ind <- subset(Seurat.external.ind, subset = (orig.ident %in% Included.Orig.ident))
sort(table(Seurat.external.ind$orig.ident))
intersect(unique(Seurat.external.ind$orig.ident), Excluded.orig.ident)

#### Recluster again
## use the bath tools for RNA contamination detection and doublet removal
setwd(".../Step10.PublicDataIntegration")
reg <- makeExperimentRegistry("ReclusterMajorCluster.Registry",
                              packages = c("tidyverse", "Seurat", "data.table", "harmony"), 
                              seed = 20221005)

reg <- loadRegistry("ReclusterMajorCluster.Registry", writeable = TRUE)

## Define instance
# define instance (data preparation)
subsample <- function(data, job, i) {
  MajorCluster <- c("Ex", "In", "Ast", "Olig", "OPC", "Mic", "Vascular")
  Cluster <- MajorCluster[i]
  data <- readRDS(data)
  data <- subset(data, MajorCluster == Cluster)
  data <- DietSeurat(data, assays = c("RNA", "decontX"), counts = TRUE)
  DefaultAssay(data) <- "RNA"
  list(data = data, Cluster = Cluster)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "ReclusterMajorCluster", 
           data = Dir, 
           fun = subsample) 

# Define tasks
HarmoneyRecluster <- function(data, job, instance) { # instance is inherited from the last step
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
                 npcs = 50, 
                 pc.genes = data@var.genes,
                 verbose = TRUE)
  
  
  ## integration based on individual
  data <- data %>% RunHarmony("orig.ident", 
                              plot_convergence = TRUE,
                              max.iter.harmony = 50)
  
  
  ## Identify clusters
  data <- data %>% 
    RunUMAP(reduction = "harmony", dims = 1:50) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
    FindClusters(resolution = 1.0) %>% 
    identity()
  
  
  # Save the individual Seurat file
  SaveDir <- ".../Step10.PublicDataIntegration/ClusterCleaning"
  
  saveRDS(data, paste0(SaveDir, "/", Cluster, ".rds"))
}

addAlgorithm(name = "HarmoneyRecluster", fun = HarmoneyRecluster)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design <- CJ(i = 1:7)

# final problem design
pdes <- list(ReclusterMajorCluster = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(HarmoneyRecluster = data.table())

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

## The above cleaning and reclustering process will be done multiple times until
## no mixed cluster population can be identified

##### Now, we are going to merge all cleaned subclusters and do the clustering again
SeuratObjs <- list.files("ClusterCleaning", pattern = ".rds")

# load the files
Seurat.list <- lapply(SeuratObjs, readRDS)

# Merge into one single Seurat object
Seurat.external <- merge(x = Seurat.list[[1]], y = unlist(lapply(2:7, function(x){Seurat.list[[x]]})))
rm(Seurat.list)

# further remove of contaminated cells
Seurat.external <- subset(Seurat.external, decontX_contamination < 0.4 & percent.mt <5)
Seurat.external <- DietSeurat(Seurat.external, assays = c("RNA", "decontX"), counts = TRUE)

# run harmony for data integration
options(repr.plot.height = 2.5, repr.plot.width = 6)

DefaultAssay(Seurat.external) <- "RNA"
Seurat.external <- FindVariableFeatures(Seurat.external, 
                                        selection.method = "vst", 
                                        nfeatures = 3000,
                                        verbose = F)
Seurat.external <- ScaleData(Seurat.external, 
                             vars.to.regress = c("nFeature_RNA", "percent.mt", "decontX_contamination"), 
                             verbose = FALSE)

Seurat.external <- RunPCA(Seurat.external, 
                          npcs = 50, 
                          pc.genes = Seurat.external@var.genes,
                          verbose = TRUE)


## integration based on individual
Seurat.external <- Seurat.external %>% RunHarmony("orig.ident", 
                                                  plot_convergence = TRUE,
                                                  max.iter.harmony = 50)


## Identify clusters
Seurat.external <- Seurat.external %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.25) %>% 
  identity()


setwd(".../Step10.PublicDataIntegration")
saveRDS(Seurat.external, "Seurat.external.final.rds")

## end of the code
















