library("harmony")
library("Seurat")
library("tidyverse")
library("data.table")

setwd(".../Step10.PublicDataIntegration")

# load the Seurat object
Seurat.external <- readRDS("Seurat.external.rds")

# run harmony for data integration
options(repr.plot.height = 2.5, repr.plot.width = 6)
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
  FindClusters(resolution = 0.3) %>% 
  identity()


### QC
# Dimplot by individual
DimPlot(Seurat.external.ind, label = T) &
  # split.by = "orig.ident", ncol = 15) &
  NoLegend()

# Featureplot
FeaturePlot(Seurat.external.ind, features = c("nFeature_RNA", "decontX_contamination", "percent.mt"))

## Orig.ident by cluster
FreqTable <- Seurat.external.ind@meta.data %>% 
  select(orig.ident, seurat_clusters) %>% 
  group_by(orig.ident, seurat_clusters) %>%
  summarise(NumberByCluster = n()) %>%
  group_by(orig.ident) %>%
  mutate(NumberByID = sum(NumberByCluster)) %>%
  mutate(ClusterPerc = NumberByCluster/NumberByID * 100)

ggplot(aes(x = orig.ident, y = ClusterPerc, fill = seurat_clusters), data = FreqTable) +
  geom_bar(stat = "identity")

# Dimplot by individual
IDs <- split(unique(Seurat.external.ind$orig.ident), ceiling(seq_along(unique(Seurat.external.ind$orig.ident))/4))

for (i in 1:length(IDs)){
  jpeg(paste0("Plots.Ident/orig.ident", i, ".jpeg"), 
       width = 12.0, height = 4.0, units = "in", res = 600)
  
  p <- DimPlot(subset(Seurat.external.ind, orig.ident %in% unlist(IDs[i])), label = F, 
               split.by = "orig.ident",
               ncol = 4) &
    NoLegend()
  
  print(p)
  dev.off()
}

## Save the data
saveRDS(Seurat.external, "Seurat.external.rds")
