library("Seurat")
library("tidyverse")

# Set the working directory
setwd("path/to/project")

SeuratObj <- readRDS("Step3.Clustering/SeuratObj.rds")

DimPlot(SeuratObj, label = T)
FeaturePlot(SeuratObj, features = c("nFeature_RNA", "percent.mt", "decontX_contamination"))

# marker gene expression
Idents(SeuratObj) <- SeuratObj$seurat_clusters
NormalizeData(SeuratObj, assay = "RNA")
DotPlot(SeuratObj, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                                      "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                                      "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                                      "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                                      "PDGFRB"), 
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

## Rename the subclusters
SeuratObj$Subcluster <- NA
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 0] <- "Olig"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 1] <- "Ex1"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 3] <- "Ex2"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 6] <- "Ex3"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 7] <- "Ex4"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 9] <- "Ex5"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 12] <- "Ex6"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 13] <- "Ex7"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 14] <- "Ex8"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 15] <- "Ex9"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 16] <- "Ex10"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 18] <- "Ex11"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 22] <- "Ex12"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 20] <- "In1"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 2] <- "Ast1"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 24] <- "Ast2"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 8] <- "In1"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 10] <- "In2"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 11] <- "In3"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 17] <- "In4"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 19] <- "In5"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 23] <- "In6"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 5] <- "OPC"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 4] <- "Mic"
SeuratObj$Subcluster[SeuratObj$seurat_clusters == 21] <- "Vasculature"


## Define Major clusters
SeuratObj$MajorCluster[grepl("Ex", SeuratObj$Subcluster, fixed = TRUE)] <- "Ex"
SeuratObj$MajorCluster[grepl("In", SeuratObj$Subcluster, fixed = TRUE)] <- "In"
SeuratObj$MajorCluster[grepl("Olig", SeuratObj$Subcluster, fixed = TRUE)] <- "Olig"
SeuratObj$MajorCluster[grepl("Mic", SeuratObj$Subcluster, fixed = TRUE)] <- "Mic"
SeuratObj$MajorCluster[grepl("Ast", SeuratObj$Subcluster, fixed = TRUE)] <- "Ast"
SeuratObj$MajorCluster[grepl("OPC", SeuratObj$Subcluster, fixed = TRUE)] <- "OPC"
SeuratObj$MajorCluster[grepl("Vasculature", SeuratObj$Subcluster, fixed = TRUE)] <- "Vascu"

# set the subcluser as the Ident
Idents(SeuratObj) <- SeuratObj$MajorCluster
DimPlot(SeuratObj, label = T)

################################## QC Metrics ##########################
# Plot1: Dimention plot
# Plot1a. Combined Dimention plot by subclusters
Idents(SeuratObj) <- SeuratObj$Subcluster
jpeg("QC/DimPlotSubcluster.jpeg", width = 10.0, height = 10.0, unit = 'in', res = 600)

DimPlot(SeuratObj, 
        label = T, 
        label.size = 8.0,
        raster = FALSE) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(size = 25, color = 'black'),
        axis.text = element_text(size = 25, color = 'black'),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm")) &
  NoLegend()

dev.off()


## Subclusters by APOE genotype, Sex, or AD diagnosis
# By APOE
jpeg("QC/DimPlotSubcluster.By.APOE.jpeg", width = 12, height = 4, unit = 'in', res = 600)
DimPlot(SeuratObj, 
        split.by = "apoe",
        label = FALSE, 
        # label.size = 8.0,
        raster = FALSE,
        pt.size = 0.1) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 20, color = 'black'),
        axis.line = element_line(linewidth = 1.0),
        axis.ticks = element_line(linewidth = 1.0),
        axis.ticks.length = unit(0.25, "cm")) &
  NoLegend()

dev.off()

# By Sex
jpeg("QC/DimPlotSubcluster.By.Sex.jpeg", width = 8, height = 4, unit = 'in', res = 600)
DimPlot(SeuratObj, 
        split.by = "sex",
        label = FALSE, 
        # label.size = 8.0,
        raster = FALSE,
        pt.size = 0.1) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 20, color = 'black'),
        axis.line = element_line(linewidth = 1.0),
        axis.ticks = element_line(linewidth = 1.0),
        axis.ticks.length = unit(0.25, "cm")) &
  NoLegend()

dev.off()

# By bg
jpeg("QC/DimPlotSubcluster.By.Dx.jpeg", width = 8, height = 4, unit = 'in', res = 600)
DimPlot(SeuratObj, 
        split.by = "bg",
        label = FALSE, 
        # label.size = 8.0,
        raster = FALSE,
        pt.size = 0.1) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 20, color = 'black'),
        axis.line = element_line(linewidth = 1.0),
        axis.ticks = element_line(linewidth = 1.0),
        axis.ticks.length = unit(0.25, "cm")) &
  NoLegend()

dev.off()

# Plot1b. Combined Dimention plot by orig.ident
IDs <- paste0("S", 1:length(unique(SeuratObj$orig.ident)))
IDchuncks <- split(IDs, ceiling(seq_along(1:length(IDs))/4))

# make the figure look better
IDchuncks$`15` <- c("S55", "S56", "S57", "S58")

for (i in 1:length(IDchuncks)){
  dat <- subset(SeuratObj, orig.ident %in% unlist(IDchuncks[i]))
  Idents(dat) <- dat$seurat_clusters
  
  jpeg(paste0("QC/SubclusterByID/DimPlo.Chunck", i,".jpeg"), width = 16, height = 5, unit = 'in', res = 600)
  
  p <- DimPlot(dat, split.by = 'orig.ident', label = T, ncol = 4,
               # label.size = 8.0,
               raster = FALSE) +
    xlab("UMAP-1") +
    ylab("UMAP-2") +
    theme(axis.title = element_text(size = 20, color = 'black'),
          axis.text = element_text(size = 20, color = 'black'),
          axis.line = element_line(linewidth = 1.0),
          axis.ticks = element_line(linewidth = 1.0),
          axis.ticks.length = unit(0.25, "cm")) &
    NoLegend()
  
  print(p)
  dev.off()
}


## Major cluster
Idents(SeuratObj) <- SeuratObj$MajorCluster
jpeg("QC/DimPlotMajorCluster.jpeg", 
     width = 10.0, 
     height = 10.0, 
     unit = 'in',
     res = 600)

DimPlot(SeuratObj, 
        label = T, 
        label.size = 8.0,
        raster = FALSE) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(size = 25, color = 'black'),
        axis.text = element_text(size = 25, color = 'black'),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm")) &
  NoLegend()

dev.off()

# Plot1b. Combined Dimention plot by orig.ident
IDs <- paste0("S", 1:length(unique(SeuratObj$orig.ident)))
IDchuncks <- split(IDs, ceiling(seq_along(1:length(IDs))/4))

# make the figure look better
IDchuncks$`15` <- c("S55", "S56", "S57", "S58")

for (i in 1:length(IDchuncks)){
  dat <- subset(SeuratObj, orig.ident %in% unlist(IDchuncks[i]))
  Idents(dat) <- dat$seurat_clusters
  
  jpeg(paste0("QC/MajorClusterByID/DimPlo.Chunck", i,".jpeg"), 
       width = 16, 
       height = 5,
       unit = 'in',
       res = 600)
  
  p <- DimPlot(dat, split.by = 'orig.ident', label = T, ncol = 4,
               # label.size = 8.0,
               raster = FALSE) +
    xlab("UMAP-1") +
    ylab("UMAP-2") +
    theme(axis.title = element_text(size = 20, color = 'black'),
          axis.text = element_text(size = 20, color = 'black'),
          axis.line = element_line(linewidth = 1.0),
          axis.ticks = element_line(linewidth = 1.0),
          axis.ticks.length = unit(0.25, "cm")) &
    NoLegend()
  
  print(p)
  dev.off()
}

## Plot 2a: Violin Plot
FeatureToPlot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination")
for (i in 1:length(FeatureToPlot)){
  p <- VlnPlot(SeuratObj,
               features = FeatureToPlot[i],
               group.by = "orig.ident",
               pt.size = 0,
               ncol = 1) +
    NoLegend()
  p <- p + geom_jitter(data = p$data, 
                       position = position_jitterdodge(jitter.width = 8.0, dodge.width = 0.9),
                       size = 0.1, alpha = 0.01)
  ggsave(paste0("QC/QC.", FeatureToPlot[i], ".pdf"), width = 20, height = 4.0, plot = p)
}


## Plot 2b: Feature plot
jpeg("QC/QC.Features.jpeg", width = 12.0, height = 8.0, unit = "in", res = 600)
FeaturePlot(SeuratObj, 
            features = c("nFeature_RNA", "decontX_contamination", "percent.mt", "percent.rb"),
            cols = c("lightgrey", "#ef4b86"),
            ncol = 2) &
  xlab("UMAP-1") &
  ylab("UMAP-2")

dev.off()

## Plot 3: Dotplot to show markers
## Dot plot to show marker gene expression
# marker gene expression
SeuratObj$Subcluster <- factor(SeuratObj$Subcluster, 
                                     levels = c("Ex1", "Ex2", "Ex3", "Ex4", "Ex5", "Ex6", "Ex7", "Ex8",
                                                "Ex9", "Ex10", "Ex11", "Ex12", "Ex13", "In1", "In2", 
                                                "In3","In4", "In5", "In6", "Olig", "OPC", "Ast1","Ast2",
                                                "Mic", "Vasculature"))


Idents(SeuratObj) <- SeuratObj$Subcluster
NormalizeData(SeuratObj, assay = "RNA")

pdf("MarkerPlot.pdf", width = 12.0, height = 6.0)
DotPlot(SeuratObj, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7",
                                      "GAD1", "GAD2", "GFAP", "AQP4", "GJA1",
                                      "ALDH1L1", "MOBP", "MAG", "MOG", "VCAN", 
                                      "PDGFRA", "CSPG4", "BCAN", "CSF1R", "CD74", 
                                      "C3","CLDN5", "FLT1", "PDGFRB"),
        cols = c("lightgrey", "#c04b57"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

dev.off()

## use decontX
NormalizeData(SeuratObj, assay = "decontX")

pdf("MarkerPlot.decontX.pdf", width = 12.0, height = 6.0)
DotPlot(SeuratObj, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7",
                                      "GAD1", "GAD2", "GFAP", "AQP4", "GJA1",
                                      "ALDH1L1", "MOBP", "MAG", "MOG", "VCAN", 
                                      "PDGFRA", "CSPG4", "BCAN", "CSF1R", "CD74", 
                                      "C3","CLDN5", "FLT1", "PDGFRB"),
        cols = c("lightgrey", "#c04b57"),
        assay = 'decontX') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

dev.off()


# Save the data
saveRDS(SeuratObj, "SeuratObj.rds")









