library("Seurat")
library("tidyverse")
library("ggpubr")
library("purrr")
library("biomaRt")


############# step1: load the Single cell experiment objects #############
setwd("path/to/project")
SeuratObj <- readRDS("Step1.DataIntegration/merged.Seurat.rds")

setwd(".../Step2.QC")

## QC for doublet identification
FeatureToPlot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination")
for (i in 1:length(FeatureToPlot)){
  p <- VlnPlot(SeuratObj,
               features = FeatureToPlot[i],
               split.by = "DFClass", 
               group.by = "orig.ident",
               pt.size = 0,
               ncol = 1) 
  p <- p + geom_jitter(data = p$data, 
                       position = position_jitterdodge(jitter.width = 8.0, dodge.width = 0.9),
                       size = 0.1, alpha = 0.01)
  ggsave(paste0("DoubletQC.", FeatureToPlot[i], ".pdf"), width = 20, height = 4.0, plot = p)
}


## subset the singlet and remove cell with decontX_contamination >= 0.4 
SeuratObj <- subset(SeuratObj, DFClass == "Singlet" & decontX_contamination < 0.4)

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
  ggsave(paste0("QC.", FeatureToPlot[i], ".pdf"), width = 20, height = 4.0, plot = p)
}

# feature scatter plots
pdf("QC.FeatureScatter.pdf", width = 6.0, height = 6)
FeatureScatter(SeuratObj, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
dev.off()

######################### step 2b: Filtering: remove cells with less than 200 features ######################## 
# Check the most abundant features (plot the median percentage of the features)
AbundantMat <- SeuratObj@assays$RNA@counts
AbundantMat <- Matrix::t(Matrix::t(AbundantMat)/Matrix::colSums(AbundantMat)) * 100
# saveRDS(AbundantMat, "AbundantMat.rds")

HigestExpressed <- order(apply(AbundantMat, 1, median), decreasing = T)[20:1]
# saveRDS(HigestExpressed, "HigestExpressed.rds")

# boxplot show the median
par(mar = c(6, 10, 2, 1))
pdf("FeatureAbundance.pdf", width = 6.0, height = 4.0)

boxplot(t(as.matrix(AbundantMat[HigestExpressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

dev.off()

## MALAT1 is the most abundant gene, check the expression of MALAT1 in each 
VlnPlot(SeuratObj,
        features = "MALAT1",
        group.by = "orig.ident",
        pt.size = 0.01,
        ncol = 1) 


# optional: filter genes
# Filter MALAT1: be noted that Seurat autmoatically udpate n_Feature and n_count 
# after the removal of the gene
SeuratObj <- SeuratObj[!grepl("MALAT1", rownames(SeuratObj)), ]

# Filter Mitocondrial
SeuratObj <- SeuratObj[!grepl("^MT-", rownames(SeuratObj)), ]

# Filter Ribossomal gene (optional if that is a problem on your data) SeuratObj
# SeuratObj <- SeuratObj[ ! grepl('^RP[SL]', rownames(SeuratObj)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
SeuratObj <- SeuratObj[!grepl("^HB[^(P)]", rownames(SeuratObj)), ]

## update the gene detection rate
SeuratObj$cngeneson <- scale(SeuratObj$nFeature_RNA)


######################## step 2d: Sample sex information confirmation ########################
# Get the chorosome information of each gene
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# fetch chromosome info plus some other annotations
genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                                 "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
                                  useCache = F))

genes.table <- read.csv("genes.table.csv")
chrY.gene <- genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
SelectY <- intersect(rownames(SeuratObj), chrY.gene)

# calculate the percentage of Y chrosome gene
SeuratObj$pct_chrY <- colSums(SeuratObj@assays$RNA@counts[SelectY, ])/colSums(SeuratObj@assays$RNA@counts)
FeatureScatter(SeuratObj, feature1 = "XIST", feature2 = "pct_chrY")


# Feature plot to show sex gene information
FeatureToPlot <- c("XIST", "pct_chrY")
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
  ggsave(paste0("QC.", FeatureToPlot[i], ".pdf"), width = 20, height = 6.0, plot = p)
}

######################## step 2e: Calculate cell-cycle scores ########################
# Before running CellCycleScoring the data need to be normalized and
# logtransformed.
SeuratObj <- NormalizeData(SeuratObj)
SeuratObj <- CellCycleScoring(object = SeuratObj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

# Feature plot to show cell cycle information
FeatureToPlot <- c("S.Score", "G2M.Score")
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
  ggsave(paste0("QC.", FeatureToPlot[i], ".pdf"), width = 20, height = 6.0, plot = p)
}

## save the Seurat object for clustering analysis
saveRDS(SeuratObj, "SeuratObj.QC.rds")
