library("Seurat")
library("tidyverse")
library("ggpubr")
library("purrr")
library("biomaRt")


############# step1: load the Single cell experiment objects #############
setwd(".../Step10.PublicDataIntegration")
Seurat.external.filter <- readRDS("Seurat.external.filter.rds")

## QC for doublet identification
FeatureToPlot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination")
for (i in 1:length(FeatureToPlot)){
  p <- VlnPlot(Seurat.external.filter,
               features = FeatureToPlot[i],
               split.by = "orig.ident", 
               group.by = "Study",
               pt.size = 0,
               ncol = 1) & NoLegend()
  p <- p + geom_jitter(data = p$data, 
                       position = position_jitterdodge(jitter.width = 8.0, dodge.width = 0.9),
                       size = 0.1, alpha = 0.01)
  ggsave(paste0("DoubletQC.", FeatureToPlot[i], ".pdf"), width = 25, height = 4.0, plot = p)
}


######################### step 2b: Filtering: remove cells with less than 200 features ######################## 
# Check the most abundant features (plot the median percentage of the features)
AbundantMat <- Seurat.external.filter@assays$RNA@counts
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

# Filter Mitocondrial
Seurat.external.filter <- Seurat.external.filter[!grepl("^MT-", rownames(Seurat.external.filter)), ]
Seurat.external.filter <- Seurat.external.filter[!grepl("^HB[^(P)]", rownames(Seurat.external.filter)), ]

## update the gene detection rate
Seurat.external.filter$cngeneson <- scale(Seurat.external.filter$nFeature_RNA)


######################## step 2d: Sample sex information confirmation ########################
# Get the chorosome information of each gene
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# fetch chromosome info plus some other annotations
genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                                 "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
                                  useCache = F))

genes.table <- read.csv(".../Step10.PublicDataIntegration/genes.table.csv")
chrY.gene <- genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
SelectY <- intersect(rownames(Seurat.external.filter), chrY.gene)

# calculate the percentage of Y chrosome gene
Seurat.external.filter$pct_chrY <- colSums(Seurat.external.filter@assays$RNA@counts[SelectY, ])/colSums(Seurat.external.filter@assays$RNA@counts)
FeatureScatter(Seurat.external.filter, feature1 = "XIST", feature2 = "pct_chrY")


# Feature plot to show sex gene information
FeatureToPlot <- c("XIST", "pct_chrY")
for (i in 1:length(FeatureToPlot)){
  p <- VlnPlot(Seurat.external.filter,
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
Seurat.external.filter <- NormalizeData(Seurat.external.filter)
Seurat.external.filter <- CellCycleScoring(object = Seurat.external.filter, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

# Feature plot to show cell cycle information
FeatureToPlot <- c("S.Score", "G2M.Score")
for (i in 1:length(FeatureToPlot)){
  p <- VlnPlot(Seurat.external.filter,
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

saveRDS(Seurat.external.filter, "Seurat.external.rds")
