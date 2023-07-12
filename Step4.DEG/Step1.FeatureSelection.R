library("Seurat")
library("biomaRt")

setwd("path/to/project")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl",
                host = 'www.ensembl.org')

ProteinEncodingGenes <- biomaRt::getBM(attributes = c("hgnc_symbol", "transcript_biotype"), 
                                       filters = "transcript_biotype",
                                       values = "protein_coding", 
                                       mart = mart)

saveRDS(ProteinEncodingGenes, "Step4.DEG/ProteinEncodingGenes.rds")

# load the data
SeuratObj <- readRDS("Step3.Clustering/SeuratObj.rds")

# Calculate the percentage of cells expressing the genes
DefaultAssay(SeuratObj) <- "RNA"

MajorCellType <- unique(SeuratObj$MajorCluster)

FeatureSelect <- list()

for (i in 1:length(MajorCellType)){
  SeuratObj <- subset(SeuratObj, MajorCluster == MajorCellType[i])
  ExpMat <- GetAssayData(SeuratObj, assay = "RNA", slot = "counts")      
  PercMat <- as.matrix(rowMeans(ExpMat > 0))*100
  colnames(PercMat) <- "Percentage"

  ## include genes expressed by expressed by at least 10% of the cells
  FeaturesIncluded <- rownames(subset(as.data.frame(PercMat), Percentage >= 10))
  
  ## Use only protein-encoding genes
  FeaturesIncluded <- FeaturesIncluded[FeaturesIncluded %in% ProteinEncodingGenes$hgnc_symbol]
  FeatureSelect[[i]] <- FeaturesIncluded
  
}

names(FeatureSelect) <- MajorCellType

# save the data
saveRDS(FeatureSelect, "Step4.DEG/FeatureSelect.rds")





