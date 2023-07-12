library('SingleR')
library('Seurat')
library("readxl")
library('tidyverse')
library('data.table')
library("scuttle")

"""
All the reference datasets can be downloaded as described in their original
publications

"""

## 1. load the Seurat object
setwd("path/to/project")
scaObjs <- readRDS(scaObjs, "Step3.Clustering/SeuratObj.rds")

Idents(scaObjs) <- scaObjs$MajorCluster

# for each of the dataset, including the inquery and reference dataset: 1) make sure the default
# assay is RNA; 2) convert the dataset to the SingleCellExperiment object; 3) Do log-norm transformation
# 1.  set the default assay
DefaultAssay(scaObjs) <- "RNA"

# 2. convert the seurat object to SingleCellExperiment object
ApoeSigexp <- as.SingleCellExperiment(scaObjs)


# 3. Normalize the data and perform log-transformation for 
ApoeSigexp <- logNormCounts(ApoeSigexp)


########################### with Jakel et al., 2019 as the reference dataset #######################
setwd(".../Castelo-BrancoPaper")
# load the expression matrix
Jakel.Mat <- fread('GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt', sep = '\t')
Jakel.Mat <- column_to_rownames(Jakel.Mat, 'V1')
Jakel.Mat <- as.matrix(Jakel.Mat)
Jakel.Meta <- fread('GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt', sep = '\t')
Jakel.Meta <- Jakel.Meta %>% column_to_rownames("Detected")

## Recoding SnRNA_MS_meta$Celltypes
Jakel.Meta$Celltypes <- fct_recode(Jakel.Meta$Celltypes,
                                   "As" = "Astrocytes",
                                   "As" = "Astrocytes2",
                                   "End" = "Endothelial_cells1",
                                   "End" = "Endothelial_cells2",
                                   "Immue" = "Immune_cells",
                                   "Mac" = "Macrophages",
                                   "Mic/Mac" = "Microglia_Macrophages",
                                   "Neuron" = "Neuron1",
                                   "Neuron" = "Neuron2",
                                   "Neuron" = "Neuron3",
                                   "Neuron" = "Neuron4",
                                   "Neuron" = "Neuron5",
                                   "Oligo" = "Oligo1",
                                   "Oligo" = "Oligo2",
                                   "Oligo" = "Oligo3",
                                   "Oligo" = "Oligo4",
                                   "Oligo" = "Oligo5",
                                   "Oligo" = "Oligo6",
                                   "OPC" = "OPCs",
                                   "Smooth" = "Vasc_smooth_muscle"
)


# construct the Seurat object 
Jakel.Seurat <- CreateSeuratObject(counts = Jakel.Mat, project = "JakelMS")

# Add the metadata to the seurat object
Jakel.Seurat <- AddMetaData(object = Jakel.Seurat, metadata = Jakel.Meta)
DefaultAssay(Jakel.Seurat) <- "RNA"

# convert the seurat object to SingleCellExperiment object
Jakel.Sigexp <- as.SingleCellExperiment(Jakel.Seurat)

# 3. Normalize the data and perform log-transformation for 
Jakel.Sigexp <- logNormCounts(Jakel.Sigexp)


########################################## Annotation #################################
Jakel.Pred.single <- SingleR(test=ApoeSigexp, ref=Jakel.Sigexp, labels=Jakel.Sigexp$Celltypes, 
                             de.method="wilcox")


Jakel.Pred.cluster <- SingleR(test=ApoeSigexp, ref=Jakel.Sigexp, labels=Jakel.Sigexp$Celltypes, 
                              method = "cluster", 
                              clusters = ApoeSigexp$MajorCluster, 
                              de.method="wilcox")

setwd("Step3.Clustering")
saveRDS(Jakel.Pred.single, "Jakel.Pred.single.rds")
saveRDS(Jakel.Pred.cluster, "Jakel.Pred.cluster.rds")


########################### with Lake et al., 2018 as the reference dataset #######################
setwd(".../LakePaper")

# load the expression matrix
Lake.Mat <- fread('GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt', sep = '\t')
Lake.Mat <- column_to_rownames(Lake.Mat, 'V1')
Lake.Mat <- as.matrix(Lake.Mat)

# creat the metadata
Lake.Meta <- as.data.frame(colnames(Lake.Mat)) 
Barcode <- Lake.Meta$`colnames(Lake.Mat)`

# Extract the cell type information 
Lake.Meta <- as.data.frame(str_split_fixed(Lake.Meta$`colnames(Lake.Mat)`, "_", 3))
rownames(Lake.Meta) <- Barcode
colnames(Lake.Meta) <- c("SubCellTypes", "Region", "Barcode")

# creat the varible: Major cell type
Lake.Meta <- Lake.Meta %>%
  mutate(MajorCluster = case_when(str_detect(SubCellTypes, "Ex*") ~ "Ex",
                                  str_detect(SubCellTypes, "In*") ~ "In",
                                  str_detect(SubCellTypes, "Ast*") ~ "Ast",
                                  str_detect(SubCellTypes, "End*") ~ "End",
                                  str_detect(SubCellTypes, "Mic*") ~ "Mic",
                                  str_detect(SubCellTypes, "OPC*") ~ "OPC",
                                  str_detect(SubCellTypes, "Oli*") ~ "Olig",
                                  str_detect(SubCellTypes, "Per*") ~ "Per"))



# construct the Seurat object 
Lake.Seurat <- CreateSeuratObject(counts = Lake.Mat, project = "Lake")

# Add the metadata to the seurat object
Lake.Seurat <- AddMetaData(object = Lake.Seurat, metadata = Lake.Meta)
DefaultAssay(Lake.Seurat) <- "RNA"

# convert the seurat object to SingleCellExperiment object
Lake.Sigexp <- as.SingleCellExperiment(Lake.Seurat)

# 3. Normalize the data and perform log-transformation for 
Lake.Sigexp <- logNormCounts(Lake.Sigexp)


########################################## Annotation #################################
setwd("Step3.Clustering")

# mapping by single cell
Lake.Pred.single <- SingleR(test=ApoeSigexp, ref=Lake.Sigexp, labels=Lake.Sigexp$MajorCluster, 
                            de.method="wilcox")

# Mapping by cluster
Lake.Pred.cluster <- SingleR(test=ApoeSigexp, ref=Lake.Sigexp, labels=Lake.Sigexp$MajorCluster, 
                             method = "cluster", 
                             clusters = ApoeSigexp$MajorCluster, 
                             de.method="wilcox")


saveRDS(Lake.Pred.single , "Lake.Pred.single.rds")
saveRDS(Lake.Pred.cluster, "Lake.Pred.cluster.rds")


#################### with Leng et al., 2021 as the reference dataset ##############
# load the data
setwd('.../KampmannPaper')
LengSigexp <- readRDS("sce.EC.scAlign.assigned.rds")

# 3. Normalize the data and perform log-transformation for 
LengSigexp <- logNormCounts(LengSigexp)

########################################## Annotation #################################
Leng.Pred.single <- SingleR(test=ApoeSigexp, ref=LengSigexp, labels=LengSigexp$clusterCellType, 
                            de.method="wilcox")


Leng.Pred.cluster <- SingleR(test=ApoeSigexp, ref=LengSigexp, labels=LengSigexp$clusterCellType, 
                             method = "cluster", 
                             clusters = ApoeSigexp$MajorCluster, 
                             de.method="wilcox")

setwd("Step3.Clustering")
saveRDS(Leng.Pred.single, "Leng.Pred.single.rds")
saveRDS(Leng.Pred.cluster, "Leng.Pred.cluster.rds")

#################### with Mathys et al., 2019 as the reference dataset ##############
setwd(".../TsaiPaper")
TsaiSeurat <- readRDS("TsaiSeurat.rds")

# 1. set the default assay
DefaultAssay(TsaiSeurat) <- "RNA"

# 2. convert the seurat object to SingleCellExperiment object
TsaiSigexp <- as.SingleCellExperiment(TsaiSeurat)

# 3. Normalize the data and perform log-transformation for 
TsaiSigexp <- logNormCounts(TsaiSigexp)

###################### Annotation with Mathys et al., 2019  dataset ######################
Mathys.Pred.single <- SingleR(test=ApoeSigexp, ref=TsaiSigexp, labels=TsaiSigexp$broad.cell.type, 
                              de.method="wilcox")


Mathys.Pred.cluster <- SingleR(test=ApoeSigexp, ref=TsaiSigexp, labels=TsaiSigexp$broad.cell.type, 
                               method = "cluster", 
                               clusters = ApoeSigexp$MajorCluster, 
                               de.method="wilcox")

setwd("Step3.Clustering")
saveRDS(Mathys.Pred.single, "Mathys.Pred.single.rds")
saveRDS(Mathys.Pred.cluster, "Mathys.Pred.cluster.rds")

