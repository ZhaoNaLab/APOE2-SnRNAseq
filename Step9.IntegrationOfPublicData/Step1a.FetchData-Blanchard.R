##########################################################################################################################################################
## Note: the data used for this analysis were downloaded from
## https://www.synapse.org/#!Synapse:syn21682120
## load required libraries
library("tidyverse")
library("synapser")
library("biomaRt")
library("ggbeeswarm")
library("magrittr")
library("SingleCellExperiment")

## download the Normalised, covariates and surrogate variable adjusted residual expression 
synapser::synLogin(email = "YourEmail@gmail.com", password = "YourPassword")

############################ get the matrix data from synapse #####################
WD <- ".../Step9.IntegrationOfPublicData/Blanchard"
setwd(WD)
## Get the raw count matrix
# metadata from synapse
Raw.count <- synapser::synGet("syn45340944", downloadLocation = paste0(WD, "/RawDat"))

## Fetch the gene name of the data
GeneNames <- synapser::synGet("syn45339741", downloadLocation = paste0(WD, "/RawDat"))
GeneNames <- read.csv("RawDat/gene_names.csv") 
GeneNames$V2 <- make.unique(GeneNames$V2)
GeneNames %<>% column_to_rownames("V2")
GeneNames$Features <- rownames(GeneNames)
GeneNames$V1 <- NULL

## Fetch the Column Metadata and modify it 
ColumnMetadata <- synapser::synGet("syn45339740", downloadLocation = paste0(WD, "/RawDat"))
ColumnMetadata <- read.csv("RawDat/column_metadata.csv", row.names = 1)
Metadata.Ind <- ColumnMetadata %>% remove_rownames() %>% unique()

## get the ROSMAP metadata
Metadat.individual <- synapser::synGet("syn3191087", downloadLocation = WD)
ROSMAP_clinical <- read.csv("ROSMAP_clinical.csv")

## add more information to the indvidual metadata
Metadata.Ind %<>% left_join(ROSMAP_clinical %>% 
                              dplyr::select(projid, braaksc, ceradsc, cogdx, individualID),
                            by = "projid")


write.csv(Metadata.Ind, "Metadata.Ind.csv")

## add more information to the cell metadata
ColumnMetadata.updated <- ColumnMetadata %>% rownames_to_column("Cells") %>%
  left_join(ROSMAP_clinical %>% 
              dplyr::select(projid, braaksc, ceradsc, cogdx, individualID),
            by = "projid") %>%
  column_to_rownames("Cells")

write.csv(ColumnMetadata.updated, "ColumnMetadata.updated.csv")


## construct the SCE object
# read the Matrix
Mat <- Matrix::readMM("RawDat/raw_counts.mtx")
colnames(Mat) <- rownames(ColumnMetadata.updated)
rownames(Mat) <- GeneNames$V2

## Construct the SCE object
sce.all <- SingleCellExperiment(list(counts = Mat), 
                                colData = ColumnMetadata.updated,
                                rowData = GeneNames)

saveRDS(sce.all, "sce.all.rds")
## end of the code
