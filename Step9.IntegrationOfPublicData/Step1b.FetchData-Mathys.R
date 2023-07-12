##########################################################################################################################################################
## Note: the data used for this analysis were downloaded from
## https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn18485175
## load required libraries
library(tidyverse)
library(synapser)
library(biomaRt)
library(ggbeeswarm)
library(Matrix)
library(magrittr)

## download the Normalised, covariates and surrogate variable adjusted residual expression 
synapser::synLogin(email = "YourEmail@gmail.com", password = "YourPassword")

## use the Cell ranger output
SynapseID <- c("syn18686383", "syn18687958", "syn18687959")
Dir <- ".../Step9.IntegrationOfPublicData/Mathys/Raw"
for (i in 1:length(SynapseID)){
  Dat <- synapser::synGet(SynapseID[i], downloadLocation = Dir)
}

# ## in the bash
# mv CellRangerOutput_barcodes.tsv barcodes.tsv
# mv CellRangerOutput_genes.tsv genes.tsv
# mv " CellRangerOutput_matrix.mtx" matrix.mtx

setwd(".../Step9.IntegrationOfPublicData/Mathys")
## Get the metadta
Metadat.Biospecimen <- synapser::synGet("syn18642936")
Metadat.Biospecimen <- read.csv(Metadat.Biospecimen$path)

Metadat.individual <- synapser::synGet("syn3191087")
Metadat.individual <- read.csv(Metadat.individual$path)

## merge the two datasets
Metadata <- Metadat.Biospecimen %>% 
  dplyr::select(projid, tissue, tissueWeight, individualID) %>%
  left_join(Metadat.individual %>% 
              dplyr::select(individualID, Study, msex, educ, race, 
                            apoe_genotype, age_death, pmi, braaksc, ceradsc,
                            cogdx), by = "individualID") %>%
  arrange(age_death)

dim(Metadata)


## Update the AD Dx information based on the metadata published in the Nature paper
Meta.paper <- read.csv("MetafromPaper.csv") %>%
  arrange(age_death)

## According the age mathing
Metadata$PathDx <- NA
Metadata$PathDx[1:36] <- Meta.paper$pathologic.diagnosis.of.AD[1:36]

## Further complement the Dx information based on the unique information
## of Education and cogdx
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 16 & Metadata$cogdx == 2 & Metadata$msex == 0] <- "NO"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 20 & Metadata$cogdx == 4 & Metadata$msex == 0] <- "YES"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 23 & Metadata$cogdx == 1 & Metadata$msex == 0] <- "NO"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 21 & Metadata$cogdx == 1 & Metadata$msex == 1] <- "NO"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 18 & Metadata$cogdx == 2 & Metadata$msex == 0] <- "NO"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 18 & Metadata$cogdx == 1 & Metadata$msex == 0] <- "NO"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 21 & Metadata$cogdx == 4 & Metadata$msex == 0] <- "YES"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 21 & Metadata$cogdx == 4 & Metadata$msex == 0] <- "YES"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 24 & Metadata$cogdx == 1 & Metadata$msex == 0] <- "NO"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 24 & Metadata$cogdx == 4 & Metadata$msex == 0] <- "YES"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 18 & Metadata$cogdx == 4 & Metadata$braaksc == 5] <- "YES"
Metadata$PathDx[Metadata$age_death == "90+" & Metadata$educ == 18 & Metadata$cogdx == 4 & Metadata$braaksc == 2] <- "NO"

## save the Metadata
write.csv(Metadata, "Metadata.csv", row.names = FALSE)

## Construct SingleCellExperiment
Counts <- readMM("RawDat/notfiltered_count_matrix.mtx")
RowMeatada <- read.delim("RawDat/notfiltered_gene_row_names.txt",
                         header = F) 

RowMeatada$V2 <- make.unique(RowMeatada$V2)
RowMeatada %<>% column_to_rownames("V2")
RowMeatada$V1 <- NULL
RowMeatada$Features <- rownames(RowMeatada)

# modify the column metadata
ColMetadata <- read.delim("RawDat/notfiltered_column_metadata.txt") %>%
  left_join(Metadata, by = "projid") %>%
  column_to_rownames("TAG")

## add rownames and column names to the count matrix
rownames(Counts) <- rownames(RowMeatada)
colnames(Counts) <- rownames(ColMetadata)

# Create a SingleCellExperiment object and run decontX
sce.all <- SingleCellExperiment(list(counts = Counts),
                                colData = ColMetadata,
                                rowData = RowMeatada)


saveRDS(sce.all, "sce.all.rds")


