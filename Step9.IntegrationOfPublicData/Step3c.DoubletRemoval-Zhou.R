library("Seurat")
library("tidyverse")
library("magrittr")

setwd(".../Step9.IntegrationOfPublicData/Zhou/SeuratObj.Ind")
SeuratObjs <- list.files()

# load the files
Seurat.list <- lapply(SeuratObjs, readRDS)

# Merge into one single Seurat object
Seurat.Zhou <- merge(x = Seurat.list[[1]], y = unlist(lapply(2:length(Seurat.list), function(x){Seurat.list[[x]]})))
rm(Seurat.list)

# add the individual information to the metadata
Meta.cell <- Seurat.Zhou@meta.data %>% 
  rownames_to_column("Cells") 


## load the metadata
Metadata.Zhou <- read.csv(".../Step9.IntegrationOfPublicData/Zhou/Metadata.csv", row.names = 1)
Metadata.Zhou %<>% rownames_to_column("orig.ident")

# Rename some varible
names(Metadata.Zhou)[names(Metadata.Zhou) == "msex"] <- "Sex"
names(Metadata.Zhou)[names(Metadata.Zhou) == "age_death"] <- "Age"
names(Metadata.Zhou)[names(Metadata.Zhou) == "pmi"] <- "PMI"
names(Metadata.Zhou)[names(Metadata.Zhou) == "apoe_genotype"] <- "Apoe"
names(Metadata.Zhou)[names(Metadata.Zhou) == "braaksc"] <- "Braak"
names(Metadata.Zhou)[names(Metadata.Zhou) == "ceradsc"] <- "CERAD"

## recode some varibles
Metadata.Zhou$Age[Metadata.Zhou$Age == "90+"] <- 90
Metadata.Zhou$Age <- round(as.numeric(Metadata.Zhou$Age), 1)

Metadata.Zhou %<>% mutate(Sex = dplyr::recode(as.character(Sex), "0" = "F", "1" = "M"))

## scale the age
Metadata.Zhou$zAge <- scale(Metadata.Zhou$Age)

# add sample information to the metadata
Meta.cell <- Meta.cell %>% left_join(Metadata.Zhou, by = "orig.ident") %>%
  column_to_rownames("Cells")


## Add the metadata to the Seurat object
Seurat.Zhou <- AddMetaData(Seurat.Zhou, metadata = Meta.cell)

# remove some varibles
Seurat.Zhou$nFeature_decontX <- NULL
Seurat.Zhou$nCount_decontX <- NULL
Seurat.Zhou$nCount_originalexp  <- NULL
Seurat.Zhou$nFeature_originalexp <- NULL

## make the new varible cngeneson (scaled cell detection rate)
Seurat.Zhou$cngeneson <- scale(Seurat.Zhou$nFeature_RNA)

## Check the doublets
VlnPlot(Seurat.Zhou, features = c("nFeature_RNA", "percent.mt"),
        group.by = "DFClass", pt.size = 0)

setwd(".../Step9.IntegrationOfPublicData/Zhou")
# saveRDS(Seurat.Zhou, "Seurat.Zhou.rds")
# 

## subset the singlet
Seurat.Zhou <- subset(Seurat.Zhou, DFClass == "Singlet")
saveRDS(Seurat.Zhou, "Seurat.Zhou.Single.rds")

