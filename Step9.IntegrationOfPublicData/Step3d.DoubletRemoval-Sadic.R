library("Seurat")
library("tidyverse")

setwd(".../Step9.IntegrationOfPublicData/Sadick/SeuratObj.Ind")
SeuratObjs <- list.files()

# load the files
Seurat.list <- lapply(SeuratObjs, readRDS)

# Merge into one single Seurat object
Seurat.Sadick <- merge(x = Seurat.list[[1]], y = unlist(lapply(2:length(Seurat.list), function(x){Seurat.list[[x]]})))
rm(Seurat.list)

# add the individual information to the metadata
Meta.cell <- Seurat.Sadick@meta.data %>% 
  rownames_to_column("Cells") 

# recode orig.ident
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor1.1"] <- "D1.1"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor1.2"] <- "D1.2"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor10"] <- "D10"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor11"] <- "D11"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor12"] <- "D12"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor13"] <- "D13"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor15"] <- "D15"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor16"] <- "D16"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor17"] <- "D17"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor2"] <- "D2"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor3"] <- "D3"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor4"] <- "D4"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor5"] <- "D5"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor6"] <- "D6"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor7"] <- "D7"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor8"] <- "D8"
Meta.cell$orig.ident[Meta.cell$orig.ident == "Donor9"] <- "D9"

## load the metadata
Metadata.Sadick <- read.csv(".../Step9.IntegrationOfPublicData/Sadick/Metadata.Ind.csv")
Metadata.Sadick %<>% dplyr::select(DONOR_NUMBER, DISDX, AGE, SEX, APOE, 
                                   RIN, PMI, A, B, C, THAL, BRAAK, CERAD)
names(Metadata.Sadick) <- c("orig.ident", "Dx", "Age", "Sex", "Apoe", "RIN", "PMI",
                            "A", "B", "C", "Thal", "Braak", "CERAD")

# Scale the age
Metadata.Sadick$zAge <- scale(Metadata.Sadick$Age)
Metadata.Sadick$orig.ident[Metadata.Sadick$orig.ident == "D1_replicate1"] <- 'D1.1'
Metadata.Sadick$orig.ident[Metadata.Sadick$orig.ident == "D1_replicate2"] <- 'D1.2'


# add sample information to the metadata
Meta.cell <- Meta.cell %>% left_join(Metadata.Sadick, by = "orig.ident") %>%
  column_to_rownames("Cells")


## Add the metadata to the Seurat object
Seurat.Sadick <- AddMetaData(Seurat.Sadick, metadata = Meta.cell)

# remove some varibles
Seurat.Sadick$nFeature_decontX <- NULL
Seurat.Sadick$nCount_decontX <- NULL
Seurat.Sadick$nCount_originalexp  <- NULL
Seurat.Sadick$nFeature_originalexp <- NULL

## make the new varible cngeneson (scaled cell detection rate)
Seurat.Sadick$cngeneson <- scale(Seurat.Sadick$nFeature_RNA)

setwd(".../Step9.IntegrationOfPublicData/Sadick")
saveRDS(Seurat.Sadick, "Seurat.Sadick.rds")

## Check the doublets
VlnPlot(Seurat.Sadick, features = c("nFeature_RNA", "percent.mt"),
        group.by = "DFClass", pt.size = 0)



## subset the singlet
Seurat.Sadick <- subset(Seurat.Sadick, DFClass == "Singlet")
saveRDS(Seurat.Sadick, "Seurat.Sadick.Single.rds")

