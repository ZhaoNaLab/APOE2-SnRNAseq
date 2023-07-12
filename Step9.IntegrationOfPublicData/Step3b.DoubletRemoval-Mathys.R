library("Seurat")
library("tidyverse")

setwd(".../Step9.IntegrationOfPublicData/Mathys/SeuratObj.Ind")
SeuratObjs <- list.files()

# load the files
Seurat.list <- lapply(SeuratObjs, readRDS)

# Merge into one single Seurat object
Seurat.Mathys <- merge(x = Seurat.list[[1]], y = unlist(lapply(2:length(Seurat.list), function(x){Seurat.list[[x]]})))
rm(Seurat.list)

# add the individual information to the metadata
Meta.cell <- Seurat.Mathys@meta.data %>% 
  rownames_to_column("Cells") 

# use Individual ID as orig.ident
Meta.cell$orig.ident <- Meta.cell$individualID

## load the metadata
Metadata.Mathys <- read.csv(".../Step9.IntegrationOfPublicData/Mathys/Metadata.csv", row.names = 1)

# Rename some varible
names(Metadata.Mathys)[names(Metadata.Mathys) == "msex"] <- "Sex"
names(Metadata.Mathys)[names(Metadata.Mathys) == "age_death"] <- "Age"
names(Metadata.Mathys)[names(Metadata.Mathys) == "pmi"] <- "PMI"
names(Metadata.Mathys)[names(Metadata.Mathys) == "AD"] <- "Dx"
names(Metadata.Mathys)[names(Metadata.Mathys) == "apoe_genotype"] <- "Apoe"
names(Metadata.Mathys)[names(Metadata.Mathys) == "braaksc"] <- "Braak"
names(Metadata.Mathys)[names(Metadata.Mathys) == "ceradsc"] <- "CERAD"

## recode some varibles
Metadata.Mathys$Age[Metadata.Mathys$Age == "90+"] <- 90
Metadata.Mathys$Age <- round(as.numeric(Metadata.Mathys$Age), 1)

Metadata.Mathys %<>% 
  mutate(Dx = dplyr::recode(PathDx, "NO" = "Ctrl", "YES" = "AD")) %>%
  mutate(Sex = dplyr::recode(as.character(Sex), "0" = "F", "1" = "M")) %>%
  dplyr::select(-PathDx)

## scale the age
Metadata.Mathys$zAge <- scale(Metadata.Mathys$Age)


# add sample information to the metadata
Meta.cell <- Meta.cell %>% left_join(Metadata.Mathys, by = "individualID") %>%
  column_to_rownames("Cells")


## Add the metadata to the Seurat object
Seurat.Mathys <- AddMetaData(Seurat.Mathys, metadata = Meta.cell)

# remove some varibles
Seurat.Mathys$nFeature_decontX <- NULL
Seurat.Mathys$nCount_decontX <- NULL
Seurat.Mathys$nCount_originalexp  <- NULL
Seurat.Mathys$nFeature_originalexp <- NULL

## make the new varible cngeneson (scaled cell detection rate)
Seurat.Mathys$cngeneson <- scale(Seurat.Mathys$nFeature_RNA)

## Check the doublets
VlnPlot(Seurat.Mathys, features = c("nFeature_RNA", "percent.mt"),
        group.by = "DFClass", pt.size = 0)

setwd(".../Step9.IntegrationOfPublicData/Mathys")

## subset the singlet
Seurat.Mathys <- subset(Seurat.Mathys, DFClass == "Singlet")
saveRDS(Seurat.Mathys, "Seurat.Mathys.Single.rds")


