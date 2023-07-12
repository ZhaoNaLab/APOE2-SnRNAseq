library("Seurat")
library("tidyverse")

setwd(".../Step9.IntegrationOfPublicData/Blanchard/SeuratObj.Ind")
SeuratObjs <- list.files()

# load the files
Seurat.list <- lapply(SeuratObjs, readRDS)

# Merge into one single Seurat object
Seurat.Blanchard <- merge(x = Seurat.list[[1]], y = unlist(lapply(2:length(Seurat.list), function(x){Seurat.list[[x]]})))
rm(Seurat.list)

# add the individual information to the metadata
Meta.cell <- Seurat.Blanchard@meta.data %>% 
  rownames_to_column("Cells") %>%
  dplyr::select(-ends_with("y")) %>%
  dplyr::select(-c(AD, msex))


colnames(Meta.cell) <- gsub("\\.x", "",colnames(Meta.cell))


# Rename some varible
names(Meta.cell)[names(Meta.cell) == "age_death"] <- "Age"
names(Meta.cell)[names(Meta.cell) == "pmi"] <- "PMI"
names(Meta.cell)[names(Meta.cell) == "AD"] <- "Dx"
names(Meta.cell)[names(Meta.cell) == "apoe_genotype"] <- "Apoe"
names(Meta.cell)[names(Meta.cell) == "braaksc"] <- "Braak"
names(Meta.cell)[names(Meta.cell) == "ceradsc"] <- "CERAD"

## recode some varibles

# add sample information to the metadata
Meta.cell %<>% column_to_rownames("Cells") %>%
  subset(., select = which(!duplicated(names(.))))


## Add the metadata to the Seurat object
Seurat.Blanchard <- AddMetaData(Seurat.Blanchard, metadata = Meta.cell)

# remove some varibles
Seurat.Blanchard$nFeature_decontX <- NULL
Seurat.Blanchard$nCount_decontX <- NULL
Seurat.Blanchard$nCount_originalexp  <- NULL
Seurat.Blanchard$nFeature_originalexp <- NULL

## make the new varible cngeneson (scaled cell detection rate)
Seurat.Blanchard$cngeneson <- scale(Seurat.Blanchard$nFeature_RNA)

## Check the doublets
VlnPlot(Seurat.Blanchard, features = c("nCount_RNA", "percent.mt"),
        group.by = "DFClass", pt.size = 0)

setwd(".../Step9.IntegrationOfPublicData/Blanchard")

## subset the singlet
Seurat.Blanchard <- subset(Seurat.Blanchard, DFClass == "Singlet")
saveRDS(Seurat.Blanchard, "Seurat.Blanchard.Single.rds")



