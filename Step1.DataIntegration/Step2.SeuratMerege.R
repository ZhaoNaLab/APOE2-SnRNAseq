# Import necessary libraries
library("Seurat")
library("tidyverse")

# Set working directory to the location where individual Seurat objects are saved
setwd(".../Step1.DataIntegration/SeuratObj.Ind")

# Create a list of Seurat object files in the current directory
SeuratObjs <- list.files()

# Load the Seurat objects into a list
Seurat.list <- lapply(SeuratObjs, readRDS)

# Merge all Seurat objects into one
merged.Seurat <- merge(x = Seurat.list[[1]], y = unlist(lapply(2:56, function(x){Seurat.list[[x]]})))
rm(Seurat.list) # Clear the Seurat.list from memory to save space

# Add individual information to the metadata
Meta.cell <- merged.Seurat@meta.data %>% 
  rownames_to_column("Cells") 

# Update orig.ident with simplified identifiers
Orig.Ident <- gsub("R", "", paste0("S", str_split_fixed(Meta.cell$Cells, "_", n = 3)[, 2]))
Meta.cell$orig.ident <- Orig.Ident

# Load metadata
Metadata <- read.csv(".../Metadata.Ind.csv")
Metadata$zAge <- scale(Metadata$Age) # Scale the age data

# Subset the merged.Seurat
merged.Seurat$orig.ident <- Orig.Ident
merged.Seurat <- subset(merged.Seurat, orig.ident %in% Metadata$orig.ident)

# Join the metadata with sample information
Meta.cell <- Meta.cell %>% left_join(Metadata, by = "orig.ident") %>%
  column_to_rownames("Cells")

# Add the updated metadata to the merged.Seurat object
merged.Seurat <- AddMetaData(merged.Seurat, metadata = Meta.cell)

# Remove unnecessary variables
merged.Seurat$nFeature_decontX <- NULL
merged.Seurat$nCount_decontX <- NULL

# Create a new variable, 'cngeneson', as a scaled version of the cell detection rate
merged.Seurat$cngeneson <- scale(merged.Seurat$nFeature_RNA)

# Reorder the metadata variables
merged.Seurat@meta.data <- merged.Seurat@meta.data %>% dplyr::select(c(1, 4:5, 2, 8, 9, 10, 7, 3, 6, 11:29))

# Change working directory to save the merged object
setwd(".../Step1.DataIntegration")

# Save the Seurat object containing both singlets and doublets
saveRDS(merged.Seurat, "merged.Seurat.rds")
