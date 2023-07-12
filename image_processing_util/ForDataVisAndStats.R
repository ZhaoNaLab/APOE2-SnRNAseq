library("tidyverse")

## load the plaque ROI data
setwd("/research/labs/moleneurosci/bug/m198507/Projects/snRNAseq/Step19.Validation/IFImageProcessing/Results")
Plaque.ROI.df <- read.csv("PlaqueROIInfo.csv")
Plaque.ROI.df <- Plaque.ROI.df %>%
  arrange(-Area)