library("Seurat")
library("harmony")
library("tidyverse")
library("magrittr")

setwd(".../Step9.IntegrationOfPublicData")

## Blanchard
Seurat.Blanchard <- readRDS("Blanchard/Seurat.Blanchard.Single.rds")
Seurat.Blanchard$Study <- "Blanchard et al"
Seurat.Blanchard <- subset(Seurat.Blanchard, orig.ident != "SeuratProject")

## Make sure the Dx and APOE genotype information are consistent with the original publication
Blanchard.Meta.Ind <- Seurat.Blanchard@meta.data %>%
  dplyr::select(orig.ident, projid, Apoe, Dx) %>%
  distinct()

table(Blanchard.Meta.Ind$Dx, Blanchard.Meta.Ind$Apoe)

## Mathys
Seurat.Mathys <- readRDS("Mathys/Seurat.Mathys.Single.rds")
Seurat.Mathys$Study <- "Mathys et al"

## Make sure the Dx and APOE genotype information are consistent with the original publication
Mathys.Meta.Ind <- Seurat.Mathys@meta.data %>%
  dplyr::select(orig.ident, projid, Apoe, Dx) %>%
  distinct()

table(Mathys.Meta.Ind$Dx)
table(Mathys.Meta.Ind$Dx, Mathys.Meta.Ind$Apoe)


## Sadick
Seurat.Sadick <- readRDS("Sadick/Seurat.Sadick.Single.rds")
Seurat.Sadick$Study <- "Sadick et al" 


## Make sure the Dx and APOE genotype information are consistent with the original publication
Sadick.Meta.Ind <- Seurat.Sadick@meta.data %>%
  dplyr::select(orig.ident, Apoe, Dx) %>%
  distinct()

table(Sadick.Meta.Ind$Dx)
table(Sadick.Meta.Ind$Dx, Sadick.Meta.Ind$Apoe)


## add additinal varibles to the metadata
Seurat.Sadick$individualID <- NA
Seurat.Sadick$projid <- NA

## Zhou
Seurat.Zhou <- readRDS("Zhou/Seurat.Zhou.Single.rds")
Seurat.Zhou$Study <- "Zhou et al" 
Seurat.Zhou@meta.data$Sex[Seurat.Zhou@meta.data$Sex == "female "] <- "F"
Seurat.Zhou@meta.data$Sex[Seurat.Zhou@meta.data$Sex == "male "] <- "M"
Seurat.Zhou$orig.ident <- Seurat.Zhou$individualID

## Make sure the Dx and APOE genotype information are consistent with the original publication
Zhou.Meta.Ind <- Seurat.Zhou@meta.data %>%
  dplyr::select(orig.ident, Apoe, Dx) %>%
  distinct()

table(Zhou.Meta.Ind$Dx)
table(Zhou.Meta.Ind$Dx, Zhou.Meta.Ind$Apoe)


## Exclude AD-R62H individuals from the downstream analysis
Seurat.Zhou <- subset(Seurat.Zhou, Dx != "AD-R62H")
Seurat.Zhou@meta.data$Dx[Seurat.Zhou@meta.data$Dx == "AD-CV"] <- "AD"


#### Get the common features
Common.features <- purrr::reduce(list(rownames(Seurat.Blanchard), rownames(Seurat.Mathys),
                                      rownames(Seurat.Sadick), rownames(Seurat.Zhou)),
                                 intersect)

### Get the common metadata varibles
Common.varibles <- purrr::reduce(list(names(Seurat.Blanchard@meta.data), names(Seurat.Mathys@meta.data),
                                      names(Seurat.Sadick@meta.data), names(Seurat.Zhou@meta.data)),
                                 intersect)


### Subset the Seurate objects based on common features
Seurat.Blanchard <- subset(Seurat.Blanchard, features = Common.features)
Seurat.Blanchard@meta.data <- Seurat.Blanchard@meta.data %>% dplyr::select(all_of(Common.varibles))

Seurat.Mathys <- subset(Seurat.Mathys, features = Common.features)
Seurat.Mathys@meta.data <- Seurat.Mathys@meta.data %>% dplyr::select(all_of(Common.varibles))

Seurat.Sadick <- subset(Seurat.Sadick, features = Common.features)
Seurat.Sadick@meta.data <- Seurat.Sadick@meta.data %>% dplyr::select(all_of(Common.varibles))

Seurat.Zhou <- subset(Seurat.Zhou, features = Common.features)
Seurat.Zhou@meta.data <- Seurat.Zhou@meta.data %>% dplyr::select(all_of(Common.varibles))


## merge the four Seuratobj
Seurat.external <- merge(x = Seurat.Blanchard, y = c(Seurat.Mathys, Seurat.Sadick, Seurat.Zhou))
Seurat.external$Age[Seurat.external$individualID == "R2880377"] <- 96


## Add new sample ID
Seurat.external@meta.data$NewID <- NA
Seurat.external@meta.data$NewID <- paste0(Seurat.external@meta.data$orig.ident, ".", Seurat.external@meta.data$Study)
Meta.external <- Seurat.external@meta.data 
Meta.external.ind <- Meta.external %>% dplyr::select(NewID, orig.ident, Sex,
                                                     Age, Dx, Apoe, Braak, CERAD,
                                                     zAge, Study, individualID) %>%
  remove_rownames() %>% unique()

## remove duplicated columns
Meta.external.cell <- Seurat.external@meta.data
Columnnames <- colnames(Meta.external.cell)
Columnnames <- gsub("\\.x|\\.y", "", Columnnames)
colnames(Meta.external.cell) <- Columnnames

# update the metadat
Meta.external.cell <- Meta.external.cell %>% subset(., select = which(!duplicated(names(.))))
# Seurat.external@meta.data <- Seurat.external@meta.data[, c(1:23, 34:38, 40:41, 43:48)]

## Select samples 
xtabs(~ Dx + Apoe + Study, data = Meta.external.ind)

# select one of the duplicated sample
# Remove subjects with No APOE genoype
Seurat.external.filter <- subset(Seurat.external,
                                 NewID %in%
                                   c("D1.2.Sadick et al",
                                      "D10.Sadick et al",
                                      "D11.Sadick et al",
                                      "D12.Sadick et al",
                                      "D13.Sadick et al",
                                      "D15.Sadick et al",
                                      "D16.Sadick et al",
                                      "D17.Sadick et al",
                                      "D2.Sadick et al",
                                      "D3.Sadick et al",
                                      "D4.Sadick et al",
                                      "D5.Sadick et al",
                                      "D6.Sadick et al",
                                      "D7.Sadick et al",
                                      "D8.Sadick et al",
                                      "D9.Sadick et al",
                                      "R1067972.Mathys et al",
                                      "R1154454.Blanchard et al",
                                      "R1375133.Zhou et al",
                                      "R1538032.Blanchard et al",
                                      "R1845714.Blanchard et al",
                                      "R1924801.Blanchard et al",
                                      "R2008064.Zhou et al",
                                      "R2081705.Mathys et al",
                                      "R2144127.Blanchard et al",
                                      "R2294544.Mathys et al",
                                      "R2474257.Blanchard et al",
                                      "R2735120.Zhou et al",
                                      "R2880377.Mathys et al",
                                      "R2895885.Mathys et al",
                                      "R3067449.Blanchard et al",
                                      "R3086211.Mathys et al",
                                      "R3177264.Mathys et al",
                                      "R3209518.Blanchard et al",
                                      "R3368249.Zhou et al",
                                      "R3405776.Mathys et al",
                                      "R3442506.Zhou et al",
                                      "R3744330.Blanchard et al",
                                      "R3840906.Blanchard et al",
                                      "R3884524.Mathys et al",
                                      "R3900996.Mathys et al",
                                      "R3997006.Blanchard et al",
                                      "R4012015.Mathys et al",
                                      "R4042599.Mathys et al",
                                      "R4146432.Mathys et al",
                                      "R4258320.Mathys et al",
                                      "R4262244.Blanchard et al",
                                      "R4379962.Blanchard et al",
                                      "R4415805.Mathys et al",
                                      "R4439627.Blanchard et al",
                                      "R4482444.Mathys et al",
                                      "R4567280.Mathys et al",
                                      "R4689636.Mathys et al",
                                      "R4728676.Mathys et al",
                                      "R4739508.Mathys et al",
                                      "R4841941.Zhou et al",
                                      "R5061712.Zhou et al",
                                      "R5138383.Zhou et al",
                                      "R5184427.Blanchard et al",
                                      "R5196723.Mathys et al",
                                      "R5385855.Blanchard et al",
                                      "R5394614.Zhou et al",
                                      "R5525186.Zhou et al",
                                      "R5741580.Blanchard et al",
                                      "R5816648.Mathys et al",
                                      "R5885245.Mathys et al",
                                      "R6086139.Blanchard et al",
                                      "R6114572.Mathys et al",
                                      "R6176158.Mathys et al",
                                      "R6231758.Zhou et al",
                                      "R6267541.Mathys et al",
                                      "R6411801.Zhou et al",
                                      "R6528749.Zhou et al",
                                      "R6636386.Blanchard et al",
                                      "R6665276.Blanchard et al",
                                      "R6911631.Blanchard et al",
                                      "R7039412.Blanchard et al",
                                      "R7066784.Blanchard et al",
                                      "R7160627.Mathys et al",
                                      "R7288382.Mathys et al",
                                      "R7423003.Zhou et al",
                                      "R7583108.Blanchard et al",
                                      "R7698313.Blanchard et al",
                                      "R7721691.Mathys et al",
                                      "R7770387.Blanchard et al",
                                      "R7912121.Blanchard et al",
                                      "R7944883.Blanchard et al",
                                      "R8399817.Zhou et al",
                                      "R8451530.Mathys et al",
                                      "R8598847.Mathys et al",
                                      "R8608442.Blanchard et al",
                                      "R8629052.Mathys et al",
                                      "R8725848.Mathys et al",
                                      "R8744945.Mathys et al",
                                      "R8983137.Blanchard et al",
                                      "R9094222.Mathys et al",
                                      "R9113571.Mathys et al",
                                      "R9176076.Blanchard et al",
                                      "R9255163.Zhou et al",
                                      "R9307768.Blanchard et al",
                                      "R9426782.Mathys et al",
                                      "R9453274.Zhou et al",
                                      "R9500594.Mathys et al",
                                      "R9551808.Mathys et al",
                                      "R9679238.Zhou et al"))


# Make new orig.ident
Meta.external.ind$orig.ident <- paste0("S", 1:nrow(Meta.external.ind))
Meta.filter <- Seurat.external.filter@meta.data %>% rownames_to_column("CellID") %>%
  dplyr::select(-orig.ident) %>%
  left_join(Meta.external.ind %>% dplyr::select(NewID, orig.ident), by = "NewID") %>%
  column_to_rownames("CellID") %>% dplyr::select(orig.ident, everything())

# update the metadata
Seurat.external.filter <- AddMetaData(Seurat.external.filter, Meta.filter)
Seurat.external.filter$NewID <- NULL
Meta.external.filer <- Seurat.external.filter@meta.data 
Meta.external.ind.filter <- Meta.external.filer %>% dplyr::select(orig.ident, Sex,
                                                     Age, Dx, Apoe, Braak, CERAD,
                                                     zAge, Study, individualID, projid) %>%
  remove_rownames() %>% unique()



# Save the data
saveRDS(Seurat.external.filter, ".../Step10.PublicDataIntegration/Seurat.external.filter.rds")
write.csv(Meta.external.ind.filter, ".../Step10.PublicDataIntegration/Meta.external.ind.filter.csv")















