library("tidyverse")
library("clusterProfiler")
library('data.table')
library("BiocParallel")
library("openxlsx")

# prepare from GMT files download from gprofile
ReadGMT_GSEA <- function (gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  fetchName <- function(x){
    unlist(x)[2]
  }
  names(pathways) <- sapply(pathwayLines, fetchName)
  
  temp <- list()
  for (i in 1:length(pathways)){
    temp[[i]] <- as.data.frame(cbind(GO_ID = names(pathways)[i],
                                     gene = pathways [[i]]))
  }
  gmt <- do.call(rbind, temp)
  return(gmt)
}

# load the GMT files
# test with GSEA function 
setwd("/research/labs/moleneurosci/bug/data/ZH_processing_data/SnRNA/GMT/gProfiler")
BPGMT <- ReadGMT_GSEA('hsapiens.GO_BP.name.gmt')
KeggGMT <- read.csv('KEGGMig.csv', row.names = 1)


# load the DEG datasets
setwd("path/to/project")

# load the DEG files
DEGFileCombind <- readRDS("Step7.SubclusterDEG/DEG.AllinOne.rds")
SubCluser <- names(DEGFileCombind)

# loop the pathway analysis
for (i in 1:length(SubCluser)){
  # what is the major cluster
  MajorCluster <-  gsub("[^a-zA-Z]", "", SubCluser[i])
  SubType <- SubCluser[i]
  
  dat <- DEGFileCombind[[i]] %>% as.data.frame() %>%
    dplyr::filter(!(contrast == "sexfemale")) %>%
    dplyr::filter(!(contrast == "bgAD" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(contrast = case_when(str_detect(Group, "RefE2") & contrast == "apoeE3"  ~ "E3vsE2",
                                str_detect(Group, "RefE2") & contrast == "apoeE4"  ~ "E4vsE2",
                                str_detect(Group, "RefE3") & contrast == "apoeE2"  ~ "E2vsE3",
                                str_detect(Group, "RefE3") & contrast == "apoeE4"  ~ "E4vsE3",
                                str_detect(Group, "RefE4") & contrast == "apoeE2"  ~ "E2vsE4",
                                str_detect(Group, "RefE4") & contrast == "apoeE3"  ~ "E3vsE4",
                                str_detect(contrast, "AD")   ~ "ADvsControl",
                                TRUE ~ contrast)) %>%
    dplyr::filter(!(contrast %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "Comb") ~ "Combined",
                             str_detect(Group, "Ctrl") ~ "Control",
                             str_detect(Group, "AD") ~ "AD",
                             TRUE ~ Group))

  Dat.split <- split(dat, list(dat$Group,dat$contrast), drop=TRUE)
  Dat.split <- Dat.split[intersect(c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3",
                                     "Combined.E2vsE4", "Control.E2vsE3", "Control.E4vsE3" ,
                                     "Control.E2vsE4", "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4",
                                     "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl"),
                                   names(Dat.split))]


  for (j in 1:length(Dat.split)){
    tryCatch({
      Dat <- Dat.split[[j]]
      Dat <- setorder(Dat, -z)
      Dat$index <- 1:nrow(Dat)
      DEGList <- with(Dat[, c('primerid', 'z')], setNames(z, primerid))
      DEGList <- DEGList[!is.na(DEGList)]
      DEGList <- DEGList[!is.infinite(DEGList)]
      
      dir.create(paste0(".../Step8.SubclusterGSEA/GSEA.significant/", MajorCluster, "/", SubType), 
                 recursive = T, 
                 showWarnings = F)
      
      Dir <- paste0(".../Step8.SubclusterGSEA/GSEA.significant/", MajorCluster, "/", SubType) 
      
      # perform GSEA analysis
      setwd(Dir)
      set.seed(2022110)
      GSEA.BP <- GSEA(DEGList, TERM2GENE =BPGMT, verbose = FALSE,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      BPPARAM=MulticoreParam(workers = 4))
      GSEA.BP.Table <- GSEA.BP@result
      
      saveRDS(GSEA.BP, paste0("BP.", names(Dat.split)[j], ".rds"))
      write.csv(GSEA.BP.Table, paste0("BP.Table.", names(Dat.split)[j], ".csv"))
      
      set.seed(2022110)
      GSEA.KEGG <- GSEA(DEGList, TERM2GENE = KeggGMT, verbose = FALSE,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        BPPARAM=MulticoreParam(workers = 4))
      GSEA.KEGG.Table <- GSEA.KEGG@result
      
      saveRDS(GSEA.KEGG, paste0("KEGG.", names(Dat.split)[j], ".rds"))
      write.csv(GSEA.KEGG.Table, paste0("KEGG.Table.", names(Dat.split)[j], ".csv"))
    }, error = function(x){})
  }
}

### Combine the GSEA results into excel files for publications
## Load the files
RootDir <- "..."
SubCluser <- names(DEGFileCombind)

## BP
for (i in 1:length(SubCluser)){
  MajorCluster <-  gsub("[^a-zA-Z]", "", SubCluser[i])
  SubType <- SubCluser[i]
  
  dir.create(paste0(".../Step8.SubclusterGSEA/GSEA.significant/", MajorCluster, "/", SubType),
             showWarnings = FALSE)
  Dir <- paste0(".../Step8.SubclusterGSEA/GSEA.significant/", MajorCluster, "/", SubType) 
  setwd(Dir)
  
  BP.files <- list.files(pattern = "*BP.Table*")
  BP.df <- lapply(BP.files, read.csv, row.names = 1)
  names(BP.df) <- gsub("BP.Table.|.csv", "", BP.files)
  
  BP.df <- BP.df[intersect(c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", "Combined.E2vsE4",
                             "Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                             "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                             "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl"), 
                           names(BP.df))]
  
  ## Put all the dataframes into an Excel workbook
  work_book <- createWorkbook()
  for (j in 1:length(BP.df)){
    addWorksheet(work_book, sheetName = names(BP.df)[j])
    writeData(work_book, names(BP.df)[j], BP.df[[j]], startCol = 1, startRow = 1)
  }
  
  dir.create(paste0(".../Step8.SubclusterGSEA/SuppementaryDataSheet/GSEA.significant/", MajorCluster))
  setwd(paste0(".../Step8.SubclusterGSEA/SuppementaryDataSheet/GSEA.significant/", MajorCluster))
  saveWorkbook(work_book,
               file = paste0(SubType, ".BP.xlsx"),
               overwrite = TRUE)
}

## KEGG
for (i in 1:length(SubCluser)){
  MajorCluster <-  gsub("[^a-zA-Z]", "", SubCluser[i])
  SubType <- SubCluser[i]
  
  Dir <- paste0(".../Step8.SubclusterGSEA/GSEA.significant/", MajorCluster, "/", SubType) 
  setwd(Dir)
  
  KEGG.files <- list.files(pattern = "*KEGG.Table*")
  KEGG.df <- lapply(KEGG.files, read.csv, row.names = 1)
  names(KEGG.df) <- gsub("KEGG.Table.|.csv", "", KEGG.files)
  
  KEGG.df <- KEGG.df[intersect(c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", "Combined.E2vsE4",
                                 "Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                                 "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                                 "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl"), 
                               names(KEGG.df))]
  
  ## Put all the dataframes into an Excel workbook
  work_book <- createWorkbook()
  for (j in 1:length(KEGG.df)){
    addWorksheet(work_book, sheetName = names(KEGG.df)[j])
    writeData(work_book, names(KEGG.df)[j], KEGG.df[[j]], startCol = 1, startRow = 1)
  }
  
  setwd(paste0(".../Step8.SubclusterGSEA/SuppementaryDataSheet/GSEA.significant/", MajorCluster))
  saveWorkbook(work_book,
               file = paste0(SubType, ".KEGG.xlsx"),
               overwrite = TRUE)
}

###################################### Do not use the thresphold #######################
for (i in 1:length(SubCluser)){
  # what is the major cluster
  MajorCluster <-  gsub("[^a-zA-Z]", "", SubCluser[i])
  SubType <- SubCluser[i]
  
  dat <- DEGFileCombind[[i]] %>% as.data.frame() %>% 
    dplyr::filter(!(contrast == "sexfemale")) %>%
    dplyr::filter(!(contrast == "bgAD" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(contrast = case_when(str_detect(Group, "RefE2") & contrast == "apoeE3"  ~ "E3vsE2",
                                str_detect(Group, "RefE2") & contrast == "apoeE4"  ~ "E4vsE2",
                                str_detect(Group, "RefE3") & contrast == "apoeE2"  ~ "E2vsE3",
                                str_detect(Group, "RefE3") & contrast == "apoeE4"  ~ "E4vsE3",
                                str_detect(Group, "RefE4") & contrast == "apoeE2"  ~ "E2vsE4",
                                str_detect(Group, "RefE4") & contrast == "apoeE3"  ~ "E3vsE4",
                                str_detect(contrast, "AD")   ~ "ADvsControl",
                                TRUE ~ contrast)) %>%
    dplyr::filter(!(contrast %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "Comb") ~ "Combined", 
                             str_detect(Group, "Ctrl") ~ "Control",
                             str_detect(Group, "AD") ~ "AD",
                             TRUE ~ Group))
  
  Dat.split <- split(dat, list(dat$Group,dat$contrast), drop=TRUE)
  Dat.split <- Dat.split[intersect(c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", 
                                     "Combined.E2vsE4", "Control.E2vsE3", "Control.E4vsE3" , 
                                     "Control.E2vsE4", "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                                     "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl"),
                                   names(Dat.split))]
  
  
  for (j in 1:length(Dat.split)){
    tryCatch({
      Dat <- Dat.split[[j]]
      Dat <- setorder(Dat, -z)
      Dat$index <- 1:nrow(Dat)
      DEGList <- with(Dat[, c('primerid', 'z')], setNames(z, primerid))
      DEGList <- DEGList[!is.na(DEGList)]
      DEGList <- DEGList[!is.infinite(DEGList)]
      
      dir.create(paste0(".../Step8.SubclusterGSEA/GSEA.All/", MajorCluster, "/", SubType), 
                 recursive = T, 
                 showWarnings = F)
      
      Dir <- paste0(".../Step8.SubclusterGSEA/GSEA.All/", MajorCluster, "/", SubType) 
      
      # perform GSEA analysis
      setwd(Dir)
      set.seed(2022110)
      GSEA.BP <- GSEA(DEGList, TERM2GENE =BPGMT, verbose = FALSE,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1.0,
                      BPPARAM=MulticoreParam(workers = 4))
      GSEA.BP.Table <- GSEA.BP@result
      
      saveRDS(GSEA.BP, paste0("BP.", names(Dat.split)[j], ".rds"))
      write.csv(GSEA.BP.Table, paste0("BP.Table.", names(Dat.split)[j], ".csv"))
      
      set.seed(2022110)
      GSEA.KEGG <- GSEA(DEGList, TERM2GENE = KeggGMT, verbose = FALSE,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1.0,
                        BPPARAM=MulticoreParam(workers = 4))
      GSEA.KEGG.Table <- GSEA.KEGG@result
      
      saveRDS(GSEA.KEGG, paste0("KEGG.", names(Dat.split)[j], ".rds"))
      write.csv(GSEA.KEGG.Table, paste0("KEGG.Table.", names(Dat.split)[j], ".csv"))
    }, error = function(x){})
  }
}


