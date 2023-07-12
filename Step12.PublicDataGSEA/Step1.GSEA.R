library("tidyverse")
library("clusterProfiler")
library('data.table')
library("BiocParallel")
library("openxlsx")

######################################### only significant ones ########################
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
setwd(".../GMT/gProfiler")
BPGMT <- ReadGMT_GSEA('hsapiens.GO_BP.name.gmt')
KeggGMT <- read.csv('KEGGMig.csv', row.names = 1)


# load the DEG datasets
setwd(".../Step11.PublicDataDEG")

# load the DEG files
DEGFileCombind <- readRDS("DEG.AllinOne.rds")
Cluser <- names(DEGFileCombind)

for (i in 1:length(Cluser)){
  dat <- DEGFileCombind[[i]] %>% as.data.frame() %>% 
    dplyr::filter(!(contrast == "SexF")) %>%
    dplyr::filter(!(contrast == "DxAD" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(contrast = case_when(str_detect(Group, "RefE2") & contrast == "ApoeE3"  ~ "E3vsE2",
                                str_detect(Group, "RefE2") & contrast == "ApoeE4"  ~ "E4vsE2",
                                str_detect(Group, "RefE3") & contrast == "ApoeE2"  ~ "E2vsE3",
                                str_detect(Group, "RefE3") & contrast == "ApoeE4"  ~ "E4vsE3",
                                str_detect(Group, "RefE4") & contrast == "ApoeE2"  ~ "E2vsE4",
                                str_detect(Group, "RefE4") & contrast == "ApoeE3"  ~ "E3vsE4",
                                str_detect(contrast, "AD")   ~ "ADvsControl",
                                TRUE ~ contrast)) %>%
    filter(!(contrast %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "Comb") ~ "Combined", 
                             str_detect(Group, "Ctrl") ~ "Control",
                             str_detect(Group, "AD") ~ "AD",
                             TRUE ~ Group))
  
  Dat.split <- split(dat, list(dat$Group,dat$contrast), drop=TRUE)
  Dat.split <- Dat.split[c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", 
                           "Combined.E2vsE4", "Control.E2vsE3", "Control.E4vsE3" , 
                           "Control.E2vsE4", "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                           "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  for (j in 1:length(Dat.split)){
    Dat <- Dat.split[[j]]
    Dat <- setorder(Dat, -z)
    Dat$index <- 1:nrow(Dat)
    DEGList <- with(Dat[, c('primerid', 'z')], setNames(z, primerid))
    DEGList <- DEGList[!is.na(DEGList)]
    dir.create(paste0(".../Step11.PublicDataDEG/GSEA.significant/", Cluser[i]), 
               recursive = T, 
               showWarnings = F)
    
    Dir <- paste0(".../Step11.PublicDataDEG/GSEA.significant/", Cluser[i]) 
    
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
  }
}

### Combine the GSEA results into excel files for publications
## Load the files
RootDir <- ".../Step11.PublicDataDEG"
Cluser <- c("Ex", "In", "Olig", "OPC", "Ast", "Mic", "Vascular")


## BP
for (i in 1:length(Cluser)){
  setwd(paste0(RootDir, "/Step11.PublicDataDEG/GSEA.significant/", Cluser[i]))
  
  BP.files <- list.files(pattern = "*BP.Table*")
  BP.df <- lapply(BP.files, read.csv, row.names = 1)
  names(BP.df) <- gsub("BP.Table.|.csv", "", BP.files)
  
  BP.df <- BP.df[c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", "Combined.E2vsE4",
                   "Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                   "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                   "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  ## Put all the dataframes into an Excel workbook
  work_book <- createWorkbook()
  for (j in 1:length(BP.df)){
    addWorksheet(work_book, sheetName = names(BP.df)[j])
    writeData(work_book, names(BP.df)[j], BP.df[[j]], startCol = 1, startRow = 1)
  }
  
  setwd(RootDir)
  saveWorkbook(work_book,
               file = paste0("Step11.PublicDataDEG/SuppementaryDataSheet/",Cluser[i], ".BP.Significant.xlsx"),
               overwrite = TRUE)
  
}

## KEGG
for (i in 1:length(Cluser)){
  setwd(paste0(RootDir, "/Step11.PublicDataDEG/GSEA.significant/", Cluser[i]))
  
  KEGG.files <- list.files(pattern = "*KEGG.Table*")
  KEGG.df <- lapply(KEGG.files, read.csv, row.names = 1)
  names(KEGG.df) <- gsub("KEGG.Table.|.csv", "", KEGG.files)
  
  KEGG.df <- KEGG.df[c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", "Combined.E2vsE4",
                       "Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                       "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                       "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  ## Put all the dataframes into an Excel workbook
  work_book <- createWorkbook()
  for (j in 1:length(KEGG.df)){
    addWorksheet(work_book, sheetName = names(KEGG.df)[j])
    writeData(work_book, names(KEGG.df)[j], KEGG.df[[j]], startCol = 1, startRow = 1)
  }
  
  setwd(RootDir)
  saveWorkbook(work_book,
               file = paste0("Step11.PublicDataDEG/SuppementaryDataSheet/",Cluser[i], ".KEGG.Significant.xlsx"),
               overwrite = TRUE)
  
}


######################################### All Pathway ########################
for (i in 1:length(Cluser)){
  dat <- DEGFileCombind[[i]] %>% as.data.frame() %>% 
    dplyr::filter(!(contrast == "SexF")) %>%
    dplyr::filter(!(contrast == "DxAD" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(contrast = case_when(str_detect(Group, "RefE2") & contrast == "ApoeE3"  ~ "E3vsE2",
                                str_detect(Group, "RefE2") & contrast == "ApoeE4"  ~ "E4vsE2",
                                str_detect(Group, "RefE3") & contrast == "ApoeE2"  ~ "E2vsE3",
                                str_detect(Group, "RefE3") & contrast == "ApoeE4"  ~ "E4vsE3",
                                str_detect(Group, "RefE4") & contrast == "ApoeE2"  ~ "E2vsE4",
                                str_detect(Group, "RefE4") & contrast == "ApoeE3"  ~ "E3vsE4",
                                str_detect(contrast, "AD")   ~ "ADvsControl",
                                TRUE ~ contrast)) %>%
    filter(!(contrast %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "Comb") ~ "Combined", 
                             str_detect(Group, "Ctrl") ~ "Control",
                             str_detect(Group, "AD") ~ "AD",
                             TRUE ~ Group))
  
  Dat.split <- split(dat, list(dat$Group,dat$contrast), drop=TRUE)
  Dat.split <- Dat.split[c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", 
                           "Combined.E2vsE4", "Control.E2vsE3", "Control.E4vsE3" , 
                           "Control.E2vsE4", "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                           "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  for (j in 1:length(Dat.split)){
    Dat <- Dat.split[[j]]
    Dat <- setorder(Dat, -z)
    Dat$index <- 1:nrow(Dat)
    DEGList <- with(Dat[, c('primerid', 'z')], setNames(z, primerid))
    DEGList <- DEGList[!is.na(DEGList)]
    dir.create(paste0(".../Step11.PublicDataDEG/GSEA.All/", Cluser[i]), 
               recursive = T, 
               showWarnings = F)
    
    Dir <- paste0(".../Step11.PublicDataDEG/GSEA.All/", Cluser[i]) 
    
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
  }
}

### Combine the GSEA results into excel files for publications
## BP
for (i in 1:length(Cluser)){
  setwd(paste0(RootDir, "/Step11.PublicDataDEG/GSEA.All/", Cluser[i]))
  
  BP.files <- list.files(pattern = "*BP.Table*")
  BP.df <- lapply(BP.files, read.csv, row.names = 1)
  names(BP.df) <- gsub("BP.Table.|.csv", "", BP.files)
  
  BP.df <- BP.df[c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", "Combined.E2vsE4",
                   "Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                   "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                   "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  ## Put all the dataframes into an Excel workbook
  work_book <- createWorkbook()
  for (j in 1:length(BP.df)){
    addWorksheet(work_book, sheetName = names(BP.df)[j])
    writeData(work_book, names(BP.df)[j], BP.df[[j]], startCol = 1, startRow = 1)
  }
  
  setwd(RootDir)
  saveWorkbook(work_book,
               file = paste0("Step11.PublicDataDEG/SuppementaryDataSheet/",Cluser[i], ".BP.All.xlsx"),
               overwrite = TRUE)
  
}

## KEGG
for (i in 1:length(Cluser)){
  setwd(paste0(RootDir, "/Step11.PublicDataDEG/GSEA.All/", Cluser[i]))
  
  KEGG.files <- list.files(pattern = "*KEGG.Table*")
  KEGG.df <- lapply(KEGG.files, read.csv, row.names = 1)
  names(KEGG.df) <- gsub("KEGG.Table.|.csv", "", KEGG.files)
  
  KEGG.df <- KEGG.df[c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", "Combined.E2vsE4",
                       "Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                       "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                       "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  ## Put all the dataframes into an Excel workbook
  work_book <- createWorkbook()
  for (j in 1:length(KEGG.df)){
    addWorksheet(work_book, sheetName = names(KEGG.df)[j])
    writeData(work_book, names(KEGG.df)[j], KEGG.df[[j]], startCol = 1, startRow = 1)
  }
  
  setwd(RootDir)
  saveWorkbook(work_book,
               file = paste0("Step11.PublicDataDEG/SuppementaryDataSheet/",Cluser[i], ".KEGG.All.xlsx"),
               overwrite = TRUE)
  
}


### Combine the GSEA results into excel files for publication
## Create a workbook for the final Excel file
final_workbook <- createWorkbook()

## BP
for (i in 1:length(Cluser)){
  setwd(paste0(RootDir, "/Step11.PublicDataDEG/GSEA.All/", Cluser[i]))
  
  BP.files <- list.files(pattern = "*BP.Table*")
  BP.df <- lapply(BP.files, read.csv, row.names = 1)
  names(BP.df) <- gsub("BP.Table.|.csv", "", BP.files)
  
  BP.df <- BP.df[c("Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                   "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                   "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  names(BP.df) <- c("E2 vs E3 in Ctrl", "E4 vs E3 in Ctrl",
                    "E2 vs E4 in Ctrl", "E2 vs E3 in AD",
                    "E4 vs E3 in AD", "E2 vs E4 in AD",
                    "AD vs Ctrl in E2", "AD vs Ctrl in E3",
                    "AD vs Ctrl in E4")
  
  ## Add a new sheet for the current Cluser
  addWorksheet(final_workbook, sheetName = Cluser[i])
  
  ## Write the data horizontally with one blank column in between
  start_col <- 1
  for (j in 1:length(BP.df)){
    # Write dataset name
    writeData(final_workbook, Cluser[i], data.frame(names(BP.df)[j]), startCol = start_col, startRow = 1, colNames = FALSE)
    # Write dataset content
    writeData(final_workbook, Cluser[i], BP.df[[j]], startCol = start_col, startRow = 2)
    start_col <- start_col + ncol(BP.df[[j]]) + 1
  }
}

## Save the final workbook with all Cluser and their data
setwd(RootDir)
saveWorkbook(final_workbook,
             file = "Step11.PublicDataDEG/SuppementaryDataSheet/All_Cluser_BP_All.xlsx",
             overwrite = TRUE)

