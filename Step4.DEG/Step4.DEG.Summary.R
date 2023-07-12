library("tidyverse")
library("data.table")
library("openxlsx")
library("magrittr")

Root.Dir <- ".../Step4.DEG"
CellType <- c("Ex", "In", "Olig", "OPC", "Ast", "Mic", "Vascu")
DEG.comp <- list()

for (i in 1:length(CellType)){
  setwd(paste0(Root.Dir, "/", CellType[i]))
  FileList <- list.files(pattern = ".csv")
  
  tempDf <- list()
  
  for (j in 1:length(FileList)){
    DEG.file <- read.csv(FileList[j], row.names = 1)
    setDT(DEG.file)[, fdr:= p.adjust(`Pr..Chisq.`, 'fdr'), by = "contrast"]
    DEG.file$DEG <- NA
    DEG.file$DEG[DEG.file$fdr < 0.05 & DEG.file$coef >= 0.1] <- "up"
    DEG.file$DEG[DEG.file$fdr < 0.05 & DEG.file$coef <= -0.1] <- "down"
    DEG.file <- arrange(DEG.file, fdr)
    DEG.file$Group <- gsub(".csv", "", FileList[j])
    tempDf[[j]] <- DEG.file
  }
  
  DEG.comp[[i]] <- do.call(rbind, tempDf)
}

names(DEG.comp) <- CellType

## Make sure that every cell types have all the comparisions
lapply(DEG.comp, function(x){table(x$Group, x$contrast)})

## Save the combined DEG results
saveRDS(DEG.comp, "DEG.AllinOne.rds")

### Save the files for publishing
### Omit Sex effect in this study 
for (i in 1:length(DEG.comp)){
  Dat <- DEG.comp[[i]] %>%
    filter(!(contrast == "sexfemale")) %>%
    mutate(contrast = case_when(str_detect(Group, "RefE2") & contrast == "apoeE3"  ~ "E3vsE2",
                                str_detect(Group, "RefE2") & contrast == "apoeE4"  ~ "E4vsE2",
                                str_detect(Group, "RefE3") & contrast == "apoeE2"  ~ "E2vsE3",
                                str_detect(Group, "RefE3") & contrast == "apoeE4"  ~ "E4vsE3",
                                str_detect(Group, "RefE4") & contrast == "apoeE2"  ~ "E2vsE4",
                                str_detect(Group, "RefE4") & contrast == "apoeE3"  ~ "E3vsE4",
                                str_detect(contrast, "AD")  ~ "ADvsControl",
                                TRUE ~ contrast)) %>%
    filter(!(contrast %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    filter(!(contrast == "ADvsControl" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "Comb") ~ "Combined", 
                             str_detect(Group, "Ctrl") ~ "Control",
                             str_detect(Group, "AD") ~ "AD",
                             TRUE ~ Group))
  
  Dat.split <- split(Dat, list(Dat$Group,Dat$contrast), drop=TRUE)
  Dat.split <- Dat.split[c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", "Combined.E2vsE4",
                           "Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                           "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                           "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  ## Check the dimentions of the datasets
  lapply(Dat.split, dim)
  
  ## Put all the dataframes into an Excel workbook
  work_book <- createWorkbook()
  for (j in 1:length(Dat.split)){
    addWorksheet(work_book, sheetName = names(Dat.split)[j])
    writeData(work_book, names(Dat.split)[j], Dat.split[[j]], startCol = 1, startRow = 1)
  }
  
  saveWorkbook(work_book,
               file = paste0("SuppementaryDataSheet/", names(DEG.comp)[i], ".MAST.DEG.xlsx"),
               overwrite = TRUE)
}


# It is probably better to put all DEG result in one excel file?
### Save the files for uploading when publishing the manuscript
### Omit Sex effect in this study 

## Create a workbook for the final Excel file
final_workbook <- createWorkbook()

for (i in 1:length(DEG.comp)){
  Dat <- DEG.comp[[i]] %>%
    filter(!(contrast == "sexfemale")) %>%
    mutate(contrast = case_when(str_detect(Group, "RefE2") & contrast == "apoeE3"  ~ "E3vsE2",
                                str_detect(Group, "RefE2") & contrast == "apoeE4"  ~ "E4vsE2",
                                str_detect(Group, "RefE3") & contrast == "apoeE2"  ~ "E2vsE3",
                                str_detect(Group, "RefE3") & contrast == "apoeE4"  ~ "E4vsE3",
                                str_detect(Group, "RefE4") & contrast == "apoeE2"  ~ "E2vsE4",
                                str_detect(Group, "RefE4") & contrast == "apoeE3"  ~ "E3vsE4",
                                str_detect(contrast, "AD")  ~ "ADvsControl",
                                TRUE ~ contrast)) %>%
    filter(!(contrast %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    filter(!(contrast == "ADvsControl" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "Comb") ~ "Combined", 
                             str_detect(Group, "Ctrl") ~ "Control",
                             str_detect(Group, "AD") ~ "AD",
                             TRUE ~ Group))
  
  Dat.split <- split(Dat, list(Dat$Group, Dat$contrast), drop=TRUE)
  Dat.split <- Dat.split[c("Control.E2vsE3", "Control.E4vsE3" , "Control.E2vsE4", 
                           "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                           "E2.ADvsControl", "E3.ADvsControl", "E4.ADvsControl")]
  
  names(Dat.split) <- c("E2 vs E3 in Ctrl", "E4 vs E3 in Ctrl",
                        "E2 vs E4 in Ctrl", "E2 vs E3 in AD",
                        "E4 vs E3 in AD", "E2 vs E4 in AD",
                        "AD vs Ctrl in E2", "AD vs Ctrl in E3",
                        "AD vs Ctrl in E4")
  
  ## Add a new sheet for the current cell type
  addWorksheet(final_workbook, sheetName = names(DEG.comp)[i])
  
  ## Write the data horizontally with one blank column in between
  start_col <- 1
  for (j in 1:length(Dat.split)){
    # Write dataset name
    writeData(final_workbook, names(DEG.comp)[i], data.frame(names(Dat.split)[j]), startCol = start_col, startRow = 1, colNames = FALSE)
    # Write dataset content
    writeData(final_workbook, names(DEG.comp)[i], Dat.split[[j]], startCol = start_col, startRow = 2)
    start_col <- start_col + ncol(Dat.split[[j]]) + 1
  }
}

## Save the final workbook with all cell types and their data
saveWorkbook(final_workbook,
             file = "SuppementaryDataSheet/All_CellTypes_MAST.DEG.xlsx",
             overwrite = TRUE)
