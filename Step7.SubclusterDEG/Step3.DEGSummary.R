library("tidyverse")
library("data.table")
library("openxlsx")
library("magrittr")

Root.Dir <- ".../Step7.SubclusterDEG/DEGData"
setwd(Root.Dir)

CellType <- list.files()
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
setwd(".../Step7.SubclusterDEG/")
saveRDS(DEG.comp, "DEG.AllinOne.rds")

### Save the files for uploading when publishing the manuscript
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
                                str_detect(contrast, "AD")  ~ "ADvsCtrl",
                                TRUE ~ contrast)) %>%
    filter(!(contrast %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    filter(!(contrast == "ADvsCtrl" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "Comb") ~ "Combined", 
                             str_detect(Group, "Ctrl") ~ "Ctrl",
                             str_detect(Group, "AD") ~ "AD",
                             TRUE ~ Group))
  
  Dat.split <- split(Dat, list(Dat$Group,Dat$contrast), drop=TRUE)
  Dat.split <- Dat.split[intersect(c("Combined.ADvsCtrl", "Combined.E2vsE3","Combined.E4vsE3", "Combined.E2vsE4",
                                     "Ctrl.E2vsE3", "Ctrl.E4vsE3" , "Ctrl.E2vsE4", 
                                     "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
                                     "E2.ADvsCtrl", "E3.ADvsCtrl", "E4.ADvsCtrl"), 
                                   names(Dat.split))]
  
 
  ## Check the dimentions of the datasets
  lapply(Dat.split, dim)
  
  ## Put all the dataframes into an Excel workbook
  work_book <- createWorkbook()
  for (j in 1:length(Dat.split)){
    addWorksheet(work_book, names(Dat.split)[j])
    writeData(work_book, names(Dat.split)[j], Dat.split[[j]], startCol = 1, startRow = 1)
  }
  
  saveWorkbook(work_book,
               file = paste0("SuppementaryDataSheet/", names(DEG.comp)[i], ".MAST.DEG.xlsx"),
               overwrite = TRUE)
}


## Summarize the DEG number and plot the results
DEG.num.Sum <- list()

for (i in 1:length(DEG.comp)){
  Dat <- setDT(DEG.comp[[i]])[, .(Number = .N), 
                              by = c("Group", "contrast", "DEG")] %>%
    drop_na() %>%
    filter(!(contrast == "sexfemale" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(Comparison = case_when(str_detect(Group, "RefE2") & contrast == "apoeE3"  ~ "E3vsE2",
                                  str_detect(Group, "RefE2") & contrast == "apoeE4"  ~ "E4vsE2",
                                  str_detect(Group, "RefE3") & contrast == "apoeE2"  ~ "E2vsE3",
                                  str_detect(Group, "RefE3") & contrast == "apoeE4"  ~ "E4vsE3",
                                  str_detect(Group, "RefE4") & contrast == "apoeE2"  ~ "E2vsE4",
                                  str_detect(Group, "RefE4") & contrast == "apoeE3"  ~ "E3vsE4",
                                  str_detect(contrast, "AD")  ~ "ADvsCtrl",
                                  str_detect(contrast, "sexfemale")  ~ "FemalevsMale",
                                  TRUE ~ contrast)) %>%
    filter(!(Comparison %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    filter(!(Comparison == "ADvsCtrl" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "AD")  ~ "AD",
                             str_detect(Group, "Comb")  ~ "Combined",
                             str_detect(Group, "Ctrl")  ~ "Ctrl",
                             TRUE ~ Group))
  
  DEG.num.Sum[[i]] <- Dat 
  names(DEG.num.Sum)[i] <- names(DEG.comp)[i]
}

## Save the data: put all data into the excel files
## Put all the dataframes into an Excel workbook
work_book <- createWorkbook()
for (i in 1:length(DEG.num.Sum)){
  addWorksheet(work_book, sheetName = names(DEG.num.Sum)[i])
  writeData(work_book, names(DEG.num.Sum)[i], DEG.num.Sum[[i]], startCol = 1, startRow = 1)
}

saveWorkbook(work_book, "SuppementaryDataSheet/DEG.Number.Summary.xlsx", overwrite = TRUE)












