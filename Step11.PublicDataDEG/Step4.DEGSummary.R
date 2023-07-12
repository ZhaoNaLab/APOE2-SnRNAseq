library("tidyverse")
library("data.table")
library("openxlsx")
library("magrittr")

Root.Dir <- ".../Step11.PublicDataDEG"
CellType <- c("Ex", "In", "Olig", "OPC", "Ast", "Mic", "Vascular")
DEG.comp <- list()

for (i in 1:length(CellType)){
  setwd(paste0(Root.Dir, "/", CellType[i]))
  FileList.Comb.RefE2 <- list.files(pattern = "Comb.RefE2.*.csv")
  FileList.Comb.RefE2.df <- lapply(FileList.Comb.RefE2, read.csv, row.names = 1)
  FileList.Comb.RefE2.df <- do.call(rbind, FileList.Comb.RefE2.df)
  
  FileList.Comb.RefE3 <- list.files(pattern = "Comb.RefE3.*.csv")
  FileList.Comb.RefE3.df <- lapply(FileList.Comb.RefE3, read.csv, row.names = 1)
  FileList.Comb.RefE3.df <- do.call(rbind, FileList.Comb.RefE3.df)
  
  FileList.Comb.RefE4 <- list.files(pattern = "Comb.RefE4.*.csv")
  FileList.Comb.RefE4.df <- lapply(FileList.Comb.RefE4, read.csv, row.names = 1)
  FileList.Comb.RefE4.df <- do.call(rbind, FileList.Comb.RefE4.df)
  
  FileList.Ctrl.RefE2 <- list.files(pattern = "Ctrl.*.RefE2.csv")
  FileList.Ctrl.RefE2.df <- lapply(FileList.Ctrl.RefE2, read.csv, row.names = 1)
  FileList.Ctrl.RefE2.df <- do.call(rbind, FileList.Ctrl.RefE2.df)
  
  FileList.Ctrl.RefE3 <- list.files(pattern = "Ctrl.*.RefE3.csv")
  FileList.Ctrl.RefE3.df <- lapply(FileList.Ctrl.RefE3, read.csv, row.names = 1)
  FileList.Ctrl.RefE3.df <- do.call(rbind, FileList.Ctrl.RefE3.df)
  
  FileList.Ctrl.RefE4 <- list.files(pattern = "Ctrl.*.RefE4.csv")
  FileList.Ctrl.RefE4.df <- lapply(FileList.Ctrl.RefE4, read.csv, row.names = 1)
  FileList.Ctrl.RefE4.df <- do.call(rbind, FileList.Ctrl.RefE4.df)
  
  FileList.AD.RefE2 <- list.files(pattern = "AD.*.RefE2.csv")
  FileList.AD.RefE2.df <- lapply(FileList.AD.RefE2, read.csv, row.names = 1)
  FileList.AD.RefE2.df <- do.call(rbind, FileList.AD.RefE2.df)
  
  FileList.AD.RefE3 <- list.files(pattern = "AD.*.RefE3.csv")
  FileList.AD.RefE3.df <- lapply(FileList.AD.RefE3, read.csv, row.names = 1)
  FileList.AD.RefE3.df <- do.call(rbind, FileList.AD.RefE3.df)
  
  FileList.AD.RefE4 <- list.files(pattern = "AD.*.RefE4.csv")
  FileList.AD.RefE4.df <- lapply(FileList.AD.RefE4, read.csv, row.names = 1)
  FileList.AD.RefE4.df <- do.call(rbind, FileList.AD.RefE4.df)
  
  FileList.E2 <- list.files(pattern = "^E2.*.csv")
  FileList.E2.df <- lapply(FileList.E2, read.csv, row.names = 1)
  FileList.E2.df <- do.call(rbind, FileList.E2.df)
  
  FileList.E3 <- list.files(pattern = "^E3.*.csv")
  FileList.E3.df <- lapply(FileList.E3, read.csv, row.names = 1)
  FileList.E3.df <- do.call(rbind, FileList.E3.df)
  
  FileList.E4 <- list.files(pattern = "^E4.*.csv")
  FileList.E4.df <- lapply(FileList.E4, read.csv, row.names = 1)
  FileList.E4.df <- do.call(rbind, FileList.E4.df)
  
  FileListCombined <- list(FileList.Comb.RefE2.df, FileList.Comb.RefE3.df, FileList.Comb.RefE4.df,
                           FileList.Ctrl.RefE2.df, FileList.Ctrl.RefE3.df, FileList.Ctrl.RefE4.df,
                           FileList.AD.RefE2.df, FileList.AD.RefE3.df, FileList.AD.RefE4.df,
                           FileList.E2.df, FileList.E3.df, FileList.E4.df)
  
  FileNames <- c("Comb.RefE2", "Comb.RefE3", "Comb.RefE4",
                 "Ctrl.RefE2", "Ctrl.RefE3", "Ctrl.RefE4",
                 "AD.RefE2", "AD.RefE3", "AD.RefE4",
                 "E2", "E3", "E4")
  
  tempDf <- list()
  
  for (j in 1:length(FileNames)){
    DEG.file <- FileListCombined[[j]]
    setDT(DEG.file)[, fdr:= p.adjust(`Pr..Chisq.`, 'fdr'), by = "contrast"]
    DEG.file$DEG <- NA
    DEG.file$DEG[DEG.file$fdr < 0.05 & DEG.file$coef >= 0.1] <- "up"
    DEG.file$DEG[DEG.file$fdr < 0.05 & DEG.file$coef <= -0.1] <- "down"
    DEG.file <- arrange(DEG.file, fdr)
    DEG.file$Group <- FileNames[j]
    tempDf[[j]] <- DEG.file
  }
  
  DEG.comp[[i]] <- do.call(rbind, tempDf) %>% unique()
}

names(DEG.comp) <- CellType

## Make sure that every cell types have all the comparisions
lapply(DEG.comp, function(x){table(x$Group, x$contrast)})

## Save the combined DEG results
setwd(Root.Dir)
saveRDS(DEG.comp, "DEG.AllinOne.rds")

### Save the files for uploading when publishing the manuscript
### Omit Sex effect in this study 
for (i in 1:length(DEG.comp)){
  Dat <- DEG.comp[[i]] %>%
    filter(!(contrast == "SexF")) %>%
    filter(!(contrast == "DxAD" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(contrast = case_when(str_detect(Group, "RefE2") & contrast == "ApoeE3"  ~ "E3vsE2",
                                str_detect(Group, "RefE2") & contrast == "ApoeE4"  ~ "E4vsE2",
                                str_detect(Group, "RefE3") & contrast == "ApoeE2"  ~ "E2vsE3",
                                str_detect(Group, "RefE3") & contrast == "ApoeE4"  ~ "E4vsE3",
                                str_detect(Group, "RefE4") & contrast == "ApoeE2"  ~ "E2vsE4",
                                str_detect(Group, "RefE4") & contrast == "ApoeE3"  ~ "E3vsE4",
                                str_detect(contrast, "AD")  ~ "ADvsControl",
                                TRUE ~ contrast)) %>%
    filter(!(contrast %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "Comb") ~ "Combined", 
                             str_detect(Group, "Ctrl") ~ "Control",
                             str_detect(Group, "AD") ~ "AD",
                             TRUE ~ Group))
  
  Dat.split <- split(Dat, list(Dat$Group,Dat$contrast), drop=TRUE)
  Dat.split <- Dat.split[c("Combined.ADvsControl", "Combined.E2vsE3","Combined.E4vsE3", 
                           "Combined.E2vsE4", "Control.E2vsE3", "Control.E4vsE3" , 
                           "Control.E2vsE4", "AD.E2vsE3", "AD.E4vsE3", "AD.E2vsE4", 
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


## Summarize the DEG number and plot the results
DEG.num.Sum <- list()

for (i in 1:length(DEG.comp)){
  Dat <- setDT(DEG.comp[[i]])[, .(Number = .N), 
                              by = c("Group", "contrast", "DEG")] %>%
    drop_na() %>%
    filter(!(contrast == "SexF" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    filter(!(contrast == "DxAD" & Group %in% c("Comb.RefE2", "Comb.RefE4"))) %>%
    mutate(Comparison = case_when(str_detect(Group, "RefE2") & contrast == "ApoeE3"  ~ "E3vsE2",
                                  str_detect(Group, "RefE2") & contrast == "ApoeE4"  ~ "E4vsE2",
                                  str_detect(Group, "RefE3") & contrast == "ApoeE2"  ~ "E2vsE3",
                                  str_detect(Group, "RefE3") & contrast == "ApoeE4"  ~ "E4vsE3",
                                  str_detect(Group, "RefE4") & contrast == "ApoeE2"  ~ "E2vsE4",
                                  str_detect(Group, "RefE4") & contrast == "ApoeE3"  ~ "E3vsE4",
                                  str_detect(contrast, "AD")  ~ "ADvsControl",
                                  str_detect(contrast, "SexF")  ~ "FemalevsMale",
                                  TRUE ~ contrast)) %>%
    filter(!(Comparison %in% c("E3vsE2", "E4vsE2", "E3vsE4"))) %>%
    mutate(Group = case_when(str_detect(Group, "AD")  ~ "AD",
                             str_detect(Group, "Comb")  ~ "Combined",
                             str_detect(Group, "Ctrl")  ~ "Control",
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

