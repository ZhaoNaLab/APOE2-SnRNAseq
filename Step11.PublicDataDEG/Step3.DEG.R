# https://mllg.github.io/batchtools/articles/batchtools.html
library("batchtools")
setwd(".../Step11.PublicDataDEG")

reg <- makeExperimentRegistry("MAST.DEG.Registry",
                              packages = c("MAST", "data.table", "Seurat", "lme4"),
                              seed = 22112019)

reg <- loadRegistry("MAST.DEG.Registry", writeable = TRUE)


# Step1: function to subset data
# sca file has been filtered and were siplit into chunks
scaObjs <- readRDS("scaObjs.rds")

# define instance (data preparation)
subsample <- function(data, job, strata, APOEgenotype, Ref, DX, i, j) { # i cell types, j Chun
  Cluster <- names(data)[i]
  if (strata == "comb"){
    data = data[[i]]
    chunk_number <- 20   
    Chunks <- split(rownames(data),             # Applying split() function
                    cut(seq_along(rownames(data)),
                        chunk_number,
                        labels = FALSE))
    
    data = data[unique(c("RPL29", Chunks[[j]])), ]
    
  } else if (strata == "Apoe"){
    data = subset(data[[i]], Apoe == APOEgenotype)
    chunk_number <- 20   
    Chunks <- split(rownames(data),             # Applying split() function
                    cut(seq_along(rownames(data)),
                        chunk_number,
                        labels = FALSE))
    
    data = data[unique(c("RPL29", Chunks[[j]])), ]
  } else {
    data = subset(data[[i]], Dx == DX)
    chunk_number <- 20   
    Chunks <- split(rownames(data),             # Applying split() function
                    cut(seq_along(rownames(data)),
                        chunk_number,
                        labels = FALSE))
    
    data = data[unique(c("RPL29", Chunks[[j]])), ]
  } 
  
  list(data = data, 
       strata = strata, 
       APOEgenotype = APOEgenotype, 
       Ref = Ref, 
       DX = DX, 
       Cluster = Cluster,
       Chunk = j)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "GetSCAs", data = scaObjs, fun = subsample) 

# create an algorithm which applies a support vector machine
MASTDEG <- function(data, job, instance) {# instance is inherited from the last step
  setwd(".../Step11.PublicDataDEG")
  if (instance$strata == "comb"){
    if(instance$Ref == "E2"){
      data = instance$data
      Cluster = instance$Cluster
      Chunk = instance$Chunk
      data@colData$Apoe <- relevel(as.factor(data@colData$Apoe), ref = "E2")
      data@colData$Dx <- relevel(as.factor(data@colData$Dx), ref = "Ctrl")
      data@colData$Sex <- relevel(as.factor(data@colData$Sex), ref = "M")
      zlmCond <- zlm(~ Dx + Apoe + Sex + zAge + cngeneson + percent.mt + decontX_contamination + (1|Study/orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     force = TRUE,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('DxAD', 'SexF', 'ApoeE3', 'ApoeE4'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluster, '/Comb.RefE2.', Chunk, '.csv'))
    } else if (instance$Ref == "E3") {
      data = instance$data
      Cluster = instance$Cluster
      Chunk = instance$Chunk
      data@colData$Apoe <- relevel(as.factor(data@colData$Apoe), ref = "E3")
      data@colData$Dx <- relevel(as.factor(data@colData$Dx), ref = "Ctrl")
      data@colData$Sex <- relevel(as.factor(data@colData$Sex), ref = "M")
      zlmCond <- zlm(~ Dx + Apoe + Sex + zAge + cngeneson + percent.mt + decontX_contamination + (1|Study/orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     force = TRUE,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('DxAD', 'SexF', 'ApoeE2', 'ApoeE4'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluster, '/Comb.RefE3.', Chunk, '.csv'))
    } else {
      data = instance$data
      Cluster = instance$Cluster
      Chunk = instance$Chunk
      data@colData$Apoe <- relevel(as.factor(data@colData$Apoe), ref = "E4")
      data@colData$Dx <- relevel(as.factor(data@colData$Dx), ref = "Ctrl")
      data@colData$Sex <- relevel(as.factor(data@colData$Sex), ref = "M")
      zlmCond <- zlm(~ Dx + Apoe + Sex + zAge + cngeneson + percent.mt + decontX_contamination + (1|Study/orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     force = TRUE,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('DxAD', 'SexF', 'ApoeE2', 'ApoeE3'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      
      write.csv(fcHurdle, paste0(Cluster, '/Comb.RefE4.', Chunk, '.csv'))
    }
  } else if (instance$strata == "Apoe"){
    data = instance$data
    Cluster = instance$Cluster
    APOEgenotype = instance$APOEgenotype
    Chunk = instance$Chunk
    data@colData$Dx <- relevel(as.factor(data@colData$Dx), ref = "Ctrl")
    data@colData$Sex <- relevel(as.factor(data@colData$Sex), ref = "M")
    zlmCond <- zlm(~ Dx + Sex + zAge + cngeneson + percent.mt + decontX_contamination + (1|Study/orig.ident),
                   method = 'glmer',
                   ebayes = FALSE,
                   fitArgsD = list(nAGQ = 0),
                   parallel = 4,
                   force = TRUE,
                   sca = data)
    
    summaryCond <- summary(zlmCond, doLRT= c('DxAD', 'SexF'))
    
    # print out the data.table
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                      summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
    
    write.csv(fcHurdle, paste0(Cluster, '/', APOEgenotype, '.', Chunk, '.csv'))
  } else {
    if (instance$Ref == "E2") {
      data = instance$data
      Cluster = instance$Cluster
      DX = instance$DX
      Chunk = instance$Chunk
      data@colData$Apoe <- relevel(as.factor(data@colData$Apoe), ref = "E2")
      
      zlmCond <- zlm(~ Apoe + zAge + cngeneson + percent.mt + decontX_contamination + (1|Study/orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     force = TRUE,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('ApoeE3', 'ApoeE4'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluster, "/", DX, '.', Chunk, '.RefE2.csv'))
    } else if (instance$Ref == "E3") {
      data = instance$data
      Cluster = instance$Cluster
      DX = instance$DX
      Chunk = instance$Chunk
      data@colData$Apoe <- relevel(as.factor(data@colData$Apoe), ref = "E3")
      
      zlmCond <- zlm(~ Apoe + zAge + cngeneson + percent.mt + decontX_contamination + (1|Study/orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     force = TRUE,
                     parallel = 4,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('ApoeE2', 'ApoeE4'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluster, "/", DX, '.', Chunk, '.RefE3.csv'))
    } else {
      data = instance$data
      Cluster = instance$Cluster
      DX = instance$DX
      Chunk = instance$Chunk
      data@colData$Apoe <- relevel(as.factor(data@colData$Apoe), ref = "E4")
      
      zlmCond <- zlm(~ Apoe + zAge + cngeneson + percent.mt + decontX_contamination + (1|Study/orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     force = TRUE,
                     parallel = 4,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('ApoeE2', 'ApoeE3'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluster, "/", DX, '.', Chunk, '.RefE4.csv'))
    }
  }
}

addAlgorithm(name = "MASTDEG", fun = MASTDEG)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design0 <- CJ(strata = "comb",
                   Ref = c("E2", "E3", "E4"),
                   i = 1:6,
                   j = 1:20)

prob.design1 <- CJ(strata = "Apoe",
                   APOEgenotype = c("E2", "E3", "E4"),
                   i = 1:6,
                   j = 1:20)

prob.design2 <- CJ(strata = "Dx",
                   DX = c("Ctrl", "AD"),
                   Ref = c("E2", "E3", "E4"),
                   i = 1:6,
                   j = 1:20)


prob.design <- rbind(prob.design0, prob.design1, prob.design2, fill = TRUE)

# convert NA to NULL
prob.design[is.na(prob.design)] <- ""


# final problem design
pdes <- list(GetSCAs = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(MASTDEG =data.table())

addExperiments(pdes, ades) # pdes set parameter to the data part; ades set parameters to the instance parts

# summarize the experiment
summarizeExperiments()

# Submitting and Collecting Results
options(batchtools.progress = FALSE)
submitJobs() 
# getStatus()
# getErrorMessages()
# getJobTable()
# findExpired()














