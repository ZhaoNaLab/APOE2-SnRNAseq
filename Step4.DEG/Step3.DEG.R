# https://mllg.github.io/batchtools/articles/batchtools.html
library("batchtools")
setwd(".../Step4.DEG")

reg <- makeExperimentRegistry("MAST.DEG.Registry",
                              packages = c("MAST", "data.table", "Seurat", "lme4"),
                              seed = 22112019)

reg <- loadRegistry("MAST.DEG.Registry", writeable = TRUE)


# Step1: function to subset data
# sca file has been filtered and were siplit into chunks
scaObjs <- readRDS('scaObjs.rds')

# define instance (data preparation)
subsample <- function(data, job, strata, APOEgenotype, Ref, Dx, i) {
  
  Cluster <- names(data)[i]
  if (strata == "comb"){
    data = data[[i]]
  } else if (strata == "apoe"){
    data = subset(data[[i]], apoe == APOEgenotype)
  } else {
    data = subset(data[[i]], bg == Dx)
    } 
  
  list(data = data, 
       strata = strata, 
       APOEgenotype = APOEgenotype, 
       Ref = Ref, 
       Dx = Dx, 
       Cluster = Cluster)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "GetSCAs", data = scaObjs, fun = subsample) 

# create an algorithm which applies a support vector machine
MASTDEG <- function(data, job, instance) {# instance is inherited from the last step
  setwd("/path/to/save")
  if (instance$strata == "comb"){
    if(instance$Ref == "E2"){
      data = instance$data
      Cluser = instance$Cluster
      data@colData$apoe <- relevel(as.factor(data@colData$apoe), ref = "E2")
      data@colData$bg <- relevel(as.factor(data@colData$bg), ref = "Ctrl")
      data@colData$sex <- relevel(as.factor(data@colData$sex), ref = "male")
      zlmCond <- zlm(~ bg + apoe + sex + zAge + cngeneson + VaD + percent.mt + Seq_sat + Bulk_RIN + decontX_contamination + (1|orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('bgAD', 'sexfemale', 'apoeE3', 'apoeE4'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluser, '/Comb.RefE2.csv'))
    } else if (instance$Ref == "E3") {
      data = instance$data
      Cluser = instance$Cluster
      data@colData$apoe <- relevel(as.factor(data@colData$apoe), ref = "E3")
      data@colData$bg <- relevel(as.factor(data@colData$bg), ref = "Ctrl")
      data@colData$sex <- relevel(as.factor(data@colData$sex), ref = "male")
      zlmCond <- zlm(~ bg + apoe + sex + zAge + cngeneson + VaD + percent.mt + Seq_sat + Bulk_RIN + decontX_contamination + (1|orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('bgAD', 'sexfemale', 'apoeE2', 'apoeE4'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluser, '/Comb.RefE3.csv'))
    } else {
      data = instance$data
      Cluser = instance$Cluster
      data@colData$apoe <- relevel(as.factor(data@colData$apoe), ref = "E4")
      data@colData$bg <- relevel(as.factor(data@colData$bg), ref = "Ctrl")
      data@colData$sex <- relevel(as.factor(data@colData$sex), ref = "male")
      zlmCond <- zlm(~ bg + apoe + sex + zAge + cngeneson + VaD + percent.mt + Seq_sat + Bulk_RIN + decontX_contamination + (1|orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('bgAD', 'sexfemale', 'apoeE2', 'apoeE3'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      
      write.csv(fcHurdle, paste0(Cluser, '/Comb.RefE4.csv'))
    }
  } else if (instance$strata == "apoe"){
    data = instance$data
    Cluser = instance$Cluster
    APOEgenotype = instance$APOEgenotype
    data@colData$bg <- relevel(as.factor(data@colData$bg), ref = "Ctrl")
    data@colData$sex <- relevel(as.factor(data@colData$sex), ref = "male")
    zlmCond <- zlm(~ bg + sex + zAge + cngeneson + VaD + percent.mt + Seq_sat + Bulk_RIN + decontX_contamination + (1|orig.ident),
                   method = 'glmer',
                   ebayes = FALSE,
                   fitArgsD = list(nAGQ = 0),
                   parallel = 4,
                   sca = data)
    
    summaryCond <- summary(zlmCond, doLRT= c('bgAD', 'sexfemale'))
    
    # print out the data.table
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                      summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
    
    write.csv(fcHurdle, paste0(Cluser, '/', APOEgenotype, '.csv'))
  } else {
    if (instance$Ref == "E2") {
      data = instance$data
      Cluser = instance$Cluster
      Dx = instance$Dx
      data@colData$apoe <- relevel(as.factor(data@colData$apoe), ref = "E2")
      
      zlmCond <- zlm(~ apoe + zAge + cngeneson + VaD + percent.mt + Seq_sat + Bulk_RIN + decontX_contamination + (1|orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('apoeE3', 'apoeE4'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      bg <- unique(as.vector(data@colData$bg))
      
      write.csv(fcHurdle, paste0(Cluser, "/", Dx, '.RefE2.csv'))
    } else if (instance$Ref == "E3") {
      data = instance$data
      Cluser = instance$Cluster
      Dx = instance$Dx
      data@colData$apoe <- relevel(as.factor(data@colData$apoe), ref = "E3")
      
      zlmCond <- zlm(~ apoe + zAge + cngeneson + VaD + percent.mt + Seq_sat + Bulk_RIN + decontX_contamination + (1|orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('apoeE2', 'apoeE4'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluser, "/", Dx, '.RefE3.csv'))
    } else {
      data = instance$data
      Cluser = instance$Cluster
      Dx = instance$Dx
      data@colData$apoe <- relevel(as.factor(data@colData$apoe), ref = "E4")
      
      zlmCond <- zlm(~ apoe + zAge + cngeneson + VaD + percent.mt + Seq_sat + Bulk_RIN + decontX_contamination + (1|orig.ident),
                     method = 'glmer',
                     ebayes = FALSE,
                     fitArgsD = list(nAGQ = 0),
                     parallel = 4,
                     sca = data)
      
      summaryCond <- summary(zlmCond, doLRT= c('apoeE2', 'apoeE3'))
      
      # print out the data.table
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                        summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
      
      write.csv(fcHurdle, paste0(Cluser, "/", Dx, '.RefE4.csv'))
    }
  }
}

addAlgorithm(name = "MASTDEG", fun = MASTDEG)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design0 <- CJ(strata = "comb",
                   Ref = c("E2", "E3", "E4"),
                   i = 1:7)

prob.design1 <- CJ(strata = "apoe",
                   APOEgenotype = c("E2", "E3", "E4"),
                   i = 1:7)

prob.design2 <- CJ(strata = "bg",
                   Dx = c("Ctrl", "AD"),
                   Ref = c("E2", "E3", "E4"),
                   i = 1:7)


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











