##########################################################################################################################################################
## Note: the data used for this analysis were downloaded from
## https://www.synapse.org/#!Synapse:syn21682120
## load required libraries
setwd(".../Step9.IntegrationOfPublicData/Zhou")
library(tidyverse)
library(synapser)
library(biomaRt)
library(ggbeeswarm)

## download the Normalised, covariates and surrogate variable adjusted residual expression 
synapser::synLogin(email = "YourEmail@gmail.com", password = "YourPassword")

############################ get the matrix data from synapse #####################
SynapseID <- c("syn21190625", "syn21190626", "syn21190635", "syn21190852", 
               "syn21190854", "syn21190861", "syn21190757", "syn21190758", 
               "syn21190759", "syn21190766", "syn21190768", "syn21190779", 
               "syn21190094", "syn21190095", "syn21190142", "syn21190788", 
               "syn21190787", "syn21190791", "syn21190781", "syn21190782", 
               "syn21190785", "syn21190871", "syn21190872", "syn21190873", 
               "syn21190867", "syn21190868", "syn21190870", "syn21190752", 
               "syn21190754", "syn21190756", "syn21190642", "syn21190643",
               "syn21190647", "syn21190691", "syn21190692", "syn21190693", 
               "syn21190609", "syn21190611", "syn21190612", "syn21190874", 
               "syn21190875", "syn21190913", "syn21190835", "syn21190836", 
               "syn21190842", "syn21190845", "syn21190846", "syn21190851", 
               "syn21190694", "syn21190695", "syn21190696", "syn21190862", 
               "syn21190863", "syn21190866", "syn21190702", "syn21190703", 
               "syn21190921", "syn21190793", "syn21190796", "syn21190834", 
               "syn21190613", "syn21190615", "syn21190624", "syn21190592", 
               "syn21190593", "syn21190604", "syn21190503", "syn21190504", 
               "syn21190510", "syn21190511", "syn21190512", "syn21190527", 
               "syn21190417", "syn21190425", "syn21190490", "syn21190497", 
               "syn21190498", "syn21190500", "syn21190144", "syn21190145", 
               "syn21190146", "syn21190539", "syn21190540", "syn21190544", 
               "syn21190494", "syn21190495", "syn21190496", "syn21190333", 
               "syn21190338", "syn21190394", "syn21190257", "syn21190262",
               "syn21190314", "syn21190528", "syn21190530", "syn21190538")

for (i in 1:length(SynapseID)){
  Dat <- synapser::synGet(SynapseID[i])
  file.copy(from = Dat$path,
            to = paste0(".../Step9.IntegrationOfPublicData/Zhou/RawDat/", 
                        grep("*.gz", unlist(str_split(Dat$path, "\\/", n=Inf)), value = TRUE)))
}


## Get the metadata
## metadata from the paper
setwd(".../Step9.IntegrationOfPublicData/Zhou")
Metadata.paper <- read.csv("MetafromPaper.csv")

# metadata from synapse
Meta.Biospecimen <- synapser::synGet("syn21738956")
Meta.Biospecimen <- read.table(Meta.Biospecimen$path, stringsAsFactors = F, header=T,
                               sep = ",")


## individual
Meta.individual <- synapser::synGet("syn3191087")
Meta.individual <- read.csv(Meta.individual$path)

## Merge the two metadata
Metadata <- Meta.Biospecimen %>% 
  left_join(Meta.individual %>% 
              dplyr::select(projid, Study, individualID), 
            by = "individualID") %>%
  left_join(Metadata.paper, by = c("projid" = "SampleID")) %>%
  dplyr::select(specimenID, individualID,  tissue, projid, Study, Sex, Age, Apoe, Trem2,
         `cogdx...CDR`, braaksc, ceradsc, Education, pmi, amyloid, tangles, SampleSource)

## rename "cogdx...CDR"
names(Metadata)[names(Metadata) == "cogdx...CDR"] <- "Dx"
Metadata$Dx[Metadata$SampleSource == "AD - Rush"] <- "AD-CV"
Metadata$Dx[Metadata$SampleSource == "Control - Rush"] <- "Ctrl"
Metadata$Dx[Metadata$SampleSource == "TREM2 - R62H - Rush"] <- "AD-R62H"

# save the data
write.csv(Metadata, "Metadata.csv", row.names = FALSE)


## Organize the data
setwd(".../Step9.IntegrationOfPublicData/Zhou")
files <- list.files("RawDat")

# get all samples
Samples <- str_split_fixed(files, "\\_", n = Inf)[, 1]

## load the metadata
Metadata <- read.csv("Metadata.csv")

## Check whether the data is complete
setdiff(unique(Samples), Metadata$specimenID)
setdiff(Metadata$specimenID, unique(Samples))

# ## In bash
# for i in "AD1"  "AD10" "AD11" "AD12" "AD13" "AD2"  "AD3"  "AD5"  "AD7" "AD8"  "AD9"  "C1" "C11"  "C12"  "C2"   "C3"   "C4"   "C5"   "C6"   "C7" "C8"   "C9"   "P10"  "P11"  "P12"  "P13"  "P2"   "P3"   "P5"   "P6"   "P7"   "P9"  
# do 
# mkdir $i
# find . \( -name $i\_\*.gz \) -exec mv -t $i {} +
#   done
# 
# ## Edit the file names
# for dir in .../Step9.IntegrationOfPublicData/Zhou/RawDat/*
#   do
# cd $dir
# for i in *.gz
# do 
# newName=$(echo "$i" | tr -cd "barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz")
# mv $i $newName
# done
# done































