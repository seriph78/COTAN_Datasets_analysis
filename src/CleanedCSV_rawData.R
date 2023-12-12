# Code that produce the csv data from the cleaned data to be used as 
# input for celltypist

library(ggplot2)
library(tibble)
library(zeallot)
library(Seurat)
library(COTAN)
options(parallelly.fork.enable = TRUE)

# CD14
outDir <- "Data/CD14Cleaned/"
sampleCondition <- "CD14_Monocytes"
cd14Obj <- readRDS(file = file.path(outDir, paste0(sampleCondition, ".cotan.RDS")))
dataset <- Read10X(file.path(outDir, "/OrigialDatahg19"))
write.csv(dataset[,getCells(cd14Obj)],file = paste0(outDir,sampleCondition, "_cleaned.csv"))

# CD14 raw not cleaned
outDir <- "Data/CD14Cleaned/"
dataset <- Read10X(file.path(outDir, "/OrigialDatahg19"))
write.csv(dataset,file = paste0(outDir, "CD14_raw_NOT_cleaned.csv"))

        
# Mouse Brain Le Manno - Loom file E15.0
outDir <- "Data/MouseCortexFromLoom/"
fb150Obj <- readRDS("Data/MouseCortexFromLoom/SourceData/e15.0_ForebrainDorsal.cotan.RDS")
fb150Obj.cl <- readRDS("Data/MouseCortexFromLoom/e15.0_ForebrainDorsal.cotan.RDS")
sampleCondition <- getMetadataElement(fb150Obj, datasetTags()[["cond"]])

write.csv(getRawData(fb150Obj)[,getCells(fb150Obj.cl)],file = paste0(outDir,sampleCondition, "_cleaned.csv"))

# Mouse Brain Le Manno - Loom file E17.5
outDir <- "Data/MouseCortexFromLoom/"
fb175Obj <- readRDS("Data/MouseCortexFromLoom/SourceData/e17.5_ForebrainDorsal.cotan.RDS")
fb175Obj.cl <- readRDS("Data/MouseCortexFromLoom/e17.5_ForebrainDorsal.cotan.RDS")

sampleCondition <- getMetadataElement(fb175Obj, datasetTags()[["cond"]])

write.csv(getRawData(fb175Obj)[,getCells(fb175Obj.cl)],file = paste0(outDir,sampleCondition, "_cleaned.csv"))

# Mouse Brain Le Manno - Loom file E13.5
outDir <- "Data/MouseCortexFromLoom/"
fb135Obj <- readRDS("Data/MouseCortexFromLoom/SourceData/e13.5_ForebrainDorsal.cotan.RDS")
fb135Obj.cl <- readRDS("Data/MouseCortexFromLoom/e13.5_ForebrainDorsal.cotan.RDS")

sampleCondition <- getMetadataElement(fb135Obj, datasetTags()[["cond"]])

write.csv(getRawData(fb135Obj)[,getCells(fb135Obj.cl)],file = paste0(outDir,sampleCondition, "_cleaned.csv"))

# Cortical cells DGE E17.5 
outDir <- "Data/Yuzwa_MouseCortex/"
cc175Obj <- readRDS("Data/Yuzwa_MouseCortex/CorticalCells_GSM2861514_E175.cotan.RDS")
sampleCondition <- getMetadataElement(cc175Obj, datasetTags()[["cond"]])
dataset <- read.csv(file.path(outDir <- "Data/Yuzwa_MouseCortex/"
                              , "GSM2861514_E175_Only_Cortical_Cells_DGE.txt"),
                    header = TRUE, sep = "\t", strip.white = TRUE,row.names = 1)

write.csv(dataset[,getCells(cc175Obj)],file = paste0(outDir,sampleCondition, "_cleaned.csv"))


# Cortical cells DGE E13.5 
outDir <- "Data/Yuzwa_MouseCortex/"
cc135Obj <- readRDS("Data/Yuzwa_MouseCortex/CorticalCells_GSM2861511_E135.cotan.RDS")
sampleCondition <- getMetadataElement(cc135Obj, datasetTags()[["cond"]])
dataset <- read.csv(file.path(outDir <- "Data/Yuzwa_MouseCortex/"
                              , "GSM2861511_E135_Only_Cortical_Cells_DGE.txt"),
                    header = TRUE, sep = "\t", strip.white = TRUE,row.names = 1)



write.csv(dataset[,getCells(cc135Obj)],file = paste0(outDir,sampleCondition, "_cleaned.csv"))


