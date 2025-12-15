
library(zeallot)

sourceDir <- file.path(".", "src")

if (!dir.exists(sourceDir)) {
  sourceDir <- "."
}

source(file.path(sourceDir, "ClusterSplitting.R"))

printGenes <- function(x) cat(paste0(x, collapse = "\n"))

extractResInfo <- function(resData) {
  c(NumClCells, GDI.df, Clusters, Genes.for.clust,
    NumCellsForAnal, Genes.for.anal, DEA.df, PVal.df,
    LR.model, Validation.assessment, Significant.genes) %<-% resData

  quantiles <- quantile(GDI.df$GDI,
                        probs = 1.0 - c(0.05, 0.03, 0.02, 0.01, 0.005, 0.0025, 0.001))
  posSortedGDIGenes <- order(GDI.df$GDI, decreasing = TRUE)

  numTopGenes <- 200L
  posTopGDIGenes <- posSortedGDIGenes[seq_len(numTopGenes)]
  posTopDEAGenes <- order(DEA.df[, 1], decreasing = T)[seq_len(numTopGenes)]

  stopifnot(attr(Validation.assessment$class,"measure") == "Misclassification Error")
  stopifnot(attr(Validation.assessment$aus,"measure")   == "AUS")

  numLRSignGenes <- length(Significant.genes$names)

  genesFrGDIInDEA <- sum(rownames(GDI.df)[posTopGDIGenes] %in%
                           rownames(DEA.df)[posTopDEAGenes]) / numTopGenes

  genesFrLRSInGDI <- sum(Significant.genes$names %in%
                           rownames(GDI.df)[posTopGDIGenes]) / numLRSignGenes

  genesFrLRSInDEA <- sum(Significant.genes$names %in%
                           rownames(DEA.df)[posTopDEAGenes]) / numLRSignGenes

  res <- list(
    "NumClCells" = NumClCells,
    "NumClGenes" = nrow(GDI.df),
    "NumCellsCl_0" = sum(Clusters == "0"),
    "NumCellsCl_1" = sum(Clusters == "1"),
    "NumGenesForClust" = length(Genes.for.clust),
    "NumCellsForAnal" = NumCellsForAnal,
    "NumGenesForAnal" = length(Genes.for.anal),
    "NumGenesWithGDIAbove_1.30" = sum(GDI.df$GDI > 1.30),
    "NumGenesWithGDIAbove_1.33" = sum(GDI.df$GDI > 1.33),
    "NumGenesWithGDIAbove_1.37" = sum(GDI.df$GDI > 1.37),
    "NumGenesWithGDIAbove_1.40" = sum(GDI.df$GDI > 1.40),
    "NumGenesWithGDIAbove_1.43" = sum(GDI.df$GDI > 1.43),
    "NumGenesWithGDIAbove_1.50" = sum(GDI.df$GDI > 1.50),
    "PercGenesWithGDIAbove_1.30" = sum(GDI.df$GDI > 1.30) / nrow(GDI.df) * 100.0,
    "PercGenesWithGDIAbove_1.33" = sum(GDI.df$GDI > 1.33) / nrow(GDI.df) * 100.0,
    "PercGenesWithGDIAbove_1.37" = sum(GDI.df$GDI > 1.37) / nrow(GDI.df) * 100.0,
    "PercGenesWithGDIAbove_1.40" = sum(GDI.df$GDI > 1.40) / nrow(GDI.df) * 100.0,
    "PercGenesWithGDIAbove_1.43" = sum(GDI.df$GDI > 1.43) / nrow(GDI.df) * 100.0,
    "PercGenesWithGDIAbove_1.46" = sum(GDI.df$GDI > 1.46) / nrow(GDI.df) * 100.0,
    "PercGenesWithGDIAbove_1.50" = sum(GDI.df$GDI > 1.50) / nrow(GDI.df) * 100.0,
    "GDIQuantile_95"    = quantiles[[1L]],
    "GDIQuantile_97"    = quantiles[[2L]],
    "GDIQuantile_98"    = quantiles[[3L]],
    "GDIQuantile_99"    = quantiles[[4L]],
    "GDIQuantile_99.5"  = quantiles[[5L]],
    "GDIQuantile_99.75" = quantiles[[6L]],
    "GDIQuantile_99.9"  = quantiles[[7L]],
    "Top1stGDIGene"  = GDI.df$GDI[posSortedGDIGenes[1]],
    "Top2ndGDIGene"  = GDI.df$GDI[posSortedGDIGenes[2]],
    "Top5thGDIGene"  = GDI.df$GDI[posSortedGDIGenes[5]],
    "Top10thGDIGene" = GDI.df$GDI[posSortedGDIGenes[10]],
    "Top20thGDIGene" = GDI.df$GDI[posSortedGDIGenes[20]],
    "Top30thGDIGene" = GDI.df$GDI[posSortedGDIGenes[30]],
    "Top50thGDIGene" = GDI.df$GDI[posSortedGDIGenes[50]],
    "MisclassificationError" = Validation.assessment$class[[1L]],
    "AreaUnderROCCurve" = Validation.assessment$auc[[1L]],
    "NumTopGenesConsidered" = numTopGenes,
    "NumLRSignificantGenes" = numLRSignGenes,
    "FractionTopGDIInHighDEAGenes" = genesFrGDIInDEA,
    "FractionLRSignInTopGDIGenes"  = genesFrLRSInGDI,
    "FractionLRSignInHighDEAGenes" = genesFrLRSInDEA,
    "NumGenesWithPValBelow5Perc" =  sum(PVal.df[, 1] < 0.05)
  )

  return(res)
}

loopOnAllClusters <- function(obj, datasetName, clusters, outDir) {
  if (!is.factor(clusters)) {
    clusters <- factor(clusters)
  }

  summaryDF <- data.frame()
  for (clName in levels(clusters)) {
    testDataName <- paste0(datasetName, "_CL-", clName)

    clSize <- sum(clusters == clName)
    if (clSize < 100L) {
      logThis(paste0("Skipped analysis of ", testDataName,
                     ": it has only ", clSize, " cells"), appendLF = TRUE)
      next
    }

    message("*", appendLF = FALSE)

    analysisData <-
      splitClusterAndRelatedData(obj,
                                 clToUse = clName,
                                 allClusters = clusters,
                                 geneRatioForClustering = 0.5,
                                 train.test.partition = 0.8)

    logThis(paste0("Saving analysis of ", testDataName))

    saveRDS(analysisData, file.path(outDir,
              paste0(testDataName, "_AnalysisData.RDS")))

    message("+", appendLF = FALSE)

    plot(analysisData$LR.model)

    summaryData <- extractResInfo(analysisData)
    summaryDF <- rbind(summaryDF, summaryData)
    rownames(summaryDF)[nrow(summaryDF)] <- testDataName

    message(" ", appendLF = FALSE)
  }

  return(summaryDF)
}

runAllAnalysis <- function(datasetNames) {
  allSummariesDF <- data.frame()
  for (datasetName in datasetNames) {
    logThis(paste0("Analysis of ", datasetName))

    obj <- readRDS(file.path(dataDir, paste0(datasetName, ".cotan.RDS")))
    for(clstName in getClusterizations(obj)) {
      logThis(paste0("Checking clusterization ", clstName))
      clusters <- getClusters(obj, clstName)

      summaryDF <- loopOnAllClusters(obj, datasetName = datasetName,
                                     clusters = clusters, outDir = outDir)

      allSummariesDF <- rbind(allSummariesDF, summaryDF)
    }
  }
  return(allSummariesDF)
}


loopOnAllAnalysisData <- function(outDir, whichDFs) {
  summaryDF <- data.frame()

  for (file in list.files(outDir, pattern = "_AnalysisData.RDS")) {
    if (stringr::str_detect(file, "e15.0") && whichDFs == "Train") {
      next
    }
    if (!stringr::str_detect(file, "e15.0") && whichDFs == "Test") {
      next
    }
    analysisData <- readRDS(file = file.path(outDir, file))

    summaryData <- extractResInfo(analysisData)
    summaryDF <- rbind(summaryDF, summaryData)
    rownames(summaryDF)[nrow(summaryDF)] <- stringr::str_remove(file, "_AnalysisData.RDS")
  }

  return(summaryDF)
}

dataDir <- file.path(".", "Data", "DataForGDITh")

if (!dir.exists(dataDir)) {
  dataDir <- "."
}

outDir <- file.path(".", "Analysis")
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

winCOTANDir <- file.path("~", "work", "COTAN",
                         "ClustersStatiticalAnalysis", "Analysis")

setLoggingLevel(2L)
setLoggingFile(file.path(outDir, "TestAnalysis.log"))

options(parallelly.fork.enable = TRUE)

# ------

trainDatasetNames <- c("CD14_Monocytes",
                       "CorticalCells_GSM2861511_E135",
                       "e15.5_ForebrainDorsal",
                       "PBMC3")

trainSummariesDF <- runAllAnalysis(datasetNames)

write.csv(trainSummariesDF, file = file.path(outDir, "TrainAnalysisSummaries.csv"),
          row.names = TRUE)

trainSummariesDF <- read.csv(file = file.path(outDir, "TrainAnalysisSummaries.csv"),
                             row.names = 1L, header = TRUE)

View(trainSummariesDF)

# -------------

testDatasetNames <- c("e15.0_ForebrainDorsal", "e15.0_ForebrainDorsal.2.3ris")

testSummariesDF <- runAllAnalysis(testDatasetNames)

write.csv(testSummariesDF, file = file.path(outDir, "TestAnalysisSummaries.csv"),
          row.names = TRUE)

testSummariesDF <- read.csv(file = file.path(outDir, "TestAnalysisSummaries.csv"),
                            row.names = 1L, header = TRUE)

View(testSummariesDF)

goodSplit2 <- (testSummariesDF$NumCellsCl_1 / testSummariesDF$NumCellsCl_0 > 0.33)

View(testSummariesDF[!goodSplit2, ])


# -------------------

goodSplit <- (trainSummariesDF$NumCellsCl_1 / trainSummariesDF$NumCellsCl_0 > 0.33)

filteredSummariesDF <- trainSummariesDF[goodSplit, ]

View(filteredSummariesDF)

colFilteredSummDF <- filteredSummariesDF[, colnames(filteredSummariesDF)[c(1,3:4,6,18:19,21,24:25)]]

colFilteredSummDF$ClusterWasUniform <- 0.5

View(colFilteredSummDF)

surelyNotUniform <- colFilteredSummDF$MisclassificationError < 0.2 & colFilteredSummDF$AreaUnderROCCurve > 0.81

colFilteredSummDF$ClusterWasUniform[surelyNotUniform] <- 0.0

surelyUniform <- colFilteredSummDF$MisclassificationError > 0.377 | colFilteredSummDF$AreaUnderROCCurve < 0.55

colFilteredSummDF$ClusterWasUniform[surelyUniform] <- 1.0

sum(colFilteredSummDF$ClusterWasUniform == 0.5)

colFilteredSummDF["PBMC3_CL-9",                         "ClusterWasUniform"] <- 0.75
colFilteredSummDF["CorticalCells_GSM2861511_E135_CL-0", "ClusterWasUniform"] <- 0.75
colFilteredSummDF["CorticalCells_GSM2861511_E135_CL-2", "ClusterWasUniform"] <- 0.75
colFilteredSummDF["PBMC3_CL-17",                        "ClusterWasUniform"] <- 0.75
colFilteredSummDF["e15.5_ForebrainDorsal_CL-4",         "ClusterWasUniform"] <- 0.75
colFilteredSummDF["PBMC3_CL-10",                        "ClusterWasUniform"] <- 0.75
colFilteredSummDF["PBMC3_CL-21",                        "ClusterWasUniform"] <- 0.50
colFilteredSummDF["e15.5_ForebrainDorsal_CL-7",         "ClusterWasUniform"] <- 0.25
colFilteredSummDF["PBMC3_CL-20",                        "ClusterWasUniform"] <- 0.25
colFilteredSummDF["e15.5_ForebrainDorsal_CL-10",        "ClusterWasUniform"] <- 0.75
colFilteredSummDF["PBMC3_CL-0",                         "ClusterWasUniform"] <- 0.25
colFilteredSummDF["PBMC3_CL-19",                        "ClusterWasUniform"] <- 0.25
colFilteredSummDF["PBMC3_CL-15",                        "ClusterWasUniform"] <- 0.25
colFilteredSummDF["PBMC3_CL-13",                        "ClusterWasUniform"] <- 0.75
colFilteredSummDF["e15.5_ForebrainDorsal_CL-3",         "ClusterWasUniform"] <- 0.75
colFilteredSummDF["CD14_Monocytes_CL-0",                "ClusterWasUniform"] <- 0.25
colFilteredSummDF["PBMC3_CL-12",                        "ClusterWasUniform"] <- 0.25
colFilteredSummDF["PBMC3_CL-8",                         "ClusterWasUniform"] <- 0.75

write.csv(colFilteredSummDF,
          file = file.path(outDir, "OurAnalysisOnSubsetSummaries.csv"),
          row.names = TRUE, col.names = TRUE)

colFilteredSummDF <- read.csv(
          file = file.path(outDir, "OurAnalysisOnSubsetSummaries.csv"),
          row.names = 1L, header = TRUE)

View(colFilteredSummDF)

#-------------

trainSummariesDF <- loopOnAllAnalysisData(outDir, whichDFs = "Train")

goodSplit <- (trainSummariesDF$NumCellsCl_1 / trainSummariesDF$NumCellsCl_0 > 0.33)

trainSummariesDF$ClusterWasUniform <- NA

trainSummariesDF[rownames(colFilteredSummDF), "ClusterWasUniform"] <-
  colFilteredSummDF[, "ClusterWasUniform"]

sum(is.na(trainSummariesDF$ClusterWasUniform) & !goodSplit)

trainSummariesDF[!goodSplit, "ClusterWasUniform"] <- 0.0

sum(is.na(trainSummariesDF$ClusterWasUniform))

write.csv(trainSummariesDF, file = file.path(outDir, "TrainAnalysisSummaries.csv"),
          row.names = TRUE)

trainSummariesDF <- read.csv(file = file.path(outDir, "TrainAnalysisSummaries.csv"),
                           row.names = 1L, header = TRUE)

View(trainSummariesDF)

#trainSummariesDF_Backup <- trainSummariesDF
trainSummariesDF <- trainSummariesDF_Backup
trainSummariesDF <- trainSummariesDF[-1:-9, ]

#-------------

plot(x = allSummariesDF$NumGenesWithGDIAbove_1.30, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$NumGenesWithGDIAbove_1.33, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$NumGenesWithGDIAbove_1.37, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$NumGenesWithGDIAbove_1.40, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$NumGenesWithGDIAbove_1.43, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$NumGenesWithGDIAbove_1.50, y = allSummariesDF$ClusterWasUniform)

plot(x = log10(allSummariesDF$PercGenesWithGDIAbove_1.30), y = allSummariesDF$ClusterWasUniform)
plot(x = log10(allSummariesDF$PercGenesWithGDIAbove_1.33), y = allSummariesDF$ClusterWasUniform)
plot(x = log10(allSummariesDF$PercGenesWithGDIAbove_1.37), y = allSummariesDF$ClusterWasUniform)
plot(x = log10(allSummariesDF$PercGenesWithGDIAbove_1.40), y = allSummariesDF$ClusterWasUniform)
plot(x = log10(allSummariesDF$PercGenesWithGDIAbove_1.43), y = allSummariesDF$ClusterWasUniform)
plot(x = log10(allSummariesDF$PercGenesWithGDIAbove_1.50), y = allSummariesDF$ClusterWasUniform)

plot(x = allSummariesDF$GDIQuantile_99,    y = log10(allSummariesDF$NumClCells))
plot(x = allSummariesDF$GDIQuantile_99 * log10(allSummariesDF$NumClCells), y = allSummariesDF$ClusterWasUniform)

plot(x = allSummariesDF$GDIQuantile_99,    y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$GDIQuantile_99.5,  y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$GDIQuantile_99.75, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$GDIQuantile_99.9,  y = allSummariesDF$ClusterWasUniform)

plot(x = allSummariesDF$Top1stGDIGene,  y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$Top2ndGDIGene,  y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$Top5thGDIGene,  y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$Top10thGDIGene, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$Top20thGDIGene, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$Top30thGDIGene, y = allSummariesDF$ClusterWasUniform)
plot(x = allSummariesDF$Top50thGDIGene, y = allSummariesDF$ClusterWasUniform)


gdiFilteredSummaryDF <- allSummariesDF[allSummariesDF$Top1stGDIGene < 1.5, ]

plot(x = gdiFilteredSummaryDF$GDIQuantile_99.9,    y = gdiFilteredSummaryDF$ClusterWasUniform)

plot(x = allSummariesDF$Top1stGDIGene,  y = allSummariesDF$AreaUnderROCCurve)
plot(x = allSummariesDF$Top5thGDIGene,  y = allSummariesDF$AreaUnderROCCurve)
plot(x = allSummariesDF$GDIQuantile_99,  y = allSummariesDF$AreaUnderROCCurve)
plot(x = allSummariesDF$GDIQuantile_99.9,  y = allSummariesDF$AreaUnderROCCurve)

plot(x = allSummariesDF$GDIQuantile_99,
     y = allSummariesDF$AreaUnderROCCurve / (1 + allSummariesDF$ClusterWasUniform))
plot(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.37, 3e-2)),
     y = allSummariesDF$AreaUnderROCCurve / (1 + allSummariesDF$ClusterWasUniform))

plot(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.33, 1e-3)), y = allSummariesDF$AreaUnderROCCurve)
plot(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.37, 3e-2)), y = allSummariesDF$AreaUnderROCCurve)
plot(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.40, 1e-3)), y = allSummariesDF$AreaUnderROCCurve)
plot(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.43, 1e-3)), y = allSummariesDF$AreaUnderROCCurve)
plot(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.50, 1e-3)), y = allSummariesDF$AreaUnderROCCurve)

plot(x = log10(pmax(allSummariesDF$NumGenesWithGDIAbove_1.40, 1)), y = allSummariesDF$AreaUnderROCCurve)

stats::cor(x = allSummariesDF$Top1stGDIGene,  y = allSummariesDF$AreaUnderROCCurve)
stats::cor(x = allSummariesDF$GDIQuantile_99.5,  y = allSummariesDF$AreaUnderROCCurve)
stats::cor(x = allSummariesDF$PercGenesWithGDIAbove_1.37, y = allSummariesDF$AreaUnderROCCurve)
stats::cor(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.33, 1e-3)), y = allSummariesDF$AreaUnderROCCurve)
stats::cor(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.37, 3e-2)), y = allSummariesDF$AreaUnderROCCurve)
stats::cor(x = log10(pmax(allSummariesDF$PercGenesWithGDIAbove_1.40, 1e-3)), y = allSummariesDF$AreaUnderROCCurve)
stats::cor(x = log10(pmax(allSummariesDF$NumGenesWithGDIAbove_1.40, 1)), y = allSummariesDF$AreaUnderROCCurve)

# ------

testSummariesDF$ClusterWasUniform <- NA

surelyNotUniform2 <- goodSplit2 & testSummariesDF$MisclassificationError < 0.2 & testSummariesDF$AreaUnderROCCurve > 0.81

testSummariesDF$ClusterWasUniform[surelyNotUniform2] <- 0.0

testSummariesDF$ClusterWasUniform[!goodSplit2] <- 0.0

surelyUniform2 <- goodSplit2 & (testSummariesDF$MisclassificationError > 0.377 | testSummariesDF$AreaUnderROCCurve < 0.55)

testSummariesDF$ClusterWasUniform[surelyUniform2] <- 1.0

sum(is.na(testSummariesDF$ClusterWasUniform))

View(testSummariesDF)

write.csv(testSummariesDF, file = file.path(outDir, "TestAnalysisSummaries.csv"),
          row.names = TRUE)

allSummariesDF <- rbind(trainSummariesDF, testSummariesDF)

sum(is.na(allSummariesDF$ClusterWasUniform))

write.csv(allSummariesDF, file = file.path(outDir, "AllAnalysisSummaries.csv"),
          row.names = TRUE)

# ------

extractGDIDFFromAllAnalysisData <- function(outDir, whichDFs) {
  for (file in list.files(outDir, pattern = "_AnalysisData.RDS")) {
    if (stringr::str_detect(file, "e15.0") && whichDFs == "Train") {
      next
    }
    if (!stringr::str_detect(file, "e15.0") && whichDFs == "Test") {
      next
    }
    analysisData <- readRDS(file = file.path(outDir, file))

    gdiDF <- analysisData[["GDI.df"]]

    fileName <-
      file.path(outDir, paste0(stringr::str_remove(file, "_AnalysisData.RDS"), "_GDI.csv"))

    write.csv(gdiDF, file = fileName, row.names = TRUE)
  }
}

extractGDIDFFromAllAnalysisData(outDir = outDir, whichDFs = "All")
