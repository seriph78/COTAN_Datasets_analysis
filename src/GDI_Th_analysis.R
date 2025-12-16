
sourceDir <- file.path(".", "src")

if (!dir.exists(sourceDir)) {
  sourceDir <- "."
}

source(file.path(sourceDir, "ClusterSplitting.R"))

extractResInfo <- function(resData) {
  c(NumClCells, GDI.df, Clusters, Genes.for.clust,
    NumCellsForAnal, Genes.for.anal, DEA.df, PVal.df,
    LR.model, Validation.assessment, Significant.genes, Fit.plot) %<-% resData

  quantiles <- quantile(GDI.df$GDI, probs = 1.0 - c(0.01, 0.005))
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
    "NumGenesWithGDIAbove_1.33" = sum(GDI.df$GDI > 1.33),
    "NumGenesWithGDIAbove_1.37" = sum(GDI.df$GDI > 1.37),
    "NumGenesWithGDIAbove_1.40" = sum(GDI.df$GDI > 1.40),
    "NumGenesWithGDIAbove_1.43" = sum(GDI.df$GDI > 1.43),
    "GDIQuantile_99"   = quantiles[[1L]],
    "GDIQuantile_99.5" = quantiles[[2L]],
    "Top10thGDIGene" = GDI.df$GDI[posSortedGDIGenes[10]],
    "Top20thGDIGene" = GDI.df$GDI[posSortedGDIGenes[10]],
    "Top30thGDIGene" = GDI.df$GDI[posSortedGDIGenes[10]],
    "Top50thGDIGene" = GDI.df$GDI[posSortedGDIGenes[10]],
    "MisclassificationError" = Validation.assessment$class[[1L]],
    "AreaUnderROCCurve" = Validation.assessment$auc[[1L]],
    "NumTopGenesConsidered" = numTopGenes,
    "NumLRSignificantGenes" = numLRSignGenes,
    "FractionTopGDIInHighDEAGenes" = genesFrGDIInDEA,
    "FractionLRSignInTopGDIGenes"  = genesFrLRSInGDI,
    "FractionLRSignInHighDEAGenes" = genesFrLRSInDEA,
    "NumGenesWithPValBelow5Perc" =  sum(PVal.df[, 1] < 0.05)
  )

#  print(res)

  return(res)
}

  # print(paste0(attr(Validation.assessment$class,"measure"),": ",
  #              Validation.assessment$class[[1L]]))
  #
  # print(paste0(attr(Validation.assessment$auc,"measure"),": ",
  #              Validation.assessment$auc[[1L]]))
  #
  # print("Number of genes with GDI above")
  # print(paste0(1.33, ": ", sum(GDI.df$GDI > 1.33)))
  # print(paste0(1.37, ": ", sum(GDI.df$GDI > 1.37)))
  # print(paste0(1.40, ": ", sum(GDI.df$GDI > 1.40)))
  # print(paste0(1.43, ": ", sum(GDI.df$GDI > 1.43)))
  #
  # print("GDI quantiles:")
  # print()
  #
  #
  # print("n-th GDI value")
  # print(paste0(10, ": ", GDI.df$GDI[posSortedGDIGenes[10]]))
  # print(paste0(20, ": ", GDI.df$GDI[posSortedGDIGenes[20]]))
  # print(paste0(30, ": ", GDI.df$GDI[posSortedGDIGenes[30]]))
  # print(paste0(50, ": ", GDI.df$GDI[posSortedGDIGenes[50]]))
  #
  #
  # print(paste0("Fraction of the top ", numTopGenes, " GDI genes among the highest ",
  #              numTopGenes, " DEA genes: ",
  #              (sum(rownames(GDI.df)[posTopGDIGenes] %in% rownames(DEA.df)[posTopDEAGenes]) /
  #                 numTopGenes)))
  #
  # print(paste0("Fraction among the ", length(Significant.genes$names),
  #              " LR significant genes in top ", numTopGenes, " GDI genes: ",
  #              (sum(Significant.genes$names %in% rownames(GDI.df)[posTopGDIGenes]) /
  #                 length(Significant.genes$names))))
  #
  # print(paste0("Fraction among the ", length(Significant.genes$names),
  #              " LR significant genes in highest ", numTopGenes, " DEA genes: ",
  #       (sum(Significant.genes$names %in% rownames(DEA.df)[posTopDEAGenes]) /
  #         length(Significant.genes$names))))
  #
  # print(paste0("Number of low p-value genes [5%] out of ",
  #              nrow(PVal.df), ": ", sum(PVal.df[, 1] < 0.05)))
  #
  # tryCatch({
  #   plot(Fit.plot)
  # },
  # error = function(err) {
  #   print(err)
  # })

dataDir <- file.path(".", "Data", "DataForGDITh")

if (!dir.exists(dataDir)) {
  dataDir <- "."
}

outDir <- file.path(".", "Analisys")
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

setLoggingLevel(2L)
setLoggingFile(file.path(outDir, "TestAnalisys.log"))

options(parallelly.fork.enable = TRUE)

# ------

obj.CD14 <- readRDS(file.path(dataDir, "CD14_Monocytes.cotan.RDS"))

clusters.CD14 <- getClusters(obj.CD14, "FromSeurat1.520")

test0.CD14 <- splitClusterAndRelatedData(obj.CD14,
                                         clToUse = "0",
                                         allClusters = clusters.CD14,
                                         geneRatioForClustering = 0.5,
                                         train.test.partition = 0.8)


test1.CD14 <- splitClusterAndRelatedData(obj.CD14,
                                         clToUse = "1",
                                         allClusters = clusters.CD14,
                                         geneRatioForClustering = 0.5,
                                         train.test.partition = 0.8)


# ------

obj.PBMC3 <- readRDS(file.path(dataDir, "PBMC3.cotan.RDS"))

clusters.PBMC3 <- getClusters(obj.PBMC3, clName = "FromSeurat1.520")

test0.PBMC3 <- splitClusterAndRelatedData(obj.PBMC3,
                                          clToUse = "0",
                                          allClusters = clusters.PBMC3,
                                          geneRatioForClustering = 0.5,
                                          train.test.partition = 0.8)

test1.PBMC3 <- splitClusterAndRelatedData(obj.PBMC3,
                                          clToUse = "1",
                                          allClusters = clusters.PBMC3,
                                          geneRatioForClustering = 0.5,
                                          train.test.partition = 0.8)

test2.PBMC3 <- splitClusterAndRelatedData(obj.PBMC3,
                                          clToUse = "2",
                                          allClusters = clusters.PBMC3,
                                          geneRatioForClustering = 0.5,
                                          train.test.partition = 0.8)

# ------

obj.CC <- readRDS(file.path(dataDir, "CorticalCells_GSM2861511_E135.cotan.RDS"))

clusters.CC <- getClusters(obj.CC, clName = "FromSeurat1.520")

test0.CC <- splitClusterAndRelatedData(obj.CC,
                                       clToUse = "0",
                                       allClusters = clusters.CC,
                                       geneRatioForClustering = 0.5,
                                       train.test.partition = 0.8)

test1.CC <- splitClusterAndRelatedData(obj.CC,
                                       clToUse = "1",
                                       allClusters = clusters.CC,
                                       geneRatioForClustering = 0.5,
                                       train.test.partition = 0.8)

# ------

allResData <- list("test0.CD14"  = test0.CD14,
                   "test1.CD14"  = test1.CD14,
                   "test0.PBMC3" = test0.PBMC3,
                   "test1.PBMC3" = test1.PBMC3,
                   "test2.PBMC3" = test2.PBMC3,
                   "test0.CC"    = test0.CC,
                   "test1.CC"    = test1.CC)

resDF <- data.frame()

for (resDataName in names(allResData)) {
  message("* ", appendLF = FALSE)
  row <- extractResInfo(allResData[[resDataName]])
  resDF <- rbind(resDF, row)
  rownames(resDF)[nrow(resDF)] <- resDataName
}

View(resDF)
