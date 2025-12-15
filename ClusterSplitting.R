
library(Seurat)
library(COTAN)

splitCluster <- function(rawData) {
  srat <- CreateSeuratObject(counts = rawData,
                             project = paste0("cluster_splitting"),
                             min.cells = 3L, min.features = 4L)

  srat <- NormalizeData(srat)
  srat <- FindVariableFeatures(srat, selection.method = "vst",
                               nfeatures = 2000L)
  srat <- ScaleData(srat, features = rownames(srat))

  maxRows <- nrow(srat@meta.data) - 1L
  srat <- RunPCA(srat, features = VariableFeatures(object = srat),
                 npcs = min(50L, maxRows))

  srat <- FindNeighbors(srat, dims = 1L:min(25L, maxRows))

  lowResolution <- 0.1
  resolution <- 0.8
  highResolution <- 2.5
  clusters <- NULL

  iter <- 1L
  maxIterations <- 50L
  repeat {
    srat <- FindClusters(srat, resolution = resolution, algorithm = 2L)
    clusters <- factor(srat[["seurat_clusters", drop = TRUE]])
    numClusters <- nlevels(clusters)

    message("Found ", numClusters, " clusters at resolution ", resolution)

    if (numClusters == 2L) {
      break
    }

    if (iter >= maxIterations) {
      warning("Solution not found: max iterations reached")
      break
    }

    resolution <- 0.75 * resolution + 0.25 *
      ifelse(numClusters > 2L, lowResolution, highResolution)

    iter <- iter + 1L
  }

  return(clusters)
}

library(glmnet)
library(caret)

LR_function <- function(obj, clusters,
                        train.test.partition = 0.8) {

  data <- getLogNormData(obj)

  row_stdev <- apply(data, 1, sd, na.rm=TRUE)
  row_stdev <- row_stdev[order(row_stdev,decreasing = T)]

  genes.to.keep <- names(row_stdev)[seq_len(min(1000, nrow(data)))]

  data.small <- data[rownames(data) %in% genes.to.keep, ]

  data.small <- t(as.matrix(data.small))
  data.small <- scale(data.small,center = TRUE,scale = TRUE)

  Cl.code <- set_names(as.data.frame(as.integer(factorToVector(clusters))), "Cl.code")

  Cl.code$cell.Names <- getCells(obj)

  # Split the data into training and test set
  set.seed(123)
  training.samples <- Cl.code[,"Cl.code"] %>%
    createDataPartition(p = train.test.partition, list = FALSE)

  train.cl.code <- Cl.code[ training.samples, ]
  test.cl.code  <- Cl.code[-training.samples, ]

  train.data <- data.small[ training.samples, ]
  test.data  <- data.small[-training.samples, ]

  stopifnot(setequal(unique(Cl.code$Cl.code), 0L:1L))

  fit <- cv.glmnet(train.data,
                   train.cl.code$Cl.code,
                   family= "binomial",
                   type.measure = "class",
                   alpha = (1 - 1e-3))

  assessment <- assess.glmnet(fit,
                              newx = test.data,
                              newy = test.cl.code$Cl.code,
                              family= "binomial")

  coef <- as.matrix(coef(fit, s = "lambda.1se"))
  coef_names <- as.data.frame(coef)
  coef_names$names <- rownames(coef_names)
  rownames(coef_names) <- rownames(coef)

  sign.genes <- coef_names[coef_names$s1 != 0,]

  return(list("LR.model" = fit,
              "Validation.assessment" = assessment,
              "Significant.genes" = sign.genes))
}


splitClusterAndRelatedData <- function(objCOTAN,
                                       allClusters, clToUse,
                                       geneRatioForClustering,
                                       train.test.partition) {
  cellsToUse <- names(allClusters)[allClusters == clToUse]

  clObj <- dropGenesCells(objCOTAN, cells = getCells(objCOTAN)[!getCells(objCOTAN) %in% cellsToUse])
  clObj <- proceedToCoex(clObj, calcCoex = TRUE, cores = 6L, saveObj = FALSE, outDir = outDir)

  gdiDF <- calculateGDI(clObj, rowsFraction = 0.05)

  numGenes <- getNumGenes(clObj)
  numCells <- getNumCells(clObj)

  set.seed(137)
  genesPosForClustering <- seq_len(numGenes)
  genesPosForAnalysis <- genesPosForClustering

  if (geneRatioForClustering > 0.0 && geneRatioForClustering < 1.0) {
    genesPosForClustering <- sample(genesPosForClustering,
                                    size = ceiling(numGenes * geneRatioForClustering),
                                    replace = FALSE)
    genesPosForAnalysis <- genesPosForAnalysis[-genesPosForClustering]
  }
  logThis(paste("Using", length(genesPosForClustering), "genes out of",
                numGenes, "for cluster splitting"), logLevel = 1L)

  clObjForSplit <- dropGenesCells(clObj, genes = getGenes(clObj)[-genesPosForClustering])
  clObjForSplit <- proceedToCoex(clObjForSplit, calcCoex = FALSE,
                                 cores = 6L, saveObj = FALSE, outDir = outDir)

  clusters <- splitCluster(rawData = getRawData(clObjForSplit))

  splitRes <- list("NumCells" = numCells,
                   "GDI.df"   = gdiDF,
                   "Clusters" = clusters,
                   "Genes.for.clust" = getGenes(clObjForSplit))

  # ------

  clObjForAnal <- dropGenesCells(clObj, genes = getGenes(clObj)[-genesPosForAnalysis])
  clObjForAnal <- clean(clObjForAnal)
  if (!identical(getCells(clObjForSplit), getCells(clObjForAnal))) {
    # clean might drop some new cells in different ways
    cellsToDrop <- getCells(clObjForAnal)[!getCells(clObjForAnal) %in% getCells(clObjForSplit)]
    clObjForAnal <- dropGenesCells(clObjForAnal, cells = cellsToDrop)
  }
  clObjForAnal <- proceedToCoex(clObjForAnal, calcCoex = FALSE,
                                cores = 6L, saveObj = FALSE, outDir = outDir)

  coexDF <- DEAOnClusters(clObjForAnal, clusters = clusters)

  pValDF <- pValueFromDEA(coexDF = coexDF, numCells = getNumCells(clObjForAnal),
                          adjustmentMethod = "BH")

  deaRes <- list("NumCellsForAnal" = getNumCells(clObjForAnal),
                 "Genes.for.anal" = getGenes(clObjForAnal),
                 "DEA.df" = coexDF,
                 "PVal.df" = pValDF)

  lrRes <- LR_function(clObjForAnal, clusters = clusters,
                       train.test.partition = train.test.partition)

  return(c(splitRes, deaRes, lrRes))
}


