
raw <- Read10X("Data/CD14Cleaned/OrigialDatahg19/")

obj <- COTAN(raw = raw)
obj <- initializeMetaDataset(obj,
                             GEO = "_",
                             sequencingMethod = "10X",
                             sampleCondition = "CD14")

obj <- clean(obj)

c(pcaCellsPlot, pcaCellsData, genesPlot,
  UDEPlot, nuPlot, zoomedNuPlot) %<-% cleanPlots(obj)
options(parallelly.fork.enable = TRUE)

cells_to_rem <- rownames(pcaCellsData)[pcaCellsData[["groups"]] == "B"]

obj <- dropGenesCells(obj, cells = cells_to_rem)

obj <- addElementToMetaDataset(obj, "Num drop B group", 1)

obj <- clean(obj)
c(pcaCellsPlot, pcaCellsData, genesPlot,
  UDEPlot, nuPlot, zoomedNuPlot) %<-% cleanPlots(obj)

advChecker <- new("AdvancedGDIUniformityCheck")
obj <- proceedToCoex(obj, calcCoex = FALSE,
                     optimizeForSpeed = TRUE, cores = 3L, deviceStr = "cuda",
                     saveObj = FALSE)

clusters <- cellsUniformClustering(obj,checker = advChecker,cores = 3,saveObj = TRUE,initialResolution = 0.5, outDir =  "Temp/")

obj <- addClusterization(obj, clName = "SplitCD14Monocytes_new",clusters = clusters$clusters)
table(clusters$clusters)

obj <- dropGenesCells(obj,cells = names(clusters$clusters[clusters$clusters == "-1"]))
obj <- clean(obj)
obj <- proceedToCoex(obj,calcCoex = F,cores = 3)

clusters.merged <- mergeUniformCellsClusters(objCOTAN = obj,checker = advChecker,
                                             clusters = getClusterizationData(obj,clName = "MergedCD14Monocytes_new")[[1]],
                                             cores = 3,saveObj = T,batchSize = 7,
                                             outDir =  "Temp/")

obj <- addClusterization(obj, clName = "MergedCD14Monocytes_new",clusters = clusters.merged$clusters)


SeuratObjAz_orig <- CreateSeuratObject(counts = getRawData(obj),
                                     assay = "RNA",
                                     meta.data = getMetadataCells(obj)[,c("CL_SplitCD14Monocytes_new","CL_MergedCD14Monocytes_new")],
                                     project = "CD14",
                                     min.cells = 3,
                                     min.features = 200)

SeuratObjAz_orig <- RunAzimuth(SeuratObjAz_orig, reference = "pbmcref")

SeuratObjAz_orig <- NormalizeData(SeuratObjAz_orig)
SeuratObjAz_orig <- ScaleData(SeuratObjAz_orig)
SeuratObjAz_orig <- FindVariableFeatures(SeuratObjAz_orig, selection.method = "vst", nfeatures = 2000)

SeuratObjAz_orig <- RunPCA(SeuratObjAz_orig, verbose = FALSE)
SeuratObjAz_orig <- RunUMAP(SeuratObjAz_orig, dims = 1:20,
                       verbose = FALSE)
SeuratObjAz_orig <- SetIdent(SeuratObjAz_orig,value = "CL_MergedCD14Monocytes_new")

DimPlot(SeuratObjAz_orig,label = T)

table(SeuratObjAz_orig@meta.data$CL_MergedCD14Monocytes_new,SeuratObjAz_orig@meta.data$predicted.celltype.l1)


objClMarckers <- findClustersMarkers(obj,clName = "MergedCD14Monocytes_new" ,n = 1000,adjustmentMethod = "BH")
head(objClMarckers)
write.csv(objClMarckers,file = "Temp/ClMarkers.csv")


for (cl in unique(objClMarckers$CL)) {
  subset1 <- objClMarckers[objClMarckers$CL == cl &
                            objClMarckers$adjPVal < 0.01 &
                            objClMarckers$logFoldCh > 0,]$Gene
  
  subset2 <- objClMarckers[objClMarckers$CL == cl 
                             ,c("Gene","logFoldCh") ]
  
  write.csv(subset1,file = paste0("Temp/CD14_cl",cl,"_ORA.txt"),row.names = F,quote = F)
  
  write.csv(subset2,file = paste0("Temp/CD14_cl",cl,"_GSEA.csv"),row.names = F)
  
}

write.csv(getGenes(obj),file = paste0("Temp/CD14_background_ORA.txt"),row.names = F,quote = F)
