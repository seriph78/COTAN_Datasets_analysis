library(Seurat)

seurat.clustering.COTAN.obj.creation <- function(raw.data,
                                                 outDir,
                                                 resolution = 0.5,dims = 20,
                                                 GEO,
                                                 sampleCondition,
                                                 sequencingMethod){
  
  pbmc <- CreateSeuratObject(counts = raw.data, project = "temp", min.cells = 3, min.features = 100)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- FindNeighbors(pbmc, dims = 1:dims)
  pbmc <- FindClusters(pbmc, resolution = resolution)
  
  dataset <- automaticCOTANObjectCreation(raw = raw.data,
                                          GEO = GEO,
                                          sequencingMethod = sequencingMethod,
                                          sampleCondition = sampleCondition,cores = 5,
                                          outDir = outDir
                                          )
  dataset <- clean(dataset)
  
  dataset <- addClusterization(dataset,clName = paste0("FromSeurat",resolution,dims),
                               clusters = pbmc$seurat_clusters)
  
  return(dataset)
}


datasets <- c("Data/CD14Cleaned/CD14_Monocytes.cotan.RDS",
              "Data/Yuzwa_MouseCortex/CorticalCells_GSM2861511_E135.cotan.RDS",
              "Data/MouseCortexFromLoom/e15.5_ForebrainDorsal.cotan.RDS",
              "Data/PBMC3/filtered/PBMC3.cotan.RDS")

outDir <- "Data/DataForGDITh/"

for (d in "Data/MouseCortexFromLoom/e15.0_ForebrainDorsal.cotan.RDS") {
  data <- readRDS(d)
  
  raw <- getRawData(data)
  
  resolution = 1.5
  dims = 20
  GEO = getMetadataDataset(data)[1,2]
  sampleCondition = getMetadataDataset(data)[3,2]
  sequencingMethod = getMetadataDataset(data)[2,2]
  
  rm(data)
  
  obj <- seurat.clustering.COTAN.obj.creation(raw.data = raw,outDir = outDir,
                                              resolution = resolution,
                                              dims = dims,
                                              GEO =GEO,
                                              sampleCondition=sampleCondition,
                                              sequencingMethod=sequencingMethod)
  saveRDS(obj,paste0(outDir,sampleCondition,".COTAN.RDS"))
  
}

devtools::load_all("../COTAN/")
options(parallelly.fork.enable = TRUE)

file = "e15.0_ForebrainDorsal.cotan.RDS"

for (file in list.files("Data/DataForGDITh/", pattern = ".cotan.RDS")) {
  
  obj <- readRDS(paste0(outDir,file))
  
  cond <- getMetadataElement(obj, datasetTags()[["cond"]])
  
  outDirCond <- file.path(outDir, cond)
  if (!dir.exists(outDirCond)) {
    dir.create(outDirCond)
  }
  
  clList <- toClustersList(getClusters(obj,clName = "FromSeurat1.520"))
  for (cl in names(clList)) {
    print(paste0("Cluster ", cl))
    
    checkClusterUniformity(obj,clusterName = cl,cells = clList[[cl]],GDIThreshold = 1.37,
                           ratioAboveThreshold = 0.005,cores = 5,saveObj = TRUE,
                           outDir = outDirCond)
    
  }
  

}

clusterSeurat <- function(raw.matrix, genesToKeep=rownames(raw.matrix), res=1){
  raw.matrix <- raw.matrix[genesToKeep,]
  pbmc <- CreateSeuratObject(counts = raw.matrix, project = "temp", min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- FindNeighbors(pbmc, dims = 1:20)
  pbmc <- FindClusters(pbmc, resolution = res)
  
  return(pbmc)
}

library(ggplot2)
plots.coex <- function(obj){
  coex <- getGenesCoex(obj)
  coex <- coex[lower.tri(coex, diag = FALSE)]
  coex.array <- as.vector(coex)
  coex.array <- as.data.frame(coex.array)
  ggplot(coex.array, aes(coex.array)) + stat_ecdf(geom = "point")
  }

cells <- names(getClusters(CD14_Monocytes.cotan)[getClusters(CD14_Monocytes.cotan)== 1])

cells <- names(getClusters(PBMC3.cotan)[getClusters(PBMC3.cotan)== 0])
obj_cl0_PBMC <- automaticCOTANObjectCreation(raw = getRawData(PBMC3.cotan)[,cells],
                                        GEO = "",sequencingMethod = "" ,sampleCondition = "PBMC_Cl_0",
                                        calcCoex = TRUE,cores = 5,saveObj = FALSE)

GDI_cl0_10 <- calculateGDI(obj_cl0_PBMC,rowsFraction = 0.1)
GDI_cl0_10$From <- "cl0_10_PBMC"

GDI_cl0_5 <- calculateGDI(obj_cl0_PBMC,rowsFraction = 0.05)
GDI_cl0_5$From <- "cl0_5_PBMC"
GDI <- rbind(GDI, GDI_cl0_10,GDI_cl0_5)

obj_cl0 <- automaticCOTANObjectCreation(raw = getRawData(CD14_Monocytes.cotan)[,cells],
                                    GEO = "",sequencingMethod = "" ,sampleCondition = "CD14Cl_0",
                                    calcCoex = TRUE,cores = 5,saveObj = FALSE)

obj_cl1 <- automaticCOTANObjectCreation(raw = getRawData(CD14_Monocytes.cotan)[,cells],
                                        GEO = "",sequencingMethod = "" ,sampleCondition = "CD14Cl_1",
                                        calcCoex = TRUE,cores = 5,saveObj = FALSE)

GDI_cl0_10 <- calculateGDI(obj_cl0,rowsFraction = 0.1)
GDI_cl0_10$From <- "cl0_10_CD14"
GDI_cl1_10 <- calculateGDI(obj_cl1,rowsFraction = 0.1)
GDI_cl1_10$From <- "cl1_10_CD14"

GDI_cl0_5 <- calculateGDI(obj_cl0,rowsFraction = 0.05)
GDI_cl0_5$From <- "cl0_5_CD14"

GDI_cl1_5 <- calculateGDI(obj_cl1,rowsFraction = 0.05)
GDI_cl1_5$From <- "cl1_5_CD14"
GDI <- rbind(GDI_cl0_10,GDI_cl0_5,GDI_cl1_10,GDI_cl1_5)


seurat.obj <- clusterSeurat(getRawData(obj),res = 0.7)

obj <- addClusterization(obj,clName = "FromSeurat_0.7",clusters = seurat.obj$seurat_clusters)


