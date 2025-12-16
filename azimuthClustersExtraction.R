library(stringr)

for (fl in list.files("Data/AzimuthLabellingMatrixes/")) {
  print(fl)
  sample.name <- str_split(fl,pattern = "_",simplify = T)[1] 
  level <- paste(str_split(fl,pattern = "_",simplify = T)[c(2,3)],collapse = "_")
  
  if (sample.name %in% c("PBMC1","PBMC2","PBMC3","PBMC4")) {
    data.in <- read.csv(paste0("Data/AzimuthLabellingMatrixes/",fl),row.names = 1)
    
    labels <- as.data.frame(colnames(data.in)[max.col(data.in)])
    labels$cell <- rownames(data.in) 
    colnames(labels)[1] <-"cluster"
    labels <- labels[,c("cell","cluster")]
    if (!exists(paste0("Data/",sample.name,"/azimuth/default/"))) {
      dir.create(paste0("Data/",sample.name,"/azimuth/default/"))
    }
    write.csv(labels,file = paste0("Data/",sample.name,"/azimuth/default/","clustering_labels.csv"),row.names = F)
    write.csv(labels,file = paste0("Data/",sample.name,"/azimuth/default/",level,".csv"))
  }
  
  
}

advChecker <- new("AdvancedGDIUniformityCheck")
advChecker <- shiftCheckerThresholds(advChecker,shift = 0.3)
clMergeMore <- mergeUniformCellsClusters(PBMC1.cotan,checkers = "AdvancedGDIUniformityCheck",cores = 3,saveObj = F)

