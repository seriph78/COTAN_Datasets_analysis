e15 <- readRDS()

radial.glia.cells <- rownames(e15@metaCells[e15@metaCells$Class == "Radial glia",])
automaticCOTANObjectCreation(raw = getRawData(e15)[radial.glia.cells],
                             GEO = getMetadataDataset(e15)[1,2],
                             sequencingMethod = getMetadataDataset(e15)[2,2],
                             sampleCondition = "E15.0_RadialGlia",
                             cores = 15,saveObj = T,
                             outDir = "Data/MouseCortexFromLoom/RadialGliaE15.0/")