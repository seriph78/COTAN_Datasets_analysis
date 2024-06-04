correlation_pvalues <- function(data,int.genes, n.cells){
  data <- t(as.matrix(data)[rownames(data) %in% int.genes,])
  
  data.cor <- fastCor(data,upperTri = T,verbose = T,optBLAS = T)
  
  data.cor <- Matrix::forceSymmetric(data.cor, uplo = "L")
  diag(data.cor) <- 1
  
  temp <- ((data.cor)**2)*n.cells
  p_values <- pchisq(as.matrix(temp), df = 1L, lower.tail = FALSE)
  #diag(p_values) <- 1.0
  
  return(list("data.cor"= data.cor,"p_values"=p_values))
  
}

correlation_plot <- function(data.cor.big, genesList, title){
  #check the presence fo all genes and drop eventually absences
  for (na in names(genesList)) {
  if(!all(genesList[[na]] %in% rownames(data.cor.big))){
      genesList[[na]] <- genesList[[na]][-which(!genesList[[na]] %in% rownames(data.cor.big))]
    }
  }
  
  data.cor <- data.cor.big[c(genesList$NPGs,
                             genesList$hk,
                             genesList$PNGs),
                           c(genesList$NPGs,
                             genesList$hk,
                             genesList$PNGs)]
  
  diag(data.cor) <- 0
  
  f1 = colorRamp2(seq(-0.5,0.5, length = 3), c("#DC0000B2", "white","#3C5488B2" ))
  
  split.genes <- factor(c(rep("NPGs",length(genesList[["NPGs"]])),
                          rep("HK",length(genesList[["hk"]])),
                          rep("PNGs",length(genesList[["PNGs"]]))),
                        levels = c("NPGs","HK","PNGs"))
  
  htmp <- Heatmap(as.matrix(data.cor),
                  #width = ncol(seurat.corMat)*unit(2.5, "mm"), 
                  height = nrow(data.cor)*unit(3, "mm"),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  col = f1,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 11),
                  column_names_gp  = gpar(fontsize = 11),
                  column_split = split.genes,
                  row_split = split.genes,
                  cluster_row_slices = FALSE, 
                  cluster_column_slices = FALSE,
                  heatmap_legend_param = list(
                    title = title, at = c(-0.5, 0, 0.5),direction = "horizontal",
                    labels = c("-0.5", "0", "0.5")
                  )
  )
  return(htmp)
}
