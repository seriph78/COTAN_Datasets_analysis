
library(tree)

allSummariesDF <- read.csv("Data/DataForGDITh/AllAnalysisSummaries.csv",row.names = 1)

allSummariesDF$ClusterWasUniformLabel <- factor(allSummariesDF$ClusterWasUniformLabel)

set.seed(42)

train <- sample(1:nrow(allSummariesDF), nrow(allSummariesDF)*0.8)

#allSummariesDF$MisclassificationError
#allSummariesDF$AreaUnderROCCurve
#allSummariesDF$ClusterWasUniform
#allSummariesDF$NumGenesWithPValBelow5Perc

colnames(allSummariesDF)

#allSummariesDF_Backup <- allSummariesDF
trainSummariesDF <- trainSummariesDF_Backup

sum(is.na(trainSummariesDF$ClusterWasUniform))

collapsedSummariesDF <- trainSummariesDF
collapsedSummariesDF$ClusterWasUniform <-
  ifelse(trainSummariesDF$ClusterWasUniform == 0.5, NA,
         round(roundSummariesDF$ClusterWasUniform))

collapsedSummariesDF <- collapsedSummariesDF[!is.na(collapsedSummariesDF$ClusterWasUniform), ]

allSummariesDF <- allSummariesDF[-1:-9, ]
collapsedSummariesDF <- collapsedSummariesDF[-1:-9, ]



#tree.GDI.WithCD14.MiscErr     <- tree(formula = (1 * (MisclassificationError)) ~
#tree.GDI.WithCD14.AreaROC     <- tree(formula = (2 * (1-AreaUnderROCCurve)) ~
#tree.GDI.WithCD14.ClIsUni     <- tree(formula = (1 * (ClusterWasUniform)) ~
#tree.GDI.WithCD14.ClIsUniC    <- tree(formula = (1 * (ClusterWasUniform)) ~
#tree.GDI.WithCD14.NumSignG    <- tree(formula = (1 * (NumGenesWithPValBelow5Perc)) ~
#tree.GDI.WithoutCD14.MiscErr  <- tree(formula = (1 * (MisclassificationError)) ~
#tree.GDI.WithoutCD14.AreaROC  <- tree(formula = (2 * (1-AreaUnderROCCurve)) ~
#tree.GDI.WithoutCD14.ClIsUni  <- tree(formula = (1 * (ClusterWasUniform)) ~
#tree.GDI.WithoutCD14.ClIsUniC <- tree(formula = (1 * (ClusterWasUniform)) ~
#tree.GDI.WithoutCD14.NumSignG <- tree(formula = (1 * (NumGenesWithPValBelow5Perc)) ~
#                   NumClCells +
tree.GDI <- tree(formula = ClusterWasUniformLabel ~  
                   NumGenesWithGDIAbove_1.50 +
                   PercGenesWithGDIAbove_1.33 +
                   PercGenesWithGDIAbove_1.37 +
                   PercGenesWithGDIAbove_1.40 +
                   PercGenesWithGDIAbove_1.43 +
                   PercGenesWithGDIAbove_1.46 +
                   PercGenesWithGDIAbove_1.50 +
                   GDIQuantile_95    +
                   GDIQuantile_97    +
                   GDIQuantile_98    +
                   GDIQuantile_99    +
                   GDIQuantile_99.5  +
                   GDIQuantile_99.75 +
                   GDIQuantile_99.9  +
                   Top1stGDIGene  +
                   Top2ndGDIGene  +
                   Top5thGDIGene  +
                   Top10thGDIGene +
                   Top20thGDIGene +
                   Top30thGDIGene +
                   #NumClCells +
                   Top50thGDIGene,
                 data = allSummariesDF,
                 subset = train)
#                 data = collapsedSummariesDF)

#tree.GDI <- tree.GDI.WithCD14.MiscErr
#tree.GDI <- tree.GDI.WithCD14.AreaROC
#tree.GDI <- tree.GDI.WithCD14.ClIsUni
#tree.GDI <- tree.GDI.WithCD14.ClIsUniC
#tree.GDI <- tree.GDI.WithCD14.NumSignG
#
#tree.GDI <- tree.GDI.WithoutCD14.MiscErr
#tree.GDI <- tree.GDI.WithoutCD14.AreaROC
#tree.GDI <- tree.GDI.WithoutCD14.ClIsUni
#tree.GDI <- tree.GDI.WithoutCD14.ClIsUniC
#tree.GDI <- tree.GDI.WithoutCD14.NumSignG

summary(tree.GDI)

plot(tree.GDI)
text(tree.GDI, pos = NULL)

cv.tree.GDI <- cv.tree(tree.GDI)
plot(cv.tree.GDI$size, cv.tree.GDI$dev, type = "b")

prune.GDI.tree <- prune.tree(tree.GDI, best = 4)
plot(prune.GDI.tree)
text(prune.GDI.tree, pretty = 0)

saveRDS(prune.GDI.tree,"Data/DataForGDITh/PrunedBestModel.RDS")


GDI.test <- allSummariesDF[-train, ]
tree.pred <- predict(prune.GDI.tree, GDI.test,
                       type = "class")
table(tree.pred, GDI.test$ClusterWasUniformLabel)
(15+3)/(15+1+2+3)

tree_frame <- prune.GDI.tree$frame
leaf_nodes <- tree_frame[tree_frame$var == "<leaf>", ]

plot(prune.GDI.tree)
text(prune.GDI.tree, pretty = 0, cex = 0.8)

for (i in 1:nrow(leaf_nodes)) {
  node_index <- as.numeric(row.names(leaf_nodes))[i]
  node_label <- paste("n=", leaf_nodes$n[i], sep="")
  
  # Manually add text to the plot at the position of the leaf nodes
  text(prune.GDI.tree$frame$yval[node_index], node_index, labels = node_label, pos = 4, cex = 0.8)
}


yhat <- predict(tree.GDI, newdata = allSummariesDF[-train, ])
GDI.test <- allSummariesDF[-train, "ClusterWasUniformLabel"]
plot(yhat, GDI.test)
abline(0, 1)
mean((yhat - GDI.test)^2)

yhat <- predict(tree.GDI, newdata = allSummariesDF[-train, ])
GDI.test <- allSummariesDF[-train, "ClusterWasUniform"]
plot(yhat, GDI.test)
abline(0, 1)
mean((yhat - GDI.test)^2)


yhat <- predict(prune.GDI.tree, newdata = allSummariesDF[-train, ])
GDI.test <- allSummariesDF[-train, "MisclassificationError"]
plot(yhat, GDI.test)
abline(0, 1)
mean((yhat - GDI.test)^2)

######## Random forest
library(randomForest)

rf.GDI <- randomForest(ClusterWasUniform ~ 
                       NumGenesWithGDIAbove_1.50 +
                         PercGenesWithGDIAbove_1.33 +
                         PercGenesWithGDIAbove_1.37 +
                         PercGenesWithGDIAbove_1.40 +
                         PercGenesWithGDIAbove_1.43 +
                         PercGenesWithGDIAbove_1.46 +
                         PercGenesWithGDIAbove_1.50 +
                         GDIQuantile_95    +
                         GDIQuantile_97    +
                         GDIQuantile_98    +
                         GDIQuantile_99    +
                         GDIQuantile_99.5  +
                         GDIQuantile_99.75 +
                         GDIQuantile_99.9  +
                         Top1stGDIGene  +
                         Top2ndGDIGene  +
                         Top5thGDIGene  +
                         Top10thGDIGene +
                         Top20thGDIGene +
                         Top30thGDIGene +
                         #NumClCells +
                         Top50thGDIGene, 
                          data = allSummariesDF,
                          subset = train, mtry = 6, importance = TRUE)



yhat.rf <- predict(rf.GDI, newdata = allSummariesDF[-train, ])
mean((yhat.rf - boston.test)^2)