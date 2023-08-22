genesList <- list(
  "NPGs"= 
    c("Nes", "Vim", "Sox2", "Sox1", "Notch1", "Hes1", "Hes5", "Pax6"),
  "PNGs"= 
    c("Map2", "Tubb3", "Neurod1", "Nefm", "Nefl", "Dcx", "Tbr1"),
  "hk"= 
    c("Calm1", "Cox6b1", "Ppia", "Rpl18", "Cox7c", "Erh", "H3f3a",
      "Taf1", "Taf2", "Gapdh", "Actb", "Golph3", "Zfr", "Sub1",
      "Tars", "Amacr"),
  "layers" = 
    c("Reln","Lhx5","Cux1","Satb2","Tle1","Mef2c","Rorb","Sox5","Bcl11b","Fezf2","Foxp2")
)


table.tot.hk <- NA
table.tot.neural <- NA
for (file in list.files("CoexData/")) {
  corr <- readRDS(paste0("CoexData/",file))
  if(str_detect(file,pattern = "Cotan")){
    code <- "COTAN coex"
  }else if(str_detect(file,pattern = "CorrSCT")){
    code <- "Seurat SCT corr."
  }else if(str_detect(file,pattern = "Corr")){
    code <- "Seurat corr."
  }
  
  table.hk <- as.data.frame(corr[,genesList$hk])
  table.hk$Gene1 <- rownames(table.hk)
  
  table.hk <- pivot_longer(as.data.frame(table.hk),cols = c(1:(ncol(table.hk)-1)),names_to = "Gene2")
  table.hk$Method <- code

  table.tot.hk <- rbind(table.tot.hk,table.hk)
  table.tot.hk <- table.tot.hk[! table.tot.hk$Gene1 == table.tot.hk$Gene2,]
 
  #Not hk
  table.neural <- as.data.frame(corr[c(genesList$NPGs,genesList$PNGs),c(genesList$NPGs,genesList$PNGs)])
  table.neural$Gene1 <- rownames(table.neural)
  
  table.neural <- pivot_longer(as.data.frame(table.neural),cols = c(1:(ncol(table.neural)-1)),names_to = "Gene2")
  table.neural$Method <- code
  
  table.tot.neural <- rbind(table.tot.neural,table.neural)
  table.tot.neural <- table.tot.neural[! table.tot.neural$Gene1 == table.tot.neural$Gene2,]
  
   
}




library(ggplot2)
table.tot.hk <- table.tot.hk[2:nrow(table.tot.hk),]
table.tot.hk$GeneType <- "Constitutive" 
table.tot.neural <- table.tot.neural[2:nrow(table.tot.neural),]
table.tot.neural$GeneType <- "Neural"

table.tot <- rbind(table.tot.hk,table.tot.neural) 
table.tot <- table.tot[!table.tot$Gene1 == table.tot$Gene2,]

identical(table.tot[table.tot$Method == "COTAN coex",]$Gene1,
          table.tot[table.tot$Method == "Seurat corr.",]$Gene1 )
identical(table.tot[table.tot$Method == "COTAN coex",]$Gene2,
          table.tot[table.tot$Method == "Seurat corr.",]$Gene2 )
cotanVSSeurat <- merge(table.tot[table.tot$Method == "COTAN coex",],
                       table.tot[table.tot$Method == "Seurat corr.",],
                       by = c("Gene1","Gene2","GeneType"))
cotanVSSeurat$DifAbs <- cotanVSSeurat$ValAbs.x - cotanVSSeurat$ValAbs.y 
cotanVSSeurat.HK <- cotanVSSeurat[cotanVSSeurat$GeneType == "Constitutive",]
cotanVSSeurat.HK$DifScaled <- scale(cotanVSSeurat.HK$DifAbs, center = F)

cotanVSSeurat.Neuro <- cotanVSSeurat[cotanVSSeurat$GeneType == "Neural",]
cotanVSSeurat.Neuro$DifScaled <- scale(cotanVSSeurat.Neuro$DifAbs, center = F)

pl1 <- ggplot(cotanVSSeurat.HK, aes(y=DifAbs,x=GeneType)) + geom_boxplot()+
  geom_point(position = "jitter", size= 0.1,alpha = 0.3)+ylim(-0.4,0.3)
pl2 <- ggplot(cotanVSSeurat.Neuro, aes(y=DifAbs,x=GeneType)) + geom_boxplot()+
  geom_point(position = "jitter", size= 0.1,alpha = 0.3)+ylim(-0.4,0.3)
plot_grid(pl1,pl2)

library(gghalves)
library(ggstatsplot)

table.tot$SqrValAbs <- table.tot$ValAbs**2


ggplot(table.tot,aes(x=Method,y=ValAbs, fill=Method)) +
geom_half_violin(alpha=0.8) +
  geom_point(position = "jitter", size= 0.1,alpha = 0.3)+ 
  geom_half_boxplot(width=0.1, alpha=0.8,side = "r") +
  theme_bw()+   theme(legend.position="none")+
  facet_grid(GeneType ~ . , space = "free_y", scales = "free")+
  scale_fill_manual(values = c("#8856A7","#EDF8FB","#B3CDE3" ))

ggbetweenstats(
  data = table.tot.neural[!table.tot.neural$Method == "Seurat corr.",],
  x = Method,
  y = value
)

ggplot(table.tot,aes(x=Method,y=log(ValAbs+0.01), fill=GeneType)) +
  geom_boxplot(alpha=0.8) +
  geom_hline(yintercept = -3.39)+
  geom_hline(yintercept = -1.86)+
  #geom_point(position = "jitter", size= 0.1,alpha = 0.3)+ 
  #geom_half_boxplot(width=0.1, alpha=0.8,side = "r") +
  theme_bw()+   #theme(legend.position="none")+
  scale_fill_manual(values = c("#8856A7","#EDF8FB","#B3CDE3" ))

ggplot(table.tot, aes(x=Method,y=(abs(value)+0.01), fill=GeneType)) +
  introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  # scale_x_discrete(name = "Condition", labels = c("Non-word", "Word")) +
  # scale_y_continuous(name = "Reaction time (ms)",
  #                    breaks = seq(200, 800, 100), 
  #                    limits = c(200, 800)) +
  scale_fill_brewer(palette = "Dark2", name = "Language group") +
  theme_minimal()

table.tot$LogValAbs <- log(abs(table.tot$value)+0.01)
table.tot$ValAbs <- abs(table.tot$value)

seur <- ggbetweenstats(
  data = table.tot[table.tot$Method == "Seurat corr.",],
  x = GeneType,type = "p",
  y = LogValAbs,mean.plotting = TRUE,mean.ci = TRUE, results.subtitle = F,
  title = "Seurat corr"
)

SCT <- ggbetweenstats(
  data = table.tot[table.tot$Method == "Seurat SCT corr.",],
  x = GeneType,type = "p",
  y = LogValAbs,mean.plotting = TRUE,mean.ci = TRUE, results.subtitle = F,
  title = "Seurat corr"
)

cotan <- ggbetweenstats(
  data = table.tot[table.tot$Method == "COTAN coex",],
  x = GeneType,type = "p",
  y = LogValAbs,results.subtitle = FALSE,title = "COTAN coex"
    
)

library("gridExtra")
library(cowplot)
grid.arrange(cotan,seur,SCT,ncol=2, nrow =2)
plot_grid(cotan,SCT, )

ggdraw() +
  cotan + SCT + seur 
  
table.tot.neural$ValAbs <- abs(table.tot.neural$value)
ggbetweenstats(
  data = table.tot.neural,
  x = Method,
  y = LogValAbs,
  type = "bayes",
  title = "Bayesian Test",
  package = "ggsci",
  palette = "nrc_npg"
)


grouped_gghistostats(
  data            = table.tot.neural,
  x               = ValAbs,
  test.value      = 0,
  grouping.var    = Method,
  plotgrid.args   = list(nrow = 1),
  annotation.args = list(tag_levels = "i")
)
