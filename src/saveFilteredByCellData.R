
library(Seurat)
library(stringr)

for (dataset in c("PBMC1", "PBMC2", "PBMC3", "PBMC4")) {
print(dataset)
raw <- Read10X(paste0("Data/",dataset,"/raw/10X/"))

cells.filter <- read.csv(paste0("Data/",dataset,"/filtered/10X/barcodes.tsv"),header = F)

cells.filter <- paste0(cells.filter$V1,"-1")

raw <- raw$`Gene Expression`[,cells.filter]

write.csv(raw,paste0("Data/",dataset,"/filtered/Filtered_only_cells.csv"))
}

