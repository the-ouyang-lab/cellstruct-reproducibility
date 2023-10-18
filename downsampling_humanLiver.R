#!/usr/bin/env Rscript

# #running this script with following command in bash:
# for x in 5000 8000 10000 20000 30000 40000 50000 60000 70000; do for y in {1..10};
# do Rscript downsampling_humanLiver.R $x $y 20 ;done;done

devtools::load_all(".")
library(cellstruct)

loc <- "~/"
seuName <- "humanLiver"
seu <- readRDS(paste0(loc,"/",seuName,".rds"))
clusterVar <- "celltype.l2"
ref.proj <- "refDR"

#for very large object, just to save space
seu@assays[["RNA"]]@counts <- seu@assays[["RNA"]]@counts[1:20,]
seu@assays[["RNA"]]@data <- seu@assays[["RNA"]]@data[1:20,]
gc()

args = commandArgs(trailingOnly=TRUE)
cells_size <- as.numeric(args[1])
seed <- as.numeric(args[2])
num_cores <- as.numeric(args[3])

set.seed(seed)
new_seu <- seu[, sample(colnames(seu), size = cells_size, replace=F)]
rm(seu)
pca <- new_seu@reductions[[ref.proj]]@cell.embeddings
clusID <- data.frame(FetchData(object = new_seu, vars = clusterVar))
write.table(cbind(pca,clusID), file = paste0(loc,seuName,"/downsampling/pca_embed/cells_",as.character(cells_size/1000),"K.",seed,".txt"), sep = "\t")

metric <- c("euclidean","cosine")
umap_neighbors <- c(15,30,50)
min_dist <- c(0.05, 0.1, 0.3, 0.5)

for (x in metric){
  for (y in umap_neighbors){
    for (z in min_dist){
      singleTuneUMAP(new_seu, x, y, z, ref.proj,
                     paste0(loc,seuName,"/downsampling/tuneUMAP/cells_",as.character(cells_size/1000),"K.",seed),
                     clusterVar, nCores = num_cores)
    }
  }
}
