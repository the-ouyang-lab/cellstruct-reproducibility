library(Seurat)
library(SeuratDisk)

loc <- "~/"
dataset <- "wilk_2020_processed.HDF5"
seu <- LoadH5Seurat(paste0(loc,dataset))

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu))
ElbowPlot(seu, ndims=50)
num_PC <- 50
seu <- RunUMAP(seu, reduction = "pca", dims = 1:num_PC)
seu <- RunTSNE(seu, reduction = "pca", dims = 1:num_PC)

library(reticulate)
library(SeuratWrappers)
os <- import("os")
py_run_string("r.os.environ['OMP_NUM_THREADS'] = '20'")
# Setup python anndata for downstream DFL/DiffMap/PAGA
sc <- import("scanpy", convert = FALSE)
ad <- import("anndata", convert = FALSE)

adata <- sc$AnnData(
  X   = np_array(t(GetAssayData(seu)[VariableFeatures(seu),]), dtype="float32"),
  obs = seu@meta.data[, c("seurat_clusters")],
  var = data.frame(geneName = VariableFeatures(seu)))
adata$obsm$update(X_pca = Embeddings(seu, "pca"))
adata$obsm$update(X_umap = Embeddings(seu, "umap"))
adata$var_names <- VariableFeatures(seu)
adata$obs_names <- colnames(seu)
sc$pp$neighbors(adata, n_neighbors = as.integer(30), n_pcs = as.integer(num_PC))

# FDL
sc$tl$draw_graph(adata, layout = "fa", init_pos = "X_umap")
oupDR <- py_to_r(adata$obsm['X_draw_graph_fa'])
rownames(oupDR) <- colnames(seu)
colnames(oupDR) <- c("FDL_1","FDL_2")
oupDR = oupDR / 10^(floor(log10(diff(range(oupDR))))-1)
seu[["fdl"]] <- CreateDimReducObject(embeddings = oupDR, key = "FDL_", assay = DefaultAssay(seu))

ref.proj <- "pca"
target.proj <- c("umap","tsne","fdl")
clusterVar <- "cell.type.fine"
seuName <- paste0(loc,"wilk_2020")
seu <- run_cellstruct(seu, seuName, ref.proj, target.proj, clusterVar, nCores = 50)

saveRDS(seu, paste0(seuName,".rds"))


#for Wilk2020
seu <- readRDS(paste0(seuName,".rds"))
set.seed(42) #start here
all.coi <- colnames(seu)[which(seu$cell.type.fine=="CD14 Monocyte"&seu$Status=="COVID")]
coi <- all.coi[sample(c(1:length(all.coi)), size = 1657)] # randomly selected 20% cells (saved in "tmp20Percent" metadata as "COVID-CD14Mono20Percent")
heatmapGS(seu, ref.proj, target.proj, coi = coi) #running it inwardly within the function
library(ggsci)
plot.embed1 #Figure S8B
clrs <- c("darkred","orchid","tomato","darkorange","lightskyblue","dodgerblue","slateblue","mediumblue","navy",
          "mediumorchid4","hotpink","cyan","plum1","purple","indianred","green","red","grey80","peru","lightcoral")
DimPlot(seu, reduction = "umap", group.by = clusterVar, cols = clrs, shuffle = T)#Figure 2A left
DimPlot(seu, reduction = "umap", group.by = "Status", shuffle = T) #Figure 2A middle

# #metadata "randomCells" - randomly selected 1000 cells, along with IDs that are interested
# seu$randomCells <- rep(NA,ncol(seu))
# seu$randomCells[randomCell] <- "Selected"
# seu$randomCells[randomCell[which(randomCell.celltype=="CD14 Monocyte")]] <- "Other CD14 Mono"
# seu$randomCells[randomCell[randOrder[726:772]]] <- "Interested CD14 Mono"


default.full <- seu
seu <- subset(default.full, subset = cell.type.fine == "CD14 Monocyte" & Status == "COVID")
#need to rereun runUMAP, runTSNE, and get FDL again
ref.proj <- "pca"
target.proj <- c("umap","tsne","fdl")
ref <- seu@reductions[[ref.proj]]@cell.embeddings
for(proj in target.proj){
  data <- seu@reductions[[proj]]@cell.embeddings
  seu <- AddMetaData(object = seu,
                     metadata = calcGS(ref, data, nCores = 50),
                     col.name = paste0(toupper(proj),"_gs"))
}
seu[["defaultUMAP"]] <- CreateDimReducObject(embeddings = seu@reductions[["umap"]]@cell.embeddings, key = "defaultUMAP_", assay = DefaultAssay(seu))
#tuneUMAP(seu, seuName, "pca", paste0(loc,"/wilk_2020_figures/tuneUMAP/full"),clusterVar, nCores = 50)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:num_PC, metric = "euclidean", n.neighbors = 30, min.dist = 0.5) # tuned UMAP hyperparameters
#umap - tunedUMAP, defaultUMAP - umap with default parameters of runUMAP
seu$defaultUMAP_gs <- seu$UMAP_gs
data <- seu@reductions[["umap"]]@cell.embeddings
seu$UMAP_gs <- calcGS(ref, data, nCores = 50)

coi <- colnames(seu)[which(!is.na(seu$tmp20Percent))] #same 20% of COVID-CD14Mono
heatmapGS(seu, ref.proj, target.proj, coi = coi) #running it inwardly within the function
plot.embed1 + Heatmap(exprs$tmp20Percent.3, name = "Group",col = pal_jco("default")(4), width = unit(5, "mm")) +
  Heatmap(exprs$admission.ventilated, name = "Severity",col = c("purple3","orange","springgreen3"), width = unit(5, "mm")) #Figure 2C

# #metadata "tmp20Percent.3" - randomly selected 1657 CD14 Monocyte, COVID cells and labeled them
# seu$tmp20Percent.3[which(colnames(seu) %in% coi[cellOrder[1:94]])] <- "GroupB"
# seu$tmp20Percent.3[which(colnames(seu) %in% coi[cellOrder[95:671]])] <- "GroupA"
# seu$tmp20Percent.3[which(colnames(seu) %in% coi[cellOrder[672:825]])] <- "GroupC"
# seu$tmp20Percent.3[which(colnames(seu) %in% coi[cellOrder[826:1657]])] <- "GroupD"


plot.list <- list()
projections <- c("defaultUMAP","umap")
for (p in 1:2){
  plot.list[[p]] <- DimPlot(seu, reduction = projections[p], group.by = "tmp20Percent.3", order = T, cols = pal_jco("default")(4), na.value = "grey90")
}
plot_grid(plotlist = plot.list, nrow = 1) # Figure 2D left, middle
DimPlot(seu, reduction = "umap", group.by = "admission.ventilated",
        shuffle = T, cols = c("purple3","orange","springgreen3")) #Figure 2D right

plot.list <- list()
exprs <- data.frame(FetchData(seu, vars = c("defaultUMAP_gs","UMAP_gs")))
for (p in 1:2){
  plot.list[[p]] <- FeaturePlot(seu, reduction = projections[p], features = colnames(exprs)[p], order = T) +
    ggtitle(paste0(colnames(exprs)[p],": ",round(mean(exprs[,p]),4))) +
    scale_color_distiller(palette = "RdYlBu", limits = c(0,1))
}
plot_grid(plotlist = plot.list, nrow = 1) #Figure 2B
saveRDS(seu, paste0(seuName,"/subsetCD14MonoCOVID.rds"))


#statistics on 2x2 contigency table
df <- data.frame("GroupA" = c(77,93,407), "GroupB" = c(0,53,41), "GroupC" = c(0,84,70),
                 "GroupD" = c(321,219,292), row.names = c("Floor.NonVent","ICU.NonVent","ICU.Vent"))
library(stats)
chisq <- chisq.test(df)
library(corrplot)
corrplot(chisq$residuals, method = "color", is.cor = FALSE, addCoef.col = 'grey50',
         cl.align.text = 'l', cl.offset = 0.25, col = rev(COL2('RdBu', 200))) #Figure 2E


#differential expression analysis
Idents(seu) <- seu$admission.ventilated
CD14MonoCOVID1657.admVent.markers <- FindAllMarkers(subset(seu, subset = tmp20Percent.3 %in% c("GroupA","GroupB","GroupC","GroupD")),
                                only.pos = T)
write.table(CD14MonoCOVID1657.admVent.markers, file = paste0(seuName,"/CD14MonoCOVID1657.admVent.markers.txt"), sep = "\t", row.names = F)

Idents(seu) <- seu$tmp20Percent.3
CD14MonoCOVID1657.group.markers <- FindAllMarkers(subset(seu, subset = tmp20Percent.3 %in% c("GroupA","GroupB","GroupC","GroupD")),
                                                  only.pos = T)
write.table(CD14MonoCOVID1657.group.markers, file = paste0(seuName,"/CD14MonoCOVID1657.group.markers.txt"), sep = "\t", row.names = F)
#refer to enrichR.R for plotting Figures 2F-G

#generate vlnplot of UMAP_gs sorted by variance
Idents(default.full) <- default.full$cell.type.fine
var.mat <- c()
clusName <- sort(unique(default.full$cell.type.fine))
for (c in clusName){
  var.mat <- rbind(var.mat,
                   cbind.data.frame(c,var(default.full$UMAP_gs[which(default.full$cell.type.fine==c)])))
}
orderClusName <- var.mat[,1][order(var.mat[,2], decreasing = T)]
default.full$ordered.cell.type.fine <- factor(default.full$cell.type.fine, levels = orderClusName)
VlnPlot(default.full, features = "UMAP_gs", group.by = "ordered.cell.type.fine") + NoLegend() +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) #Figure S8A


#validation with wilk et al 2020
DimPlot(seu, reduction = "umap", group.by = "tmp20Percent.3", order = T, split.by = "sample",
        cols = pal_jco("default")(4), na.value = "grey90") #Figure S9A

seu.exprs <- as.data.frame(FetchData(seu, vars = c("sample","tmp20Percent.3")))
plot_data <- seu.exprs[which(!is.na(seu.exprs$tmp20Percent.3)),]
colnames(plot_data)[2] <- "Group"
library(dplyr)
plot_data <- plot_data %>%
  group_by(sample, Group) %>%
  summarise(Count=n()) %>% as.data.frame()
# Grouped
ggplot(plot_data, aes(fill=sample, y=Count, x=Group)) + geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("C1 A"="darkred","C1 B"="lightcoral","C2"="navy","C3"="seagreen","C4"="orange",
                                "C5"="darkgoldenrod","C6"="plum3","C7"="grey")) #Figure S9B
