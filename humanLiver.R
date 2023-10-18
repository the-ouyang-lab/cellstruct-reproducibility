library(cellstruct)
library(dplyr)
library(cowplot)

loc <- "~/"
seuName <- paste0(loc,"humanLiver")
seu <- readRDS(paste0(seuName,".rds"))
clusterVar <- "celltype.l2"

num_PC <- 50
seu <- RunUMAP(seu, reduction = "refDR", dims = 1:num_PC)
seu <- RunTSNE(seu, reduction = "refDR", dims = 1:num_PC)

library(reticulate)
library(SeuratWrappers)
os <- import("os")
py_run_string("r.os.environ['OMP_NUM_THREADS'] = '20'")
sc <- import("scanpy", convert = FALSE)
ad <- import("anndata", convert = FALSE)

adata <- sc$AnnData(
  X   = np_array(t(GetAssayData(seu)[VariableFeatures(seu),]), dtype="float32"),
  obs = seu@meta.data[, clusterVar],
  var = data.frame(geneName = VariableFeatures(seu)))
adata$obsm$update(X_pca = Embeddings(seu, "refDR"))
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

ref.proj <- "refDR"
target.proj <- c("umap","tsne","fdl")
#run cellstruct on default hyperparameters
# (need to change the metadata stored for UMAP metric scores to "defaultUMAP...", avoiding replacement by the scores for tuned UMAP later)
seu <- run_cellstruct(seu, seuName, ref.proj, target.proj, clusterVar, nCores = 50)

tuned.full <- tuneUMAP(seu, seuName, "refDR", paste0(seuName,"/humanLiver"),
                    clusterVar, nCores = 50) #tuneUMAP returns seurat object, after running cellstruct on UMAP with max meanGS

cols = c("bisque4","lightskyblue","dodgerblue","orange","navy","mediumblue","hotpink","peru","grey","palegreen",
         "limegreen","olivedrab","maroon","magenta","orchid","seagreen","royalblue","yellow","indianred","violet",
         "mistyrose","khaki1","moccasin") #color scheme used in cell type annotation

## default UMAP embedding and score are stored in tuned.full (defaultUMAP)

#four cells selected for Figure S3
final.coi <- c("aizarani_v1_CD34pos315_5_17","macparland_v1_P3TLH_CGTGAGCTCCTAGAAC_1", #Hepatocyte
               "aizarani_v1_CD34pos_10_129","macparland_v1_P4TLH_CTTCTCTCAGGCTCAC_1")
final.list <- list()
target.proj <- c("umap","defaultUMAP","tsne","fdl")
for (p in 1:length(target.proj)){
  plot.list <- list()
  markedCells.loc <- as.data.frame(tuned.full@reductions[[target.proj[p]]]@cell.embeddings[match(final.coi,rownames(embed1)),])
  colnames(markedCells.loc) <- c("dim1","dim2")
  for (c in 1:length(final.coi)){
    plot.list[[c]] <- dimredGS(tuned.full, "refDR", target.proj[p], final.coi[c]) +
      geom_point(data = markedCells.loc, aes(x=dim1, y=dim2), shape=4, size=1, stroke=1, colour="black")
  }
  final.list <- c(final.list, plot.list)
}
plot_grid(plotlist = final.list, nrow = 4) #Figure S3B
heatmapGS(tuned.full, "refDR", c("umap","defaultUMAP","tsne","fdl"), coi = final.coi,
          target = final.coi, dist_percentile = 1) #Figure S3A


df <- data.frame(FetchData(tuned.full, vars = c("UMAP_gs","TSNE_gs","FDL_gs","UMAP_ls","TSNE_ls","FDL_ls",
                                                "UMAP_gc","TSNE_gc","FDL_gc",clusterVar)))
Idents(tuned.full) <- tuned.full$celltype.l2
var.mat <- c()
clusName <- sort(unique(df$celltype.l2))
for (c in clusName){
  var.mat <- rbind(var.mat, cbind.data.frame(c,var(df$UMAP_gs[which(df$celltype.l2==c)])))
}
orderClusName <- var.mat[,1][order(var.mat[,2], decreasing = T)]
tuned.full$ordered.celltype.l2 <- factor(tuned.full$celltype.l2, levels = orderClusName)
VlnPlot(tuned.full, features = "UMAP_gs", group.by = "ordered.celltype.l2") + NoLegend() +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) #Figure S4


#all.coi <- rownames(df)[which(df$celltype.l2=="Hepatocyte")] #randomly select 284 cells
all.coi <- rownames(df)[which(df$celltype.l2=="Cholangiocyte")] #randomly select 363 cells
set.seed(42)
coi <- all.coi[sample(c(1:length(all.coi)), size = 363)]
heatmapGS(tuned.full, "refDR", c("umap","tsne","fdl"), coi = coi) #14x5 inches (Figure 1F, S6A)
#for manual annotation of different groups of cells
# Heatmap(as.matrix(prop.embed1[cellOrder[246:284],]), name = lgd_name, col = col_fun, column_title = toupper(ref.proj),
#         bottom_annotation = HeatmapAnnotation(Target_cluster=randomCell.celltype), column_order = randOrder, cluster_rows = F) +
#   rowAnnotation(rn = anno_text(cellOrder[246:284]))
outlier.cells <- read.table(paste0(loc,"humanLiver/noWeight/cholangiocytesTop78.txt"),sep = "\t", header = T)
coi <- unique(c(coi,outlier.cells$outlierCells)) #428 cholangiocytes plot in Figure S6A
cholan.coi <- c("aizarani_v1_P310_15_163","aizarani_v1_Epcampos315_4_143","aizarani_v1_P308_31_92") #the three cells marked  out in Figure S6C
dimredGS(tuned.full, "refDR", "umap", cholan.coi[1]) #4x13.5 inches for 3 cells

#metadata "HepaGroup", "CholanGroup"
#4x5 inches (14x3.5 inches for 3 panels) - Figure S5, S6B (changing the group.by variable)
DimPlot(tuned.full, reduction = "umap", group.by = "HepaGroup", order = T, cols = c(pal_jco("default")(4)[c(1,4)],"orange"), na.value = "grey90")

heatmapGC(tuned.full, "refDR", c("umap","tsne","fdl")) #Figure 1D

############ for LS analysis using boundary cells (Figure S1) ######################################################
ref <- tuned.full@reductions[["refDR"]]@cell.embeddings
add.data <- merge(ref, df, by = 'row.names', all = T)
CToi <- c("CD4 T", "CD8 T")
Tcell.idx <- which(add.data$celltype.l2 %in% CToi)
embed1 <- add.data[,2:51]
rownames(embed1) <- add.data$Row.names

num_neighbor = 30
nearest.embed1 <- nn2(embed1, k = num_neighbor+1)
Tcell.nn <- nearest.embed1$nn.idx[Tcell.idx,2:ncol(nearest.embed1$nn.idx)]
Tcell.nn.celltype <- t(apply(Tcell.nn, 1, function (x) add.data$celltype.l2[x]))
pca.freq.table <- matrix(, nrow = length(Tcell.idx), ncol = 2, dimnames = list(add.data$Row.names[Tcell.idx],CToi))
pca.freq.table[,1] <- rowSums(Tcell.nn.celltype == CToi[1])
pca.freq.table[,2] <- rowSums(Tcell.nn.celltype == CToi[2])
boundaryCell <- rownames(pca.freq.table)[which(pca.freq.table[,1]>=11&pca.freq.table[,1]<=19&
                                                 pca.freq.table[,2]>=11&pca.freq.table[,2]<=19&(pca.freq.table[,1]+pca.freq.table[,2])==num_neighbor)]
nonBoundaryCell <- rownames(pca.freq.table)[which(!rownames(pca.freq.table) %in% boundaryCell)]

embed2.list <- list()
nearest.embed2.list <- list()
target.proj <- c("umap","tsne","fdl")
for (p in 1:length(target.proj)){
  embed2.list[[p]] <- tuned.full@reductions[[target.proj[p]]]@cell.embeddings[,1:2]
  nearest.embed2.list[[p]] <- nn2(embed2.list[[p]], k = num_neighbor+1)
}
coi_idx <- match(boundaryCell,rownames(embed2.list[[1]]))
embed2.neighbor_idx.list <- vector("list",length(embed2.list))
for (n in 1:length(coi_idx)){
  for (p in 1:length(embed2.list)){
    embed2.neighbor_idx.list[[p]][[n]] <- nearest.embed2.list[[p]]$nn.idx[coi_idx[n],2:ncol(nearest.embed2.list[[p]]$nn.idx)]
  }
}

embed2.adata <- add.data[match(rownames(embed2.list[[1]]),add.data$Row.names),]
embed2.freq.list <- list()
for (p in 1:length(embed2.list)){
  tmp.nn.celltype <- t(apply(matrix(unlist(embed2.neighbor_idx.list[[p]]), ncol = num_neighbor, byrow = TRUE), 1, function (x) embed2.adata$celltype.l2[x]))
  freq.table <- matrix(, nrow = length(boundaryCell), ncol = 3, dimnames = list(boundaryCell,c(CToi,"proportion")))
  freq.table[,1] <- rowSums(tmp.nn.celltype == CToi[1])
  freq.table[,2] <- rowSums(tmp.nn.celltype == CToi[2])
  freq.table[,3] <- rowSums(tmp.nn.celltype == embed2.adata$celltype.l2[match(boundaryCell,embed2.adata$Row.names)])/num_neighbor
  embed2.freq.list[[p]] <- freq.table[,3]
}
plot.list <- list()
for(p in 1:length(target.proj)){
  proj <- target.proj[p]
  tuned.full <- AddMetaData(object = tuned.full, metadata = rep(NA,ncol(tuned.full)),
                     col.name = paste0(toupper(proj),"_nnProp"))
  tuned.full@meta.data[paste0(toupper(proj),"_nnProp")][match(boundaryCell,colnames(tuned.full)),] <- embed2.freq.list[[p]]
  plot.list[[p]] <- FeaturePlot(tuned.full, reduction = proj, features = paste0(toupper(proj),"_nnProp"), order = T)+
                                scale_color_distiller(palette = "RdYlBu", limits = c(0,1), na.value = "grey90")
}
plot_grid(plotlist = plot.list, nrow = 1) #Figure S1E

tuned.full$PCA_nnProp <- rep(NA,ncol(tuned.full))
tuned.full$PCA_nnProp[match(add.data$Row.names[Tcell.idx],colnames(tuned.full))] <- rowSums(Tcell.nn.celltype == add.data$celltype.l2[Tcell.idx])/num_neighbor
tuned.full$PCA_nnProp[match(nonBoundaryCell,colnames(tuned.full))] <- NA
tmp.exprs <- na.omit(data.frame(FetchData(tuned.full, vars = paste0(c("PCA",toupper(target.proj)),"_nnProp"))))
tmp.exprs$PCA_nnProp <- as.character(round(tmp.exprs$PCA_nnProp,4))
ggplot(reshape::melt(tmp.exprs), aes(x=PCA_nnProp, y=value, fill=variable)) + geom_boxplot(width=0.4) #Figure S1F


## finding the rank of target embedding neighbors relative to ref neighbor
pca_nn = 150
pca <- embed1[match(rownames(embed2.list[[1]]),rownames(embed1)),]
nearest150.pca <- nn2(pca, k = pca_nn+1)
rank.list <- list()
for (p in 1:length(target.proj)){
  rank.list[[p]] <- matrix(, nrow = nrow(pca), ncol = num_neighbor)
  for (n in 1:nrow(pca)){
    rank.list[[p]][n,] <- match(nearest.embed2.list[[p]]$nn.idx[n,2:(num_neighbor+1)], nearest150.pca$nn.idx[n,2:(pca_nn+1)])
  }
}
data <- cbind.data.frame(UMAP=rowSums(!is.na(rank.list[[1]])),TSNE=rowSums(!is.na(rank.list[[2]])),FDL=rowSums(!is.na(rank.list[[3]])),
              celltype=embed2.adata$celltype.l2)
ggplot(reshape::melt(data[1:3]), aes(x=variable,y=value)) + geom_violin() + geom_boxplot(width=0.1) + xlab("") +
  ylab("Count of refDR neighbors") #i.e. within top 150 refDR neighbors (4x6 inches) - Figure S1C

tmp.mat1 <- reshape::melt(data) %>% group_by(celltype,variable) %>%
  dplyr::summarise(median_count=median(value)) %>%
  as.data.frame()
tmp.mat2 <- reshape::melt(df[,c(4:6,10)]) %>% group_by(celltype.l2,variable) %>%
  dplyr::summarise(median_LS=median(value)) %>%
  as.data.frame()
final.mat <- cbind.data.frame(tmp.mat1,tmp.mat2)[,c(1:3,6)]
ggplot(final.mat, aes(x=median_count,y=median_LS, color=variable))+
     geom_point(alpha=0.5) + stat_cor(method = "pearson", label.x = 0, label.y = c(1,0.95,0.9)) #Figure S1D
#################################################################################

###################### for stability analysis of GS ####################################
ref <- tuned.full@reductions[["refDR"]]@cell.embeddings
proj <- c("umap","tsne","fdl")
for (p in proj){
  print(paste0(p,date()))
  data <- tuned.full@reductions[[p]]@cell.embeddings
  for (k in c(1:10)){
    print(k)
    tuned.full <- AddMetaData(object = tuned.full,
                              metadata = calcGS(ref, data, num_waypoint = k*1000, nCores = 40),
                              col.name = paste0(toupper(p),"_gs.",k,"K"))
  }
}
gs.exprs <- as.data.frame(FetchData(tuned.full, vars = c(paste0("UMAP_gs.",seq(1,10),"K"),
                              paste0("TSNE_gs.",seq(1,10),"K"),paste0("FDL_gs.",seq(1,10),"K"))))
tmp_data <- melt(gs.exprs[,21:30]) #changing the columns to plot UMAP, tSNE, and FDL respectively
ggplot(tmp_data, aes(x=variable,y=value)) + geom_violin() + geom_boxplot(width=0.1) + ylim(0,1) #Figure S7A
##########################################################################################

######### downsampling analysis ##################################################
#refer to downsampling_humanLiver.R for downsampling dataset and tuning them
#tuningUMAP for complete dataset
metric <- c("euclidean","cosine")
umap_neighbors <- c(15,30,50)
min_dist <- c(0.05, 0.1, 0.3, 0.5)
for (x in metric){
  for (y in umap_neighbors){
    for (z in min_dist){
      singleTuneUMAP(tuned.full, x, y, z, "refDR", paste0(loc,"/downsampling/tuneUMAP/full"), clusterVar, nCores = 50)
    }
  }
}
#refer to compiling_downsampling.R for compiling the mean GS across all cells in a UMAP embedding
#generated under a combination of hyperparameters
#refer to plot_downsamplingFigure.R for plotting Figure S7B
##############################################################################################################
saveRDS(tuned.full, paste0(seuName,".rds"))

