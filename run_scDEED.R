#downloaded "Simulated Data-Done" from scDEED github
#most scripts are adopted from Simulation.Rmd, downloaded from their github

loc <- "~/"

#```{r Functions}
find_neighbors = function(embeddings, neighbors){
  original_neighbors = matrix(ncol = neighbors, nrow = dim(embeddings)[1])
  rownames(original_neighbors) = rownames(embeddings)
  dist = as.matrix(dist(embeddings))
  for (i in 1:dim(embeddings)[1]){
    order = order(dist[i, ], decreasing = F)
    original_neighbors[i, ] = rownames(embeddings)[order[2:(neighbors+1)]]
  }
  return(original_neighbors)}

cellclass_ids = function(data, ident){
  Idents(data) = ident
  classes = unique(data@active.ident)
  results = list()
  for (i in 1:length(classes)){
    results[[classes[i]]] = which(data@active.ident == classes[i])
  }
  return(results)
}

find_cluster_neighbors = function(embeddings,cellclass_ids, neighbors){

  cluster_centers = matrix(nrow = length(cellclass_ids), ncol = dim(embeddings)[2])
  for (i in 1:length(cellclass_ids)){
    dat = as.matrix(embeddings[cellclass_ids[[i]], ])
    cluster_centers[i, ] = colMedians(dat)
  }
  rownames(cluster_centers) = names(cellclass_ids)
  original_neighbors = matrix(ncol = neighbors, nrow = length(cellclass_ids))

  cluster_centers = as.matrix(dist(cluster_centers))
  for (i in 1:nrow(original_neighbors)){
    order = order(cluster_centers[i, ], decreasing = F)
    original_neighbors[i, ] = colnames(cluster_centers)[order[2:(neighbors+1)]]
  }
  return(original_neighbors)}


KNN_overlap = function(data, neighbors, reduction){
  new_embeddings = Embeddings(data, reduction)
  pca_embeddings = Embeddings(data, 'pca')
  new_embeddings = find_neighbors(new_embeddings, neighbors)
  pca_embeddings = find_neighbors(pca_embeddings, neighbors)
  res = c()
  for (i in 1:nrow(new_embeddings)){
    res[i] = length(intersect(new_embeddings[i, ], pca_embeddings[i, ]))
  }
  return(list(mean = mean(res), list = res))
}

KNC_overlap = function(data, neighbors, reduction, cellclass_ids){
  new_embeddings = Embeddings(data, reduction)
  pca_embeddings = Embeddings(data, 'pca')
  new_embeddings = find_cluster_neighbors(new_embeddings, cellclass_ids, neighbors)
  pca_embeddings = find_cluster_neighbors(pca_embeddings, cellclass_ids, neighbors)
  res = c()
  for (i in 1:nrow(new_embeddings)){
    res[i] = length(intersect(new_embeddings[i, ], pca_embeddings[i, ]))
  }
  return(list(mean = mean(res), list = res))
}

#```{r Get data}
library(scDesign2)
data_mat <- readRDS(system.file("extdata", "mouse_sie_10x.rds", package = "scDesign2"))
nonspikes <- which(!grepl("ercc", rownames(data_mat), ignore.case = TRUE))
print(paste("number of spike-ins:", nrow(data_mat)-length(nonspikes)))
#> [1] "number of spike-ins: 9"
data_mat <- data_mat[nonspikes, ,drop = FALSE]
# explore basic structure of data -------------------------------------------------------
dim(data_mat)
celltypes = colnames(data_mat)


KN_results = data.frame('embeddingScore' = c(), 'TrustworthyCells' = c(),'DubiousCells' = c(), 'NeitherCells' = c(), 'KNN' = c(), 'KNC' = c(), 'iteration' = c())
all_counts = data.frame('embeddingScore' = c(), 'TrustworthyCells' = c(),'DubiousCells' = c(), 'NeitherCells' = c())
for (i in 1:20){
  print(paste0("Simulation:",i))
  simulated <- readRDS(paste0(loc,"Simulated Data-Done/simulated_data_",i,".Rds"))
  simulated$celltypes <- celltypes

  simulated <- NormalizeData(simulated, verbose = FALSE) #LogNormalize
  simulated <- FindVariableFeatures(simulated, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  simulated <- ScaleData(simulated)
  simulated <- RunPCA(simulated)
  simu_umap.1 = scDEED(simulated, num_pc = 12, use_method = "umap",visualization = TRUE, min.dist = 0.3,
                     n_neighbors = 30, optimize_neib = F, optimize_min = F)
  simu_tsne.1 = scDEED(simulated, num_pc = 12, use_method = "tsne",visualization = TRUE, perplexity = 30)

  results <- c()
  counts <- c()
  x_list <- c("scDEED-tSNE","scDEED-UMAP","cellstruct-tSNE","cellstruct-UMAP")
  y_list <- list("Trustworthy cells", c("Dubious cells","Neither cells"))
  names <- c("Trustworthy","Non-trustworthy")
  reduction <- c("umap","tsne") # 2, 1
  for (x in 1:4){
      y = 1 #subset only the trustworthy cells, then rerun scDEED
      tmp <- x_list[x]
      if (x == 1){
        tmp.obj <- subset(simu_tsne.1, subset = scDEED.classify.tSNE %in% y_list[[y]])
        tmp.counts <- noquote(cbind(ES=x_list[x],t(table(factor(simu_tsne.1$scDEED.classify.tSNE, level=unlist(y_list))))))
      }else if (x == 2){
        tmp.obj <- subset(simu_umap.1$object, subset = scDEED.classify.UMAP %in% y_list[[y]])
        tmp.counts <- noquote(cbind(ES=x_list[x],t(table(factor(simu_umap.1$object$scDEED.classify.UMAP, level=unlist(y_list))))))
      }else if (x == 3){
        tmp.obj <- subset(simu_tsne.1, subset = cellstruct.classify.tSNE %in% y_list[[y]])
        tmp.counts <- noquote(cbind(ES=x_list[x],t(table(factor(simu_tsne.1$cellstruct.classify.tSNE, level=unlist(y_list))))))
      }else{
        tmp.obj <- subset(simu_umap.1$object, subset = cellstruct.classify.UMAP %in% y_list[[y]])
        tmp.counts <- noquote(cbind(ES=x_list[x],t(table(factor(simu_umap.1$object$cellstruct.classify.UMAP, level=unlist(y_list))))))
      }


      tmp.obj <- FindVariableFeatures(tmp.obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      tmp.obj <- ScaleData(tmp.obj)
      tmp.obj <- RunPCA(tmp.obj)
      if (x %in% c(1,3)){
        tmp.obj <- RunTSNE(tmp.obj, dims = 1:12, seed.use = 100, perplexity=30) #default tSNE parameters
        tmp.obj = scDEED(tmp.obj, num_pc = 12, use_method = "tsne",visualization = TRUE, perplexity = 30)
        if (x == 1){
          tmp <- as.data.frame(cbind(tmp, t(table(factor(tmp.obj$scDEED.classify.tSNE, level=unlist(y_list))))))
        }else if (x == 3){
          tmp <- as.data.frame(cbind(tmp, t(table(factor(tmp.obj$cellstruct.classify.tSNE, level=unlist(y_list))))))
        }
      }else if (x %in% c(2,4)){
        tmp.obj <- RunUMAP(tmp.obj, dims = 1:12, seed.use = 100, min.dist = 0.3, n.neighbors = 30) #default UMAP parameters
        simu_umap.2 = scDEED(tmp.obj, num_pc = 12, use_method = "umap",visualization = TRUE, min.dist = 0.3,
                           n_neighbors = 30, optimize_neib = F, optimize_min = F)
        tmp.obj <- simu_umap.2$object
        if (x == 2){
          tmp <- as.data.frame(cbind(tmp, t(table(factor(tmp.obj$scDEED.classify.UMAP, level=unlist(y_list))))))
        }else if (x == 4){
          tmp <- as.data.frame(cbind(tmp, t(table(factor(tmp.obj$cellstruct.classify.UMAP, level=unlist(y_list))))))
        }
      }


      knn <- KNN_overlap(tmp.obj, 10, reduction[(x %% 2)+1])$mean
      cellclass.ids = cellclass_ids(tmp.obj, 'celltypes')
      knc <- KNC_overlap(tmp.obj, 4, reduction[(x %% 2)+1], cellclass.ids)$mean
      # tmp <- cbind.data.frame(tmp, names[y], ncol(tmp.obj), knn, knc)
      tmp <- cbind.data.frame(tmp, knn, knc)
      results <- rbind(results, tmp)
      counts <- rbind(counts, tmp.counts)
  }
  colnames(results) = c('embeddingScore', 'TrustworthyCells','DubiousCells', 'NeitherCells', 'KNN','KNC')
  counts <- as.data.frame(counts)

  results$iteration = i
  counts$iteration = i
  KN_results = rbind(KN_results,results)
  all_counts = rbind(all_counts, counts)
}
#write.table(all_counts, file=paste0(loc,"simulatedData/cellCounts.txt"), row.names = F, sep = "\t", quote = F)
#write.table(KN_results, file=paste0(loc,"simulatedData/trustworthySubset.KN_results.txt"), row.names = F, sep = "\t", quote = F)

plot_results <- KN_results#[which(KN_results$classification=="Trustworthy"),]
plot_results[c('method', 'DR')] <- str_split_fixed(plot_results$embeddingScore, '-', 2)
ggplot(plot_results, aes(x = KNN, y = KNC, shape = DR, col = method))+ geom_point() +
  scale_color_manual(values = c("scDEED"="navy","cellstruct"="darkred")) +
  ggtitle("Biological quantification of trustworthy cells") #Figure S11

summary(KN_results[which(KN_results$embeddingScore == "cellstruct-UMAP"), "KNC"])
summary(KN_results[which(KN_results$embeddingScore == "scDEED-UMAP"), "KNC"])
t.test(KN_results[which(KN_results$embeddingScore == "cellstruct-UMAP"), "KNC"],
       KN_results[which(KN_results$embeddingScore == "scDEED-UMAP"), "KNC"],
       paired = T)


## for Figure S10 (analysis with simulated dataset 1)
#simulated data 1 (from paste0(loc,"Simulated Data-Done/simulated_data_1.Rds"))
simulated <- simu_umap.1$object
clrs <- c("magenta", "darkgreen","lightgreen","khaki4","khaki1","steelblue","lightblue","navy",
          "cyan","tan1","grey","lightcoral","red","darkred","darkgoldenrod")
DimPlot(simulated, reduction = "umap", group.by = "celltypes", label = T, shuffle = T,cols = clrs) #Figure S10A
Idents(simulated) <- simulated$celltypes
exprs <- as.data.frame(FetchData(simulated, vars = c("UMAP_gs","rho_UMAP","celltypes")))

plot.list <- list()
for (p in 1:2){
  plot.list[[p]] <- FeaturePlot(simulated, reduction = "umap", features = colnames(exprs)[p], order = T) +
    ggtitle(paste0(colnames(exprs)[p],": ",round(mean(exprs[,p]),4))) +
    scale_color_distiller(palette = "RdYlBu", limits = c(0,1))
}
plot_grid(plotlist = plot.list, nrow = 1) #Figure S10B
p1 <- DimPlot(simulated, reduction = "umap", group.by = "cellstruct.classify.UMAP", cols = c("blue","seagreen"))
p2 <- DimPlot(simulated, reduction = "umap", group.by = "scDEED.classify.UMAP", cols = c("red","blue","seagreen"))
plot_grid(p1,p2, nrow = 1) # Figure S10F

#generating boxplot side-by-side
ggplot(melt(exprs), aes(x=celltypes,y=value,fill=variable)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust=1)) #Figure S10C

df <- melt(exprs[,c(1,3)]) %>% group_by(celltypes,variable) %>%
        dplyr::summarise(median_count=median(value)) %>%
        as.data.frame()

heatmapGS(simulated, "pca", "umap", colnames(simulated)[which(simulated$celltypes=="Enterocyte.Mature.Proximal")]) #Figure S10D
heatmapGS(simulated, "pca", "umap", colnames(simulated)[which(simulated$celltypes=="Stem")]) #Figure S10E

#saveRDS(simulated, paste0(loc,"simulatedData/sim1.rds"))

