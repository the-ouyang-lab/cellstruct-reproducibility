library(ggplot2)
library(cowplot)
library(enrichR)
setEnrichrSite("Enrichr") # Human genes
dbs <- "Reactome_2022"

loc <- "~/wilk2020/"
tmp.markers <- read.table(paste0(loc,"CD14MonoCOVID1657.group.markers.txt"), sep="\t", header = T)
groupClus <- c("GroupA","GroupB","GroupC")

# tmp.markers <- read.table(paste0(loc,"CD14MonoCOVID1657.admVent.markers.txt"), sep="\t", header = T)
# groupClus <- unique(tmp.markers$cluster)
enriched.list <- list()
for (i in 1:length(groupClus)){
  print(groupClus[i])
  enriched.list[[i]] <- enrichr(tmp.markers[which(tmp.markers$p_val_adj<0.05&tmp.markers$avg_log2FC>0&tmp.markers$cluster==groupClus[i]),7], dbs)
}

plot.list <- list()
for (i in 1:length(groupClus)){
  # plot.list[[i]] <- plotEnrich(enriched.list[[i]][[1]][which(enriched.list[[i]][[1]]$Adjusted.P.value<0.05),], showTerms = 10, numChar = 40, y = "Count", orderBy = "P.value")+
  #   ggtitle(groupClus[i])
  enriched.list[[i]][[1]] <- cbind(enriched.list[[i]][[1]], truncTerm=stringr::str_trunc(enriched.list[[i]][[1]]$Term, 40),
                                   log10AdjP=-1*log10(enriched.list[[i]][[1]]$Adjusted.P.value),
                                   log10CombinedScore=log10(enriched.list[[i]][[1]]$Combined.Score))
  plot.list[[i]] <- ggplot(enriched.list[[i]][[1]][which(enriched.list[[i]][[1]]$Adjusted.P.value<0.05),][c(1:10),],
                           aes(x=reorder(truncTerm,log10AdjP),y=log10CombinedScore, fill=log10AdjP))+geom_bar(stat="identity") +
    scale_fill_gradient(low="blue",high="red") + coord_flip() +
     ggtitle(groupClus[i]) + ylab("log10(Combined score)") + xlab("") +theme_bw() #3x6 inches

}
plot_grid(plotlist = plot.list, ncol=3) #Figures 2F, 2G

