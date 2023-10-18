
loc <- "~/"
seuName <- "humanLiver"
main <- "cells_70K"
seed <- seq(1,10)
for (s in seed){
  print(s)
  filename <- paste0(main,".",s)
  out.mat <- c()
  for (x in c("euclidean","cosine")){
    for (y in c(15,30,50)){
      for (z in c(0.05, 0.1, 0.3, 0.5)){
        data <- read.table(paste0(loc,seuName,"/downsampling/tuneUMAP/",filename,".",x,"_",y,"_",z,".txt"), sep = "\t", header = T)
        out.mat <- rbind(out.mat, cbind.data.frame(metric=x, n.neighbors=y, min.dist=z,
                                                  mean_GS=mean(data$gs), mean_LS=mean(data$ls), mean_GSLS=mean(c(data$gs,data$ls))))
      }
    }
  }
  write.table(out.mat, paste0(loc,seuName,"/downsampling/results/",main,"_param.",s,".txt"), sep = "\t", row.names = F)
}

