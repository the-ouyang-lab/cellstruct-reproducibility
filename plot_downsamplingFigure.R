
loc <- "~/"
seuName <- "humanLiver"

plot_cells <- c(5,8,seq(10,70,10))
col <- 4
out.mat <- c()
for (i in plot_cells){
  for (j in 1:10){
    data <- read.table(paste0(loc,seuName,"/downsampling/results/cells_",i,"K_param.",j,".txt"), sep = "\t", header = T)
    out.mat <- rbind(out.mat,
                     cbind.data.frame(parameters=paste(data[,1], data[,2], data[,3],sep=":"), variable=i, value=data[,col]))
  }
}

#reading full dataset ########################################################################
data <- read.table(paste0(loc,seuName,"/",seuName,".unweighted_param.txt"), sep = "\t", header = T)
out.mat <- rbind(out.mat,
                 cbind.data.frame(parameters=paste(data[,1], data[,2], data[,3],sep=":"), variable=79.5, value=data[,col]))
####################################################################################################

final <- data_summary(out.mat, varname="value", groupnames=c("parameters", "variable"))
#final$variable <- as.numeric(as.vector(final$variable))
pdf(paste0(loc,seuName,"_downsamplingCells.pdf"), w=9, h=5)
ggplot(final, aes(x=variable, y=value, colour=parameters, group=parameters)) + geom_line() +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(0.05)) +
  xlab("Number of cells (in K)") + ylab("Mean of GS")
dev.off()

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

