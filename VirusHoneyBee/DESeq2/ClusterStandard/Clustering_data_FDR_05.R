library(rtracklayer)
library(Rsamtools)
library(grid)
library(GenomicAlignments)
library(ggplot2)
library(GGally)
library(edgeR)
library(stringr)
library(EDASeq)
library(dplyr)
library(matrixStats)
library(gridExtra)
library(reshape2)
library(scales)
library(bigPint)
library(matrixStats)

getPCP <- function(nC){
  
  set.seed(1)
  colList = scales::hue_pal()(nC+1)
  colList[2] = "#E9AA0D"
  colList[3] = "#EA502F"
  k = cutree(hc, k=nC)
  
  yMin = min(sigDatas[,1:nColumns])
  yMax = max(sigDatas[,1:nColumns])
  
  ###########################
  
  sbsDF <- data.frame()
  for (i in 1:nC){
    x = as.data.frame(sigDatas[which(k==i),])
    xNames = rownames(x)
    xFDR = metrics[which(metrics$ID %in% xNames),]$padj
    sbsDF = rbind(sbsDF, data.frame(Cluster = paste("Cluster", i), FDR = xFDR))
  }
  
  plot_clusters = lapply(1:nC, function(i){
    j = rev(order(table(k)))[i]
    x = as.data.frame(sigDatas[which(k==j),])
    nGenes = nrow(x)
    x$cluster = "color"
    x$cluster2 = factor(x$cluster)
    xNames = rownames(x)
    xFDR = metrics[which(metrics$ID %in% xNames),]$padj
    scatMatMetrics = list()
    scatMatMetrics[[currPair]] = metrics[which(metrics$ID %in% xNames),]
    scatMatMetrics[[currPair]]$FDR = 10e-10
    scatMatMetrics[[currPair]]$ID = as.factor(as.character(scatMatMetrics[[currPair]]$ID))
    
    plotDatas = datas[, c(ncol(datas), 1:ncol(datas)-1)]
    plotData = data[,c(9,1:6)]
    plotData[,c(2:7)] <- log(plotData[,c(2:7)] +1)
    
    fileName = paste(getwd(), "/", outDir, "/", currPair, "_SM_", nC, "_", i, ".jpg", sep="")
    ret <- plotDEG(data = plotDatas, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "padj", threshVal = 0.05, degPointColor = colList[i+1], fileName=fileName)
    
    x$ID = xNames
    saveRDS(xNames, file=paste0("Sig_", nC,"_", i,".Rds"))
    
    pcpDat <- melt(x[,c(1:(nColumns+1))], id.vars="ID")
    colnames(pcpDat) <- c("ID", "Sample", "Count")
    boxDat$Sample <- as.character(boxDat$Sample)
    pcpDat$Sample <- as.character(pcpDat$Sample)
    
    p = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i+1], alpha=0.15) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Standardized Count") + theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1))
    
    fileName = paste(getwd(), "/", outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
    jpeg(fileName)
    plot(p)
    invisible(dev.off())
    p
  })
  
  ggBP = ggplot(sbsDF, aes(x=Cluster, y=FDR)) +
    stat_boxplot(geom ='errorbar') + 
    geom_boxplot(outlier.shape=NA, aes(fill=Cluster), alpha = 0.3) +
    geom_point(aes(fill=Cluster), shape=21, position=position_jitter(width=0.3), alpha=0.5) +
    scale_fill_manual(values=colList[c(2:length(colList), 1)])
  jpeg(file = paste(getwd(), "/", outDir, "/", currPair, "_boxplot_", nC, ".jpg", sep=""), width=1000, height=700)
  ggBP
  invisible(dev.off())

   jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_", nC, ".jpg", sep=""), width=1000, height=700)
  p = do.call("grid.arrange", c(plot_clusters, ncol=ceiling(nC/2)))
   invisible(dev.off())
}
  
i=1
metricsAll <- readRDS("../dataMetrics.Rds")
pairs <- names(metricsAll)
currPair <- pairs[i]

pair1 <- strsplit(currPair, "_")[[1]][1]
pair2 <- strsplit(currPair, "_")[[1]][2]

metrics <- metricsAll[[currPair]]
data <- as.data.frame(readRDS("../data.Rds"))
data <- data[,which(sapply(colnames(data), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(pair1, pair2))]
data<- cbind(ID = rownames(data), data)
data$ID <- as.character(data$ID)
nms <- colnames(data[-1])
nColumns <- length(data)-1

logData = data
logData[,-1] <- log(data[,-1]+1)

RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
data = data %>% mutate(mean= rowMeans(as.matrix(.[nms])), stdev=rowSds(as.matrix(.[nms])))
data$mean <- as.numeric(data$mean)
rownames(data)=data_Rownames
data$ID <- data_Rownames

datas <- t(apply(as.matrix(data[,1:nColumns]), 1, scale))
datas <- as.data.frame(datas)
colnames(datas) <- colnames(data[,1:nColumns])
datas$ID <- rownames(datas)
# Indices of the NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(datas$N.1))
# Set these filtered values that have all same values for samples to 0
datas[nID,1:nColumns] <- 0

# Combine the filtered and remaining data
boxDat <- melt(datas, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

sigDatas = datas[which(metricsAll[["C_T"]]$padj<0.05),]

saveRDS(sigDatas$ID, file="Sig.Rds")

dendo = sigDatas
rownames(dendo) = NULL
d = dist(as.matrix(dendo))
hc = hclust(d, method="ward.D")

plotName = currPair
outDir = "Clustering_data_FDR_05"

fileName = paste(getwd(), "/", outDir, "/", currPair, "_dendogram.jpg", sep="")
jpeg(fileName)
plot(hc, main="data Dendogram", xlab=NA, sub=NA)
invisible(dev.off())

getPCP(4)

