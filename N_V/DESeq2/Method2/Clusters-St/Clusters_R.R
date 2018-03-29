library(bigPint)
library(data.table)
library(EDASeq)
library(plyr)
library(gridExtra)

outDir = getwd()
data <- readRDS("../../../data/data.Rds")
data <- data[,c(1:6, 13:18, 7:12, 19:24)]
setDT(data, keep.rownames = TRUE)[]
colnames(data)[1] = "ID"
data <- as.data.frame(data)

metrics <- readRDS("../dataMetrics.Rds")[["C_R"]]
sigMets = metrics[which(metrics$padj<0.05),]
sigR <- sigMets[which(sigMets$log2FoldChange<0),]


# File output information
plotName = "R"
logSoy = data
logSoy[,-1] <- log(data[,-1]+1)







RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
data = mutate(data, mean = (NC.1+NC.2+NC.3+NC.4+NC.5+NC.6+NR.1+NR.2+NR.3+NR.4+NR.5+NR.6+VC.1+VC.2+VC.3+VC.4+VC.5+VC.6+VR.1+VR.2+VR.3+VR.4+VR.5+VR.6)/6, stdev = RowSD(cbind(NC.1,NC.2,NC.3,NC.4,NC.5,NC.6,NR.1,NR.2,NR.3,NR.4,NR.5,NR.6,VC.1,VC.2,VC.3,VC.4,VC.5,VC.6,VR.1,VR.2,VR.3,VR.4,VR.5,VR.6)))
rownames(data)=data_Rownames
data$ID <- data_Rownames

dataqps <- t(apply(as.matrix(data[,1:24]), 1, scale))
dataqps <- as.data.frame(dataqps)
colnames(dataqps) <- colnames(data[,1:24])
dataqps$ID <- rownames(dataqps)

# Comine the filtered and remaining data
fulls <- dataqps
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

# Indices of the 1311 NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(dataqps$NC.1))
# Set these filtered values that have all same values for samples to 0
dataqps[nID,1:24] <- 0

#####################################################

yMin = min(dataqps[,1:24])
yMax = max(dataqps[,1:24])

x = as.data.frame(dataqps[which(dataqps$ID %in% sigR$ID),])
x$cluster = "color"
x$cluster2 = factor(x$cluster)
xNames = rownames(x)
xSig = x
xSigNames = rownames(xSig)
nGenes = nrow(xSig)
xSig$ID = xSigNames


dendo = xSig[,1:24] # or dataqps? (If do fulls, then have NAs introduced by conversion)
rownames(dendo) = NULL
d = dist(as.matrix(dendo))
hc = hclust(d, method="ward.D")

fileName = paste(outDir, "/", plotName, "_dendogram.jpg", sep="")
jpeg(fileName)
plot(hc, main="data Dendogram", xlab=NA, sub=NA)
invisible(dev.off())





getPCP <- function(nC){
  
  set.seed(1)
  colList = scales::seq_gradient_pal("purple4", "purple", "Lab")(seq(0,1,length.out=nC))

  k = cutree(hc, k=nC)
  
  plot_clusters = lapply(1:nC, function(j){
    i = rev(order(table(k)))[j]
    x = as.data.frame(xSig[,1:24][which(k==i),])
    nGenes = nrow(x)
    x$cluster = "color"
    x$cluster2 = factor(x$cluster)
    xNames = rownames(x)
    x$ID = xNames
    xSigNames = rownames(x)
    saveRDS(xSigNames, file=paste0(outDir, "/Sig_", nC, "_", j, ".Rds"))
    
    
    
    # scatMatMetrics = list()
    # scatMatMetrics[["N_V"]] = metrics[which(metrics$ID %in% x$ID),]
    # scatMatMetrics[["N_V"]]$FDR = 10e-10
    # scatMatMetrics[["N_V"]]$ID = as.factor(as.character(scatMatMetrics[["N_V"]]$ID))
    # 
    # fileName = paste(getwd(), "/", plotName, "_Sig_SM_Orig_", nC, "_", j, ".jpg", sep="")
    # ret <- plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "padj", threshVal = 0.05, degPointColor = colList[j], fileName=fileName)
    # jpeg(fileName, height=700, width=700)
    # ret[[plotName]] + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Cluster ", j, " Significant Original Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14), axis.text=element_text(size=14), axis.title=element_text(size=18), strip.text = element_text(size = 14))
    # invisible(dev.off())
    
    
    
    pcpDat <- melt(x[,c(1:24,27)], id.vars="ID")
    colnames(pcpDat) <- c("ID", "Sample", "Count")
    boxDat$Sample <- as.character(boxDat$Sample)
    pcpDat$Sample <- as.character(pcpDat$Sample)
    
    boxDat$Sample <- factor(boxDat$Sample, levels = unique(boxDat$Sample),ordered = TRUE)
    
    p = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[j], alpha=0.5) + ylab("Standardized Count") + ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=32), axis.text=element_text(size=16), axis.title=element_text(size=32))
    
    fileName = paste(outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
    jpeg(fileName, width=650, height=500)
    plot(p)
    invisible(dev.off())
    p
  })
  jpeg(file = paste(outDir, "/", plotName, "_", nC, ".jpg", sep=""), width=1800, height=1800)
  # We allow up to 4 plots in each column
  if (nC %in% c(1,2,3,4)){
    p = do.call("grid.arrange", c(plot_clusters, ncol=2))
  }
  else{
    p = do.call("grid.arrange", c(plot_clusters, ncol=ceiling(nC/4)))
  }
  invisible(dev.off())
}

for (i in 2:6){
  getPCP(i)
}
