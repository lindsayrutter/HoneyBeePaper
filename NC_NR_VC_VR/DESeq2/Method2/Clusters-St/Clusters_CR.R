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
library(data.table)
library(dplyr)

outDir = getwd()
data <- readRDS("../../../data/data.Rds")
data <- data[,c(1:6, 13:18, 7:12, 19:24)]
setDT(data, keep.rownames = TRUE)[]
colnames(data)[1] = "ID"
data <- as.data.frame(data)

metrics <- readRDS("../dataMetrics.Rds")[["C_R"]]
sigMets = metrics[which(metrics$padj<0.05),]
sigC <- sigMets[which(sigMets$log2FoldChange<0),]
sigR <- sigMets[which(sigMets$log2FoldChange>0),]


# File output information
plotName = "C_R"
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

yMin = min(dataqps[,1:24])
yMax = max(dataqps[,1:24])

#####################################################

colPurple = scales::seq_gradient_pal("purple4", "purple", "Lab")(seq(0,1,length.out=9))
colOrange = scales::seq_gradient_pal("orangered4", "darkorange2", "Lab")(seq(0,1,length.out=9))
colList = c(colOrange[5], colPurple[5])
Type = c("R Diet", "C Diet")

###########################

sigIDs = list(sigR$ID, sigC$ID)

plot_clustersSig = lapply(1:2, function(i){ 
  x = as.data.frame(dataqps[which(dataqps$ID %in% sigIDs[[i]]),])
  x$cluster = "color"
  x$cluster2 = factor(x$cluster)
  xNames = rownames(x)
  metricFDR = metrics[which(as.character(metrics$ID) %in% xNames),]
  sigID = metricFDR[metricFDR$padj<0.05,]$ID
  xSig = x[which(rownames(x) %in% sigID),]
  xSigNames = rownames(xSig)
  nGenes = nrow(xSig)
  
  xSig$ID = xSigNames
  pcpDat <- melt(xSig[,c(1:25)], id.vars="ID")
  colnames(pcpDat) <- c("ID", "Sample", "Count")
  pcpDat$Sample <- as.character(pcpDat$Sample)
  
  pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i], alpha=0.1) + ylab("Standardized Count") + ggtitle(paste("Significant Genes for ",  Type[i] ," (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=28), axis.text=element_text(size=28), axis.title=element_text(size=28))
  
  fileName = paste(outDir, "/", plotName, "_Sig_", i, ".jpg", sep="")
  jpeg(fileName)
  plot(pSig)
  invisible(dev.off())
  pSig
})

jpeg(file = paste(outDir, "/", plotName, "_Sig.jpg", sep=""), width=1300, height=700)
# We allow up to 4 plots in each column
p = do.call("grid.arrange", c(plot_clustersSig, ncol=2))
invisible(dev.off())

