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

data <- as.data.frame(readRDS("../../../data/data.Rds"))
metricsAll <- readRDS("../dataMetrics.Rds")

data$ID <- rownames(data)
data <- data[,c(7,1:6)]

# Filter, normalize, and standardize the data so each gene has mean=0 and stdev=1
res <- filterStandardizeKL(data)
# Fitered data standardized
datas <- res[["fulls"]]

# Indices of the NAN rows
nID <- which(is.nan(datas$C.1))
# Set these filtered values that have all same values for samples to 0
datas[nID,1:6] <- 0

yMin = min(datas[,1:6])
yMax = max(datas[,1:6])

datas <- datas[,c(7,1:6)]
nGenes <- 213

metricsAll[["C_T"]]$ID = as.factor(as.character(metricsAll[["C_T"]]$ID))
metricsAll[["C_T"]]$FDR = 1
metricsAll[["C_T"]]$FDR[1:nGenes] = 0.01

ret <- plotDEG(data = datas, dataMetrics = metricsAll, option="scatterPoints", threshVar = "FDR", threshVal = 0.05) # degPointColor = colList[j]
ret[["C_T"]] + xlab("Standardized Count") + ylab("Standardized Count") + ggtitle(paste("Galbraith EdgeR-Btwn Top Significant Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=11), axis.text=element_text(size=11), axis.title=element_text(size=12), strip.text = element_text(size = 10))

filterStandardizeKL <- function(data){
  
  data_Rownames <- data$ID
  data = data[,-1]
  rownames(data) <- data_Rownames
  data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
  data = as.data.frame(data)
  # Add mean and standard deviation for each row/gene
  data = mutate(data, mean = (C.1+C.2+C.3+T.1+T.2+T.3)/6, stdev = RowSD(cbind(C.1,C.2,C.3,T.1,T.2,T.3)))
  rownames(data)=data_Rownames
  data$ID <- data_Rownames
  
  dataqps <- t(apply(as.matrix(data[,1:6]), 1, scale))
  dataqps <- as.data.frame(dataqps)
  colnames(dataqps) <- colnames(data[,1:6])
  dataqps$ID <- rownames(dataqps)
  
  # Combine the filtered and remaining data
  fulls <- dataqps
  boxDat <- melt(fulls, id.vars="ID")
  colnames(boxDat) <- c("ID", "Sample", "Count")
  
  # Indices of the 775 NAN rows. They had stdev=0 in the filt data
  nID <- which(is.nan(dataqps$K.1))
  # Set these filtered values that have all same values for samples to 0
  dataqps[nID,1:6] <- 0
  
  # Return several parameters
  list(datas = dataqps, fulls = fulls)
}

# This function calculates the standard deviation of each row in a data frame
RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
