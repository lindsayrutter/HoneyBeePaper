# Run in /Users/suzu/HoneyBeePaper/VirusHoneyBee/DESeq2

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

data <- as.data.frame(readRDS("../../Galbraith/DESeq2/data.Rds"))
dataMetrics <- readRDS("../../Galbraith/DESeq2/dataMetrics.Rds")

logData = data
logData[,-1] <- log(data[,-1]+1)

# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
datas = as.data.frame(t(apply(as.matrix(data), 1, scale)))
datas$ID = as.character(rownames(datas))
datas = datas[,c(7,1:6)]
colnames(datas) = c("ID", colnames(data))
nID = which(is.nan(datas[,2]))
datas[nID,2:length(datas)] = 0

plotClusters(data=datas, dataMetrics = dataMetrics, threshVar = "padj")

saveRDS(sigDatas$ID, file="Sig.Rds")


