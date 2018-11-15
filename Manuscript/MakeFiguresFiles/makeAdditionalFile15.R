library(bigPint)
library(EDASeq)

# Read in data and dataMetrics
data <- as.data.frame(readRDS("../../Galbraith/DESeq2/data.Rds"))
dataMetrics <- readRDS("../../Galbraith/DESeq2/dataMetrics.Rds")

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

# Define color vector
colList = scales::hue_pal()(5)
colList[2] = "#E9AA0D"
colList[3] = "#EA502F"
colList = colList[2:5]

# Plot clusters as parallel coordinate lines
# Set verbose=TRUE to save images and .rds files of gene IDs to directory "ClusterFiles"
ret <- plotClusters(data=datas, dataMetrics = dataMetrics, threshVar = "padj", clusterAllData = FALSE, colList = colList, yAxisLabel = "Standardized count", outDir = "ClusterFiles/Galbraith", saveFile = FALSE, verbose = TRUE)
plot(ret[["N_V_4"]])
