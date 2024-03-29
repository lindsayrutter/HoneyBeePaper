library(EDASeq)
library(bigPint)

# Read in data and dataMetrics
data <- readRDS("../../N_V/DESeq2/data.Rds")
dataMetrics <- readRDS("../../N_V/DESeq2/dataMetrics.Rds")

# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
datas = as.data.frame(t(apply(as.matrix(data), 1, scale)))
datas$ID = as.character(rownames(datas))
datas = datas[,c(25,1:24)]
colnames(datas) = c("ID", colnames(data))
nID = which(is.nan(datas[,2]))
datas[nID,2:length(datas)] = 0

# Define color vector
myCols = scales::hue_pal()(5)
colList = c(myCols[1], "#C11B8D", myCols[2:4])
colList = colList[2:5]

# Plot clusters as parallel coordinate lines
# Set verbose=TRUE to save images and .rds files of gene IDs to directory "ClusterFiles"
ret <- plotClusters(data=datas, dataMetrics = dataMetrics, threshVar = "padj", clusterAllData = FALSE, colList = colList, lineAlpha = 1, lineSize = 0.4, yAxisLabel = "Standardized count", vxAxis = TRUE, outDir = "ClusterFiles/N_V", saveFile = FALSE, verbose = TRUE)
plot(ret[["N_V_4"]])
