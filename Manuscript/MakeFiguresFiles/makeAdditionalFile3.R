library(EDASeq)
library(bigPint)

# Read in data and dataMetrics
data <- readRDS("../../C_R/DESeq2/data.Rds")
dataMetrics <- readRDS("../../C_R/DESeq2/dataMetrics.Rds")

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
myCols = scales::hue_pal()(6)

# Plot clusters as parallel coordinate lines
ret <- plotClusters(data=datas, dataMetrics = dataMetrics, nC=6, threshVar = "padj", clusterAllData = FALSE, colList = colList, lineSize = 0.1, yAxisLabel = "Standardized count", vxAxis = TRUE, saveFile = FALSE)
plot(ret[["C_R_6"]])
