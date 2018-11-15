library(bigPint)
library(EDASeq)

# Read in data and dataMetrics
data <- as.data.frame(readRDS("../../NC_NR_VC_VR/data/data.Rds"))

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

# Rename columns to fit bigPint specifications
colnames(datas) <- c("ID", "N.C.1", "N.C.2", "N.C.3", "N.C.4", "N.C.5", "N.C.6", "N.R.1", "N.R.2", "N.R.3", "N.R.4", "N.R.5", "N.R.6", "V.C.1", "V.C.2", "V.C.3", "V.C.4", "V.C.5", "V.C.6", "V.R.1", "V.R.2", "V.R.3", "V.R.4", "V.R.5", "V.R.6")

# Define color vector
colList = c("#00A600FF", "#7570B3", "#f97976", "#0066FFFF")

# Read in geneList for tolerance genes
geneList <- readRDS("../../ResistanceTolerance/resistance.rds")

# Plot clusters as parallel coordinate lines
# Set verbose=TRUE to save images and .rds files of gene IDs to directory "ClusterFiles"
ret <- plotClusters(data=datas, geneList=geneList, clusterAllData = FALSE, yAxisLabel = "Standardized count", colList = colList, saveFile = FALSE, vxAxis = TRUE, lineAlpha = 1, lineSize = 0.3)
plot(ret[["N_V_4"]])
