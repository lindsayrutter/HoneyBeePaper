library(bigPint)
library(EDASeq)

# Read in data and dataMetrics
data <- as.data.frame(readRDS("../../Galbraith/DESeq2/data.Rds"))

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

# Read in geneList for Cluster 1 created in makeFigure3.R
geneList <- readRDS("ClusterFiles/Galbraith/N_V_4_1.rds")

# Plot clusters as parallel coordinate lines
# Set verbose=TRUE to save images and .rds files of gene IDs to directory "ClusterFiles"
ret <- plotSM(data=datas, geneList = geneList, saveFile = FALSE, pointColor = "#E9AA0D")
ret[["N_V"]] + labs(x = "Standardized count", y = "Standardized count", title = paste0("Cluster 1 (n = ", length(geneList), ")")) + theme_gray() + theme(plot.title = element_text(hjust = 0.5, size = 13), axis.title = element_text(size = 13)) 

