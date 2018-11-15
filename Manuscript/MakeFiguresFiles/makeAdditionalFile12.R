library(bigPint)
library(EDASeq)
library(dplyr)

# Read in data and dataMetrics
data <- as.data.frame(readRDS("../../N_V/DESeq2/data.Rds"))

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
datas = datas %>% select(ID, N.4, N.5, N.6, V.4, V.5, V.6)

# Read in geneList for Cluster 1 created in makeFigure4.R
geneList <- readRDS("ClusterFiles/N_V/N_V_4_1.rds")

# Plot cluster as scatterplot matrix
ret <- plotSM(data=datas, geneList = geneList, saveFile = FALSE, pointColor = "#C11B8D")
ret[["N_V"]] + labs(x = "Standardized count", y = "Standardized count", title = paste0("Cluster 1 (n = ", length(geneList), "): Replicates 4, 5, 6")) + theme_gray() + theme(plot.title = element_text(hjust = 0.5, size = 13), axis.title = element_text(size = 13)) 
