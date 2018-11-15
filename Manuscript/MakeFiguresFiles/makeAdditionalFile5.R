library(EDASeq)
library(bigPint)

# Read in data and dataMetrics
data <- readRDS("../../Galbraith/DESeq2/data.Rds")
data <- as.data.frame(data)
data <- log(data+1)
data$ID <- rownames(data)
data <- data[,c(7, 1:6)]
dataMetrics <- readRDS("../../Galbraith/DESeq2/dataMetrics.Rds")

# Read in geneList for Cluster 1 created in makeFigure3.R
geneList <- readRDS("ClusterFiles/Galbraith/N_V_4_1.rds")

# Create interactive litre plot
# Inside the application, set the "Metric" field to "padj" and set the "Metric order" to increasing. AdditionalFile5.png was created by collecting nine of the first DEGs.
plotLitreApp(data = data, dataMetrics = dataMetrics, geneList = geneList, pointColor = "#E9AA0D")
