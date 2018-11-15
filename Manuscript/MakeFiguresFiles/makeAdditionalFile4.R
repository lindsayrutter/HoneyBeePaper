library(EDASeq)
library(bigPint)

# Read in data and dataMetrics
data <- readRDS("../../N_V/DESeq2/data.Rds")
data <- as.data.frame(data)
data <- log(data+1)
data$ID <- rownames(data)
data <- data[,c(25, 1:24)]
dataMetrics <- readRDS("../../N_V/DESeq2/dataMetrics.Rds")

# Read in geneList for Cluster 1 created in makeFigure4.R
geneList <- readRDS("ClusterFiles/N_V/N_V_4_1.rds")

# Create interactive litre plot
# Inside the application, set the "Metric" field to "padj" and set the "Metric order" to increasing. AdditionalFile4.png was created by collecting nine of the first DEGs.
plotLitreApp(data = data, dataMetrics = dataMetrics, geneList = geneList, pointColor = "#C11B8D")
