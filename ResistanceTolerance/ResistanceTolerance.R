library(readr)
source("functionRT.R")

# tolerance <- readRDS("intNCVC.Rds") #122
# resistance <- readRDS("VCminInt.Rds") #125
# virus <- readRDS("../N_V/DESeq2/RD_VIRUS_TOTAL.Rds") #43

#data <- as.data.frame(readRDS("../N_V/data/data.Rds"))
data <- as.data.frame(readRDS("../NC_NR_VC_VR/data/data.Rds"))


# Next lines to standardize
RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
data = data %>% mutate(mean= rowMeans(as.matrix(.[nms])), stdev=rowSds(as.matrix(.[nms])))
data$mean <- as.numeric(data$mean)
rownames(data)=data_Rownames
data$ID <- data_Rownames
datas <- t(apply(as.matrix(data[,1:nColumns]), 1, scale))
datas <- as.data.frame(datas)
colnames(datas) <- colnames(data[,1:nColumns])
datas$ID <- rownames(datas)
# Indices of the NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(datas$N.1))
# Set these filtered values that have all same values for samples to 0
datas[nID,1:nColumns] <- 0


# Next lines for logging
# data <- log(data+1)

# Another approach?
# data <- as.matrix(data)
# coldata = data.frame(row.names = colnames(data), treatment = unlist(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1])))
# dds = DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ treatment)
# dds <- DESeq(dds)

type = "Raw"
functionRT(data=data, type=type)

type = "BLN"
dataBTN <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
dataBTN = as.data.frame(dataBTN)
functionRT(data=dataBTN, type=type)

type = "Log"
dataLog <- log(data+1)
functionRT(data=dataLog, type=type)

