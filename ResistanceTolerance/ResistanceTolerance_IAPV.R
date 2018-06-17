library(readr)
library(ggplot2)
library(EDASeq)
library(dplyr)
library(genefilter)
library(nlme)
source("functionRT_IAPV.R")

data <- as.data.frame(readRDS("../NC_NR_VC_VR/data/data.Rds"))

RawPVal = data.frame()
type = "Raw"
RawPVal <- functionRT(data=data, type=type, RawPVal=RawPVal)

type = "BLN"
dataBTN <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
dataBTN = as.data.frame(dataBTN)
RawPVal <- functionRT(data=dataBTN, type=type, RawPVal=RawPVal)

type = "Log"
dataLog <- log(data+1)
RawPVal <- functionRT(data=dataLog, type=type, RawPVal=RawPVal)

type = "Standardize"
# Next lines to standardize
RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
colData <- colnames(data)
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
rowNamesDat = rownames(data)
# Add mean and standard deviation for each row/gene
data = data %>% mutate(mean= rowMeans(as.matrix(data)), stdev=rowSds(as.matrix(data)))
data$mean <- as.numeric(data$mean)
datas <- t(apply(as.matrix(data[,1:(ncol(data)-2)]), 1, scale))
datas <- as.data.frame(datas)
colnames(datas) <- colData
# Indices of the NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(datas$NC.1))
# Set these filtered values that have all same values for samples to 0
datas[nID,] <- 0
datas <- datas[,1:24]
rownames(datas) <- rowNamesDat
RawPVal_IAPV <- functionRT(data=datas, type=type, RawPVal=RawPVal)
