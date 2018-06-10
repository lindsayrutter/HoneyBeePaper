library(readr)
library(ggplot2)
library(EDASeq)
library(dplyr)
library(genefilter)
library(nlme)
source("functionRT_Mort.R")

data <- as.data.frame(readRDS("../NC_NR_VC_VR/data/data.Rds"))

type = "Raw"
functionRT(data=data, type=type)

type = "BLN"
dataBTN <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
dataBTN = as.data.frame(dataBTN)
functionRT(data=dataBTN, type=type)

type = "Log"
dataLog <- log(data+1)
functionRT(data=dataLog, type=type)

type = "Standardize"
# Next lines to standardize
RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
colData <- colnames(data)
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
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
functionRT(data=datas, type=type)
