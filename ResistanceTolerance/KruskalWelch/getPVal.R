library(readr)
library(ggplot2)
library(EDASeq)
library(dplyr)
library(genefilter)
library(nlme)
library(onewaytests)
source("functionRT_IAPV.R")

data <- as.data.frame(readRDS("../../NC_NR_VC_VR/data/data.Rds"))

PVal <- data.frame()
PVal <- functionRT(data=data, type="IAPV", colNum=7, PVal=PVal)
PVal <- functionRT(data=data, type="SBV", colNum=6, PVal=PVal)
PVal <- functionRT(data=data, type="Mortality", colNum=5, PVal=PVal)

saveRDS(PVal, "PVal.Rds")
