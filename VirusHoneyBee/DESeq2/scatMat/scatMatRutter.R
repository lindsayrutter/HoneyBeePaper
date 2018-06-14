library(bigPint)

data <- as.data.frame(readRDS("../../../N_V/DESeq2/data.Rds"))
logData <- log(data+1)
logData$ID <- rownames(logData)
logData <- logData[,c(25,1:24)]

logDataSub <- logData[,c(1,2:4, 14:16)]
ret <- plotScatterStatic(data = logDataSub, pointSize=0.1, option="point")

logDataSub <- logData[,c(1,5:7, 17:19)]
ret <- plotScatterStatic(data = logDataSub, pointSize=0.1, option="point")

logDataSub <- logData[,c(1,8:10, 20:22)]
ret <- plotScatterStatic(data = logDataSub, pointSize=0.1, option="point")

logDataSub <- logData[,c(1,11:13, 23:25)]
ret <- plotScatterStatic(data = logDataSub, pointSize=0.1, option="point")
