library(bigPint)

data <- as.data.frame(readRDS("../data.Rds"))
logData <- log(data+1)
logData$ID <- rownames(logData)
logData <- logData[,c(7,1:6)]

# scatMatGalbraith.png
ret <- plotScatterStatic(data = logData, pointSize=0.1, option="point")

ret[[1]] + xlab("Logged Count") + ylab("Logged Count")

plotScatterStatic(data = logData, pointSize=0.5, option="orthogonal")

plotScatterStatic(data = logData, pointSize=0.5, piLevel =0.9999, option="prediction")
