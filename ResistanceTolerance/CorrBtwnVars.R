library(GGally)

Variables <- read_csv("~/HoneyBeePaper/Variables.csv")
ColNames <- colnames(Variables)
attributes(Variables) <- NULL

output <- matrix(unlist(Variables), ncol = 9)
output <- as.data.frame(output)
colnames(output) <- ColNames

output[,1] <- as.integer(output[,1])
output[,2] <- as.character(output[,2])
output[,3] <- as.character(output[,3])
output[,4] <- as.integer(output[,4])
output[,5] <- as.numeric(output[,5])
output[,6] <- as.numeric(output[,6])
output[,7] <- as.numeric(output[,7])
output[,8] <- as.numeric(output[,8])
output[,9] <- as.numeric(output[,9])

ggpairs(output[,c(5:7)])
