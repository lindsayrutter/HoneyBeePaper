library("merTools")
library("readr")

d <- read_delim("GSE65659_AntiviralResponseReadCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

data <- as.data.frame(d)
rownames(data) <- data[,1]
data <- data[,2:7]
colnames(data) <- c("N.1", "N.2", "N.3", "V.1", "V.2", "V.3")

saveRDS(data, "data.Rds")
