library("merTools")
library("readr")

d <- read_delim("GSE65659_AntiviralResponseReadCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

data <- as.data.frame(d)
rownames(data) <- data[,1]
data <- data[,2:7]
colnames(data) <- c("C.1", "C.2", "C.3", "T.1", "T.2", "T.3")

saveRDS(data, "data.Rds")
