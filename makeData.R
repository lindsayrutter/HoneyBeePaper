library("merTools")

d <- read_delim("~/VirusHoneyBee/GSE65659_AntiviralResponseReadCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

data <- as.data.frame(d)
data <- data[,1:7]
colnames(data) <- c("ID", "C.1", "C.2", "C.3", "T.1", "T.2", "T.3")

saveRDS(data, "data/data.Rds")
