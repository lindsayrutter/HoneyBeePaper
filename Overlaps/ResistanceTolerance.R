library(readr)

# tolerance <- readRDS("intNCVC.Rds") #122
# resistance <- readRDS("VCminInt.Rds") #125
# virus <- readRDS("../N_V/DESeq2/RD_VIRUS_TOTAL.Rds") #43

#data <- as.data.frame(readRDS("../N_V/data/data.Rds"))
data <- as.data.frame(readRDS("../NC_NR_VC_VR/data/data.Rds"))

# Normalize for sequencing depth and other distributional differences between lanes
# data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
# data = as.data.frame(data)

# Another approach?
# data <- as.matrix(data)
# coldata = data.frame(row.names = colnames(data), treatment = unlist(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1])))
# dds = DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ treatment)
# dds <- DESeq(dds)

rsq <- function (x, y) cor(x, y) ^ 2

Variables <- read_csv("~/HoneyBeePaper/Variables.csv")
varIAPV <- as.data.frame(Variables[,7])[,1]

RSQ_virus_4 = list()
for (i in 1:4){
  temp <- readRDS(paste0("../N_V/DESeq2/ClusterStandard/Sig_4_", i, ".Rds"))
  temp2 <- data[which(rownames(data) %in% temp),]
  
  temp3 = c()
  for (j in 1:nrow(temp2)){
    temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
  }
  RSQ_virus_4[[i]] <- temp3
}

