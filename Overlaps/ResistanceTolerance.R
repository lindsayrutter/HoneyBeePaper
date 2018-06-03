library(readr)

tolerance <- readRDS("intNCVC.Rds") #122
resistance <- readRDS("VCminInt.Rds") #125
virus <- readRDS("../N_V/DESeq2/RD_VIRUS_TOTAL.Rds") #43
virus_4_1 <- readRDS("../N_V/DESeq2/ClusterStandard/Sig_4_1.Rds") #28




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



data_virus_4_1 <- data[which(rownames(data) %in% virus_4_1),]

Variables <- read_csv("~/HoneyBeePaper/Variables.csv")
View(Variables)
varIAPV <- Variables[,c(2,7)]

R2_row1 = as.matrix(cbind(as.numeric(data_virus_4_1[1,]), as.numeric(unlist(varIAPV[,2]))))
corr(R2_row1) # from boot package


