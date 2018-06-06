library(readr)

# tolerance <- readRDS("intNCVC.Rds") #122
# resistance <- readRDS("VCminInt.Rds") #125
# virus <- readRDS("../N_V/DESeq2/RD_VIRUS_TOTAL.Rds") #43

#data <- as.data.frame(readRDS("../N_V/data/data.Rds"))
data <- as.data.frame(readRDS("../NC_NR_VC_VR/data/data.Rds"))


# Next lines to standardize
RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
data = data %>% mutate(mean= rowMeans(as.matrix(.[nms])), stdev=rowSds(as.matrix(.[nms])))
data$mean <- as.numeric(data$mean)
rownames(data)=data_Rownames
data$ID <- data_Rownames

datas <- t(apply(as.matrix(data[,1:nColumns]), 1, scale))
datas <- as.data.frame(datas)
colnames(datas) <- colnames(data[,1:nColumns])
datas$ID <- rownames(datas)
# Indices of the NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(datas$N.1))
# Set these filtered values that have all same values for samples to 0
datas[nID,1:nColumns] <- 0

# Next two lines for BLN
# data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
# data = as.data.frame(data)

# Next lines for logging
# data <- log(data+1)

# Another approach?
# data <- as.matrix(data)
# coldata = data.frame(row.names = colnames(data), treatment = unlist(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1])))
# dds = DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ treatment)
# dds <- DESeq(dds)

rsq <- function (x, y) cor(x, y) ^ 2

Variables <- read_csv("~/HoneyBeePaper/Variables.csv")
varIAPV <- as.data.frame(Variables[,7])[,1]

boxVirus <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(boxVirus) <- c("Group", "R2")
for (i in 1:4){
  temp <- readRDS(paste0("../N_V/DESeq2/ClusterStandard/Sig_4_", i, ".Rds"))
  temp2 <- data[which(rownames(data) %in% temp),]
  temp3 = c()
  for (j in 1:nrow(temp2)){
    temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
  }
  boxVirus <- rbind(boxVirus, data.frame(Group = rep(paste0("virus",i), length(temp3)), R2 = temp3))
}

boxTolerance <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(boxTolerance) <- c("Group", "R2")
for (i in 1:4){
  temp <- readRDS(paste0("../ResistanceTolerance/Clustering_Tolerance/Sig_4_", i, ".Rds"))
  temp2 <- data[which(rownames(data) %in% temp),]
  temp3 = c()
  for (j in 1:nrow(temp2)){
    temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
  }
  boxTolerance <- rbind(boxTolerance, data.frame(Group = rep(paste0("tolerance",i), length(temp3)), R2 = temp3))
}

boxResistance <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(boxResistance) <- c("Group", "R2")
for (i in 1:4){
  temp <- readRDS(paste0("../ResistanceTolerance/Clustering_Resistance/Sig_4_", i, ".Rds"))
  temp2 <- data[which(rownames(data) %in% temp),]
  temp3 = c()
  for (j in 1:nrow(temp2)){
    temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
  }
  boxResistance <- rbind(boxResistance, data.frame(Group = rep(paste0("resistance",i), length(temp3)), R2 = temp3))
}

boxData <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(boxData) <- c("Group", "R2")
for (i in 1:nrow(data)){
  temp3[i] <- rsq(as.numeric(data[i,]), varIAPV)
}
temp3 <- temp3[-which(is.na(temp3))]
boxData <- rbind(boxData, data.frame(Group = rep(paste0("data"), length(temp3)), R2 = temp3))

plotVirus = rbind(boxVirus, boxData)
ggplot(plotVirus, aes(x=Group, y=R2)) + geom_boxplot(fill="palegreen2")

plotTolerance = rbind(boxTolerance, boxData)
ggplot(plotTolerance, aes(x=Group, y=R2)) + geom_boxplot(fill="palegreen2")

plotResistance = rbind(boxResistance, boxData)
ggplot(plotResistance, aes(x=Group, y=R2)) + geom_boxplot(fill="palegreen2")
