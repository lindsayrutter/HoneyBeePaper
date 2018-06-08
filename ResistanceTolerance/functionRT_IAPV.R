functionRT <- function(data, type){
  rsq <- function (x, y) cor(x, y) ^ 2
  
  Variables <- read_csv("~/HoneyBeePaper/Variables.csv")
  varIAPV <- as.data.frame(Variables[,7])[,1]
  
  boxVirus <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxVirus) <- c("Group", "R2_IAPV")
  for (i in 1:4){
    temp <- readRDS(paste0("../N_V/DESeq2/ClusterStandard/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    for (j in 1:nrow(temp2)){
      temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
    }
    boxVirus <- rbind(boxVirus, data.frame(Group = rep(paste0("virus",i), length(temp3)), R2_IAPV = temp3))
  }
  
  boxTolerance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxTolerance) <- c("Group", "R2_IAPV")
  for (i in 1:4){
    temp <- readRDS(paste0("../ResistanceTolerance/Clustering_Tolerance/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    for (j in 1:nrow(temp2)){
      temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
    }
    boxTolerance <- rbind(boxTolerance, data.frame(Group = rep(paste0("tolerance",i), length(temp3)), R2_IAPV = temp3))
  }
  
  boxResistance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxResistance) <- c("Group", "R2_IAPV")
  for (i in 1:4){
    temp <- readRDS(paste0("../ResistanceTolerance/Clustering_Resistance/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    for (j in 1:nrow(temp2)){
      temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
    }
    boxResistance <- rbind(boxResistance, data.frame(Group = rep(paste0("resistance",i), length(temp3)), R2_IAPV = temp3))
  }
  
  boxData <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxData) <- c("Group", "R2_IAPV")
  for (i in 1:nrow(data)){
    temp3[i] <- rsq(as.numeric(data[i,]), varIAPV)
  }
  temp3 <- temp3[-which(is.na(temp3))]
  boxData <- rbind(boxData, data.frame(Group = rep(paste0("data"), length(temp3)), R2_IAPV = temp3))
  
  plotVirus = rbind(boxVirus, boxData)
  png(paste0('Virus_', type, '.jpg'))
  print({
    ggplot(plotVirus, aes(x=Group, y=R2_IAPV)) + geom_boxplot(fill="palegreen2")
  })
  dev.off()  
  
  plotTolerance = rbind(boxTolerance, boxData)
  png(paste0('Tolerance_', type, '.jpg'))
  print({
    ggplot(plotTolerance, aes(x=Group, y=R2_IAPV)) + geom_boxplot(fill="palegreen2")
  })
  dev.off()  
  
  plotResistance = rbind(boxResistance, boxData)
  png(paste0('Resistance_', type, '.jpg'))
  print({
    ggplot(plotResistance, aes(x=Group, y=R2_IAPV)) + geom_boxplot(fill="palegreen2")
  })
  dev.off()  
}