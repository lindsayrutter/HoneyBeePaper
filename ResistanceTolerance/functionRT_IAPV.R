functionRT <- function(data, type, RawPVal){
  rsq <- function (x, y) cor(x, y) ^ 2
  
  Variables <- read_csv("~/HoneyBeePaper/Variables.csv")
  varIAPV <- as.data.frame(Variables[,7])[,1]
  
  boxVirus <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxVirus) <- c("Group", "R2_IAPV")
  allVirus <- c()
  for (i in 1:4){
    temp <- readRDS(paste0("../N_V/DESeq2/ClusterStandard/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    allVirus <- c(allVirus, temp)
    if (nrow(temp2)>0){
      for (j in 1:nrow(temp2)){
        temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
      } 
    }
    boxVirus <- rbind(boxVirus, data.frame(Group = rep(paste0("virus",i), length(temp3)), R2_IAPV = temp3))
  }
  
  boxTolerance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxTolerance) <- c("Group", "R2_IAPV")
  allTolerance <- c()
  for (i in 1:4){
    temp <- readRDS(paste0("../ResistanceTolerance/Clustering_Tolerance/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    allTolerance <- c(allTolerance, temp)
    if(nrow(temp2)>0){
      for (j in 1:nrow(temp2)){
        temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
      } 
    }
    boxTolerance <- rbind(boxTolerance, data.frame(Group = rep(paste0("tolerance",i), length(temp3)), R2_IAPV = temp3))
  }
  
  boxResistance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxResistance) <- c("Group", "R2_IAPV")
  allResistance <- c()
  for (i in 1:4){
    temp <- readRDS(paste0("../ResistanceTolerance/Clustering_Resistance/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    allResistance <- c(allResistance, temp)
    if(nrow(temp2)>0){
      for (j in 1:nrow(temp2)){
        temp3[j] <- rsq(as.numeric(temp2[j,]), varIAPV)
      }
    }
    boxResistance <- rbind(boxResistance, data.frame(Group = rep(paste0("resistance",i), length(temp3)), R2_IAPV = temp3))
  }
  
  boxDataVirus <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxDataVirus) <- c("Group", "R2_IAPV")
  dataMinVirus <- data[-which(rownames(data) %in% allVirus),]
  temp3 = c()
  if (nrow(dataMinVirus)>0){
    for (i in 1:nrow(dataMinVirus)){
      temp3[i] <- rsq(as.numeric(dataMinVirus[i,]), varIAPV)
    } 
  }
  temp3 <- temp3[-which(is.na(temp3))]
  boxDataVirus <- rbind(boxDataVirus, data.frame(Group = rep(paste0("data"), length(temp3)), R2_IAPV = temp3))
  
  boxDataTolerance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxDataTolerance) <- c("Group", "R2_IAPV")
  dataMinTolerance <- data[-which(rownames(data) %in% allTolerance),]
  temp3 = c()
  if (nrow(dataMinVirus)>0){
    for (i in 1:nrow(dataMinTolerance)){
      temp3[i] <- rsq(as.numeric(dataMinTolerance[i,]), varIAPV)
    }
  }
  temp3 <- temp3[-which(is.na(temp3))]
  boxDataTolerance <- rbind(boxDataTolerance, data.frame(Group = rep(paste0("data"), length(temp3)), R2_IAPV = temp3))
  
  boxDataResistance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxDataResistance) <- c("Group", "R2_IAPV")
  dataMinResistance <- data[-which(rownames(data) %in% allResistance),]
  temp3 = c()
  if (nrow(dataMinVirus)>0){
    for (i in 1:nrow(dataMinResistance)){
      temp3[i] <- rsq(as.numeric(dataMinResistance[i,]), varIAPV)
    }
  }
  temp3 <- temp3[-which(is.na(temp3))]
  boxDataResistance <- rbind(boxDataResistance, data.frame(Group = rep(paste0("data"), length(temp3)), R2_IAPV = temp3))
  
  
  plotVirus = rbind(boxVirus, boxDataVirus)
  png(paste0('Virus_IAPV_', type, '.jpg'))
  print({
    # Test populations are not identical not assuming normality or equal variance
    kruskal.test(R2_IAPV ~ Group, data = plotVirus)
    output <- pairwise.wilcox.test(plotVirus$R2_IAPV, plotVirus$Group, p.adjust.method = "BH")
    RawPVal <- rbind(RawPVal, data.frame(Method = rep(type,4), Type = rep("virus", 4), PVal = c(signif(output[[3]][4],3), signif(output[[3]][8],3), signif(output[[3]][12],3), signif(output[[3]][16],3))))
    label1 = paste0(as.character(length(which(plotVirus$Group=="virus1"))), "\n", as.character(signif(output[[3]][4],3)))
    label2 = paste0(as.character(length(which(plotVirus$Group=="virus2"))), "\n", as.character(signif(output[[3]][8],3)))
    label3 = paste0(as.character(length(which(plotVirus$Group=="virus3"))), "\n", as.character(signif(output[[3]][12],3)))
    label4 = paste0(as.character(length(which(plotVirus$Group=="virus4"))), "\n", as.character(signif(output[[3]][16],3)))
    labelDF = data.frame(plot.labels=c("virus1","virus2","virus3","virus4","data"), labels = c(label1,label2,label3,label4,length(which(plotVirus$Group=="data"))), V1 = rep(0.5,5))
    ggplot(plotVirus, aes(x=Group, y=R2_IAPV)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, aes(x = plot.labels, y = V1, label = labels)) + ylab("R2 with IAPV titers") + theme_gray()
  })
  dev.off()  
  
  plotTolerance = rbind(boxTolerance, boxDataTolerance)
  png(paste0('Tolerance_IAPV_', type, '.jpg'))
  print({
    # Test populations are not identical not assuming normality or equal variance
    kruskal.test(R2_IAPV ~ Group, data = plotTolerance)
    output <- pairwise.wilcox.test(plotTolerance$R2_IAPV, plotTolerance$Group, p.adjust.method = "BH")
    RawPVal <- rbind(RawPVal, data.frame(Method = rep(type,4), Type = rep("tolerance", 4), PVal = c(signif(output[[3]][4],3), signif(output[[3]][8],3), signif(output[[3]][12],3), signif(output[[3]][16],3))))
    label1 = paste0(as.character(length(which(plotTolerance$Group=="tolerance1"))), "\n", as.character(signif(output[[3]][4],3)))
    label2 = paste0(as.character(length(which(plotTolerance$Group=="tolerance2"))), "\n", as.character(signif(output[[3]][8],3)))
    label3 = paste0(as.character(length(which(plotTolerance$Group=="tolerance3"))), "\n", as.character(signif(output[[3]][12],3)))
    label4 = paste0(as.character(length(which(plotTolerance$Group=="tolerance4"))), "\n", as.character(signif(output[[3]][16],3)))
    labelDF = data.frame(plot.labels=c("tolerance1","tolerance2","tolerance3","tolerance4","data"), labels = c(label1,label2,label3,label4,length(which(plotTolerance$Group=="data"))), V1 = rep(0.5,5))
    ggplot(plotTolerance, aes(x=Group, y=R2_IAPV)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, aes(x = plot.labels, y = V1, label = labels)) + ylab("R2 with IAPV titers") + theme_gray()
  })
  dev.off()  
  
  plotResistance = rbind(boxResistance, boxDataResistance)
  png(paste0('Resistance_IAPV_', type, '.jpg'))
  print({
    # Test populations are not identical not assuming normality or equal variance
    kruskal.test(R2_IAPV ~ Group, data = plotResistance)
    output <- pairwise.wilcox.test(plotResistance$R2_IAPV, plotResistance$Group, p.adjust.method = "BH")
    RawPVal <- rbind(RawPVal, data.frame(Method = rep(type,4), Type = rep("resistance", 4), PVal = c(signif(output[[3]][4],3), signif(output[[3]][8],3), signif(output[[3]][12],3), signif(output[[3]][16],3))))
    label1 = paste0(as.character(length(which(plotResistance$Group=="resistance1"))), "\n", as.character(signif(output[[3]][4],3)))
    label2 = paste0(as.character(length(which(plotResistance$Group=="resistance2"))), "\n", as.character(signif(output[[3]][8],3)))
    label3 = paste0(as.character(length(which(plotResistance$Group=="resistance3"))), "\n", as.character(signif(output[[3]][12],3)))
    label4 = paste0(as.character(length(which(plotResistance$Group=="resistance4"))), "\n", as.character(signif(output[[3]][16],3)))
    labelDF = data.frame(plot.labels=c("resistance1","resistance2","resistance3","resistance4","data"), labels = c(label1,label2,label3,label4,length(which(plotResistance$Group=="data"))), V1 = rep(0.5,5))
    ggplot(plotResistance, aes(x=Group, y=R2_IAPV)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, aes(x = plot.labels, y = V1, label = labels)) + ylab("R2 with IAPV titers") + theme_gray()
  })
  dev.off()  
return(RawPVal)
  }
