functionRT <- function(data, type){
  rsq <- function (x, y) cor(x, y) ^ 2
  
  Variables <- read_csv("~/HoneyBeePaper/Variables.csv")
  varSBV <- as.data.frame(Variables[,6])[,1]
  
  boxVirus <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxVirus) <- c("Group", "R2_SBV")
  allVirus <- c()
  for (i in 1:4){
    temp <- readRDS(paste0("../N_V/DESeq2/ClusterStandard/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    allVirus <- c(allVirus, temp)
    for (j in 1:nrow(temp2)){
      temp3[j] <- rsq(as.numeric(temp2[j,]), varSBV)
    }
    boxVirus <- rbind(boxVirus, data.frame(Group = rep(paste0("virus",i), length(temp3)), R2_SBV = temp3))
  }
  
  boxTolerance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxTolerance) <- c("Group", "R2_SBV")
  allTolerance <- c()
  for (i in 1:4){
    temp <- readRDS(paste0("../ResistanceTolerance/Clustering_Tolerance/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    allTolerance <- c(allTolerance, temp)
    for (j in 1:nrow(temp2)){
      temp3[j] <- rsq(as.numeric(temp2[j,]), varSBV)
    }
    boxTolerance <- rbind(boxTolerance, data.frame(Group = rep(paste0("tolerance",i), length(temp3)), R2_SBV = temp3))
  }
  
  boxResistance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxResistance) <- c("Group", "R2_SBV")
  allResistance <- c()
  for (i in 1:4){
    temp <- readRDS(paste0("../ResistanceTolerance/Clustering_Resistance/Sig_4_", i, ".Rds"))
    temp2 <- data[which(rownames(data) %in% temp),]
    temp3 = c()
    allResistance <- c(allResistance, temp)
    for (j in 1:nrow(temp2)){
      temp3[j] <- rsq(as.numeric(temp2[j,]), varSBV)
    }
    boxResistance <- rbind(boxResistance, data.frame(Group = rep(paste0("resistance",i), length(temp3)), R2_SBV = temp3))
  }
  
  boxDataVirus <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxDataVirus) <- c("Group", "R2_SBV")
  dataMinVirus <- data[-which(rownames(data) %in% allVirus),]
  temp3 = c()
  for (i in 1:nrow(dataMinVirus)){
    temp3[i] <- rsq(as.numeric(dataMinVirus[i,]), varSBV)
  }
  temp3 <- temp3[-which(is.na(temp3))]
  boxDataVirus <- rbind(boxDataVirus, data.frame(Group = rep(paste0("data"), length(temp3)), R2_SBV = temp3))
  
  boxDataTolerance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxDataTolerance) <- c("Group", "R2_SBV")
  dataMinTolerance <- data[-which(rownames(data) %in% allTolerance),]
  temp3 = c()
  for (i in 1:nrow(dataMinTolerance)){
    temp3[i] <- rsq(as.numeric(dataMinTolerance[i,]), varSBV)
  }
  temp3 <- temp3[-which(is.na(temp3))]
  boxDataTolerance <- rbind(boxDataTolerance, data.frame(Group = rep(paste0("data"), length(temp3)), R2_SBV = temp3))
  
  boxDataResistance <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(boxDataResistance) <- c("Group", "R2_SBV")
  dataMinResistance <- data[-which(rownames(data) %in% allResistance),]
  temp3 = c()
  for (i in 1:nrow(dataMinResistance)){
    temp3[i] <- rsq(as.numeric(dataMinResistance[i,]), varSBV)
  }
  temp3 <- temp3[-which(is.na(temp3))]
  boxDataResistance <- rbind(boxDataResistance, data.frame(Group = rep(paste0("data"), length(temp3)), R2_SBV = temp3))
  
  
  plotVirus = rbind(boxVirus, boxDataVirus)
  png(paste0('Virus_SBV_', type, '.jpg'))
  print({
    # Test populations are not identical not assuming normality or equal variance
    kruskal.test(R2_SBV ~ Group, data = plotVirus)
    output <- pairwise.wilcox.test(plotVirus$R2_SBV, plotVirus$Group, p.adjust.method = "BH")
    label1 = paste0(as.character(length(which(plotVirus$Group=="virus1"))), "\n", as.character(signif(output[[3]][4],3)))
    label2 = paste0(as.character(length(which(plotVirus$Group=="virus2"))), "\n", as.character(signif(output[[3]][8],3)))
    label3 = paste0(as.character(length(which(plotVirus$Group=="virus3"))), "\n", as.character(signif(output[[3]][12],3)))
    label4 = paste0(as.character(length(which(plotVirus$Group=="virus4"))), "\n", as.character(signif(output[[3]][16],3)))
    labelDF = data.frame(plot.labels=c("virus1","virus2","virus3","virus4","data"), labels = c(label1,label2,label3,label4,length(which(plotVirus$Group=="data"))), V1 = rep(0.5,5))
    ggplot(plotVirus, aes(x=Group, y=R2_SBV)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, aes(x = plot.labels, y = V1, label = labels)) + ylab("R2 with SBV titers") + theme_gray()
  })
  dev.off()  
  
  plotTolerance = rbind(boxTolerance, boxDataTolerance)
  png(paste0('Tolerance_SBV_', type, '.jpg'))
  print({
    # Test populations are not identical not assuming normality or equal variance
    kruskal.test(R2_SBV ~ Group, data = plotTolerance)
    output <- pairwise.wilcox.test(plotTolerance$R2_SBV, plotTolerance$Group, p.adjust.method = "BH")
    label1 = paste0(as.character(length(which(plotTolerance$Group=="tolerance1"))), "\n", as.character(signif(output[[3]][4],3)))
    label2 = paste0(as.character(length(which(plotTolerance$Group=="tolerance2"))), "\n", as.character(signif(output[[3]][8],3)))
    label3 = paste0(as.character(length(which(plotTolerance$Group=="tolerance3"))), "\n", as.character(signif(output[[3]][12],3)))
    label4 = paste0(as.character(length(which(plotTolerance$Group=="tolerance4"))), "\n", as.character(signif(output[[3]][16],3)))
    labelDF = data.frame(plot.labels=c("tolerance1","tolerance2","tolerance3","tolerance4","data"), labels = c(label1,label2,label3,label4,length(which(plotTolerance$Group=="data"))), V1 = rep(0.5,5))
    ggplot(plotTolerance, aes(x=Group, y=R2_SBV)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, aes(x = plot.labels, y = V1, label = labels)) + ylab("R2 with SBV titers") + theme_gray()
  })
  dev.off()  
  
  plotResistance = rbind(boxResistance, boxDataResistance)
  png(paste0('Resistance_SBV_', type, '.jpg'))
  print({
    # Test populations are not identical not assuming normality or equal variance
    kruskal.test(R2_SBV ~ Group, data = plotResistance)
    output <- pairwise.wilcox.test(plotResistance$R2_SBV, plotResistance$Group, p.adjust.method = "BH")
    label1 = paste0(as.character(length(which(plotResistance$Group=="resistance1"))), "\n", as.character(signif(output[[3]][4],3)))
    label2 = paste0(as.character(length(which(plotResistance$Group=="resistance2"))), "\n", as.character(signif(output[[3]][8],3)))
    label3 = paste0(as.character(length(which(plotResistance$Group=="resistance3"))), "\n", as.character(signif(output[[3]][12],3)))
    label4 = paste0(as.character(length(which(plotResistance$Group=="resistance4"))), "\n", as.character(signif(output[[3]][16],3)))
    labelDF = data.frame(plot.labels=c("resistance1","resistance2","resistance3","resistance4","data"), labels = c(label1,label2,label3,label4,length(which(plotResistance$Group=="data"))), V1 = rep(0.5,5))
    ggplot(plotResistance, aes(x=Group, y=R2_SBV)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, aes(x = plot.labels, y = V1, label = labels)) + ylab("R2 with SBV titers") + theme_gray()
  })
  dev.off()  
}
