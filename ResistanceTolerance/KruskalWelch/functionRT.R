functionRT <- function(data, type, colNum, PVal){
  rsq <- function (x, y) cor(x, y) ^ 2
  
  Variables <- read_csv("~/HoneyBeePaper/Variables.csv")
  myVar <- as.data.frame(Variables[,colNum])[,1]
  
  readFile <- c("../../N_V/DESeq2/ClusterStandard/Sig_4_", "../../ResistanceTolerance/Clustering_Tolerance/Sig_4_", "../../ResistanceTolerance/Clustering_Resistance/Sig_4_")
  strVar <- c("virus", "tolerance", "resistance")
  
  getR2 <- function(rf, sv){
    dfR2 <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(dfR2) <- c("Group", "R2")
    allR2 <- c()
    for (i in 1:4){
      temp <- readRDS(paste0(rf, i, ".Rds")) # different
      temp2 <- data[which(rownames(data) %in% temp),]
      temp3 = c()
      allR2 <- c(allR2, temp)
      if (nrow(temp2)>0){
        for (j in 1:nrow(temp2)){
          temp3[j] <- rsq(as.numeric(temp2[j,]), myVar)
        } 
      }
      dfR2 <- rbind(dfR2, data.frame(Group = rep(paste0("Cluster",i), length(temp3)), R2 = temp3))
    }
    dfR2Box <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(dfR2Box) <- c("Group", "R2")
    dataMin <- data[-which(rownames(data) %in% allR2),]
    temp3 = c()
    if (nrow(dataMin)>0){
      for (i in 1:nrow(dataMin)){
        temp3[i] <- rsq(as.numeric(dataMin[i,]), myVar)
      } 
    }
    temp3 <- temp3[-which(is.na(temp3))]
    dfR2Box <- rbind(dfR2Box, data.frame(Group = rep(paste0("data"), length(temp3)), R2 = temp3))
    list(dfR2, dfR2Box)
  }
  
  virusR2 <- getR2(readFile[1], strVar[1])
  plotVirus <- rbind(virusR2[[1]], virusR2[[2]])
  toleranceR2 <- getR2(readFile[2], strVar[2])
  plotTolerance <- rbind(toleranceR2[[1]], toleranceR2[[2]])
  resistanceR2 <- getR2(readFile[3], strVar[3])
  plotResistance <- rbind(resistanceR2[[1]], resistanceR2[[2]])
  
  makePlot <- function(inputR2, degGroup){
    png(paste0(degGroup, type, '.jpg'))
    print({
      # Kruskal test
      kruskal.test(R2 ~ Group, data = inputR2)
      output <- pairwise.wilcox.test(inputR2$R2, inputR2$Group, p.adjust.method = "BH")
      # Welch test
      welch1 = welch.test(R2 ~ Group, data = inputR2[which(inputR2$Group %in% c("Cluster1", "data")),])
      welch2 = welch.test(R2 ~ Group, data = inputR2[which(inputR2$Group %in% c("Cluster2", "data")),])
      welch3 = welch.test(R2 ~ Group, data = inputR2[which(inputR2$Group %in% c("Cluster3", "data")),])
      welch4 = welch.test(R2 ~ Group, data = inputR2[which(inputR2$Group %in% c("Cluster4", "data")),])
      
      PVal <- rbind(PVal, data.frame(Method = rep(degGroup,4), Type = rep(type, 4), Cluster = c(1,2,3,4), PKruskal = c(signif(output[[3]][4],3), signif(output[[3]][8],3), signif(output[[3]][12],3), signif(output[[3]][16],3)), PWelch = c(welch1[[3]], welch2[[3]], welch3[[3]], welch4[[3]])))
      
      label1 = paste0(as.character(length(which(inputR2$Group=="Cluster1"))), "\n", as.character(signif(output[[3]][4],3)))
      label2 = paste0(as.character(length(which(inputR2$Group=="Cluster2"))), "\n", as.character(signif(output[[3]][8],3)))
      label3 = paste0(as.character(length(which(inputR2$Group=="Cluster3"))), "\n", as.character(signif(output[[3]][12],3)))
      label4 = paste0(as.character(length(which(inputR2$Group=="Cluster4"))), "\n", as.character(signif(output[[3]][16],3)))
      labelDF = data.frame(plot.labels=c("Cluster1","Cluster2","Cluster3","Cluster4","data"), labels = c(label1,label2,label3,label4,length(which(inputR2$Group=="data"))), V1 = rep(0.5,5))
      ggplot(inputR2, aes(x=Group, y=R2)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, size=6, aes(x = plot.labels, y = V1, label = labels)) + ylab(paste0("R2 with", type)) + theme_gray() + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.position="none")
    })
    dev.off()
    PVal
  }
  PVal <- makePlot(plotVirus, "virus")
  PVal <- makePlot(plotTolerance, "tolerance")
  PVal <- makePlot(plotResistance, "resistance")
  PVal
}

