library(readxl)
library(readr)
library(reshape)

noninfected <- readRDS("../N_V/DESeq2/RDN.Rds")
infected <- readRDS("../N_V/DESeq2/RDV.Rds")
chestnut <- readRDS("../C_R/DESeq2/RDC.Rds")
rockrose <- readRDS("../C_R/DESeq2/RDR.Rds")
tolerance <- readRDS("../ResistanceTolerance/tolerance.Rds")
resistance <- readRDS("../ResistanceTolerance/resistance.Rds")

data <- readRDS("../N_V/data/data.Rds")
rutterAll <- rownames(data)

geneTable <- read_csv("am.gene_info.csv", col_names=FALSE)
colnames(geneTable)[2] = "Entrez"
colnames(geneTable)[6] = "BeeBase"
geneTable <- geneTable[which(geneTable$BeeBase %in% grep("GB", geneTable$BeeBase, value=TRUE)),]
geneTable$BeeBase <- unlist(lapply(geneTable$BeeBase, function (x) unlist(strsplit(x, "[:]"))[2]))
# Only 8764 of the 15,314 Rutter BeeBase IDs are in the conversion table
geneTable <- geneTable[which(geneTable$BeeBase %in% rutterAll),]

johnson <- read_excel("Table_S1.xlsx")
colJohnson <- unlist(as.data.frame(johnson[2,]))
attributes(colJohnson) <- NULL
colnames(johnson) <- colJohnson
johnson <- johnson[3:nrow(johnson),]
colnames(johnson)[1] = "Entrez"

###########################################################
# Boxplot background Data
backBox <- johnson[which(johnson$Entrez %in% geneTable$Entrez),2:ncol(johnson)]

Antenna <- as.numeric(unlist(backBox[,grep("Antenna", colnames(backBox), ignore.case=TRUE)]))
Ganglia <- as.numeric(unlist(backBox[,grep("Ganglia", colnames(backBox), ignore.case=TRUE)]))
Hypopharyngeal <- as.numeric(unlist(backBox[,grep("Hypopharyngeal", colnames(backBox), ignore.case=TRUE)]))
Mandibular <- as.numeric(unlist(backBox[,grep("Mandibular", colnames(backBox), ignore.case=TRUE)]))
Midgut <- as.numeric(unlist(backBox[,grep("Midgut", colnames(backBox), ignore.case=TRUE)]))
Malpighian <- as.numeric(unlist(backBox[,grep("Malpighian", colnames(backBox), ignore.case=TRUE)]))
Muscle <- as.numeric(unlist(backBox[,grep("Muscle", colnames(backBox), ignore.case=TRUE)]))
Nasonov <- as.numeric(unlist(backBox[,grep("Nasonov", colnames(backBox), ignore.case=TRUE)]))
Sting <- as.numeric(unlist(backBox[,grep("Sting", colnames(backBox), ignore.case=TRUE)]))
Brain <- as.numeric(unlist(backBox[,grep("Brain", colnames(backBox), ignore.case=TRUE)]))

backBox2 <- data.frame(Antenna=Antenna, Ganglia=Ganglia, Hypopharyngeal=Hypopharyngeal, Mandibular=Mandibular, Midgut=Midgut, Malpighian=Malpighian, Muscle=Muscle, Nasonov=Nasonov, Sting=Sting, Brain=Brain)

backBox <- melt(backBox2)
colnames(backBox) <- c("Tissue", "Count")
backBox$Count <- log(backBox$Count + 1)
backBox$Cluster <- "allData"

png(paste0('All.jpg'))
print({
  ggplot(backBox, aes(x = Tissue, y = Count)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90))
})
dev.off()

###########################################################

totalClusterBox = data.frame()
i=1
clusterVec = c("noninfected", "infected", "chestnut", "rockrose", "tolerance", "resistance")
for (cluster in clusterVec){
  geneCluster <- geneTable[which(geneTable$BeeBase %in% get(cluster)),]
  colnames(geneCluster)[2] = "Entrez"
  johnsonCluster <- as.data.frame(johnson[which(johnson$Entrez %in% geneCluster$Entrez),])
  johnsonCluster[,2:ncol(johnsonCluster)] <- as.data.frame(sapply(johnsonCluster[,2:ncol(johnsonCluster)], as.numeric))
  
  # Boxplot foreground cluster
  clusterBox <- johnsonCluster[,2:ncol(johnsonCluster)]
  clusterBox <- clusterBox[,grep("Forager", colnames(clusterBox), ignore.case=TRUE)]
  
  Antenna <- as.numeric(unlist(clusterBox[,grep("Antenna", colnames(clusterBox), ignore.case=TRUE)]))
  Ganglia <- as.numeric(unlist(clusterBox[,grep("Ganglia", colnames(clusterBox), ignore.case=TRUE)]))
  Hypopharyngeal <- as.numeric(unlist(clusterBox[,grep("Hypopharyngeal", colnames(clusterBox), ignore.case=TRUE)]))
  Mandibular <- as.numeric(unlist(clusterBox[,grep("Mandibular", colnames(clusterBox), ignore.case=TRUE)]))
  Midgut <- as.numeric(unlist(clusterBox[,grep("Midgut", colnames(clusterBox), ignore.case=TRUE)]))
  Malpighian <- as.numeric(unlist(clusterBox[,grep("Malpighian", colnames(clusterBox), ignore.case=TRUE)]))
  Muscle <- as.numeric(unlist(clusterBox[,grep("Muscle", colnames(clusterBox), ignore.case=TRUE)]))
  Nasonov <- as.numeric(unlist(clusterBox[,grep("Nasonov", colnames(clusterBox), ignore.case=TRUE)]))
  Sting <- as.numeric(unlist(clusterBox[,grep("Sting", colnames(clusterBox), ignore.case=TRUE)]))
  Brain <- as.numeric(unlist(clusterBox[,grep("Brain", colnames(clusterBox), ignore.case=TRUE)]))
  
  clusterBox2 <- data.frame(Antenna=Antenna, Ganglia=Ganglia, Hypopharyngeal=Hypopharyngeal, Mandibular=Mandibular, Midgut=Midgut, Malpighian=Malpighian, Muscle=Muscle, Nasonov=Nasonov, Sting=Sting, Brain=Brain)
  
  clusterBox <- melt(clusterBox2)
  colnames(clusterBox) <- c("Tissue", "Count")
  clusterBox$Cluster <- c(clusterVec[i])
  clusterBox$Count <- log(clusterBox$Count + 1)
  totalClusterBox = rbind(totalClusterBox, clusterBox)
  
  i=i+1
}

TissueVec <- c("Antenna", "Ganglia", "Hypopharyngeal", "Mandibular", "Midgut", "Malpighian", "Muscle", "Nasonov", "Sting", "Brain")
heatmapDat <- data.frame()

i=1
for (tissue in TissueVec){
  tissueBox = data.frame()
  tissueBox = rbind(tissueBox, backBox[which(backBox$Tissue==tissue),])
  tissueBox = rbind(tissueBox, totalClusterBox[which(totalClusterBox$Tissue==tissue),])
  
  tissueBox$Cluster = as.factor(tissueBox$Cluster)
  kruskal.test(Count ~ Cluster, data = tissueBox)
  output <- pairwise.wilcox.test(tissueBox$Count, tissueBox$Cluster, p.adjust.method = "BH")
  heatmapDat <- rbind(heatmapDat, data.frame(levels(tissueBox$Cluster)[-1], TissueVec[i], output[[3]][1:6]))
  
  png(paste0(TissueVec[i], '_All.jpg'))
  print({
    ggplot(tissueBox, aes(x = Cluster, y = Count)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90)) +labs(title=TissueVec[i]) 
  })
  dev.off()
  i=i+1
}

colnames(heatmapDat) <- c("Group", "Tissue", "pvalue")
h <- ggplot(heatmapDat, aes(Group, Tissue )) + geom_tile(aes(fill = pvalue), color = "white") + scale_fill_gradient(low = "black", high = "white", trans = "log") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 12), plot.title = element_text(size=16), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 90, hjust = 1)) + labs(fill = "kruskalPvalue") + ggtitle("All")

png(paste0('Heatmap_All.jpg'))
print({
  h 
})
dev.off()
