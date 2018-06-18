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
# Only 8764 of the 15,314 Rutter BeeBase IDs are in the conversion table (convert between EntrezID and BeeBase)
geneTable <- geneTable[which(geneTable$BeeBase %in% rutterAll),]

johnson <- read_excel("Table_S1.xlsx")
colJohnson <- unlist(as.data.frame(johnson[2,]))
attributes(colJohnson) <- NULL
colnames(johnson) <- colJohnson
johnson <- johnson[3:nrow(johnson),]
colnames(johnson)[1] = "Entrez"
write.csv(johnson, "johnson.csv")

###########################################################
# Boxplot background Data
# Only 8739 Johnson Entrez is in geneTable Entrez
backBox <- johnson[which(johnson$Entrez %in% geneTable$Entrez),2:ncol(johnson)]

# All 6 values for each tissue type for each gene
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
backBox$Group <- "allData"

png(paste0('All.jpg'))
print({
  ggplot(backBox, aes(x = Tissue, y = Count)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90))
})
dev.off()

###########################################################

totalGroupBox = data.frame()
i=1
GroupVec = c("noninfected", "infected", "chestnut", "rockrose", "tolerance", "resistance")
for (Group in GroupVec){
  geneGroup <- geneTable[which(geneTable$BeeBase %in% get(Group)),]
  colnames(geneGroup)[2] = "Entrez"
  johnsonGroup <- as.data.frame(johnson[which(johnson$Entrez %in% geneGroup$Entrez),])
  johnsonGroup[,2:ncol(johnsonGroup)] <- as.data.frame(sapply(johnsonGroup[,2:ncol(johnsonGroup)], as.numeric))
  
  # Boxplot foreground Group
  GroupBox <- johnsonGroup[,2:ncol(johnsonGroup)]
  GroupBox <- GroupBox[,grep("Forager", colnames(GroupBox), ignore.case=TRUE)]
  
  Antenna <- as.numeric(unlist(GroupBox[,grep("Antenna", colnames(GroupBox), ignore.case=TRUE)]))
  Ganglia <- as.numeric(unlist(GroupBox[,grep("Ganglia", colnames(GroupBox), ignore.case=TRUE)]))
  Hypopharyngeal <- as.numeric(unlist(GroupBox[,grep("Hypopharyngeal", colnames(GroupBox), ignore.case=TRUE)]))
  Mandibular <- as.numeric(unlist(GroupBox[,grep("Mandibular", colnames(GroupBox), ignore.case=TRUE)]))
  Midgut <- as.numeric(unlist(GroupBox[,grep("Midgut", colnames(GroupBox), ignore.case=TRUE)]))
  Malpighian <- as.numeric(unlist(GroupBox[,grep("Malpighian", colnames(GroupBox), ignore.case=TRUE)]))
  Muscle <- as.numeric(unlist(GroupBox[,grep("Muscle", colnames(GroupBox), ignore.case=TRUE)]))
  Nasonov <- as.numeric(unlist(GroupBox[,grep("Nasonov", colnames(GroupBox), ignore.case=TRUE)]))
  Sting <- as.numeric(unlist(GroupBox[,grep("Sting", colnames(GroupBox), ignore.case=TRUE)]))
  Brain <- as.numeric(unlist(GroupBox[,grep("Brain", colnames(GroupBox), ignore.case=TRUE)]))
  
  GroupBox2 <- data.frame(Antenna=Antenna, Ganglia=Ganglia, Hypopharyngeal=Hypopharyngeal, Mandibular=Mandibular, Midgut=Midgut, Malpighian=Malpighian, Muscle=Muscle, Nasonov=Nasonov, Sting=Sting, Brain=Brain)
  
  GroupBox <- melt(GroupBox2)
  colnames(GroupBox) <- c("Tissue", "Count")
  GroupBox$Group <- c(GroupVec[i])
  GroupBox$Count <- log(GroupBox$Count + 1)
  totalGroupBox = rbind(totalGroupBox, GroupBox)
  
  i=i+1
}

TissueVec <- c("Antenna", "Ganglia", "Hypopharyngeal", "Mandibular", "Midgut", "Malpighian", "Muscle", "Nasonov", "Sting", "Brain")
heatmapDat <- data.frame()

i=1
for (tissue in TissueVec){
  tissueBox = data.frame()
  tissueBox = rbind(tissueBox, backBox[which(backBox$Tissue==tissue),])
  tissueBox = rbind(tissueBox, totalGroupBox[which(totalGroupBox$Tissue==tissue),])
  
  tissueBox$Group = as.factor(tissueBox$Group)
  #kruskal.test(Count ~ Group, data = tissueBox)
  output <- pairwise.wilcox.test(tissueBox$Count, tissueBox$Group, p.adjust.method = "BH") #alternative="greater"
  heatmapDat <- rbind(heatmapDat, data.frame(levels(tissueBox$Group)[-1], TissueVec[i], output[[3]][1:6]))
  
  png(paste0(TissueVec[i], '_All.jpg'))
  print({
    ggplot(tissueBox, aes(x = Group, y = Count)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90)) +labs(title=TissueVec[i]) 
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
