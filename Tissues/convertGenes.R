library(readxl)
library(readr)
library(reshape)

data <- readRDS("../N_V/data/data.Rds")
rutterAll <- rownames(data)

geneTable <- read_csv("am.gene_info.csv", col_names=FALSE)
colnames(geneTable)[2] = "Entrez"
colnames(geneTable)[6] = "BeeBase"
geneTable <- geneTable[which(geneTable$BeeBase %in% grep("GB", geneTable$BeeBase, value=TRUE)),]
geneTable$BeeBase <- unlist(lapply(geneTable$BeeBase, function (x) unlist(strsplit(x, "[:]"))[2]))
# Only 8764 of the 15,314 Rutter BeeBase IDs are in the conversion table
geneTable <- geneTable[which(geneTable$BeeBase %in% rutterAll),]

i=1
cluster <- readRDS(paste0("../ResistanceTolerance/Clustering_Tolerance/Sig_4_", i, ".Rds"))
geneCluster <- geneTable[which(geneTable$BeeBase %in% cluster),]
colnames(geneCluster)[2] = "Entrez"

johnson <- read_excel("Table_S1.xlsx")
colJohnson <- unlist(as.data.frame(johnson[2,]))
attributes(colJohnson) <- NULL
colnames(johnson) <- colJohnson
johnson <- johnson[3:nrow(johnson),]
colnames(johnson)[1] = "Entrez"

johnsonCluster <- as.data.frame(johnson[which(johnson$Entrez %in% geneCluster$Entrez),])
johnsonCluster[,2:ncol(johnsonCluster)] <- as.data.frame(sapply(johnsonCluster[,2:ncol(johnsonCluster)], as.numeric))

###########################################################
# Boxplot background Data
backBoxForager <- johnson[which(johnson$Entrez %in% geneTable$Entrez),2:ncol(backBoxForager)]
backBoxForager <- backBoxForager[,grep("Forager", colnames(backBoxForager), ignore.case=TRUE)]

Antenna <- as.numeric(unlist(backBoxForager[,grep("Antenna", colnames(backBoxForager), ignore.case=TRUE)]))
Ganglia <- as.numeric(unlist(backBoxForager[,grep("Ganglia", colnames(backBoxForager), ignore.case=TRUE)]))
Hypopharyngeal <- as.numeric(unlist(backBoxForager[,grep("Hypopharyngeal", colnames(backBoxForager), ignore.case=TRUE)]))
Mandibular <- as.numeric(unlist(backBoxForager[,grep("Mandibular", colnames(backBoxForager), ignore.case=TRUE)]))
Midgut <- as.numeric(unlist(backBoxForager[,grep("Midgut", colnames(backBoxForager), ignore.case=TRUE)]))
Malpighian <- as.numeric(unlist(backBoxForager[,grep("Malpighian", colnames(backBoxForager), ignore.case=TRUE)]))
Muscle <- as.numeric(unlist(backBoxForager[,grep("Muscle", colnames(backBoxForager), ignore.case=TRUE)]))
Nasonov <- as.numeric(unlist(backBoxForager[,grep("Nasonov", colnames(backBoxForager), ignore.case=TRUE)]))
Sting <- as.numeric(unlist(backBoxForager[,grep("Sting", colnames(backBoxForager), ignore.case=TRUE)]))
Brain <- as.numeric(unlist(backBoxForager[,grep("Brain", colnames(backBoxForager), ignore.case=TRUE)]))

backBoxForager2 <- data.frame(Antenna=Antenna, Ganglia=Ganglia, Hypopharyngeal=Hypopharyngeal, Mandibular=Mandibular, Midgut=Midgut, Malpighian=Malpighian, Muscle=Muscle, Nasonov=Nasonov, Sting=Sting, Brain=Brain)

backBoxForager <- melt(backBoxForager2)
colnames(backBoxForager) <- c("Tissue", "Count")
backBoxForager$Count <- log(backBoxForager$Count + 1)

ggplot(backBoxForager, aes(x = Tissue, y = Count)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90))

###########################################################
# Boxplot foreground cluster
clusterBox <- johnsonCluster[which(johnson$Entrez %in% geneTable$Entrez),2:ncol(johnsonCluster)]
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
clusterBox$Count <- log(clusterBox$Count + 1)

ggplot(clusterBox, aes(x = Tissue, y = Count)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90))
