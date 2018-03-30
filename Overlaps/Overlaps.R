library(VennDiagram)
library(readr)

draw.triple.venn(65, 75, 85, 35, 15, 25, 5, c("First", "Second", "Third"))

#######################################################
# Get overlap between all genes between two studies
dataT = read.delim(file="../N_V/data/AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
dataG = read_delim("../VirusHoneyBee/data/GSE65659_AntiviralResponseReadCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
dataG <- as.data.frame(dataG)

# namesG has length 15314; namesT has length 11825
namesT = rownames(dataT)
namesG = unname(sapply(dataG$id, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1]))

N = intersect(namesG, namesT)
N = length(N)
#######################################################

############# Venn diagrams for Galbraith ############# 
####################################################### 

colList = scales::hue_pal()(4)

gDESeq = readRDS("../VirusHoneyBee/DESeq2/Method1/dataMetrics.Rds")
gDESeq = gDESeq[["C_T"]]
gD = gDESeq[which(gDESeq$padj <0.05),]$ID

gEdgeR = readRDS("../VirusHoneyBee/EdgeR/edgeR/dataMetrics.Rds")
gEdgeR = gEdgeR[["C_T"]]
gE = gEdgeR[which(gEdgeR$FDR <0.05),]$ID

gEdgeRB = readRDS("../VirusHoneyBee/EdgeR/edgeR-btwnLane/dataMetrics.Rds")
gEdgeRB = gEdgeRB[["C_T"]]
gEB = gEdgeRB[which(gEdgeRB$FDR <0.05),]$ID

gLimma = readRDS("../VirusHoneyBee/LimmaVoom/dataMetrics.Rds")
gLimma = gLimma[["C_T"]]
gL = gLimma[which(gLimma$adj.P.Val <0.05),]$ID

intgD = length(gD)
intgE = length(gE)
intgEB = length(gEB)
intgL = length(gL)

gDE = intersect(gD, gE)
gDEB = intersect(gD, gEB)
gDL = intersect(gD, gL)
gEEB = intersect(gE, gEB)
gEL = intersect(gE, gL)
gEBL = intersect(gEB, gL)
gDEEB = intersect(intersect(gD, gE), gEB)
gDEL = intersect(intersect(gD, gE), gL)
gDEBL = intersect(intersect(gD, gEB), gL)
gEEBL = intersect(intersect(gE, gEB), gL)
gDEEBL = intersect(intersect(gE, gEB), intersect(gD, gL))

intgDE = length(gDE)
intgDEB = length(gDEB)
intgDL = length(gDL)
intgEEB = length(gEEB)
intgEL = length(gEL)
intgEBL = length(gEBL)
intgDEEB = length(gDEEB)
intgDEL = length(gDEL)
intgDEBL = length(gDEBL)
intgEEBL = length(gEEBL)
intgDEEBL = length(gDEEBL)

fileName = paste(getwd(), "/Venn_Galbraith.jpg", sep="")
jpeg(fileName)
draw.quad.venn(add.title = "Test", area1=intgD, area2=intgE, area3=intgEB, area4=intgL, n12=intgDE, n13=intgDEB, n14=intgDL, n23=intgEEB, n24=intgEL, n34=intgEBL, n123=intgDEEB, n124=intgDEL, n134=intgDEBL, n234=intgEEBL, n1234=intgDEEBL, c("g-DESeq2", "g-EdgeR", "g-EdgeR-btwn", "g-Limma"), cat.col = colList)
invisible(dev.off())

############### Venn diagrams for Toth ############## 
##################################################### 

tDESeq = readRDS("../N_V/DESeq2/Method1/dataMetrics.Rds")
tDESeq = tDESeq[["N_V"]]
tD = tDESeq[which(tDESeq$padj <0.05),]$ID

tEdgeR = readRDS("../N_V/EdgeR/edgeR/dataMetrics.Rds")
tEdgeR = tEdgeR[["N_V"]]
tE = tEdgeR[which(tEdgeR$FDR <0.05),]$ID

tEdgeRB = readRDS("../N_V/EdgeR/edgeR-btwnLane/dataMetrics.Rds")
tEdgeRB = tEdgeRB[["N_V"]]
tEB = tEdgeRB[which(tEdgeRB$FDR <0.05),]$ID

tLimma = readRDS("../N_V/LimmaVoom/dataMetrics.Rds")
tLimma = tLimma[["N_V"]]
tL = tLimma[which(tLimma$adj.P.Val <0.05),]$ID

inttD = length(tD)
inttE = length(tE)
inttEB = length(tEB)
inttL = length(tL)

tDE = intersect(tD, tE)
tDEB = intersect(tD, tEB)
tDL = intersect(tD, tL)
tEEB = intersect(tE, tEB)
tEL = intersect(tE, tL)
tEBL = intersect(tEB, tL)
tDEEB = intersect(intersect(tD, tE), tEB)
tDEL = intersect(intersect(tD, tE), tL)
tDEBL = intersect(intersect(tD, tEB), tL)
tEEBL = intersect(intersect(tE, tEB), tL)
tDEEBL = intersect(intersect(tE, tEB), intersect(tD, tL))

inttDE = length(tDE)
inttDEB = length(tDEB)
inttDL = length(tDL)
inttEEB = length(tEEB)
inttEL = length(tEL)
inttEBL = length(tEBL)
inttDEEB = length(tDEEB)
inttDEL = length(tDEL)
inttDEBL = length(tDEBL)
inttEEBL = length(tEEBL)
inttDEEBL = length(tDEEBL)

fileName = paste(getwd(), "/Venn_Toth.jpg", sep="")
jpeg(fileName)
draw.quad.venn(area1=inttD, area2=inttE, area3=inttEB, area4=inttL, n12=inttDE, n13=inttDEB, n14=inttDL, n23=inttEEB, n24=inttEL, n34=inttEBL, n123=inttDEEB, n124=inttDEL, n134=inttDEBL, n234=inttEEBL, n1234=inttDEEBL, c("t-DESeq2", "t-EdgeR", "t-EdgeR-bwn", "t-Limma"), cat.col = colList)
invisible(dev.off())

fileName = paste(getwd(), "/Venn_Toth_NoLimma.jpg", sep="")
jpeg(fileName)
draw.triple.venn(area1=inttD, area2=inttE, area3=inttEB, n12=inttDE, n13=inttDEB, n23=inttEEB, n123=inttDEEB, c("t-DESeq2", "t-EdgeR", "t-EdgeR-bwn"), cat.col = colList[1:3])
invisible(dev.off())

############# Venn diagrams for Galbraith versus Toth per Method ############# 
##############################################################################

gD2 = sapply(gD, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
fileName = paste(getwd(), "/Venn_DESeq2.jpg", sep="")
jpeg(fileName)
DcrossArea = length(intersect(gD2, tD))
draw.pairwise.venn(area1=intgD, area2=inttD, cross.area = DcrossArea, category=c("g-DESeq2", "t-DESeq2"))
invisible(dev.off())

# Overlap between T and G for DESeq2 is significant at p-value < 2.2e-16
# 11825 | 16 | 992 | 27
fisher.test(matrix(c(N, intgD-DcrossArea, inttD-DcrossArea, DcrossArea), nrow=2), alternative="greater")

gE2 = sapply(gE, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
fileName = paste(getwd(), "/Venn_EdgeR.jpg", sep="")
jpeg(fileName)
EcrossArea = length(intersect(gE2, tE))
draw.pairwise.venn(area1=intgE, area2=inttE, cross.area = length(intersect(gE2, tE)), category=c("g-EdgeR", "t-EdgeR"))
invisible(dev.off())

# Overlap between T and G for EdgeR is significant at p-value = 9.501e-10
# 11825 | 9 | 653 | 11
fisher.test(matrix(c(N, intgE-EcrossArea, inttE-EcrossArea, EcrossArea), nrow=2), alternative="greater")

gEB2 = sapply(gEB, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
fileName = paste(getwd(), "/Venn_EdgeR-btwn.jpg", sep="")
jpeg(fileName)
EBcrossArea = length(intersect(gEB2, tEB))
draw.pairwise.venn(area1=intgEB, area2=inttEB, cross.area = length(intersect(gEB2, tEB)), category=c("g-EdgeR-btwn", "t-EdgeR-btwn"))
invisible(dev.off())

# Overlap between T and G for EdgeR-btwn is significant at p-value = 4.32e-05
# 11825 | 7 | 3317 | 13
fisher.test(matrix(c(N, intgEB-EBcrossArea, inttEB-EBcrossArea, EBcrossArea), nrow=2), alternative="greater")


