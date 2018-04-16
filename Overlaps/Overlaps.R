library(VennDiagram)
library(readr)
library(readxl)

#######################################################
# Read in DEG files
RDC <- readRDS("../C_R/DESeq2/RDC.Rds")
RDR <- readRDS("../C_R/DESeq2/RDR.Rds")
RD_DIET_TOTAL <- readRDS("../C_R/DESeq2/RD_DIET_TOTAL.Rds")
REC <- readRDS("../C_R/EdgeR/edgeR/REC.Rds")
RER <- readRDS("../C_R/EdgeR/edgeR/RER.Rds")
RE_DIET_TOTAL <- readRDS("../C_R/EdgeR/edgeR/RE_DIET_TOTAL.Rds")
RLC <- readRDS("../C_R/LimmaVoom/RLC.Rds")
RLR <- readRDS("../C_R/LimmaVoom/RLR.Rds")
RL_DIET_TOTAL <- readRDS("../C_R/LimmaVoom/RL_DIET_TOTAL.Rds")
RDV <- readRDS("../N_V/DESeq2/RDV.Rds")
RDN <- readRDS("../N_V/DESeq2/RDN.Rds")
RD_VIRUS_TOTAL <- readRDS("../N_V/DESeq2/RD_VIRUS_TOTAL.Rds")
REV <- readRDS("../N_V/EdgeR/edgeR/REV.Rds")
REN <- readRDS("../N_V/EdgeR/edgeR/REN.Rds")
RE_VIRUS_TOTAL <- readRDS("../N_V/EdgeR/edgeR/RE_VIRUS_TOTAL.Rds")
GDV <- readRDS("../VirusHoneyBee/DESeq2/GDV.Rds")
GDC <- readRDS("../VirusHoneyBee/DESeq2/GDC.Rds")
GD_TOTAL <- readRDS("../VirusHoneyBee/DESeq2/GD_TOTAL.Rds")
GEV <- readRDS("../VirusHoneyBee/EdgeR/edgeR/GEV.Rds")
GEC <- readRDS("../VirusHoneyBee/EdgeR/edgeR/GEC.Rds")
GE_TOTAL <- readRDS("../VirusHoneyBee/EdgeR/edgeR/GE_TOTAL.Rds")
GLV <- readRDS("../VirusHoneyBee/LimmaVoom/GLV.Rds")
GLC <- readRDS("../VirusHoneyBee/LimmaVoom/GLC.Rds")
GL_TOTAL <- readRDS("../VirusHoneyBee/LimmaVoom/GL_TOTAL.Rds")

#######################################################
# Get overlap between all genes between two studies
# namesG has length 15314; namesR has length 15314
# Overlap has 11825

dataR = read.delim(file="../N_V/data/AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
dataG = read_delim("../VirusHoneyBee/data/GSE65659_AntiviralResponseReadCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
dataG <- as.data.frame(dataG)

namesR = rownames(dataR)
namesG = unname(sapply(dataG$id, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1]))

N = intersect(namesG, namesR)
N = length(N)

############# Venn diagrams for Galbraith (Total) ############# 

int_GD_GE_TOTAL = intersect(GD_TOTAL, GE_TOTAL)
int_GD_GL_TOTAL = intersect(GD_TOTAL, GL_TOTAL)
int_GE_GL_TOTAL = intersect(GE_TOTAL, GL_TOTAL)
int_GD_GE_GL_TOTAL = intersect(intersect(GD_TOTAL, GE_TOTAL), GL_TOTAL)

fileName = paste(getwd(), "/Venn_Galbraith_Total.jpg", sep="")
jpeg(fileName)
draw.triple.venn(area1=length(GD_TOTAL), area2=length(GE_TOTAL), area3=length(GL_TOTAL), n12=length(int_GD_GE_TOTAL), n13=length(int_GD_GL_TOTAL), n23=length(int_GE_GL_TOTAL), n123=length(int_GD_GE_GL_TOTAL), c("DESeq2", "EdgeR", "Limma"))
invisible(dev.off())

############# Venn diagrams for Galbraith (Virus) ############# 

int_GDV_GEV = intersect(GDV, GEV)
int_GDV_GLV = intersect(GDV, GLV)
int_GEV_GLV = intersect(GEV, GLV)
int_GDV_GEV_GLV = intersect(intersect(GDV, GEV), GLV)

fileName = paste(getwd(), "/Venn_Galbraith_Virus.jpg", sep="")
jpeg(fileName)
draw.triple.venn(area1=length(GDV), area2=length(GEV), area3=length(GLV), n12=length(int_GDV_GEV), n13=length(int_GDV_GLV), n23=length(int_GEV_GLV), n123=length(int_GDV_GEV_GLV), c("DESeq2", "EdgeR", "Limma"))
invisible(dev.off())

############# Venn diagrams for Galbraith (Control) ############# 

int_GDC_GEC = intersect(GDC, GEC)
int_GDC_GLC = intersect(GDC, GLC)
int_GEC_GLC = intersect(GEC, GLC)
int_GDC_GEC_GLC = intersect(intersect(GDC, GEC), GLC)

fileName = paste(getwd(), "/Venn_Galbraith_Control.jpg", sep="")
jpeg(fileName)
draw.triple.venn(area1=length(GDC), area2=length(GEC), area3=length(GLC), n12=length(int_GDC_GEC), n13=length(int_GDC_GLC), n23=length(int_GEC_GLC), n123=length(int_GDC_GEC_GLC), c("DESeq2", "EdgeR", "Limma"))
invisible(dev.off())

############### Venn diagrams for Rutter (Diet Total) ############## 

int_RD_RE_DIET_TOTAL = intersect(RD_DIET_TOTAL, RE_DIET_TOTAL)
int_RD_RL_DIET_TOTAL = intersect(RD_DIET_TOTAL, RL_DIET_TOTAL)
int_RE_RL_DIET_TOTAL = intersect(RE_DIET_TOTAL, RL_DIET_TOTAL)
int_RD_RE_RL_DIET_TOTAL = intersect(intersect(RD_DIET_TOTAL, RE_DIET_TOTAL), RL_DIET_TOTAL)

fileName = paste(getwd(), "/Venn_Rutter_Diet_Total.jpg", sep="")
jpeg(fileName)
draw.triple.venn(area1=length(RD_DIET_TOTAL), area2=length(RE_DIET_TOTAL), area3=length(RL_DIET_TOTAL), n12=length(int_RD_RE_DIET_TOTAL), n13=length(int_RD_RL_DIET_TOTAL), n23=length(int_RE_RL_DIET_TOTAL), n123=length(int_RD_RE_RL_DIET_TOTAL), c("DESeq2", "EdgeR", "Limma"))
invisible(dev.off())

############### Venn diagrams for Rutter (Diet Castanea) ############## 

int_RDC_REC = intersect(RDC, REC)
int_RDC_RLC = intersect(RDC, RLC)
int_REC_RLC = intersect(REC, RLC)
int_RDC_REC_RLC = intersect(intersect(RDC, REC), RLC)

fileName = paste(getwd(), "/Venn_Rutter_Diet_Castanea.jpg", sep="")
jpeg(fileName)
draw.triple.venn(area1=length(RDC), area2=length(REC), area3=length(RLC), n12=length(int_RDC_REC), n13=length(int_RDC_RLC), n23=length(int_REC_RLC), n123=length(int_RDC_REC_RLC), c("DESeq2", "EdgeR", "Limma"))
invisible(dev.off())

############### Venn diagrams for Rutter (Diet Rockrose) ############## 

int_RDR_RER = intersect(RDR, RER)
int_RDR_RLR = intersect(RDR, RLR)
int_RER_RLR = intersect(RER, RLR)
int_RDR_RER_RLR = intersect(intersect(RDR, RER), RLR)

fileName = paste(getwd(), "/Venn_Rutter_Diet_Rockrose.jpg", sep="")
jpeg(fileName)
draw.triple.venn(area1=length(RDR), area2=length(RER), area3=length(RLR), n12=length(int_RDR_RER), n13=length(int_RDR_RLR), n23=length(int_RER_RLR), n123=length(int_RDR_RER_RLR), c("DESeq2", "EdgeR", "Limma"))
invisible(dev.off())

############### Rutter Diet (overlaps between up and down) ############## 

# No cases where upregulated in one pipeline (DESeq, EdgeR, Limma) and downregulated in another
length(intersect(RDC, RER))
length(intersect(RDC, RLR))
length(intersect(REC, RDR))
length(intersect(REC, RLR))
length(intersect(RLC, RDR))
length(intersect(RLC, RER))



############### Venn diagrams for Rutter (Virus Total) ############## 

int_RD_RE_VIRUS_TOTAL = intersect(RD_VIRUS_TOTAL, RE_VIRUS_TOTAL)

fileName = paste(getwd(), "/Venn_Rutter_Virus_Total.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(RD_VIRUS_TOTAL), area2=length(RE_VIRUS_TOTAL), cross.area=length(int_RD_RE_VIRUS_TOTAL), c("DESeq2", "EdgeR"))
invisible(dev.off())

############### Venn diagrams for Rutter (Virus Infected) ############## 

int_RDV_REV = intersect(RDV, REV)

fileName = paste(getwd(), "/Venn_Rutter_Virus_Infected.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(RDV), area2=length(REV), cross.area=length(int_RDV_REV), c("DESeq2", "EdgeR"))
invisible(dev.off())

############### Venn diagrams for Rutter (Virus Healthy) ############## 

int_RDN_REN = intersect(RDN, REN)

fileName = paste(getwd(), "/Venn_Rutter_Virus_Healthy.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(RDN), area2=length(REN), cross.area=length(int_RDN_REN), c("DESeq2", "EdgeR"))
invisible(dev.off())

############### Rutter Virus (overlaps between up and down) ############## 

# No cases where upregulated in one pipeline (DESeq, EdgeR, Limma) and downregulated in another
length(intersect(RDV, REN))
length(intersect(REV, RDN))












############### Venn diagrams for Rutter and Galbraith, DEL (Total) ##############

int_GD_RD_TOTAL = intersect(GD_TOTAL, RD_VIRUS_TOTAL)
int_GD_RE_TOTAL = intersect(GD_TOTAL, RE_VIRUS_TOTAL)
int_GE_RD_TOTAL = intersect(GE_TOTAL, RD_VIRUS_TOTAL)
int_GE_RE_TOTAL = intersect(GE_TOTAL, RE_VIRUS_TOTAL)
int_GL_RD_TOTAL = intersect(GL_TOTAL, RD_VIRUS_TOTAL)
int_GL_RE_TOTAL = intersect(GL_TOTAL, RE_VIRUS_TOTAL)

int_GD_GE_RD_TOTAL = Reduce(intersect, list(GD_TOTAL, GE_TOTAL, RD_VIRUS_TOTAL))
int_GD_GE_RE_TOTAL = Reduce(intersect, list(GD_TOTAL, GE_TOTAL, RE_VIRUS_TOTAL))
int_GD_GL_RD_TOTAL = Reduce(intersect, list(GD_TOTAL, GL_TOTAL, RD_VIRUS_TOTAL))
int_GD_GL_RE_TOTAL = Reduce(intersect, list(GD_TOTAL, GL_TOTAL, RE_VIRUS_TOTAL))
int_GD_RD_RE_TOTAL = Reduce(intersect, list(GD_TOTAL, RD_VIRUS_TOTAL, RE_VIRUS_TOTAL))
int_GE_GL_RD_TOTAL = Reduce(intersect, list(GE_TOTAL, GL_TOTAL, RD_VIRUS_TOTAL))
int_GE_GL_RE_TOTAL = Reduce(intersect, list(GE_TOTAL, GL_TOTAL, RE_VIRUS_TOTAL))
int_GE_RD_RE_TOTAL = Reduce(intersect, list(GE_TOTAL, RD_VIRUS_TOTAL, RE_VIRUS_TOTAL))
int_GL_RD_RE_TOTAL = Reduce(intersect, list(GL_TOTAL, RD_VIRUS_TOTAL, RE_VIRUS_TOTAL))

int_GD_GE_GL_RD_TOTAL = Reduce(intersect, list(GD_TOTAL, GE_TOTAL, GL_TOTAL, RD_VIRUS_TOTAL))
int_GD_GE_GL_RE_TOTAL = Reduce(intersect, list(GD_TOTAL, GE_TOTAL, GL_TOTAL, RE_VIRUS_TOTAL))
int_GD_GE_RD_RE_TOTAL = Reduce(intersect, list(GD_TOTAL, GE_TOTAL, RD_VIRUS_TOTAL, RE_VIRUS_TOTAL))
int_GD_GL_RD_RE_TOTAL = Reduce(intersect, list(GD_TOTAL, GL_TOTAL, RD_VIRUS_TOTAL, RE_VIRUS_TOTAL))
int_GE_GL_RD_RE_TOTAL = Reduce(intersect, list(GE_TOTAL, GL_TOTAL, RD_VIRUS_TOTAL, RE_VIRUS_TOTAL))
int_GD_GE_GL_RD_RE_TOTAL = Reduce(intersect, list(GD_TOTAL, GE_TOTAL, GL_TOTAL, RD_VIRUS_TOTAL, RE_VIRUS_TOTAL))

# Reference five-set diagram
venn.plot <- draw.quintuple.venn(
  area1 = length(GD_TOTAL),
  area2 = length(GE_TOTAL),
  area3 = length(GL_TOTAL),
  area4 = length(RD_VIRUS_TOTAL),
  area5 = length(RE_VIRUS_TOTAL),
  n12 = length(int_GD_GE_TOTAL),
  n13 = length(int_GD_GL_TOTAL),
  n14 = length(int_GD_RD_TOTAL),
  n15 = length(int_GD_RE_TOTAL),
  n23 = length(int_GE_GL_TOTAL),
  n24 = length(int_GE_RD_TOTAL),
  n25 = length(int_GE_RE_TOTAL),
  n34 = length(int_GL_RD_TOTAL),
  n35 = length(int_GL_RE_TOTAL),
  n45 = length(int_RD_RE_VIRUS_TOTAL),
  n123 = length(int_GD_GE_GL_TOTAL),
  n124 = length(int_GD_GE_RD_TOTAL),
  n125 = length(int_GD_GE_RE_TOTAL),
  n134 = length(int_GD_GL_RD_TOTAL),
  n135 = length(int_GD_GL_RE_TOTAL),
  n145 = length(int_GD_RD_RE_TOTAL),
  n234 = length(int_GE_GL_RD_TOTAL),
  n235 = length(int_GE_GL_RE_TOTAL),
  n245 = length(int_GE_RD_RE_TOTAL),
  n345 = length(int_GL_RD_RE_TOTAL),
  n1234 = length(int_GD_GE_GL_RD_TOTAL),
  n1235 = length(int_GD_GE_GL_RE_TOTAL),
  n1245 = length(int_GD_GE_RD_RE_TOTAL),
  n1345 = length(int_GD_GL_RD_RE_TOTAL),
  n2345 = length(int_GE_GL_RD_RE_TOTAL),
  n12345 = length(int_GD_GE_GL_RD_RE_TOTAL),
  category = c("", "", " ", " ", " "),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
);

tiff(filename = "Quintuple_Venn_diagram.jpg", compression = "lzw");
grid.draw(venn.plot);
dev.off();










############### Venn diagrams for Supplemental Galbraith versus Toth per Method ############## 

# Keep 753 DEGs from Galbraith paper
gDSupp = read_excel("../VirusHoneyBee/GalbraithSuppFile.xlsx")
gDSupp = as.data.frame(gDSupp)

keepColNms <- gDSupp[3,]
keepColNms <- as.data.frame(keepColNms)
colnames(keepColNms) <- NULL
keepColNms <- unlist(keepColNms)

keepDEGs <- gDSupp[4:(4+753-1),]
keepDEGs <- as.data.frame(keepDEGs)
colnames(keepDEGs) <- keepColNms

# Extract 753 IDs from Galbraith paper
gDSupp = keepDEGs$id
gDSupp = unname(sapply(gDSupp, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1]))
gtDSupp = intersect(gDSupp, GD_TOTAL)

fileName = paste(getwd(), "/Venn_Galbraith_GalbraithSupp_Virus.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GD_TOTAL), area2=length(gDSupp), cross.area = length(gtDSupp), category=c("Galbraith Us", "Galbraith Supp"))
invisible(dev.off())

# Extract 607 Upregeulated IDs from Galbraith paper
gDSuppUp = keepDEGs[which(keepDEGs$`up/down`=="+"),]$id
gDSuppUp = unname(sapply(gDSuppUp, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1]))
gtDSuppUpp = intersect(gDSuppUp, GDV)

fileName = paste(getwd(), "/Venn_Galbraith_GalbraithSupp_VirusUp.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GDV), area2=length(gDSuppUp), cross.area = length(gtDSuppUpp), category=c("Galbraith Us", "Galbraith Supp"))
invisible(dev.off())

# Extract 146 Upregeulated IDs from Galbraith paper
gDSuppDown = keepDEGs[which(keepDEGs$`up/down`=="-"),]$id
gDSuppDown = unname(sapply(gDSuppDown, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1]))
gtDSuppDown = intersect(gDSuppDown, GDC)

fileName = paste(getwd(), "/Venn_Galbraith_GalbraithSupp_VirusDown.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GDC), area2=length(gDSuppDown), cross.area = length(gtDSuppDown), category=c("Galbraith Us", "Galbraith Supp"))
invisible(dev.off())

# No cases where upregulated in one study (Galbriath, Rutter) and downregulated in another
length(intersect(gDSuppUp, GDC))
length(intersect(gDSuppDown, GDV))


############### Venn diagrams for DESeq (Total) Rutter vs. Galbraith ##############

fileName = paste(getwd(), "/Venn_GR_DESeq_Total.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GD_TOTAL), area2=length(RD_VIRUS_TOTAL), cross.area = length(intersect(RD_VIRUS_TOTAL, GD_TOTAL)), category=c("Galbraith", "Rutter"))
invisible(dev.off())

############### Venn diagrams for DESeq (Virus) Rutter vs. Galbraith ##############

fileName = paste(getwd(), "/Venn_GR_DESeq_Virus.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GDV), area2=length(RDV), cross.area = length(intersect(RDV, GDV)), category=c("Galbraith", "Rutter"))
invisible(dev.off())

############### Venn diagrams for DESeq (Healthy) Rutter vs. Galbraith ##############

fileName = paste(getwd(), "/Venn_GR_DESeq_Healthy.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GDC), area2=length(RDN), cross.area = length(intersect(RDN, GDC)), category=c("Galbraith", "Rutter"))
invisible(dev.off())

############### Rutter vs. Galbriath DESeq (overlaps between up and down) ############## 

length(intersect(RDV, GDC))
length(intersect(RDN, GDV)) # One overlap! ("GB51305")







############### Venn diagrams for EdgeR (Total) Rutter vs. Galbraith ##############

fileName = paste(getwd(), "/Venn_GR_EdgeR_Total.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GE_TOTAL), area2=length(RE_VIRUS_TOTAL), cross.area = length(intersect(RE_VIRUS_TOTAL, GE_TOTAL)), category=c("Galbraith", "Rutter"))
invisible(dev.off())

############### Venn diagrams for EdgeR (Virus) Rutter vs. Galbraith ##############

fileName = paste(getwd(), "/Venn_GR_EdgeR_Virus.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GEV), area2=length(REV), cross.area = length(intersect(REV, GEV)), category=c("Galbraith", "Rutter"))
invisible(dev.off())

############### Venn diagrams for EdgeR (Healthy) Rutter vs. Galbraith ##############

fileName = paste(getwd(), "/Venn_GR_EdgeR_Healthy.jpg", sep="")
jpeg(fileName)
draw.pairwise.venn(area1=length(GEC), area2=length(REN), cross.area = length(intersect(REN, GEC)), category=c("Galbraith", "Rutter"))
invisible(dev.off())

############### Rutter vs. Galbriath DESeq (overlaps between up and down) ############## 

length(intersect(REV, GDC))
length(intersect(GEV, REN)) # One overlap! ("GB51305")




###################### Save GO Terms ##################### 
##########################################################
  
# Save GO 
gD = sapply(gD, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gE = sapply(gE, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gEB = sapply(gEB, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gL = sapply(gL, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gDE = sapply(gDE, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gDEB = sapply(gDEB, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gDL = sapply(gDL, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gEEB = sapply(gEEB, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gEL = sapply(gEL, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gEBL = sapply(gEBL, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gDEEB = sapply(gDEEB, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gDEL = sapply(gDEL, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gDEBL = sapply(gDEBL, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gEEBL = sapply(gEEBL, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])
gDEEBL = sapply(gDEEBL, function(x) strsplit(strsplit(x, "[|]")[[1]][3], "[-]")[[1]][1])

saveRDS(unname(gD), "GO/gD.rds")
saveRDS(unname(gE), "GO/gE.rds")
saveRDS(unname(gEB), "GO/gEB.rds")
saveRDS(unname(gL), "GO/gL.rds")
saveRDS(unname(gDE), "GO/gDE.rds")
saveRDS(unname(gDEB), "GO/gDEB.rds")
saveRDS(unname(gDL), "GO/gDL.rds")
saveRDS(unname(gEEB), "GO/gEEB.rds")
saveRDS(unname(gEL), "GO/gEL.rds")
saveRDS(unname(gEBL), "GO/gEBL.rds")
saveRDS(unname(gDEEB), "GO/gDEEB.rds")
saveRDS(unname(gDEL), "GO/gDEL.rds")
saveRDS(unname(gDEBL), "GO/gDEBL.rds")
saveRDS(unname(gEEBL), "GO/gEEBL.rds")
saveRDS(unname(gDEEBL), "GO/gDEEBL.rds")

saveRDS(tD, "GO/tD.rds")
saveRDS(tE, "GO/tE.rds")
saveRDS(tEB, "GO/tEB.rds")
saveRDS(tL, "GO/tL.rds")
saveRDS(tDE, "GO/tDE.rds")
saveRDS(tDEB, "GO/tDEB.rds")
saveRDS(tDL, "GO/tDL.rds")
saveRDS(tEEB, "GO/tEEB.rds")
saveRDS(tEL, "GO/tEL.rds")
saveRDS(tEBL, "GO/tEBL.rds")
saveRDS(tDEEB, "GO/tDEEB.rds")
saveRDS(tDEL, "GO/tDEL.rds")
saveRDS(tDEBL, "GO/tDEBL.rds")
saveRDS(tEEBL, "GO/tEEBL.rds")
saveRDS(tDEEBL, "GO/tDEEBL.rds")

tgDInt = intersect(gD2, tD)
tgDUnion = union(gD2, tD)
saveRDS(tgDInt, "GO/tgDInt.rds")
saveRDS(tgDUnion, "GO/tgDUnion.rds")

tgEInt = intersect(gE2, tE)
tgEUnion = union(gE2, tE)
saveRDS(tgEInt, "GO/tgEInt.rds")
saveRDS(tgEUnion, "GO/tgEUnion.rds")

tgEBInt = intersect(gEB2, tEB)
tgEBUnion = union(gEB2, tEB)
saveRDS(tgEBInt, "GO/tgEBInt.rds")
saveRDS(tgEBUnion, "GO/tgEBUnion.rds")

