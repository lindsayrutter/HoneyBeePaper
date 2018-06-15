library(edgeR)
library(ggplot2)
library(GGally)
library(EDASeq)
library(utils)
library(data.table)
library(bigPint)

thisPath <- getwd()

beeCounts <- readRDS("../../data/data.Rds")
beeCounts <- cbind(ID=rownames(beeCounts), beeCounts)
beeCounts$ID <- as.character(beeCounts$ID)

plotPermutationsD(beeCounts, nPerm = 20, topThresh = 50, option="log", outDir = getwd())
