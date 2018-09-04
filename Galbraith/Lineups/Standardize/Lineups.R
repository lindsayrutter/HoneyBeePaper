library(edgeR)
library(ggplot2)
library(GGally)
library(EDASeq)
library(utils)
library(data.table)
library(bigPint)

thisPath <- getwd()

beeCounts <- readRDS("../../data/data.Rds")

plotPermutations(beeCounts, nPerm = 10, topThresh = 50, option="standardize", outDir = getwd())
