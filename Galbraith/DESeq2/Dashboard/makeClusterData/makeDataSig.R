library(bigPint)

ID <- readRDS("../../ClusterStandard/Sig.Rds")
#load("../../Bioinformatics/Pictures/FilterNotSig/soybean_ir_noFilt_metrics.rda")

metricsAll <- readRDS("../../dataMetrics.Rds")

allMetrics = metricsAll[["C_T"]]

metricsCluster <- allMetrics[which(allMetrics$ID %in% ID),]
metricsCluster$ID <- as.character(metricsCluster$ID)

metrics <- list()
metrics[["C_T"]] <- metricsCluster

save(metrics, file = "Sigmetrics.rda")
