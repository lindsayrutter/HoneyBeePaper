library(bigPint)

ID <- readRDS("../../ClusterStandard/Sig.Rds")
#load("../../Bioinformatics/Pictures/FilterNotSig/soybean_ir_noFilt_metrics.rda")

metricsAll <- readRDS("../../dataMetrics.Rds")

allMetrics = metricsAll[["N_V"]]

metricsCluster <- allMetrics[which(allMetrics$ID %in% ID),]
metricsCluster$ID <- as.character(metricsCluster$ID)

metrics <- list()
metrics[["N_V"]] <- metricsCluster

save(metrics, file = "Sigmetrics.rda")
