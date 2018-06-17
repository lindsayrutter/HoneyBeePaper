IAPV <- readRDS("IAPV/PVal_IAPV.Rds")
IAPV$Response = "IAPV"
IAPVIAPV[order(IAPV$PVal),]
SBV <- readRDS("SBV/PVal_SBV.Rds")
SBV$Response = "SBV"
Mort <- readRDS("Mort/PVal_Mort.Rds")
Mort$Response = "Mort"

IAPV_PValues <- IAPV[order(IAPV$PVal),]
SBV_PValues <- SBV[order(SBV$PVal),]
Mort_PValues <- Mort[order(Mort$PVal),]

All_PValues <- rbind(IAPV, SBV, Mort)

write.csv(IAPV_PValues, "IAPV_PValues.csv")
write.csv(SBV_PValues, "SBV_PValues.csv")
write.csv(Mort_PValues, "Mort_PValues.csv")
write.csv(All_PValues, "All_PValues.csv")
