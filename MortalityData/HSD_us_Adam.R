library(plyr)
library(ggplot2)
library(multcompView)
library(readr)
library(multcomp)

d <- read_csv("day3Mortality.csv")
d = as.data.frame(d)
attr(d, "spec") <- NULL
# Note: There are 14 NC; 15 VR, VC, NR
d = d[which(d$Treatment %in% c("None Chestnut", "None Rockrose", "Virus Chestnut", "Virus Rockrose")),]
d$Treatment

d$Treatment = as.factor(d$Treatment)
colnames(d)[3] = "Day3Mortality"
levels(d$Treatment) <- c('NC', 'NR', 'VC', 'VR')

d$Virus = sapply(d$Treatment, function(x) substr(x, 1, 1))
d$Diet = sapply(d$Treatment, function(x) substr(x, 2, 2))
d$Virus = as.factor(d$Virus)
d$Diet = as.factor(d$Diet)

mortcomp2 = lme(Day3Mortality ~ Treatment, data=d, random = ~1|Experiment) 
anova(mortcomp2)
summary(glht(mortcomp2, linfct=mcp(Treatment="Tukey")))

labelDF = data.frame(plot.labels=c("NR","VC","VR","NC"), labels = c("ab","bc","c","a"), V1 = c(0.22, 0.479, 0.793, 0.164))

plotMortality = ggplot(d, aes(x=Treatment, y=Day3Mortality)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, aes(x = plot.labels, y = V1, label = labels)) + ylab("Day 3 Mortality Rate")

###### Plot IAPV Values ########

d <- read_csv("logIAPV.csv")
d = as.data.frame(d)
attr(d, "spec") <- NULL
# Note: There are 14 NC; 15 VR, VC, NR

d = d[which(d$treatment %in% c("None Chestnut", "None Rockrose", "Virus Chestnut", "Virus Rockrose")),]

d$treatment = as.factor(d$treatment)
colnames(d)[3] = "logIAPV"
levels(d$treatment) <- c('NC', 'NR', 'VC', 'VR')
colnames(d)[1] <- "Treatment"

d$Virus = sapply(d$Treatment, function(x) substr(x, 1, 1))
d$Diet = sapply(d$Treatment, function(x) substr(x, 2, 2))
d$Virus = as.factor(d$Virus)
d$Diet = as.factor(d$Diet)

mortcomp2 = lme(logIAPV ~ Treatment, data=d, random = ~1|Experiment)
anova(mortcomp2)
summary(glht(mortcomp2, linfct=mcp(Treatment="Tukey")))


labelDF = data.frame(plot.labels=c("NR","VC","VR","NC"), labels = c("a","ab","b","a"), V1 = c(4.2, 5.3, 8.52, 3.8))

plotMortality = ggplot(d, aes(x=Treatment, y=logIAPV)) + geom_boxplot(fill="paleturquoise2") + geom_text(data = labelDF, aes(x = plot.labels, y = V1, label = labels)) + ylab("Log IAPV Titer")
