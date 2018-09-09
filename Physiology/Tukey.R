library(plyr)
library(ggplot2)
library(multcompView)
library(readr)
library(multcomp)
library(nlme)

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

mortcomp = lme(Day3Mortality ~ Diet*Virus, data=d, random = ~1|Experiment) 
mortcomp2 = lme(Day3Mortality ~ Treatment, data=d, random = ~1|Experiment) 
anova(mortcomp)
summary(glht(mortcomp2, linfct=mcp(Treatment="Tukey")))

# Mortality (All)
labelDF = data.frame(plot.labels=c("NR","VC","VR","NC"), labels = c("ab","bc","c","a"), V1 = c(0.22, 0.479, 0.793, 0.164))
plotMortality = ggplot(d, aes(x=Treatment, y=Day3Mortality)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, size = 8, aes(x = plot.labels, y = V1, label = labels)) + ylab("Day 3 Mortality Rate") + theme_gray() + ylim(0,0.8) + theme(text=element_text(size=20))

# Mortality (N vs V)
labelDF = data.frame(plot.labels=c("N","V"), labels = c("a","b"), V1 = c(0.22, 0.52))
plotMortality2 = ggplot(d, aes(x=Virus, y=Day3Mortality)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, size = 8, aes(x = plot.labels, y = V1, label = labels)) + theme_gray() + ylim(0,0.8) +theme(axis.title.y=element_blank(), text=element_text(size=20))
# summary(glht(mortcomp, linfct=mcp(Virus="Tukey")))

# Mortality (R vs C)
labelDF = data.frame(plot.labels=c("C","R"), labels = c("a","b"), V1 = c(0.48, 0.54))
plotMortality3 = ggplot(d, aes(x=Diet, y=Day3Mortality)) + geom_boxplot(fill="palegreen2") + geom_text(data = labelDF, size = 8, aes(x = plot.labels, y = V1, label = labels)) + theme_gray() + ylim(0,0.8)+theme(axis.title.y=element_blank(), text=element_text(size=20))
# summary(glht(mortcomp, linfct=mcp(Diet="Tukey")))

###### Plot IAPV Values ########

d <- read_delim("logIAPV.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
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

mortcomp = lme(logIAPV ~ Diet*Virus, data=d, random = ~1|Experiment) 
mortcomp2 = lme(logIAPV ~ Treatment, data=d, random = ~1|Experiment)
anova(mortcomp)
summary(glht(mortcomp2, linfct=mcp(Treatment="Tukey")))

# IAPV (All)
labelDF = data.frame(plot.labels=c("NR","VC","VR","NC"), labels = c("a","ab","b","a"), V1 = c(4.25, 5.35, 8.57, 3.85))
plotIAPV = ggplot(d, aes(x=Treatment, y=logIAPV)) + geom_boxplot(fill="paleturquoise2") + geom_text(data = labelDF, size = 8, aes(x = plot.labels, y = V1, label = labels)) + ylab("Log IAPV Titer") +theme_gray() + theme(text=element_text(size=20)) + ylim(2,9)

# IAPV (N vs V)
labelDF = data.frame(plot.labels=c("N","V"), labels = c("a","b"), V1 = c(4.3, 8.7))
plotIAPV2 = ggplot(d, aes(x=Virus, y=logIAPV)) + geom_boxplot(fill="paleturquoise2") + geom_text(data = labelDF, size = 8, aes(x = plot.labels, y = V1, label = labels)) +theme_gray() + theme(axis.title.y=element_blank(), text=element_text(size=20)) +ylim(2,9)

# IAPV (R vs C)
labelDF = data.frame(plot.labels=c("C","R"), labels = c("a","a"), V1 = c(4.3, 8.7))
plotIAPV3 = ggplot(d, aes(x=Diet, y=logIAPV)) + geom_boxplot(fill="paleturquoise2") + geom_text(data = labelDF, size = 8, aes(x = plot.labels, y = V1, label = labels)) +theme_gray() + theme(axis.title.y=element_blank(), text=element_text(size=20)) +ylim(2,9)


