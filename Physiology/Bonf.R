library(agricolae)

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

# Mortality (All)
model<-aov(Day3Mortality~Treatment, data=d)
out <- LSD.test(model,"Treatment", p.adj="bonferroni")
# Variation range: max and min
plot(out, ylab="Mortality")
# Used to print the actual p-values
anova(model)

# Mortality (N vs V)
model<-aov(Day3Mortality~Virus, data=d)
out <- LSD.test(model,"Virus", p.adj="bonferroni")
# Variation range: max and min
plot(out, ylab="Mortality")
# Used to print the actual p-values
anova(model)

# Mortality (R vs C)
model<-aov(Day3Mortality~Diet, data=d)
out <- LSD.test(model,"Diet", p.adj="bonferroni")
# Variation range: max and min
plot(out, ylab="Mortality")
# Used to print the actual p-values
anova(model)

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

# IAPV (All)
model<-aov(logIAPV~Treatment, data=d)
out <- LSD.test(model,"Treatment", p.adj="bonferroni")
# Variation range: max and min
plot(out, ylab="Log IAPV")
# Used to print the actual p-values
anova(model)

# IAPV (N vs V)
model<-aov(logIAPV~Virus, data=d)
out <- LSD.test(model,"Virus", p.adj="bonferroni")
# Variation range: max and min
plot(out, ylab="Log IAPV")
# Used to print the actual p-values
anova(model)

# IAPV (R vs C)
model<-aov(logIAPV~Diet, data=d)
out <- LSD.test(model,"Diet", p.adj="bonferroni")
# Variation range: max and min
plot(out, ylab="Log IAPV")
# Used to print the actual p-values
anova(model)

