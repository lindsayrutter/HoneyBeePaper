library(plyr)
library(ggplot2)
library(multcompView)
library(readr)

d <- read_csv("day3Mortality.csv")
d = as.data.frame(d)
attr(d, "spec") <- NULL
# Note: There are 14 NC; 15 VR, VC, NR
d = d[which(d$Treatment %in% c("None Chestnut", "None Rockrose", "Virus Chestnut", "Virus Rockrose")),]
d$Treatment

d$Treatment = as.factor(d$Treatment)
colnames(d)[2] = "Day3Mortality"
levels(d$Treatment) <- c('NC', 'NR', 'VC', 'VR')

a <- aov(Day3Mortality~Treatment, data=d)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)

generate_label_df1 <- function(HSD, flev){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[flev]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$Day3Mortality)) + 0.05)
  boxplot.df[2,]$V1 = 0.22
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}

plotMortality = ggplot(d, aes(x=Treatment, y=Day3Mortality)) + geom_boxplot(fill="palegreen2") + geom_text(data = generate_label_df1(tHSD, 'Treatment'), aes(x = plot.labels, y = V1, label = labels)) + ylab("Day 3 Mortality Rate")

d <- read_csv("logIAPV.csv")
d = as.data.frame(d)
attr(d, "spec") <- NULL
# Note: There are 14 NC; 15 VR, VC, NR

d = d[which(d$treatment %in% c("None Chestnut", "None Rockrose", "Virus Chestnut", "Virus Rockrose")),]

d$treatment = as.factor(d$treatment)
colnames(d)[2] = "logIAPV"
levels(d$treatment) <- c('NC', 'NR', 'VC', 'VR')
colnames(d)[1] <- "Treatment"

a <- aov(logIAPV~Treatment, data=d)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)

generate_label_df2 <- function(HSD, flev){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[flev]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$logIAPV)) + 0.25)
  #boxplot.df[2,]$V1 = 0.22
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}

plotIAPV= ggplot(d, aes(x=Treatment, y=logIAPV)) + geom_boxplot(fill="paleturquoise2") + geom_text(data = generate_label_df2(tHSD, 'Treatment'), aes(x = plot.labels, y = V1, label = labels)) + ylab("Log IAPV Titer")

