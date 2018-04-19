library(agricolae)
data(sweetpotato)
# Analysis of variance model
model<-aov(yield~virus, data=sweetpotato)
out <- HSD.test(model,"virus", group=TRUE,console=TRUE,
                main="Yield of sweetpotato\nDealt with different virus")
#stargraph
# Variation range: max and min
plot(out)
#endgraph
out<-HSD.test(model,"virus", group=FALSE)
print(out$comparison)
# Old version HSD.test()
df<-df.residual(model)
MSerror<-deviance(model)/df
with(sweetpotato,HSD.test(yield,virus,df,MSerror, group=TRUE,console=TRUE, main="Yield of sweetpotato. Dealt with different virus"))






library(plyr)
library(ggplot2)
library(multcompView)

set.seed(0)
lev <- gl(3, 10)
y <- c(rnorm(10), rnorm(10) + 0.1, rnorm(10) + 3)
d <- data.frame(lev=lev, y=y)

a <- aov(y~lev, data=d)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)

p_base <- ggplot(d, aes(x=lev, y=y)) + geom_boxplot() + geom_text(data = generate_label_df(tHSD, 'lev'), aes(x = plot.labels, y = V1, label = labels))

generate_label_df <- function(HSD, flev){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[flev]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$y)) + 0.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}
