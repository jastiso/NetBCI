# Energy
library(R.matlab)
library(ggplot2)
library(dplyr)
library(coin)
library(lmPerm)
library(car)
library(aplpack)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
setwd("/Users/stiso/Documents/R/NetBCI/")
bands = c('theta', 'alpha', 'beta', 'gamma')
nSubj = 6

# load data
E = readMat(paste('data/E.mat'))
E_corr = readMat(paste('data/E_corr.mat'))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), E = numeric(0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.high[,b_idx]
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.other[,b_idx]
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.low[,b_idx]
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
mypallette = colorRampPalette(brewer.pal(4, 'Blues'))
plot = ggplot(plot_data, aes(x = cond, y = E, col = cond, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = wes_palette("Royal1")) + 
  labs(x = 'Subgraph', y = 'Energy')  + theme_minimal()
ggsave(paste('E', '.png', sep = ''))


# scatter plot
E_plot_data = data.frame(BE = E_corr$E.exp.corr[,2], E = E_corr$E.exp.corr[,1], band = E_corr$band.order, subj = E_corr$subj.order)
E_plot_data$band = as.factor(E_plot_data$band)
E_plot_data$subj = as.factor(E_plot_data$subj)
scatterplot = ggplot(E_plot_data, aes(x = log(BE), y = E, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'ENergy') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = wes_palette("Royal1"))
ggsave(paste('dot', '_e', '.png', sep = ''))

# stats
# anova
fit1 = lm(E ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(E ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(E ~ cond, data = plot_data)
summary(stat1)
# corr
stat2 = lmp(E ~ log(BE), data = E_plot_data)
summary(stat2)
