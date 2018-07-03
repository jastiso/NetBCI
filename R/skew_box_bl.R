# Skew
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
S = readMat(paste('data/S_bl.mat'))
S_corr = readMat(paste('data/S_corr_bl.mat'))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), S = numeric(0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.high[,b_idx]
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.other[,b_idx]
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.low[,b_idx]
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
mypallette = colorRampPalette(brewer.pal(4, 'Blues'))
plot = ggplot(plot_data, aes(x = cond, y = S, col = cond, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = wes_palette("Moonrise1")) + 
  labs(x = 'Subgraph', y = 'Energy')  + theme_minimal()
ggsave(paste('S', '_bl.png', sep = ''))


# scatter plot
S_plot_data = data.frame(BE = S_corr$S.exp.corr[,2], S = S_corr$S.exp.corr[,1], band = S_corr$band.order, subj = S_corr$subj.order)
S_plot_data$band = as.factor(S_plot_data$band)
S_plot_data$subj = as.factor(S_plot_data$subj)
scatterplot = ggplot(S_plot_data, aes(x = log(BE), y = S, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'ENergy') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = wes_palette("Moonrise1"))
ggsave(paste('dot', '_s', '_bl.png', sep = ''))

# stats
# anova
fit1 = lm(S ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(S ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(S ~ cond, data = plot_data)
summary(stat1)
# corr
stat2 = lmp(S ~ log(BE), data = S_plot_data)
summary(stat2)
