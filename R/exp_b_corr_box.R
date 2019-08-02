# Behavior and BE correlations

library(ggplot2)
library(dplyr)
library(coin)
library(lmPerm)
library(car)
library(aplpack)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
library(coin)
library(perm)
setwd("/Users/stiso/Documents/R/NetBCI/")
bands = c( 'alpha', 'beta', 'low_gamma', 'gamma')
nSubj = 20
sens = 'grad'

# load data
max = readMat(paste('data/gc/', sens, '/max_exp_b_corr.mat', sep = ''))
min = readMat(paste('data/gc/', sens, '/min_exp_b_corr.mat', sep = ''))


data = data.frame(BE = character(length = 0), cond = character(length = 0), band = character(length = 0), corr = numeric(0), stringsAsFactors = FALSE)
for (b in c(1:length(bands))){
  #max
  data[(nrow(data)+1):(nrow(data)+nSubj),1] = rep('max', times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),2] = rep('task', times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),3] = rep(bands[b], times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),4] = max$max.corrs[,b]
  
  #min
  data[(nrow(data)+1):(nrow(data)+nSubj),1] = rep('min', times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),2] = rep('task', times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),3] = rep(bands[b], times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),4] = min$min.corrs[,b]
  
 
}

plot = ggplot(data, aes(x = BE, y = corr, fill = band))
plot + geom_boxplot(notch = FALSE, lwd = 1) +
  scale_fill_manual(values = brewer.pal(4,'Greys')) + 
  labs(x = 'Subgraph', y = 'Corr')  + theme_minimal()
ggsave(paste('b_exp_corr', '.png', sep = ''))

# stats
# anova
fit1 = lmp(corr ~ BE + band, data = data)
summary(fit1)
Anova(fit1)
fit2 = lm(corr ~ band+BE, data = data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(corr ~ BE, data = data)
summary(stat1)
# for each condition
bl = filter(data, cond == "bl")
task = filter(data, cond == 'task')
stat = lmp(corr ~ BE + band, data = bl)
summary(stat)
stat = lmp(corr ~ BE + band, data = task)
Anova(stat)
# specififc ttest
stat = permTS(filter(data, BE == 'min')$corr, filter(data, BE == 'min_bl')$corr)
stat
stat = permTS(filter(data, BE == 'max')$corr, filter(data, BE == 'max_bl')$corr)
stat

