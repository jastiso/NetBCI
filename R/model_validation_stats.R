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
library(ez)
library(plyr)
setwd("/Users/stiso/Documents/R/NetBCI/")
bands = c('alpha', 'beta', 'low_gamma')
nSubj = 20
sens = 'grad'

raw_data = readMat(paste('data/wpli/model_validation.mat', sep = ''))
data = data.frame(corr = c(raw_data$empirical[1,], raw_data$null[1,]), model = c(rep('emp', times = nSubj), rep('null', times = nSubj)))

plot = ggplot(data, aes(x = model, y = corr, fill = model) )
plot + geom_violin(trim = FALSE, position = position_dodge(0.75)) + 
  geom_boxplot(width = 0.15,
               position = position_dodge(0.75)
  ) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Temporal Predictor', y = 'Correlation (Pearsons R)')  + theme_minimal() + scale_fill_manual(values = c(rgb(81/255,184/255,161/255), 'grey'))
ggsave(paste(sens, '_model_validation', '.pdf', sep = ''))

stat = permTS(raw_data$empirical[1,], raw_data$null[1,])
stat




################
# only 3rd SG
#################
raw_data = readMat(paste('data/wpli/grad/model_validation1.mat', sep = ''))
data1 = data.frame(corr = c(raw_data$high[,1], raw_data$low[,1],raw_data$high.null[,1], raw_data$low.null[,1]), 
                  model = c(rep('emp', times = nSubj*2), rep('null', times = nSubj*2)),
                  sg = c(rep('high', times = nSubj), rep('low', times = nSubj), rep('high', times = nSubj), rep('low', times = nSubj)))

plot = ggplot(data1, aes(x = sg, y = corr, fill = model) )
plot + geom_violin(trim = TRUE, position = position_dodge(0.75)) + 
  geom_boxplot(width = 0.15,
               position = position_dodge(0.75)
  ) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, position = position_dodge(0.75)) + 
  labs(x = 'Temporal Predictor', y = 'Correlation (Pearsons R)')  + theme_minimal() + scale_fill_manual(values = c(rgb(81/255,184/255,161/255), 'grey'))
ggsave(paste(sens, '_model_validation1', '.pdf', sep = ''))

stat = t.test(raw_data$high[,1], raw_data$high.null[,1], paired=TRUE)
stat

stat = t.test(raw_data$low[,1], raw_data$low.null[,1], paired=TRUE)
stat
