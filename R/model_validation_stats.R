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
  labs(x = 'Temporal Predictor', y = 'Correlation (Pearsons R)')  + theme_minimal()
ggsave(paste(sens, '_model_validation', '.png', sep = ''))

