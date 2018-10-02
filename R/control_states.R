# Regional expression of control states

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
bands = c('alpha', 'beta', 'low_gamma', 'gamma')
nSubj = 20
sens = 'grad'

raw_data = readMat(paste('data/wpli/', sens, '/control_state.mat', sep = ''))
data = data.frame(band = unlist(raw_data$band.ord), region = unlist(raw_data$region.ord), high = raw_data$h.region[,1], low = raw_data$l.region[1,])

plot = ggplot(data, aes(x = region, y = high, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_high_control_state', '.png', sep = ''))

plot = ggplot(data, aes(x = region, y = low, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_low_control_state', '.png', sep = ''))
