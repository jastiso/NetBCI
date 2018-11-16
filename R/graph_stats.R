# Look at some graph metrics

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

raw_data = readMat(paste('data/wpli/', sens, '/graph_stats.mat', sep = ''))
data = data.frame(band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), high = raw_data$high.eff[,1], high2 = raw_data$high2.eff[,1],
                  high3 = raw_data$high3.eff[,1], low = raw_data$low.eff[,1], zero = raw_data$zero.eff[,1], high_rm = raw_data$high.eff.rm[,1], high2_rm = raw_data$high2.eff.rm[,1],
                  high3_rm = raw_data$high3.eff.rm[,1], low_rm = raw_data$low.eff[,1], zero_rm = raw_data$zero.eff.rm[,1])

# reformat by condition
data_high = select(data, band, subj, high)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high'] = 'eff'

data_high2 = select(data, band, subj, high2)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2'] = 'eff'

data_high3 = select(data, band, subj, high3)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3'] = 'eff'

data_low = select(data, band, subj, low)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low'] = 'eff'

data_zero = select(data, band, subj, zero)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero'] = 'eff'

data2 = rbind(data_high, data_low, data_high2, data_high3, data_zero)

data_alpha = filter(data2, band == "alpha")
data_beta = filter(data2, band == "beta")
data_gamma = filter(data2, band == "low_gamma")

plot = ggplot(data2, aes(x = cond, y = eff, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Local Eff')  + theme_minimal()
ggsave(paste(sens, '_eff', '.png', sep = ''))


# stats
beta_low2 = t.test(filter(data, band == 'beta')$high2,filter(data, band == 'beta')$low, paired=TRUE)
beta_low2

beta_zero2 = t.test(filter(data, band == 'beta')$high2,filter(data, band == 'beta')$zero, paired=TRUE)
beta_zero2

beta_low3 = t.test(filter(data, band == 'beta')$high3,filter(data, band == 'beta')$low, paired=TRUE)
beta_low3

beta_zero3 = t.test(filter(data, band == 'beta')$high3,filter(data, band == 'beta')$zero, paired=TRUE)
beta_zero3

alpha_low = t.test(filter(data, band == 'alpha')$high,filter(data, band == 'alpha')$low, paired=TRUE)
alpha_low

alpha_zero = t.test(filter(data, band == 'alpha')$high,filter(data, band == 'alpha')$zero, paired=TRUE)
alpha_zero
