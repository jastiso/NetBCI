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
                  high3_rm = raw_data$high3.eff.rm[,1], low_rm = raw_data$low.eff[,1], zero_rm = raw_data$zero.eff.rm[,1], high_nc = raw_data$high.nc[,1], high2_nc = raw_data$high2.nc[,1],
                  high3_nc = raw_data$high3.nc[,1], low_nc = raw_data$low.nc[,1], zero_nc = raw_data$zero.nc[,1], high_pc = raw_data$high.pc[,1], high2_pc = raw_data$high2.pc[,1],
                  high3_pc = raw_data$high3.pc[,1], low_pc = raw_data$low.pc[,1], zero_pc = raw_data$zero.pc[,1])

####################

# Efficiency

####################


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






####################

# Number of communities

####################


# reformat by condition
data_high = select(data, band, subj, high_nc)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high_nc'] = 'num_comm'

data_high2 = select(data, band, subj, high2_nc)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2_nc'] = 'num_comm'

data_high3 = select(data, band, subj, high3_nc)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3_nc'] = 'num_comm'

data_low = select(data, band, subj, low_nc)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low_nc'] = 'num_comm'

data_zero = select(data, band, subj, zero_nc)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero_nc'] = 'num_comm'

data2 = rbind(data_high, data_low, data_high2, data_high3, data_zero)

data_alpha = filter(data2, band == "alpha")
data_beta = filter(data2, band == "beta")
data_gamma = filter(data2, band == "low_gamma")

plot = ggplot(data2, aes(x = cond, y = num_comm, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Number of Communities')  + theme_minimal()
ggsave(paste(sens, '_num_comm', '.png', sep = ''))


# stats
beta_low2_nc = t.test(filter(data, band == 'beta')$high2_nc,filter(data, band == 'beta')$low_nc, paired=TRUE)
beta_low2_nc

beta_zero2_nc = t.test(filter(data, band == 'beta')$high2_nc,filter(data, band == 'beta')$zero_nc, paired=TRUE)
beta_zero2_nc

beta_low3_nc = t.test(filter(data, band == 'beta')$high3_nc,filter(data, band == 'beta')$low_nc, paired=TRUE)
beta_low3_nc

beta_zero3_nc = t.test(filter(data, band == 'beta')$high3_nc,filter(data, band == 'beta')$zero_nc, paired=TRUE)
beta_zero3_nc

alpha_low_nc = t.test(filter(data, band == 'alpha')$high_nc,filter(data, band == 'alpha')$low_nc, paired=TRUE)
alpha_low_nc

alpha_zero_nc = t.test(filter(data, band == 'alpha')$high_nc,filter(data, band == 'alpha')$zero_nc, paired=TRUE)
alpha_zero_nc



####################

# Participation Coefficient

####################


# reformat by condition
data_high = select(data, band, subj, high_pc)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high_pc'] = 'pc'

data_high2 = select(data, band, subj, high2_pc)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2_pc'] = 'pc'

data_high3 = select(data, band, subj, high3_pc)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3_pc'] = 'pc'

data_low = select(data, band, subj, low_pc)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low_pc'] = 'pc'

data_zero = select(data, band, subj, zero_pc)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero_pc'] = 'pc'

data2 = rbind(data_high, data_low, data_high2, data_high3, data_zero)

data_alpha = filter(data2, band == "alpha")
data_beta = filter(data2, band == "beta")
data_gamma = filter(data2, band == "low_gamma")


plot = ggplot(data2, aes(x = cond, y = pc, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Participation Coefficient')  + theme_minimal()
ggsave(paste(sens, '_pc', '.png', sep = ''))


# stats
beta_low2_pc = t.test(filter(data, band == 'beta')$high2_pc,filter(data, band == 'beta')$low_pc, paired=TRUE)
beta_low2_pc

beta_zero2_pc = t.test(filter(data, band == 'beta')$high2_pc,filter(data, band == 'beta')$zero_pc, paired=TRUE)
beta_zero2_pc

beta_low3_pc = t.test(filter(data, band == 'beta')$high3_pc,filter(data, band == 'beta')$low_pc, paired=TRUE)
beta_low3_pc

beta_zero3_pc = t.test(filter(data, band == 'beta')$high3_pc,filter(data, band == 'beta')$zero_pc, paired=TRUE)
beta_zero3_pc

alpha_low_pc = t.test(filter(data, band == 'alpha')$high_pc,filter(data, band == 'alpha')$low_pc, paired=TRUE)
alpha_low_pc

alpha_zero_pc = t.test(filter(data, band == 'alpha')$high_pc,filter(data, band == 'alpha')$zero_pc, paired=TRUE)
alpha_zero_pc
