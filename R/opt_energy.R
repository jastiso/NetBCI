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

raw_data = readMat(paste('data/wpli/', sens, '/opt_energy.mat', sep = ''))
data = data.frame(band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), high = raw_data$u.high[,1], high2 = raw_data$u.high2[,1],
                  high3 = raw_data$u.high3[,1], low = raw_data$u.low[,1],zero = raw_data$u.zero[,1], slope = raw_data$slope[1,], imp = raw_data$imp[1,], 
                  sess_diff = raw_data$sess.diff[1,])
data = mutate(data, diff = log10(low)-log10(high2))

# reformat by condition
data_high = select(data, band, subj, slope, diff, high)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high'] = 'opt_u'

data_high2 = select(data, band, subj,slope, diff, high2)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2'] = 'opt_u'

data_high3 = select(data, band, subj,slope,diff, high3)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3'] = 'opt_u'

data_low = select(data, band, subj, slope,diff, low)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low'] = 'opt_u'

data_zero = select(data, band, subj, slope,diff, zero)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero'] = 'opt_u'

data2 = rbind(data_low, data_high2)

data_alpha = filter(data2, band == "alpha")
data_beta = filter(data2, band == "beta")
data_gamma = filter(data2, band == "low_gamma")

plot = ggplot(data2, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u', '.png', sep = ''))

plot = ggplot(data_alpha, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_alpha', '.png', sep = ''))

plot = ggplot(data_beta, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_beta', '.png', sep = ''))

plot = ggplot(data_gamma, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_gamma', '.png', sep = ''))


## Stats
beta_low2 = t.test(filter(data, band == 'beta')$high2,filter(data, band == 'beta')$low, paired=TRUE)
beta_low2

beta_zero2 = t.test(filter(data, band == 'beta')$high2,filter(data, band == 'beta')$zero, paired=TRUE)
beta_zero2

beta_low3 = t.test(filter(data, band == 'beta')$high3,filter(data, band == 'beta')$low, paired=TRUE)
beta_low3

beta_zero3 = t.test(filter(data, band == 'beta')$high3,filter(data, band == 'beta')$zero, paired=TRUE)
beta_zero3



# corr with slope
scatterplot = ggplot(data_beta, aes(x = diff, y = slope)) 
  scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
    theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u.png')

corr_beta_high2 = cor.test(log(filter(data, band == 'beta')$high2),filter(data, band == 'beta')$slope)
corr_beta_high2

corr_beta_low = cor.test(log(filter(data, band == 'beta')$low),filter(data, band == 'beta')$slope)
corr_beta_low

corr_beta_diff = cor.test(filter(data, band == 'beta')$diff,filter(data, band == 'beta')$slope, method = "spearman")
corr_beta_diff

fit1 = lm(slope ~ diff, filter(data, band == 'beta'))
summary(fit1)

################################################

# Different Control Set

################################################



raw_data = readMat(paste('data/wpli/', sens, '/opt_energy_control.mat', sep = ''))
data_control = data.frame(band = unlist(raw_data$band.order.c), subj = unlist(raw_data$subj.order.c), high = raw_data$u.high.c[,1], high2 = raw_data$u.high2.c[,1],
                  high3 = raw_data$u.high3.c[,1], low = raw_data$u.low.c[,1],zero = raw_data$u.zero.c[,1])

# reformat by condition
data_high_control = select(data_control, band, subj, high)
data_high_control$cond = 'high'
names(data_high_control)[names(data_high_control)=='high'] = 'opt_u'

data_high2_control = select(data_control, band, subj, high2)
data_high2_control$cond = 'high2'
names(data_high2_control)[names(data_high2_control)=='high2'] = 'opt_u'

data_high3_control = select(data_control, band, subj, high3)
data_high3_control$cond = 'high3'
names(data_high3_control)[names(data_high3_control)=='high3'] = 'opt_u'

data_low_control = select(data_control, band, subj, low)
data_low_control$cond = 'low'
names(data_low_control)[names(data_low_control)=='low'] = 'opt_u'

data_zero_control = select(data_control, band, subj, zero)
data_zero_control$cond = 'zero'
names(data_zero_control)[names(data_zero_control)=='zero'] = 'opt_u'

data2_control = rbind(data_high_control, data_low_control, data_high2_control, data_high3_control, data_zero_control)

data_alpha_control = filter(data2_control, band == "alpha")
data_beta_control = filter(data2_control, band == "beta")
data_gamma_control = filter(data2_control, band == "low_gamma")

plot = ggplot(data2_control, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_control', '.png', sep = ''))

plot = ggplot(data_alpha_control, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_alpha_control', '.png', sep = ''))

plot = ggplot(data_beta_control, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_beta_control', '.png', sep = ''))

plot = ggplot(data_gamma_control, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_gamma_control', '.png', sep = ''))


## Stats
beta_low2_control = t.test(filter(data_control, band == 'beta')$high2,filter(data_control, band == 'beta')$low, paired=TRUE)
beta_low2_control

beta_zero2_control = t.test(filter(data_control, band == 'beta')$high2,filter(data_control, band == 'beta')$zero, paired=TRUE)
beta_zero2_control

beta_low3_control = t.test(filter(data_control, band == 'beta')$high3,filter(data_control, band == 'beta')$low, paired=TRUE)
beta_low3_control

beta_zero3_control = t.test(filter(data_control, band == 'beta')$high3,filter(data_control, band == 'beta')$zero, paired=TRUE)
beta_zero3_control



################################################

# UPR

################################################



raw_data = readMat(paste('data/wpli/', sens, '/opt_energy_pr.mat', sep = ''))
data_pr = data.frame(band = unlist(raw_data$band.order.pr), subj = unlist(raw_data$subj.order.pr), high = raw_data$u.high.pr[,1], high2 = raw_data$u.high2.pr[,1],
                          high3 = raw_data$u.high3.pr[,1], low = raw_data$u.low.pr[,1],zero = raw_data$u.zero.pr[,1])

# reformat by condition
data_high_pr = select(data_pr, band, subj, high)
data_high_pr$cond = 'high'
names(data_high_pr)[names(data_high_pr)=='high'] = 'opt_u'

data_high2_pr = select(data_pr, band, subj, high2)
data_high2_pr$cond = 'high2'
names(data_high2_pr)[names(data_high2_pr)=='high2'] = 'opt_u'

data_high3_pr = select(data_pr, band, subj, high3)
data_high3_pr$cond = 'high3'
names(data_high3_pr)[names(data_high3_pr)=='high3'] = 'opt_u'

data_low_pr = select(data_pr, band, subj, low)
data_low_pr$cond = 'low'
names(data_low_pr)[names(data_low_pr)=='low'] = 'opt_u'

data_zero_pr = select(data_pr, band, subj, zero)
data_zero_pr$cond = 'zero'
names(data_zero_pr)[names(data_zero_pr)=='zero'] = 'opt_u'

data2_pr = rbind(data_high_pr, data_low_pr, data_high2_pr, data_high3_pr, data_zero_pr)

data_alpha_pr = filter(data2_pr, band == "alpha")
data_beta_pr = filter(data2_pr, band == "beta")
data_gamma_pr = filter(data2_pr, band == "low_gamma")

plot = ggplot(data2_pr, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_pr', '.png', sep = ''))

plot = ggplot(data_alpha_pr, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_alpha_pr', '.png', sep = ''))

plot = ggplot(data_beta_pr, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_beta_pr', '.png', sep = ''))

plot = ggplot(data_gamma_pr, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_gamma_pr', '.png', sep = ''))


## Stats
beta_low2_pr = t.test(filter(data_pr, band == 'beta')$high2,filter(data_pr, band == 'beta')$low, paired=TRUE)
beta_low2_pr

beta_zero2_pr = t.test(filter(data_pr, band == 'beta')$high2,filter(data_pr, band == 'beta')$zero, paired=TRUE)
beta_zero2_pr

beta_low3_pr = t.test(filter(data_pr, band == 'beta')$high3,filter(data_pr, band == 'beta')$low, paired=TRUE)
beta_low3_pr

beta_zero3_pr = t.test(filter(data_pr, band == 'beta')$high3,filter(data_pr, band == 'beta')$zero, paired=TRUE)
beta_zero3_pr

