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
                  high3 = raw_data$u.high3[,1], low = raw_data$u.low[,1],zero = raw_data$u.zero[,1], slope = raw_data$slope[1,], fin = raw_data$fin[1,], 
                  sess_diff = raw_data$sess.diff[1,], max = raw_data$maxi[1,], high2_m = raw_data$u.high2.m[,1], high3_m = raw_data$u.high3.m[,1], 
                  low_m = raw_data$u.low.m[,1], high_a = raw_data$u.high.a[,1], high2_a = raw_data$u.high2.a[,1],
                  high3_a = raw_data$u.high3.a[,1], low_a = raw_data$u.low.a[,1])
data = mutate(data, diff = log10(low)-log10(high2))
data = mutate(data, diff3 = log10(low)-log10(high3))
data = mutate(data, diffm = log10(low_m)-log10(high2_m))
data = mutate(data, diffm3 = log10(low_m)-log10(high3_m))
data = mutate(data, diffa = log10(low_a)-log10(high2_a))
data = mutate(data, diffa3 = log10(low_a)-log10(high3_a))

# reformat by condition
data_high = select(data, band, subj, slope, fin, max, sess_diff, diff, diff3, high)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high'] = 'opt_u'

data_high2 = select(data, band, subj, slope, fin, max, sess_diff, diff, diff3, high2)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2'] = 'opt_u'

data_high3 = select(data, band, subj, slope, fin, max, sess_diff, diff, diff3, high3)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3'] = 'opt_u'

data_low = select(data, band, subj, slope, fin, max, sess_diff, diff, diff3, low)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low'] = 'opt_u'

data_zero = select(data, band, subj, slope, fin, max, sess_diff, diff, diff3, zero)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero'] = 'opt_u'

data2 = rbind(data_low, data_high2, data_high3, data_high, data_zero)

data_high2_m = select(data, band, subj, slope, fin, max, sess_diff, diffm, diffm3, high2_m)
data_high2_m$cond = 'high2'
names(data_high2_m)[names(data_high2_m)=='high2_m'] = 'opt_u'

data_high3_m = select(data, band, subj, slope, fin, max, sess_diff, diffm, diffm3, high3_m)
data_high3_m$cond = 'high3'
names(data_high3_m)[names(data_high3_m)=='high3_m'] = 'opt_u'

data_low_m = select(data, band, subj, slope, fin, max, sess_diff, diffm, diffm3, low_m)
data_low_m$cond = 'low'
names(data_low_m)[names(data_low_m)=='low_m'] = 'opt_u'

data_m = rbind(data_low_m, data_high2_m, data_high3_m)

data_alpha = filter(data2, band == "alpha")
data_beta = filter(data2, band == "beta")
data_gamma = filter(data2, band == "low_gamma")

data_alpha_m = filter(data_m, band == "alpha")
data_beta_m = filter(data_m, band == "beta")
data_gamma_m = filter(data_m, band == "low_gamma")

plot = ggplot(data2, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u', '.png', sep = ''))

plot = ggplot(data_beta, aes(x = cond, y = log(opt_u), fill = cond ))
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binwidth = 1, binaxis='y', stackdir='center', dotsize=.5) +
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal() + geom_line(aes(group = subj), alpha = 0.6)
ggsave(paste(sens, '_opt_u_line', '.png', sep = ''))

plot = ggplot(data_m, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_m', '.png', sep = ''))

plot = ggplot(data_alpha_m, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_alpha_m', '.png', sep = ''))

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
beta_low2 = t.test(log(filter(data, band == 'beta')$high2),log(filter(data, band == 'beta')$low), paired=TRUE)
beta_low2

beta_zero2 = t.test(log(filter(data, band == 'beta')$high2),log(filter(data, band == 'beta')$zero), paired=TRUE)
beta_zero2

beta_low3 = t.test(filter(data, band == 'beta')$high3,filter(data, band == 'beta')$low, paired=TRUE)
beta_low3

beta_zero3 = t.test(log(filter(data, band == 'beta')$high3),log(filter(data, band == 'beta')$zero), paired=TRUE)
beta_zero3

# other bands
alpha_low2 = t.test(log(filter(data, band == 'alpha')$high2),log(filter(data, band == 'alpha')$low), paired=TRUE)
alpha_low2

alpha_m = t.test(log(filter(data, band == 'alpha')$high2_m),log(filter(data, band == 'alpha')$low_m), paired=TRUE)
alpha_m

gamma_low2 = t.test(log(filter(data, band == 'low_gamma')$high3),log(filter(data, band == 'low_gamma')$low), paired=TRUE)
gamma_low2

gamma_low2 = t.test(log(filter(data, band == 'low_gamma')$high3_m),log(filter(data, band == 'low_gamma')$low_m), paired=TRUE)
gamma_low2

# maintainance
beta_low2_m = t.test(log(filter(data, band == 'beta')$high2_m),log(filter(data, band == 'beta')$low_m), paired=TRUE)
beta_low2_m

beta_low3_m = t.test(log(filter(data, band == 'beta')$high3_m),log(filter(data, band == 'beta')$low_m), paired=TRUE)
beta_low3_m



# corr with slope
scatterplot = ggplot(data_beta, aes(x = diff, y = slope)) 
  scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
    theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u.png')

scatterplot = ggplot(data_beta, aes(x = diff3, y = slope)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high3)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u3.png')

scatterplot = ggplot(data_beta, aes(x = log(opt_u), y = slope, color = cond)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'opt_u', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u_high3.png')

# corr with max
scatterplot = ggplot(data_beta, aes(x = diff, y = max)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'max') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('max_opt_u.png')

scatterplot = ggplot(data_beta, aes(x = log(opt_u), y = max, color = cond)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'max') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('max_opt_u_cond.png')

# corr with final
scatterplot = ggplot(data_beta, aes(x = diff, y = fin)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'fin') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('fin_opt_u.png')

scatterplot = ggplot(data_beta, aes(x = log(opt_u), y = fin, color = cond)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(opt_u)', y = 'fin') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('fin_opt_u_cond.png')

corr_beta_high2 = cor.test(log(filter(data, band == 'beta')$high2),filter(data, band == 'beta')$slope)
corr_beta_high2

corr_beta_high3 = cor.test(log(filter(data, band == 'beta')$high3),filter(data, band == 'beta')$slope)
corr_beta_high3

corr_beta_low = cor.test(log(filter(data, band == 'beta')$low),filter(data, band == 'beta')$slope)
corr_beta_low

corr_beta_diff = cor.test(filter(data, band == 'beta')$diff,filter(data, band == 'beta')$slope, method = "pearson")
corr_beta_diff

corr_beta_diff = cor.test(filter(data, band == 'beta')$diff,filter(data, band == 'beta')$max)
corr_beta_diff

fit1 = lm(slope ~ diff, filter(data, band == 'beta'))
summary(fit1)

fit3 = lm(max ~ diff, filter(data, band == 'beta'))
summary(fit3)

fit4 = lm(fin ~ diff, filter(data, band == 'beta'))
summary(fit4)


# maintainance
scatterplot = ggplot(data_beta_m, aes(x = diffm, y = slope)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u_m.png')

scatterplot = ggplot(data_m, aes(x = log(opt_u), y = slope, color = cond)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u_m.png')

scatterplot = ggplot(data_beta_m, aes(x = log(opt_u), y = slope, color = cond)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'opt_u', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u_cond_m.png')

plot = ggplot(data_beta_m, aes(x = cond, y = log(opt_u), fill = cond ))
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  geom_dotplot(binwidth = 1, binaxis='y', stackdir='center', dotsize=.5) +
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal() + geom_line(aes(group = subj), alpha = 0.6)
ggsave(paste(sens, '_opt_u_m', '.png', sep = ''))


##############################################

# Attention Networks

#############################################

# reformat by condition
data_high_a = select(data, band, subj, slope, fin, max, diffa, diffa3, high_a)
data_high_a$cond = 'high'
names(data_high_a)[names(data_high_a)=='high_a'] = 'opt_u'

data_high2_a = select(data, band, subj, slope, fin, max, diffa, diffa3, high2_a)
data_high2_a$cond = 'high2'
names(data_high2_a)[names(data_high2_a)=='high2_a'] = 'opt_u'

data_high3_a = select(data, band, subj, slope, fin, max, diffa, diffa3, high3_a)
data_high3_a$cond = 'high3'
names(data_high3_a)[names(data_high3_a)=='high3_a'] = 'opt_u'

data_low_a = select(data, band, subj, slope, fin, max, diffa, diffa3, low_a)
data_low_a$cond = 'low'
names(data_low_a)[names(data_low_a)=='low_a'] = 'opt_u'


data_a = rbind(data_low_a, data_high2_a, data_high3_a, data_high_a)


plot = ggplot(data_a, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_a', '.png', sep = ''))

data_alpha_a = filter(data_a, band == "alpha")
data_beta_a = filter(data_a, band == "beta")
data_gamma_a = filter(data_a, band == "low_gamma")

# Stats
alpha_a = t.test(log(filter(data, band == 'alpha')$high2_a),log(filter(data, band == 'alpha')$low_a), paired=TRUE)
alpha_a

beta_a = t.test(log(filter(data, band == 'beta')$high2_a),log(filter(data, band == 'beta')$low_a), paired=TRUE)
beta_a

gamma_a = t.test(log(filter(data, band == 'low_gamma')$high2_a),log(filter(data, band == 'low_gamma')$low_a), paired=TRUE)
gamma_a


# corr with slope
scatterplot = ggplot(data_alpha_a, aes(x = diffa, y = slope)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_diff_a.png')

scatterplot = ggplot(data, aes(x = log(low_a), y = slope, color = band)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u_a.png')






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
                          high3 = raw_data$u.high3.pr[,1], low = raw_data$u.low.pr[,1],zero = raw_data$u.zero.pr[,1], slope = raw_data$slope[1,], fin = raw_data$fin[1,], 
                          max = raw_data$maxi[1,])
data_pr = mutate(data_pr, diff = log10(low)-log10(high2))
data_pr = mutate(data_pr, diff3 = log10(low)-log10(high3))

# reformat by condition
data_high_pr = select(data_pr, band, subj, slope, fin, max, diff, high)
data_high_pr$cond = 'high'
names(data_high_pr)[names(data_high_pr)=='high'] = 'opt_u'

data_high2_pr = select(data_pr, band, subj, slope, fin, max, diff, high2)
data_high2_pr$cond = 'high2'
names(data_high2_pr)[names(data_high2_pr)=='high2'] = 'opt_u'

data_high3_pr = select(data_pr, band, subj, slope, fin, max, diff, high3)
data_high3_pr$cond = 'high3'
names(data_high3_pr)[names(data_high3_pr)=='high3'] = 'opt_u'

data_low_pr = select(data_pr, band, subj, slope, fin, max, diff, low)
data_low_pr$cond = 'low'
names(data_low_pr)[names(data_low_pr)=='low'] = 'opt_u'

data_zero_pr = select(data_pr, band, subj, slope, fin, max, diff, zero)
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

## Relationship to slope

scatterplot = ggplot(data_beta_pr, aes(x = diff, y = slope)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u_pr.png')

corr_beta_high2 = cor.test(log(filter(data_pr, band == 'beta')$high2),filter(data_pr, band == 'beta')$slope)
corr_beta_high2

# combined plot
data_comb = data_frame(diff = c(data$diff, data_pr$diff), slope = c(data$slope, data_pr$slope), model = c(rep('emp', times = 60), rep('upr', times = 60)), band = c(data$band, data_pr$band) )
data_comb = filter(data_comb, band == "2")

scatterplot = ggplot(data_comb, aes(x = diff, y = slope, color = model)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise2",4))
ggsave('slope_opt_u_comb.png')
