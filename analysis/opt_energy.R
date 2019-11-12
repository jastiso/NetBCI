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
dens = readMat(paste('data/wpli/', sens, '/density.mat', sep = ''))
data = data.frame(band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), high = raw_data$u.high[,1], high2 = raw_data$u.high2[,1],
                  high3 = raw_data$u.high3[,1], low = raw_data$u.low[,1], slope = raw_data$slope[1,], 
                   high_a = raw_data$u.high.a[,1], high2_a = raw_data$u.high2.a[,1], high3_a = raw_data$u.high3.a[,1], low_a = raw_data$u.low.a[,1])
data = mutate(data, diff3 = log10(low)-log10(high3))
data = mutate(data, diffa3 = log10(low_a)-log10(high3_a))
data = mutate(data, diffa23 = log10(high2_a)-log10(high3_a))
data = mutate(data, diffa13 = log10(high_a)-log10(high3_a))
data = mutate(data, diffa2l = log10(low_a)-log10(high2_a))
data = mutate(data, diffa1l = log10(low_a)-log10(high_a))

# add density 
dens_data = data.frame(band = unlist(dens$band.order), subj = unlist(dens$subj.order), high3 = dens$d.high3[,1], low = dens$d.low[,1])
dens_data = mutate(data, dens_diff = low-high3)
# check order is correct - want both to be false
any(dens_data$band != data$band)
any(dens_data$subj != data$subj)

data$dens_diff = dens_data$dens_diff

# reformat by condition
data_high = select(data, band, subj, slope, diff3, high)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high'] = 'opt_u'

data_high2 = select(data, band, subj, slope, diff3, high2)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2'] = 'opt_u'

data_high3 = select(data, band, subj, slope, diff3, high3)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3'] = 'opt_u'

data_low = select(data, band, subj, slope, diff3, low)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low'] = 'opt_u'

data2 = rbind(data_low, data_high2, data_high3, data_high)


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
beta_low3 = t.test(log(filter(data, band == 'beta')$high3),log(filter(data, band == 'beta')$low), paired=TRUE)
beta_low3

alpha_low3 = t.test(log(filter(data, band == 'alpha')$high3),log(filter(data, band == 'alpha')$low), paired=TRUE)
alpha_low3

gamma_low3 = t.test(log(filter(data, band == 'low_gamma')$high3),log(filter(data, band == 'low_gamma')$low), paired=TRUE)
gamma_low3




# corr with slope
scatterplot = ggplot(data_beta, aes(x = diff3, y = slope)) 
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

# corr stats

corr_beta_diff = cor.test(filter(data, band == 'beta')$diff3,filter(data, band == 'beta')$slope)
corr_beta_diff

corr_alpha_diff = cor.test(filter(data, band == 'alpha')$diff3,filter(data, band == 'alpha')$slope)
corr_alpha_diff

corr_alpha_diff = cor.test(filter(data, band == 'low_gamma')$diff3,filter(data, band == 'low_gamma')$slope)
corr_alpha_diff



##############################################

# Attention Networks

#############################################

# reformat by condition
data_high_a = select(data, band, subj, slope, diffa23, diffa13, diffa3, high_a)
data_high_a$cond = 'high'
names(data_high_a)[names(data_high_a)=='high_a'] = 'opt_u'

data_high2_a = select(data, band, subj, slope, diffa23, diffa13, diffa3, high2_a)
data_high2_a$cond = 'high2'
names(data_high2_a)[names(data_high2_a)=='high2_a'] = 'opt_u'

data_high3_a = select(data, band, subj, slope, diffa23, diffa13, diffa3, high3_a)
data_high3_a$cond = 'high3'
names(data_high3_a)[names(data_high3_a)=='high3_a'] = 'opt_u'

data_low_a = select(data, band, subj, slope, diffa23, diffa13, diffa3, low_a)
data_low_a$cond = 'low'
names(data_low_a)[names(data_low_a)=='low_a'] = 'opt_u'

data_a = rbind(data_high3_a, data_low_a)


plot = ggplot(data_a, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_violin(aes(fill = band), trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(aes(fill = band), width = 0.15,
               position = position_dodge(0.9)
  ) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Optimal Energy')  + theme_minimal()
ggsave(paste(sens, '_opt_u_a', '.png', sep = ''))

data_alpha_a = filter(data_a, band == "alpha")
data_beta_a = filter(data_a, band == "beta")
data_gamma_a = filter(data_a, band == "low_gamma")

# Stats
alpha_a = t.test(log(filter(data, band == 'alpha')$high3_a),log(filter(data, band == 'alpha')$low_a), paired=TRUE)
alpha_a

beta_a = t.test(log(filter(data, band == 'beta')$high3_a),log(filter(data, band == 'beta')$low_a), paired=TRUE)
beta_a

gamma_a = t.test(log(filter(data, band == 'low_gamma')$high3_a),log(filter(data, band == 'low_gamma')$low_a), paired=TRUE)
gamma_a


# corr with slope
scatterplot = ggplot(data_beta_a, aes(x = diffa3, y = slope)) 
scatterplot+ geom_smooth(method="lm", color = rgb(81/255,184/255,161/255)) + geom_point(size = 6, color = rgb(81/255,184/255,161/255)) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() + scale_color_manual(values =  rgb(81/255,184/255,161/255))
ggsave('slope_diff_a.png')

scatterplot = ggplot(data_beta_a, aes(x = diff2a3, y = slope)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_diff_a2.png')

scatterplot = ggplot(data_beta_a, aes(x = diffa, y = slope)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_diff23_a.png')

scatterplot = ggplot(data, aes(x = (diffa3), y = slope, color = band)) 
scatterplot+ geom_smooth(method="lm", size = 4) + geom_point(size = 4) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('slope_opt_u_a.png')

scatterplot = ggplot(data_beta_a, aes(x = log(opt_u), y = slope, color = cond)) 
scatterplot+ geom_smooth(method="lm", size = 3) + geom_point(size = 4) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() + #scale_color_manual(values =  wes_palette("Moonrise1",4))
 #scale_color_manual(values = c(rgb(53/255,57/255,61/255), rgb(169/255,225/255,186/255)))
  #scale_color_manual(values = c(rgb(0,169/255,208/255), rgb(215/255,190/255,123/255)))
  scale_color_manual(values = c(rgb(215/255,190/255,123/255), rgb(33/255,67/255,104/255)))
ggsave('slope_opt_u_interaction.svg')

# individual SGs
corr_beta_high3 = cor.test(log(filter(data, band == 'beta')$high3_a),filter(data, band == 'beta')$slope)
corr_beta_high3

corr_beta_low = cor.test(log(filter(data, band == 'beta')$low_a),filter(data, band == 'beta')$slope)
corr_beta_low

#difference
corr_beta_diff = cor.test(filter(data, band == 'beta')$diffa3,filter(data, band == 'beta')$slope, method = "pearson")
corr_beta_diff

# other high SG
corr_beta_diff23 = cor.test(filter(data, band == 'beta')$diffa23,filter(data, band == 'beta')$slope)
corr_beta_diff23

corr_beta_diff13 = cor.test(filter(data, band == 'beta')$diffa13,filter(data, band == 'beta')$slope)
corr_beta_diff13

corr_beta_diff2l = cor.test(filter(data, band == 'beta')$diffa2l,filter(data, band == 'beta')$slope)
corr_beta_diff2l

corr_beta_diff1l = cor.test(filter(data, band == 'beta')$diffa1l,filter(data, band == 'beta')$slope)
corr_beta_diff1l

# controlling for density
fit1 = lm(slope ~ diffa3 + dens_diff, filter(data, band == 'beta'))
summary(fit1)



################################################

# Different Control Set

################################################



raw_data = readMat(paste('data/wpli/', sens, '/opt_energy_control.mat', sep = ''))
data_control = data.frame(slope = raw_data$slope[1,], band = unlist(raw_data$band.order.c), subj = unlist(raw_data$subj.order.c),
                  high3 = raw_data$u.high3.c[,1], low = raw_data$u.low.c[,1])
data_control = mutate(data_control, diff3 = log10(low)-log10(high3))


# reformat by condition
data_diff_control = select(data_control, band, subj, diff3, slope)
data_diff_control$cond = 'diff'
names(data_diff_control)[names(data_diff_control)=='diff3'] = 'opt_u'

data_high3_control = select(data_control, band, subj, high3, slope)
data_high3_control$cond = 'high3'
names(data_high3_control)[names(data_high3_control)=='high3'] = 'opt_u'

data_low_control = select(data_control, band, subj, low, slope)
data_low_control$cond = 'low'
names(data_low_control)[names(data_low_control)=='low'] = 'opt_u'

data2_control = rbind(data_low_control, data_high3_control, data_diff_control)

data_alpha_control = filter(data2_control, band == "alpha")
data_beta_control = filter(data2_control, band == "beta")
data_gamma_control = filter(data2_control, band == "low_gamma")

plot = ggplot(data2_control, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_violin(aes(fill = band), trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(aes(fill = band), width = 0.15,
               position = position_dodge(0.9)
  ) +
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


## Relationship to slope

scatterplot = ggplot(data_diff_control, aes(x = opt_u, y = slope, color = band)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('slope_opt_u_cont.png')

corr_beta_high3 = cor.test(filter(data_control, band == 'beta')$diff3,filter(data_control, band == 'beta')$slope)
corr_beta_high3

################################################

# UPR

################################################

# here, we are only loading data for attentional target states

raw_data = readMat(paste('data/wpli/', sens, '/opt_energy_pr.mat', sep = ''))
data_pr = data.frame(band = unlist(raw_data$band.order.pr), subj = unlist(raw_data$subj.order.pr),
                          high3 = raw_data$u.high3.pr[,1], low = raw_data$u.low.pr[,1], slope = raw_data$slope[1,])
data_pr = mutate(data_pr, diff3 = log10(low)-log10(high3))

# reformat by condition
data_high3_pr = select(data_pr, band, subj, slope, diff3, high3)
data_high3_pr$cond = 'high3'
names(data_high3_pr)[names(data_high3_pr)=='high3'] = 'opt_u'

data_low_pr = select(data_pr, band, subj, slope, diff3, low)
data_low_pr$cond = 'low'
names(data_low_pr)[names(data_low_pr)=='low'] = 'opt_u'

data2_pr = rbind(data_low_pr, data_high3_pr)

data_alpha_pr = filter(data2_pr, band == "alpha")
data_beta_pr = filter(data2_pr, band == "beta")
data_gamma_pr = filter(data2_pr, band == "low_gamma")

plot = ggplot(data2_pr, aes(x = cond, y = log(opt_u), fill = band) )
plot + geom_violin(aes(fill = band), trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(aes(fill = band), width = 0.15,
               position = position_dodge(0.9)
  ) + 
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


## Relationship to slope

scatterplot = ggplot(data_beta_pr, aes(x = diff3, y = slope)) 
scatterplot+ geom_smooth(method="lm", color = "black") + geom_point(size = 6, color ="grey") + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('slope_opt_u_pr.png')

corr_beta_high3 = cor.test((filter(data_pr, band == 'beta')$diff3),filter(data_pr, band == 'beta')$slope)
corr_beta_high3

# combined plot
data_comb = data_frame(diff = c(data$diffa3, data_pr$diff3, data_control$diff3, data$diff3), slope = c(data$slope, data_pr$slope, data_control$slope, data$slope), model = c(rep('emp', times = 60), rep('upr', times = 60), rep('cont', times = 60), rep('state', times = 60)), band = c(data$band, data_pr$band, data_control$band, data$band) )
data_comb = filter(data_comb, band == "2")

plot_data = filter(data_comb, model != "cont", model != "state")
plot_data$model = as.factor(plot_data$model)
plot_data$model <- relevel(plot_data$model, "upr")
scatterplot = ggplot(plot_data, aes(x = diff, y = slope, color = model)) 
scatterplot+ geom_smooth(method="lm", size=3) + geom_point(size = 4) + 
  labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() + #scale_color_manual(values = c(rgb(215/255,190/255,123/255), rgb(33/255,67/255,104/255)))
  scale_color_manual(values = c('grey', rgb(81/255,184/255,161/255)))
  #scale_color_manual(values =  rev(wes_palette("GrandBudapest1",4)))
ggsave('slope_opt_u_comb.svg')



########################

# State controls

########################

raw_data = readMat(paste('data/wpli/', sens, '/opt_energy_state.mat', sep = ''))
data_state = data.frame(slope = raw_data$slope[1,], band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), high3_mag = raw_data$u.high3.mag[,1], 
                        low_mag = raw_data$u.low.mag[,1], high3_2 = raw_data$u.high3.2[,1], low_2 = raw_data$u.low.2[,1],high3_cent = raw_data$u.high3.centered[,1], low_cent = raw_data$u.low.centered[,1],
                        high3_012 = raw_data$u.high3.012[,1], low_012 = raw_data$u.low.012[,1], high3_inv = raw_data$u.high3.inv[,1], low_inv = raw_data$u.low.inv[,1])
data_state = mutate(data_state, diff_mag = log10(low_mag)-log10(high3_mag))
data_state = mutate(data_state, diff_2 = log10(low_2)-log10(high3_2))
data_state = mutate(data_state, diff_cent = log10(low_cent)-log10(high3_cent))
data_state = mutate(data_state, diff_012 = log10(low_012)-log10(high3_012))
data_state = mutate(data_state, diff_inv = log10(low_inv)-log10(high3_inv))

## Relationship to slope

scatterplot = ggplot(filter(data_state, band == "beta"), aes(x = diff_cent, y = slope)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('slope_opt_state_control.png')

corr_beta_mag = cor.test(filter(data_state, band == 'beta')$diff_mag,filter(data_state, band == 'beta')$slope)
corr_beta_mag

corr_beta_2 = cor.test(filter(data_state, band == 'beta')$diff_2,filter(data_state, band == 'beta')$slope)
corr_beta_2

corr_beta_cent = cor.test(filter(data_state, band == 'beta')$diff_cent,filter(data_state, band == 'beta')$slope)
corr_beta_cent

corr_beta_012 = cor.test(filter(data_state, band == 'beta')$diff_012,filter(data_state, band == 'beta')$slope)
corr_beta_012

corr_beta_inv = cor.test(filter(data_state, band == 'beta')$diff_inv,filter(data_state, band == 'beta')$slope)
corr_beta_inv


########################

# parameters

########################

raw_data = readMat(paste('data/wpli/', sens, '/opt_energy_params.mat', sep = ''))
data_params = data.frame(slope = raw_data$slope[1,], band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), high3_1 = raw_data$u.high3.1[,1], 
                        low_1 = raw_data$u.low.1[,1], high3_2 = raw_data$u.high3.2[,1], low_2 = raw_data$u.low.2[,1])
data_params = mutate(data_params, diff_1 = log10(low_1)-log10(high3_1))
data_params = mutate(data_params, diff_2 = log10(low_2)-log10(high3_2))

## Relationship to slope

scatterplot = ggplot(data_params, aes(x = diff_1, y = slope, color = band)) 
scatterplot+ geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'log(low)-log(high)', y = 'slope') +
  theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('slope_opt_state_params.png')

corr_beta_1 = cor.test(filter(data_params, band == 'beta')$diff_1,filter(data_params, band == 'beta')$slope)
corr_beta_1

corr_beta_2 = cor.test(filter(data_params, band == 'beta')$diff_2,filter(data_params, band == 'beta')$slope)
corr_beta_2


########################

# magnitude

########################

raw_data = readMat(paste('data/wpli/', sens, '/opt_energy_mag.mat', sep = ''))
data_mag = data.frame(slope = raw_data$slope[1,], subj = unlist(raw_data$subj.order), sim = raw_data$sim.order[1,],
                         high3 = raw_data$u.high3[,1], low = raw_data$u.low[,1])
data_mag = mutate(data_mag, diff = log10(low)-log10(high3))

nSim = 250

## Relationship to slope
stat = list()
for (i in 1:nSim){
  curr = dplyr::filter(data_mag, sim == i)
  pval = cor.test(curr$diff, curr$slope)
  stat[i] = (pval$p.value)
}
df = data.frame(stat = unlist(stat))  

plot = ggplot(df, aes(x = stat)) 
plot+ geom_histogram(binwidth=.01) + geom_vline(xintercept = (0.0103), color = "red", size=1) +
  theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('slope_opt_state_mag.pdf')

sum(df$stat <= 0.0103)

###################
# Normality Tests
####################

par(mfrow=c(1,3))
qqPlot(unique(filter(data, band == 'alpha')$diff3))
qqPlot(unique(filter(data, band == 'beta')$diff3))
qqPlot(unique(filter(data, band == 'low_gamma')$diff3))

par(mfrow=c(1,3))
qqPlot(unique(filter(data, band == 'alpha')$diffa3))
qqPlot(unique(filter(data, band == 'beta')$diffa3))
qqPlot(unique(filter(data, band == 'low_gamma')$diffa3))

par(mfrow=c(2,2))
qqPlot(unique(filter(data, band == 'alpha')$diffa23))
qqPlot(unique(filter(data, band == 'beta')$diffa13))
qqPlot(unique(filter(data, band == 'low_gamma')$high3_a))
qqPlot(unique(filter(data, band == 'low_gamma')$low_a))


par(mfrow=c(2,3))
qqPlot(unique(filter(data_state, band == 'alpha')$diff_mag))
qqPlot(unique(filter(data_state, band == 'alpha')$diff_2))
qqPlot(unique(filter(data_state, band == 'alpha')$diff_cent))
qqPlot(unique(filter(data_state, band == 'alpha')$diff_012))
qqPlot(unique(filter(data_state, band == 'alpha')$diff_inv))



