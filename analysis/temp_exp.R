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

raw_data = readMat(paste('data/wpli/', sens, '/temp_exp.mat', sep = ''))
dens = readMat(paste('data/wpli/', sens, '/density.mat', sep = ''))
dens_data = data.frame(band = unlist(dens$band.order), subj = unlist(dens$subj.order), high = dens$d.high[,1], high2 = dens$d.high2[,1], high3 = dens$d.high3[,1], low = dens$d.low[,1])
data = data.frame(band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), high_e = raw_data$e.high[1,], high2_e = raw_data$e.high2[1,],
                  high3_e = raw_data$e.high3[1,], low_e = raw_data$e.low[1,],zero_e = raw_data$e.zero[1,], high_p = raw_data$p.high[1,], high2_p = raw_data$p.high2[1,],
                  high3_p = raw_data$p.high3[1,], low_p = raw_data$p.low[1,],zero_p = raw_data$p.zero[1,], slope = raw_data$slope[1,])
# check order is correct - want both to be false
any(dens_data$band != data$band)
any(dens_data$subj != data$subj)

data_high = select(data, band, subj, slope, high_e, high_p)
data_high$cond = 'high'
data_high$dens = dens_data$high
names(data_high)[names(data_high)=='high_e'] = 'e'
names(data_high)[names(data_high)=='high_p'] = 'p'

data_high2 = select(data, band, subj, slope, high2_e, high2_p)
data_high2$cond = 'high2'
data_high2$dens = dens_data$high2
names(data_high2)[names(data_high2)=='high2_e'] = 'e'
names(data_high2)[names(data_high2)=='high2_p'] = 'p'

data_high3 = select(data, band, subj, slope, high3_e, high3_p)
data_high3$cond = 'high3'
data_high3$dens = dens_data$high3
names(data_high3)[names(data_high3)=='high3_e'] = 'e'
names(data_high3)[names(data_high3)=='high3_p'] = 'p'

data_low = select(data, band, subj, slope, low_e, low_p)
data_low$cond = 'low'
data_low$dens = dens_data$low
names(data_low)[names(data_low)=='low_e'] = 'e'
names(data_low)[names(data_low)=='low_p'] = 'p'

data2 = rbind(data_high, data_high2, data_high3, data_low)
data2 = mutate(data2,hl = cond == c('high','high2','high3'))

data_alpha = filter(data2, band == "alpha")
data_beta = filter(data2, band == "beta")
data_gamma = filter(data2, band == "low_gamma")


## Plots

plot = ggplot(data2, aes(x = band, y = log(e), fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Moonrise3",5)) + 
  labs(x = 'Loading', y = 'Energy')  + theme_minimal()
ggsave(paste(sens, '_energy', '.png', sep = ''))

plot = ggplot(data2, aes(x = cond, y = p) )
plot + geom_violin(aes(fill = band), trim = FALSE, position = position_dodge(0.75)) + 
  geom_boxplot(aes(fill = band), width = 0.15,
    position = position_dodge(0.75)
  ) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",3)) + 
  #scale_fill_manual(values = wes_palette("Set2")) + 
  labs(x = 'Loading', y = 'Peak')  + theme_minimal()
ggsave(paste(sens, '_peak.pdf', sep = ''))


##
plot = ggplot(data2, aes(x = dens, y = p, color = cond) )
plot + geom_point(size = 4) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'dens', y = 'peak')  + theme_minimal()

## Stats
fit = lmp(p ~ cond + band + dens, data2)
summary(fit)
anova(fit)

fit_beta_e = lm(e~cond, data_beta)
anova(fit_beta_e)
fit_beta_p = lm(p~cond, data_beta)
anova(fit_beta_p)

fit_alpha_e = lm(e~cond, data_alpha)
anova(fit_alpha_e)
fit_alpha_p = lm(p~cond, data_alpha)
anova(fit_alpha_p)

fit_gamma_e = lm(e~cond, data_gamma)
anova(fit_gamma_e)
fit_gamma_p = lm(p~cond, data_gamma)
anova(fit_gamma_p)

fit_p = lm(p~cond+band, data2)
anova(fit_p)
fit_e = lm(e~cond+band, data2)
anova(fit_e)

# is high later than low?
fit_p = lm(p~hl+band, data2)
anova(fit_p)
summary(fit_p)

# ttest
beta_peak = t.test(log(filter(data, band == 'beta')$high2_p),log(filter(data, band == 'beta')$high3_p), paired=TRUE)
beta_peak

alpha_peak = t.test(log(filter(data, band == 'alpha')$high2_p),log(filter(data, band == 'alpha')$high3_p), paired=TRUE)
alpha_peak

gamma_peak = t.test(log(filter(data, band == 'low_gamma')$high2_p),log(filter(data, band == 'low_gamma')$high3_p), paired=TRUE)
gamma_peak



#######################

#UPR

#######################
raw_data = readMat(paste('data/wpli/', sens, '/temp_exp_pr.mat', sep = ''))
data_pr = data.frame(band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), high_e = raw_data$e.high[1,], high2_e = raw_data$e.high2[1,],
                  high3_e = raw_data$e.high3[1,], low_e = raw_data$e.low[1,],zero_e = raw_data$e.zero[1,], high_p = raw_data$p.high[1,], high2_p = raw_data$p.high2[1,],
                  high3_p = raw_data$p.high3[1,], low_p = raw_data$p.low[1,],zero_p = raw_data$p.zero[1,], slope = raw_data$slope[1,])


data_pr_high = select(data_pr, band, subj, slope, high_e, high_p)
data_pr_high$cond = 'high'
names(data_pr_high)[names(data_pr_high)=='high_e'] = 'e'
names(data_pr_high)[names(data_pr_high)=='high_p'] = 'p'

data_pr_high2 = select(data_pr, band, subj, slope, high2_e, high2_p)
data_pr_high2$cond = 'high2'
names(data_pr_high2)[names(data_pr_high2)=='high2_e'] = 'e'
names(data_pr_high2)[names(data_pr_high2)=='high2_p'] = 'p'

data_pr_high3 = select(data_pr, band, subj, slope, high3_e, high3_p)
data_pr_high3$cond = 'high3'
names(data_pr_high3)[names(data_pr_high3)=='high3_e'] = 'e'
names(data_pr_high3)[names(data_pr_high3)=='high3_p'] = 'p'

data_pr_low = select(data_pr, band, subj, slope, low_e, low_p)
data_pr_low$cond = 'low'
names(data_pr_low)[names(data_pr_low)=='low_e'] = 'e'
names(data_pr_low)[names(data_pr_low)=='low_p'] = 'p'

data_pr2 = rbind(data_pr_high, data_pr_high2, data_pr_high3, data_pr_low)
data_pr2 = mutate(data_pr2,hl = cond == c('high','high2','high3'))

data_pr_alpha = filter(data_pr2, band == "alpha")
data_pr_beta = filter(data_pr2, band == "beta")
data_pr_gamma = filter(data_pr2, band == "low_gamma")


# stats
fit_pr = lmp(p ~ cond + band, data_pr2)
summary(fit_pr)
anova(fit_pr)

plot = ggplot(data_pr2, aes(x = cond, y = p) )
plot + geom_violin(aes(fill = band), trim = FALSE, position = position_dodge(0.75)) + 
  geom_boxplot(aes(fill = band), width = 0.15,
               position = position_dodge(0.75)
  ) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = brewer.pal(3,"Greys")) + 
  #scale_fill_manual(values = wes_palette("Set2")) + 
  labs(x = 'Loading', y = 'Peak')  + theme_minimal()
ggsave(paste(sens, '_peak_upr.pdf', sep = ''))

