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
data = data.frame(band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), high_e = raw_data$e.high[1,], high2_e = raw_data$e.high2[1,],
                  high3_e = raw_data$e.high3[1,], low_e = raw_data$e.low[1,],zero_e = raw_data$e.zero[1,], high_p = raw_data$p.high[1,], high2_p = raw_data$p.high2[1,],
                  high3_p = raw_data$p.high[1,], low_p = raw_data$p.low[1,],zero_p = raw_data$p.zero[1,], high_m = raw_data$m.high[1,], high2_m = raw_data$m.high2[1,],
                  high3_m = raw_data$m.high3[1,], low_m = raw_data$m.low[1,],zero_m = raw_data$m.zero[1,], slope = raw_data$slope[1,])

data_high = select(data, band, subj, slope, high_e, high_p, high_m)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high_e'] = 'e'
names(data_high)[names(data_high)=='high_p'] = 'p'
names(data_high)[names(data_high)=='high_m'] = 'm'

data_high2 = select(data, band, subj, slope, high2_e, high2_p, high2_m)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2_e'] = 'e'
names(data_high2)[names(data_high2)=='high2_p'] = 'p'
names(data_high2)[names(data_high2)=='high2_m'] = 'm'

data_high3 = select(data, band, subj, slope, high3_e, high3_p, high3_m)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3_e'] = 'e'
names(data_high3)[names(data_high3)=='high3_p'] = 'p'
names(data_high3)[names(data_high3)=='high3_m'] = 'm'

data_low = select(data, band, subj, slope, low_e, low_p, low_m)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low_e'] = 'e'
names(data_low)[names(data_low)=='low_p'] = 'p'
names(data_low)[names(data_low)=='low_m'] = 'm'

data_zero = select(data, band, subj, slope, zero_e, zero_p, zero_m)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero_e'] = 'e'
names(data_zero)[names(data_zero)=='zero_p'] = 'p'
names(data_zero)[names(data_zero)=='zero_m'] = 'm'

data2 = rbind(data_high, data_high2, data_high3, data_low, data_zero)
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
plot + geom_violin(aes(fill = band), trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(aes(fill = band), width = 0.15,
    position = position_dodge(0.9)
  ) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",3)) + 
  #scale_fill_manual(values = wes_palette("Set2")) + 
  labs(x = 'Loading', y = 'Peak')  + theme_minimal()
ggsave(paste(sens, '_peak.pdf', sep = ''))

plot = ggplot(data2, aes(x = cond, y = m, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Loading', y = 'Max')  + theme_minimal()
ggsave(paste(sens, '_max', '.png', sep = ''))




## Stats
fit_beta_e = lm(e~cond, data_beta)
anova(fit_beta_e)
fit_beta_p = lm(p~cond, data_beta)
anova(fit_beta_p)
fit_beta_m = lm(m~cond, data_beta)
anova(fit_beta_m)

fit_alpha_e = lm(e~cond, data_alpha)
anova(fit_alpha_e)
fit_alpha_p = lm(p~cond, data_alpha)
anova(fit_alpha_p)
fit_alpha_m = lm(m~cond, data_alpha)
anova(fit_alpha_m)

fit_gamma_e = lm(e~cond, data_gamma)
anova(fit_gamma_e)
fit_gamma_p = lm(p~cond, data_gamma)
anova(fit_gamma_p)
fit_gamma_m = lm(m~cond, data_gamma)
anova(fit_gamma_m)

fit_p = lm(p~cond+band, data2)
anova(fit_p)
fit_e = lm(e~cond+band, data2)
anova(fit_e)
fit_m = lm(m~cond+band, data2)
anova(fit_m)

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
