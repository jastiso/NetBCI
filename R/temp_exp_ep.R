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
data = data.frame(band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), slope = raw_data$slope[1,], high_e = raw_data$e.high[1,], high2_e = raw_data$e.high2[1,],
                  high3_e = raw_data$e.high3[1,], low_e = raw_data$e.low[1,], high_p = raw_data$p.high[1,], high2_p = raw_data$p.high2[1,],
                  high3_p = raw_data$p.high3[1,], low_p = raw_data$p.low[1,])

data_high = select(data, band, subj, slope, high_e, high_p)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high_e'] = 'e'
names(data_high)[names(data_high)=='high_p'] = 'p'

data_high2 = select(data, band, subj, slope, high2_e, high2_p)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2_e'] = 'e'
names(data_high2)[names(data_high2)=='high2_p'] = 'p'

data_high3 = select(data, band, subj, slope, high3_e, high3_p)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3_e'] = 'e'
names(data_high3)[names(data_high3)=='high3_p'] = 'p'

data_low = select(data, band, subj, slope, low_e, low_p)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low_e'] = 'e'
names(data_low)[names(data_low)=='low_p'] = 'p'

data2 = rbind(data_high, data_high2, data_high3, data_low)
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

data2 = mutate(data2, label = paste(data2$band,data2$cond))
plot = ggplot(data2, aes(x = p, y = label, color = band) )
plot + geom_point( color = "black",size=6) + 
  geom_point(size=5) + 
  #geom_dotplot(binaxis='y', stackdir=, 'center', dotsize=.5)
  scale_color_manual(values = wes_palette("Royal1",3)) + 
  #scale_fill_manual(values = "black") + 
  labs(x = 'Loading', y = 'Peak')  + theme_minimal()
ggsave(paste(sens, '_peak3.svg', sep = ''))





## Stats
fitp = lm(p ~ cond + band, data2)
summary(fitp)
anova(fitp)

fite = lm(e ~ cond + band, data2)
summary(fite)
anova(fite)

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

fit_gamma_e = lm(e~cond, data_gamma)
anova(fit_gamma_e)
fit_gamma_p = lm(p~cond, data_gamma)

fit_p = lm(p~cond+band, data2)
anova(fit_p)
fit_e = lm(e~cond+band, data2)
anova(fit_e)

# corr with slope
beta_corr = cor.test((filter(data, band == 'beta')$high3_p),(filter(data, band == 'beta')$slope))
beta_corr

alpha_corr = cor.test((filter(data, band == 'alpha')$high3_p),(filter(data, band == 'alpha')$slope))
alpha_corr

gamma_corr = cor.test((filter(data, band == 'low_gamma')$high3_p),(filter(data, band == 'low_gamma')$slope))
gamma_corr

# ttest
beta_peak = t.test((filter(data, band == 'beta')$high2_p),log(filter(data, band == 'beta')$high3_p), paired=TRUE)
beta_peak

alpha_peak = t.test(log(filter(data, band == 'alpha')$high2_p),log(filter(data, band == 'alpha')$high3_p), paired=TRUE)
alpha_peak

gamma_peak = t.test(log(filter(data, band == 'low_gamma')$high2_p),log(filter(data, band == 'low_gamma')$high3_p), paired=TRUE)
gamma_peak



#######################

#UPR

#######################
raw_data = readMat(paste('data/wpli/', sens, '/temp_exp_pr.mat', sep = ''))
data_pr = data.frame(band = unlist(raw_data$band.order), subj = unlist(raw_data$subj.order), slope = raw_data$slope[1,], high_e = raw_data$e.high[1,], high2_e = raw_data$e.high2[1,],
                  high3_e = raw_data$e.high3[1,], low_e = raw_data$e.low[1,], high_p = raw_data$p.high[1,], high2_p = raw_data$p.high2[1,],
                  high3_p = raw_data$p.high3[1,], low_p = raw_data$p.low[1,])


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

data_pr_alpha = filter(data_pr2, band == "alpha")
data_pr_beta = filter(data_pr2, band == "beta")
data_pr_gamma = filter(data_pr2, band == "low_gamma")


# stats
fit_pr = lm(p ~ cond + band, data_pr2)
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

