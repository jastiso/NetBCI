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



## stats

# reformat
data_high = select(data, -low)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high'] = 'state'
data_low = select(data, -high)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low'] = 'state'
data2 = rbind(data_high, data_low)

data_alpha = filter(data2, band == 'alpha')
data_beta = filter(data2, band == 'beta')
data_low_gamma = filter(data2, band == 'low_gamma')
data_gamma = filter(data2, band == 'gamma')

# plot sub data
plot = ggplot(data_beta, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'beta_control_state', '.png', sep = ''))


plot = ggplot(data_alpha, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'alpha_control_state', '.png', sep = ''))

plot = ggplot(data_low_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'low_gamma_control_state', '.png', sep = ''))

plot = ggplot(data_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'gamma_control_state', '.png', sep = ''))



# beta should have more motor exp in left
stat_beta = t.test(filter(data, region == 'Left_parietal', band == 'beta')$high,filter(data, region == 'Left_parietal', band == 'beta')$low, paired=TRUE)
stat_beta
control_beta = t.test(filter(data, region == 'Right_parietal', band == 'beta')$high,filter(data, region == 'Right_parietal', band == 'beta')$low, paired=TRUE)
control_beta

# alpha should have less motor exp in left
stat_alpha = t.test(filter(data, region == 'Left_parietal', band == 'alpha')$high,filter(data, region == 'Left_parietal', band == 'alpha')$low, paired=TRUE)
stat_alpha
control_alpha = t.test(filter(data, region == 'Right_parietal', band == 'alpha')$high,filter(data, region == 'Right_parietal', band == 'alpha')$low, paired=TRUE)
control_alpha

# gamma should have lateralized frontal exp
fit_low_gamma = lm(state ~ region + cond, data_low_gamma)
anova(fit_low_gamma)

fit_gamma = lm(state ~ region + cond, data_gamma)
anova(fit_gamma)





#####################################################################

## Uniform Phase Randomization

####################################################################

raw_data = readMat(paste('data/wpli/', sens, '/control_state_pr.mat', sep = ''))
pr_data = data.frame(band = unlist(raw_data$band.ord), region = unlist(raw_data$region.ord), high = raw_data$h.region[,1], low = raw_data$l.region[1,])

plot = ggplot(pr_data, aes(x = region, y = high, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_high_control_state', '.png', sep = ''))

plot = ggplot(pr_data, aes(x = region, y = low, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_low_control_state', '.png', sep = ''))






## stats

# reformat
data_high = select(pr_data, -low)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high'] = 'state'
data_low = select(pr_data, -high)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low'] = 'state'
data2 = rbind(data_high, data_low)

data_alpha = filter(data2, band == 'alpha')
data_beta = filter(data2, band == 'beta')
data_low_gamma = filter(data2, band == 'low_gamma')
data_gamma = filter(data2, band == 'gamma')

# plot sub data
plot = ggplot(data_beta, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_beta_control_state', '.png', sep = ''))


plot = ggplot(data_alpha, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_alpha_control_state', '.png', sep = ''))

plot = ggplot(data_low_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_low_gamma_control_state', '.png', sep = ''))

plot = ggplot(data_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_gamma_control_state', '.png', sep = ''))



# beta should have more motor exp in left
stat_beta = t.test(filter(pr_data, region == 'Left_parietal', band == 'beta')$high,filter(pr_data, region == 'Left_parietal', band == 'beta')$low, paired=TRUE)
stat_beta
control_beta = t.test(filter(pr_data, region == 'Right_parietal', band == 'beta')$high,filter(pr_data, region == 'Right_parietal', band == 'beta')$low, paired=TRUE)
control_beta

# alpha should have less motor exp in left
stat_alpha = t.test(filter(pr_data, region == 'Left_parietal', band == 'alpha')$high,filter(pr_data, region == 'Left_parietal', band == 'alpha')$low, paired=TRUE)
stat_alpha
control_alpha = t.test(filter(pr_data, region == 'Right_parietal', band == 'alpha')$high,filter(pr_data, region == 'Right_parietal', band == 'alpha')$low, paired=TRUE)
control_alpha

# gamma should have lateralized frontal exp
fit_low_gamma = lm(state ~ region + cond, data_low_gamma)
anova(fit_low_gamma)

fit_gamma = lm(state ~ region + cond, data_gamma)
anova(fit_gamma)
