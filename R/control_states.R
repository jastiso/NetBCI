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
bands = c('alpha', 'beta', 'low_gamma')
nSubj = 20
sens = 'grad'

raw_data = readMat(paste('data/wpli/', sens, '/control_state.mat', sep = ''))
data = data.frame(band = unlist(raw_data$band.ord), region = unlist(raw_data$region.ord), high = raw_data$h.region[1,], high2 = raw_data$h2.region[1,],
                  high3 = raw_data$h3.region[1,], high_ev2 = raw_data$h.region2[1,], low = raw_data$l.region[1,], low_ev2 = raw_data$l.region2[1,], 
                  zero = raw_data$z.region[1,], zero_ev2 = raw_data$z.region2[1,], small_ev = raw_data$small.region[1,], high2_ev2 = raw_data$h2.region2[1,], 
                  high3_ev2 = raw_data$h3.region2[1,] )

plot = ggplot(data, aes(x = region, y = high, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_high_control_state', '.png', sep = ''))

plot = ggplot(data, aes(x = region, y = low, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_low_control_state', '.png', sep = ''))

plot = ggplot(data, aes(x = region, y = high2, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_high2_control_state', '.png', sep = ''))

plot = ggplot(data, aes(x = region, y = high3, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_high3_control_state', '.png', sep = ''))

## stats

# reformat
data_high = select(data, band, region, high)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high'] = 'state'

data_high2 = select(data, band, region, high2)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2'] = 'state'

data_high3 = select(data, band, region, high3)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3'] = 'state'

data_low = select(data, band, region, low)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low'] = 'state'

data_zero = select(data, band, region, zero)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero'] = 'state'

data_small = select(data, band, region, small_ev)
data_small$cond = 'small'
names(data_small)[names(data_small)=='small_ev'] = 'state'

data2 = rbind(data_high, data_low, data_high2, data_high3, data_zero, data_small)

data_alpha = filter(data2, band == 'alpha')
data_beta = filter(data2, band == 'beta')
data_low_gamma = filter(data2, band == 'low_gamma')


# plot sub data
plot = ggplot(data_beta, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'beta_control_state', '.png', sep = ''))

plot = ggplot(filter(data_beta, region == c('Left_motor', "Right_motor")), aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'beta_control_state_region', '.png', sep = ''))

plot = ggplot(data_alpha, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'alpha_control_state', '.png', sep = ''))

plot = ggplot(filter(data_alpha, region == c('Left_motor', "Right_motor")), aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'alpha_control_state_region', '.png', sep = ''))

plot = ggplot(data_low_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, 'low_gamma_control_state', '.png', sep = ''))


# beta should have low motor exp in right

right_beta = t.test(filter(data, region == 'Right_motor', band == 'beta')$high,filter(data, region == 'Right_motor', band == 'beta')$low, paired=TRUE)
right_beta
right_beta_z = t.test(filter(data, region == 'Right_motor', band == 'beta')$high,filter(data, region == 'Right_motor', band == 'beta')$zero, paired=TRUE)
right_beta_z
high_beta = t.test(filter(data, region == 'Left_motor', band == 'beta')$high,filter(data, region == 'Right_motor', band == 'beta')$high, paired = TRUE)
high_beta
low_beta = t.test(filter(data, region == 'Left_motor', band == 'beta')$low,filter(data, region == 'Right_motor', band == 'beta')$low, paired = TRUE)
low_beta

# beta should have more motor exp in right

left_beta = t.test(filter(data, region == 'Left_motor', band == 'beta')$high,filter(data, region == 'Left_motor', band == 'beta')$low, paired = TRUE)
left_beta
left_beta2 = t.test(filter(data, region == 'Left_motor', band == 'beta')$high2,filter(data, region == 'Left_motor', band == 'beta')$low, paired = TRUE)
left_beta2
left_beta3 = t.test(filter(data, region == 'Left_motor', band == 'beta')$high3,filter(data, region == 'Left_motor', band == 'beta')$low, paired = TRUE)
left_beta3
left_beta_z = t.test(filter(data, region == 'Left_motor', band == 'beta')$high,filter(data, region == 'Left_motor', band == 'beta')$zero, paired=TRUE)
left_beta_z
left_beta_z2 = t.test(filter(data, region == 'Left_motor', band == 'beta')$high2,filter(data, region == 'Left_motor', band == 'beta')$zero, paired=TRUE)
left_beta_z2
left_beta_z3 = t.test(filter(data, region == 'Left_motor', band == 'beta')$high3,filter(data, region == 'Left_motor', band == 'beta')$zero, paired=TRUE)
left_beta_z3
high_beta = t.test(filter(data, region == 'Left_motor', band == 'beta')$high,filter(data, region == 'Right_motor', band == 'beta')$high, paired = TRUE)
high_beta
low_beta = t.test(filter(data, region == 'Left_motor', band == 'beta')$low,filter(data, region == 'Right_motor', band == 'beta')$low, paired = TRUE)
low_beta

# alpha should have less motor exp in left
left_alpha = t.test(filter(data, region == 'Left_motor', band == 'alpha')$high,filter(data, region == 'Left_motor', band == 'alpha')$low, paired=TRUE)
left_alpha
right_alpha = t.test(filter(data, region == 'Right_motor', band == 'alpha')$high,filter(data, region == 'Right_motor', band == 'alpha')$low, paired=TRUE)
right_alpha

# gamma should have ???
left_low_gamma = t.test(filter(data, region == 'Left_motor', band == 'low_gamma')$high,filter(data, region == 'Left_motor', band == 'low_gamma')$low, paired=TRUE)
left_low_gamma
right_low_gamma = t.test(filter(data, region == 'Right_motor', band == 'low_gamma')$high,filter(data, region == 'Right_motor', band == 'low_gamma')$low, paired=TRUE)
right_low_gamma
frontal_low_gamma = t.test(filter(data, region == 'Left_frontal', band == 'low_gamma')$high,filter(data, region == 'Left_frontal', band == 'low_gamma')$low, paired=TRUE)
frontal_low_gamma
frontal_low_gamma = t.test(filter(data, region == 'Right_frontal', band == 'low_gamma')$high,filter(data, region == 'Right_frontal', band == 'low_gamma')$low, paired=TRUE)
frontal_low_gamma




#####################################################################

## Uniform Phase Randomization

####################################################################

raw_data = readMat(paste('data/wpli/', sens, '/control_state_pr.mat', sep = ''))
pr_data = data.frame(band = unlist(raw_data$band.ord), region = unlist(raw_data$region.ord), high = raw_data$h.region[1,], high2 = raw_data$h2.region[1,],
                     high3 = raw_data$h3.region[1,], high_ev2 = raw_data$h.region2[1,], low = raw_data$l.region[1,], low_ev2 = raw_data$l.region2[1,],
                     zero = raw_data$z.region[1,], zero_ev2 = raw_data$z.region2[1,])

plot = ggplot(pr_data, aes(x = region, y = high, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = brewer.pal(4, "Greys")) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_high_control_state', '.png', sep = ''))

plot = ggplot(pr_data, aes(x = region, y = low, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = brewer.pal(4, "Greys")) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_low_control_state', '.png', sep = ''))






## stats

# reformat
pr_data_high = select(pr_data, band, region, high)
pr_data_high$cond = 'high'
names(pr_data_high)[names(pr_data_high)=='high'] = 'state'

pr_data_high2 = select(pr_data, band, region, high2)
pr_data_high2$cond = 'high2'
names(pr_data_high2)[names(pr_data_high2)=='high2'] = 'state'

pr_data_high3 = select(pr_data, band, region, high3)
pr_data_high3$cond = 'high3'
names(pr_data_high3)[names(pr_data_high3)=='high3'] = 'state'

pr_data_low = select(pr_data, band, region, low)
pr_data_low$cond = 'low'
names(pr_data_low)[names(pr_data_low)=='low'] = 'state'

pr_data_zero = select(pr_data, band, region, zero)
pr_data_zero$cond = 'zero'
names(pr_data_zero)[names(pr_data_zero)=='zero'] = 'state'

pr_data2 = rbind(pr_data_high, pr_data_low, pr_data_high2, pr_data_high3, pr_data_zero)

pr_data_alpha = filter(pr_data2, band == 'alpha')
pr_data_beta = filter(pr_data2, band == 'beta')
pr_data_low_gamma = filter(pr_data2, band == 'low_gamma')

# plot sub data
plot = ggplot(pr_data_beta, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_beta_control_state', '.png', sep = ''))


plot = ggplot(pr_data_alpha, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_alpha_control_state', '.png', sep = ''))

plot = ggplot(pr_data_low_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_pr_low_gamma_control_state', '.png', sep = ''))




# beta should have more motor exp in left
left_beta = t.test(filter(pr_data, region == 'Left_motor', band == 'beta')$high,filter(pr_data, region == 'Left_motor', band == 'beta')$low, paired=TRUE)
left_beta
right_beta = t.test(filter(pr_data, region == 'Right_motor', band == 'beta')$high,filter(pr_data, region == 'Right_motor', band == 'beta')$low, paired=TRUE)
right_beta
right_beta_z = t.test(filter(pr_data, region == 'Right_motor', band == 'beta')$high,filter(pr_data, region == 'Right_motor', band == 'beta')$zero, paired=TRUE)
right_beta_z
high_beta = t.test(filter(pr_data, region == 'Right_motor', band == 'beta')$high,filter(pr_data, region == 'Left_motor', band == 'beta')$high, paired=TRUE)
high_beta

# alpha should have less motor exp in left
left_alpha = t.test(filter(pr_data, region == 'Left_motor', band == 'alpha')$high,filter(pr_data, region == 'Left_motor', band == 'alpha')$low, paired=TRUE)
left_alpha
right_alpha = t.test(filter(pr_data, region == 'Right_motor', band == 'alpha')$high,filter(pr_data, region == 'Right_motor', band == 'alpha')$low, paired=TRUE)
right_alpha

# gamma should ??
fit_low_gamma = lm(state ~ region + cond, data_low_gamma)
anova(fit_low_gamma)


# cross model test
null_beta = t.test(filter(pr_data, region == 'Left_motor', band == 'beta')$high,filter(data, region == 'Left_motor', band == 'beta')$high, paired=TRUE)
null_beta
null_beta_low = t.test(filter(pr_data, region == 'Right_motor', band == 'beta')$low,filter(data, region == 'Right_motor', band == 'beta')$low, paired=TRUE)
null_beta_low
null_beta_zero = t.test(filter(pr_data, region == 'Right_motor', band == 'beta')$low,filter(data, region == 'Right_motor', band == 'beta')$zero, paired=TRUE)
null_beta_zero
null_alpha = t.test(filter(pr_data, region == 'Left_motor', band == 'alpha')$high,filter(data, region == 'Left_motor', band == 'alpha')$high, paired=TRUE)
null_alpha







#####################################################################

## Eigenvector2

####################################################################

# reformat
data_high = select(data, band, region, high_ev2)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high_ev2'] = 'state'

data_high2 = select(data, band, region, high2_ev2)
data_high2$cond = 'high2'
names(data_high2)[names(data_high2)=='high2_ev2'] = 'state'

data_high3 = select(data, band, region, high3_ev2)
data_high3$cond = 'high3'
names(data_high3)[names(data_high3)=='high3_ev2'] = 'state'

data_low = select(data,  band, region, low_ev2)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low_ev2'] = 'state'

data_zero = select(data,  band, region, zero_ev2)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero_ev2'] = 'state'

data2 = rbind(data_high, data_low, data_zero, data_high2, data_high3)

data_alpha = filter(data2, band == 'alpha')
data_beta = filter(data2, band == 'beta')
data_low_gamma = filter(data2, band == 'low_gamma')

# plot sub data
plot = ggplot(data_beta, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_beta_ev2', '.png', sep = ''))


plot = ggplot(data_alpha, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_alpha_ev2', '.png', sep = ''))

plot = ggplot(data_low_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_low_gamma_ev2', '.png', sep = ''))


# alpha should have less motor exp in left
stat_alpha = t.test(filter(data, region == 'Left_motor', band == 'alpha')$high_ev2,filter(data, region == 'Left_motor', band == 'alpha')$low_ev2, paired=TRUE)
stat_alpha
stat_alpha_z = t.test(filter(data, region == 'Left_motor', band == 'alpha')$high_ev2,filter(data, region == 'Left_motor', band == 'alpha')$zero_ev2, paired=TRUE)
stat_alpha_z
# beta 
stat_beta = t.test(filter(data, region == 'Left_motor', band == 'beta')$high3_ev2,filter(data, region == 'Left_motor', band == 'beta')$low_ev2, paired=TRUE)
stat_beta
stat_beta_z = t.test(filter(data, region == 'Left_motor', band == 'beta')$high3_ev2,filter(data, region == 'Left_motor', band == 'beta')$zero_ev2, paired=TRUE)
stat_beta_z
stat_gamma = t.test(filter(data, region == 'Left_motor', band == 'low_gamma')$high_ev2,filter(data, region == 'Left_motor', band == 'low_gamma')$low_ev2, paired=TRUE)
stat_gamma
stat_gamma_z = t.test(filter(data, region == 'Left_motor', band == 'low_gamma')$high_ev2,filter(data, region == 'Left_motor', band == 'low_gamma')$zero_ev2, paired=TRUE)
stat_gamma_z







#####################################################################

## Eigenvector2 - phase randomized

####################################################################
# reformat
pr_data_high = select(pr_data, band, region, high_ev2)
pr_data_high$cond = 'high'
names(pr_data_high)[names(pr_data_high)=='high_ev2'] = 'state'

pr_data_low = select(pr_data,  band, region, low_ev2)
pr_data_low$cond = 'low'
names(pr_data_low)[names(pr_data_low)=='low_ev2'] = 'state'

pr_data_zero = select(pr_data,  band, region, zero_ev2)
pr_data_zero$cond = 'zero'
names(pr_data_zero)[names(pr_data_zero)=='zero_ev2'] = 'state'

pr_data2 = rbind(pr_data_high, pr_data_low, pr_data_zero)

pr_data_alpha = filter(pr_data2, band == 'alpha')
pr_data_beta = filter(pr_data2, band == 'beta')
pr_data_low_gamma = filter(pr_data2, band == 'low_gamma')

# plot sub pr_data
plot = ggplot(pr_data_beta, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_beta_ev2', '.png', sep = ''))


plot = ggplot(pr_data_alpha, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_alpha_ev2', '.png', sep = ''))

plot = ggplot(pr_data_low_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_low_gamma_ev2', '.png', sep = ''))


# alpha should have less motor exp in left
stat_alpha = t.test(filter(pr_data, region == 'Left_motor', band == 'alpha')$high_ev2,filter(pr_data, region == 'Left_motor', band == 'alpha')$low_ev2, paired=TRUE)
stat_alpha
stat_alpha_zero = t.test(filter(pr_data, region == 'Left_motor', band == 'alpha')$high_ev2,filter(pr_data, region == 'Left_motor', band == 'alpha')$zero_ev2, paired=TRUE)
stat_alpha_zero
# beta more
stat_beta = t.test(filter(pr_data, region == 'Left_motor', band == 'beta')$high,filter(pr_data, region == 'Left_motor', band == 'beta')$low, paired=TRUE)
stat_beta








#####################################################################

## Different Control Set

####################################################################


raw_data = readMat(paste('data/wpli/', sens, '/cont_control_state.mat', sep = ''))
cont_data = data.frame(band = unlist(raw_data$band.ord), region = unlist(raw_data$region.ord), high = raw_data$cont.h.region[1,], low = raw_data$cont.l.region[1,],
                       high_ev2 = raw_data$cont.h.region2[1,], low_ev2 = raw_data$cont.l.region2[1,],
                       zero_ev2 = raw_data$cont.z.region2[1,])

# reformat
data_high = select(cont_data, band, region, high_ev2)
data_high$cond = 'high'
names(data_high)[names(data_high)=='high_ev2'] = 'state'

data_low = select(cont_data,  band, region, low_ev2)
data_low$cond = 'low'
names(data_low)[names(data_low)=='low_ev2'] = 'state'

data_zero = select(cont_data,  band, region, zero_ev2)
data_zero$cond = 'zero'
names(data_zero)[names(data_zero)=='zero_ev2'] = 'state'

data2 = rbind(data_high, data_low, data_zero)

data_alpha = filter(data2, band == 'alpha')
data_beta = filter(data2, band == 'beta')
data_low_gamma = filter(data2, band == 'low_gamma')

# plot sub data
plot = ggplot(data_beta, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_null_beta_ev2', '.png', sep = ''))


plot = ggplot(data_alpha, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_null_alpha_ev2', '.png', sep = ''))

plot = ggplot(data_low_gamma, aes(x = region, y = state, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'State Value')  + theme_minimal()
ggsave(paste(sens, '_null_low_gamma_ev2', '.png', sep = ''))




# beta should have more motor exp in left
control_beta = t.test(filter(cont_data, region == 'Right_motor', band == 'beta')$high_ev2,filter(cont_data, region == 'Right_motor', band == 'beta')$low_ev2, paired=TRUE)
control_beta


# alpha should have less motor exp in left
stat_alpha = t.test(filter(cont_data, region == 'Right_motor', band == 'alpha')$high_ev2,filter(cont_data, region == 'Right_motor', band == 'alpha')$low_ev2, paired=TRUE)
stat_alpha
stat_alpha_z = t.test(filter(cont_data, region == 'Right_motor', band == 'alpha')$high_ev2,filter(cont_data, region == 'Right_motor', band == 'alpha')$zero_ev2, paired=TRUE)
stat_alpha_z
