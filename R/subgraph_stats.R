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

raw_data = readMat(paste('data/wpli/subgraph_stats.mat', sep = ''))
nObs = length(raw_data$Vertex.wi)
#data = data_frame(slope = raw_data$slope[,1], imp = raw_data$imp[,1], difference = raw_data$diff[,1], noise_sg = as.factor(raw_data$noise.sg[,1]), 
#                  n = raw_data$sg.n[,1], subj = as.factor(unlist(raw_data$subjects)), band = unlist(raw_data$band), cond = unlist(raw_data$cond), behave_coeff = raw_data$bc[,1], temp_energy = raw_data$temp.energy[,1], 
#                  coeff_med = raw_data$coeff.med[,1], coeff_autocorr = raw_data$coeff.autocorr[,1], coeff_skew = raw_data$coeff.skew[,1], coeff_min = raw_data$coeff.min[,1], 
#                  coeff_quant = raw_data$coeff.quant[,1], lf_wi = raw_data$Left.frontal.wi[,1], lf_bw = raw_data$Left.frontal.bw[,1],
#                  lo_wi = raw_data$Left.occipital.wi[,1], lo_bw = raw_data$Left.occipital.bw[,1], lp_wi = raw_data$Left.parietal.wi[,1], lp_bw = raw_data$Left.parietal.bw[,1],
#                  lt_wi = raw_data$Left.temporal.wi[,1], lt_bw = raw_data$Left.temporal.bw[,1], v_wi = raw_data$Vertex.wi[,1], v_bw = raw_data$Vertex.bw[,1], rf_wi = raw_data$Right.frontal.wi[,1], 
#                  rf_bw = raw_data$Right.frontal.bw[,1], ro_wi = raw_data$Right.occipital.wi[,1], ro_bw = raw_data$Right.occipital.bw[,1], rp_wi = raw_data$Right.parietal.wi[,1],
#                  rp_bw = raw_data$Right.parietal.bw[,1], rt_wi = raw_data$Right.temporal.wi[,1], rt_bw = raw_data$Right.temporal.bw[,1], rfv = raw_data$rfv[,1],
#                  lfv = raw_data$lfv[,1])

data = data_frame(slope = raw_data$slope[,1], imp = raw_data$imp[,1], difference = raw_data$diff[,1], n = raw_data$sg.n[,1], 
                  subj = as.factor(unlist(raw_data$subjects)), band = unlist(raw_data$band), cond = unlist(raw_data$cond), behave_coeff = raw_data$bc[,1], 
                  temp_energy = raw_data$temp.energy[,1], coeff_med = raw_data$coeff.med[,1], coeff_autocorr = raw_data$coeff.autocorr[,1], 
                  coeff_skew = raw_data$coeff.skew[,1], coeff_min = raw_data$coeff.min[,1], coeff_quant = raw_data$coeff.quant[,1], 
                  lf = raw_data$Left.frontal[,1],lo = raw_data$Left.occipital[,1],lm = raw_data$Left.motor[,1], rf = raw_data$Right.frontal[,1], 
                  ro = raw_data$Right.occipital[,1], rm = raw_data$Right.motor[,1] )

# Add fields
#data = dplyr::mutate(group_by(data,subj,band), noise_sg = temp_energy == max(temp_energy))
#data = dplyr::mutate(data, find_noise_bc = (noise_sg == TRUE & cond == 'high'))





# Scatter Plots - Temporal Expression and Coefficients
scatterplot = ggplot(data, aes(x = log(temp_energy), y = coeff_skew, color = cond))
scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Coefficient Skew') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('energy_skew_subj.png')                  

scatterplot = ggplot(data, aes(x = coeff_skew, y = coeff_med, color = cond))
scatterplot + geom_point(size = 6) + labs(x = 'Coeff_Skew', y = 'Coefficient Median') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('skew_med_subj.png')

scatterplot = ggplot(data, aes(x = coeff_skew, y = coeff_min, color = cond))
scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Coefficient Min') + 
   theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('energy_min_subj.png')

scatterplot = ggplot(data, aes(x = coeff_skew, y = coeff_quant, color = cond))
scatterplot + geom_point(size = 6) + labs(x = 'skewness', y = 'Coefficient Quantile') + 
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('skew_quant_subj.png')

scatterplot = ggplot(data, aes(x = log(temp_energy), y = coeff_autocorr, color = cond))
scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Autocorrelation') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('energy_ac_subj.png') 



# Regional Strength
# reformat data
data_rm = select(data, slope, subj, band, cond, rm)
data_rm$region = 'right_motor'
names(data_rm)[names(data_rm)=='rm'] = 'strength'

data_lm = select(data, slope, subj, band, cond, lm)
data_lm$region = 'left_motor'
names(data_lm)[names(data_lm)=='lm'] = 'strength'

data_ro = select(data, slope, subj, band, cond, ro)
data_ro$region = 'right_occipital'
names(data_ro)[names(data_ro)=='ro'] = 'strength'

data_lo = select(data, slope, subj, band, cond, lo)
data_lo$region = 'left_occipital'
names(data_lo)[names(data_lo)=='lo'] = 'strength'

data_rf = select(data, slope, subj, band, cond, rf)
data_rf$region = 'right_frontal'
names(data_rf)[names(data_rf)=='rf'] = 'strength'

data_lf = select(data, slope, subj, band, cond, lf)
data_lf$region = 'left_frontal'
names(data_lf)[names(data_lf)=='lf'] = 'strength'

data2 = rbind(data_rm, data_lm, data_ro, data_lo, data_rf, data_lf)

data_alpha = filter(data2, band == 'alpha')
data_beta = filter(data2, band == 'beta')
data_low_gamma = filter(data2, band == 'low_gamma')


plot = ggplot(data_alpha, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_alpha_str', '.png', sep = ''))

plot = ggplot(data_beta, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_beta_str', '.png', sep = ''))

plot = ggplot(data_low_gamma, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_gamma_str', '.png', sep = ''))



# stats - across regions
# alpha
fit_alpha_high = lm(strength ~ region, filter(data_alpha, cond == "high"))
anova(fit_alpha_high)
summary(fit_alpha_high)

fit_alpha_high2 = lm(strength ~ region, filter(data_alpha, cond == "high2"))
anova(fit_alpha_high2)
summary(fit_alpha_high2)

fit_alpha_high3 = lm(strength ~ region, filter(data_alpha, cond == "high3"))
anova(fit_alpha_high3)
summary(fit_alpha_high3)

fit_alpha_low = lm(strength ~ region, filter(data_alpha, cond == "low"))
anova(fit_alpha_low)
summary(fit_alpha_low)

fit_alpha_zero = lm(strength ~ region, filter(data_alpha, cond == "zero"))
anova(fit_alpha_zero)
summary(fit_alpha_zero)

# beta
fit_beta_high = lm(strength ~ region, filter(data_beta, cond == "high"))
anova(fit_beta_high)
summary(fit_beta_high)

fit_beta_high2 = lm(strength ~ region, filter(data_beta, cond == "high2"))
anova(fit_beta_high2)
summary(fit_beta_high2)

fit_beta_high3 = lm(strength ~ region, filter(data_beta, cond == "high3"))
anova(fit_beta_high3)
summary(fit_beta_high3)

fit_beta_low = lm(strength ~ region, filter(data_beta, cond == "low"))
anova(fit_beta_low)
summary(fit_beta_low)

fit_beta_zero = lm(strength ~ region, filter(data_beta, cond == "zero"))
anova(fit_beta_zero)
summary(fit_beta_zero)

# gamma
fit_low_gamma_high = lm(strength ~ region, filter(data_low_gamma, cond == "high"))
anova(fit_low_gamma_high)
summary(fit_low_gamma_high)

fit_low_gamma_high2 = lm(strength ~ region, filter(data_low_gamma, cond == "high2"))
anova(fit_low_gamma_high2)
summary(fit_low_gamma_high2)

fit_low_gamma_high3 = lm(strength ~ region, filter(data_low_gamma, cond == "high3"))
anova(fit_low_gamma_high3)
summary(fit_low_gamma_high3)

fit_low_gamma_low = lm(strength ~ region, filter(data_low_gamma, cond == "low"))
anova(fit_low_gamma_low)
summary(fit_low_gamma_low)

fit_low_gamma_zero = lm(strength ~ region, filter(data_low_gamma, cond == "zero"))
anova(fit_low_gamma_zero)
summary(fit_low_gamma_zero)

# make plots of only regions with significant regional differences
curr_alpha = filter(data_alpha, cond == "high")
plot = ggplot(curr_alpha, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_alpha_str_sub', '.png', sep = ''))

curr_beta = filter(data_beta, cond == c("high2", "high3", "low"))
plot = ggplot(curr_beta, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_beta_str_sub', '.png', sep = ''))

curr_gamma = filter(data_low_gamma, cond == c("high", "high3", "low"))
plot = ggplot(curr_gamma, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_gamma_str_sub', '.png', sep = ''))




# Relationship between regional expression and BC
scatterplot = ggplot(data, aes(x = lf, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(data, aes(x = rf, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(data, aes(x = lo, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(data, aes(x = ro, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(data, aes(x = lm, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(data, aes(x = rm, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')


#scatterplot = ggplot(filter(data, cond == 'high'), aes(x = lf_bw, y = slope, color = band))
#scatterplot + geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
#   theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
#ggsave('region_slope.png')





# Data without noise SG
#data_clean = filter(data, noise_sg == FALSE)
#scatterplot = ggplot(data_clean, aes(x = log(temp_energy), y = coeff_skew, color = subj))
#scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Coefficient Skew') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
#ggsave('energy_skew_subj_clean.png') 


#scatterplot = ggplot(data_clean, aes(x = behave_coeff, y = coeff_med, color = band))
#scatterplot + geom_point(size = 6) + labs(x = 'Behavioral Coeff', y = 'Coefficient Med') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
#ggsave('energy_bc_band_clean.png') 

#scatterplot = ggplot(data_clean, aes(x = rf_bw, y = behave_coeff, color = band))
#scatterplot + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))




#####################################################################

## Uniform Phase Randomization

####################################################################


raw_data = readMat(paste('data/wpli/pr_subgraph_stats.mat', sep = ''))
#pr_data = data_frame(slope = raw_data$slope[,1], imp = raw_data$imp[,1], difference = raw_data$diff[,1], noise_sg = as.factor(raw_data$noise.sg[,1]), 
#                  n = raw_data$sg.n[,1], subj = as.factor(unlist(raw_data$subjects)), band = unlist(raw_data$band), cond = unlist(raw_data$cond), behave_coeff = raw_data$bc[,1], temp_energy = raw_data$temp.energy[,1], 
#                  coeff_med = raw_data$coeff.med[,1], coeff_autocorr = raw_data$coeff.autocorr[,1], coeff_skew = raw_data$coeff.skew[,1], coeff_min = raw_data$coeff.min[,1], 
#                  coeff_quant = raw_data$coeff.quant[,1], lf_wi = raw_data$Left.frontal.wi[,1], lf_bw = raw_data$Left.frontal.bw[,1],
#                  lo_wi = raw_data$Left.occipital.wi[,1], lo_bw = raw_data$Left.occipital.bw[,1], lp_wi = raw_data$Left.parietal.wi[,1], lp_bw = raw_data$Left.parietal.bw[,1],
#                  lt_wi = raw_data$Left.temporal.wi[,1], lt_bw = raw_data$Left.temporal.bw[,1], v_wi = raw_data$Vertex.wi[,1], v_bw = raw_data$Vertex.bw[,1], rf_wi = raw_data$Right.frontal.wi[,1], 
#                  rf_bw = raw_data$Right.frontal.bw[,1], ro_wi = raw_data$Right.occipital.wi[,1], ro_bw = raw_data$Right.occipital.bw[,1], rp_wi = raw_data$Right.parietal.wi[,1],
#                  rp_bw = raw_data$Right.parietal.bw[,1], rt_wi = raw_data$Right.temporal.wi[,1], rt_bw = raw_data$Right.temporal.bw[,1], rfv = raw_data$rfv[,1],
#                  lfv = raw_data$lfv[,1])
pr_data = data_frame(slope = raw_data$slope[,1], imp = raw_data$imp[,1], difference = raw_data$diff[,1], n = raw_data$sg.n[,1], 
                  subj = as.factor(unlist(raw_data$subjects)), band = unlist(raw_data$band), cond = unlist(raw_data$cond), behave_coeff = raw_data$bc[,1], 
                  temp_energy = raw_data$temp.energy[,1], coeff_med = raw_data$coeff.med[,1], coeff_autocorr = raw_data$coeff.autocorr[,1], 
                  coeff_skew = raw_data$coeff.skew[,1], coeff_min = raw_data$coeff.min[,1], coeff_quant = raw_data$coeff.quant[,1], 
                  lf = raw_data$Left.frontal[,1],lo = raw_data$Left.occipital[,1],lm = raw_data$Left.motor[,1], rf = raw_data$Right.frontal[,1], 
                  ro = raw_data$Right.occipital[,1], rm = raw_data$Right.motor[,1] )


# Regional Strength
# reformat data
pr_data_rm = select(pr_data, slope, subj, band, cond, rm)
pr_data_rm$region = 'right_motor'
names(pr_data_rm)[names(pr_data_rm)=='rm'] = 'strength'

pr_data_lm = select(pr_data, slope, subj, band, cond, lm)
pr_data_lm$region = 'left_motor'
names(pr_data_lm)[names(pr_data_lm)=='lm'] = 'strength'

pr_data_ro = select(pr_data, slope, subj, band, cond, ro)
pr_data_ro$region = 'right_occipital'
names(pr_data_ro)[names(pr_data_ro)=='ro'] = 'strength'

pr_data_lo = select(pr_data, slope, subj, band, cond, lo)
pr_data_lo$region = 'left_occipital'
names(pr_data_lo)[names(pr_data_lo)=='lo'] = 'strength'

pr_data_rf = select(pr_data, slope, subj, band, cond, rf)
pr_data_rf$region = 'right_frontal'
names(pr_data_rf)[names(pr_data_rf)=='rf'] = 'strength'

pr_data_lf = select(pr_data, slope, subj, band, cond, lf)
pr_data_lf$region = 'left_frontal'
names(pr_data_lf)[names(pr_data_lf)=='lf'] = 'strength'

pr_data2 = rbind(pr_data_rm, pr_data_lm, pr_data_ro, pr_data_lo, pr_data_rf, pr_data_lf)

pr_data_alpha = filter(pr_data2, band == 'alpha')
pr_data_beta = filter(pr_data2, band == 'beta')
pr_data_low_gamma = filter(pr_data2, band == 'low_gamma')


plot = ggplot(pr_data_alpha, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_alpha_str_pr', '.png', sep = ''))

plot = ggplot(pr_data_beta, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_beta_str_pr', '.png', sep = ''))

plot = ggplot(pr_data_low_gamma, aes(x = region, y = strength, fill = cond) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  #scale_fill_manual(values = wes_palette("Moonrise1",4)) + 
  labs(x = 'Region', y = 'Strength')  + theme_minimal()
ggsave(paste(sens, '_gamma_str_pr', '.png', sep = ''))



# stats
# alpha
fit_alpha_high = lm(strength ~ region, filter(pr_data_alpha, cond == "high"))
anova(fit_alpha_high)
summary(fit_alpha_high)

fit_alpha_high2 = lm(strength ~ region, filter(pr_data_alpha, cond == "high2"))
anova(fit_alpha_high2)
summary(fit_alpha_high2)

fit_alpha_high3 = lm(strength ~ region, filter(pr_data_alpha, cond == "high3"))
anova(fit_alpha_high3)
summary(fit_alpha_high3)

fit_alpha_low = lm(strength ~ region, filter(pr_data_alpha, cond == "low"))
anova(fit_alpha_low)
summary(fit_alpha_low)

fit_alpha_zero = lm(strength ~ region, filter(pr_data_alpha, cond == "zero"))
anova(fit_alpha_zero)
summary(fit_alpha_zero)

# beta
fit_beta_high = lm(strength ~ region, filter(pr_data_beta, cond == "high"))
anova(fit_beta_high)
summary(fit_beta_high)

fit_beta_high2 = lm(strength ~ region, filter(pr_data_beta, cond == "high2"))
anova(fit_beta_high2)
summary(fit_beta_high2)

fit_beta_high3 = lm(strength ~ region, filter(pr_data_beta, cond == "high3"))
anova(fit_beta_high3)
summary(fit_beta_high3)

fit_beta_low = lm(strength ~ region, filter(pr_data_beta, cond == "low"))
anova(fit_beta_low)
summary(fit_beta_low)

fit_beta_zero = lm(strength ~ region, filter(pr_data_beta, cond == "zero"))
anova(fit_beta_zero)
summary(fit_beta_zero)

# gamma
fit_low_gamma_high = lm(strength ~ region, filter(pr_data_low_gamma, cond == "high"))
anova(fit_low_gamma_high)
summary(fit_low_gamma_high)

fit_low_gamma_high2 = lm(strength ~ region, filter(pr_data_low_gamma, cond == "high2"))
anova(fit_low_gamma_high2)
summary(fit_low_gamma_high2)

fit_low_gamma_high3 = lm(strength ~ region, filter(pr_data_low_gamma, cond == "high3"))
anova(fit_low_gamma_high3)
summary(fit_low_gamma_high3)

fit_low_gamma_low = lm(strength ~ region, filter(pr_data_low_gamma, cond == "low"))
anova(fit_low_gamma_low)
summary(fit_low_gamma_low)

fit_low_gamma_zero = lm(strength ~ region, filter(pr_data_low_gamma, cond == "zero"))
anova(fit_low_gamma_zero)
summary(fit_low_gamma_zero)



# Relationship between regional expression and BC
scatterplot = ggplot(pr_data, aes(x = lf, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(pr_data, aes(x = rf, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(pr_data, aes(x = lo, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(pr_data, aes(x = ro, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(pr_data, aes(x = lm, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')

scatterplot = ggplot(pr_data, aes(x = rm, y = slope, color = band, group = band))
scatterplot + geom_point(size = 6) + labs(x = 'Strength', y = 'Slope') + 
  geom_smooth(method="lm") + theme_minimal() + scale_color_manual(values =  wes_palette("Royal1",4))
ggsave('lf_bc.png')


#scatterplot = ggplot(filter(data, cond == 'high'), aes(x = lf_bw, y = slope, color = band))
#scatterplot + geom_smooth(method="lm") + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
#   theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
#ggsave('region_slope.png')





# pr_data without noise SG
#pr_data_clean = filter(pr_data, noise_sg == 0)
#scatterplot = ggplot(pr_data_clean, aes(x = log(temp_energy), y = coeff_skew, color = subj))
#scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Coefficient Skew') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
#ggsave('pr_energy_skew_subj_clean.png') 


#scatterplot = ggplot(pr_data_clean, aes(x = behave_coeff, y = coeff_med, color = band))
#scatterplot + geom_point(size = 6) + labs(x = 'Behavioral Coeff', y = 'Coefficient Med') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
#ggsave('pr_energy_bc_band_clean.png') 

#scatterplot = ggplot(pr_data_clean, aes(x = rf_bw, y = behave_coeff, color = band))
#scatterplot + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))




# Comparing null models
stat_alpha = t.test(filter(data, cond == 'high', band == 'alpha')$lm, filter(pr_data, cond == 'high', band == 'alpha')$lm, paired=TRUE)
stat_alpha
