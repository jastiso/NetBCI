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

raw_data = readMat(paste('data/wpli/subgraph_stats.mat', sep = ''))
nObs = length(raw_data$Vertex.wi)
data = data_frame(noise_sg = as.factor(raw_data$noise.sg[,1]), n = raw_data$sg.n[,1], subj = as.factor(unlist(raw_data$subjects)), band = unlist(raw_data$band), cond = unlist(raw_data$cond), behave_coeff = raw_data$bc[,1], temp_energy = raw_data$temp.energy[,1], 
                  coeff_med = raw_data$coeff.med[,1], coeff_autocorr = raw_data$coeff.autocorr[,1], coeff_skew = raw_data$coeff.skew[,1], coeff_min = raw_data$coeff.min[,1], coeff_quant = raw_data$coeff.quant[,1], lf_wi = raw_data$Left.frontal.wi[,1], lf_bw = raw_data$Left.frontal.bw[,1],
                  lo_wi = raw_data$Left.occipital.wi[,1], lo_bw = raw_data$Left.occipital.bw[,1], lp_wi = raw_data$Left.parietal.wi[,1], lp_bw = raw_data$Left.parietal.bw[,1],
                  lt_wi = raw_data$Left.temporal.wi[,1], lt_bw = raw_data$Left.temporal.bw[,1], v_wi = raw_data$Vertex.wi[,1], v_bw = raw_data$Vertex.bw[,1], rf_wi = raw_data$Right.frontal.wi[,1], 
                  rf_bw = raw_data$Right.frontal.bw[,1], ro_wi = raw_data$Right.occipital.wi[,1], ro_bw = raw_data$Right.occipital.bw[,1], rp_wi = raw_data$Right.parietal.wi[,1],
                  rp_bw = raw_data$Right.parietal.bw[,1], rt_wi = raw_data$Right.temporal.wi[,1], rt_bw = raw_data$Right.temporal.bw[,1], rfv = raw_data$rfv[,1],
                  lfv = raw_data$lfv[,1])

# Add fields
#data = dplyr::mutate(group_by(data,subj,band), noise_sg = temp_energy == max(temp_energy))
#data = dplyr::mutate(data, find_noise_bc = (noise_sg == TRUE & cond == 'high'))





# Scatter Plots - Temporal Expression and Coefficients
scatterplot = ggplot(data, aes(x = log(temp_energy), y = coeff_skew, color = noise_sg))
scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Coefficient Skew') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('energy_skew_subj.png')                  

scatterplot = ggplot(data, aes(x = coeff_skew, y = coeff_med, color = noise_sg))
scatterplot + geom_point(size = 6) + labs(x = 'Coeff_Skew', y = 'Coefficient Median') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('skew_med_subj.png')

scatterplot = ggplot(data, aes(x = coeff_skew, y = coeff_min, color = noise_sg))
scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Coefficient Min') + 
   theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('energy_min_subj.png')

scatterplot = ggplot(data, aes(x = coeff_skew, y = coeff_quant, color = noise_sg))
scatterplot + geom_point(size = 6) + labs(x = 'skewness', y = 'Coefficient Quantile') + 
  theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('skew_quant_subj.png')

scatterplot = ggplot(data, aes(x = log(temp_energy), y = coeff_autocorr, color = noise_sg))
scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Autocorrelation') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('energy_ac_subj.png') 






# Relationship between regional expression and BC
scatterplot = ggplot(data, aes(x = lf_wi, y = behave_coeff, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'Mean Within Region Edge', y = 'Behavioral Coefficient') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('lf_wi_bc.png')

scatterplot = ggplot(data, aes(x = v_bw, y = behave_coeff, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('lf_wi_bc.png')

scatterplot = ggplot(data, aes(x = rf_wi, y = behave_coeff, color = cond))
scatterplot + geom_point(size = 6) + labs(x = 'Mean Within Region Edge', y = 'Behavioral Coefficient') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('rf_wi_bc.png')

scatterplot = ggplot(data, aes(x = rf_bw, y = behave_coeff, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('rf_wi_bc.png')

scatterplot = ggplot(data, aes(x = lfv, y = behave_coeff, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('lfv_bc.png')

scatterplot = ggplot(data, aes(x = rfv, y = behave_coeff, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('rfv_bc.png')





# Data without noise SG
data_clean = filter(data, noise_sg == FALSE)
scatterplot = ggplot(data_clean, aes(x = log(temp_energy), y = coeff_skew, color = subj))
scatterplot + geom_point(size = 6) + labs(x = 'log(Temporal Energy)', y = 'Coefficient Skew') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('energy_skew_subj_clean.png') 


scatterplot = ggplot(data_clean, aes(x = behave_coeff, y = coeff_med, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'Behavioral Coeff', y = 'Coefficient Med') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
ggsave('energy_bc_band_clean.png') 

scatterplot = ggplot(data_clean, aes(x = rf_bw, y = behave_coeff, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'Mean Between Region Edge', y = 'Behavioral Coefficient') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  wes_palette("Moonrise1",4))
