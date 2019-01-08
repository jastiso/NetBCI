# Behavior and BE correlations

library(GGally)
library(ggplot2)
library(dplyr)
library(coin)
library(lmPerm)
library(car)
library(aplpack)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
library(coin)
library(perm)
library(lme4)
setwd("/Users/stiso/Documents/R/NetBCI/")
bands = c( 'alpha', 'beta', 'low_gamma')
nSubj = 20
sens = 'grad'


#############################################################################################

# effects with learning rate

#############################################################################################
beh = readMat(paste('data/wpli/', sens, '/exp_beh_cor.mat', sep = ''))

corr_data = data.frame(band = character(length = 0), n_zero = numeric(0), max = numeric(0), max2 = numeric(0), max3 = numeric(0), max4 = numeric(0), 
                       min = numeric(0), sd = numeric(0), sum = numeric(0), slope = numeric(0), subj = numeric(0), stringsAsFactors = FALSE)
for (b in c(1:length(bands))){

  corr_data[(nrow(corr_data)+1):(nrow(corr_data)+nSubj),1] = rep(bands[b], times = nSubj)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),2] = beh$num.zero[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),3] = beh$max.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),4] = beh$max.exp2[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),5] = beh$max.exp3[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),6] = beh$max.exp4[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),7] = beh$min.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),8] = beh$sd.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),9] = beh$sum.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),10] = t(beh$betas)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),11] = t(beh$Subj)

}
corr_data$subj = as.factor(corr_data$subj)
# save p_values
stats = data.frame(p_val = rep(0, times = 12), b_val = rep(0, times = 12), pred = c('mean', 'sd', 'high', 'high2', 'high3', 'low', 
                                                                                    'mean', 'sd', 'high', 'high2', 'high3', 'low'), 
                   model = c(rep('emp', times = 6), rep('UPR', times = 6)))


# scatter plot for all bands - mean
plot = ggplot(corr_data, aes(x = slope, y = sum, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Mean')  + theme_minimal()
ggsave(paste('sum_bc_corr.png', sep = ''))

# scatter plot for all bands - std
plot = ggplot(corr_data, aes(x = slope, y = sd, color = band))
plot  +  geom_point(aes(color = band), size = 4) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  geom_smooth(method="lm", size = 3, color = "black") +
  labs(x = 'Slope', y = 'Standard Deviation')  + theme_minimal()
ggsave(paste('sd_bc_corr.png', sep = ''))


# scatter plot for all bands - individual subgraphs
plot = ggplot(corr_data, aes(x = slope, y = max, col = band))
plot +geom_smooth(method="lm", color = 'black',size = 3) +  geom_point(size = 4) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max BC')  + theme_minimal()
ggsave(paste('max_bc_corr.png', sep = ''))


# scatter plot for all bands - second highest
plot = ggplot(corr_data, aes(x = slope, y = max2, col = band))
plot +geom_smooth(method="lm", color = "black", size = 3) +  geom_point(size = 4) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max2 BC')  + theme_minimal()
ggsave(paste('max_bc_corr2.png', sep = ''))

# scatter plot for all bands - third highest
plot = ggplot(corr_data, aes(x = slope, y = max3, col = band))
plot +geom_smooth(method="lm", color = "black", size = 3) +  geom_point(size = 4) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max3 BC')  + theme_minimal()
ggsave(paste('max_bc_corr3.png', sep = ''))

# scatter plot for all bands - 4th highest
plot = ggplot(corr_data, aes(x = slope, y = max4, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max4 BC')  + theme_minimal()
ggsave(paste('max_bc_corr4.png', sep = ''))

# scatter plot for all bands - smallest nonzero
plot = ggplot(corr_data, aes(x = slope, y = min, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max BC')  + theme_minimal()
ggsave(paste('min_bc_corr.pdf', sep = ''))

# scatter plot for all bands - number of 0 graphs
plot = ggplot(corr_data, aes(x = slope, y = n_zero, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'N Zero')  + theme_minimal()
ggsave(paste('nzero_corr.png', sep = ''))


# stats
# mean
fit1 = lm(slope ~ sum + band, data=  corr_data)
summary(fit1)
Anova(fit1)

stats$b_val[1] = fit1$coefficients[2]
stats$p_val[1] = summary(fit1)$coefficients[2,4]

# sd
fit2 = lm(slope ~ sd + band, data=  corr_data)
summary(fit2)
Anova(fit1)

stats$b_val[2] = fit2$coefficients[2]
stats$p_val[2] = summary(fit2)$coefficients[2,4]

# max
fit3 = lm(slope ~ max + band, data=  corr_data)
summary(fit3)
Anova(fit3)

stats$b_val[3] = fit3$coefficients[2]
stats$p_val[3] = summary(fit3)$coefficients[2,4]

# second max
fit4 = lm(slope ~ max2 + band, data=  corr_data)
summary(fit4)
Anova(fit4)

stats$b_val[4] = fit4$coefficients[2]
stats$p_val[4] = summary(fit4)$coefficients[2,4]

# third max
fit5 = lm(slope ~ max3 + band, data=  corr_data)
summary(fit5)
Anova(fit5)

stats$b_val[5] = fit5$coefficients[2]
stats$p_val[5] = summary(fit5)$coefficients[2,4]


# min
fit6 = lm(slope ~ min + band, data=  corr_data)
summary(fit6)
Anova(fit6)

stats$b_val[6] = fit6$coefficients[2]
stats$p_val[6] = summary(fit)$coefficients[2,4]




#############################################################################################

# effects with slope - PR

#############################################################################################
beh_pr = readMat(paste('data/wpli/', sens, '/exp_beh_cor_pr.mat', sep = ''))

corr_data_pr = data.frame(band = character(length = 0), max = numeric(0), max2 = numeric(0), max3 = numeric(0), min = numeric(0), sd = numeric(0), 
                          sum = numeric(0), slope = numeric(0), subj = numeric(0), stringsAsFactors = FALSE)
for (b in c(1:length(bands))){
  
  corr_data_pr[(nrow(corr_data_pr)+1):(nrow(corr_data_pr)+nSubj),1] = rep(bands[b], times = nSubj)
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),2] = beh_pr$max.exp.null[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),3] = beh_pr$max.exp.null2[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),4] = beh_pr$max.exp.null3[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),5] = beh_pr$min.exp.null[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),6] = beh_pr$sd.exp.null[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),7] = beh_pr$sum.exp.null[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),8] = t(beh_pr$betas)
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),9] = t(beh_pr$Subj)
  
}
corr_data_pr$subj = as.factor(corr_data_pr$subj)

# scatter plot for all bands - mean
plot = ggplot(corr_data_pr, aes(x = slope, y = sum, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = rev(brewer.pal(4,'Greys'))) +
  labs(x = 'Slope', y = 'Mean BC')  + theme_minimal()
ggsave(paste('mean_bc_corr3_pr.png', sep = ''))

# scatter plot for all bands - std
plot = ggplot(corr_data_pr, aes(x = slope, y = sd, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = rev(brewer.pal(4,'Greys'))) +
  labs(x = 'Slope', y = 'SD BC')  + theme_minimal()
ggsave(paste('std_bc_corr_pr.pdf', sep = ''))

# scatter plot for all bands - max
plot = ggplot(corr_data_pr, aes(x = slope, y = max, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = rev(brewer.pal(4,'Greys'))) +
  labs(x = 'Slope', y = 'Max BC')  + theme_minimal()
ggsave(paste('max_bc_corr_pr.pdf', sep = ''))

# scatter plot for all bands - second highest
plot = ggplot(corr_data_pr, aes(x = slope, y = max2, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = rev(brewer.pal(4,'Greys'))) +
  labs(x = 'Slope', y = 'Max2 BC')  + theme_minimal()
ggsave(paste('max_bc_corr2_pr.pdf', sep = ''))

# scatter plot for all bands - third highest
plot = ggplot(corr_data_pr, aes(x = slope, y = max3, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = rev(brewer.pal(4,'Greys'))) +
  labs(x = 'Slope', y = 'Max3 BC')  + theme_minimal()
ggsave(paste('max_bc_corr3_pr.pdf', sep = ''))

# scatter plot for all bands - smallest nonzero
plot = ggplot(corr_data_pr, aes(x = slope, y = min, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = rev(brewer.pal(4,'Greys'))) +
  labs(x = 'Slope', y = 'Min BC')  + theme_minimal()
ggsave(paste('min_bc_corr_pr.pdf', sep = ''))

# stats
# mean
fit7 = lm(slope ~ sum + band, data=  corr_data_pr)
summary(fit7)
Anova(fit7)

stats$b_val[7] = fit7$coefficients[2]
stats$p_val[7] = summary(fit7)$coefficients[2,4]

# sd
fit8 = lm(slope ~ sd + band, data=  corr_data_pr)
summary(fit8)
Anova(fit8)

stats$b_val[8] = fit8$coefficients[2]
stats$p_val[8] = summary(fit8)$coefficients[2,4]

#  max
fit9 = lm(slope ~ max + band, data=  corr_data_pr)
summary(fit9)
Anova(fit9)

stats$b_val[9] = fit9$coefficients[2]
stats$p_val[9] = summary(fit9)$coefficients[2,4]

# 2nd max
fit10 = lm(slope ~ max2 + band, data=  corr_data_pr)
summary(fit10)
Anova(fit10)

stats$b_val[10] = fit10$coefficients[2]
stats$p_val[10] = summary(fit10)$coefficients[2,4]

# third
fit11 = lm(slope ~ max3 + band, data=  corr_data_pr)
summary(fit11)
Anova(fit11)

stats$b_val[11] = fit11$coefficients[2]
stats$p_val[11] = summary(fit11)$coefficients[2,4]

# min
fit12 = lm(slope ~ min + band, data=  corr_data_pr)
summary(fit12)
Anova(fit12)

stats$b_val[12] = fit12$coefficients[2]
stats$p_val[12] = summary(fit12)$coefficients[2,4]





#####################################

# Combined Data for Plots

####################################

corr_data$band = lapply(corr_data$band, paste, '_emp', sep = '')
corr_data_pr$band = lapply(corr_data_pr$band, paste, '_upr', sep = '')

data_cmb = data.frame(max3 = c(corr_data$max3, corr_data_pr$max3), max2 = c(corr_data$max2, corr_data_pr$max2), sd = c(corr_data$sd, corr_data_pr$sd), 
                      slope = c(corr_data$slope, corr_data_pr$slope), band = c(unlist(corr_data$band), unlist(corr_data_pr$band)), 
                      model = c(rep('emp', times = 60), rep('upr', times = 60)))
data_cmb$band = as.factor(data_cmb$band)
levels(data_cmb$band)
# emp, upr, lines
#tmp = c(wes_palette('Royal1',3), (brewer.pal(3,'Greys')), 'black', 'white')
tmp = c(rgb(81/255,184/255,161/255), rgb(81/255,184/255,161/255), rgb(81/255,184/255,161/255), 'grey', 'grey', 'grey', rgb(81/255,184/255,161/255), 'white')
# alpha & beta band,line_emp, gamma band,  line_upr
color_vect = c(tmp[1], tmp[4], tmp[2], tmp[5], tmp[7], tmp[3], tmp[7], tmp[8])
color_vect

plot = ggplot(data_cmb, aes(x = slope, y = max3, color = band, group = model))
plot  +  geom_point(size = 4) + 
  scale_color_manual(values = color_vect) +
  #scale_fill_manual(values = c('blue', 'black')) +
  geom_smooth(aes(color = model), method="lm", size = 3) +
  labs(x = 'Slope', y = 'Max3')  + theme_minimal() 
ggsave(paste('max3_bc_corr_cmb.png', sep = ''))

plot = ggplot(data_cmb, aes(x = slope, y = max2, color = band, group = model))
plot  +  geom_point(size = 4) + 
  scale_color_manual(values = color_vect) +
  #scale_fill_manual(values = c('blue', 'black')) +
  geom_smooth(aes(color = model), method="lm", size = 3) +
  labs(x = 'Slope', y = 'Max2')  + theme_minimal() 
ggsave(paste('max2_bc_corr_cmb.png', sep = ''))

plot = ggplot(data_cmb, aes(x = slope, y = sd, color = band, group = model))
plot  +  geom_point(size = 4) + 
  scale_color_manual(values = color_vect) +
  #scale_fill_manual(values = c('blue', 'black')) +
  geom_smooth(aes(color = model), method="lm", size = 3) +
  labs(x = 'Slope', y = 'sd')  + theme_minimal() 
ggsave(paste('sd_bc_corr_cmb.png', sep = ''))





# bar plots
plot = ggplot(stats, aes(x = pred, y = p_val, fill = model, group = model))
plot + geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c(rgb(81/255,184/255,161/255), 'grey')) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.008, linetype = "dashed", color = "red") +
  labs(x = 'Predictor', y = 'P-value')  + theme_minimal() 
  ggsave(paste('pred_cmb_p.pdf', sep = ''))
  
  plot = ggplot(stats, aes(x = pred, y = b_val, fill = model, group = model))
  plot + geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c(rgb(81/255,184/255,161/255), 'grey')) +
    labs(x = 'Predictor', y = 'Coefficient')  + theme_minimal() 
  ggsave(paste('pred_cmb_b.pdf', sep = ''))
