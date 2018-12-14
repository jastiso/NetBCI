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

# load data
max = readMat(paste('data/wpli/', sens, '/max_exp_b_corr.mat', sep = ''))
min = readMat(paste('data/wpli/', sens, '/min_exp_b_corr.mat', sep = ''))


data = data.frame(BE = character(length = 0), cond = character(length = 0), band = character(length = 0), corr = numeric(0), stringsAsFactors = FALSE)
for (b in c(1:length(bands))){
  #max
  data[(nrow(data)+1):(nrow(data)+nSubj),1] = rep('max', times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),2] = rep('task', times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),3] = rep(bands[b], times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),4] = max$max.corrs[,b]
  
  #min
  data[(nrow(data)+1):(nrow(data)+nSubj),1] = rep('min', times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),2] = rep('task', times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),3] = rep(bands[b], times = nSubj)
  data[(nrow(data)-nSubj+1):nrow(data),4] = min$min.corrs[,b]
  
  
}

plot = ggplot(data, aes(x = BE, y = corr, fill = band))
plot + geom_boxplot(notch = FALSE, lwd = 1) +
  scale_fill_manual(values = brewer.pal(4,'Greys')) + 
  labs(x = 'Subgraph', y = 'Corr')  + theme_minimal()
ggsave(paste('b_exp_corr', '.png', sep = ''))

# stats
# anova
fit1 = lmp(corr ~ BE + band, data = data)
summary(fit1)
Anova(fit1)
fit2 = lm(corr ~ band+BE, data = data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(corr ~ BE, data = data)
summary(stat1)
# for each condition
max = filter(data, BE == "max")
min = filter(data, BE == 'min')
stat = lmp(corr ~ BE + band, data = bl)
summary(stat)
stat = lmp(corr ~ BE + band, data = task)
Anova(stat)
# specififc ttest
stat = permTS(min$corr, max$corr)
stat
stat = permTS(min$corr)
stat



#############################################################################################

# effects with slope

#############################################################################################
beh = readMat(paste('data/wpli/', sens, '/exp_beh_cor.mat', sep = ''))

corr_data = data.frame(band = character(length = 0), n_zero = numeric(0), max = numeric(0), max2 = numeric(0), max3 = numeric(0), max4 = numeric(0), 
                       min = numeric(0), min2 = numeric(0), sd = numeric(0), sum = numeric(0), slope = numeric(0), 
                       fin = numeric(0), maximum = numeric(0), subj = numeric(0), diff_exp = numeric(0), stringsAsFactors = FALSE)
for (b in c(1:length(bands))){
  
  corr_data[(nrow(corr_data)+1):(nrow(corr_data)+nSubj),1] = rep(bands[b], times = nSubj)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),2] = beh$num.zero[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),3] = beh$max.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),4] = beh$max.exp2[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),5] = beh$max.exp3[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),6] = beh$max.exp4[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),7] = beh$min.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),8] = beh$min.exp2[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),9] = beh$sd.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),10] = beh$sum.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),11] = t(beh$betas)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),12] = t(beh$final)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),13] = t(beh$maximum)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),14] = t(beh$Subj)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),15] = beh$diff.exp[,b]
  
}
corr_data$subj = as.factor(corr_data$subj)

# scatter plot for all bands - mean
plot = ggplot(corr_data, aes(x = slope, y = sum, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max BC')  + theme_minimal()
ggsave(paste('sum_bc_corr.png', sep = ''))

# scatter plot for all bands - std
plot = ggplot(corr_data, aes(x = slope, y = sd, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max BC')  + theme_minimal()
ggsave(paste('sd_bc_corr.png', sep = ''))

# scatter plot for all bands - diff
plot = ggplot(corr_data, aes(x = slope, y = diff_exp, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max BC')  + theme_minimal()
ggsave(paste('diff_bc_corr.png', sep = ''))

# scatter plot for all bands
plot = ggplot(corr_data, aes(x = slope, y = max, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max BC')  + theme_minimal()
ggsave(paste('max_bc_corr.png', sep = ''))

plot = ggplot(corr_data, aes(x = slope, y = sum, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Mean BC')  + theme_minimal()
ggsave(paste('sum_bc_corr.png', sep = ''))

plot = ggplot(corr_data, aes(x = slope, y = sum, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Std BC')  + theme_minimal()
ggsave(paste('std_bc_corr.png', sep = ''))

# scatter plot for all bands - second highest
plot = ggplot(corr_data, aes(x = slope, y = max2, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Max2 BC')  + theme_minimal()
ggsave(paste('max_bc_corr2.png', sep = ''))

# scatter plot for all bands - third highest
plot = ggplot(corr_data, aes(x = slope, y = max3, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
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
ggsave(paste('min_bc_corr.png', sep = ''))

# scatter plot for all bands - 2nd smallest nonzero
plot = ggplot(corr_data, aes(x = slope, y = min2, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'Min 2 BC')  + theme_minimal()
ggsave(paste('min2_bc_corr.png', sep = ''))

# scatter plot for all bands - number of 0 graphs
plot = ggplot(corr_data, aes(x = slope, y = n_zero, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = wes_palette('Royal1',4)) +
  labs(x = 'Slope', y = 'N Zero')  + theme_minimal()
ggsave(paste('nzero_corr.png', sep = ''))

# do you want to include gamma?
# corr_data = filter(corr_data, band != "low_gamma")

# stats

# mean
fit1 = lm(slope ~ sum + band, data=  corr_data)
summary(fit1)
Anova(fit1)

fit_perm = lmp(slope ~ sum + band, data=  corr_data)
summary(fit_perm)
anova(fit_perm)

# sd
fit1 = lm(slope ~ sd + band, data=  corr_data)
summary(fit1)
Anova(fit1)

fit_perm = lmp(slope ~ sd + band, data=  corr_data)
summary(fit_perm)
anova(fit_perm)

# max
fit1 = lm(slope ~ max + band, data=  corr_data)
summary(fit1)
Anova(fit1)

fit1 = lmer(slope ~ max + band + (1 | subj), data=  corr_data)
summary(fit1)
Anova(fit1)

fit_perm = lmp(slope ~ max + band, data=  corr_data)
summary(fit_perm)
anova(fit_perm)

# second max
fit2 = lm(slope ~ max2 + band, data=  corr_data)
summary(fit2)
Anova(fit2)

fit2_perm = lmp(slope ~ max2 + band, data=  corr_data)
summary(fit2_perm)
anova(fit2_perm)

# third max
fit3 = lm(slope ~ max3 + band, data=  corr_data)
summary(fit3)
Anova(fit3)

fit3_perm = lmp(slope ~ max3 + band, data=  corr_data)
summary(fit3_perm)
anova(fit3_perm)

# 4 max
fit4 = lm(slope ~ max4 + band, data=  corr_data)
summary(fit4)
Anova(fit4)

fit4_perm = lmp(slope ~ max4 + band, data=  corr_data)
summary(fit4_perm)
anova(fit4_perm)

# min
fit_min = lm(slope ~ min + band, data=  corr_data)
summary(fit_min)
Anova(fit_min)

fit_min_perm = lmp(slope ~ min + band, data=  corr_data)
summary(fit_min_perm)
anova(fit_min_perm)


# min
fit_min2 = lm(slope ~ min2 + band, data=  corr_data)
summary(fit_min2)
Anova(fit_min2)

fit_min_perm2 = lmp(slope ~ min2 + band, data=  corr_data)
summary(fit_min_perm2)
anova(fit_min_perm2)





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


# scatter plot for all bands
plot = ggplot(corr_data_pr, aes(x = slope, y = max, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = brewer.pal(4,'Greys')) +
  labs(x = 'Slope', y = 'Max BC')  + theme_minimal()
ggsave(paste('max_bc_corr_pr.png', sep = ''))

# scatter plot for all bands - second highest
plot = ggplot(corr_data_pr, aes(x = slope, y = max2, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = brewer.pal(4,'Greys')) +
  labs(x = 'Slope', y = 'Max2 BC')  + theme_minimal()
ggsave(paste('max_bc_corr2_pr.png', sep = ''))

# scatter plot for all bands - third highest
plot = ggplot(corr_data_pr, aes(x = slope, y = max3, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = brewer.pal(4,'Greys')) +
  labs(x = 'Slope', y = 'Max3 BC')  + theme_minimal()
ggsave(paste('max_bc_corr3_pr.png', sep = ''))

# scatter plot for all bands - smallest nonzero
plot = ggplot(corr_data_pr, aes(x = slope, y = min, col = band, group = band))
plot +geom_smooth(method="lm") +  geom_point(size = 6) + 
  scale_color_manual(values = brewer.pal(4,'Greys')) +
  labs(x = 'Slope', y = 'Min BC')  + theme_minimal()
ggsave(paste('min_bc_corr_pr.png', sep = ''))

# stats
# max
fit1 = lm(slope ~ max + band, data=  corr_data_pr)
summary(fit1)
Anova(fit1)

fit_perm = lmp(slope ~ max + band, data=  corr_data_pr)
summary(fit_perm)
anova(fit_perm)

# second max
fit2 = lm(slope ~ max2 + band, data=  corr_data_pr)
summary(fit2)
Anova(fit2)

fit2_perm = lmp(slope ~ max2 + band, data=  corr_data_pr)
summary(fit2_perm)
anova(fit2_perm)

# third max
fit3 = lm(slope ~ max3 + band, data=  corr_data_pr)
summary(fit3)
Anova(fit3)

fit3_perm = lmp(slope ~ max3 + band, data=  corr_data_pr)
summary(fit3_perm)
anova(fit3_perm)

# min
fit4 = lm(slope ~ min + band, data=  corr_data_pr)
summary(fit4)
Anova(fit4)

fit4_perm = lmp(slope ~ min + band, data=  corr_data_pr)
summary(fit4_perm)
anova(fit4_perm)

# mean
fit1 = lm(slope ~ sum + band, data=  corr_data_pr)
summary(fit1)
Anova(fit1)

fit_perm = lmp(slope ~ sum + band, data=  corr_data_pr)
summary(fit_perm)
anova(fit_perm)

# sd
fit1 = lm(slope ~ sd + band, data=  corr_data_pr)
summary(fit1)
Anova(fit1)

fit_perm = lmp(slope ~ sd + band, data=  corr_data_pr)
summary(fit_perm)
anova(fit_perm)
