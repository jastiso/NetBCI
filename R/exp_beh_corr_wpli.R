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
library(lm.beta)
setwd("/Users/stiso/Documents/R/NetBCI/")
bands = c( 'alpha', 'beta', 'low_gamma')
nSubj = 20
sens = 'grad'


#############################################################################################

# effects with learning rate

#############################################################################################
beh = readMat(paste('data/wpli/', sens, '/exp_beh_cor.mat', sep = ''))

corr_data = data.frame(band = character(length = 0), max = numeric(0), max2 = numeric(0), max3 = numeric(0), max4 = numeric(0), 
                       min = numeric(0), min2 = numeric(0), slope = numeric(0), subj = numeric(0), stringsAsFactors = FALSE)
for (b in c(1:length(bands))){

  corr_data[(nrow(corr_data)+1):(nrow(corr_data)+nSubj),1] = rep(bands[b], times = nSubj)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),2] = beh$max.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),3] = beh$max.exp2[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),4] = beh$max.exp3[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),5] = beh$max.exp4[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),6] = beh$min.exp[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),7] = beh$min.exp2[,b]
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),8] = t(beh$betas)
  corr_data[(nrow(corr_data)-nSubj+1):nrow(corr_data),9] = t(beh$Subj)

}
corr_data$subj = as.factor(corr_data$subj)
# save p_values
stats = data.frame(p_val = rep(0, times = 12), b_val = rep(0, times = 12), pred = c('high', 'high2', 'high3', 'high4', 'low', 'low2', 
                                                                                    'high', 'high2', 'high3', 'high4', 'low', 'low2'), 
                   model = c(rep('emp', times = 6), rep('UPR', times = 6)))



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




# stats

# max
fit1 = lmp(slope ~ max + band, data=  corr_data)
summary(fit1)

stats$b_val[1] = fit1$coefficients[2]
stats$p_val[1] = summary(fit1)$coefficients[2,3]

# second max
fit2 = lmp(slope ~ max2 + band, data=  corr_data)
summary(fit2)

stats$b_val[2] = fit2$coefficients[2]
stats$p_val[2] = summary(fit2)$coefficients[2,3]

# third max
fit3 = lmp(slope ~ max3 + band, data=  corr_data)
summary(fit3)

stats$b_val[3] = fit3$coefficients[2]
stats$p_val[3] = summary(fit3)$coefficients[2,3]

# 4th max
fit4 = lmp(slope ~ max4 + band, data=  corr_data)
summary(fit4)

stats$b_val[4] = fit4$coefficients[2]
stats$p_val[4] = summary(fit4)$coefficients[2,3]


# min
fit5 = lmp(slope ~ min + band, data=  corr_data)
summary(fit5)

stats$b_val[5] = fit5$coefficients[2]
stats$p_val[5] = summary(fit5)$coefficients[2,3]

# 2nd min
fit6 = lmp(slope ~ min2 + band, data=  corr_data)
summary(fit6)

stats$b_val[6] = fit6$coefficients[2]
stats$p_val[6] = summary(fit6)$coefficients[2,3]






#############################################################################################

# effects with slope - PR

#############################################################################################
beh_pr = readMat(paste('data/wpli/', sens, '/exp_beh_cor_pr.mat', sep = ''))

corr_data_pr = data.frame(band = character(length = 0), max = numeric(0), max2 = numeric(0), max3 = numeric(0), max4 = numeric(0), min = numeric(0), min2 = numeric(0), 
                          slope = numeric(0), subj = numeric(0), stringsAsFactors = FALSE)
for (b in c(1:length(bands))){
  
  corr_data_pr[(nrow(corr_data_pr)+1):(nrow(corr_data_pr)+nSubj),1] = rep(bands[b], times = nSubj)
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),2] = beh_pr$max.exp.null[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),3] = beh_pr$max.exp.null2[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),4] = beh_pr$max.exp.null3[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),5] = beh_pr$max.exp.null4[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),6] = beh_pr$min.exp.null[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),7] = beh_pr$min.exp.null2[,b]
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),8] = t(beh_pr$betas)
  corr_data_pr[(nrow(corr_data_pr)-nSubj+1):nrow(corr_data_pr),9] = t(beh_pr$Subj)
  
}
corr_data_pr$subj = as.factor(corr_data_pr$subj)



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
#  max
fit7 = lmp(slope ~ max + band, data=  corr_data_pr)
summary(fit7)

stats$b_val[7] = fit7$coefficients[2]
stats$p_val[7] = summary(fit7)$coefficients[2,3]

# 2nd max
fit8 = lmp(slope ~ max2 + band, data=  corr_data_pr)
summary(fit8)

stats$b_val[8] = fit8$coefficients[2]
stats$p_val[8] = summary(fit8)$coefficients[2,3]

# third
fit9 = lmp(slope ~ max3 + band, data=  corr_data_pr)
summary(fit9)

stats$b_val[9] = fit9$coefficients[2]
stats$p_val[9] = summary(fit9)$coefficients[2,3]

# fourth
fit10 = lmp(slope ~ max4 + band, data=  corr_data_pr)
summary(fit10)

stats$b_val[10] = fit10$coefficients[2]
stats$p_val[10] = summary(fit10)$coefficients[2,3]

# min
fit11 = lmp(slope ~ min + band, data=  corr_data_pr)
summary(fit11)

stats$b_val[11] = fit11$coefficients[2]
stats$p_val[11] = summary(fit11)$coefficients[2,3]

# min2
fit12 = lmp(slope ~ min2 + band, data=  corr_data_pr)
summary(fit12)

stats$b_val[12] = fit12$coefficients[2]
stats$p_val[12] = summary(fit12)$coefficients[2,3]





#########################
# Bootstapping
#########################

# Bootstrap erro bars
bootdata=list()
nBoot = 500
for (i in 1:nBoot) {
  #curr=sample(seq(1,20), 20, replace=TRUE)
  curr=sample(seq(1,nrow(corr_data)), nrow(corr_data), replace=TRUE)
  curr_data = data.frame()
  for (j in 1:length(curr)){
    #curr_data = rbind(curr_data,dplyr::filter(corr_data, subj == curr[j]))
    curr_data = rbind(curr_data,corr_data[curr[j],])
  }
  bootdata[[i]] = curr_data
}


# Bootstrap for null
bootdata_upr=list()
for (i in 1:nBoot) {
  #curr=sample(seq(1,20), 20, replace=TRUE)
  curr=sample(seq(1,nrow(corr_data_pr)), nrow(corr_data_pr), replace=TRUE)
  curr_data = data.frame(corr_data_pr[0,])
  for (j in 1:length(curr)){
    #curr_data = rbind(curr_data,dplyr::filter(corr_data_pr, subj == curr[j]))
    curr_data = rbind(curr_data,corr_data_pr[curr[j],])
  }
  bootdata_upr[[i]] = curr_data
}

# get betas and pvals
boot_stats = data.frame(p_val = rep(0, times = 12*nBoot), b_val = rep(0, times = 12*nBoot), 
                        pred = rep(c('high', 'high2', 'high3', 'high4', 'low', 'low2', 'high', 'high2', 'high3', 'high4', 'low', 'low2'), times = nBoot), 
                        model = c(rep('emp', times = 6*nBoot), rep('UPR', times = 6*nBoot)))
cnt = 1
for (i in 1:nBoot){
 
  
  #1
  curr_1 = lmp(slope ~ max + band, data=  bootdata[[i]])
  boot_stats$b_val[cnt] = (curr_1)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_1)$coefficients[2,3]
  cnt = cnt + 1
  
  #2
  curr_2 = lmp(slope ~ max2 + band, data=  bootdata[[i]])
  boot_stats$b_val[cnt] = (curr_2)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_2)$coefficients[2,3]
  cnt = cnt + 1
  
  #3
  curr_3 = lmp(slope ~ max3 + band, data=  bootdata[[i]])
  boot_stats$b_val[cnt] = (curr_3)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_3)$coefficients[2,3]
  cnt = cnt + 1
  
  #4
  curr_4 = lmp(slope ~ max4 + band, data=  bootdata[[i]])
  boot_stats$b_val[cnt] = (curr_4)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_4)$coefficients[2,3]
  cnt = cnt + 1
  
  #low
  curr_min = lmp(slope ~ min + band, data=  bootdata[[i]])
  boot_stats$b_val[cnt] = (curr_min)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_min)$coefficients[2,3]
  cnt = cnt + 1
  
  #l2
  curr_min2 = lmp(slope ~ min2 + band, data=  bootdata[[i]])
  boot_stats$b_val[cnt] = (curr_min2)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_min2)$coefficients[2,3]
  cnt = cnt + 1
}
# for for null data
for (i in 1:nBoot){

  #1
  curr_1 = lmp(slope ~ max + band, data=  bootdata_upr[[i]])
  boot_stats$b_val[cnt] = (curr_1)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_1)$coefficients[2,3]
  cnt = cnt + 1
  
  #2
  curr_2 = lmp(slope ~ max2 + band, data=  bootdata_upr[[i]])
  boot_stats$b_val[cnt] = (curr_2)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_2)$coefficients[2,3]
  cnt = cnt + 1
  
  #3
  curr_3 = lmp(slope ~ max3 + band, data=  bootdata_upr[[i]])
  boot_stats$b_val[cnt] = (curr_3)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_3)$coefficients[2,3]
  cnt = cnt + 1
  
  #4
  curr_4 = lmp(slope ~ max4 + band, data=  bootdata_upr[[i]])
  boot_stats$b_val[cnt] = (curr_4)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_4)$coefficients[2,3]
  cnt = cnt + 1
  
  #min
  curr_min = lmp(slope ~ min + band, data=  bootdata_upr[[i]])
  boot_stats$b_val[cnt] = (curr_min)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_min)$coefficients[2,3]
  cnt = cnt + 1
  
  #l2
  curr_min2 = lmp(slope ~ min2 + band, data=  bootdata_upr[[i]])
  boot_stats$b_val[cnt] = (curr_min2)$coefficients[2]
  boot_stats$p_val[cnt] = summary(curr_min2)$coefficients[2,3]
  cnt = cnt + 1
}

# summarize

boot_stats$p_val = log10(boot_stats$p_val)
boot_stats$p_val[boot_stats$p_val == -Inf] = -16
bootstrap_summary = boot_stats %>%
  group_by(pred, model) %>%
  dplyr::summarise(mean_p = median(p_val), sd_p = sd(p_val)/sqrt(nBoot), mean_b = mean(b_val), sd_b = sd(b_val)/sqrt(nBoot))

p<-ggplot(boot_stats, aes(x=log(p_val), fill = pred, color=model)) + 
  geom_histogram()
p




#####################################

# Combined Data for Plots

####################################


# bar plots
stats_rank = merge(stats,bootstrap_summary)



plot = ggplot(stats_rank, aes(x = pred, y = log10(p_val), fill = model, group = model))
plot + geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c(rgb(81/255,184/255,161/255), 'grey')) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed", color = "black") +
  geom_hline(yintercept = log10(0.008), linetype = "dashed", color = "red") +
  geom_errorbar(aes(ymin = mean_p - sd_p, ymax = (mean_p) + (sd_p)), width=0.2, position=position_dodge(.9)) +
  labs(x = 'Predictor', y = 'P-value')  + theme_minimal() 
  ggsave(paste('pred_cmb_p.pdf', sep = ''))
  
plot = ggplot(stats_rank, aes(x = pred, y = b_val, fill = model, group = model))
  plot + geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c(rgb(81/255,184/255,161/255), 'grey')) +
    geom_errorbar(aes(ymin = mean_b - sd_b, ymax = mean_b + sd_b), width=0.2, position=position_dodge(.9)) +
    labs(x = 'Predictor', y = 'Coefficient')  + theme_minimal() 
  ggsave(paste('pred_cmb_b.pdf', sep = ''))

  
  
##################
  # Normality test
##################


qqPlot(unique(corr_data$slope))

par(mfrow=c(3,2))
qqPlot(unique(corr_data$max))
qqPlot(unique(corr_data$max2))
qqPlot(unique(corr_data$max3))
qqPlot(unique(corr_data$max4))
qqPlot(unique(corr_data$min))
qqPlot(unique(corr_data$min2))

par(mfrow=c(3,2))
qqPlot(unique(corr_data_pr$max))
qqPlot(unique(corr_data_pr$max2))
qqPlot(unique(corr_data_pr$max3))
qqPlot(unique(corr_data_pr$max4))
qqPlot(unique(corr_data_pr$min))
qqPlot(unique(corr_data_pr$min2))
  
  
  