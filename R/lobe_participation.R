# Regional expression

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
bands = c('alpha', 'beta', 'low_gamma','gamma')
nSubj = 20
sens = 'grad'

h = readMat(paste('data/gc/', sens, '/high_participation.mat', sep = ''))
l = readMat(paste('data/gc/', sens, '/low_participation.mat', sep = ''))

high = data.frame(subj = character(len = length(h$band.h)), band = character(len = length(h$band.h)), region = character(len = length(h$band.h)), part = numeric(length(h$band.h)), slope = numeric(length(h$band.h)))
high$band = unlist(h$band.h)
high$subj = unlist(h$subj.h)
high$region = unlist(h$region.h)
high$part = unlist(h$part.h)
high$slope = unlist(h$slope)

low = data.frame(subj = character(len = length(h$band.h)), band = character(len = length(h$band.h)), region = character(len = length(h$band.h)), part = numeric(length(h$band.h)), slope = numeric(length(h$band.h)))
low$band = unlist(l$band.l)
low$subj = unlist(l$subj.l)
low$region = unlist(l$region.l)
low$part = unlist(l$part.l)
low$slope = unlist(l$slope)

# plot
plot = ggplot(high, aes(x = region, y = part, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = rev(brewer.pal(4,'Blues'))) + 
  labs(x = 'Region', y = 'Participation')  + theme_minimal()
ggsave(paste(sens, '_high_participation', '.png', sep = ''))

# plot
plot = ggplot(low, aes(x = region, y = part, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = rev(brewer.pal(4,'Greens'))) + 
  labs(x = 'Region', y = 'Exp')  + theme_minimal()
ggsave(paste(sens, '_low_participation', '.png', sep = ''))


#####################

# Time Shift

######################

h_ts = readMat(paste('data/gc/', sens, '/ts_high_participation.mat', sep = ''))
l_ts = readMat(paste('data/gc/', sens, '/ts_low_participation.mat', sep = ''))

high_ts = data.frame(subj = character(len = length(h_ts$band.h)), band = character(len = length(h_ts$band.h)), region = character(len = length(h_ts$band.h)), part = numeric(length(h_ts$band.h)), slope = numeric(length(h$band.h)))
high_ts$band = unlist(h_ts$band.h)
high_ts$subj = unlist(h_ts$subj.h)
high_ts$region = unlist(h_ts$region.h)
high_ts$part = unlist(h_ts$part.h)
high_ts$slope = unlist(h_ts$slope)

low_ts = data.frame(subj = character(len = length(h_ts$band.h)), band = character(len = length(h_ts$band.h)), region = character(len = length(h_ts$band.h)), part = numeric(length(h_ts$band.h)), slope = numeric(length(h$band.h)))
low_ts$band = unlist(l_ts$band.l)
low_ts$subj = unlist(l_ts$subj.l)
low_ts$region = unlist(l_ts$region.l)
low_ts$part = unlist(l_ts$part.l)
low_ts$slope = unlist(l_ts$slope)

# plot
plot = ggplot(high_ts, aes(x = region, y = (part), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = rev(brewer.pal(3,'Greys'))) + 
  labs(x = 'Region', y = 'Exp')  + theme_minimal()
ggsave(paste(sens, '_high_ts_part', '.png', sep = ''))

# plot
plot = ggplot(low_ts, aes(x = region, y = (part), fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = rev(brewer.pal(3,'Greys'))) + 
  labs(x = 'Region', y = 'Exp')  + theme_minimal()
ggsave(paste(sens, '_low_ts_part', '.png', sep = ''))



############## Correlation
stats_exp = data.frame(region = character(length = 0), band = character(length = 0), r = numeric(0), p = numeric(0), stringsAsFactors = FALSE)
stats_exp_low = data.frame(region = character(length = 0), band = character(length = 0), r = numeric(0), p = numeric(0), stringsAsFactors = FALSE)

for (i in unique(high$region)){
  for (j in unique(high$band)){
    # high
    curr = filter(high, region == i, band == j)
    plot = ggplot(curr, aes(x = (part), y = slope))  
    plot + geom_point(size = 6, color = "slateblue") + labs(x = 'Participation', y = 'Slope') + 
      geom_smooth(method="lm", color = 'black') + theme_minimal()
    ggsave(paste(sens, '_high_part_', i, '_', j, '.png', sep = ''))
    
    curr_stat = cor.test(curr$slope, (curr$part), method = 'spearman')
    stats_exp[nrow(stats_exp)+1,1] = i
    stats_exp[nrow(stats_exp),2] = j
    stats_exp[nrow(stats_exp),3] = curr_stat$estimate
    stats_exp[nrow(stats_exp),4] = curr_stat$p.value
    
    # low
    curr_l = filter(low, region == i, band == j)
    plot = ggplot(curr_l, aes(x = (part), y = slope))  
    plot + geom_point(size = 6, color = "thistle") + labs(x = 'Participation', y = 'Slope') + 
      geom_smooth(method="lm", color = 'black') + theme_minimal() 
    ggsave(paste(sens, '_low_part_', i, '_', j, '.png', sep = ''))
    
    curr_stat_l = cor.test(curr_l$slope, (curr_l$part), method = 'spearman')
    stats_exp_low[nrow(stats_exp_low)+1,1] = i
    stats_exp_low[nrow(stats_exp_low),2] = j
    stats_exp_low[nrow(stats_exp_low),3] = curr_stat_l$estimate
    stats_exp_low[nrow(stats_exp_low),4] = curr_stat_l$p.value
  }
}



# stats
stat = lm(slope ~ part*band*region, data=high)
Anova(stat)
summary(stat)

stat = lm(slope ~ part*band*region, data=low)
Anova(stat)
summary(stat)

#################### 

#combined stats

####################
all = high
all$model = rep('aligned', times = nrow(high))
high_ts$model = rep('shifted', times = nrow(high_ts))
all[(nrow(high)+1):(nrow(high)+nrow(high_ts)),] = high_ts
stat = lmp(part ~ region + model + band, data=all)
Anova(stat)
anova(stat)
summary(stat)
perm = ezPerm(data = all, dv = log(exp), wid = .(subj), within = NULL, between = .(model, region, band))
perm


all_low = low
all_low$model = rep('aligned', times = nrow(low))
low_ts$model = rep('shifted', times = nrow(low_ts))
all_low[(nrow(low)+1):(nrow(low)+nrow(low_ts)),] = low_ts
stat = lmp(log(exp) ~ region + model + band, data=all_low)
Anova(stat)
anova(stat)
perm = ezPerm(data = all_low, dv = log(exp), wid = .(subj), within = NULL, between = .(model, region, band))
perm

###############################

# Compare for Band

################################
curr_band = 'alpha'
dodge <- position_dodge(width = 0.4)

beta_H = high
beta_H$model = character(dim(high)[1])
beta_H$model = rep('aligned', times = dim(high)[1])
high_ts$model = rep('shifted', times = dim(high)[1])
beta_H[(nrow(beta_H)+1):(nrow(beta_H)+dim(high_ts)[1]),] = high_ts
beta_H = filter(beta_H, band == curr_band)

#plot
plot = ggplot(beta_H, aes(x = region, y = part, fill = model) )
plot + geom_violin(position = dodge, trim = FALSE) + geom_boxplot(position = dodge, width = 0.2, outlier.color = NA) +
  #geom_dotplot(binaxis='y', stackdir='center', stackgroups=TRUE, dotsize=.3)
  scale_fill_manual(values = (brewer.pal(3,'YlGn'))) + 
  labs(x = 'Region', y = 'Exp')  + theme_minimal()
ggsave(paste(sens, '_', curr_band, '_part', '.png', sep = ''))

stat1 = lmp(part ~ as.factor(region) + as.factor(model), data = beta_H)
stat2 = lmp(part ~ region, data = beta_H)
summary(stat1)
Anova(stat1)
anova(stat1)
anova(stat1, stat2, test='Chisq')


beta_L = low
beta_L$model = character(dim(high)[1])
beta_L$model = rep('aligned', times = dim(high)[1])
low_ts$model = rep('shifted', times = dim(high)[1])
beta_L[(nrow(beta_L)+1):(nrow(beta_L)+dim(high_ts)[1]),] = low_ts
beta_L = filter(beta_L, band == curr_band)

#plot
plot = ggplot(beta_L, aes(x = region, y = part, fill = model) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = (brewer.pal(3,'YlGn'))) + 
  labs(x = 'Region', y = 'Exp')  + theme_minimal()
ggsave(paste(sens, '_', curr_band, '_L_part', '.png', sep = ''))

beta_L$subj = as.factor(beta_L$subj)
beta_L$model = as.factor(beta_L$model)
beta_L$region = as.factor(beta_L$region)
beta_L$exp = as.numeric(beta_L$exp)

stat1 = lmp(log(exp) ~ region + model, data = beta_L)
stat2 = lmp(log(exp) ~ region, data = beta_L)
summary(stat1)
Anova(stat1)
anova(stat1)
anova(stat1, stat2, test='Chisq')


perm_data = ddply(.data = beta_L, .variables = .( subj, model, region ), .fun = function(x){
  to_return = data.frame(exp = mean(x$exp))
  return(to_return)
})
perm = ezPerm(data = perm_data, dv = log(exp), wid = .(subj), within = NULL, between = .(model, region))
perm





