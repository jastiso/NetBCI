## Get SG data for wpli
library(R.matlab)
library(ggplot2)
library(dplyr)
library(coin)
library(lmPerm)
library(lmerTest)
library(car)
library(aplpack)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
setwd("/Users/stiso/Documents/R/NetBCI/")
bands = c('alpha', 'beta', 'low_gamma', 'gamma')
nSubj = 20
sens = 'grad'

#############################################################################################

# Skew

#############################################################################################
S = readMat(paste('data/wpli/', sens, '/Ssg.mat', sep = ''))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), S = numeric(0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.high[,b_idx]
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.other[,b_idx]
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.low[,b_idx]
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
mypallette = colorRampPalette(brewer.pal(4, 'Blues'))
plot = ggplot(plot_data, aes(x = band, y = S, col = cond, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = rev(brewer.pal(4,'Blues'))) + 
  labs(x = 'Subgraph', y = 'Skew')  + theme_minimal()
ggsave(paste(sens, '_Ssg_wpli', '.png', sep = ''))



# stats
# anova
fit1 = lm(S ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(S ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# anova beta
tmp = filter(plot_data, band == 'beta')
statb1 = lmp(S ~ cond, data = tmp)
summary(statb1)
# post hoc
stat1 = lmp(S ~ cond, data = plot_data)
summary(stat1)



#############################################################################################

# Median

#############################################################################################
# load data
St = readMat(paste('data/wpli/', sens, '/St.mat', sep = ''))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), Med = numeric(0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = St$St.high[,b_idx]
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = St$St.other[,b_idx]
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = St$St.low[,b_idx]
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
plot = ggplot(plot_data, aes(x = band, y = Med, col = cond, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = rev(brewer.pal(4,"Greens"))) + 
  labs(x = 'Subgraph', y = 'Median')  + theme_minimal()
ggsave(paste(sens, '_Med_wpli', '.png', sep = ''))




# stats
# anova
fit1 = lm(Med ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(Med ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# anova beta
tmp = filter(plot_data, band == 'beta')
Medatb1 = lmp(Med ~ cond, data = tmp)
summary(Medatb1)
# post hoc
stat1 = lmp(Med ~ cond, data = plot_data)
summary(stat1)

