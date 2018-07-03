# Compare Statistics Across Null Models
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
setwd("/Users/stiso/Documents/R/NetBCI/")
bands = c('alpha', 'beta', 'low_gamma')
nSubj = 19
sens = 'grad'

#############################################################################################

# Skew

#############################################################################################
S_ts = readMat(paste('data/gc/', sens, '/S_ts.mat', sep = ''))
S = readMat(paste('data/gc/', sens, '/S.mat', sep = ''))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), S = numeric(0), model = character(length = 0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.high[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.other[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S$S.low[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  # time shifted
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S_ts$S.high[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S_ts$S.other[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = S_ts$S.low[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
mypallette = colorRampPalette(brewer.pal(4, 'Blues'))
plot = ggplot(plot_data, aes(x = band, y = S, col = cond, fill = model) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = rev(brewer.pal(4,'Pastel1'))) + 
  labs(x = 'Subgraph', y = 'Skew')  + theme_minimal()
ggsave(paste(sens, '_S_ts_al', '.png', sep = ''))


# stats
# anova
fit1 = lm(S ~ cond + band + model, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(S ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(S ~ cond + band + model, data = plot_data)
summary(stat1)
# anova beta
tmp = filter(plot_data, band == 'beta')
statb1 = lmp(S ~ model, data = tmp)
summary(statb1)



#############################################################################################

# Mean

#############################################################################################
M_ts = readMat(paste('data/gc/', sens, '/M_ts.mat', sep = ''))
M = readMat(paste('data/gc/', sens, '/M.mat', sep = ''))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), M = numeric(0), model = character(length = 0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M$M.high[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M$M.other[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M$M.low[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  # time shifted
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M_ts$M.high[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M_ts$M.other[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M_ts$M.low[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
mypallette = colorRampPalette(brewer.pal(4, 'Blues'))
plot = ggplot(plot_data, aes(x = band, y = M, col = cond, fill = model) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = rev(brewer.pal(4,'Pastel2'))) + 
  labs(x = 'Subgraph', y = 'Mean')  + theme_minimal()
ggsave(paste(sens, '_M_ts_al', '.png', sep = ''))


# stats
# anova
fit1 = lm(M ~ cond + band + model, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(M ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(M ~ cond + band + model, data = plot_data)
summary(stat1)
# anova beta
tmp = filter(plot_data, band == 'beta')
statb1 = lmp(M ~ model, data = tmp)
summary(statb1)


#############################################################################################

# Energy

#############################################################################################
E_ts = readMat(paste('data/gc/', sens, '/E_ts.mat', sep = ''))
E = readMat(paste('data/gc/', sens, '/E.mat', sep = ''))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), E = numeric(0), model = character(length = 0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.high[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.other[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.low[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  # time shifted
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E_ts$E.high[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E_ts$E.other[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E_ts$E.low[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
mypallette = colorRampPalette(brewer.pal(4, 'Blues'))
plot = ggplot(plot_data, aes(x = band, y = E, col = cond, fill = model) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = rev(brewer.pal(4,'Dark2'))) + 
  labs(x = 'Subgraph', y = 'Energy')  + theme_minimal()
ggsave(paste(sens, '_E_ts_al', '.png', sep = ''))


# stats
# anova
fit1 = lm(E ~ cond + band + model, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(E ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(E ~ cond + band + model, data = plot_data)
summary(stat1)
# anova beta
tmp = filter(plot_data, band == 'beta')
statb1 = lmp(E ~ model, data = tmp)
summary(statb1)


#############################################################################################

# Entropy

#############################################################################################
H_ts = readMat(paste('data/gc/', sens, '/H_ts.mat', sep = ''))
H = readMat(paste('data/gc/', sens, '/H.mat', sep = ''))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), H = numeric(0), model = character(length = 0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H$H.high[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H$H.other[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H$H.low[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('aligned', times = nSubj)
  cnt = cnt + nSubj
  
  # time shifted
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H_ts$H.high[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H_ts$H.other[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H_ts$H.low[,b_idx]
  plot_data[cnt:(cnt+nSubj-1),4] = rep('shifted', times = nSubj)
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
mypallette = colorRampPalette(brewer.pal(4, 'Blues'))
plot = ggplot(plot_data, aes(x = band, y = H, col = cond, fill = model) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = rev(brewer.pal(4,'Accent'))) + 
  labs(x = 'Subgraph', y = 'Entropy')  + theme_minimal()
ggsave(paste(sens, '_H_ts_al', '.png', sep = ''))


# stats
# anova
fit1 = lm(H ~ cond + band + model, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(H ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(H ~ cond + band + model, data = plot_data)
summary(stat1)
# anova beta
tmp = filter(plot_data, band == 'beta')
statb1 = lmp(H ~ model, data = tmp)
summary(statb1)



