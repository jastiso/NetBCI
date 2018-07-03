## Get SG data for AC
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
bands = c('theta', 'alpha', 'low_gamma')
nSubj = 16

#############################################################################################

# Skew

#############################################################################################
S = readMat(paste('data/ac/S_bl.mat'))
S_corr = readMat(paste('data/ac/S_corr_bl.mat'))
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
  scale_fill_manual(values = brewer.pal(3,'Greys')) + 
  labs(x = 'Subgraph', y = 'Skew')  + theme_minimal()
ggsave(paste('S_ac_bl', '.png', sep = ''))


# scatter plot
#S_plot_data = data.frame(BE = S_corr$S.exp.corr[,2], S = S_corr$S.exp.corr[,1], band = S_corr$band.order, subj = S_corr$subj.order)
#S_plot_data$band = as.factor(S_plot_data$band)
#S_plot_data$subj = as.factor(S_plot_data$subj)
#scatterplot = ggplot(S_plot_data, aes(x = log(BE), y = S, color = band))
#scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'ENergy') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = wes_palette("Moonrise1"))
#ggsave(paste('ac_dot', '_s', '.png', sep = ''))

# stats
# anova
fit1 = lm(S ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(S ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(S ~ cond, data = plot_data)
summary(stat1)
# corr
#stat2 = lmp(S ~ log(BE), data = S_plot_data)
#summary(stat2)


#############################################################################################

# Energy

#############################################################################################

# load data
E = readMat(paste('data/ac/E_bl.mat'))
E_corr = readMat(paste('data/ac/E_corr_bl.mat'))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), E = numeric(0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.high[,b_idx]
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.other[,b_idx]
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = E$E.low[,b_idx]
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
plot = ggplot(plot_data, aes(x = band, y = E, col = cond, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values =  brewer.pal(3,'Greys')) + 
  labs(x = 'Subgraph', y = 'Energy')  + theme_minimal()
ggsave(paste('E_ac_bl', '.png', sep = ''))


# scatter plot
#E_plot_data = data.frame(BE = E_corr$E.exp.corr[,2], E = E_corr$E.exp.corr[,1], band = E_corr$band.order, subj = E_corr$subj.order)
#E_plot_data$band = as.factor(E_plot_data$band)
#E_plot_data$subj = as.factor(E_plot_data$subj)
#scatterplot = ggplot(E_plot_data, aes(x = log(BE), y = E, color = band))
#scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'ENergy') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = wes_palette("Royal1"))
#ggsave(paste('dot', '_e', '.png', sep = ''))

# stats
# anova
fit1 = lm(E ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(E ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(E ~ cond, data = plot_data)
summary(stat1)
# corr
#stat2 = lmp(E ~ log(BE), data = E_plot_data)
#summary(stat2)


#############################################################################################

# Entropy

#############################################################################################
# load data
H = readMat(paste('data/ac/H_bl.mat'))
H_corr = readMat(paste('data/ac/H_corr_bl.mat'))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), H = numeric(0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H$H.high[,b_idx]
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H$H.other[,b_idx]
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = H$H.low[,b_idx]
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
plot = ggplot(plot_data, aes(x = band, y = H, col = cond, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values = brewer.pal(3,"Greys")) + 
  labs(x = 'Subgraph', y = 'Shannon Entropy')  + theme_minimal()
ggsave(paste('H_ac_bl', '.png', sep = ''))


# scatter plot
#H_plot_data = data.frame(BE = H_corr$H.exp.corr[,2], H = H_corr$H.exp.corr[,1], band = H_corr$band.order, subj = H_corr$subj.order)
#H_plot_data$band = as.factor(H_plot_data$band)
#H_plot_data$subj = as.factor(H_plot_data$subj)
#scatterplot = ggplot(H_plot_data, aes(x = log(BE), y = H, color = band))
#scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'Shannon Entropy') + 
#  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = wes_palette("GrandBudapest2"))
#ggsave(paste('dot', '_h', '.png', sep = ''))

# stats
# anova
fit1 = lm(H ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(H ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# post hoc
stat1 = lmp(H ~ cond, data = plot_data)
summary(stat1)
# corr
#stat2 = lmp(H ~ log(BE), data = H_plot_data)
#summary(stat2)



# scatter plot
H_plot_data$E = E_plot_data$E
s = ggplot(H_plot_data, aes(x = E, y = H, color = log(BE)))
s + geom_point(size = 6) + labs(x = 'Energy', y = 'Shannon Entropy') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_continuous(colors = mypallette)




