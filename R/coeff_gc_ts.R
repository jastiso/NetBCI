## Get SG data for gC, time shifted null model
library(R.matlab)
library(ggplot2)
library(dplyr)
library(coin)
library(lmPerm)
library(lme4)
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
S = readMat(paste('data/gc/', sens, '/S_ts.mat', sep = ''))
S_corr = readMat(paste('data/gc/', sens, '/S_corr_ts.mat', sep = ''))
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
  scale_fill_manual(values = rev(brewer.pal(4,'Greys'))) + 
  labs(x = 'Subgraph', y = 'Skew')  + theme_minimal()
ggsave(paste(sens, '_S_gc_ts', '.png', sep = ''))


# scatter plot
S_plot_data = data.frame(BE = S_corr$S.exp.corr[,2], S = S_corr$S.exp.corr[,1], band = S_corr$sbo, subj = S_corr$sso)
S_plot_data$band = as.factor(S_plot_data$band)
S_plot_data$subj = as.factor(S_plot_data$subj)
scatterplot = ggplot(S_plot_data, aes(x = log(BE), y = S, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'Skew') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = rev(brewer.pal(4,'Greys')))
ggsave(paste(sens, '_gc_dot_s_ts', '.png', sep = ''))

# beta
S_beta = filter(S_plot_data, band == '2') # beta is level 2
scatterplot = ggplot(S_beta, aes(x = log(BE), y = S, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'Skew') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = rev(brewer.pal(4,'Greys')))
ggsave(paste(sens, '_gc_beta_s_ts.png', sep = ''))


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
# corr
stat2 = lmp(S ~ log(BE) + band, data = S_plot_data)
summary(stat2)
# corr beta
statb = lme4::lmer(S ~ log(BE) + (1|subj), data = S_beta)
Anova(statb)

#############################################################################################

# Mean

#############################################################################################

# load data
M = readMat(paste('data/gc/', sens, '/M_ts.mat', sep = ''))
M_corr = readMat(paste('data/gc/', sens, '/M_corr_ts.mat', sep = ''))
plot_data = data.frame(cond = character(length = 0), band = character(length = 0), M = numeric(0), stringsAsFactors = FALSE)
cnt = 1
b_idx = 1;
for (b in bands){
  plot_data[cnt:(cnt+nSubj-1),1] = rep('high BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M$M.high[,b_idx]
  cnt = cnt + nSubj
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('other BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M$M.other[,b_idx]
  cnt = cnt + nSubj
  
  
  plot_data[cnt:(cnt+nSubj-1),1] = rep('low BE', times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),2] = rep(b, times = nSubj)
  plot_data[cnt:(cnt+nSubj-1),3] = M$M.low[,b_idx]
  cnt = cnt + nSubj
  
  b_idx = b_idx + 1
}
plot_data$cond = factor(plot_data$cond, levels(factor(plot_data$cond))[1:3])
plot_data$band = as.factor(plot_data$band)


# box plot
plot = ggplot(plot_data, aes(x = band, y = M, col = cond, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_color_manual(values = c(rgb(0,0,0, alpha = 1), rgb(.3,.3,.3, alpha = 1), rgb(.6,.6,.6, alpha = 1))) + 
  scale_fill_manual(values =  rev(brewer.pal(4,'Greys'))) + 
  labs(x = 'Subgraph', y = 'Energy')  + theme_minimal()
ggsave(paste(sens, '_M_gc_ts.png', sep = ''))


# scatter plot
M_plot_data = data.frame(BE = M_corr$M.exp.corr[,2], M = M_corr$M.exp.corr[,1], band = M_corr$mbo, subj = M_corr$mso)
M_plot_data$band = as.factor(M_plot_data$band)
M_plot_data$subj = as.factor(M_plot_data$subj)
scatterplot = ggplot(M_plot_data, aes(x = log(BE), y = M, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'Mean') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  rev(brewer.pal(4,'Greys')))
ggsave(paste(sens, '_dot_m_gc_ts.png', sep = ''))

# beta
M_beta = filter(M_plot_data, band == '2') # beta is level 2
M_plot_data = data.frame(BE = M_corr$M.exp.corr[,2], M = M_corr$M.exp.corr[,1], band = M_corr$mbo, subj = M_corr$mso)
M_plot_data$band = as.factor(M_plot_data$band)
M_plot_data$subj = as.factor(M_plot_data$subj)
scatterplot = ggplot(M_beta, aes(x = log(BE), y = M, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'Mean') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  rev(brewer.pal(4,'Greys')))
ggsave(paste(sens, '_beta_m_gc_ts.png', sep = ''))


# stats
# anova
fit1 = lm(M ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(M ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# anova beta
tmp = filter(plot_data, band == 'beta')
statb1 = lmp(M ~ cond, data = tmp)
summary(statb1)
# post hoc
stat1 = lmp(M ~ cond, data = plot_data)
summary(stat1)
# corr
stat2 = lmp(M ~ log(BE), data = M_plot_data)
summary(stat2)
# beta
statb = lme4::lmer(M ~ log(BE) + (1|subj), data = M_beta)
Anova(statb)


#############################################################################################

# Energy

#############################################################################################

# load data
E = readMat(paste('data/gc/', sens, '/E_ts.mat', sep = ''))
E_corr = readMat(paste('data/gc/', sens, '/E_corr_ts.mat', sep = ''))
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
  scale_fill_manual(values =  rev(brewer.pal(4,'Greys'))) + 
  labs(x = 'Subgraph', y = 'Energy')  + theme_minimal()
ggsave(paste(sens, '_E_gc_ts', '.png', sep = ''))


# scatter plot
E_plot_data = data.frame(BE = E_corr$E.exp.corr[,2], E = E_corr$E.exp.corr[,1], band = E_corr$ebo, subj = E_corr$eso)
E_plot_data$band = as.factor(E_plot_data$band)
E_plot_data$subj = as.factor(E_plot_data$subj)
scatterplot = ggplot(E_plot_data, aes(x = log(BE), y = E, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'ENergy') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  rev(brewer.pal(4,'Greys')))
ggsave(paste(sens, '_dot_e_gc_ts', '.png', sep = ''))

# beta
E_beta = filter(E_plot_data, band == '2') # beta is level 2
scatterplot = ggplot(E_beta, aes(x = log(BE), y = E, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'ENergy') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values =  rev(brewer.pal(4,'Greys')))
ggsave(paste(sens, '_beta_e_gc_ts', '.png', sep = ''))


# stats
# anova
fit1 = lm(E ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(E ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# anova beta
tmp = filter(plot_data, band == 'beta')
statb1 = lmp(E ~ cond, data = tmp)
summary(statb1)
# post hoc
stat1 = lmp(E ~ cond, data = plot_data)
summary(stat1)
# corr
stat2 = lmp(E ~ log(BE), data = E_plot_data)
summary(stat2)
# corr
statb = lme4::lmer(E ~ log(BE) + (1|subj), data = E_beta)
Anova(statb)

#############################################################################################

# Entropy

#############################################################################################
# load data
H = readMat(paste('data/gc/', sens, '/H_ts.mat', sep = ''))
H_corr = readMat(paste('data/gc/', sens, '/H_corr_ts.mat', sep = ''))
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
  scale_fill_manual(values = rev(brewer.pal(4,'Greys'))) + 
  labs(x = 'Subgraph', y = 'Shannon Entropy')  + theme_minimal()
ggsave(paste(sens, '_H_gc_ts.png', sep = ''))


# scatter plot
H_plot_data = data.frame(BE = H_corr$H.exp.corr[,2], H = H_corr$H.exp.corr[,1], band = H_corr$hbo, subj = H_corr$hso)
H_plot_data$band = as.factor(H_plot_data$band)
H_plot_data$subj = as.factor(H_plot_data$subj)
scatterplot = ggplot(H_plot_data, aes(x = log(BE), y = H, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'Shannon Entropy') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = rev(brewer.pal(4,'Greys')))
ggsave(paste(sens, '_dot_h_gc_ts', '.png', sep = ''))

# beta
H_beta = filter(H_plot_data, band == '2') # beta is level 2
scatterplot = ggplot(H_beta, aes(x = log(BE), y = H, color = band))
scatterplot + geom_point(size = 6) + labs(x = 'log(Behavioral Expression)', y = 'Shannon Entropy') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() + scale_color_manual(values = rev(brewer.pal(4,'Greys')))
ggsave(paste(sens, '_beta_h_gc_ts', '.png', sep = ''))

# stats
# anova
fit1 = lm(H ~ cond + band, data=  plot_data)
summary(fit1)
Anova(fit1)
fit2 = lm(H ~ band*cond, data = plot_data)
summary(fit2)
Anova(fit2)
# anova beta
tmp = filter(plot_data, band == 'beta')
statb1 = lmp(H ~ cond, data = tmp)
summary(statb1)
# post hoc
stat1 = lmp(H ~ cond, data = plot_data)
summary(stat1)
# corr
stat2 = lmp(H ~ log(BE), data = H_plot_data)
summary(stat2)
# corr beta
statb = lme4::lmer(H ~ log(BE) + (1|subj), data = H_beta)
Anova(statb)

# scatter plot
H_plot_data$E = E_plot_data$E
s = ggplot(H_plot_data, aes(x = E, y = H, color = log(BE)))
s + geom_point(size = 6) + labs(x = 'Energy', y = 'Shannon Entropy') + 
  geom_smooth(method="lm", color = 'black') + theme_minimal() #+ scale_color_continuous(colors = mypallette)




