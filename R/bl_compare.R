# Compare BE for BL and task

# Energy
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
bands = c('theta', 'alpha', 'beta', 'gamma')
nSubj = 6

E_bl = readMat(paste('data/E_corr_bl.mat'))
E = readMat(paste('data/E_corr.mat'))
E_plot_data = data.frame(BE = E$E.exp.corr[,2], E = E$E.exp.corr[,1], band = E$band.order, subj = E$subj.order)
E_plot_data_bl = data.frame(BE_bl = E_bl$E.exp.corr[,2],  E_bl = E_bl$E.exp.corr[,1], band = E_bl$band.order, subj = E_bl$subj.order)
E_plot_data$band = as.factor(E_plot_data$band)
E_plot_data$subj = as.factor(E_plot_data$subj)
E_plot_data_bl$band = as.factor(E_plot_data_bl$band)
E_plot_data_bl$subj = as.factor(E_plot_data_bl$subj)

plot = ggplot(E_plot_data, aes(x = band, y = BE, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_fill_manual(values = wes_palette("Zissou1")) + 
  labs(x = 'Graph', y = 'BE')  + theme_minimal()
ggsave(paste('BE', '.png', sep = ''))

plot = ggplot(E_plot_data_bl, aes(x = band, y = BE_bl, fill = band) )
plot + geom_boxplot(notch = FALSE, lwd = 1) + 
  scale_fill_manual(values = wes_palette("Zissou1")) + 
  labs(x = 'Graph', y = 'BE')  + theme_minimal() + ylim(1, 15)
ggsave(paste('BE_bl', '.png', sep = ''))

stat1 = t.test(E_plot_data$BE, E_plot_data_bl$BE_bl)
stat1
data_np = data.frame(E = c(E_plot_data$BE, E_plot_data_bl$BE_bl), Label = c(rep("Task", times = length(E_plot_data$BE)), rep("BL", times = length(E_plot_data_bl$BE_bl))))
stat_np = lmp(E ~ Label, data_np)
summary(stat_np)
