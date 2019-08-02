# Consistency QQ plot

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
bands = c('alpha', 'beta', 'low_gamma')
nSubj = 20
sens = 'grad'

raw_data = readMat('data/wpli/consistency.mat')

data = data.frame(band = unlist(raw_data$band), model = unlist(raw_data$model), rank = unlist(raw_data$cond), consistency = raw_data$data[1,])

###################
# Normality Tests
####################

par(mfrow=c(1,1))
qqPlot(data$consistency)
