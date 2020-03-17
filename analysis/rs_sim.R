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

raw_data = read.csv(paste('data/wpli/', sens, '/rs_sim.csv', sep = ''))

plot = ggplot(raw_data, aes(x = sg, y = sim, fill = band) )
plot + geom_violin(aes(fill = band), trim = FALSE, position = position_dodge(0.75)) + 
  geom_boxplot(aes(fill = band), width = 0.15,
               position = position_dodge(0.75)
  ) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)
  scale_fill_manual(values = wes_palette("Royal1",3)) + 
  #scale_fill_manual(values = wes_palette("Set2")) + 
  labs(x = 'Loading', y = 'Cosine Similarity to Resting State')  + theme_minimal()
ggsave(paste(sens, '_rs_sim.pdf', sep = ''))

# stats
fit = lmp(sim ~ sg + band, raw_data)
summary(fit)
Anova(fit)

# display medians
data = dplyr::summarise(group_by(raw_data, sg, band), med = median(sim))
data
