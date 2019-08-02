# Behavior - replicating MCs paper
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
library(viridis)
setwd("/Users/stiso/Documents/R/NetBCI/")

raw_data = readMat("/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/0_Behavior/behavior_updated_trials.mat")
raw_data = raw_data$behavior.updated[,,1]$BCI[,,1]$Perf[,,1]
data = data_frame(subj = numeric(0), session = numeric(0), test = numeric(0), score = numeric(0), avg_score = numeric(0))
for (sess in names(raw_data)){
  curr = raw_data[[sess]][,,1]
  for (s in seq(1,length(curr$Runs))){
    run_data = unlist(curr$Runs[s])
    tmp_data = data_frame(subj = numeric(length(run_data)), session = numeric(length(run_data)), test = numeric(length(run_data)), score = numeric(length(run_data)),
                          avg_score = numeric(length(run_data)))
    tmp_data$subj = rep(s, length(run_data))
    tmp_data$avg_score = rep(unlist(curr$Avg[s]), length(run_data))
    tmp_data$session = rep(sess, times = length(run_data))
    tmp_data$score = run_data
    tmp_data$test = seq(1,length(run_data))
    data = rbind(data, tmp_data)
  }
}
data$test = factor(data$test)
data$subj = facotr(data$subj)
data$time = paste(data$session, data$test)
data$sess_cont = data$session
data$sess_cont[data$session == "Sess1"] = 1
data$sess_cont[data$session == "Sess2"] = 2
data$sess_cont[data$session == "Sess3"] = 3
data$sess_cont[data$session == "Sess4"] = 4
data$sess_cont = as.numeric(data$sess_cont)

data_avg = dplyr::summarise(group_by(data,subj,session), mean = mean(score), n = n())

plot = ggplot(data, aes(x = session, y = score, fill = session, color = test) )
plot + 
  geom_violin(, trim = TRUE, position = position_dodge(0.75)) + 
  geom_boxplot(width = 0.15,
               position = position_dodge(0.75)
  ) +
  #geom_line(aes(group = subj, color = subj), alpha = 0.2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.2, width = 0.15,
               position = position_dodge(0.75)) +
  #scale_fill_manual(values = wes_palette("Royal2",4)) + 
  scale_fill_manual(values = brewer.pal(4, "Greys")) + 
  #scale_fill_manual(values = c(rgb(137/255,157/255,164/255), rgb(247/255,251/255,254/255), rgb(149/255,192/255,76/255), rgb(63/255,67/255,71/255) )) +
  #scale_fill_viridis(alpha = .5) +
  scale_color_manual(values = rep("black", times = 6)) + 
  labs(x = 'Time', y = 'BCI Score')  + theme_minimal()
ggsave('behavior.pdf')

avg_data = filter(data, test == 1)
plot = ggplot(avg_data, aes(x = session, y = avg_score, fill = as.factor(sess_cont)) )
plot + 
  geom_violin( trim = TRUE, position = position_dodge(0.75)) + 
  geom_boxplot(width = 0.15,
               position = position_dodge(0.75)
  ) +
  geom_line(aes(group = subj), alpha = 0.2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, width = 0.15,
               position = position_dodge(0.75)) +
  #scale_fill_manual(values = wes_palette("Royal2",4)) + 
  scale_fill_manual(values = brewer.pal(4, "Greys")) + 
  #scale_fill_manual(values = c(rgb(149/255,192/255,76/255), rgb(247/255,251/255,254/255), rgb(137/255,157/255,164/255), rgb(53/255,57/255,61/255) )) +
  #scale_fill_viridis(alpha=1/2) +
  scale_color_manual(values = rep("grey", times = 6)) + 
  labs(x = 'Time', y = 'BCI Score')  + theme_minimal()
ggsave('behavior_avg.pdf')

# Stats
fit = with(data_avg, aov(mean ~ session + Error(as.factor(subj))))
summary(fit)



##################
# Normality test
##################


qqPlot(unique(corr_data$slope))

par(mfrow=c(2,2))
qqPlot(filter(avg_data, session == "Sess1")$score)
qqPlot(filter(avg_data, session == "Sess2")$score)
qqPlot(filter(avg_data, session == "Sess3")$score)
qqPlot(filter(avg_data, session == "Sess4")$score)

