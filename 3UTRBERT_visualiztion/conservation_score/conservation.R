data <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/conservation_score/1000_conservation_20_window_all_small_range.csv", header=TRUE)


right_data_all <- data[1:2955,2]
left_data_all <- data[2956:3321,2]


right_data <- right_data_all[right_data_all <= 2]
left_data <- left_data_all[left_data_all <= 2]

right_data <- right_data[right_data >= -1]
left_data <- left_data[left_data >= -1]

#install.package(“ggplot2”)
#install.packages('mixtools')
library(ggplot2)
library(tidyverse)
library(mixtools)
theme_set(theme_classic())

p<-ggplot(data, aes(x = conservation_score))

p + geom_density(aes(color = left_or_right)) + xlim(-1,2)+

  scale_color_manual(values=c("#999999", "#E69F00"), labels=c("left", "right")) + theme(legend.position=c(0.1,0.9))+
  labs(title="Conservation Score Density Plot",x="Conservation scores", y = "Density", color='')+
  theme(legend.text = element_text(size = 14))+ geom_vline(xintercept = mean(right_data,na.rm=TRUE), col = "#E69F00", lwd = 0.6, linetype="dashed")+
  geom_vline(xintercept = mean(left_data,na.rm=TRUE), col = "#999999", lwd = 0.6, linetype="dashed")

ks.test(right_data_all, left_data_all, alternative="two.sided")

data2 <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/conservation_score/1000_conservation_no_window.csv", header=TRUE)
pos_data_all2 <- data2[1:4195,2]
neg_data_all2 <- data2[4195:360573,2]


#pos_data <- pos_data_all[pos_data_all <= 0.5]
#neg_data <- neg_data_all[neg_data_all <= 0.5]

#install.package(“ggplot2”)
#install.packages('mixtools')
library(ggplot2)
library(tidyverse)
library(mixtools)
theme_set(theme_classic())

p<-ggplot(data2, aes(x = attention_score))

p + geom_density(aes(color = pos_or_neg)) + xlim(-0.1,0.75)+ ylim(0,4) + 
  
  scale_color_manual(values=c("#999999", "#E69F00"), labels=c("Negative", "Positive")) + theme(legend.position=c(0.1,0.9))+
  labs(title="Conservation Score Density Plot",x="Conservation scores", y = "Density", color='')+
  theme(legend.text = element_text(size = 14))+ geom_vline(xintercept = mean(pos_data_all2,na.rm=TRUE), col = "#E69F00", lwd = 0.6, linetype="dashed")+
  geom_vline(xintercept = mean(neg_data_all2,na.rm=TRUE), col = "#999999", lwd = 0.6, linetype="dashed")

ks.test(pos_data_all2, neg_data_all2, alternative="two.sided")