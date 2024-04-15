data <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/no_window_merged_no_overlap.csv", header=TRUE)
data2 <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/window20_merged_no_overlap.csv", header=TRUE)

pos_data_all <- data[1:2129,3]
neg_data_all <- data[2130:1749908,3]


pos_data <- pos_data_all[pos_data_all <= 0.5]
neg_data <- neg_data_all[neg_data_all <= 0.5]

#install.packages('mixtools')
library(ggplot2)
library(tidyverse)
library(mixtools)
theme_set(theme_classic())
library(ggpubr)
library(cowplot)




p1<-ggplot(data[which(data$attention_score <= 0.5),], aes(x = attention_score))
p1 <- p1 + geom_density(aes(color = pos_or_neg, x=attention_score, fill=pos_or_neg), lwd=0.8, show.legend = F) +
  scale_fill_manual(values=alpha(c("white", "#95cc51"), .1))+
  geom_vline(xintercept = mean(pos_data,na.rm=TRUE), col = "#95cc51", lwd = 0.3, linetype="dashed")+
  geom_vline(xintercept = mean(neg_data,na.rm=TRUE), col = "#3c7781", lwd = 0.3, linetype="dashed")+
  
  scale_color_manual(values=c("#3c7781", "#95cc51")) + 
  labs(title="Attention Density Plot",x="Attention scores", y = "Density", color='')+
  theme(text = element_text(size = 14))+
  annotate("rect", xmin = 0, xmax = 0.019, ymin =2.87, ymax = 3,alpha = .23,fill = "white", color="#3c7781")+
  annotate("text", x = 0.07-0.01, y = 2.935, label = "Negative sites",
           color="Black",size = 5, angle=0 ) +
  
  annotate("rect", xmin = 0, xmax = 0.019, ymin =2.72, ymax = 2.85,alpha = .23,fill = "#95cc51", color="#95cc51")+
  annotate("text", x = 0.0685-0.01, y = 2.785, label = "Positive sites",
           color="Black",size = 5, angle=0 ) +
  geom_segment(aes(x=0, xend=0.019, y=2.64, yend=2.64), color="#3c7781",
               linetype="dashed")+
  geom_segment(aes(x=0, xend=0.019, y=2.69, yend=2.69), color="#95cc51",
               linetype="dashed") +
  annotate("text", x = 0.0782-0.01, y = 2.665, label = "Average attention",
           color="Black",size = 5, angle=0 )

ggsave(p1,file = "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/m6a_attention_score_no_window.png", width = 10.5, height = 6, type = "cairo", dpi = 800) 



p<-ggplot(data2[which(data2$attention_score <= 0.5),], aes(x = attention_score))
p <- p + geom_density(aes(color = pos_or_neg, x=attention_score, fill=pos_or_neg), lwd=0.8, show.legend = F) +
  scale_fill_manual(values=alpha(c("white", "#95cc51"), .1))+
  geom_vline(xintercept = 0.305, col = "#95cc51", lwd = 0.3, linetype="dashed")+
  geom_vline(xintercept = mean(neg_data,na.rm=TRUE), col = "#3c7781", lwd = 0.3, linetype="dashed")+
  
  scale_color_manual(values=c("#3c7781", "#95cc51")) + 
  labs(title="Attention Density Plot",x="Attention scores", y = "Density", color='')+
  theme(text = element_text(size = 14))  +
  annotate("rect", xmin = 0, xmax = 0.019, ymin =2.87, ymax = 3,alpha = .23,fill = "white", color="#3c7781")+
  annotate("text", x = 0.07-0.01, y = 2.935, label = "Negative sites",
           color="Black",size = 5, angle=0 ) +
  
  annotate("rect", xmin = 0, xmax = 0.019, ymin =2.72, ymax = 2.85,alpha = .23,fill = "#95cc51", color="#95cc51")+
  annotate("text", x = 0.0685-0.01, y = 2.785, label = "Positive sites",
           color="Black",size = 5, angle=0 ) +
  geom_segment(aes(x=0, xend=0.019, y=2.64, yend=2.64), color="#3c7781",
                                                             linetype="dashed")+
  geom_segment(aes(x=0, xend=0.019, y=2.69, yend=2.69), color="#95cc51",
               linetype="dashed") +
  annotate("text", x = 0.0782-0.01, y = 2.665, label = "Average attention",
           color="Black",size = 5, angle=0 )
  

ggsave(p,file = "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/m6a_attention_score_20_window.png", width = 10.5, height = 6, type = "cairo", dpi = 1200) 

# #install.packages("plotmm")
# #install.packages("patchwork")
# #install.packages("EMCluster")
# #install.packages("flexmix")
# library(plotmm)
# library(patchwork)
# library(EMCluster)
# library(flexmix)
# 
# plot_mm(out,2)
# 
# 
# 
# 
# 
# 
# 
# 
# pos_data_all <- data[1:4829,3]
# neg_data_all <- data[4830:1937394,3]
# 
# pos_data <- pos_data_all[pos_data_all <= 0.5]
# neg_data <- neg_data_all[neg_data_all <= 0.5]
# 
# 
# 
# ks.test(pos_data, neg_data, alternative="two.sided")
