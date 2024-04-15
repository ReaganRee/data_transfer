data <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/window20_merged_no_overlap.csv", header=TRUE)
library(plotmm)
library(patchwork)
library(EMCluster)
library(flexmix)
library(ggplot2)
library(tidyverse)
library(mixtools)
library(AdaptGauss)

pos_data_all <- data[1:2129,3]
neg_data_all <- data[2130:1749908,3]


# pos_data <- pos_data_all[pos_data_all <= 0.5] * 2
# neg_data <- neg_data_all[neg_data_all <= 0.5] * 2
pos_data <- pos_data_all[pos_data_all <= 0.5] 
neg_data <- neg_data_all[neg_data_all <= 0.5] 

set.seed(100)
d = c(pos_data)
d <- d[!is.na(d)]
out<-normalmixEM(d, epsilon=0.16, arbvar=TRUE, fast=TRUE)
print(out)

theme_set(theme_classic())

p<-ggplot(data[which(data$attention_score <= 0.5),], aes(x = attention_score))
p <- p + geom_density(aes(color = pos_or_neg, x=attention_score, fill=pos_or_neg), lwd=0.8, show.legend = F) +
  scale_fill_manual(values=alpha(c("#B8b5b2", "#d2ffcf"), .1))+
  geom_vline(xintercept = mean(pos_data,na.rm=TRUE), col = "#99c391", lwd = 0.3, linetype="dashed")+
  geom_vline(xintercept = mean(neg_data,na.rm=TRUE), col = "#747474", lwd = 0.3, linetype="dashed")+

  scale_color_manual(values=c("#747474", "#99c391")) + 
  labs(title="Attention Density Plot",x="Attention scores", y = "Density", color='')+
  theme(legend.text = element_text(size = 14))+
  stat_function(geom = "line", fun = plot_mix_comps_normal, # here is the function
                args = list(out$mu[1], out$sigma[1], lam = out$lambda[1]),
                colour = "#eda3a0", lwd = 0.8)+
  stat_function(geom = "line", fun = plot_mix_comps_normal, # here again as k = 2
                args = list(out$mu[2], out$sigma[2], lam = out$lambda[2]),
                colour = "#366aa6", lwd = 0.8) + scale_x_continuous(expand = c(0, 0.0005)) + 
  theme(plot.margin=margin(t = 40,  # Top margin
                           r = 40,  # Right margin
                           b = 40,  # Bottom margin
                           l = 40), panel.spacing.y = unit(5, "lines"))+
  annotate("rect", xmin = 0.02, xmax = 0.035, ymin =2.9, ymax = 3,alpha = .23,fill = "#B8b5b2", color="#747474")+
  annotate("text", x = 0.058, y = 2.95, label = "Negative sites",
           color="Black",size = 2.2, angle=0 ) +
  
  annotate("rect", xmin = 0.02, xmax = 0.035, ymin =2.77, ymax = 2.87,alpha = .23,fill = "#d2ffcf", color="#99c391")+
  annotate("text", x = 0.0565, y = 2.82, label = "Positive sites",
           color="Black",size = 2.2, angle=0 ) +
  
  annotate("rect", xmin = 0.02, xmax = 0.035, ymin =2.64, ymax = 2.74, color="#eda3a0", fill=NA)+
  annotate("text", x = 0.0565, y = 2.69, label = "Fake positive",
           color="Black",size = 2.2, angle=0 ) +
  
  annotate("rect", xmin = 0.02, xmax = 0.035, ymin =2.51, ymax = 2.61, color="#366aa6", fill=NA)+
  annotate("text", x = 0.056, y = 2.56, label = "True positive",
           color="Black",size = 2.2, angle=0 )+
  geom_segment(aes(x=0.02, xend=0.035, y=2.4, yend=2.4), color="#747474",
               linetype="dashed")+
  geom_segment(aes(x=0.02, xend=0.035, y=2.45, yend=2.45), color="#99c391",
              linetype="dashed") +
  annotate("text", x = 0.06, y = 2.425, label = "Average attention",
           color="Black",size = 2.2, angle=0 )

ggsave(p,file = "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/m6a_attention_score.png", width = 10.5, height = 6, type = "cairo", dpi = 800) 







#print(Intersect2Mixtures(0.2185390, 0.07763533, 0.5525921, 0.4034567, 0.06057512, 0.4474079))
# pos_data_all <- data[1:4829,3]
# neg_data_all <- data[4830:1937394,3]
# 
# pos_data <- pos_data_all[pos_data_all <= 0.5]
# neg_data <- neg_data_all[neg_data_all <= 0.5]
# 
# 
# 
# ks.test(pos_data, neg_data, alternative="two.sided")
