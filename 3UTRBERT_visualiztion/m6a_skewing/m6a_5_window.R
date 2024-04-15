data <- read.csv("/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/window15_merged.csv", header=TRUE)


#install.package(“ggplot2”)
library(ggplot2)
p<-ggplot(data, aes(x = attention_score))

p + geom_density(adjust = 0.000001)
p + geom_density(aes(color = pos_or_neg))
p + geom_density(aes(fill = pos_or_neg), alpha=0.4)



p + geom_density(aes(fill = pos_or_neg), alpha=0.4) + xlim(-0.0005,100)
