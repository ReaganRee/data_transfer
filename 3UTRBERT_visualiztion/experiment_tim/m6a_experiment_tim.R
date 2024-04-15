install.packages("hrbrthemes")
library(ggplot2)
library(dplyr)
library(patchwork) # To display 2 charts together
library(hrbrthemes)
# draw all data points > 0.003 that has not been normalized
data <- read.csv("//Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_tim_pos_threshold.csv", header=TRUE)

# p<-ggplot(data, aes(x = Relative_distance, y=Frequency))

# p + geom_point() + geom_smooth()

# draw all data points
data2 <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_tim_pos_all.csv", header=TRUE)
#p2<-ggplot(data2, aes(x = Relative_distance, y=Frequency))

#p2 + geom_point() + geom_smooth()

# draw normalized result
data3 <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_tim_pos_final_normalized.csv", header=TRUE)
p3<-ggplot()


p3 + scale_y_continuous(
  
  # Features of the first axis
  name = "Frequency",
  
  # Add a second axis and specify its features
  sec.axis = sec_axis( trans=~.*0.002, name="Frequency")
) +
  theme_ipsum()+
  geom_point(shape=18, data = data3, aes(x = Relative_distance, y=Frequency*500), color = "#fd7f6f") + geom_smooth(data = data3, aes(x = Relative_distance, y=Frequency*500), color = "#fd7f6f", fill="#fd7f6f") +
  geom_point(shape=16, data = data2, aes(x = Relative_distance, y=Frequency), color = "#7eb0d5") + geom_smooth(data = data2, aes(x = Relative_distance, y=Frequency), color = "#7eb0d5", fill = "#7eb0d5") + 
  geom_point(shape=17, data = data, aes(x = Relative_distance, y=Frequency), color = "#bd7ebe") + geom_smooth(data = data, aes(x = Relative_distance, y=Frequency), color = "#bd7ebe", fill = "#bd7ebe") +
  theme(axis.title.y.right = element_text(angle =270, colour="#fd7f6f")) + theme(axis.text.y.right = element_text(colour="#fd7f6f")) + geom_vline(alpha=0.6, xintercept = 0.3, linetype = "dashed") +
  geom_vline(alpha=0.6, xintercept = 0.7, linetype = "dashed") + annotate("rect", xmin = 0, xmax = 0.3, ymin =0, ymax = 500,alpha = .1,fill = "blue") + labs(x="Relative DIstance", y="Frequency")