# Load dataset from github
#mydata <- read.table("/Users/reagan/Desktop/3UTRBERT_visualiztion/small_graph_under_landscape/small_graph.csv", header=TRUE)
mydata <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/small_graph_under_landscape/GRWD1_HepG2/small_graph_averaged.csv", header=TRUE)
# Make the histogram
library(ggplot2)
library(ggalt)


#install.packages("ggalt")

p <- ggplot(data = mydata,
            mapping = aes(
              x = X.1,
              y = X0))
p +  theme_classic() + theme(panel.border = element_rect(color = "black",
                                                              fill = NA,
                                                              size = 2)) +
  geom_xspline( spline_shape = -0.7) +
  geom_vline(alpha=1.3, xintercept = 0, linetype = "dotted") + labs(x = "Distance from peak center (bp)", y = "Avg attention") 
