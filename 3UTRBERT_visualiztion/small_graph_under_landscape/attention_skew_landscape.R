# Libraries
library(ggplot2)
library(hrbrthemes)

# Dummy data
data <- read.csv("/Users/reagan/Downloads/RBM15_sameRBP.csv")

# Chart
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)



# With transparency (right)
ggplot(data=data, aes(x=attention, group=label, fill=label)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()