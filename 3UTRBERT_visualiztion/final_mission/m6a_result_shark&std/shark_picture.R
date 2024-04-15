A <- read.csv("/Users/reagan/Desktop/final_mission/m6a_result_shark&std/mcc.csv", header = T)
library(ggplot2)
library(ggsci)
#install.packages("ggsci")
#BiocManager::install("rphylopic")
library(tidyverse)
library(rphylopic)

rayurl = "/Users/reagan/Downloads/icons/rna_icon_3mer.png"
library(grid)
raylogo1 = readPNG(getURLContent(rayurl), native = T)
im1 <- rasterGrob(raylogo1, interpolate=TRUE)

rayur2 = "/Users/reagan/Downloads/icons/rna_icon_4mer.png"
raylogo2 = readPNG(getURLContent(rayur2), native = T)
im2 <- rasterGrob(raylogo2, interpolate=TRUE)


rayur3 = "/Users/reagan/Downloads/icons/rna_icon_5mer.png"
raylogo3 = readPNG(getURLContent(rayur3), native = T)
im3 <- rasterGrob(raylogo3, interpolate=TRUE)

rayur4 = "/Users/reagan/Downloads/icons/rna_icon_6mer.png"
raylogo4 = readPNG(getURLContent(rayur4), native = T)
im4 <- rasterGrob(raylogo4, interpolate=TRUE)

rayur5 = "/Users/reagan/Downloads/icons/rna_icon_deepm6aseq.png"
raylogo5 = readPNG(getURLContent(rayur5), native = T)
im5 <- rasterGrob(raylogo5, interpolate=TRUE)

rayur6 = "/Users/reagan/Downloads/icons/rna_icon_iMRM.png"
raylogo6 = readPNG(getURLContent(rayur6), native = T)
im6 <- rasterGrob(raylogo6, interpolate=TRUE)

rayur7 = "/Users/reagan/Downloads/icons/rna_icon_SCRAMP.png"
raylogo7 = readPNG(getURLContent(rayur7), native = T)
im7 <- rasterGrob(raylogo7, interpolate=TRUE)

rayur8 = "/Users/reagan/Downloads/icons/rna_icon_WHISTLE.png"
raylogo8 = readPNG(getURLContent(rayur8), native = T)
im8 <- rasterGrob(raylogo8, interpolate=TRUE)



ggplot()+
  geom_point(shape=15, size=10, data = A, aes(x = cell_line, y=X3mer), color="#AE2012", alpha=0.85)+
  geom_point(shape=15, size=10, data = A, aes(x = cell_line, y=X4mer), color="#005F73", alpha=0.85)+
  geom_point(shape=15, size=10, data = A, aes(x = cell_line, y=X5mer), color="#0A9396", alpha=0.85)+
  geom_point(shape=15, size=10, data = A, aes(x = cell_line, y=X6mer), color="#94D2BD", alpha=0.85)+
  geom_point(shape=15, size=10, data = A, aes(x = cell_line, y=DeepM6ASeq), color="#EE9B00", alpha=0.85)+
  geom_point(shape=15, size=10, data = A, aes(x = cell_line, y=iMRM), color="#CA6702", alpha=0.85)+
  geom_point(shape=15, size=10, data = A, aes(x = cell_line, y=SCRAMP), color="#91b12f", alpha=0.85)+
  geom_point(shape=15, size=10, data = A, aes(x = cell_line, y=WHISTLE), color="#9B5DE5", alpha=0.85)+
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line=element_line(colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10,angle = 90,
                                   vjust = 0.5,hjust = 0.8, 
                                   color = 'black'),
        axis.text.y = element_text(size = 10, color = 'black'),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = "NA")+
  labs(y=expression('MCC'))+
  ylim(0.5, 1.16)+
  geom_hline(yintercept = 1.02, cex=0.5)+
  geom_hline(yintercept = 1.13, cex=0.5)+
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5), color="grey", linetype=2,cex=0.5)+
  annotate(geom = 'text', label="MCC", x=5.5, y=1.157, size=10)+
  annotation_custom(im1, xmin=0.5, xmax=1.5, ymin=0.6, ymax=0.9)+
  annotation_custom(im1, xmin=1.5, xmax=2.5, ymin=0.6, ymax=0.9)+
  annotation_custom(im1, xmin=2.5, xmax=3.5, ymin=0.6, ymax=0.9)+
  annotation_custom(im2, xmin=3.5, xmax=4.5, ymin=0.65, ymax=1)+
  annotation_custom(im3, xmin=4.5, xmax=5.5, ymin=0.65, ymax=1)+
  annotation_custom(im3, xmin=5.5, xmax=6.5, ymin=0.65, ymax=1)+
  annotation_custom(im4, xmin=6.5, xmax=7.5, ymin=0.65, ymax=1)+
  annotation_custom(im4, xmin=7.5, xmax=8.5, ymin=0.65, ymax=1)


