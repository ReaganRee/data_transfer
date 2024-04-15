library(BiocManager)
#install("motifStack")
library("motifStack")

motifs<-importMatrix(dir(path ="/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/motifs_PWM", full.names = TRUE))

set.seed(10)
pfms <- motifs
domain <- read.csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/domain_withGC.csv", header = TRUE, sep=",")



hc <- clusterMotifs(pfms)
library(ade4)
phylog <- hclust2phylog(hc)
## reorder the pfms by the order of hclust
leaves <- names(phylog$leaves)
print(leaves)
pfms <- pfms[leaves]


## extract the motif signatures
motifSig <- motifSignature(pfms, phylog, cutoffPval = 0.0001, min.freq=1)
sig <- signatures(motifSig)
gpCol <- sigColor(motifSig)


library(RColorBrewer)
color <- brewer.pal(10, "Set3")
## plot logo stack with radial style

pfms <- DNAmotifToRNAmotif(pfms)



pfms <- lapply(pfms, trimMotif, t=0.4)



#pfms <- trimMotif(pfms, t=0.4)
#motifStack(pfms, layout="radialPhylog",  #circle is the most inner circle radius,
#           circle=4.8, cleaves = 0.3, 
#           clabel.leaves = 0.5,  #clabel.labels is the evolution line and RBP label size
#           col.bg=rep(color, each=26), col.bg.alpha=0.3, 
#           col.leaves=rep(color, each=26),
#           col.inner.label.circle=rep(color, each=26), 
#           inner.label.circle.width=0.15, #里面环的宽度
#           col.outer.label.circle=rep(color, each=26), 
#           outer.label.circle.width=2, #外面环的宽度
#           circle.motif=5.2, motifScale="logarithmic",
#           angle=360)



RRM_color_list <- list()
KH_color_list <- list()
PUF_color_list <- list()
CSD_color_list <- list()
counter = 1

family_dict <- c("FMR"="#E5FFCC", "SRP"="#CCFFFF", "SSB"="#E5CCFF", "HNR"="#FFCCFF", "NXF"="#FFE5CC", "PUM"="#00CCCC", "SRS"="#FD7272",
                 "CEL"="#E28EFF", "PTB"="#8EC7FF", "RBM"="#8EFFF1", "YBX"="#E5FF8E", "TIA"="#FFDD6C","MBN"="#FFBB6C", "DAZ" = "#E4CFCC", 
                 "ELA" = "#8EC7AB", "NOV" = "#9EE9FF", "PCB" = "#A5BF3E", "FUS" = "blue", "SFP" = "red", "KHS" = 'brown')
tree_color_list <- list()


none_color = "#EFEFEF"
for (item in pfms) {
  tree_color_list[counter] <- family_dict[substr(strsplit(item$name, split="_")[[1]][[1]],1,3)]
  i = match(item$name, domain$name)
  if (domain$domain[i] == "KH"){
    KH_color_list[counter] <- "#003333"
    CSD_color_list[counter] <- none_color
    RRM_color_list[counter] <- none_color
    PUF_color_list[counter] <- none_color
  }
  if (domain$domain[i] == "RRM"){
    
    KH_color_list[counter] <- none_color
    CSD_color_list[counter] <- none_color
    RRM_color_list[counter] <- "#009999"
    PUF_color_list[counter] <- none_color
  }
  if (domain$domain[i] == "Other"){
    
    KH_color_list[counter] <- none_color
    CSD_color_list[counter] <- none_color
    RRM_color_list[counter] <- none_color
    PUF_color_list[counter] <- none_color
  }
  if (domain$domain[i] == "PUF"){

    KH_color_list[counter] <- none_color
    CSD_color_list[counter] <- none_color
    RRM_color_list[counter] <- none_color
    PUF_color_list[counter] <- "brown"
  }
  if (domain$domain[i] == "CSD"){

    KH_color_list[counter] <- none_color
    CSD_color_list[counter] <- "#F69156"
    RRM_color_list[counter] <- none_color
    PUF_color_list[counter] <- none_color
  }
  
  #domain_list[counter] <- strsplit(item$name, split="_")[[1]][[2]]
  
  #name_list[counter] <- strsplit(item$name, split="_")[[1]][[1]]
  counter <- counter + 1
}
RRM_color_list <- trimws(RRM_color_list)
CSD_color_list <- trimws(CSD_color_list)
KH_color_list <- trimws(KH_color_list)
PUF_color_list <- trimws(PUF_color_list)
tree_color_list <- trimws(tree_color_list)
#

#for (i in 1:75){
#  pfms[[i]]$color["A"] <- "#81B046"
#  pfms[[i]]$color["G"] <- "#EEAD55"
#  pfms[[i]]$color["C"] <- "#354973"
#  pfms[[i]]$color["U"] <- "#CC2829"
#}
new_tree_color_list <- c(rep("#009999", each=56), rep("red", each=4), rep("#003333", each=15),rep("#FFCCFF", each=8), 
                         rep("brown", each=7),rep("blue", each=12), rep("#FC9C60", each=13), rep("dark green", each=2),
                         rep("#9EE9FF", each=10))


motifCircos(phylog=phylog, pfms=DNAmotifAlignment(pfms),motifScale="linear", angle=360, r.rings=c(2, 0.11,0.13,0.13,0.16,0.15,0.22,0.17), 
            col.rings=list(rep("white", 128), 
                           KH_color_list, 
                           rep("white", 128),
                           PUF_color_list,
                           rep("white", 128),
                           CSD_color_list,
                           rep("white", 128),
                           RRM_color_list),
            col.tree.bg=new_tree_color_list, col.tree.bg.alpha=.3,
            #col.leaves.bg=tree_color_list,col.leaves.bg.alpha=.3,
            R=3.5,
            r.pfms = 1.5,
            cleaves=.1, r.tree=1.5, r.leaves=.1, clabel.leaves=1.2,
            #col.inner.label.circle = color_list,
            plotIndex=FALSE, IndexCex=.7)


legend(-3.4, 2.4, legend=c("RRM","CSD", "KH", "PUF", "Other"),
       fill=c("#009999", "#F69156", "#003333", "brown","#EFEFEF"),
       border="white", lty=NULL, bty = "n")



                                                                                                                                 
                                                                                                                                       
                                                                                                                                         
motifCircos(phylog=phylog, pfms=DNAmotifAlignment(pfms),motifScale="logarithmic", angle=360,# r.rings=c(0.35), 
            #col.rings=list(rep("white", 75)), 
            col.tree.bg=new_tree_color_list, col.tree.bg.alpha=0.3,
            #col.leaves.bg=tree_color_list,col.leaves.bg.alpha=.3,
            R=3.5,
            r.pfms = 1.5,
            cleaves=0.1, r.tree=0.9, r.leaves=0.5, clabel.leaves=0.5,
            #col.inner.label.circle = color_list,
            plotIndex=FALSE, IndexCex=.7)