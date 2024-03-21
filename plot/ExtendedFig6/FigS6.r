#
# Author Yuxiang Li
#

### plot the genome-wide D1 vs D2 DML distribution
library("tidyverse")

dat_CG = read.table("dml_M_D1vsD2_delta10percent1e-5.mergeSS2DS.test.bedgraph.new.txt",header=F)
dat_CG$DML = ifelse(dat_CG$V4>0,"hyper","hypo")
dat_CA = read.table("dml.D1vsD2.mCA_delta10percent1e-4.mergeSS2DS.txt.bedgraph.new.txt",header=F)
dat_CA$DML = ifelse(dat_CA$V4>0,"hyper","hypo")
chr_lab <- read.table("chr_label.txt")

p1.1 <- ggplot(dat_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.1, alpha=0.3)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="red","hypo"="blue4"))

p1.2 <- p1.1 +theme_classic(base_size = 5)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=1,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,expand = c(0, 0), limits = c(0,2725537669))+ geom_hline(yintercept = 0, size=0.3)

library(gridExtra)

p2.1 <- ggplot(dat_CA,aes(x=V2,y=V4,color=DML))+geom_point(size=0.1, alpha=0.3)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="red","hypo"="blue4"))

p2.2 <- p2.1 +theme_classic(base_size = 5)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=1,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,labels = chr_lab$V1,expand = c(0, 0), limits = c(0,2725537669))+ geom_hline(yintercept = 0, size=0.3)

p <- grid.arrange(p1.2, p2.2, nrow = 2)
ggsave("figre1.13.png",plot=p, device=png,width=24, height=6,units="cm",dpi=600) 


### plot the  cell-type-specific DHZ

#Drd1
library("tidyverse")

dat_M_CG = read.table("dml_M_D1vsD2_delta20percent1e-4.mergeSS2DS.txt.bedgraph",header=F)
dat_M_CG$DML = ifelse(dat_M_CG$V4>0,"hyper","hypo")
dat_M_CA = read.table("dml.D1vsD2.mCA_delta10percent1e-4.mergeSS2DS.txt.bedgraph",header=F)
dat_M_CA$DML = ifelse(dat_M_CA$V4>0,"hyper","hypo")


chr_lab <- read.table("chr_label.txt")

p1.1 <- ggplot(dat_M_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="red","hypo"="blue"))

p1.2 <- p1.1 +theme_classic(base_size = 5,base_line_size = 0.2)+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,expand = c(0, 0), limits = c(35000000,55000000))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")

library(gridExtra)

p2.1 <- ggplot(dat_M_CA,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="red","hypo"="blue"))

p2.2 <- p2.1 +theme_classic(base_size = 5,base_line_size = 0.2)+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,labels = chr_lab$V1,expand = c(0, 0), limits = c(35000000,55000000))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")


p <- grid.arrange(p1.2, p2.2, nrow = 2, padding = unit(c(0,0,0,0), "cm"))
ggsave("figre1.D1D2.Drd1.test25.png",plot=p, device=png,width=10, height=10,units="cm",dpi=600) 
ggsave("figre1.D1D2.Drd1.test29.svg",plot=p, device=svg,width=8, height=4,units="cm") 
ggsave("figre1.D1D2.Drd1.test30.svg",plot=p, device=svg,width=8, height=4.5,units="cm") 




 

