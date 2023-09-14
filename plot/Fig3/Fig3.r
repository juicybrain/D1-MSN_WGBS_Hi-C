####
#
#  Author Yuxiang Li
#
###  PCA analysis of the DMR methylation dynamics of all D1-MSN samples 

set.seed(123)
library(dplyr)

wd="D:/D1_WGBS/"
dat_DMR_FM_merge = read.table(paste0(wd, "DMR_meCG.txt"), header=T)
colnames(dat_DMR_FM_merge) = c("pos", "Male_WT_Rep2","Male_WT_Rep3","Male_WT_Rep4","Male_WT_Rep5","Male_WT_Rep6","Male_WT_Rep1","Male_KO_Rep3","Male_KO_Rep4","Male_KO_Rep5","Male_KO_Rep1","Male_KO_Rep2","Female_WT_Rep1","Female_WT_Rep2","Female_WT_Rep3","Female_WT_Rep4","Female_WT_Rep5","Female_KO_Rep3","Female_KO_Rep4","Female_KO_Rep5","Female_KO_Rep6","Female_KO_Rep1","Female_KO_Rep2","M","F","WT", "KO","chr")

dat_all = dat_DMR_FM_merge

sample_group= sapply(strsplit(colnames(dat_all[,2:23]),split="_Re"),"[[",1)

df2 = data.frame(Group=sample_group)
Group=factor(df2$Group,levels=c("Male_WT","Female_WT","Male_KO", "Female_KO"))
col= list(Group=c("Female_KO" ="darkmagenta", "Male_KO"="darkorange2","Female_WT" = "purple", "Male_WT"="orange" ))

dat2 <- scale(t(dat_all[, 2:23]),scale=F, center=T)
dat3 <- t(dat2)

inputMatrix = dat3
probesetvar = apply(inputMatrix,1,var)
ord = order(probesetvar,decreasing=TRUE)
pca = prcomp(t(inputMatrix[ord,]),scale=F)
ss = summary(pca)
pca_dat <- pca$x[,c(1,2)]
pca_dat <- as.data.frame(pca_dat)
pca_dat$name = row.names(pca_dat)

pca_dat$group <- sapply(strsplit(pca_dat$name,split="_Rep"),"[[",1)

g1_DMR <- ggplot(pca_dat,aes(x=pca_dat$PC1,y=pca_dat$PC2,label=pca_dat$name)) + geom_point(aes(color=group,shape=group)) + geom_text(aes(label=""),hjust=-0.5,vjust=-0.5,size=3)+scale_color_manual(values = c("Male_KO"="orange","Male_WT"="orange","Female_WT"="purple","Female_KO"="purple")) + scale_shape_manual(values=c(13,1,13,1))+xlab("PC1")+ylab("PC2")+
  #coord_cartesian(xlim=c(-20,20),ylim = c(-20, 20),expand=F)+
  xlab("PC1 ( 20.3% )")+ylab("PC2 ( 17.8% )") + stat_ellipse(geom = "polygon", linetype = "dotted",size=0.3,alpha=0.05,show.legend=F,aes(color=group),level=0.8)
# + stat_ellipse(type = "norm", linetype = "dotted",size=0.5,alpha=0.5,show.legend=T,aes(color=group),level=0.7)

g1.1_DMR = g1_DMR + theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                                 axis.text.x=element_blank(),
                                                 axis.text.y=element_blank(),
                                                # axis.text.x =element_text(size = 8,angle = 0,vjust=0.5),
                           legend.position ="right",legend.title=element_text(face='italic',size=5),
                           panel.grid.major=element_line(colour=NA),
                   #        legend.title = element_text(size=6),
                           legend.text = element_text(size=5),
                           legend.key.size = unit(0, 'cm'),
                           legend.key.height = unit(0.2, 'cm'),
                           legend.key.width = unit(0, 'cm'),
                                                 panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())
ggsave("D1_DMR.PCA.svg",plot=g1.1_DMR,width = 6,height = 4,unit="cm")

rm(list=ls())

###  Plot the enrichment of TF in DMR by motif analysis 

set.seed(2)
library(tidyverse)
library(reshape2)


dat_tf<- read.table("CG_DMR_motif.txt",header=T)
dat_tf$gene = factor(dat_tf$gene, levels = c("AP-1","Jun-AP1","Atf1","Atf3","JunB","BATF","Bach2","Fosl2","Fra1","Fra2", "Mef2a","Mef2b","Mef2c","Sox2","Sox6","Sox9","Sox10","Sox15","Sox17","RFX","Rfx2","Rfx6","Stat3","STAT5", "ETS","Ets1-distal","Foxa2","Six1","Lhx2","Dlx3","REST-NRSF","NFE2L2","Tbox:Smad","p300"))
dat_tf$Group = factor(dat_tf$Group, levels=c("KO_hypo",  "KO_hyper", "WT_hypo",  "WT_hyper", "F_hypo",   "F_hyper",  "M_hypo",   "M_hyper"))


p_tf2 <- ggplot(dat_tf, aes(x=Group,y=gene, size= logP )) +
  geom_point(aes(fill=DMR), color="black", pch=21, stroke=0.5 ) + 
  scale_size(range=c(1,2)) +
  xlab("") + 
  coord_flip() +
  scale_fill_manual(values=c("hyper"="red", "hypo"="blue"))



p_tf.2 = p_tf2 + theme_bw(base_size = 5)+ theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 90, vjust=0),
                                                      #axis.text.x=element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                    axis.text.y =element_text(size = 5, vjust=0),
                                                     legend.position ="top",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0.2, 'cm'),
                                                     legend.key.height = unit(0.2, 'cm'),
                                                     legend.key.width = unit(0.2, 'cm'),
                                                     legend.box = "horizontal"
                                                ) 
ggsave("CG.TF.svg",plot=p_tf.2, device="svg",width=8, height=4.5,unit="cm") 
rm(list=ls())

###  plot genome-wide DML distribution


library("tidyverse")

dat_M_CG = read.table("dml_D1_Smth_M_KOvsCon_5percent1e-4.txt.bedgraph.txt",header=F)
dat_M_CG$DML = ifelse(dat_M_CG$V4>0,"hyper","hypo")
dat_F_CG = read.table("dml_D1_Smth_F_KOvsCon_5percentP1e-4.mergedSS2DS.txt.bedgraph.txt",header=F)
dat_F_CG$DML = ifelse(dat_F_CG$V4>0,"hyper","hypo")
dat_Ctl_CG = read.table("dml_D1_Smth_Con_MvsF_5percent1e-4.mergeSS2DS.txt.bedgraph.txt",header=F)
dat_Ctl_CG$DML = ifelse(dat_Ctl_CG$V4>0,"hyper","hypo")
dat_Ko_CG = read.table("dml_D1_Smth_KO_Male_vs_Female_5percent1e-4.mergeSS2DS.txt.bedgraph.txt",header=F)
dat_Ko_CG$DML = ifelse(dat_Ko_CG$V4>0,"hyper","hypo")

chr_lab <- read.table("chr_label.txt")
p1.1 <- ggplot(dat_M_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="orange2","hypo"="orange2"))

p1.2 <- p1.1 +theme_classic(base_size = 5,base_line_size = 0.2)+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,expand = c(0, 0), limits = c(0,2462745373))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")

library(gridExtra)

p2.1 <- ggplot(dat_F_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="darkorchid","hypo"="darkorchid"))

p2.2 <- p2.1 +theme_classic(base_size = 5,base_line_size = 0.2)+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,labels = chr_lab$V1,expand = c(0, 0), limits = c(0,2462745373))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")


p3.1 <- ggplot(dat_Ctl_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="darkgreen","hypo"="darkgreen"))

p3.2 <- p3.1 +theme_classic(base_size = 5,base_line_size = 0.2)+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,labels = chr_lab$V1,expand = c(0, 0), limits = c(0,2462745373))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")

p4.1 <- ggplot(dat_Ko_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="cyan3","hypo"="cyan3"))

p4.2 <- p4.1 +theme_classic(base_size = 5,base_line_size = 0.2)+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,labels = chr_lab$V1,expand = c(0, 0), limits = c(0,2462745373))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")

p <- grid.arrange(p1.2, p2.2,p3.2,p4.2, nrow = 4, padding = unit(c(0,0,0,0), "cm"))
ggsave("figre1.test25.png",plot=p, device=png,width=24, height=6,units="cm",dpi=600) 

ggsave("figre1.test29.svg",plot=p, device=svg,width=30, height=9,units="cm") 
ggsave("figre1.test30.svg",plot=p, device=svg,width=30, height=10,units="cm") 
















