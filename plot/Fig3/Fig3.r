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

###  Plot genome-wide DML distribution


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


###  Plot mCG in the hotspot on chr8  


############################ chr8 ###############################################
dat = read.table("promoter/chr8.cDMR_promoter_1.bedgraph", header=F)
dat$V5 = 1
for (i in 2:10){
  dat_tmp = read.table(paste0("promoter/chr8.cDMR_promoter_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i
  dat=rbind(dat,dat_tmp)
}

for (i in 1:50){
  dat_tmp = read.table(paste0("gene_body/chr8.cDMR_gene_body_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i+10
  dat=rbind(dat,dat_tmp)
}

for (i in 1:10){
  dat_tmp = read.table(paste0("downstream/chr8.cDMR_downstream_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i+60
  dat=rbind(dat,dat_tmp)
}

dat_s = dat[order(dat$V4, dat$V5),]

write.table(dat_s, "chr8.cDMR.scaled.txt", quote=F, row.names = F, sep="\t")

#library(reshape)
#dat_test = head(dat_s)

#dat_test = dat_test[,-c(1,2,3,6)]
#colnames(dat_test) = c("gene","region",rep("MN",6), rep("MP",5), rep("FN", 5), rep("FP",6))

dat_MN = dat_s[, 7:12]
dat_MP = dat_s[, 13:17]
dat_FN = dat_s[, 18:22]
dat_FP = dat_s[, 23:28]


for (spl in list("dat_MN", "dat_MP", "dat_FN", "dat_FP")){
      
      assign(paste0(spl, "_mean"), rowMeans(get(spl)) )
      assign(paste0(spl, "_sd"),  apply(get(spl), 1, sd) )
  }

 dat_M_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "MN_mean"= dat_MN_mean, "MN_sd" = dat_MN_sd, "MP_mean" = dat_MP_mean, "MP_sd"= dat_MP_sd)
 dat_F_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "FN_mean"= dat_FN_mean, "FN_sd" = dat_FN_sd, "FP_mean" = dat_FP_mean, "FP_sd"= dat_FP_sd)
 
 
 dat_MN_mean_plot = dat_M_mean_sd[,1:4]
 colnames(dat_MN_mean_plot) = c("gene", "region", "mean", "sd") 
 dat_MN_mean_plot$group = "MN"
 
 dat_MP_mean_plot = dat_M_mean_sd[,c(1,2,5,6)]
 colnames(dat_MP_mean_plot) = c("gene", "region", "mean", "sd")
 dat_MP_mean_plot$group = "MP"

 dat_FN_mean_plot = dat_F_mean_sd[,1:4]
 colnames(dat_FN_mean_plot) = c("gene", "region", "mean", "sd") 
 dat_FN_mean_plot$group = "FN"
 
 dat_FP_mean_plot = dat_F_mean_sd[,c(1,2,5,6)]
 colnames(dat_FP_mean_plot) = c("gene", "region", "mean", "sd")
 dat_FP_mean_plot$group = "FP"
 
 
 dat_plot_M = rbind(dat_MN_mean_plot, dat_MP_mean_plot)
 
 dat_plot_M_s = dat_plot_M[order(dat_plot_M$gene, dat_plot_M$region ),]
 
 dat_plot_M_1 = dat_plot_M_s[1:140, ]
 
 
 g1 = ggplot(dat_plot_M_1, aes(x = region, y = mean, color = group, group = group)) +
  geom_line(size=0.25) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd, fill = group),linetype=0, alpha = 0.2) +
 # geom_point(size=0.1) +
 # theme_minimal() +
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c("MN" = "blue", "MP" = "red")) +
  scale_fill_manual(values = c("MN" = "blue", "MP" = "red"))+
  scale_x_continuous(limits=c(0,70),breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) +  scale_y_continuous(limits = c(0.5, 1)) + geom_vline(xintercept = 10, linetype="dashed", color = "black", size=0.25) +
  geom_vline(xintercept = 60, linetype="dashed", color = "black", size=0.25) 
 
 g1.1 = g1  + theme_classic(base_size=5, base_line_size = 0.25)+ theme(legend.position ="right", legend.key.size = unit(0.1, "cm"))
 
 ggsave( "chr.male.8.8.svg",plot=g1.1, width = 4,height = 2, units = "cm")
 
 ################### female ##################
 
 dat_plot_F = rbind(dat_FN_mean_plot, dat_FP_mean_plot)
 
 dat_plot_F_s = dat_plot_F[order(dat_plot_F$gene, dat_plot_F$region ),]
 
 dat_plot_F_1 = dat_plot_F_s[1:140, ]
 
 
 g1 = ggplot(dat_plot_F_1, aes(x = region, y = mean, color = group, group = group)) +
  geom_line(size=0.25) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd, fill = group),linetype=0, alpha = 0.2) +
 # geom_point(size=0.1) +
 # theme_minimal() +
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c("FN" = "blue", "FP" = "red")) +
  scale_fill_manual(values = c("FN" = "blue", "FP" = "red"))+
  scale_x_continuous(limits=c(0,70),breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) +  scale_y_continuous(limits = c(0.5, 1)) + geom_vline(xintercept = 10, linetype="dashed", color = "black", size=0.25) +
  geom_vline(xintercept = 60, linetype="dashed", color = "black", size=0.25) 
 
 g1.1 = g1  + theme_classic(base_size=5, base_line_size = 0.25)+ theme(legend.position ="right", legend.key.size = unit(0.1, "cm"))
 
 ggsave( "chr.female.8.8.svg",plot=g1.1, width = 4,height = 2, units = "cm")
 
 
##################### beta regression ############################
 

dat_test = data.frame(MN_mean=dat_MN_mean_plot$mean[11:60], MP_mean=dat_MP_mean_plot$mean[11:60], pos=1:50)
library(betareg)
library(reshape)
df = melt(dat_test, id="pos")
colnames(df) = c("pos", "group", "mC")
epsilon <- 0.0001
df$mC <- ifelse(df$mC == 0, df$mC + epsilon, df$mC)
df$mC <- ifelse(df$mC == 1, df$mC - epsilon, df$mC)

# Fit the model
model_M <- betareg(mC ~ group + pos, data = df)

# Print the summary of the model
summary(model_M)


dat_test = data.frame(FN_mean=dat_FN_mean_plot$mean[11:60], FP_mean=dat_FP_mean_plot$mean[11:60], pos=1:50)
library(betareg)
library(reshape)
df = melt(dat_test, id="pos")
colnames(df) = c("pos", "group", "mC")
epsilon <- 0.0001
df$mC <- ifelse(df$mC == 0, df$mC + epsilon, df$mC)
df$mC <- ifelse(df$mC == 1, df$mC - epsilon, df$mC)

# Fit the model
model_F <- betareg(mC ~ group + pos, data = df)

# Print the summary of the model
summary(model_F)


###  Plot conversion-rate corrected mCH in the hotspot on chr8 


############################ chr8 ###############################################

library("tidyverse")
dat_CA = read.table("mCA.corrected.txt", header=T)
dat_CT = read.table("mCT.corrected.txt", header=T)
dat_CC = read.table("mCC.corrected.txt", header=T)


dat_c = dat_CA[, 51:94] + dat_CT[, 51:94] + dat_CC[,-1]

dat_s = dat_c[,c(T,F)] / dat_c[,c(F,T)]

dat_s = cbind(dat_CA[, 1:6], dat_s)

#library(reshape)
#dat_test = head(dat_s)

#dat_test = dat_test[,-c(1,2,3,6)]
#colnames(dat_test) = c("gene","region",rep("MN",6), rep("MP",5), rep("FN", 5), rep("FP",6))

dat_MN = dat_s[, 7:12]
dat_MP = dat_s[, 13:17]
dat_FN = dat_s[, 18:22]
dat_FP = dat_s[, 23:28]


for (spl in list("dat_MN", "dat_MP", "dat_FN", "dat_FP")){
      
      assign(paste0(spl, "_mean"), rowMeans(get(spl)) )
      assign(paste0(spl, "_sd"),  apply(get(spl), 1, sd) )  
      
  }

 dat_M_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "MN_mean"= dat_MN_mean, "MN_sd" = dat_MN_sd, "MP_mean" = dat_MP_mean, "MP_sd"= dat_MP_sd)
 dat_F_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "FN_mean"= dat_FN_mean, "FN_sd" = dat_FN_sd, "FP_mean" = dat_FP_mean, "FP_sd"= dat_FP_sd)
 
 
 dat_MN_mean_plot = dat_M_mean_sd[,1:4]
 colnames(dat_MN_mean_plot) = c("gene", "region", "mean", "sd") 
 dat_MN_mean_plot$group = "MN"
 
 dat_MP_mean_plot = dat_M_mean_sd[,c(1,2,5,6)]
 colnames(dat_MP_mean_plot) = c("gene", "region", "mean", "sd")
 dat_MP_mean_plot$group = "MP"

 dat_FN_mean_plot = dat_F_mean_sd[,1:4]
 colnames(dat_FN_mean_plot) = c("gene", "region", "mean", "sd") 
 dat_FN_mean_plot$group = "FN"
 
 dat_FP_mean_plot = dat_F_mean_sd[,c(1,2,5,6)]
 colnames(dat_FP_mean_plot) = c("gene", "region", "mean", "sd")
 dat_FP_mean_plot$group = "FP"
 
 
 dat_plot_M = rbind(dat_MN_mean_plot, dat_MP_mean_plot)
 
 dat_plot_M_s = dat_plot_M[order(dat_plot_M$gene, dat_plot_M$region ),]
 
 dat_plot_M_1 = dat_plot_M_s[141:280, ]
 
 
 g1 = ggplot(dat_plot_M_1, aes(x = region, y = mean, color = group, group = group)) +
  geom_line(size=0.25) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd, fill = group),linetype=0, alpha = 0.2) +
 # geom_point(size=0.1) +
 # theme_minimal() +
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c("MN" = "blue", "MP" = "red")) +
  scale_fill_manual(values = c("MN" = "blue", "MP" = "red"))+
  scale_x_continuous(limits=c(0,70),breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) +  scale_y_continuous(limits = c(0.005, 0.03)) + geom_vline(xintercept = 10, linetype="dashed", color = "black", size=0.25) +
  geom_vline(xintercept = 60, linetype="dashed", color = "black", size=0.25) 
 
 g1.1 = g1  + theme_classic(base_size=5, base_line_size = 0.25)+ theme(legend.position ="right", legend.key.size = unit(0.1, "cm"))
 
 ggsave( "chr.male.8.svg",plot=g1.1, width = 4,height = 2, units = "cm")
 
 ################### female ##################
 
 dat_plot_F = rbind(dat_FN_mean_plot, dat_FP_mean_plot)
 
 dat_plot_F_s = dat_plot_F[order(dat_plot_F$gene, dat_plot_F$region ),]
 
 dat_plot_F_1 = dat_plot_F_s[141:280, ]
 
 
 g1 = ggplot(dat_plot_F_1, aes(x = region, y = mean, color = group, group = group)) +
  geom_line(size=0.25) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd, fill = group),linetype=0, alpha = 0.2) +
 # geom_point(size=0.1) +
 # theme_minimal() +
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c("FN" = "blue", "FP" = "red")) +
  scale_fill_manual(values = c("FN" = "blue", "FP" = "red"))+
  scale_x_continuous(limits=c(0,70),breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) +  scale_y_continuous(limits = c(0.005, 0.03)) + geom_vline(xintercept = 10, linetype="dashed", color = "black", size=0.25) +
  geom_vline(xintercept = 60, linetype="dashed", color = "black", size=0.25) 
 
 g1.1 = g1  + theme_classic(base_size=5, base_line_size = 0.25)+ theme(legend.position ="right", legend.key.size = unit(0.1, "cm"))
 
 ggsave( "chr.female.8.8.svg",plot=g1.1, width = 4,height = 2, units = "cm")
 
 
 
 

#chr8
  dat_test = data.frame(MP_mean=dat_MP_mean_plot$mean[81:130], MN_mean=dat_MN_mean_plot$mean[81:130], pos=1:50)
library(betareg)
library(reshape)
df = melt(dat_test, id="pos")
colnames(df) = c("pos", "group", "mCH")
epsilon <- 0.0001
df$mCH <- ifelse(df$mCH == 0, df$mCH + epsilon, df$mCH)
df$mCH <- ifelse(df$mCH == 1, df$mCH - epsilon, df$mCH)

# Fit the model
model_M <- betareg(mCH ~ group +pos, data = df)

# Print the summary of the model
summary(model_M)

#chr8
dat_test = data.frame(FP_mean=dat_FP_mean_plot$mean[81:130], FN_mean=dat_FN_mean_plot$mean[81:130], pos=1:50)
library(betareg)
library(reshape)
df = melt(dat_test, id="pos")
colnames(df) = c("pos", "group", "mCH")
epsilon <- 0.0001
df$mCH <- ifelse(df$mCH == 0, df$mCH + epsilon, df$mCH)
df$mCH <- ifelse(df$mCH == 1, df$mCH - epsilon, df$mCH)

# Fit the model
model_F <- betareg(mCH ~ group + pos, data = df)

# Print the summary of the model
summary(model_F)











