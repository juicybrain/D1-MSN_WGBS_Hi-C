####
#
#  Author Yuxiang Li
#


### Plot the DML distribution across the whole genome
library("tidyverse")

dat_M_CG = read.table("dml_D1_Smth_M_WTvsCre_5percent1e-4.txt.bedgraph.new.txt",header=F)
dat_M_CG$DML = ifelse(dat_M_CG$V4>0,"hyper","hypo")
dat_F_CG = read.table("dml_D1_Smth_F_WTvsCre_5percentP1e-4.mergedSS2DS.txt.bedgraph.new.txt",header=F)
dat_F_CG$DML = ifelse(dat_F_CG$V4>0,"hyper","hypo")
dat_Ctl_CG = read.table("dml_D1_Smth_WT_Male_vs_Female_delta5percent1e-4.mergeSS2DS.txt.bedgraph.new.txt",header=F)
dat_Ctl_CG$DML = ifelse(dat_Ctl_CG$V4>0,"hyper","hypo")
dat_Ko_CG = read.table("dml_D1_Smth_KO_Male_vs_Female_delta5percent1e-4.mergeSS2DS.txt.bedgraph.new.txt",header=F)
dat_Ko_CG$DML = ifelse(dat_Ko_CG$V4>0,"hyper","hypo")

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
  scale_x_continuous(breaks=chr_lab$V2,expand = c(0, 0), limits = c(0,2462745373))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")

library(gridExtra)

p2.1 <- ggplot(dat_F_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="red","hypo"="blue"))

p2.2 <- p2.1 +theme_classic(base_size = 5,base_line_size = 0.2)+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,labels = chr_lab$V1,expand = c(0, 0), limits = c(0,2462745373))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")


p3.1 <- ggplot(dat_Ctl_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="red","hypo"="blue"))

p3.2 <- p3.1 +theme_classic(base_size = 5,base_line_size = 0.2)+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"),
  axis.text.x =element_blank(),
  axis.text.y =element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,labels = chr_lab$V1,expand = c(0, 0), limits = c(0,2462745373))+ geom_hline(yintercept = 0, size=0.3, linetype = "dotted")



p4.1 <- ggplot(dat_Ko_CG,aes(x=V2,y=V4,color=DML))+geom_point(size=0.03, alpha=1)+ scale_y_continuous(expand = c(0, 0), limits = c(-1,1))+ylab("")+xlab("")+scale_color_manual(values=c("hyper"="red","hypo"="blue"))

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
ggsave("figre1.test31.png",plot=p, device=png,width=24, height=6,units="cm",dpi=600) 
ggsave("figre1.test32.svg",plot=p, device=svg,width=30, height=8,units="cm") 
ggsave("figre1.test33.svg",plot=p, device=svg,width=30, height=12,units="cm") 

### Plot chr8, chr14 DHZ 

# similar code as above with adjustment of the x-axis (scale_x_continuous(limits = c(XXX,XXX)))

### plot chr8 and chr14 regional mCG level

#firstly, slicing the DHZ into 50 bins, calculate the mCG level of 50 bins along with up- and downstream 10 bins. Then perform the plot and statistic test.
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

#plot se instead of sd
for (spl in list("dat_MN", "dat_MP", "dat_FN", "dat_FP")){
      
      assign(paste0(spl, "_mean"), rowMeans(get(spl)) )
      assign(paste0(spl, "_sd"),  apply(get(spl), 1, sd) )
  }

 dat_M_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "MN_mean"= dat_MN_mean, "MN_sd" = dat_MN_sd/sqrt(6), "MP_mean" = dat_MP_mean, "MP_sd"= dat_MP_sd/sqrt(5))
 dat_F_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "FN_mean"= dat_FN_mean, "FN_sd" = dat_FN_sd/sqrt(5), "FP_mean" = dat_FP_mean, "FP_sd"= dat_FP_sd/sqrt(6))
 
 
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
 
 
 
###### precise p-value #################
 
 
#library("pspearman")
#out <- spearman.test(dat_FN_mean_plot$mean[11:60],  dat_FP_mean_plot$mean[11:60], alternative = "two.sided", approximation="t-distribution")
#out$p.value

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

 
 
 ###### wilcox.test ###################
 
wilcox.test(dat_FN_mean_plot$mean[11:60], dat_FP_mean_plot$mean[11:60],alternative = "greater", paired=T)$p.value
print("p.value < 1.75e-34") 
wilcox.test(dat_MN_mean_plot$mean[11:60], dat_MP_mean_plot$mean[11:60], alternative = "greater",paired=T)$p.value
print("p.value < 1.81e-33") 


 
 ################### chr14 ##########################
 
 
dat = read.table("promoter/chr14.cDMR_promoter_1.bedgraph", header=F)
dat$V5 = 1

for (i in 2:10){
  dat_tmp = read.table(paste0("promoter/chr14.cDMR_promoter_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i
  dat=rbind(dat,dat_tmp)
}

for (i in 1:50){
  dat_tmp = read.table(paste0("gene_body/chr14.cDMR_gene_body_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i+10
  dat=rbind(dat,dat_tmp)
}

for (i in 1:10){
  dat_tmp = read.table(paste0("downstream/chr14.cDMR_downstream_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i+60
  dat=rbind(dat,dat_tmp)
}

dat_s = dat[order(dat$V4, dat$V5),]

write.table(dat_s, "chr14.cDMR.male.scaled.txt", quote=F, row.names = F, sep="\t")


dat_MN = dat_s[, 7:12]
dat_MP = dat_s[, 13:17]
dat_FN = dat_s[, 18:22]
dat_FP = dat_s[, 23:28]


for (spl in list("dat_MN", "dat_MP", "dat_FN", "dat_FP")){
      assign(paste0(spl, "_mean"), rowMeans(get(spl)) )
      assign(paste0(spl, "_sd"),  apply(get(spl), 1, sd) )
  }

 dat_M_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "MN_mean"= dat_MN_mean, "MN_sd" = dat_MN_sd/sqrt(6), "MP_mean" = dat_MP_mean, "MP_sd"= dat_MP_sd/sqrt(5))
 dat_F_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "FN_mean"= dat_FN_mean, "FN_sd" = dat_FN_sd/sqrt(5), "FP_mean" = dat_FP_mean, "FP_sd"= dat_FP_sd/sqrt(6))
 
 
 
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


   ################ MFN ##################
 
 dat_plot_MFN = rbind(dat_MN_mean_plot, dat_FN_mean_plot)
 
 dat_plot_MFN_s = dat_plot_MFN[order(dat_plot_MFN$gene, dat_plot_MFN$region ),]
 
 dat_plot_M_1 = dat_plot_MFN_s[1:140, ]
 g1 = ggplot(dat_plot_M_1, aes(x = region, y = mean, color = group, group = group)) +
  geom_line(size=0.25) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd, fill = group),linetype=0, alpha = 0.2) +
 # geom_point(size=0.1) +
 # theme_minimal() +
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c("MN" = "orange", "FN" = "purple")) +
  scale_fill_manual(values = c("MN" = "orange", "FN" = "purple"))+
  scale_x_continuous(limits=c(0,70),breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) +  scale_y_continuous(limits = c(0.2, 1)) + geom_vline(xintercept = 10, linetype="dashed", color = "black", size=0.25) +
  geom_vline(xintercept = 60, linetype="dashed", color = "black", size=0.25) 
 
 g1.1 = g1  + theme_classic(base_size=5, base_line_size = 0.25)+ theme(legend.position ="right", legend.key.size = unit(0.1, "cm"))
 
 ggsave( "chr14.MFN.svg",plot=g1.1, width = 4,height = 2, units = "cm")
 
 
 
    ########### MFP ##############
 
 dat_plot_MFP = rbind(dat_MP_mean_plot, dat_FP_mean_plot)
 
 dat_plot_MFP_s = dat_plot_MFP[order(dat_plot_MFP$gene, dat_plot_MFP$region ),]
 
 dat_plot_M_1 = dat_plot_MFP_s[1:140, ]
 g1 = ggplot(dat_plot_M_1, aes(x = region, y = mean, color = group, group = group)) +
  geom_line(size=0.25) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd, fill = group),linetype=0, alpha = 0.1) +
 # geom_point(size=0.1) +
 # theme_minimal() +
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c("MP" = "orange", "FP" = "purple")) +
  scale_fill_manual(values = c("MP" = "orange", "FP" = "purple"))+
  scale_x_continuous(limits=c(0,70),breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) +  scale_y_continuous(limits = c(0.2, 1)) + geom_vline(xintercept = 10, linetype="dashed", color = "black", size=0.25) +
  geom_vline(xintercept = 60, linetype="dashed", color = "black", size=0.25) 
 
 g1.1 = g1  + theme_classic(base_size=5, base_line_size = 0.25)+ theme(legend.position ="right", legend.key.size = unit(0.1, "cm"))
 
 ggsave( "chr14.MFP.svg",plot=g1.1, width = 4,height = 2, units = "cm")
 
 
  
###### beta regression ################
  dat_test = data.frame(FN_mean=dat_FN_mean_plot$mean[11:60], MN_mean=dat_MN_mean_plot$mean[11:60], pos=1:50)
library(betareg)
library(reshape)
df = melt(dat_test, id="pos")
colnames(df) = c("pos", "group", "mC")
epsilon <- 0.0001
df$mC <- ifelse(df$mC == 0, df$mC + epsilon, df$mC)
df$mC <- ifelse(df$mC == 1, df$mCG - epsilon, df$mC)

# Fit the model
model_N <- betareg(mC ~ group + pos, data = df)

# Print the summary of the model
summary(model_N)


dat_test = data.frame(FP_mean=dat_FP_mean_plot$mean[11:60], MP_mean=dat_MP_mean_plot$mean[11:60], pos=1:50)
library(betareg)
library(reshape)
df = melt(dat_test, id="pos")
colnames(df) = c("pos", "group", "mC")
epsilon <- 0.0001
df$mC <- ifelse(df$mC == 0, df$mC + epsilon, df$mC)
df$mC <- ifelse(df$mC == 1, df$mC - epsilon, df$mC)

# Fit the model
model_P <- betareg(mC ~ group + pos, data = df)

# Print the summary of the model
summary(model_P)

 
t.test(dat_FN_mean_plot$mean[11:60], dat_MN_mean_plot$mean[11:60],alternative="less", paired=T)$p.value

print(p.value < 4.54e-30)

t.test(dat_FP_mean_plot$mean[11:60], dat_MP_mean_plot$mean[11:60],alternative="less", paired=T)$p.value

 
 
 
 

