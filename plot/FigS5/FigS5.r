#
#   Author Yuxiang Li
#
### plot qPCR result of D1 D2 marker genes expression in D1-tdTomato cells
library("tidyverse")

dat_plot = read.table("D1.qPCR.txt",header=T)
df_mRNA_mean = dat_plot%>%group_by(gene)%>%summarise(mRNA_mean=mean(mRNA) )
df_mRNA_se = dat_plot%>%group_by(gene)%>%summarise(mRNA_se=sd(mRNA)/sqrt(n()) )
dat_plot = merge(df_mRNA_mean,df_mRNA_se,by=c("gene"))
dat_plot$gene=factor(dat_plot$gene, levels=c("Drd1","Pdyn","Drd2","Gfap","Gpr6","Hprt1"))

# Finished bar plot
# scale_color_manual(values=c("CG"="black","nonCG"="gray"))  
p1<- ggplot(dat_plot, aes(x=gene, y=mRNA_mean, fill=gene)) + 
  geom_col(stat="identity",
           position=position_dodge(),size=0.4,width=0.4, col="black") +
  geom_errorbar(aes(ymin=mRNA_mean-mRNA_se, ymax= mRNA_mean+mRNA_se),size=0.2, width=0.2,col="black",
                 position=position_dodge(.9))+ scale_fill_manual(values = c("Drd1" ="coral", "Pdyn"="coral","Drd2" = "green", "Gprc6"="orange", "Gfap"="darkolivegreen2", "Hprt1"="cadetblue2")) + xlab("")+ylab("mRNA level (normalized to Gapdh)") + coord_cartesian(ylim = c(0,90),expand=F)+expand_limits(x =c(0.5,6.5))


p1.2 <- p1 + theme_classic(base_size = 6)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                              #  axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                                axis.text.x=element_blank(),
                                                 legend.position ="",
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm'),
                                                 #panel.background = element_rect(fill = NA,size=0.5,colour = "black"),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())

#+ scale_y_continuous(expand = c(0, 0), limits = c(0.5,1))
ggsave("D1-tdTomato.qPCR.svg",plot=p1.2, device="svg",width=4, height=3,unit="cm")


### plot the qPCR validation reslult of Tet1 KO in D1-MSN

library("tidyverse")
dat_plot = read.table("Tet1qPCR_1.txt",header=T)
df_mRNA_mean = dat_plot%>%group_by(group)%>%summarise(mRNA_mean=mean(FC) )
df_mRNA_se = dat_plot %>% group_by(group) %>% summarise(mRNA_se = sd(FC) / sqrt(n()))
dat_plot = merge(df_mRNA_mean,df_mRNA_se,by=c("group"))
dat_plot$group=factor(dat_plot$group, levels=c("M_Ctl","M_KO","F_Ctl","F_KO"))
dat_plot_se = merge(df_mRNA_mean,df_mRNA_se,by=c("group"))
dat_plot_se$group=factor(dat_plot_se$group, levels=c("M_Ctl","M_KO","F_Ctl","F_KO"))

p1<- ggplot(dat_plot, aes(x=group, y=mRNA_mean, fill=group)) + 
  geom_col(stat="identity",
           position=position_dodge(),size=0.4,width=0.4, col="black") +
  geom_errorbar(aes(ymin=mRNA_mean-mRNA_se, ymax= mRNA_mean+mRNA_se),size=0.2, width=0.2,col="black",
                 position=position_dodge(.9))+ scale_fill_manual(values = c("M_Ctl"="orange", "M_KO"="darkorange2" , "F_Ctl" = "purple", "F_KO" ="darkmagenta" ))  + xlab("")+ylab("mRNA level (fold change)") + coord_cartesian(ylim = c(0,2),expand=F)+expand_limits(x =c(0.5,4.5))

p1.2 <- p1 + theme_classic(base_size = 6)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                              #  axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                                axis.text.x=element_blank(),
                                                 legend.position ="",
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm'),
                                             #    panel.background = element_rect(fill = NA,size=0.5,colour = "black"),
                                            #     plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())

ggsave("D1.Tet1KO.svg",plot=p1.2, device="svg",width=3, height=3,unit="cm")

### plot global mCG level in all D1-MSN replicates
# only autosomes
library(tidyverse)
library(Unicode)
data_mCG = read.table("global_meCG_level.auto.txt",header=T)
data_mCG$Group = factor(data_mCG$Group,levels = c("male_WT","male_KO","female_WT","female_KO"))
italic_p <- u_char_inspect(u_char_from_name("MATHEMATICAL ITALIC SMALL P"))["Char"]
p1 <- ggplot(data_mCG, aes(x=Group, y=100*mCG_level,color=Sex)) + geom_violin()+geom_jitter(shape=16, size=0.5,position=position_jitter(0.1))+
  scale_color_manual(values = c("male"="orange","female"="purple"))+xlab("")+ylab("global meCG level (%)")+coord_cartesian(xlim=c(0.5,4.5),ylim = c(78, 81),expand=F) +     stat_summary(fun = "mean", geom = "crossbar", width = 0.5, size=0.1, colour = "black")

p1.1=p1 + theme_classic(base_size = 5)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                                 #  axis.text.x =element_text(size = 6,angle = 45,vjust=0.5),
                                                axis.text.x=element_blank(),
                                                 legend.position ="right",
                                              legend.title=element_text(size=5), panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0, 'cm'),
                                                 panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())

ggsave("global_meCG.auto.1.svg",plot=p1.1,width=4, height=4,unit="cm")

# all chromosomes
library(tidyverse)
library(Unicode)
data_mCG = read.table("global_meCG_level.txt",header=T)
data_mCG$Group = factor(data_mCG$Group,levels = c("male_WT","male_KO","female_WT","female_KO"))
italic_p <- u_char_inspect(u_char_from_name("MATHEMATICAL ITALIC SMALL P"))["Char"]
p1 <- ggplot(data_mCG, aes(x=Group, y=100*mCG_level,color=Sex)) + geom_violin()+geom_jitter(shape=16, size=0.5,position=position_jitter(0.1))+
  scale_color_manual(values = c("male"="orange","female"="purple"))+xlab("")+ylab("global meCG level (%)")+coord_cartesian(xlim=c(0.5,4.5),ylim = c(78, 81),expand=F) + stat_summary(fun = "mean", geom = "crossbar", width = 0.5, size=0.1, colour = "black")

p1.1=p1 + theme_classic(base_size = 5)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                                 #  axis.text.x =element_text(size = 6,angle = 45,vjust=0.5),
                                                axis.text.x=element_blank(),
                                                 legend.position ="right",
                                              legend.title=element_text(size=5), panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0, 'cm'),
                                                 panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())

ggsave("global_meCG.all.svg",plot=p1.1,width=4, height=4,unit="cm")





