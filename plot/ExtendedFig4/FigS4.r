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


### plot the PCA of genome-wide mCG level

## autosome mCG

set.seed(123)
library("tidyverse")
wd="D:/D1_WGBS/final version/fig1_version_4/PCA/CG_100k/"
dat_100k = read.table(paste0(wd, "D1F_WGBS_10.txt"), header=F)
colnames(dat_100k)=c( "chr", "start", "end", unlist(strsplit("D1F_WGBS_10.txt",split=".txt"))[1]

for (spl in c("D1F_WGBS_11","D1F_WGBS_12","D1F_WGBS_1","D1F_WGBS_2","D1F_WGBS_3","D1F_WGBS_4","D1F_WGBS_5","D1F_WGBS_6","D1F_WGBS_7","D1F_WGBS_8","D1M_WGBS_10","D1M_WGBS_11","D1M_WGBS_12","D1M_WGBS_1","D1M_WGBS_2","D1M_WGBS_3","D1M_WGBS_4","D1M_WGBS_5","D1M_WGBS_7","D1M_WGBS_8","D1M_WGBS_9")) {

  dat_tmp= read.table(paste0(wd,spl,".txt"),  header=F)
  colnames(dat_tmp) = c("chr", "start", "end", spl)
  dat_100k = merge(dat_100k, dat_tmp, by=c("chr", "start", "end" ),all=F)  
}

colnames(dat_100k) = c("chr", "start", "end" ,"Female_KO_Rep1","Female_KO_Rep2","Female_WT_Rep1", "Female_WT_Rep2", "Female_WT_Rep3", "Female_KO_Rep3", "Female_WT_Rep4", "Female_KO_Rep4", "Female_KO_Rep5", "Female_KO_Rep6", "Female_WT_Rep5","Male_KO_Rep1","Male_WT_Rep1","Male_KO_Rep2", "Male_WT_Rep2", "Male_WT_Rep3", "Male_WT_Rep4", "Male_WT_Rep5", "Male_WT_Rep6", "Male_KO_Rep3", "Male_KO_Rep4", "Male_KO_Rep5")

rownames(dat_100k) = paste0( dat_100k$chr, dat_100k$start, dat_100k$end)


inputMatrix = dat_100k%>%filter(chr!="chrX")%>% dplyr::select(-c(1,2,3))
probesetvar = apply(inputMatrix,1,var)
ord = order(probesetvar,decreasing=TRUE)[1:12000]
pca = prcomp(t(inputMatrix[ord,]),scale=F)
ss1 = summary(pca)
pca_dat <- pca$x[,c(1,2,3,4)]
pca_dat <- as.data.frame(pca_dat)
pca_dat$name =row.names(pca_dat)
  pca_dat$rep = sapply(strsplit(pca_dat$name,split="_"),"[[",3)

pca_dat$group <- sapply(strsplit(pca_dat$name,split="_Rep"),"[[",1)
pca_dat$TET1 <- sapply(strsplit(pca_dat$name,split="_"),"[[",2)

g1 <- ggplot(pca_dat,aes(x=pca_dat$PC1,y=pca_dat$PC2,label=pca_dat$rep)) + geom_point(aes(color=group,shape=TET1)) + geom_text(aes(label=""),hjust=1,vjust=-1,size=1)+
  #coord_cartesian(xlim=c(-20,20),ylim = c(-20, 20),expand=F) + 
  scale_color_manual(values = c("Male_KO"="orange","Male_WT"="orange","Female_WT"="purple","Female_KO"="purple")) +
  scale_shape_manual(values=c("KO"=13,"WT"=1))+xlab("PC1 ( 16.6% )")+ylab("PC2 ( 9.1% )")


g1.1 = g1 + theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
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
ggsave(paste0(wd, "D1_FM_100kPCA_mCG.svg"),plot=g1.1,width = 6,height = 4,unit="cm")

# all chrX

inputMatrix = dat_100k%>%filter(chr=="chrX")%>% select(-c(1,2,3))
probesetvar = apply(inputMatrix,1,var)
ord = order(probesetvar,decreasing=TRUE)[1:1000]
pca = prcomp(t(inputMatrix[ord,]),scale=F)
ss6 = summary(pca)
pca_dat <- pca$x[,c(1,2,3,4)]

pca_dat <- as.data.frame(pca_dat)
pca_dat$name =row.names(pca_dat)
pca_dat$rep = sapply(strsplit(pca_dat$name,split="_"),"[[",3)
pca_dat$group <- sapply(strsplit(pca_dat$name,split="_Rep"),"[[",1)
pca_dat$TET1 <- sapply(strsplit(pca_dat$name,split="_"),"[[",2)

g1 <-  ggplot(pca_dat,aes(x=pca_dat$PC1,y=pca_dat$PC2,label=pca_dat$rep)) + geom_point(aes(color=group,shape=TET1)) + geom_text(aes(label=""),hjust=1,vjust=-1,size=1)+
  #coord_cartesian(xlim=c(-20,20),ylim = c(-20, 20),expand=F) + 
  scale_color_manual(values = c("Male_KO"="orange","Male_WT"="orange","Female_WT"="purple","Female_KO"="purple")) +
  scale_shape_manual(values=c("KO"=13,"WT"=1))+xlab("PC1 ( 95% )")+ylab("PC2 ( 0.7% )")

g1.1 = g1 + theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
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
ggsave(paste0(wd, "D1_FM_100kPCA_mCG.chrX.svg"),plot=g1.1,width = 6,height = 4,unit="cm")

