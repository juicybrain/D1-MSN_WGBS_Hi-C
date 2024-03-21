#
#	Author: Yuxiang Li
#

### PC1 z-test

# male chr8
dat_MN = read.table(paste0( "M_N.0923.distalN.1M.pca.txt.bedgraph"), header=F)
dat_MP = read.table(paste0( "M_P.0923.distalN.1M.pca.txt.bedgraph"), header=F)
colnames(dat_MN) = c("chr","start","end","PC1N")
colnames(dat_MP) = c("chr","start","end","PC1P")
dat_M = merge(dat_MN,dat_MP, by=c("chr","start","end"), all=F)
write.table(dat_M, "dat_M.1M.PC1.eigenvector.txt", row.names = F, col.names = T, quote=F, sep="\t")

#dat_M=read.table("dat_M.PC1.txt", header=T)
dat_M_8_14 = dat_M[which(dat_M$chr=="chr8"),]
dif_chr8 = c(35000000,36000000,37000000,38000000,39000000,40000000,40000000,41000000,42000000,43000000,50000000,51000000,52000000,53000000)

dat_M_8_14$diff =ifelse( dat_M_8_14$chr == "chr8", ifelse(dat_M_8_14$start %in% dif_chr8,"S8","N"), ifelse(dat_M_8_14$start %in% dif_chr8, "S14", "N" ) ) 

dat_M_8_14 = dat_M_8_14[order(dat_M_8_14$diff, decreasing = F), ]
p3 = ggplot(dat_M_8_14, aes(x=PC1N, y=PC1P,color=diff))+ geom_point( size=0.3) + scale_color_manual(values = c("S8"="blue","N"="gray", "S14"="cyan"))  +
  geom_vline(xintercept=0,lty=4,col="grey",lwd=0.1) + geom_hline(yintercept = 0,lty=4,col="grey",lwd=0.1)+ scale_x_continuous(expand = c(0, 0), limits = c(-0.2,0.2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-0.2,0.2)) + geom_abline(slope = 1, col="black", lwd=0.1, lty="dashed", alpha=0.5)


p3.2 = p3 + theme_classic(base_size = 5)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 0,vjust=0.5),
                                                     axis.text.y =element_text(size = 5,angle = 0,vjust=0.5),axis.title.x = element_blank(), axis.title.y = element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="none",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0, 'cm'),
                                                     legend.key.height = unit(0.5, 'cm'),
                                                     legend.key.width = unit(0, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                    panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
                                                    plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank())
ggsave("D1_M.chr8_diff_hic_PC1.2.svg",plot=p3.2, device="svg",width=3.7, height=3.5,unit="cm") 

# female chr8

dat_MN = read.table(paste0( "F_N.0923.distalN.100k.pca.k4me3.bedgraph"), header=F)
dat_MP = read.table(paste0( "F_P.0923.distalN.100k.pca.k4me3.bedgraph"), header=F)
colnames(dat_MN) = c("chr","start","end","PC1N")
colnames(dat_MP) = c("chr","start","end","PC1P")
dat_M = merge(dat_MN,dat_MP, by=c("chr","start","end"), all=F)
write.table(dat_M, "dat_M.1M.PC1.eigenvector.txt", row.names = F, col.names = T, quote=F, sep="\t")
dat_M_8_14 = dat_M[which(dat_M$chr=="chr8"),]
dif_chr8 = c(35000000,36000000,37000000,38000000,39000000,40000000,41000000,42000000,43000000,50000000,51000000,52000000)
dat_M_8_14$diff =ifelse( dat_M_8_14$chr == "chr8", ifelse(dat_M_8_14$start %in% dif_chr8,"S8","N"), ifelse(dat_M_8_14$start %in% dif_chr8, "S14", "N" ) ) 
dat_M_8_14 = dat_M_8_14[order(dat_M_8_14$diff, decreasing = F), ]
p3 = ggplot(dat_M_8_14, aes(x=PC1N, y=PC1P,color=diff))+ geom_point( size=0.3) + scale_color_manual(values = c("S8"="blue","N"="gray", "S14"="cyan"))  +
  geom_vline(xintercept=0,lty=4,col="grey",lwd=0.1) + geom_hline(yintercept = 0,lty=4,col="grey",lwd=0.1)+ scale_x_continuous(expand = c(0, 0), limits = c(-0.2,0.2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-0.2,0.2)) + geom_abline(slope = 1, col="black", lwd=0.1, lty="dashed", alpha=0.5)
p3.2 = p3 + theme_classic(base_size = 5)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 0,vjust=0.5),
                                                     axis.text.y =element_text(size = 5,angle = 0,vjust=0.5),axis.title.x = element_blank(), axis.title.y = element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="none",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0, 'cm'),
                                                     legend.key.height = unit(0.5, 'cm'),
                                                     legend.key.width = unit(0, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                    panel.background = element_rect(fill = NA,size=0.3,colour = "black"),
                                                    plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank())

ggsave("D1_F_diff_hic_PC1.2.svg",plot=p3.2, device="svg",width=3.7, height=3.5,unit="cm") 


### loops in DIR

#male
library("tidyverse")
library(ggpubr)
dat_MNP = read.table("M.loop.in.DIR.txt", header=T)
dat_FNP = read.table("F.loop.in.DIR.txt", header=T)
 
g1 = ggpaired(dat_MNP, cond1="MN", cond2="MP", fill="condition") + scale_fill_manual(values = c("MN"="orange",  "MP" = "darkorange2")) + theme_classic(base_size = 12)+ ylim(0,150)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                              axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                              #  axis.text.x=element_blank(),
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

g2 = ggpaired(dat_FNP, cond1="FN", cond2="FP", fill="condition") + scale_fill_manual(values = c("FN"="purple",  "FP" = "darkmagenta"))  + theme_classic(base_size = 12)+ ylim(0,150)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                                axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                              #  axis.text.x=element_blank(),
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

ggsave("D1.MN_MP.svg",plot=g1, device="svg",width=6, height=6,unit="cm")
ggsave("D1.FN_FP.svg",plot=g2, device="svg",width=6, height=6,unit="cm")







