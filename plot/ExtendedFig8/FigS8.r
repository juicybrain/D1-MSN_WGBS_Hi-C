#
#   Author Yuxiang Li
#

### DMR in loop D1-MSN anchors
library("tidyverse")

dat = read.table(paste0("Loop_stats_DMR.txt"), header=T,sep="\t")


g1 = ggplot(data = dat, mapping = aes(x = factor(group,levels = c("M KO vs WT","F KO vs WT","WT M vs F","KO M vs F")), y = Quantity, fill = Loop)) + geom_bar(stat = 'identity', position = 'stack',size=0.6,width=0.6) + scale_fill_manual(values = c("in"="skyblue","out"="darkgray"))+xlab("") + scale_y_continuous(expand = c(0, 0),limits=c(0, 1000))+ylab("DMR quantity")

g1.2 <- g1 + theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 90, vjust=0.5),
                                                    #  axis.text.x=element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0.2, 'cm'),
                                                     legend.key.height = unit(0.2, 'cm'),
                                                     legend.key.width = unit(0.2, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.background = element_rect(fill = NA,size=0.5,colour = "black"),
                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank())                                               
ggsave("DMR_in_Loop.g4.DMR.svg",plot=g1.2, device="svg",width=9, height=3,unit="cm") 
