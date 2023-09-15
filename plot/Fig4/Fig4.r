#
#  Author Yuxiang Li
#

### plot the SCC score of pair-wised comparison of Hi-C samples
library("tidyverse")

df_res=data.frame("pair"=NA,"SCC"=NA,"chr"=NA)
for (chr in c(1:19, "X")){
  dat_tmp = read.table(paste0("./HiC-cor/chr",chr,".cor.matrix.hicRep.txt"), header=T)
  row.names(dat_tmp) = colnames(dat_tmp)   
  d_tmp = c()
  m_tmp = as.matrix(dat_tmp)
  
for (i in 1:dim(m_tmp)[1]){
  for(j in 1:dim(m_tmp)[1]){
    if (i>j){
     d_tmp = rbind(d_tmp, c(paste0(strsplit(rownames(m_tmp)[i],split="_R")[[1]][1],"_" , strsplit(colnames(m_tmp)[j], split = "_R")[[1]][1] ),  m_tmp[i,j], paste0("chr",chr) ) )  
    }
  }
}
  
  df_tmp=as.data.frame(d_tmp)
  colnames(df_tmp)=c("pair", "SCC", "chr")

  assign(paste0("d",chr), df_tmp)
  df_res=rbind(df_res, df_tmp)
  
}

df_backup=df_res

df_res=df_backup[-1,]

df_res$sex= paste0(sapply(strsplit(df_res$pair,split="_"),"[[",1) , sapply(strsplit(df_res$pair,split="_"),"[[",3) )

df_res$tet1= paste0(sapply(strsplit(df_res$pair,split="_"),"[[",2) ,  sapply(strsplit(df_res$pair,split="_"),"[[",4) )


df_res$sex=str_replace_all(df_res$sex,"FM","MF")
df_res$tet1=str_replace_all(df_res$tet1,"KOWT","WTKO")

df_res$chr=factor(df_res$chr,levels = c("chr1", "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX"))
p5 = ggplot(df_res, aes(x=chr, y=as.numeric(SCC))) + geom_violin(width=1,alpha=0.5, size=0.2)+
 # geom_jitter(size=0.7, alpha=1, position=position_jitter(0.3),aes(color=tet1, shape=sex))+
 geom_point( size=2, position = position_jitter(seed = 1, width = 0.3) ,alpha=1, aes(color=tet1, shape=sex) ) +
  scale_color_manual(values = c("KOKO"="blue","WTWT"="blue","WTKO"="red"))+xlab("")+
  ylab("Hi-C data Pair-wise similarity")+ ylim(0.92,1) + scale_shape_manual(values=c("FF"=4,"MM"=4,"MF"=3 ))

p5.1 = p5 + theme_classic(base_size = 6)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                               axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                               axis.line.x.bottom=element_line(size=0.2),
                                               # axis.text.x=element_blank(),
                                                 legend.position ="right",
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm'),
                                                 panel.background = element_rect(fill = NA,size=0.2,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank()) 

ggsave(filename = "./Hic.cor.svg",plot=p5.1,height = 10,width = 18,unit="cm")







