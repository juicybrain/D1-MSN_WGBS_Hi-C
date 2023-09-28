#
# Author Yuxiang Li
#
### K-means clustering of D1-MSN DMR and plot


set.seed(123)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

wd="D:/D1_WGBS/final version/fig1_version_4/DMR_cluster_08.05/"
dat_DMR_FM_merge = read.table(paste0(wd, "DMR_meCG.txt"), header=T)
colnames(dat_DMR_FM_merge) = c("pos", "Male_WT_Rep2","Male_WT_Rep3","Male_WT_Rep4","Male_WT_Rep5","Male_WT_Rep6","Male_WT_Rep1","Male_KO_Rep3","Male_KO_Rep4","Male_KO_Rep5","Male_KO_Rep1","Male_KO_Rep2","Female_WT_Rep1","Female_WT_Rep2","Female_WT_Rep3","Female_WT_Rep4","Female_WT_Rep5","Female_KO_Rep3","Female_KO_Rep4","Female_KO_Rep5","Female_KO_Rep6","Female_KO_Rep1","Female_KO_Rep2","M","F","WT", "KO","chr")
dat_all = dat_DMR_FM_merge
sample_group= sapply(strsplit(colnames(dat_all[,2:23]),split="_Re"),"[[",1)
df2 = data.frame(Group=sample_group)
Group=factor(df2$Group,levels=c("Male_WT","Female_WT","Male_KO", "Female_KO"))
col= list(Group=c("Female_KO" ="darkmagenta", "Male_KO"="darkorange2","Female_WT" = "purple", "Male_WT"="orange" ))
ha = HeatmapAnnotation(df = df2, col = col)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# 74 distinctive colors in R:
distinctive_colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
dat2 <- scale(t(dat_all[, 2:23]),scale=F, center=T)
dat3 <- t(dat2)
row.names(dat3) = dat_all$pos
write.table(dat3,"D1_DMR.new.07.19.scale.txt",quote=F,sep="\t")
dat3[is.na(dat3)] = 0


svg('D1_DMR.center.noscale_clr12_04.24.centroid.svg',width = 8,height = 8)
h1 <- Heatmap(dat3,name="normalized MeCG",
            row_km = 12,
            row_gap = unit(1, "mm"),
            column_gap = unit(1, "mm"),
            row_km_repeats = 10,
            top_annotation = ha,
            col=colorRamp2(c(-0.5,0,0.5), c("blue", "white", "red")),
            clustering_distance_rows =  "euclidean",
            clustering_method_rows = "centroid",
            show_row_names = FALSE,
            row_title = "%s",
            border = TRUE,
            show_row_dend = F,
            use_raster = TRUE,
            column_split=Group,
            cluster_column_slices = FALSE,
            column_order=c("Male_WT_Rep6","Male_WT_Rep2","Male_WT_Rep3","Male_WT_Rep4","Male_WT_Rep1","Male_WT_Rep5","Female_WT_Rep1","Female_WT_Rep4","Female_WT_Rep5","Female_WT_Rep3","Female_WT_Rep2", "Male_KO_Rep5","Male_KO_Rep4","Male_KO_Rep3","Male_KO_Rep1","Male_KO_Rep2", "Female_KO_Rep1","Female_KO_Rep2","Female_KO_Rep4","Female_KO_Rep6","Female_KO_Rep5","Female_KO_Rep3")      
)

h1
dev.off()

dat_test = dat_all[,24:27]
df4=data.frame(dat_test)
colrow2 = list("M"=c("M_hyper"= "black", "M_hypo"="sky blue"), "F"= c("F_hyper"="black","F_hypo"="sky blue"), "WT"=c("WT_hyper"="black","WT_hypo"="sky blue"), "KO"=c("KO_hyper"= "black", "KO_hypo"="sky blue")   )
df5 = data.frame(dat_all$chr)
df5$dat_all.chr = factor(df5$dat_all.chr, levels=c(paste0("chr",1:19)))

colrow3 = list("dat_all.chr"=setNames(  c("black", "darkred", "orange", "gold", "darkgreen", "cyan", "blue", "purple", "coral", "chartreuse1", "burlywood2", "cornsilk", "darkgray", "darkblue", "deeppink", "darkorange4", "cyan4", "red", "darkorchid4") , levels(df5$dat_all.chr) ))

ha3 = HeatmapAnnotation(df = df4, which = 'row', col = colrow2,na_col = "white")
ha4 = HeatmapAnnotation(df=df5, which='row',col=colrow3,na_col="white")
svg('D1_MF_DMR_clr13_scale.test_13_8_fontchange_1.centroid.center.4.svg',width = 8,height = 8)
h3=draw(h1+ha3+ha4 )
dev.off()

    rc1.list <-  row_order(h3) 
    lapply(rc1.list, function(x) length(x))  #check/confirm size gene clusters
    library(magrittr) # needed to load the pipe function '%>%'
    
            clu_df <- lapply(names(rc1.list), function(i){
            #clu_df <- lapply(c(1:13), function(i){
            out <- data.frame(DMR = rownames(dat3[rc1.list[[i]],]),
            Cluster = paste0("cluster", i),
            stringsAsFactors = FALSE)
            return(out)
            }) %>% do.call(rbind, .)
            
 write.table(clu_df,paste0("D1_FM_DMR_2022.h3.clr13.centroid.center.07.03.2023.txt"),row.names=F,quote=F,sep="\t")
 
 #### elbow test
 library("purrr")
 wss <- function(k){
   kmeans(dat3,k,nstart = 10)$tot.withinss
 }
 k.values = 1:15
 wss_values <- map_dbl(k.values, wss)
 
 g1 = plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")

 
 library(factoextra)
 fviz_nbclust(dat3, kmeans, method = "wss")
 fviz_nbclust(dat3, kmeans, method = "wss", k.max = 20, iter.max=15)
 fviz_nbclust


rm(list=ls())   

           
### Plot the survey DML/DMR in/out DHR
library("tidyverse")
dat = read.table(paste0("DHR_stats_DMR.txt"), header=T,sep="\t")
g1 = ggplot(data = dat, mapping = aes(x = factor(group,levels = c("M KO vs WT hyper","M KO vs WT hypo","F KO vs WT hyper","F KO vs WT hypo","WT M vs F hyper","WT M vs F hypo","KO M vs F hyper","KO M vs F hypo")), y = Quantity, fill = DHR)) + geom_bar(stat = 'identity', position = 'stack',size=0.6,width=0.6) + scale_fill_manual(values = c("in"="skyblue","out"="darkgray"))+xlab("") + scale_y_continuous(expand = c(0, 0),limits=c(0, 800))+ylab("DMR quantity")
g1.2 <- g1 + theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     #axis.text.x =element_text(size = 5,angle = 90, vjust=0.5),
                                                      axis.text.x=element_blank(),
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
                                                
ggsave("DML_in_DHR.g4.DMR.svg",plot=g1.2, device="svg",width=6, height=3,unit="cm") 

dat = read.table(paste0(wd,"DHR_stats_DML.1.txt"), header=T,sep="\t")
g2 = ggplot(data = dat, mapping = aes(x = factor(group,levels = c("M KO vs WT hyper","M KO vs WT hypo","F KO vs WT hyper","F KO vs WT hypo","WT M vs F hyper","WT M vs F hypo","KO M vs F hyper","KO M vs F hypo")), y = Quantity, fill = DHR)) + geom_bar(stat = 'identity', position = 'stack',size=0.6,width=0.6) + scale_fill_manual(values = c("in"="skyblue","out"="darkgray"))+xlab("") + scale_y_continuous(expand = c(0, 0),limits=c(0, 7000)) + ylab("DML quantity")
g2.2 <- g2 + theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     #axis.text.x =element_text(size = 5,angle = 90, vjust=0.5),
                                                      axis.text.x=element_blank(),
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
                                                 
ggsave("DML_in_DHR.g4.DML.svg",plot=g2.2, device="svg",width=6, height=3,unit="cm") 


library(ggpubr)
theme_set(theme_pubr())
library("gridExtra")
figure1 <- grid.arrange( g1.2, g2.2,layout_matrix = matrix(c(2,1), nrow =2))
ggsave("DML_DMR.test3.svg",plot=figure1, device="svg",width=6, height=4,unit="cm")



### plot size of DHR

M_hyper = read.table("D1_M_WTvsCre.hyper.cDMR.bed",header=F)
M_hypo = read.table("D1_M_WTvsCre.hypo.cDMR.bed",header=F)
F_hyper = read.table("D1_F_WTvsCre.hyper.cDMR.bed",header=F)
F_hypo = read.table("D1_F_WTvsCre.hypo.cDMR.bed",header=F)
Ctl_hyper = read.table("D1_WT_Male_vs_Female.hyper.cDMR.bed", header=F)
Ctl_hypo = read.table("D1_WT_Male_vs_Female.hypo.cDMR.bed", header=F)
Ko_hyper = read.table("D1_KO_Male_vs_Female.hyper.cDMR.bed", header=F)
Ko_hypo = read.table("D1_KO_Male_vs_Female.hypo.cDMR.bed", header=F)
groups=c()
size=c()

for (spl in c("M_hyper", "M_hypo", "F_hyper", "F_hypo", "Ctl_hyper", "Ctl_hypo", "Ko_hyper", "Ko_hypo")){
  assign(paste(spl,"size",sep = "_"),  get(spl)$V3- get(spl)$V2 )
  groups =c(groups, rep(spl,length(get(paste(spl,"size",sep = "_")))))
  size=c(size,unlist(get(paste(spl,"size",sep = "_"))))
}
df = data.frame(DHR_size=size, group=groups)
df$group = factor(df$group, levels=c("M_hyper", "M_hypo", "F_hyper", "F_hypo", "Ctl_hyper", "Ctl_hypo", "Ko_hyper", "Ko_hypo"))

library(ggplot2)

g3=ggplot(df, aes(x = group, fill=df$group ,y = log10(DHR_size))) +
  geom_boxplot(width=0.3) +
  stat_summary(fun = function(x) quantile(x, 0.75),
               geom = "point",
               shape = 5,
               size = 0,
               color = "black",
               show.legend = FALSE) +
  ylab("Values") +
  xlab("Group") + scale_fill_manual(values = c("M_hyper"="orange", "M_hypo"="orange", "F_hyper"="purple", "F_hypo"="purple", "Ctl_hyper"="darkgreen", "Ctl_hypo"="darkgreen", "Ko_hyper"="cyan", "Ko_hypo"="cyan"))+ ylim(c(3,7))

g3.1 = g3 + theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     #axis.text.x =element_text(size = 5,angle = 90, vjust=0.5),
                                                      axis.text.x=element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0.2, 'cm'),
                                                     legend.key.height = unit(0.2, 'cm'),
                                                     legend.key.width = unit(0.2, 'cm'),
                                                     legend.box = "top",
                                                    # panel.grid.major=element_line(colour=NA),
                                                    # panel.background = element_rect(fill = NA,size=0.5,colour = "black"),
                                                    # plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank()) 
                                                 
ggsave("DHRs.length.svg",plot=g3.1, device="svg",width=8, height=3,unit="cm") 

           
###  Calculate the genome-wide correlation of frequencies of CG and CH DMLs
library("tidyverse")

dat <- read.table("CG_CH_density.cor.f.txt",header=T)
colnames(dat) = c("chr", "S","E", "M_CG_hyper", "M_CG_hypo","F_CG_hyper", "F_CG_hypo", "Ctl_CG_hyper", "Ctl_CG_hypo", "Ko_CG_hyper", "Ko_CG_hypo","M_CH_hyper", "M_CH_hypo","F_CH_hyper", "F_CH_hypo", "Ctl_CH_hyper", "Ctl_CH_hypo", "Ko_CH_hyper", "Ko_CH_hypo")
selected_columns <- dat[, -c(1:3)]
correlation_matrix_pearson <- cor(selected_columns, method = "pearson")
                
  cor_pvalue_matrix <- function(df) {
  # Get the number of columns
  n <- ncol(df)

  # Initialize an empty matrix to store p-values
  pvalue_matrix <- matrix(NA, n, n)

  # Iterate through each pair of variables and compute p-values
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      pvalue_matrix[i, j] <- cor.test(df[[i]], df[[j]], method = "pearson")$p.value
      pvalue_matrix[j, i] <- pvalue_matrix[i, j]
    }
    pvalue_matrix[i, i] <- 1
  }
  pvalue_matrix[n, n] <- 1

  # Set column and row names
  colnames(pvalue_matrix) <- colnames(df)
  rownames(pvalue_matrix) <- colnames(df)

  return(pvalue_matrix)
}

cor_matrix = correlation_matrix_pearson
pvalue_matrix <- cor_pvalue_matrix(selected_columns)

library(ComplexHeatmap)
library(circlize)
heatmap <- Heatmap(cor_matrix,
                   name = "correlation",
                   col = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "darkred")),
                   show_column_names = TRUE,
                   show_row_names = TRUE,
                   cluster_columns = T,
                   cluster_rows = T
)

svg("CG.CH.frequency.pearson.heatmap.0722.svg")
draw(heatmap, heatmap_legend_side = "bot")
dev.off()


           
 

