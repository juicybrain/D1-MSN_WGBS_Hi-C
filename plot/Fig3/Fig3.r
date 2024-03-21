####
#
#  Author Yuxiang Li
#



### plot the global DNA methylation level
library("tidyverse")
dat_plot = read.table("global_meCG_level.txt",header=T)

df_mRNA_mean = dat_plot%>%group_by(group)%>%summarise(mRNA_mean=mean(mCG_level) )

#df_mRNA_sd = dat_plot%>%group_by(group)%>%summarise(mRNA_sd=sd(mCG_level) )
df_mRNA_se = dat_plot %>% group_by(group) %>% summarise(mRNA_se = sd(mCG_level) / sqrt(n()))

dat_plot = merge(df_mRNA_mean,df_mRNA_se,by=c("group"))
dat_plot$group=factor(dat_plot$group, levels=c("M_Ctl","M_KO","F_Ctl","F_KO"))
dat_plot_se = merge(df_mRNA_mean,df_mRNA_se,by=c("group"))
dat_plot_se$group=factor(dat_plot_se$group, levels=c("M_Ctl","M_KO","F_Ctl","F_KO"))
 
p1<- ggplot(dat_plot, aes(x=group, y=mRNA_mean, fill=group)) + 
  geom_col(stat="identity",
           position=position_dodge(),size=0.4,width=0.4, col="black") +
  geom_errorbar(aes(ymin=mRNA_mean-mRNA_se, ymax= mRNA_mean+mRNA_se),size=0.2, width=0.2,col="black",
                 position=position_dodge(.9))+ scale_fill_manual(values = c("M_Ctl"="orange", "M_KO"="orange" , "F_Ctl" = "purple", "F_KO" ="purple" ))  + xlab("")+ylab("mRNA level (fold change)") + coord_cartesian(ylim = c(0.79,0.81),expand=F)+expand_limits(x =c(0.5,4.5))


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

ggsave("D1.mCG.svg",plot=p1.2, device="svg",width=3, height=3,unit="cm")



### Upset plot show the overlap of D1-MSN DMRs
        library(regioneR)
        library(ComplexHeatmap)


        M_hyper = toGRanges("M_hyper.bed")
        M_hypo = toGRanges("M_hypo.bed")
        F_hyper = toGRanges("F_hyper.bed")
        F_hypo = toGRanges("F_hypo.bed")
        KO_hypo = toGRanges("KO_hypo.bed")
        KO_hyper = toGRanges("KO_hyper.bed")
        WT_hypo = toGRanges("WT_hypo.bed")
        WT_hyper = toGRanges("WT_hyper.bed")
        
        
        input_data=list(M_hyper=M_hyper, M_hypo=M_hypo, F_hyper=F_hyper, F_hypo=F_hypo, WT_hyper=WT_hyper, WT_hypo=WT_hypo, KO_hyper=KO_hyper, KO_hypo=KO_hypo)
        m2_intersect = make_comb_mat(input_data,mode="intersect")
         
        svg("DMR.intersect.test5.svg",width=6,height=4)
        UpSet(m2_intersect,set_order = c("M_hyper", "M_hypo", "F_hyper", "F_hypo", "WT_hyper", "WT_hypo", "KO_hyper", "KO_hypo"),comb_order = order(comb_size(m2_intersect)))
        dev.off()



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



### plot GO annotation of CG DMR

library("ChIPseeker")
library("tidyverse")
library("DOSE")
library("clusterProfiler") 
library("org.Mm.eg.db") 
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 

dat_M_hyper <- readPeakFile("dmrs_D1_M_WTvsCre.delta0.p1e-3.txt.bedgraph.hyper.bedgraph")
dat_M_hypo <- readPeakFile("dmrs_D1_M_WTvsCre.delta0.p1e-3.txt.bedgraph.hypo.bedgraph")
dat_F_hyper <- readPeakFile("dmrs_D1_F_WTvsCre_p1e-3.delta0.mergedSS2DS.txt.bedgraph.hyper.bedgraph")
dat_F_hypo <- readPeakFile("dmrs_D1_F_WTvsCre_p1e-3.delta0.mergedSS2DS.txt.bedgraph.hypo.bedgraph")
dat_WT_hyper <- readPeakFile("dmrs_D1_WT_Male_vs_Female.delta0p1e-3.mergeSS2DS.txt.bedgraph.hyper.bedgraph")
dat_WT_hypo <- readPeakFile("dmrs_D1_WT_Male_vs_Female.delta0p1e-3.mergeSS2DS.txt.bedgraph.hypo.bedgraph")
dat_KO_hyper <- readPeakFile("dmrs_D1_KO_Male_vs_Female.delta0p1e-3.mergeSS2DS.txt.bedgraph.hyper.bedgraph")
dat_KO_hypo <- readPeakFile("dmrs_D1_KO_Male_vs_Female.delta0p1e-3.mergeSS2DS.txt.bedgraph.hypo.bedgraph")

peaks <- list(M_hyper=dat_M_hyper, M_hypo=dat_M_hypo, F_hyper=dat_F_hyper, F_hypo=dat_F_hypo, WT_hyper=dat_WT_hyper, WT_hypo=dat_WT_hypo,  KO_hyper=dat_KO_hyper, KO_hypo = dat_KO_hypo )
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 100), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Mm.eg.db")
dat_gene_list = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

compGO_BP <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "BP",pvalueCutoff=0.05,qvalueCutoff=0.05)
p6.1 <- dotplot(compGO_BP,showCategory = 10, includeAll=T)+ scale_color_gradient(high="gold", low="red")
ggsave("MF.DMR.GO_BP.top10.svg",plot=p6.1, device="svg",width=8, height=10)



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




### Cluster D1MSN-DMRs

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

