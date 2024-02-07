#
#  Author Yuxiang Li (mathewlyx@gmail.com)
#
### Plot the mCG level of CG CG DMR associated with DMG in the "Amphetamine addiction pathway"

set.seed(123)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

wd="Amphetamine_addiction_genes/"
datam_DMR = read.table(paste0(wd, "Am_addiction.DMR.txt"), header=F)


colnames(datam_DMR) = c("chr", "Male_WT_Rep2","Male_WT_Rep3","Male_WT_Rep4","Male_WT_Rep5","Male_WT_Rep6","Male_WT_Rep1","Male_KO_Rep3","Male_KO_Rep4","Male_KO_Rep5","Male_KO_Rep1","Male_KO_Rep2","Female_WT_Rep1","Female_WT_Rep2","Female_WT_Rep3","Female_WT_Rep4","Female_WT_Rep5","Female_KO_Rep3","Female_KO_Rep4","Female_KO_Rep5","Female_KO_Rep6","Female_KO_Rep1","Female_KO_Rep2","M","F","WT","KO","GENE","Sex","Tet1","Sex:Tet1","cl")
dat_all = datam_DMR
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
row.names(dat3) = dat_all$chr


write.table(dat3,"D1_DMR.new.07.19.scale.txt",quote=F,sep="\t")
dat3[is.na(dat3)] = 0


dat_all$cl=factor(dat_all$cl,levels = c("A","B"))

svg('Amphetamine_gene_1.svg',width = 8,height = 8)
h1 <- Heatmap(dat3,name="normalized MeCG",
#           row_km = 2,
            row_split = dat_all$cl,
            row_gap = unit(1, "mm"),
            column_gap = unit(1, "mm"),
 #           na_col = "black",
            row_km_repeats = 10,
            top_annotation = ha,
            col=colorRamp2(c(-0.5,-0.25,0, 0.25,0.5), c("darkblue", "blue", "white", "red", "darkred")),
 #           right_annotation = ha2,
#            rowAnnotation(foo = anno_histogram(m)),
 #           clustering_distance_columns = "euclidean",
#            clustering_method_columns = "centroid",
            clustering_distance_rows =  "euclidean",
           clustering_method_rows = "centroid",
            show_row_names = FALSE,
            row_title = "%s",
            border = TRUE,
            show_row_dend = F,
            use_raster = TRUE,
#           column_split=Group,
#           column_km=4,
#           width = unit(12, "cm"), height = unit(20, "cm"),
            column_split=Group,
            cluster_column_slices = FALSE,
#            show_column_names = T,
              ## order follows pretested K-means cluster results (n=2)
            row_order = c("chr11_57176086_57176100", "chr2_174300445_174300496", "chr2_174300692_174300800", "chr6_119007467_119007573", "chr6_119013481_119013528", "chr3_136565031_136565049", "chr3_136835425_136835483", "chr14_70273363_70274681","chr12_100303350_100303476", "chr17_87448937_87448998", "chr2_92027799_92027966", "chr6_136053082_136053179", "chr6_136055946_136056048", "chr6_136119790_136120045", "chr6_136158990_136159043", "chr6_136205579_136205766", "chr6_136212625_136212970", "chr6_136227720_136227763", "chr13_53800455_53800506","chr3_126610345_126610404","chr11_6006301_6006543"),
             
             )

h1
dev.off()

ha_genes=rowAnnotation(foo = anno_mark(at = c(1:21), labels = c(dat_all$GENE)))

dat_test = dat_all[,24:27]
df4=data.frame(dat_test)
colrow2 = list("M"=c("M_hyper"= "black", "M_hypo"="grey"), "F"= c("F_hyper"="black","F_hypo"="grey"), "WT"=c("WT_hyper"="black","WT_hypo"="grey"), "KO"=c("KO_hyper"= "black", "KO_hypo"="grey")   )

ha3 = HeatmapAnnotation(df = df4, which = 'row', col = colrow2,na_col = "white")
svg('Amphetamine_addiction.genes.0725.svg',width = 8,height = 10)
  h3=draw(h1+ha3+ha_genes )
dev.off()

dat_AM_DMR = data.frame("Region" = row.names(dat3))
 dat_AM_DMR = cbind(dat_AM_DMR, dat3[,1:22])
anov_res = data.frame()

library("Publish")

dat_DMR = dat_AM_DMR
for (i in 1:21) {   
  assign(dat_DMR[i,1], as.data.frame(t(dat_DMR[i,-1]))) 
  dat_tmp = get(dat_DMR[i,1])
  colnames(dat_tmp)=c(dat_DMR[i,1])
  dat_tmp$Sex=c(rep(1,11),rep(0,11))
  dat_tmp$TET1=c(rep(1,6),rep(0,5),rep(1,5),rep(0,6) )
  df_aov_smy = c(dat_DMR[i,1],publish(summary(aov(data=dat_tmp, as.formula(paste0(dat_DMR[i,1],"~", "Sex+TET1+Sex*TET1")))),print=F))
  anov_res=rbind(anov_res,df_aov_smy)
}
colnames(anov_res)=c("loci",paste0("Df_",c("Sex","TET1","Sex:TET1")),  paste0("F-value_",c("Sex","TET1","Sex:TET1")), paste0("Pr",c("Sex","TET1","Sex:TET1")) )
rm(list=c(dat_DMR$Region))



dat_p = anov_res[8:10]
dat_p[12,3] <- 0.0001
dat_p_m =matrix(as.numeric(unlist(dat_p)), nrow=21, byrow = F )
rownames(dat_p_m) = rownames(dat3)

dat_p_m[dat_p_m > 0.05] <- NA
dat_p_m_1 = dat_p_m/0.05

g1 = Heatmap(dat_p_m_1 , name = "p-value",show_row_names = F,col=colorRamp2(c(1e-4, 1), c("red", "gold")), na_col="NA",cluster_rows = F, cluster_columns = F, use_raster = F,  show_column_names = T,row_dend_reorder = FALSE,column_dend_reorder = F,cluster_column_slices = FALSE)
svg('AM_DMR_p.pvalue.2.0725.svg',width = 12,height =6)
draw(h1+g1+ha3+ha_genes)
dev.off()

write.table(anov_res,"am_genes.anova_res.txt",quote=F,sep="\t")



### plot DHR gene

set.seed(123)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

wd="D:/D1_WGBS/final version/nature/version4/fig3/chosen_DHR_Genes/final_picked_whole_genes/"
datam_DMR = read.table(paste0(wd, "DHR.selected.genes.CG.txt"), header=F)


colnames(datam_DMR) = c("chr","start","end", "gene", "Male_WT_Rep2","Male_WT_Rep3","Male_WT_Rep4","Male_WT_Rep5","Male_WT_Rep6","Male_WT_Rep1","Male_KO_Rep3","Male_KO_Rep4","Male_KO_Rep5","Male_KO_Rep1","Male_KO_Rep2","Female_WT_Rep1","Female_WT_Rep2","Female_WT_Rep3","Female_WT_Rep4","Female_WT_Rep5","Female_KO_Rep3","Female_KO_Rep4","Female_KO_Rep5","Female_KO_Rep6","Female_KO_Rep1","Female_KO_Rep2")

dat_all = datam_DMR

sample_group= sapply(strsplit(colnames(dat_all[,5:26]),split="_Re"),"[[",1)

df2 = data.frame(Group=sample_group)
Group=factor(df2$Group,levels=c("Male_WT","Female_WT","Male_KO", "Female_KO"))
col= list(Group=c("Female_KO" ="darkmagenta", "Male_KO"="darkorange2","Female_WT" = "purple", "Male_WT"="orange" ))

ha = HeatmapAnnotation(df = df2, col = col)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# 74 distinctive colors in R:
distinctive_colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

dat2 <- scale(t(dat_all[, 5:26]),scale=F, center=T)
dat3 <- t(dat2)
row.names(dat3) = paste(dat_all$chr,dat_all$start, dat_all$gene,sep="_")


write.table(dat3,"D1_DMR.new.07.19.scale.txt",quote=F,sep="\t")
dat3[is.na(dat3)] = 0

svg('DHR_gene_1.cluster4.svg',width = 8,height = 8)
h1 <- Heatmap(dat3,name="normalized MeCG",
           row_km =4,
            row_gap = unit(1, "mm"),
            column_gap = unit(1, "mm"),
            row_km_repeats = 10,
            top_annotation = ha,
            col=colorRamp2(c(-0.1,-0.05,0,0.05,0.1), c("darkblue","blue", "white", "red","darkred")),
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

ha_genes=rowAnnotation(foo = anno_mark(at = c(1:42), labels = c(dat_all$gene)))


svg('DHR.genes.07.25.svg',width = 8,height = 10)
  h3=draw(h1+ha_genes )
dev.off()


dat_AM_DMR = data.frame("Region" = row.names(dat3))
 dat_AM_DMR = cbind(dat_AM_DMR, dat3[,1:22])
anov_res = data.frame()

library("Publish")

dat_DMR = dat_AM_DMR
for (i in 1:23) {   
  assign(dat_DMR[i,1], as.data.frame(t(dat_DMR[i,-1]))) 
  dat_tmp = get(dat_DMR[i,1])
  colnames(dat_tmp)=c(dat_DMR[i,1])
  dat_tmp$Sex=c(rep(1,11),rep(0,11))
  dat_tmp$TET1=c(rep(1,6),rep(0,5),rep(1,5),rep(0,6) )
  df_aov_smy = c(dat_DMR[i,1],publish(summary(aov(data=dat_tmp, as.formula(paste0(dat_DMR[i,1],"~", "Sex+TET1+Sex*TET1")))),print=F))
  anov_res=rbind(anov_res,df_aov_smy)
}
colnames(anov_res)=c("loci",paste0("Df_",c("Sex","TET1","Sex:TET1")),  paste0("F-value_",c("Sex","TET1","Sex:TET1")), paste0("Pr",c("Sex","TET1","Sex:TET1")) )
rm(list=c(dat_DMR$Region))

dat_p = anov_res[8:10]
dat_p[4,1] = 0.0001
dat_p[19:22, 2]=0.0001

dat_p_m =matrix(as.numeric(unlist(dat_p)), nrow=23, byrow = F )
rownames(dat_p_m) = rownames(dat3)

dat_p_m[dat_p_m > 0.05] <- NA
dat_p_m_1 = dat_p_m/0.05

g1 = Heatmap(dat_p_m_1 , name = "p-value",show_row_names = F,col=colorRamp2(c(1e-4, 1), c("red", "gold")), na_col="NA",cluster_rows = F, cluster_columns = F, use_raster = F,  show_column_names = T,row_dend_reorder = FALSE,column_dend_reorder = F,cluster_column_slices = FALSE)
svg('DHR.genes_p.pvalue.2.07.25.svg',width = 12,height =12)
draw(ha_genes+h1+g1+ha4)
dev.off()

write.table(anov_res,"DHR_genes.anova_res.txt",quote=F,sep="\t")

rm(list=ls())
### Plot the mCG level around the whole gene loci of DHR representative genes: Hcn1, Zdhhc2, and Grin2b (scale according to the gene length to 10bins upstream TSS, 50 bins for gene body, and 10 bins downstream of TTS)

library("tidyverse")

############################ male hypo ###############################################
dat = read.table("promoter/D1_M_WTvsCre.hypo.target.gene.f_promoter_1.bedgraph", header=F)
dat$V5 = 1
for (i in 2:10){
  dat_tmp = read.table(paste0("promoter/D1_M_WTvsCre.hypo.target.gene.f_promoter_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i
  dat=rbind(dat,dat_tmp)            }

for (i in 1:50){
  dat_tmp = read.table(paste0("gene_body/D1_M_WTvsCre.hypo.target.gene.f_gene_body_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i+10
  dat=rbind(dat,dat_tmp)            }

for (i in 1:10){
  dat_tmp = read.table(paste0("downstream/D1_M_WTvsCre.hypo.target.gene.f_downstream_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i+60
  dat=rbind(dat,dat_tmp)            }

dat_s = dat[order(dat$V4, dat$V5),]
write.table(dat_s, "cDMR.male.hypo.gene.scaled.txt", quote=F, row.names = F, sep="\t")

dat_MN = dat_s[, 7:12]
dat_MP = dat_s[, 13:17]
dat_FN = dat_s[, 18:22]
dat_FP = dat_s[, 23:28]

for (spl in list("dat_MN", "dat_MP", "dat_FN", "dat_FP")){     
      assign(paste0(spl, "_mean"), rowMeans(get(spl)) )
      assign(paste0(spl, "_sd"),  apply(get(spl), 1, sd) )             }

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
 
# plot all genes 
  # Create a plot for the current category
    dat_KOvsCtl_M = dat_MP - dat_MN_mean
    dat_KOvsCtl_F = dat_FP - dat_FN_mean

  for (spl in list("dat_KOvsCtl_M", "dat_KOvsCtl_F")){
      assign(paste0(spl, "_mean"), rowMeans(get(spl)) )
      assign(paste0(spl, "_sd"),  apply(get(spl), 1, sd) ) }

  dat_KOvsCtl_M_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "MvsF_Ctl_mean"= dat_KOvsCtl_M_mean, "MvsN_Ctl_sd" = dat_KOvsCtl_M_sd/sqrt(5))
  dat_KOvsCtl_F_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "MvsF_Ctl_mean"= dat_KOvsCtl_F_mean, "MvsN_Ctl_sd" = dat_KOvsCtl_F_sd/sqrt(6)) 
   
  colnames(dat_KOvsCtl_M_mean_sd) = c("gene", "region", "mean", "sd") 
  dat_KOvsCtl_M_mean_sd$group = "KOvsCtl_M"
  colnames(dat_KOvsCtl_F_mean_sd) = c("gene", "region", "mean", "sd") 
  dat_KOvsCtl_F_mean_sd$group = "KOvsCtl_F"
  dat_plot_Ctl_Ko = rbind(dat_KOvsCtl_M_mean_sd, dat_KOvsCtl_F_mean_sd )
  dat_plot_Ctl_Ko = dat_plot_Ctl_Ko[order(dat_plot_Ctl_Ko$gene, dat_plot_Ctl_Ko$region ),]
  df=dat_plot_Ctl_Ko
  write.table(df, "KOvsCon.MF.cDMRgene.txt", quote=F,sep="\t", row.names = F)

 uniq_genes <- unique(dat_plot_M_s$gene)
 for (gene in uniq_genes) {
  # Subset the data frame for the current category
  subset_df <- df[df$gene == gene, ]
  # Create a plot for the current category
  p <- ggplot(subset_df, aes(x = region, y = mean, color = group, group = group)) +
  geom_line(size=0.25) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd, fill = group), alpha = 0.3, size=0) +
  labs(x = "X-axis label", y = "Y-axis label", title = paste0(gene," Line Plot with Standard Error Shadows")) +
  scale_color_manual(values = c("KOvsCtl_M" = "orange", "KOvsCtl_F" = "purple")) +
  scale_fill_manual(values = c("KOvsCtl_M" = "orange", "KOvsCtl_F" = "purple"))+
  scale_x_continuous(breaks = c(0,10, 60,70), labels = c("-10", "TSS", "TES", "10")) + theme_classic(base_size=5, base_line_size =0.2)  
  ggsave( paste0(gene,".meCG.gene.in.KOvsCon.MF.svg"),plot=p, width = 2,height = 1)
  } 


### plot M-F
    dat_MvsF_Ctl = dat_s[, 7:12] - dat_FN_mean
    dat_MvsF_Ko =  dat_s[, 13:17] - dat_FP_mean
  for (spl in list("dat_MvsF_Ctl", "dat_MvsF_Ko")){
      assign(paste0(spl, "_mean"), rowMeans(get(spl)) )
      assign(paste0(spl, "_sd"),  apply(get(spl), 1, sd) )
  }

   dat_MvsF_ctl_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "MvsF_Ctl_mean"= dat_MvsF_Ctl_mean, "MvsN_Ctl_sd" = dat_MvsF_Ctl_sd/sqrt(6))
   dat_MvsF_Ko_mean_sd = data.frame("gene"= dat_s$V4, "region"=dat_s$V5, "MvsF_Ko_mean"= dat_MvsF_Ko_mean, "MvsN_Ko_sd" = dat_MvsF_Ko_sd/sqrt(5))   
 colnames(dat_MvsF_ctl_mean_sd) = c("gene", "region", "mean", "sd") 
 dat_MvsF_ctl_mean_sd$group = "MvsF_Ctl"
 
 colnames(dat_MvsF_Ko_mean_sd) = c("gene", "region", "mean", "sd") 
 dat_MvsF_Ko_mean_sd$group = "MvsF_Ko"
 dat_plot_Ctl_Ko = rbind(dat_MvsF_ctl_mean_sd, dat_MvsF_Ko_mean_sd )
 dat_plot_Ctl_Ko = dat_plot_Ctl_Ko[order(dat_plot_Ctl_Ko$gene, dat_plot_Ctl_Ko$region ),]
 df=dat_plot_Ctl_Ko
 write.table(df, "MvsF.cDMRgene.txt", quote=F,sep="\t", row.names = F)

 uniq_genes <- unique(dat_plot_M_s$gene)
 for (gene in uniq_genes) {
  # Subset the data frame for the current category
  subset_df <- df[df$gene == gene, ]  
  # Create a plot for the current category
  p <- ggplot(subset_df, aes(x = region, y = mean, color = group, group = group)) +
  geom_line(size=0.25) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd, fill = group), alpha = 0.3, size=0) +
  labs(x = "X-axis label", y = "Y-axis label", title = paste0(gene," Line Plot with Standard Error Shadows")) +
  scale_color_manual(values = c("MvsF_Ctl" = "blue", "MvsF_Ko" = "red")) +
  scale_fill_manual(values = c("MvsF_Ctl" = "blue", "MvsF_Ko" = "red"))+
  scale_x_continuous(breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) + theme_classic(base_size=5, base_line_size =0.2)
  ggsave( paste0(gene,".meCG.gene.in.MvsF.Ctl.svg"),plot=p, width = 2,height = 1)
 }
 
 




