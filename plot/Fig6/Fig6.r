#
#  Author Yuxiang Li (mathewlyx@gmail.com)
#
### GO and KEGG pathway analysis of sex-covergent and non-convergent DMRs

library("ChIPseeker")
library("tidyverse")
library("DOSE")
library("clusterProfiler") 
library("org.Mm.eg.db") 
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 

dat_sex_convergent <- readPeakFile("sex.consist.CGDMR.annoedbyloop.bed")
dat_sex_divergent <-  readPeakFile("sex.specific.CGDMR.annoedbyloop.bed")

peaks <- list(consis=dat_sex_convergent, diff= dat_sex_divergent )

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 100), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Mm.eg.db")
dat_gene_list = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
dat_geneName_list = lapply(peakAnnoList, function(i) as.data.frame(i)$SYMBOL)

for (i in 1:length(dat_geneName_list)) {
    cat("Table", i, "\n")
    write.table(dat_geneName_list[[i]], paste("CG_DMR_", names(dat_geneName_list)[i], ".txt"), quote=F, row.names = F, col.names = F  ) 
    cat("\n")
}

all_name= dat_geneName_list %>%flatten()%>%as_vector()%>%unique()

write.table(all_name,"Sex_convergent_and_divergent_DMR.gene.name.txt",quote=F,row.names = F)



compGO_CC <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "CC",pvalueCutoff=0.05,qvalueCutoff=0.05)
p7.1 <- dotplot(compGO_CC,showCategory = 10,includeAll=T) + scale_color_gradient(low="red", high="gold")
ggsave("sex_diff_DMR.loop3.GO.CC.p0.01.0.5.svg",plot=p7.1, device="svg",width=10, height=12)

GO_cluster_summary <- as.data.frame(compGO_CC)
ComB <- function(x){paste(bitr(unlist(strsplit(x,"/")), fromType = "ENTREZID",toType = c("ENSEMBL","SYMBOL"), OrgDb = org.Mm.eg.db)$SYMBOL,collapse = "/")}
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "Sex_convergent_and_divergent_DMR.GO.CC.txt",sep="\t",quote=F,row.names = F)

compGO_MF <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "MF",pvalueCutoff=0.05,qvalueCutoff=0.05)
p8.1 <- dotplot(compGO_MF,showCategory = 10,includeAll=T) + scale_color_gradient(low="red", high="gold")
ggsave("Sex_convergent_and_divergent_DMR.GO.MF.svg",plot=p8.1, device="svg",width=10, height=12)

GO_cluster_summary <- as.data.frame(compGO_MF)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "Sex_convergent_and_divergent_DMR.GO.MF.txt",sep="\t",quote=F,row.names = F)

compGO_BP <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "BP",pvalueCutoff=0.05,qvalueCutoff=0.05)
p9.1 <- dotplot(compGO_BP,showCategory = 10,includeAll=T) + scale_color_gradient(low="red", high="gold")
ggsave("Sex_convergent_and_divergent_DMR.GO.BP.svg",plot=p9.1, device="svg",width=6, height=8)

GO_cluster_summary <- as.data.frame(compGO_BP)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "Sex_convergent_and_divergent_DMR.GO.BP.txt",sep="\t",quote=F,row.names = F)

compReactome <- compareCluster(dat_gene_list,fun="enrichPathway",organism = "mouse", pvalueCutoff=0.05)

p10.1 <- dotplot(compReactome,showCategory = 10,includeAll=T) + scale_color_gradient(low="red", high="gold")
ggsave("Sex_convergent_and_divergent_DMR.Reactome.svg",plot=p10.1, device="svg",width=10, height=12)

compKEGG <- compareCluster(dat_gene_list,fun="enrichKEGG",organism = "mouse", pvalueCutoff=1)
p10.1 <- dotplot(compKEGG,showCategory = 15,includeAll=T) + scale_color_gradient(low="red", high="gold")
ggsave("Sex_convergent_and_divergent_DMR.KEGG.p1.svg",plot=p10.1, device="svg",width=5.2, height=6.5)

GO_cluster_summary <- as.data.frame(compKEGG)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "Sex_convergent_and_divergent_DMR.KEGG.txt",sep="\t",quote=F,row.names = F)


### Plot the mCG level of CG CG DMR associated with DMG in the "Amphetamine addiction pathway"

set.seed(123)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

datam_DMR = read.table(paste0( "DHR.synaptic_genes.CG.txt"), header=F)


colnames(datam_DMR) = c("chr","start","end", "gene", "Male_WT_Rep2","Male_WT_Rep3","Male_WT_Rep4","Male_WT_Rep5","Male_WT_Rep6","Male_WT_Rep1","Male_KO_Rep3","Male_KO_Rep4","Male_KO_Rep5","Male_KO_Rep1","Male_KO_Rep2","Female_WT_Rep1","Female_WT_Rep2","Female_WT_Rep3","Female_WT_Rep4","Female_WT_Rep5","Female_KO_Rep3","Female_KO_Rep4","Female_KO_Rep5","Female_KO_Rep6","Female_KO_Rep1","Female_KO_Rep2", "MN", "MP", "FN", "FP")

dat_all = datam_DMR

sample_group= sapply(strsplit(colnames(dat_all[,5:26]),split="_Re"),"[[",1)

df2 = data.frame(Group=sample_group)
Group=factor(df2$Group,levels=c("Male_WT","Male_KO","Female_WT","Female_KO"))
col= list(Group=c("Female_KO" ="darkmagenta", "Male_KO"="darkorange2","Female_WT" = "purple", "Male_WT"="orange" ))

ha = HeatmapAnnotation(df = df2, col = col)


library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# 74 distinctive colors in R:
distinctive_colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

dat2 <- scale(t(dat_all[, 27:30]),scale=F, center=T)
dat3 <- t(dat2)
row.names(dat3) = paste(dat_all$chr,dat_all$start, dat_all$gene,sep="_")


write.table(dat3,"D1_DMR.new.07.19.scale.txt",quote=F,sep="\t")
dat3[is.na(dat3)] = 0

svg('DHR_gene_1.1.svg',width = 8,height = 8)
h1.1 <- Heatmap(dat3,name="normalized MeCG",
           row_km =1,
            row_gap = unit(1, "mm"),
            column_gap = unit(1, "mm"),
            row_km_repeats = 10,
            col=colorRamp2(c(-0.1,-0.05,0,0.05,0.1), c("darkblue","blue", "white", "red","darkred")),
            clustering_distance_rows =  "euclidean",
           clustering_method_rows = "centroid",
            show_row_names = FALSE,
            row_title = "%s",
            border = TRUE,
            show_row_dend = F,
            use_raster = F,
            cluster_column_slices = FALSE,
          column_order=c("MN", "MP", "FN", "FP")      
)

h1.1
dev.off()

ha_genes=rowAnnotation(foo = anno_mark(at = c(1:22), labels = c(dat_all$gene)))

svg('DHR.genes.07.25.svg',width = 8,height = 10)
  h3=draw(h1.1+ha_genes  )
dev.off()

dat2 <- scale(t(dat_all[, 5:26]),scale=F, center=T)
dat3 <- t(dat2)
row.names(dat3) = paste(dat_all$chr,dat_all$start, dat_all$gene,sep="_")
write.table(dat3,"D1_DMR.new.07.19.scale.txt",quote=F,sep="\t")
dat3[is.na(dat3)] = 0

dat_AM_DMR = data.frame("Region" = row.names(dat3))
 dat_AM_DMR = cbind(dat_AM_DMR, dat3[,1:22])
anov_res = data.frame()

library("Publish")

dat_DMR = dat_AM_DMR
for (i in 1:22) {   
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
dat_p[3,1] = 0.0001
dat_p[17:20, 2]=0.0001

dat_p_m =matrix(as.numeric(unlist(dat_p)), nrow=22, byrow = F )
rownames(dat_p_m) = rownames(dat3)

dat_p_m[dat_p_m > 0.05] <- NA
dat_p_m_1 = dat_p_m/0.05
g1 = Heatmap(dat_p_m_1 , name = "p-value",show_row_names = F,col=colorRamp2(c(1e-4, 1), c("olivedrab4", "yellow")), na_col="NA",cluster_rows = F, cluster_columns = F, use_raster = F,  show_column_names = T,row_dend_reorder = FALSE,column_dend_reorder = F,cluster_column_slices = FALSE)
svg('AM_DMR_p.pvalue.2.07.25.2.svg',width = 12,height =12)
draw(ha_genes+h1.1+g1)
dev.off()
write.table(anov_res,"DHR_genes.anova_res.txt",quote=F,sep="\t")


### GO and KEGG pathway analysis of DHZ gene


library("ChIPseeker")
library("tidyverse")
library("DOSE")
library("clusterProfiler") 
library("org.Mm.eg.db") 
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 


dat_in <- readPeakFile("all.DMR.in.DHR.bed")
dat_out <- readPeakFile("all.DMR.out.DHR.bed")
dat_DHR_genes <- readPeakFile("DHR.loop.gene.bed")

peaks <- list(dmr_in=dat_in, dmr_out=dat_out, DHR=dat_DHR_genes )
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 100), verbose=FALSE,addFlankGeneInfo=TRUE,                       #flankDistance=5000,
                       annoDb="org.Mm.eg.db")


dat_geneName_anno = lapply(peakAnnoList, function(i) as.data.frame(i))
library(writexl)
write_xlsx(dat_geneName_anno, "anno.xlsx")


dat_gene = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
dat_gene_list = list("DMR_in"=unique(dat_gene$dmr_in), "DMR_out"=  unique(setdiff(dat_gene$dmr_out,  dat_gene$DHR)),  "DHR" = unique(dat_gene$DHR))
for (i in 1:length(dat_gene_list)){ print(length(dat_gene_list[[i]]))}
dat_geneName_list = lapply(peakAnnoList, function(i) as.data.frame(i)$SYMBOL)
for (i in 1:length(dat_geneName_list)) {
    cat("Table", i, "\n")
    write.table(dat_geneName_list[[i]], paste("CG_DMR_", names(dat_geneName_list)[i], ".txt"), quote=F, row.names = F, col.names = F  ) 
    cat("\n")}
all_name= dat_geneName_list %>%flatten()%>%as_vector()%>%unique()
write.table(all_name,"all.DMR.DHR.gene.name.txt",quote=F,row.names = F)

compGO_CC <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "CC",pvalueCutoff=0.05)
p7.1 <- dotplot(compGO_CC,showCategory = 5,includeAll=T)
p7.2 = p7.1 + scale_color_gradient(low="red", high="gold")
ggsave("DMR.DHR.cc.svg",plot=p7.2, device="svg",width=5, height=6)
GO_cluster_summary <- as.data.frame(compGO_CC)
ComB <- function(x){paste(bitr(unlist(strsplit(x,"/")), fromType = "ENTREZID",toType = c("ENSEMBL","SYMBOL"), OrgDb = org.Mm.eg.db)$SYMBOL,collapse = "/")}
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "DMR.DHR.GO.CC.txt",sep="\t",quote=F,row.names = F)
compGO_MF <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "MF",pvalueCutoff=0.05,qvalueCutoff=0.05)
p8.1 <- dotplot(compGO_MF,showCategory = 10,includeAll=T)+scale_color_gradient(low="red", high="gold")
ggsave("DMR.DHR.GO.MF.svg",plot=p8.1, device="svg",width=10, height=12)
GO_cluster_summary <- as.data.frame(compGO_MF)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "DMR.DHR.GO.MF.txt",sep="\t",quote=F,row.names = F)
compGO_BP <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "BP",pvalueCutoff=0.05,qvalueCutoff=0.05)
compGO_BP_simp <- simplify(compGO_BP, cutoff = 0.5, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)
p9.1 <- dotplot(compGO_BP,showCategory = 10,includeAll=T)+scale_color_gradient(low="red", high="gold")
ggsave("DMR.DHR.GO.BP.top10.longer.svg",plot=p9.1, device="svg",width=6, height=26)
GO_cluster_summary <- as.data.frame(compGO_BP)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "DMR.DHR.GO.BP.txt",sep="\t",quote=F,row.names = F)
p9.2 <- dotplot(compGO_BP_simp,showCategory = 15,includeAll=T)+scale_color_gradient(low="red", high="gold")
ggsave("DMR.DHR.GO.BP.simp.svg",plot=p9.2, device="svg",width=6, height=25)
GO_cluster_summary <- as.data.frame(compGO_BP_simp)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "DMR.DHR.GO.BP.simp.txt",sep="\t",quote=F,row.names = F)
compReactome <- compareCluster(dat_gene_list,fun="enrichPathway",organism = "mouse", pvalueCutoff=0.1)
p10.1 <- dotplot(compReactome,showCategory = 10,includeAll=T)+scale_color_gradient(low="red", high="gold")
ggsave("DMR.DHR.Reactome.svg",plot=p10.1, device="svg",width=10, height=12)
GO_cluster_summary <- as.data.frame(compReactome)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "DMR.DHR.Reactome.txt",sep="\t",quote=F,row.names = F)
compKEGG <- compareCluster(dat_gene_list,fun="enrichKEGG",organism = "mouse", pvalueCutoff=0.05)
p10.1 <- dotplot(compKEGG,showCategory = 15,includeAll=T)+scale_color_gradient(low="red", high="gold")
ggsave("DMR.DHR.KEGG.p1.svg",plot=p10.1, device="svg",width=6, height=8)
GO_cluster_summary <- as.data.frame(compKEGG)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "DMR.DHR.KEGG.txt",sep="\t",quote=F,row.names = F)


### Plot the delta mCG in HCN1 with standard error(se)

library("tidyverse")
dat = read.table("promoter/DHR_synaptic_genes.batch2_promoter_1.bedgraph", header=F)
dat$V5 = 1
for (i in 2:10){
  dat_tmp = read.table(paste0("promoter/DHR_synaptic_genes.batch2_promoter_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i
  dat=rbind(dat,dat_tmp)
}

for (i in 1:50){
  dat_tmp = read.table(paste0("gene_body/DHR_synaptic_genes.batch2_gene_body_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i+10
  dat=rbind(dat,dat_tmp)
}

for (i in 1:10){
  dat_tmp = read.table(paste0("downstream/DHR_synaptic_genes.batch2_downstream_",i,".bedgraph"),header =F)
  dat_tmp$V5 = i+60
  dat=rbind(dat,dat_tmp)
}

dat_s = na.omit(dat[order(dat$V4, dat$V5),])



write.table(dat_s, "DHR_synaptic_genes.batch2.scaled.txt", quote=F, row.names = F, sep="\t")

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
 
 ### plot KO-Con 
 
    dat_KOvsCtl_M = dat_MP - dat_MN_mean
    dat_KOvsCtl_F = dat_FP - dat_FN_mean

  for (spl in list("dat_KOvsCtl_M", "dat_KOvsCtl_F")){
      assign(paste0(spl, "_mean"), rowMeans(get(spl)) )
      assign(paste0(spl, "_sd"),  apply(get(spl), 1, sd) )
  }

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
 # geom_point(size=0.1) +
 # theme_minimal() +
  labs(x = "X-axis label", y = "Y-axis label", title = paste0(gene," Line Plot with Standard Error Shadows")) +
  scale_color_manual(values = c("KOvsCtl_M" = "orange", "KOvsCtl_F" = "purple")) +
  scale_fill_manual(values = c("KOvsCtl_M" = "orange", "KOvsCtl_F" = "purple"))+
  scale_x_continuous(breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) + theme_classic(base_size=5, base_line_size =0.2)
  
  ggsave( paste0(gene,".meCG.gene.in.KOvsCon.MF.svg"),plot=p, width = 2,height = 1)
  }
   

 # plot MvsF  
 
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
 # geom_point(size=0.1) +
 # theme_minimal() +
  labs(x = "X-axis label", y = "Y-axis label", title = paste0(gene," Line Plot with Standard Error Shadows")) +
  scale_color_manual(values = c("MvsF_Ctl" = "orange4", "MvsF_Ko" = "green")) +
  scale_fill_manual(values = c("MvsF_Ctl" = "orange4", "MvsF_Ko" = "green"))+
  scale_x_continuous(breaks = c(0,10, 60,70), labels = c("-10","TSS", "TES","10")) +scale_y_continuous(limits = c(-0.2, 0.2)) + theme_classic(base_size=5, base_line_size =0.2)
  
  ggsave( paste0(gene,".meCG.gene.in.MvsF.Ctl.svg"),plot=p, width = 2,height = 1) 
 }











