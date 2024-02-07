#
#   Autor Yuxiang Li
#
### upset plot of D1-MSNs from four comparisons (Tet1-wise: male KO vs male Con, female KO vs female Con, sex-wise: male Con vs female Con, male KO vs female KO)
       
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


        
### Plot the mCG  level of representative DMRs near CpG island

set.seed(123)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

 
datam_DMR = read.table(paste0("./DMR.selected.me.txt"), header=F)

colnames(datam_DMR) = c("chr", "Male_WT_Rep2","Male_WT_Rep3","Male_WT_Rep4","Male_WT_Rep5","Male_WT_Rep6","Male_WT_Rep1","Male_KO_Rep3","Male_KO_Rep4","Male_KO_Rep5","Male_KO_Rep1","Male_KO_Rep2","Female_WT_Rep1","Female_WT_Rep2","Female_WT_Rep3","Female_WT_Rep4","Female_WT_Rep5","Female_KO_Rep3","Female_KO_Rep4","Female_KO_Rep5","Female_KO_Rep6","Female_KO_Rep1","Female_KO_Rep2","M","F","WT","KO","GENE")

dat_all = datam_DMR

sample_group= sapply(strsplit(colnames(dat_all[,2:23]),split="_Re"),"[[",1)

df2 = data.frame(Group=sample_group)
Group=factor(df2$Group,levels=c("Male_WT","Female_WT","Male_KO", "Female_KO"))
col= list(Group=c("Female_KO" ="darkmagenta", "Male_KO"="darkorange2","Female_WT" = "purple", "Male_WT"="orange" ))

ha = HeatmapAnnotation(df = df2, col = col)


library(RColorBrewer)


dat3[is.na(dat3)] = 0

svg('CGI_gene_2.svg',width = 8,height = 8)
h1 <- Heatmap(dat3,name="normalized MeCG",
             row_km = 1,
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
            row_order = c("chr12_3235113_3235165", "chr12_3235610_3235659", "chr8_40287042_40287098", "chr10_63583673_63583799", "chr10_63583981_63584030"),
            column_order=c("Male_WT_Rep6","Male_WT_Rep2","Male_WT_Rep3","Male_WT_Rep4","Male_WT_Rep1","Male_WT_Rep5","Female_WT_Rep1","Female_WT_Rep4","Female_WT_Rep5","Female_WT_Rep3","Female_WT_Rep2", "Male_KO_Rep5","Male_KO_Rep4","Male_KO_Rep3","Male_KO_Rep1","Male_KO_Rep2", "Female_KO_Rep1","Female_KO_Rep2","Female_KO_Rep4","Female_KO_Rep6","Female_KO_Rep5","Female_KO_Rep3")      
)

h1
dev.off()

ha_genes=rowAnnotation(foo = anno_mark(at = c(1:35), labels = c(dat_all$GENE)))


dat_test = dat_all[,24:27]
df4=data.frame(dat_test)
colrow2 = list("M"=c("M_hyper"= "black", "M_hypo"="grey"), "F"= c("F_hyper"="black","F_hypo"="gray"), "WT"=c("WT_hyper"="black","WT_hypo"="grey"), "KO"=c("KO_hyper"= "black", "KO_hypo"="grey")   )
ha3 = HeatmapAnnotation(df = df4, which = 'row', col = colrow2,na_col = "white")
svg('CGI.genes.2.svg',width = 8,height = 3)
  h3=draw(h1+ha3+ha_genes )
dev.off()


library("Publish")
dat_AM_DMR = data.frame("Region" =row.names(dat3))
 dat_AM_DMR = cbind(dat_AM_DMR, dat3[,1:22])
anov_res = data.frame()

dat_DMR = dat_AM_DMR
for (i in 1:5) {   
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
dat_p[4,1] <- 0.0001
dat_p_m =matrix(as.numeric(unlist(dat_p)), nrow=5, byrow = F )
rownames(dat_p_m) = rownames(dat3)

dat_p_m[dat_p_m > 0.05] <- NA
dat_p_m_1 = dat_p_m/0.05
write.table(anov_res,"CGI.anova_res.0722.txt",quote=F,sep="\t")
g2 = Heatmap(dat_p_m_1 , name = "p-value",show_row_names = F,col=colorRamp2(c(1e-4, 1), c("olivedrab4", "yellow")), na_col="white",cluster_rows = F, cluster_columns = F, use_raster = F,  show_column_names = T,row_dend_reorder = FALSE,column_dend_reorder = F,cluster_column_slices = FALSE)
svg('Ch_DMR_p.pvalue.4.selceted.0722.2.svg',width = 8,height =4)
draw(h1+g2+ha3+ha_genes)
dev.off()

### plot interaction network between AP1, Mef2 and their target genes

library("igraph")
library("readr")
library("tidyr")
library("RColorBrewer")

df = read.table("AP1_Mef2.txt",header=T)
df = df[which(df$DMR=="M_hyper"  | df$DMR=="M_hypo" | df$DMR=="F_hyper" | df$DMR=="F_hypo" ), ]
Y = unique(df[,c(1,3)])
Y1 <- aggregate(DMR ~ gene, data = Y, FUN = function(x) paste(x, collapse = "_"))
write.table(Y1, "MF.meta.txt", quote=F,row.names = F, sep="\t")
Y1 = read.table("MF.meta.fix.txt", header=T)
X = as.data.frame(table(df[,1:2]))
X1<-subset(X,Freq>0) 

Stucont2<-graph_from_data_frame(X1, directed = FALSE, vertices = Y1)
E(Stucont2)$weight<-E(Stucont2)$Freq # Assigning edge attribute to each edge

set.seed(1)
library(RColorBrewer) # This is the color library
pal = c("pink", "steelblue2","red2", "gray", "blue3", "black", "black", "green2", "darkolivegreen")
deg = degree (Stucont2)
set.seed(9)
theta <- -pi / 3.5  # 45 degrees   //    (Stucont2) %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
# Rotating layout

svg("MF_network-1.3.fruchterman7.4.svg", width =30, height = 30)
plot(Stucont2,edge.color = "#00000080", edge.lty=2, vertex.label.cex=1, vertex.label.dist=0.4, 
     vertex.color=pal[as.numeric(as.factor(vertex_attr(Stucont2, "DMR")))],
     vertex.size = sqrt(deg)+3, edge.width=0.5,
       layout = layout.fruchterman.reingold(Stucont2) %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2))
dev.off()
##  only AP1 for M vs F 

df = read.table("AP1_Mef2.txt",header=T)
df= df[which(df$TF=="AP_1"  ), ]
df = df[which(df$DMR=="Con_hyper"  | df$DMR=="Con_hypo" | df$DMR=="KO_hyper" | df$DMR=="KO_hypo" ), ]
Y = unique(df[,c(1,3)])
Y1 <- aggregate(DMR ~ gene, data = Y, FUN = function(x) paste(x, collapse = "_"))
write.table(Y1, "Con_KO.meta.txt", quote=F,row.names = F, sep="\t")
Y1 = read.table("Con_KO.meta.fix.txt", header=T)
X = as.data.frame(table(df[,1:2]))
X1<-subset(X,Freq>0) # Delete all the edges having weight equal to 0
head(X1)
#Create an igraph object from the dataframes
library(igraph)

Stucont2<-graph_from_data_frame(X1, directed = FALSE, vertices = Y1)
E(Stucont2)$weight<-E(Stucont2)$Freq # Assigning edge attribute to each edge

set.seed(1)
library(RColorBrewer) # This is the color library

pal2 = c( "red2", "gray", "orange", "black", "blue3", "white", "pink", "gray", "steelblue2", "darkolivegreen"  )
deg = degree (Stucont2)

set.seed(9)
theta <- -pi / 3.5  # 45 degrees   //    (Stucont2) %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
# Rotating layout
svg("Con_KO_network-1.3.fruchterman7.4.svg", width =30, height = 30)
plot(Stucont2,edge.color = "#00000080", edge.lty=2, vertex.label.cex=1, vertex.label.dist=0.4, 
     vertex.color=pal2[as.numeric(as.factor(vertex_attr(Stucont2, "DMR")))],
     vertex.size = sqrt(deg)+3, edge.width=0.5,
       layout = layout.fruchterman.reingold(Stucont2) %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2))
dev.off()





                
                


