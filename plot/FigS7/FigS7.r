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
 
 
 ################################################
 ## PCA
 ################################################
 

