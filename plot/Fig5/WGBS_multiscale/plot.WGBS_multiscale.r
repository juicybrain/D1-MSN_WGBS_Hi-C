#
#  Author Yuxiang Li
#
#!/usr/bin/Rscript

dat_1 = read.table("male.mm10.1.bed.counts.bed.meCG.bed",header=F)
colnames(dat_1) = c("chr","start","end","M1","F1")

dat_all = dat_1
for (i in 2:51){
    dat_tmp = read.table(paste0("male.mm10.",i,".bed.counts.bed.meCG.bed"),header=F)
    colnames(dat_tmp) = c("chr","start","end", paste0("M",i),  paste0("F",i))
    dat_all=merge(dat_all,dat_tmp,by=c("chr","start","end"),all=T)
}

dim(dat_all)

        dat_matrix = dat_all[,-c(1,2,3)]
        dat_M=as.matrix(dat_matrix[,paste0("M",1:51)])
        dat_F=as.matrix(dat_matrix[,paste0("F",1:51)])

        library("ComplexHeatmap")
        library("circlize")
        g3 = Heatmap(dat_F,name = "multiscale_F",show_row_names = F,  col=colorRamp2(c(-0.2, -0.1, -0.01, 0, 0.01,  0.1, 0.2),   c("darkblue", "blue", "white","white", "white", "red", "red4")), na_col="black",cluster_rows = F, cluster_columns = F, use_raster = T,  show_column_names = F,row_dend_reorder = FALSE,column_dend_reorder = F,cluster_column_slices = FALSE)
        g4 = Heatmap(dat_M,name = "multiscale_M",show_row_names = F,  col=colorRamp2(c(-0.2, -0.1, -0.01, 0, 0.01,  0.1, 0.2),  c("darkblue", "blue", "white","white", "white", "red", "red4")), na_col="black",cluster_rows = F, cluster_columns = F, use_raster = T,  show_column_names = F,row_dend_reorder = FALSE,column_dend_reorder = F,cluster_column_slices = FALSE)

library("zoo")
        svg('chr8.MF.KOminusWT.multiscale.svg',width = 3, height =8)
        draw(g3+g4)
        dev.off()
