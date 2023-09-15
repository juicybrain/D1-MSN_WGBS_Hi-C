#!/usr/bin/Rscript

library(hicrep)
library(strawr)
library(reshape2)
library(ggplot2)
library("ComplexHeatmap")
library("circlize")

        args<-commandArgs(T)
        chr=args[1]
        res=as.numeric(args[2])
        med = args[3]

            for (spl in c("HiC_10_S10.allValidPairs","HiC_11_S11.allValidPairs","HiC_12_S12.allValidPairs","HiC_1_S1.allValidPairs","HiC_2_S2.allValidPairs","HiC_3_S3.allValidPairs","HiC_4_S4.allValidPairs","HiC_5_S5.allValidPairs","HiC_7_S7.allValidPairs","HiC_8_S8.allValidPairs","HiC_9_S9.allValidPairs")){ assign(paste0(spl,"_",chr), hic2mat(paste0(spl,".hic"), chromosome1 = chr, chromosome2 = chr, resol = res, method = med))}

            a=c()
            for (spl in c("HiC_10_S10.allValidPairs","HiC_11_S11.allValidPairs","HiC_12_S12.allValidPairs","HiC_1_S1.allValidPairs","HiC_2_S2.allValidPairs","HiC_3_S3.allValidPairs","HiC_4_S4.allValidPairs","HiC_5_S5.allValidPairs","HiC_7_S7.allValidPairs","HiC_8_S8.allValidPairs","HiC_9_S9.allValidPairs")){for (spl2 in c("HiC_10_S10.allValidPairs","HiC_11_S11.allValidPairs","HiC_12_S12.allValidPairs","HiC_1_S1.allValidPairs","HiC_2_S2.allValidPairs","HiC_3_S3.allValidPairs","HiC_4_S4.allValidPairs","HiC_5_S5.allValidPairs","HiC_7_S7.allValidPairs","HiC_8_S8.allValidPairs","HiC_9_S9.allValidPairs") ){ a=c(a,cor(as.vector(get(paste0(spl,"_",chr))), as.vector(get(paste0(spl2,"_",chr))), use = "complete.obs"))}}

        cor_MTX = matrix(a,nrow=11,byrow=T)
        melted_cormat <- melt(cor_MTX)

        pdf(paste0(chr,"_",med ,"_",res,"_ggplot.pdf"),width=10, height =9)
        ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + geom_tile( ) + scale_fill_gradientn(colours =c("black","blue","cyan","green","yellow","pink"), limits=c(0.9,1)) + geom_text(size=2,aes(label =round(value,3))) + theme_minimal()
        dev.off()


                    colnames(cor_MTX) = c("M_KO_Rep3",   "F_WT_Rep3",   "F_KO_Rep3",  "M_WT_Rep1",    "M_KO_Rep1",   "F_WT_Rep1",   "F_KO_Rep1",   "M_KO_Rep2",   "F_WT_Rep2",   "F_KO_Rep2",   "M_WT_Rep2" )
                    rownames(cor_MTX) = c("M_KO_Rep3",   "F_WT_Rep3",   "F_KO_Rep3", "M_WT_Rep1",    "M_KO_Rep1",   "F_WT_Rep1",   "F_KO_Rep1",   "M_KO_Rep2",   "F_WT_Rep2",   "F_KO_Rep2",   "M_WT_Rep2" )


                    write.table(cor_MTX, paste0(chr,".cor.matrix.txt"), quote=F, sep="\t", row.names=F, col.names=T)
                    min_cor = min(cor_MTX)
                    step = (1-min(cor_MTX))/5

                    svg(paste0(chr, "_", med, "_", res, "_test_complexHeatmap.svg"),width=7, height =6)
                    Heatmap(cor_MTX,   name=" r ", col=colorRamp2(c(min_cor, min_cor+step, min_cor+2*step, min_cor+3*step, min_cor+4*step, 1), c("blue","cyan","green","orange","pink","red")), show_row_names =T,show_column_names=T)
                    dev.off()
