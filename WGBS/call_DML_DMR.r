### This is an R script for identifying CG DMLs and DMRs for WGBS data comparing male Tet1 KO to male control
#
#  Author  Yuxiang Li (mathewlyx@gmail.com)
#
                
dat_all = read.table("../D1_MF_all_sample.fixed.s.autosome.DSS",header=F)
keep <- which(rowSums(dat_all[,c(4,6,8,10,12,14)]>=5)>=3 & rowSums(dat_all[,c(16,18,20,22,24)]>=5)>=3 )
dat_keep <- dat_all[keep,]
        dat_M7 = dat_keep[,c(1,2,16,15)]
        colnames(dat_M7) <- c("chr",     "pos",     "N",       "X")
        dat_M8 = dat_keep[,c(1,2,18,17)]
        colnames(dat_M8) <- c("chr",     "pos",     "N",       "X")
        dat_M9 = dat_keep[,c(1,2,20,19)]
        colnames(dat_M9) <- c("chr",     "pos",     "N",       "X")
        dat_M10 = dat_keep[,c(1,2,22,21)]
        colnames(dat_M10) <- c("chr",     "pos",     "N",       "X")
        dat_M12 = dat_keep[,c(1,2,24,23)]
        colnames(dat_M12) <- c("chr",     "pos",     "N",       "X")

        dat_M1 = dat_keep[,c(1,2,4,3)]
        colnames(dat_M1) <- c("chr",     "pos",     "N",       "X")
        dat_M2 = dat_keep[,c(1,2,6,5)]
        colnames(dat_M2) <- c("chr",     "pos",     "N",       "X")
        dat_M3= dat_keep[,c(1,2,8,7)]
        colnames(dat_M3) <- c("chr",     "pos",     "N",       "X")
        dat_M4= dat_keep[,c(1,2,10,9)]
        colnames(dat_M4) <- c("chr",     "pos",     "N",       "X")
        dat_M5 = dat_keep[,c(1,2,12,11)]
        colnames(dat_M5) <- c("chr",     "pos",     "N",       "X")
        dat_M11 = dat_keep[,c(1,2,14,13)]
        colnames(dat_M11) <- c("chr",     "pos",     "N",       "X")

BSobj = makeBSseqData( list(dat_M1,dat_M2,dat_M3,dat_M4,dat_M5,dat_M7,dat_M8,dat_M9,dat_M10,dat_M11,dat_M12),c("M1","M2","M3","M4","M5","M7","M8","M9","M10","M11","M12"))
dmlTest = DMLtest(BSobj, group2=c("M1","M2","M3","M4","M5","M11"), group1=c("M7","M8","M9","M10","M12"),smoothing=TRUE,smoothing.span=100)
save(dmlTest,file="dmlTest_M_Ko_vs_Ctl_smooth_100_rSave")
        # CG dmls used in Figure 3d for the visualization of hotspot regions
        dmls0.05 = callDML(dmlTest, delta=0.05, p.threshold=1e-4)
        write.table(dmls0.05,"dml_D1_Smth_M_KOvsCon_5percent1e-4.txt",sep="\t",row.name=FALSE,quote=F)

        # identify CG DMRs 
        dmrs0 = callDMR(dmlTest, delta=0, p.threshold=1e-3,minlen=6, minCG=3,dis.merge=100)
        write.table(dmrs0,"dmrs_D1_M_KOvsCon.delta0.p1e-3.txt",sep="\t",quote=F,row.name=FALSE)


