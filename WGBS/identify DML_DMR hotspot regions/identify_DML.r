### identify DMLs (CG coverage >10, delta beta value >0.05, p-value <0.01) for DHR identification
#
#  Author Yuxiang Li
#

dat_all = read.table("../../D1_MF_all_sample.fixed.s.autosome.DSS",header=F)
keep <- which(rowSums(dat_all[,c(4,6,8,10,12,14)]>=10)>=3 & rowSums(dat_all[,c(16,18,20,22,24)]>=10)>=3 )
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
dmlTest = DMLtest(BSobj, group2=c("M1","M2","M3","M4","M5","M11"), group1=c("M7","M8","M9","M10","M12"),smoothing=F)

        dmls0.05 = callDML(dmlTest, delta=0.05, p.threshold=1e-2)
        write.table(dmls0.05,"dml_D1_Smth_M_KOvsCon_5percent1e-2.noS.txt",sep="\t",row.name=FALSE,quote=F)

