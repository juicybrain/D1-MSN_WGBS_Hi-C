#
#   Author Yuxiang Li
#
### differential analysis of diagonal Hi-C contact (1.5Mb) comparing male KO to male Con, similar code was used to compare female KO to female Con, and male to female groups

library(edgeR)
set.seed(123)

wd="./1.5M500ks/"
files=dir(wd)[grep("1.5Mw500s.bed", dir(wd))]
dat_1.5M = read.table(paste0(wd, files[1]), header=F, skip=1)
colnames(dat_1.5M)=c( "chr","start","end", unlist(strsplit(files[1], split=".100k"))[1])

for (i in 2:length(files)) {
#for (i in 2:5) {
  dat_tmp= read.table(paste0(wd, files[i]),  header=F, skip=1)
  colnames(dat_tmp) = c("chr", "start", "end", unlist(strsplit(files[i], split=".100k"))[1])
  dat_1.5M = merge(dat_1.5M, dat_tmp, by=c("chr","start", "end"),all=F)  
}
#countData_auto= dat_1.5M[which(dat_1.5M$chr!="chrX"),]
#dat_M_kovswt=countData_auto[,c(4,5,9,10,11)]
#dat_M_kovswt= countData_auto[,c(5,7,8,10,11,13)]

dat_M_kovswt = dat_1.5M[,c(5, 8, 11, 4, 14)]

colnames(dat_M_kovswt) = c("M_ko3" ,"M_ko1","M_ko2", "M_WT1","M_WT2" )
dat_M_kovswt = apply(dat_M_kovswt,2,as.numeric)
rownames(dat_M_kovswt)=paste(dat_1.5M$chr, dat_1.5M$start, dat_1.5M$end,sep="_")
rownames(dat_M_kovswt) = rownames(dat_M_kovswt)
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}

dat_M_kovswt_fixed = delete.na(dat_M_kovswt)
colData <-  data.frame("Treatment"=c(rep("M_KO",3),rep("M_WT",2)),"Sample"= c("M_ko3" ,"M_ko1","M_ko2", "M_WT1","M_WT2" ))
y <- DGEList(dat_M_kovswt_fixed, samples=colData$Sample, group=colData$Treatment)
y <- calcNormFactors(y)

#which(is.na(dat_M_kovswt), arr.ind=TRUE)

#Estimate Error Model
design <-model.matrix(~Treatment, data=colData)
y <- estimateDisp(y, design)
#compute p-values, then output
et <- exactTest(y, pair = c(2,1))
res <- topTags(et,Inf)
tidyResult <- data.frame(Gene=rownames(res$table), res$table)
tidyResult$threshold <- ifelse(tidyResult$FDR<0.05 & abs(tidyResult$logFC)>=0.1375,ifelse(tidyResult$logFC>0.1375,"up","down"),"NOT")

dat_count = as.data.frame(dat_M_kovswt_fixed)
dat_count$Gene = row.names(dat_M_kovswt_fixed)
tidyResult_m = merge(tidyResult, dat_count,by="Gene",all=T)

write.table(tidyResult_m,file="M_kovswt_MvsF_Hi_C_res1.5Mw500s.txt",sep="\t",row.names=FALSE, quote=F)
write("Dispersion = ",stderr())
write(y$common.dispersion,stderr())

###Volcano###
library("tidyverse")

windowsFonts(Times=windowsFont("TT Times New Roman"))
p1 <- ggplot(data = tidyResult, aes(x = logFC , y = -log10(FDR), color=threshold)) +
  geom_point(alpha=1, size = 0.3) +
  xlim(c(-1, 1)) +
  ylim(c(0, 80)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-0.1375,0.1375),lty=4,col="grey",lwd=0.2)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.2)+
  theme(
    legend.position="right",
    panel.grid=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_color_manual(name = "", values = c("blue", "red", "grey"), limits = c("down", "up", "NOT"))

#+labs(x="log2 (Fold Change)",y="-log10 (P-value)",title="Drd1-WT vs M_kovswt")
p1.1 = p1 +theme_classic(base_size = 5)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 0,vjust=0.5),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0, 'cm'),
                                                     legend.key.height = unit(0.5, 'cm'),
                                                     legend.key.width = unit(0, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.grid.minor = element_blank())

ggsave("D1M_kovswt.diff.svg",plot=p1.1, device="svg",width=4.5, height=3.5,unit="cm")  


### Plot TAD corner score comparisons

library("tidyverse")
library(ggpubr)

dat_MN = read.table("MN_TAD_in_DIR.bed.txt")
colnames(dat_MN) = c("chr","MN" )
dat_MP = read.table("MP_TAD_in_DIR.bed.txt")
colnames(dat_MP) = c("chr","MP" )

dat_FN = read.table("FN_TAD_in_DIR.bed.txt")
colnames(dat_FN) = c("chr","FN" )

dat_FP = read.table("FP_TAD_in_DIR.bed.txt")


colnames(dat_FP) = c("chr","FP" )

dat = merge(dat_MN, dat_MP, by="chr",all=F)

dat = merge(dat, dat_FN, by="chr",all=F)

dat = merge(dat, dat_FP, by="chr",all=F)


t.test(dat$MN, dat$MP,paired = T)
t.test(dat$FN, dat$FP,paired = T)
 
g1 = ggpaired(dat, cond1="MN", cond2="MP", fill="condition") + scale_fill_manual(values = c("MN"="orange",  "MP" = "darkorange2")) + theme_classic(base_size = 6)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                                axis.text.x=element_blank(),
                                                 legend.position ="",
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm'),
                                                 panel.grid.minor = element_blank())

g2 = ggpaired(dat, cond1="FN", cond2="FP", fill="condition") + scale_fill_manual(values = c("FN"="purple",  "FP" = "darkmagenta"))  + theme_classic(base_size = 6)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                                axis.text.x=element_blank(),
                                                 legend.position ="",
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm'),
                                                 panel.grid.minor = element_blank())
ggsave("D1.MN_MP.svg",plot=g1, device="svg",width=6, height=6,unit="cm")
ggsave("D1.FN_FP.svg",plot=g2, device="svg",width=6, height=6,unit="cm")
