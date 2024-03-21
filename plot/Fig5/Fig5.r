#
#  Author Yuxiang Li
#

### plot the SCC score of pair-wised comparison of Hi-C samples
library("tidyverse")

df_res=data.frame("pair"=NA,"SCC"=NA,"chr"=NA)
for (chr in c(1:19, "X")){
  dat_tmp = read.table(paste0("./HiC-cor/chr",chr,".cor.matrix.hicRep.txt"), header=T)
  row.names(dat_tmp) = colnames(dat_tmp)   
  d_tmp = c()
  m_tmp = as.matrix(dat_tmp)
  
for (i in 1:dim(m_tmp)[1]){
  for(j in 1:dim(m_tmp)[1]){
    if (i>j){
     d_tmp = rbind(d_tmp, c(paste0(strsplit(rownames(m_tmp)[i],split="_R")[[1]][1],"_" , strsplit(colnames(m_tmp)[j], split = "_R")[[1]][1] ),  m_tmp[i,j], paste0("chr",chr) ) )  
    }
  }
}
  
  df_tmp=as.data.frame(d_tmp)
  colnames(df_tmp)=c("pair", "SCC", "chr")

  assign(paste0("d",chr), df_tmp)
  df_res=rbind(df_res, df_tmp)
  
}

df_backup=df_res

df_res=df_backup[-1,]

df_res$sex= paste0(sapply(strsplit(df_res$pair,split="_"),"[[",1) , sapply(strsplit(df_res$pair,split="_"),"[[",3) )

df_res$tet1= paste0(sapply(strsplit(df_res$pair,split="_"),"[[",2) ,  sapply(strsplit(df_res$pair,split="_"),"[[",4) )


df_res$sex=str_replace_all(df_res$sex,"FM","MF")
df_res$tet1=str_replace_all(df_res$tet1,"KOWT","WTKO")

df_res$chr=factor(df_res$chr,levels = c("chr1", "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX"))
p5 = ggplot(df_res, aes(x=chr, y=as.numeric(SCC))) + geom_violin(width=1,alpha=0.5, size=0.2)+
 # geom_jitter(size=0.7, alpha=1, position=position_jitter(0.3),aes(color=tet1, shape=sex))+
 geom_point( size=2, position = position_jitter(seed = 1, width = 0.3) ,alpha=1, aes(color=tet1, shape=sex) ) +
  scale_color_manual(values = c("KOKO"="blue","WTWT"="blue","WTKO"="red"))+xlab("")+
  ylab("Hi-C data Pair-wise similarity")+ ylim(0.92,1) + scale_shape_manual(values=c("FF"=4,"MM"=4,"MF"=3 ))

p5.1 = p5 + theme_classic(base_size = 6)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                                               axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                               axis.line.x.bottom=element_line(size=0.2),
                                               # axis.text.x=element_blank(),
                                                 legend.position ="right",
                                              legend.title=element_text(size=5), 
                                              panel.grid.major=element_line(colour=NA),
                                               legend.text = element_text(size=5),
                                               legend.key.size = unit(0.2, 'cm'),
                                                legend.key.height = unit(0.2, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm'),
                                                 panel.background = element_rect(fill = NA,size=0.2,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank()) 

ggsave(filename = "./Hic.cor.svg",plot=p5.1,height = 10,width = 18,unit="cm")

### plot the SCC score of pair-wised comparison of Hi-C samples at Chr8

library(ComplexHeatmap)
library("circlize")
  dat = read.table("./chr8.cor.matrix.hicRep.txt", header=T)
  rownames(dat) = colnames(dat)
	min_cor = min(dat)
        step = (1-min(dat))/5
 
        svg(paste0("chr8",  "_cor_complexHeatmap.svg"),width=7, height =6)
        Heatmap(dat,   name=" r ", col=colorRamp2(c(min_cor, min_cor+step, min_cor+2*step, min_cor+3*step, min_cor+4*step, 1), c("blue","cyan","green","orange","pink","red")), show_row_names =T,show_column_names=T)
	dev.off()


### HiC differential analysis

################ male KO vs Control ###################
library(edgeR)
set.seed(123)
#Read Data in

wd="./1.5M500ks/"
files=dir(wd)[grep("1.5Mw500s.bed", dir(wd))]
dat_1.5M = read.table(paste0(wd, files[1]), header=F, skip=1)
colnames(dat_1.5M)=c( "chr","start","end", unlist(strsplit(files[1], split=".100k"))[1])

for (i in 2:length(files)) {
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
  ylim(c(0, 60)) +
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
p1.1 = p1 +theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0, 'cm'),
                                                     legend.key.height = unit(0.5, 'cm'),
                                                     legend.key.width = unit(0, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank())

ggsave("D1M_KOvsCon.diff.svg",plot=p1.1, device="svg",width=5, height=4,unit="cm") 

######################## female KO vs Con ################################

library(edgeR)
set.seed(123)
#Read Data in

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

dat_F_kovswt = dat_1.5M[, c(7,10,13,6,9,12)]

colnames(dat_F_kovswt) = c("F_KO1" ,"F_KO2","F_KO3", "F_WT1","F_WT2","F_WT3" )
dat_F_kovswt = apply(dat_F_kovswt,2,as.numeric)
rownames(dat_F_kovswt)=paste(dat_1.5M$chr, dat_1.5M$start, dat_1.5M$end,sep="_")
#rownames(dat_F_kovswt) = rownames(dat_F_kovswt)
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}


dat_F_kovswt_fixed = delete.na(dat_F_kovswt)
colData <-  data.frame("Treatment"=c(rep("F_KO",3),rep("F_WT",3)),"Sample"= c("F_KO1" ,"F_KO2","F_KO3", "F_WT1","F_WT2","F_WT3" ))

y <- DGEList(dat_F_kovswt_fixed, samples=colData$Sample, group=colData$Treatment)
y <- calcNormFactors(y)


#Estimate Error Model
design <-model.matrix(~Treatment, data=colData)
y <- estimateDisp(y, design)


#compute p-values, then output
et <- exactTest(y, pair = c(2,1))
res <- topTags(et,Inf)
tidyResult <- data.frame(Gene=rownames(res$table), res$table)

tidyResult$threshold <- ifelse(tidyResult$FDR<0.05 & abs(tidyResult$logFC)>=0.1375,ifelse(tidyResult$logFC>0.1375,"up","down"),"NOT")
dat_count = as.data.frame(dat_F_kovswt_fixed )
dat_count$Gene = row.names(dat_F_kovswt_fixed)

tidyResult_m = merge(tidyResult, dat_count,by="Gene",all=T)


write.table(tidyResult_m,file="F_kovswt_MvsF_Hi_C_res1.5Mw500s.txt",sep="\t",row.names=FALSE, quote=F)
write("Dispersion = ",stderr())
write(y$common.dispersion,stderr())

###Volcano###
library("tidyverse")

windowsFonts(Times=windowsFont("TT Times New Roman"))
p1 <- ggplot(data = tidyResult, aes(x = logFC , y = -log10(FDR), color=threshold)) +
  geom_point(alpha=0.8, size = 0.3) +
  xlim(c(-1, 1)) +
  #ylim(c(0, 60)) +
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

p1.1 = p1 +theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0, 'cm'),
                                                     legend.key.height = unit(0.5, 'cm'),
                                                     legend.key.width = unit(0, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank())

ggsave("D1F_KOvsCon.diff.svg",plot=p1.1, device="svg",width=5, height=4,unit="cm")  

########################## Control Male vs Female ##############################


library(edgeR)
set.seed(123)

wd="/1.5M500ks/"
files=dir(wd)[grep("1.5Mw500s.bed", dir(wd))]
dat_1.5M = read.table(paste0(wd, files[1]), header=F, skip=1)
colnames(dat_1.5M)=c( "chr","start","end", unlist(strsplit(files[1], split=".100k"))[1])

for (i in 2:length(files)) {
#for (i in 2:5) {
  dat_tmp= read.table(paste0(wd, files[i]),  header=F, skip=1)
  colnames(dat_tmp) = c("chr", "start", "end", unlist(strsplit(files[i], split=".100k"))[1])
  dat_1.5M = merge(dat_1.5M, dat_tmp, by=c("chr","start", "end"),all=T)  
}

countData_auto= dat_1.5M[which(dat_1.5M$chr!="chrX"),]

#dat_WT=countData_auto[,c(4,5,9,10,11)]
dat_WT= countData_auto[,-c(1,2,3,5,7,8,10,11,13)]
dat_WT = dat_WT[,c(1,5,2,3,4)]

colnames(dat_WT) = c("WT_M1","WT_M2" ,"WT_F1","WT_F3","WT_F2")
dat_WT = apply(dat_WT,2,as.numeric)
rownames(dat_WT)=paste(countData_auto$chr, countData_auto$start,countData_auto$end,sep="_")


rownames(dat_WT) = rownames(dat_WT)
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}

dat_WT_fixed = delete.na(dat_WT)
colData <-  data.frame("Treatment"=c(rep("WT_M",2),rep("WT_F",3)),"Sample"= c("WT_M1","WT_M2" ,"WT_F1","WT_F3","WT_F2"))

 
y <- DGEList(dat_WT_fixed, samples=colData$Sample, group=colData$Treatment)
y <- calcNormFactors(y)

#which(is.na(dat_WT), arr.ind=TRUE)

#Estimate Error Model
design <-model.matrix(~Treatment, data=colData)
y <- estimateDisp(y, design)


#compute p-values, then output
et <- exactTest(y, pair = c(1,2))
res <- topTags(et,Inf)
tidyResult <- data.frame(Gene=rownames(res$table), res$table)

tidyResult$threshold <- ifelse(tidyResult$FDR<0.05 & abs(tidyResult$logFC)>=0.1375,ifelse(tidyResult$logFC>0.1375,"up","down"),"NOT")

dat_count = as.data.frame( dat_WT_fixed )
dat_count$Gene = row.names(dat_WT_fixed)
tidyResult_m = merge(tidyResult, dat_count,by="Gene",all=T)

write.table(tidyResult_m,file="WT_MvsF_Hi_C_res1.5Mw500s.allchr.txt",sep="\t",row.names=FALSE,quote=F)

write("Dispersion = ",stderr())
write(y$common.dispersion,stderr())



###Volcano###
library("tidyverse")

windowsFonts(Times=windowsFont("TT Times New Roman"))
p1 <- ggplot(data = tidyResult, aes(x = logFC , y = -log10(FDR), color=threshold)) +
  geom_point(alpha=1, size = 0.3) +
  xlim(c(-1, 1)) +
  ylim(c(0,5)) +
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

#+labs(x="log2 (Fold Change)",y="-log10 (P-value)",title="Drd1-WT vs WT")
p1.1 = p1 +theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0, 'cm'),
                                                     legend.key.height = unit(0.5, 'cm'),
                                                     legend.key.width = unit(0, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank())

ggsave("D1Con_MvsF.diff.2.allchr.svg",plot=p1.1, device="svg",width=5, height=4,unit="cm")  

#################### TET1 KO male vs TET1 KO female #########################

library(edgeR)
set.seed(123)

wd="D:/D1_WGBS/final version/fig3/hic_diff_0708/1.5M500ks/"
files=dir(wd)[grep("1.5Mw500s.bed", dir(wd))]
dat_1.5M = read.table(paste0(wd, files[1]), header=F, skip=1)
colnames(dat_1.5M)=c( "chr","start","end", unlist(strsplit(files[1], split=".100k"))[1])

for (i in 2:length(files)) {
  dat_tmp= read.table(paste0(wd, files[i]),  header=F, skip=1)
  colnames(dat_tmp) = c("chr", "start", "end", unlist(strsplit(files[i], split=".100k"))[1])
  dat_1.5M = merge(dat_1.5M, dat_tmp, by=c("chr","start", "end"),all=F)  
}

countData_auto= dat_1.5M[which(dat_1.5M$chr!="chrX"),]

#dat_ko=countData_auto[,c(4,5,9,10,11)]
dat_ko= countData_auto[,c(5,8,11,7,10,13)]
colnames(dat_ko) = c("ko_M1","ko_M2","ko_M3" ,"ko_F1","ko_F2","ko_F3")
dat_ko = apply(dat_ko,2,as.numeric)
rownames(dat_ko)=paste(countData_auto$chr, countData_auto$start,countData_auto$end,sep="_")
rownames(dat_ko) = rownames(dat_ko)
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}

dat_ko_fixed = delete.na(dat_ko)
colData <-  data.frame("Treatment"=c(rep("ko_M",3),rep("ko_F",3)),"Sample"= c("ko_M1","ko_M2","ko_M3" ,"ko_F1","ko_F2","ko_F3"))

 
 
y <- DGEList(dat_ko_fixed, samples=colData$Sample, group=colData$Treatment)
y <- calcNormFactors(y)

#which(is.na(dat_ko), arr.ind=TRUE)

#Estimate Error Model
design <-model.matrix(~Treatment, data=colData)
y <- estimateDisp(y, design)


#compute p-values, then output
et <- exactTest(y,pair=c(1,2))
res <- topTags(et,Inf)
tidyResult <- data.frame(Gene=rownames(res$table), res$table)

tidyResult$threshold <- ifelse(tidyResult$FDR<0.05 & abs(tidyResult$logFC)>=0.1375,ifelse(tidyResult$logFC>0.1375,"up","down"),"NOT")

dat_count = as.data.frame(dat_ko_fixed )
dat_count$Gene = row.names(dat_ko_fixed)

tidyResult_m = merge(tidyResult, dat_count,by="Gene",all=T)

write.table(tidyResult_m, file="ko_MvsF_Hi_C_res1.5Mw500s.txt",sep="\t",row.names=FALSE)

write("Dispersion = ",stderr())
write(y$common.dispersion,stderr())

###Volcano###
library("tidyverse")

windowsFonts(Times=windowsFont("TT Times New Roman"))
p1 <- ggplot(data = tidyResult, aes(x = logFC , y = -log10(FDR), color=threshold)) +
  geom_point(alpha=0.8, size = 0.3) +
  xlim(c(-1, 1)) +
 ylim(c(0, 5)) +
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

#+labs(x="log2 (Fold Change)",y="-log10 (P-value)",title="Drd1-WT vs KO")
p1.1 = p1 +theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 5,angle = 45,vjust=0.5),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0, 'cm'),
                                                     legend.key.height = unit(0.5, 'cm'),
                                                     legend.key.width = unit(0, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank())

ggsave("D1KO_MvsF.diff.svg",plot=p1.1, device="svg",width=5, height=4,unit="cm")  









