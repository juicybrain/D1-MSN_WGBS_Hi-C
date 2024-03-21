#
#   Author Yuxiang Li
#


### Plot DMR genomic feature 
library("tidyverse")
library("reshape")
dat = read.table("DMR.genomic_feature.input.barchart.txt", header=T)
dat_p = melt(dat,id.vars = "Group")

dat_p$Group=factor(dat_p$Group, levels=c("M","F","Con","KO"))

colors <- rainbow(7)
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2")

g1= ggplot(dat_p, aes(x = Group, y = value, fill = variable)) + 
    geom_bar(position = "fill",stat = "identity") + scale_fill_brewer(palette="Set3")

 

### Plot the length of DHR
library("tidyverse")
M_hyper = read.table("M.hyper.cDMR.bed",header=F)
M_hypo = read.table("M.hypo.cDMR.bed",header=F)

F_hyper = read.table("F.hyper.cDMR.bed",header=F)
F_hypo = read.table("F.hypo.cDMR.bed",header=F)

Ctl_hyper = read.table("Con_hyper.cDMR.bed", header=F)
Ctl_hypo = read.table("Con_hypo.cDMR.bed", header=F)

Ko_hyper = read.table("KO.hyper.cDMR.bed", header=F)
Ko_hypo = read.table("KO.hypo.cDMR.bed", header=F)

groups=c()
size=c()

for (spl in c("M_hyper", "M_hypo", "F_hyper", "F_hypo", "Ctl_hyper", "Ctl_hypo", "Ko_hyper", "Ko_hypo")){
  assign(paste(spl,"size",sep = "_"),  get(spl)$V3- get(spl)$V2 )
  groups =c(groups, rep(spl,length(get(paste(spl,"size",sep = "_")))))
  size=c(size,unlist(get(paste(spl,"size",sep = "_"))))
}


df = data.frame(DHR_size=size, group=groups)
df$group = factor(df$group, levels=c("M_hyper", "M_hypo", "F_hyper", "F_hypo", "Ctl_hyper", "Ctl_hypo", "Ko_hyper", "Ko_hypo"))

library(ggplot2)

g3=ggplot(df, aes(x = group, fill=df$group ,y = log10(DHR_size))) +
  geom_boxplot(width=0.3) +
  stat_summary(fun = function(x) quantile(x, 0.75),
               geom = "point",
               shape = 5,
               size = 0,
               color = "black",
               show.legend = FALSE) +
  ylab("Values") +
  xlab("Group") + scale_fill_manual(values = c("M_hyper"="orange", "M_hypo"="orange", "F_hyper"="purple", "F_hypo"="purple", "Ctl_hyper"="darkgreen", "Ctl_hypo"="darkgreen", "Ko_hyper"="cyan", "Ko_hypo"="cyan"))+ ylim(c(3,7))


g3.1 = g3 + theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     #axis.text.x =element_text(size = 5,angle = 90, vjust=0.5),
                                                      axis.text.x=element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0.2, 'cm'),
                                                     legend.key.height = unit(0.2, 'cm'),
                                                     legend.key.width = unit(0.2, 'cm'),
                                                     legend.box = "top",
                                                    # panel.grid.major=element_line(colour=NA),
                                                    # panel.background = element_rect(fill = NA,size=0.5,colour = "black"),
                                                    # plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank()) 
                                                 
ggsave("DHRs.length.svg",plot=g3.1, device="svg",width=8, height=3,unit="cm") 


### Plot the proportion of CG DMRs and DMLs inside of DHZs


library("tidyverse")
wd="D:/D1_WGBS/final version/fig2_version_2/final_DML_DMR_cDMR_0611/"
dat = read.table(paste0("DHR_stats_DMR.txt"), header=T,sep="\t")
g1 = ggplot(data = dat, mapping = aes(x = factor(group,levels = c("M KO vs WT hyper","M KO vs WT hypo","F KO vs WT hyper","F KO vs WT hypo","WT M vs F hyper","WT M vs F hypo","KO M vs F hyper","KO M vs F hypo")), y = Quantity, fill = DHR)) + geom_bar(stat = 'identity', position = 'stack',size=0.6,width=0.6) + scale_fill_manual(values = c("in"="skyblue","out"="darkgray"))+xlab("") + scale_y_continuous(expand = c(0, 0),limits=c(0, 800))+ylab("DMR quantity")

g1.2 <- g1 + theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     #axis.text.x =element_text(size = 5,angle = 90, vjust=0.5),
                                                      axis.text.x=element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0.2, 'cm'),
                                                     legend.key.height = unit(0.2, 'cm'),
                                                     legend.key.width = unit(0.2, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.background = element_rect(fill = NA,size=0.5,colour = "black"),
                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank()) 
                                                 
ggsave("DML_in_DHR.g4.DMR.svg",plot=g1.2, device="svg",width=6, height=3,unit="cm") 



dat = read.table(paste0(wd,"DHR_stats_DML.1.txt"), header=T,sep="\t")
g2 = ggplot(data = dat, mapping = aes(x = factor(group,levels = c("M KO vs WT hyper","M KO vs WT hypo","F KO vs WT hyper","F KO vs WT hypo","WT M vs F hyper","WT M vs F hypo","KO M vs F hyper","KO M vs F hypo")), y = Quantity, fill = DHR)) + geom_bar(stat = 'identity', position = 'stack',size=0.6,width=0.6) + scale_fill_manual(values = c("in"="skyblue","out"="darkgray"))+xlab("") + scale_y_continuous(expand = c(0, 0),limits=c(0, 7000)) + ylab("DML quantity")
g2.2 <- g2 + theme_classic(base_size = 7)+ theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.title= element_blank(),
                                                     #axis.text.x =element_text(size = 5,angle = 90, vjust=0.5),
                                                      axis.text.x=element_blank(),
                                                     #legend.position =c(0.3,0.7),
                                                     legend.position ="right",
                                                     legend.text = element_text(size=5),
                                                     legend.key.size = unit(0.2, 'cm'),
                                                     legend.key.height = unit(0.2, 'cm'),
                                                     legend.key.width = unit(0.2, 'cm'),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.background = element_rect(fill = NA,size=0.5,colour = "black"),
                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank())                                                 
ggsave("DML_in_DHR.g4.DML.svg",plot=g2.2, device="svg",width=6, height=3,unit="cm") 


### Plot the anti-correlation of CG DML and CH DML

library("tidyverse")

dat <- read.table("CG_CH_density.cor.f.txt",header=T)
colnames(dat) = c("chr", "S","E", "M_CG_hyper", "M_CG_hypo","F_CG_hyper", "F_CG_hypo", "Ctl_CG_hyper", "Ctl_CG_hypo", "Ko_CG_hyper", "Ko_CG_hypo","M_CH_hyper", "M_CH_hypo","F_CH_hyper", "F_CH_hypo", "Ctl_CH_hyper", "Ctl_CH_hypo", "Ko_CH_hyper", "Ko_CH_hypo")
selected_columns <- dat[, -c(1:3)]
correlation_matrix_pearson <- cor(selected_columns, method = "pearson")

cor_pvalue_matrix <- function(df) {
  # Get the number of columns
  n <- ncol(df)

  # Initialize an empty matrix to store p-values
  pvalue_matrix <- matrix(0, n, n)

  # Iterate through each pair of variables and compute p-values
  for (i in 1:n) {
    for (j in i:n) {
      pvalue_matrix[i, j] <- cor.test(df[[i]], df[[j]], method = "pearson")$p.value
      pvalue_matrix[j, i] <- pvalue_matrix[i, j]
    }
   # pvalue_matrix[i, i] <- 1
  }
  #pvalue_matrix[n, n] <- 1

  # Set column and row names
  colnames(pvalue_matrix) <- colnames(df)
  rownames(pvalue_matrix) <- colnames(df)

  return(pvalue_matrix)
}

cor_matrix = correlation_matrix_pearson

pvalue_matrix <- cor_pvalue_matrix(selected_columns)




color_mapping <- function(cor_matrix, pvalue_matrix, sig_level = 4e-4) {
  n <- ncol(cor_matrix)
  col_matrix <- matrix("NA", n, n)

  for (i in 1:n) {
    for (j in i:n) {
      if (pvalue_matrix[i, j] < sig_level) {
        if (cor_matrix[i, j] > 0) {
          col_matrix[i, j] <- "blue"
        } else {
          col_matrix[i, j] <- "red"
        }
      }
    }
  }

  return(col_matrix)
}

col_matrix <- color_mapping(cor_matrix, pvalue_matrix)

library(ComplexHeatmap)
library(circlize)

heatmap <- Heatmap(cor_matrix,
                   name = "correlation",
                   col = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "darkred")),
                   show_column_names = TRUE,
                   show_row_names = TRUE,
                   cluster_columns = T,
                   cluster_rows = T,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                  if (pvalue_matrix[i, j] < 4.16e-4) {
    text_color <- ifelse(cor_matrix[i, j] > 0, "blue", "red")
    grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, gp = gpar(col = text_color, fontsize = 6))
  }
}
     
)

svg("heatmap.test.pearson.03.18..svg")
draw(heatmap, heatmap_legend_side = "bot")
dev.off()











                
                


