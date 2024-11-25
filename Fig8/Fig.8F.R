

setwd("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新")

library(psych)
library(data.table)

#GSVA通路表达文件
CIBER <- fread(paste("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/19上调通路4down_hallmark.rdsExp_t.csv",sep=''),header = T)
#CIBER <- t(CIBER)
CIBER <- as.data.frame(CIBER)
rownames(CIBER) <- CIBER[,1]
CIBER <- CIBER[,-1]
dat_im <- CIBER
dat_im[] <- lapply(dat_im, as.numeric)

sample_name <- "TCGA"
data<- fread(paste('D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/TLSICDscore_ssgsea_t.csv',sep=''), header = T)
data <- as.data.frame(data)
rownames(data) <- data[,1]
dat_gene <- data[,-1,drop = FALSE] #不要删除任何维度,即使结果是一个单列或单行的数据框,也会保持为数据框的形式.

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene) ))==2))


dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]

colSums(dat_im)

data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")
# data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr") #最终横纵坐标需互换
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

write.csv(data.r,file= paste("data.r_TLSICDscore.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_TLSICDscore.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#66FFFF", "white", "#FF00FF"))(paletteLength)

test <- data.r
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.3, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.8/paletteLength, 0.6, length.out=floor(paletteLength/2)))

#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_genedown_cell_TLSICDscore.pdf",sep=''),width =3,height = 10)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()





   ####
Sel <- 'Low-low'

rf_risk <- read.csv('D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/scale_TCGA_HHLL.csv')
rf_risk <- rf_risk[,c(1,10)] #group High-high Low-low
colnames(rf_risk) <- c("SampleName","risk")
Low <- rf_risk[which(rf_risk$risk==Sel),]$SampleName

Sel <- 'High-high'

rf_risk <- read.csv('D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/scale_TCGA_HHLL.csv')
rf_risk <- rf_risk[,c(1,10)] #group High-high Low-low
colnames(rf_risk) <- c("SampleName","risk")
High <- rf_risk[which(rf_risk$risk==Sel),]$SampleName


  #### !!新加代码，赋值data_gene，否则data_gene是向量，不是数据框
data<- fread(paste('D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/TLSICDscore_ssgsea_t.csv',sep=''), header = T)
data <- as.data.frame(data)
rownames(data) <- data[,1]
dat_gene <- data[,-1,drop = FALSE] #drop = FALSE即不删除任何维度,即使结果是一个单列或单行的数据框,也会保持为数据框的形式.
 ####




  ###Low-low
dat_geneLow <- dat_gene[match(Low,rownames(dat_gene)),]  #单跑这行会报错
dat_imLow <- dat_im[match(Low,rownames(dat_im)),]

library(psych)
data.corr <- corr.test(dat_imLow, dat_geneLow, method="pearson", adjust="fdr") #横纵坐标需互换
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

write.csv(data.r,file= paste("data.r_TLSICDscore_LL.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_TLSICDscore_LL.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#66FFFF", "white", "#FF00FF"))(paletteLength)

test <- data.r
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.3, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.8/paletteLength, 0.6, length.out=floor(paletteLength/2)))

#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("TLSICDscore_LL.pdf",sep=''),width =3,height = 10)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()





     ####High-high
dat_geneHigh <- dat_gene[match(High,rownames(dat_gene)),]
dat_imHigh <- dat_im[match(High,rownames(dat_im)),]

library(psych)
data.corr <- corr.test(dat_imHigh, dat_geneHigh, method="pearson", adjust="fdr") #横纵坐标已互换
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

# na_counts <- colSums(is.na(data.r))#这段需注释？？
# data.r <- data.r[, na_counts != nrow(data.r)]#这段需注释？？
# 
# na_counts <- colSums(is.na(data.p))#这段需注释？？
# data.p <- data.p[, na_counts != nrow(data.p)]#这段需注释？？
# 


write.csv(data.r,file= paste("data.r_TLSICDscore_HH.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_TLSICDscore_HH.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#66FFFF", "white", "#FF00FF"))(paletteLength)

test <- data.r
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.3, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.8/paletteLength, 0.6, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("TLSICDscore_HH.pdf",sep=''),width =3,height = 10)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()

