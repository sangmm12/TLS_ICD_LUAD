

setwd("D:/LUAD/11.Immune/2.ImmScore_72gene_Cor")

library(psych)
library(data.table)

#三种免疫类型打分StromalScore、ImmuneScore、ESTIMATEScore
CIBER <- fread(paste("D:/LUAD/11.Immune/2.ImmScore_72gene_Cor/3ImmuScores_LUAD.csv",sep=''),header = T)
#CIBER <- t(CIBER)
CIBER <- as.data.frame(CIBER)
rownames(CIBER) <- CIBER[,1]
CIBER <- CIBER[,-1]
dat_im <- CIBER
dat_im[] <- lapply(dat_im, as.numeric)

sample_name <- "TCGA"
data<- fread(paste('D:/LUAD/11.Immune/2.ImmScore_72gene_Cor/Exp_',sample_name,'_TLS_ICD.csv',sep=''), header = T)
data <- as.data.frame(data)
rownames(data) <- data[,1]
dat_gene <- data[,-1]

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene) ))==2))


dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]


colSums(dat_im)


data.corr <- corr.test( dat_gene, dat_im,method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p.csv",sep=''),quote=F)


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
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.3, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.8/paletteLength, 0.6, length.out=floor(paletteLength/2)))

#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_genedown_cell.pdf",sep=''),width =2.5,height = 22)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()



#####################################
Sel <- 'Low-low'

rf_risk <- read.csv('D:/LUAD/11.Immune/2.ImmScore_72gene_Cor/scale_TCGA_HHLL.csv')
rf_risk <- rf_risk[,c(1,10)] #group High-high Low-low
colnames(rf_risk) <- c("SampleName","risk")
Low <- rf_risk[which(rf_risk$risk==Sel),]$SampleName

Sel <- 'High-high'

rf_risk <- read.csv('D:/LUAD/11.Immune/2.ImmScore_72gene_Cor/scale_TCGA_HHLL.csv')
rf_risk <- rf_risk[,c(1,10)] #group High-high Low-low
colnames(rf_risk) <- c("SampleName","risk")
High <- rf_risk[which(rf_risk$risk==Sel),]$SampleName


################################################
dat_geneLow <- dat_gene[match(Low,rownames(dat_gene)),]
dat_imLow <- dat_im[match(Low,rownames(dat_im)),]


library(psych)
data.corr <- corr.test( dat_geneLow, dat_imLow,method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

write.csv(data.r,file= paste("data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p.csv",sep=''),quote=F)


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
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.3, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.8/paletteLength, 0.6, length.out=floor(paletteLength/2)))

#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_genedown_Low_cell.pdf",sep=''),width =2.5,height = 22)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()



######################################
dat_geneHigh <- dat_gene[match(High,rownames(dat_gene)),]
dat_imHigh <- dat_im[match(High,rownames(dat_im)),]


library(psych)
data.corr <- corr.test( dat_geneHigh, dat_imHigh,method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

na_counts <- colSums(is.na(data.r))
data.r <- data.r[, na_counts != nrow(data.r)]

na_counts <- colSums(is.na(data.p))
data.p <- data.p[, na_counts != nrow(data.p)]

write.csv(data.r,file= paste("data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p.csv",sep=''),quote=F)


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
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.3, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.8/paletteLength, 0.6, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_genedown_High_cell.pdf",sep=''),width =2.5,height = 22)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()

