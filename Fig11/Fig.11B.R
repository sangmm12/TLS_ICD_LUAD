
##药物和免疫细胞的相关性分析
setwd("D:/LUAD/15.3.drugPredict_5up23downHHLLdiff/1.drug1_GDSC/zong/2.Immune_drug_cor")
library(data.table)

sample_name <- "TCGA"
dat <- fread(paste("D:/LUAD/15.3.drugPredict_5up23downHHLLdiff/1.drug1_GDSC/zong/2.Immune_drug_cor/out_put_",sample_name,"_DEGs_zong.csv",sep=''))
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat <- dat[,-1]

##读取上一个文件夹路径中筛选的药物和基因的data.r
datar <- fread(paste("D:/LUAD/15.3.drugPredict_5up23downHHLLdiff/1.drug1_GDSC/zong/2.Immune_drug_cor/data.r60.csv",sep=''))
drug <- colnames(datar)[-1]

dat_gene <- dat[,match(drug,colnames(dat))]

is.numeric(dat_gene)

CIBER <- fread(paste("D:/LUAD/15.3.drugPredict_5up23downHHLLdiff/1.drug1_GDSC/zong/2.Immune_drug_cor/TCGA_LUAD_Cibersort_NoTcellsCD4.csv",sep=''),header = T)
CIBER <- t(CIBER)
CIBER <- as.data.frame(CIBER)
colnames(CIBER) <- CIBER[1,]
CIBER <- CIBER[-1,]
dat_im <- CIBER
dat_im[] <- lapply(dat_im, as.numeric)
dat_im <- t(dat_im)     ##转置

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene) ))==2))


dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]

# for(i in 1:length(colnames(dat_gene)))
# {
#   dat_gene[,i] = as.numeric(unlist(dat_gene[,i]))
# }
# i=1
# nrow(dat_gene)
# ncol(dat_gene)
# 

colSums(dat_im)


library(psych)
data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
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
#myBreaks <- c(seq(-0.9, 0, length.out = floor(paletteLength/2) + 1))
myBreaks <- c(seq(-0.4, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.6/paletteLength, 0.4, length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_GDSC3drug_cell.pdf",sep=''),width =4,height = 8)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()
