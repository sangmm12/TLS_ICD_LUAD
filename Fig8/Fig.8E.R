

setwd("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新")

cancer <- "LUAD"
library(data.table)
data <-  fread(paste("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/Exp_TCGA_TLS_ICD.csv",sep=''))
data <- as.data.frame(data)
rownames(data) <- data[,1]
dat_gene <- data[,-1]


data_cell <- readRDS(paste("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/hallmark.rds",sep=''))
write.csv(data_cell,file= paste("hallmark.rdsExp.csv",sep=''),quote=F)  #??将hallmark.rds数据可视化输出为csv表，用于后续雷达图输入文件分析

DEG_BLCA <- read.csv("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/t_results_hallmark.csv",header = T)

GSVAdown <- DEG_BLCA$X[which(DEG_BLCA$t_value<= -2.58)]#只读取了下调的通路
GSVAup <- DEG_BLCA$X[which(DEG_BLCA$t_value> 2.58)]#增加一行代码，把上调的通路也读取进来
GSVAdown <- c(GSVAdown,GSVAup)#并合并到GSVAdown中
dat_GSVAdown <- data_cell[match(GSVAdown,rownames(data_cell)),]

dat_im <- t(dat_GSVAdown)

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene)))==2))

dat_im <- dat_im[match(all_name,rownames(dat_im)),]

dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]

colnames(dat_im) <- sub("HALLMARK_", "", colnames(dat_im))


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
data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")  #dat_im,dat_gene互换后，pdf的横纵坐标即可互换
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

#新加代码  #输出最后的 data.p，第一行是通路，第一列是基因
data.p_matrix <- data.corr$p
print(data.p_matrix)

#计算每列中p值小于0.05的基因数目  #将计数结果添加到data.p_matrix的最后一行
num_genes_significant <- apply(data.p_matrix < 0.05, 2, sum)
data.p_matrix_with_counts <- rbind(data.p_matrix, num_genes_significant)

 #输出最终结果
print(data.p_matrix_with_counts)
write.csv(data.p_matrix_with_counts, file = "GSVA_Gene_P_count.csv", quote = FALSE)

# paste("data.r_",cancer,".csv",sep='')
# write.csv(data.r,file= paste("data.r_xcell.csv",sep=''),quote=F)
# write.csv(data.p,file= paste("data.p_xcell.csv",sep=''),quote=F)

#截止跑到此

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
myColor <- colorRampPalette(c("#9999FF", "white", "#FF8000"))(paletteLength)

test <- data.r
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_GSVA_72gene.pdf",sep=''),width =7,height = 19.5)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat,
)

dev.off()


########################High-high Low-low

# 读入样本分组文件
sample_groups <- read.csv("scale_TCGA_HHLL.csv", header = TRUE)
# 根据分组选择样本
high_high_samples <- sample_groups$SampleName[sample_groups$group == "High-high"]
low_low_samples <- sample_groups$SampleName[sample_groups$group == "Low-low"]

# 提取对应样本的数据
dat_gene_high <- dat_gene[high_high_samples, ]
dat_gene_low <- dat_gene[low_low_samples, ]

# 提取对应样本的通路数据
dat_im_high <- dat_im[high_high_samples, ]
dat_im_low <- dat_im[low_low_samples, ]



##计算 High-high 组的相关性
data.corr_high <- corr.test(dat_gene_high, dat_im_high, method = "pearson", adjust = "fdr")
data.r_high <- data.corr_high$r
data.p_high <- data.corr_high$p
sig.mat_high <- matrix(sapply(data.p_high, getSig), nrow=nrow(data.p_high))  #新加代码，计算high组中的显著p值
str(sig.mat_high)

#新加代码  #输出最后的 data.p，第一行是通路，第一列是基因
data.p_matrix <- data.corr_high$p
print(data.p_matrix)
#计算每列中p值小于0.05的基因数目  #将计数结果添加到data.p_matrix的最后一行
num_genes_significant <- apply(data.p_matrix < 0.05, 2, sum)
data.p_matrix_with_counts <- rbind(data.p_matrix, num_genes_significant)
print(data.p_matrix_with_counts)
write.csv(data.p_matrix_with_counts, file = "GSVA_Gene_pHigh.csv", quote = FALSE)
write.csv(data.r_high,file= paste("GSVA_Gene_rHigh.csv",sep=''),quote=F)


# 计算满足【p值小于0.05且r值大于0】的基因数目
num_genes_significant_high <- apply(data.p_high < 0.05 & data.r_high > 0, 2, sum)
# 将计数结果添加到data.p_matrix的最后一行
data.p_matrix_high <- data.corr_high$p
data.p_matrix_with_counts_high <- rbind(data.p_matrix_high, num_genes_significant_high)
# 输出新表
write.csv(data.p_matrix_with_counts_high, file = "GSVA_p0.05_r0_High.csv", quote = FALSE)

# 生成 High-high 组相关性热图
pdf(paste("cor_GSVA_72gene_High-high.pdf", sep = ''), width = 7, height = 19.5)
pheatmap(data.r_high,
         color = myColor,
         breaks = myBreaks,
         clustering_method = "average", cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = sig.mat_high)
dev.off()



    ##计算 Low-low 组的相关性
data.corr_low <- corr.test(dat_gene_low, dat_im_low, method = "pearson", adjust = "fdr")
data.r_low <- data.corr_low$r
data.p_low <- data.corr_low$p
sig.mat_low <- matrix(sapply(data.p_low, getSig), nrow=nrow(data.p_low))  #新加代码，计算high组中的显著p值
str(sig.mat_low)

#输出最后的 data.p，第一行是通路，第一列是基因
data.p_matrix_low <- data.corr_low$p
print(data.p_matrix_low)
#计算每列中p值小于0.05的基因数目  #将计数结果添加到data.p_matrix的最后一行
num_genes_significant_low <- apply(data.p_matrix_low < 0.05, 2, sum)
data.p_matrix_with_counts_low <- rbind(data.p_matrix_low, num_genes_significant_low)
print(data.p_matrix_with_counts_low)
write.csv(data.p_matrix_with_counts_low, file = "GSVA_Gene_pLow.csv", quote = FALSE)
write.csv(data.r_low,file= paste("GSVA_Gene_rLow.csv",sep=''),quote=F)


# 计算满足【p值小于0.05且r值大于0】的基因数目
num_genes_significant_low <- apply(data.p_low < 0.05 & data.r_low > 0, 2, sum)
# 将计数结果添加到data.p_matrix的最后一行
data.p_matrix_low <- data.corr_low$p
data.p_matrix_with_counts_low <- rbind(data.p_matrix_low, num_genes_significant_low)
# 输出新表
write.csv(data.p_matrix_with_counts_low, file = "GSVA_p0.05_r0_Low.csv", quote = FALSE)

# 生成 Low-low 组相关性热图
pdf(paste("cor_GSVA_72gene_Low-low.pdf", sep = ''), width = 7, height = 19.5)
pheatmap(data.r_low,
         color = myColor,
         breaks = myBreaks,
         clustering_method = "average", cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = sig.mat_low)
dev.off()







#      ###
# 
# library(data.table)
# data <-  fread(paste("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/LUAD_convert_exp.txt",sep=''))
# data <- as.data.frame(data)
# 
# rownames(data) <- data[,1]
# data <- data[,-1]
# 
# genelist <- c("NR4A2","HMGB2","PIM3","ZNF331","ID2","ZFP36","NR4A3","MALAT1","LYZ","FTL",
#               "LGMN","NFKB1","JUNB","HLA-E","BCL2A1","PLA2G7","DUSP2")
# 
# dat_gene <- data[match(genelist,rownames(data)),]
# 
# dat_gene=as.data.frame(lapply(dat_gene,as.numeric),check.names=F)
# 
# rownames(dat_gene) <- genelist
# 
# dat_gene <- t(dat_gene)
# 
# 
# 
# 
# 
# data_cell <- readRDS(paste("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/hallmark.rds",sep=''))
# 
# DEG_BLCA <- read.csv("D:/1.LUAD/14.GSVA/1.GSVA/11.13改_新/t_results_hallmark.csv",header = T)
# 
# GSVAdown <- DEG_BLCA$X[which(DEG_BLCA$t_value<= -2.58)]
# 
# dat_GSVAdown <- data_cell[match(GSVAdown,rownames(data_cell)),]
# 
# dat_im <- t(dat_GSVAdown)
# 
# all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene)))==2))
# 
# dat_im <- dat_im[match(all_name,rownames(dat_im)),]
# 
# dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
# 
# colnames(dat_im) <- sub("HALLMARK_", "", colnames(dat_im))
# 
# # for(i in 1:length(colnames(dat_gene)))
# # {
# #   dat_gene[,i] = as.numeric(unlist(dat_gene[,i]))
# # }
# # i=1
# # nrow(dat_gene)
# # ncol(dat_gene)
# # 
# 
# colSums(dat_im)
# 
# 
# 
# library(psych)
# data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
# data.r <- data.corr$r  # 相关系数
# data.p <- data.corr$p  # p值
# 
# paste("data.r_",cancer,".csv",sep='')
# write.csv(data.r,file= paste("data.r_gene.csv",sep=''),quote=F)
# write.csv(data.p,file= paste("data.p_gene.csv",sep=''),quote=F)
# 
# 
# library(pheatmap)
# getSig <- function(dc) {
#   print(dc)
#   sc <- ' '
#   if (dc < 0.0001) {sc <- '****'}
#   else if (dc < 0.001){sc <- '***'}
#   else if (dc < 0.01){sc <- '**'}
#   else if (dc < 0.05) {sc <- '*'}
#   else{sc <- ''
#   }
#   return(sc)
# }
# sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
# str(sig.mat)
# 
# paletteLength <- 1000
# myColor <- colorRampPalette(c("#9999FF", "white", "#FF8000"))(paletteLength)
# 
# test <- data.r
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
# #myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
# #seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))
# 
# 
# #chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...
# 
# pdf(paste("cor_KM2_gene.pdf",sep=''),width =9,height = 8)
# pheatmap(data.r, 
#          color=myColor,
#          breaks=myBreaks,
#          clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
# )
# 
# dev.off()
# 
