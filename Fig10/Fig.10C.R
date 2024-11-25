
#引用包
library(limma)
library(ggpubr)

#目录
setwd("D:/LUAD/15.TIDE/1.TIDE")

cancer=c("LUAD")

# #准备TIDE输入文件
# #读取表达文件，并对输入文件整理
# file_path <- paste("TCGA_",cancer,"_TPM.txt", sep = "")
# rt=read.table(file_path, header=T, sep="\t", check.names=F)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
# 
# #标准化，按照行计算平均值，再拿这一行所有表达量-该平均值。
# #既基因的均值为0
# Expr <- t(apply(exp, 1, function(x)x-(mean(x)))) 
# 
# write_path <- paste(cancer,"_tcga_normalize.txt", sep = "")
# write.table(data.frame(ID=rownames(Expr),Expr,check.names = F),write_path, sep="\t", quote=F, row.names = TRUE)




#读取TIDE数据
tide_path <- paste(cancer,"_TIDE_TCGA.csv", sep = "")
tide=read.csv(tide_path, header=T, row.names=1)

#仅仅保留肿瘤样本
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,,drop=F]
#修改样本名
tide[,2] = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
#去除重复
tide1=tide[!duplicated(tide$V2),]
#整理
rownames(tide1) = tide1[,2]
tide = as.data.frame(tide1[,-2])
rownames(tide) = rownames(tide1)
colnames(tide) = "TIDE"
tide=avereps(tide)

#读取风险数据文件  ##High-high,Low-low分组
risk_path <- paste(cancer,"_riskTCGA_HHLLonly.csv", sep = "")
risk=read.csv(risk_path, header=T, row.names=1)

#合并数据
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]  #将表格里hhll组的group改名为risk
data=cbind(tide, risk)

#设置比较组
data$risk=ifelse(data$risk=="High-high", "High-high", "Low-low")
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

# 按照列排序
data <- data[order(data$risk), ]


#画图

pdf(file=paste(cancer,"_","TIDE.pdf", sep = ""), width=3.5, height=4.5)

ggviolin(data, x="risk", y="TIDE", fill = "risk", 
         xlab="", ylab="TIDE",
         palette=c("#DC0000B2","#4DBBD5B2"),  # 使用红色和绿色"#E64B35B2","#00A087B2"
         legend.title="Risk",
         add = "boxplot", add.params = list(fill="white")) + 
  stat_compare_means(comparisons = my_comparisons)

dev.off()




mycol <- c("#BD6263","#8EA325","#A9D179","#84CAC0","#F5AE6B","#BCB8D3","#4387B5")

mycol2 = c( "#E64B35B2",  ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
            "#8491B4B2", "#91D1C2B2",  "#7E6148B2")


