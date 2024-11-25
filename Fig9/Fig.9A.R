######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
scoreFile="ImmuScores.txt"       #微环境打分文件
cluFile="TCGA_HHLL.txt"        #分型的结果文件
setwd("D:/LUAD/11.Immune/1.Stromal_Immune_ESTIMATEScore")     #设置工作目录

#读取肿瘤微环境打分文件，并对数据进行整理
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
# row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
# score=avereps(score)

#读取聚类文件，并对输入文件整理
cluster=read.table(cluFile, header=F, sep="\t", row.names=1, check.names=F)

#样品取交集
sameSample=intersect(row.names(score), row.names(cluster))
score1=score[sameSample,,drop=F]
cluster1=cluster[sameSample,,drop=F]
colnames(cluster1)=c("Cluster")
data=cbind(score1, cluster1)
data$Cluster=paste0("Cluster", data$Cluster)

#设置比较组
type=levels(factor(data[,"Cluster"]))
data$Cluster=factor(data$Cluster, levels=type)
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

# #设置颜色
# bioCol=c("#00BFC4", "#F8766D","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
# bioCol=bioCol[1:length(levels(factor(data$Cluster)))]
# 
# #肿瘤微环境差异分析
# for(i in colnames(data)[1:(ncol(data)-1)]){
#   #绘制箱线图
#   boxplot=ggboxplot(data, x="Cluster", y=i, fill="Cluster",
#                     xlab="",
#                     ylab=i,
#                     legend.title="Cluster",
#                     palette=bioCol,add = "jitter"
#   )+ 
#     stat_compare_means(comparisons=my_comparisons)
#   
# 
#    #输出图形
#   pdf(file=paste0(i, ".pdf"), width=5, height=4.5)
#   print(boxplot)
#   dev.off()
# }


# 设置颜色  #High-high、Low-low
bioCol <- c( "#F8766D","#00BFC4","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol <- bioCol[1:length(levels(factor(data$Cluster)))]

# 肿瘤微环境差异分析
for(i in colnames(data)[1:(ncol(data)-1)]) {
  # 绘制箱线图
  boxplot <- ggboxplot(data, x = "Cluster", y = i, fill = "Cluster",
                       xlab = "",
                       ylab = i,
                       legend.title = "Cluster",
                       palette = bioCol
  ) + 
    stat_compare_means(comparisons = my_comparisons) +
    geom_jitter(aes(color = Cluster), width = 0.2)
  
  # 输出图形
  pdf(file = paste0(i, ".pdf"), width = 4, height = 4.5)
  print(boxplot)
  dev.off()  
}

