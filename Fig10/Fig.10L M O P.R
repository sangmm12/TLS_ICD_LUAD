######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="GSE63885_CutoffGroup.txt"    #TLSICD的ssGSVA打分cutoff分组文件文件
cliFile="GSE63885_Response.txt"            #临床数据文件
trait="Response"                    #临床性状
setwd("D:/LUAD/20.1.ImmunetherapyGSE/GSE63885_re")     #设置工作目录

#读取输入文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,,drop=F], cli[sameSample,,drop=F])

#定义临床性状的颜色
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]

#统计高低评分组病人数目
rt1=rt[,c(trait, "group")]
colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))
#计算高低评分组的百分率
df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
#百分比位置
df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("Low", "High"))



    ###新加代码，对高低风险组中的应答比例进行Fisher检验
# 提取低风险组和高风险组中的应答人数
low_responders <- df$Freq[df$group == "Low" & df$trait == levels(df$trait)[1]]
low_total <- sum(df$Freq[df$group == "Low"])
high_responders <- df$Freq[df$group == "High" & df$trait == levels(df$trait)[1]]
high_total <- sum(df$Freq[df$group == "High"])

# 进行Fisher精确检验
fisher_pvalue <- fisher.test(matrix(c(low_responders, low_total - low_responders, high_responders, high_total - high_responders), 
                                    nrow = 2))$p.value

# 绘制百分率图
p <- ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values = bioCol) +
  xlab("TLSICDscore") + ylab("Percent weight") + guides(fill = guide_legend(title = trait)) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_bw() +
  # annotate("text", x = 1.5, y = max(df$percent) + 5, label = paste0("Fisher's p = ", signif(fisher_pvalue, digits = 3)), size = 4)
  annotate("text", x = 2, y = max(df$percent) + 5, label = paste0("Fisher's p = ", signif(fisher_pvalue, digits = 3)), size = 4, hjust = 1)


# 保存图形到PDF文件
pdf(file = "barplot.GSE63885_response.pdf", width = 4, height = 5)
print(p)
dev.off()


 ##
#设置比较组
rt2=rt[,c(trait, "TLSICDscore")]
colnames(rt2)=c("trait", "TLSICDscore")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#绘制箱线图
boxplot=ggboxplot(rt2, x="trait", y="TLSICDscore", fill="trait",
                  xlab=trait,
                  ylab="TLSICDscore",
                  legend.title=trait,
                  palette=bioCol
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file="boxplot.GSE63885_Diff.pdf",width=4,height=4.5)
print(boxplot)
dev.off()

