

library(ggplot2)
library(ggpubr)
library(ggExtra)
library(data.table)
# 设置工作目录
setwd("D:\\LUAD\\11.Immune\\3.ImmScore_TLSICDscore_cor")

# 读取三种免疫细胞类型的打分文件
immu_scores <- read.csv("3ImmuScores_LUAD.csv", header = TRUE)

# 读取TLS得分文件/ICDscore/TLSICDscore
tls_scores <- read.csv("TLSICDscore_ssgsea_t.csv", header = TRUE)

# 检查数据结构
str(immu_scores)
str(tls_scores)

# 根据 SampleName 对齐两个数据框
merged_data <- merge(immu_scores, tls_scores, by = "SampleName")

# 散点图函数
scatter_plot <- function(x, y, x_label, y_label) {
  df <- data.frame(x = x, y = y)
  p <- ggplot(df, aes(x = x, y = y)) +
    xlab(x_label) + ylab(y_label) +
    geom_point() + geom_smooth(method = "lm", formula = y ~ x) +
    theme_bw() +
    stat_cor(method = "spearman", aes(x = x, y = y))
  return(p)
}


# 对StromalScore和TLSscore进行相关性分析
stromal_TLSICD_cor <- cor.test(merged_data$StromalScore, merged_data$TLSICDscore, method = "spearman")
stromal_TLSICD_cor_coef <- stromal_TLSICD_cor$estimate
stromal_TLSICD_p_value <- stromal_TLSICD_cor$p.value

# 对ImmuneScore和TLSscore进行相关性分析
immune_TLSICD_cor <- cor.test(merged_data$ImmuneScore, merged_data$TLSICDscore, method = "spearman")
immune_TLSICD_cor_coef <- immune_TLSICD_cor$estimate
immune_TLSICD_p_value <- immune_TLSICD_cor$p.value

# 对ESTIMATEScore和TLSscore进行相关性分析
estimate_TLSICD_cor <- cor.test(merged_data$ESTIMATEScore, merged_data$TLSICDscore, method = "spearman")
estimate_TLSICD_cor_coef <- estimate_TLSICD_cor$estimate
estimate_TLSICD_p_value <- estimate_TLSICD_cor$p.value
# print(estimate_tls_p_value)


# 绘制StromalScore与TLSscore的散点图和边际密度图
p_stromal_TLSICD <- scatter_plot(merged_data$StromalScore, merged_data$TLSICDscore, "StromalScore", "TLSICDscore")
p_stromal_TLSICD_density <- ggMarginal(p_stromal_TLSICD, type = "density", xparams = list(fill = "orange"), yparams = list(fill = "blue"))

# 绘制ImmuneScore与TLSscore的散点图和边际密度图
p_immune_TLSICD <- scatter_plot(merged_data$ImmuneScore, merged_data$TLSICDscore, "ImmuneScore", "TLSICDscore")
p_immune_TLSICD_density <- ggMarginal(p_immune_TLSICD, type = "density", xparams = list(fill = "orange"), yparams = list(fill = "blue"))

# 绘制ESTIMATEScore与TLSscore的散点图和边际密度图
p_estimate_TLSICD <- scatter_plot(merged_data$ESTIMATEScore, merged_data$TLSICDscore, "ESTIMATEScore", "TLSICDscore")
p_estimate_TLSICD_density <- ggMarginal(p_estimate_TLSICD, type = "density", xparams = list(fill = "orange"), yparams = list(fill = "blue"))

# 
# #出图
# pdf(file="cor.pdf",width=5,height=4.8)
# print(p1)
# dev.off()
# 
# #出图2
# pdf(file="cor.density.pdf",width=5,height=5)
# print(p2)
# dev.off()
#
# 保存图像
pdf(file = "StromalScore_TLSICDscore_cor_density.pdf", width = 5, height = 5)
# print(p_stromal_TLSICD)
print(p_stromal_TLSICD_density)
dev.off()

pdf(file = "ImmuneScore_TLSICDscore_cor_density.pdf", width = 5, height = 5)
# print(p_immune_TLSICD)
print(p_immune_TLSICD_density)
dev.off()

pdf(file = "ESTIMATEScore_TLSICDscore_cor_density.pdf", width = 5, height = 5)
# print(p_estimate_TLSICD)
print(p_estimate_TLSICD_density)
dev.off()

