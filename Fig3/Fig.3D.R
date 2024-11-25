

library(data.table)
library(survival)
library(survminer)
setwd('D:/LUAD/77.ReRiskGroup2_selected/TCGA')
data <- read.csv("scale_TCGA_HHLL.csv")  # 替换"your_data.csv"为你的CSV文件路径

# 对数据进行排序
rt <- data[order(data$Sum),]  #Sum值是两个riskscore归一化后相加的riskscore
rt <- rt[order(rt$group), ]

# 进行生存分析
diff <- survdiff(Surv(Time, Status) ~ group, data = rt)
pValue <- 1 - pchisq(diff$chisq, df = 1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific = TRUE)

fit <- survfit(Surv(Time, Status) ~ group, data = rt)

# 生成生存曲线图
surPlot <- ggsurvplot(fit,
                      data = rt,
                      legend.labs = c(unique(rt$group)),
                      legend = "top",
                      legend.title = "Risk",
                      pval = paste0("p=", pValue),
                      pval.size = 5,
                      xlab = "Time (years)",
                      break.time.by = ceiling((max(rt$Time, na.rm = TRUE))/4),
                      risk.table.title = "",
                      palette = c("red", "green"),
                      risk.table = T,
                      risk.table.height = 0.25,
                      extend = 0.05)  # 添加extend参数以延长时间范围

# 4. 保存生存曲线图为PDF
pdf(file = "TCGA_survival_HHLL.pdf", onefile = FALSE, width = 6, height = 6)
print(surPlot)
dev.off()
print(surPlot)
