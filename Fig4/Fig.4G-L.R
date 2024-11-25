

setwd('D:/LUAD/7.4.NewGSE3/4.MergeRisk_Scale/5GSE')
library(timeROC)
library(survival)

data <- read.csv("risk_scale_5GSE.NoLog_HHLL.csv")  #读取只有HHLL包含生存状态、生存时间和riskscore的CSV文件

# 提取数据中的相关列
status <- data$Status.y
full_time <- data$Time.y
predict_full <- data$Sum

# 计算ROC曲线
 ROC_rt=timeROC(T=full_time,delta=status,
                marker=predict_full,cause=1,
                weighting='aalen',
                times=c(1,2,3),ROC=TRUE)

 # 设置绘图参数
 pdf(file = "roc_5GSE.NoLog_HHLL.pdf", width = 6, height = 6)
 
 # 绘制ROC曲线
 plot(ROC_rt, time = 1, col = 'green3', title = FALSE, lwd = 2)
 plot(ROC_rt, time = 2, col = 'blue4', add = TRUE, title = FALSE, lwd = 2)
 plot(ROC_rt, time = 3, col = 'darkred', add = TRUE, title = FALSE, lwd = 2)
 
 # 添加图例
 legend('bottomright',
        c(paste0('AUC at 1 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
          paste0('AUC at 2 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
          paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[3]))),
        col = c("green3", 'blue4', 'darkred'), lwd = 2, bty = 'n')
 
 dev.off()     

