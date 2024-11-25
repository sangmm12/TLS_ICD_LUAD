
setwd('D:/LUAD/11.Immune/6.HighHigh_LowLow_link/p')

library(linkET)
#library(devtools)
#devtools::install_github ("Hy4m/linkET", force = TRUE)
library(ggplot2)
library(dplyr)
library(GSVA)
library(data.table)

#'hot-low'替换为Low-low
# Sel <- 'hot-low'
Sel <- 'Low-low'

rf_risk <- read.csv('D:/LUAD/11.Immune/6.HighHigh_LowLow_link/scale_TCGA_HHLL.csv')
rf_risk <- rf_risk[,c(1,10)]  #第10行hhll分组
colnames(rf_risk) <- c("SampleName","risk")
# hot_low <- rf_risk[which(rf_risk$risk==Sel),]$SampleName
Low_low <- rf_risk[which(rf_risk$risk==Sel),]$SampleName  

#'cold-high'替换为High-high
# Sel <- 'cold-high'
Sel <- 'High-high'

rf_risk <- read.csv('D:/LUAD/11.Immune/6.HighHigh_LowLow_link/scale_TCGA_HHLL.csv')
rf_risk <- rf_risk[,c(1,10)]
colnames(rf_risk) <- c("SampleName","risk")
# cold_high <- rf_risk[which(rf_risk$risk==Sel),]$SampleName
High_high <- rf_risk[which(rf_risk$risk==Sel),]$SampleName


#########读取CIBER算法免疫细胞exp  ##此时样本名已删除为4部分
CIBER <- fread(paste("D:/LUAD/11.Immune/6.HighHigh_LowLow_link/im_cibersort_TCGA_LUAD_NoPval.csv",sep=''),header = T)
colnames(CIBER)[1] <- 'SampleName'
CIBER <- CIBER[!duplicated(CIBER$SampleName), ] #删除重复样本名

######读取#riskscore   #选取归一化后相加的riskscore
sample_name <- "TCGA"
riskscore <- fread("D:/LUAD/11.Immune/6.HighHigh_LowLow_link/scale_TCGA_HHLL.csv",header = T)
riskscore <- as.data.frame(riskscore)
dat_riskscore <- riskscore[,c(1,13)] #选取归一化后相加的第13行riskscore
colnames(dat_riskscore)[1] <- 'SampleName'

###################ssGSVA  ##TLSICDscore
sample_name <- "TCGA"
dat <- fread(paste('D:/LUAD/11.Immune/6.HighHigh_LowLow_link/TLSICDscore_ssgsea.csv',sep=''), header = T)
dat <- as.data.frame(dat)
dat <- dat[!duplicated(dat[,c(1)]),]
rownames(dat) <- dat[,c(1)]
dat <- dat[,-c(1)]
dat <- na.omit(dat)

###新##
# 读取TLS得分  #输入文件需要改转置，并加上SampleName
tls_scores <- read.csv("D:/LUAD/11.Immune/6.HighHigh_LowLow_link/TLSscore_ssgsea_t.csv")
colnames(tls_scores) <- c("SampleName", "TLSscore") #给矩阵的变量命名
                                                    #将数据框tls_scores的列名分别设为"SampleName"和"TLSscore"。
# 读取ICD得分  #输入文件需要改转置，并加上SampleName
icd_scores <- read.csv("D:/LUAD/11.Immune/6.HighHigh_LowLow_link/ICDscore_ssgsea_t.csv")
colnames(icd_scores) <- c("SampleName", "ICDscore")


#hot_low  ##Low-low        #ssgsea可视为是TLSICDscore
######################################################
# ssgsea <- dat[,which(colnames(dat)%in%hot_low)]
ssgsea <- dat[,which(colnames(dat)%in%Low_low)]
ssgsea <-as.data.frame(t(ssgsea))
ssgsea$SampleName <- rownames(ssgsea)

ssg <- ssgsea
colnames(ssg)[1] <- "TLSICDscore"
write.csv(ssg,file= paste("TLSICDscore_LL_data.csv",sep=''),quote=F)


out <- merge(CIBER,ssgsea,by='SampleName')
out <- merge(out,dat_riskscore,by='SampleName')
out <- merge(out, tls_scores, by="SampleName") #新，合并数据
out <- merge(out, icd_scores, by="SampleName")

riskscore <- data.frame(risk = out$riskscore)
TLSICDscore <- data.frame(risk = out$TLSICD)
TLSscore <- data.frame(risk = out$TLSscore)#新，合并数据
ICDscore <- data.frame(risk = out$ICDscore)

risks <- riskscore
rownames(risks) <- out$SampleName
write.csv(risks,file= paste("riskscore_LL_data.csv",sep=''),quote=F)

#新输出
TLS_scores <- TLSscore
rownames(TLS_scores) <- out$SampleName
write.csv(TLS_scores,file= paste("TLSscore_LL_data.csv",sep=''),quote=F)

ICD_scores <- ICDscore
rownames(ICD_scores) <- out$SampleName
write.csv(ICD_scores,file= paste("ICDscore_LL_data.csv",sep=''),quote=F)


rownames(out) <- out$SampleName
out <- subset(out,select=-c(riskscore,SampleName,TLSscore,ICDscore,TLSICD)) ###out是只留免疫细胞的表达
out[] <- lapply(out, as.numeric)

data.corr <- qcorrplot(correlate(out), type = "lower", diag = FALSE)
data_cor <- data.corr$data  # 相关系数
write.csv(data_cor,file= paste("Imm.cor_LL.csv",sep=''),quote=F)


mantel1 <- mantel_test(riskscore, out,
                       spec_select = list(riskscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色


mantel2 <- mantel_test(TLSICDscore, out,
                       spec_select = list(TLSICDscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色

# mantel <- rbind(mantel1,mantel2)

#新 把TLSscore和ICDscore加入新赋值点mental3和mental4
mantel3 <-  mantel_test(TLSscore, out,
                        spec_select = list(TLSscore = 1:1)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# 对ICDscore进行Mantel测试
mantel4  <-  mantel_test(ICDscore, out,
                         spec_select = list(ICDscore = 1:1)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# 合并所有Mantel测试结果
mantel  <-  rbind(mantel1, mantel2, mantel3, mantel4)

write.csv(mantel,file='4种score.Immu_LL.csv',row.names = F)


pdf(paste('linkET_Lowlow.pdf',sep=''),width=length(colnames(CIBER))/2.6,height=length(colnames(CIBER))/2.6)
qcorrplot(correlate(out), type = "upper", diag = FALSE) +#热图绘制
  geom_square() +#热图绘制
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes里面是线条格式，data对应的是mantel test 计算结果，curvature控制线条曲率
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-0.5, 0.5)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides()函数调整标签样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()




   ###cold_high  #High-high  #注意ssgsea的赋值，不能只是一列数字，需要有样本名。
#ssgsea可视为是TLSICDscore 
ssgsea <- dat[,which(colnames(dat)%in%High_high)]
ssgsea <-as.data.frame(t(ssgsea))
ssgsea$SampleName <- rownames(ssgsea)

ssg <- ssgsea
colnames(ssg)[1] <- "ssgsea"
write.csv(ssg,file= paste("TLSICDscore_HH_data.csv",sep=''),quote=F)

out <- merge(CIBER,ssgsea,by='SampleName')
out <- merge(out,dat_riskscore,by='SampleName')
out <- merge(out, tls_scores, by="SampleName") #新，合并数据
out <- merge(out, icd_scores, by="SampleName")

riskscore <- data.frame(risk = out$riskscore)
TLSICDscore <- data.frame(risk = out$TLSICD)
TLSscore <- data.frame(risk = out$TLSscore)#新，合并数据
ICDscore <- data.frame(risk = out$ICDscore)

risks <- riskscore
rownames(risks) <- out$SampleName
write.csv(risks,file= paste("riskscore_HH_data.csv",sep=''),quote=F)

#新输出
TLS_scores <- TLSscore
rownames(TLS_scores) <- out$SampleName
write.csv(TLS_scores,file= paste("TLSscore_HH_data.csv",sep=''),quote=F)

ICD_scores <- ICDscore
rownames(ICD_scores) <- out$SampleName
write.csv(ICD_scores,file= paste("ICDscore_HH_data.csv",sep=''),quote=F)


rownames(out) <- out$SampleName
out <- subset(out,select=-c(riskscore,SampleName,TLSICD,TLSscore,ICDscore))
out[] <- lapply(out, as.numeric)
###################
data.corr <- qcorrplot(correlate(out), type = "lower", diag = FALSE)
data_cor <- data.corr$data  # 相关系数

write.csv(data_cor,file= paste("Immu.cor_HH.csv",sep=''),quote=F)



mantel1 <- mantel_test(riskscore, out,
                       spec_select = list(riskscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色


mantel2 <- mantel_test(TLSICDscore, out,
                       spec_select = list(TLSICDscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色

# mantel <- rbind(mantel1,mantel2)

#新 把TLSscore和ICDscore加入新赋值点mental3和mental4
mantel3 <-  mantel_test(TLSscore, out,
                        spec_select = list(TLSscore = 1:1)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# 对ICDscore进行Mantel测试
mantel4  <-  mantel_test(ICDscore, out,
                         spec_select = list(ICDscore = 1:1)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# 合并所有Mantel测试结果
mantel  <-  rbind(mantel1, mantel2, mantel3, mantel4)
write.csv(mantel,file='4种score.Immu_HH.csv',row.names = F)


pdf(paste('linkET_Highhigh.pdf',sep=''),width=length(colnames(CIBER))/2.6,height=length(colnames(CIBER))/2.6)
qcorrplot(correlate(out), type = "lower", diag = FALSE) +#热图绘制
  geom_square() +#热图绘制
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes里面是线条格式，data对应的是mantel test 计算结果，curvature控制线条曲率
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-0.5, 0.5)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides()函数调整标签样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()

