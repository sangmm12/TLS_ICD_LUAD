

library(pROC)
library(stringr)
library(data.table)
setwd('D:/LUAD/7.ROC_AUC_new/222.SupervisedPCA_XGboost_ICD_2_Selected')
# 假设你有两个向量，分别存储了正常组和疾病组中该基因的表达量


sample_name <- "TCGA"

# genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_down_p_HR.csv",sep=''),header = T)
# genedown <- genedat$genename
# genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_down_p_HR.csv",sep=''),header = T)
# geneup <- genedat$genename
# geneupdown <- c(geneup,genedown)

# #gene_names <- c("CYP2U1","KCTD12","AP003486.1","ZBED3","IGIP","GGTA1P","SPARCL1","WAC-AS1","LINC00847","CCR2")
# genedown <- gsub("-", "_", genedown)
# gene_names <- genedown

#       #TLS
# gene_names <- c("LAMP3", "CCL2", "CCL3","CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21","CXCL9",
#                 "CXCL10", "CXCL11","CXCR4", "BCL6", "CXCL13", "CD200", "FBLN7", "ICOS", "SGPP2", "SH2D1A",
#                 "TIGIT", "PDCD1", "CD4", "CCR5", "CXCR3", "CSF2", "IGSF6", "IL2RA", "CD38", "CD40", "CD5",
#                 "MS4A1", "SDC1", "GFI1", "IL1R1", "IL1R2", "IL10", "CCL20", "IRF4", "TRAF6", "STAT5A" ,"TNFRSF17")

      #ICD
gene_names <- c("ENTPD1", "NT5E", "CALR","HMGB1", "HSP90AA1", "ATG5", "BAX", "CASP8", "PDIA3","EIF2AK3", "PIK3CA", 
                "CXCR3","IFNA1", "IFNB1", "IL10", "IL6", "TNF", "CASP1", "IL1R1", "IL1B", "NLRP3", "P2RX7", "LY96", 
                "MYD88", "TLR4", "CD4", "CD8A", "CD8B", "FOXP3", "IFNG", "IFNGR1", "IL17A", "IL17RA", "PRF1")

auc_normal_cold <- c()
out_put_list <- c()

# gene_name = "ENTPD1"
for(gene_name in gene_names){
  
  dat <- fread(paste('D:/LUAD/7.ROC_AUC_new/222.SupervisedPCA_XGboost_ICD_2_Selected/TCGA_34ICDexp_244Low.csv',sep=''), header = T)
  # dat$V1 <- gsub("-", "_", dat$V1)
  dat <- dat[, -c(ncol(dat) - 3):-ncol(dat)]  #减去后4列，group和Time等
  temp_exp <- as.data.frame(dat)
  # temp_exp <- temp_exp[which(temp_exp[,c(1)]==gene_name),]
  temp_exp <- temp_exp[,which(colnames(temp_exp)==gene_name)]
  # temp_exp <- temp_exp[,-1]
  expc <- as.numeric(as.vector(t(temp_exp)))  #expcold = 244Low group
  
  dat <- fread(paste('D:/LUAD/7.ROC_AUC_new/222.SupervisedPCA_XGboost_ICD_2_Selected/TCGA_34ICDexp_243High.csv',sep=''), header = T)
  # dat$V1 <- gsub("-", "_", dat$V1)
  dat <- dat[, -c(ncol(dat) - 3):-ncol(dat)]  #减去后4列，group和Time等
  temp_exp <- as.data.frame(dat)
  # temp_exp <- temp_exp[which(temp_exp[,c(1)]==gene_name),]
  temp_exp <- temp_exp[,which(colnames(temp_exp)==gene_name)]
  # temp_exp <- temp_exp[,-1]
  exph <- as.numeric(as.vector(t(temp_exp)))  #exphot = 243High 组
  
  normal_expression <- expc    #normal= cold Low
  disease_expression <- exph   #disese= hot high
  
  # 合并两组表达量数据，并创建对应的类别标签 (0表示正常组，1表示疾病组)
  all_expression <- c(normal_expression, disease_expression)
  labels <- factor(c(rep(0, length(normal_expression)), rep(1, length(disease_expression))))
  
  # 计算 ROC 曲线的真阳性率和假阳性率
  roc_curve <- roc(labels, all_expression)
  
  # 绘制 ROC 曲线
  # pdf(paste('normal_cold_',gene_name,'.pdf',sep=''))
  # pdf(paste('Low_',gene_name,'_TLSmodel.pdf',sep=''))  #原代码
  pdf(paste('Low_', gene_name, '_TLSmodel.pdf', sep = ''), family = "Times", pointsize = 20)
                                                                                           #AUC字体颜色
  p <- plot(roc_curve, main = paste("ROC Curve for Gene Expression of ",gene_name,sep=''), col = "#FFC107", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
  print(p)
  dev.off()
  
  #输入dev.off()直到报错
  plot(roc_curve, main = paste("",gene_name,sep=''), col = "black", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
  eval(parse(text = paste(gene_name,'<- recordPlot()')))
  out_put_list <- append(out_put_list,c(gene_name,NULL))
  
  auc_normal_cold <- c(auc_normal_cold,as.numeric(p$auc))
 
}

# pdf("auc_normal_cold.pdf",width=40,height=60)
pdf("auc_TLSmodel.pdf",width=35,height=30, family = "Times", pointsize = 20)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=7))',sep='')))

dev.off()


###上一段代码是以low作为对照组



# auc_normal_hot <- c()
# out_put_list <- c()
# 
# for(gene_name in gene_names){
#   
#   dat <- fread(paste('D:/LUAD/7.ROC_AUC_new/222.SupervisedPCA_XGboost_ICD_2_Selected/exp_',sample_name,'_cold.csv',sep=''), header = T)
#   dat$V1 <- gsub("-", "_", dat$V1)
#   temp_exp <- as.data.frame(dat)
#   temp_exp <- temp_exp[which(temp_exp[,c(1)]==gene_name),]
#   temp_exp <- temp_exp[,-1]
#   expc <- as.numeric(as.vector(t(temp_exp)))
#   
#   dat <- fread(paste('D:/LUAD/7.ROC_AUC_new/222.SupervisedPCA_XGboost_ICD_2_Selected/exp_',sample_name,'_hot.csv',sep=''), header = T)
#   dat$V1 <- gsub("-", "_", dat$V1)
#   temp_exp <- as.data.frame(dat)
#   temp_exp <- temp_exp[which(temp_exp[,c(1)]==gene_name),]
#   temp_exp <- temp_exp[,-1]
#   exph <- as.numeric(as.vector(t(temp_exp)))
#   
#   normal_expression <- exph
#   disease_expression <- expc
#   
#   # 合并两组表达量数据，并创建对应的类别标签 (0表示正常组，1表示疾病组)
#   all_expression <- c(normal_expression, disease_expression)
#   labels <- factor(c(rep(0, length(normal_expression)), rep(1, length(disease_expression))))
#   
#   # 计算 ROC 曲线的真阳性率和假阳性率
#   roc_curve <- roc(labels, all_expression)
#   
#   # 绘制 ROC 曲线
#   pdf(paste('hot/normal_hot_',gene_name,'.pdf',sep=''))
#   p <- plot(roc_curve, main = paste("ROC Curve for Gene Expression of ",gene_name,sep=''), col = "#FFC107", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
#   print(p)
#   dev.off()
#   
#   #输入dev.off()直到报错
#   plot(roc_curve, main = paste("",gene_name,sep=''), col = "#FFC107", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
#   eval(parse(text = paste(gene_name,'<- recordPlot()')))
#   out_put_list <- append(out_put_list,c(gene_name,NULL))
#   
#   
#   auc_normal_hot <- c(auc_normal_hot,as.numeric(p$auc))
#   
# }
# 
# pdf("auc_normal_hot.pdf",width=40,height=60)
# 
# eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=8))',sep='')))
# 
# dev.off()
# 
# dat <- cbind(Gene = gene_names,auc_normal_cold,auc_normal_hot)
# 
# write.csv(dat,file=paste("auc.csv",sep=''),row.names = T)

