

setwd('D:/LUAD/17.ssGSVA/2.12')
library(data.table)

#Tumor features
genelist <- c("MKI67","ESCO2","CETN3","CDK2","CCND1","CCNE1","AURKA","AURKB","E2F1","MYBL2","BUB1","PLK1","CCNB1","MCM2","MCM6",
              "SNAI1","SNAI2","TWIST1","TWIST2","ZEB1","ZEB2","CDH2")





# col1<- rep("T_cell", times=length(genelist))  
col1<- rep("Tumor features", times=length(genelist))

col2<- genelist
immune <- data.frame(col1,col2)
colnames(immune)  <- c("gs_name","gene_symbol")


gcSample = split(immune$gene_symbol,
                 immune$gs_name)

#读取exp1
exp1 <-  fread(paste("D:/LUAD/17.ssGSVA/2.12/LUAD_Count_Exp.txt",sep=''))
exp1 <- as.data.frame(exp1)

#AI改的代码：将第一列基因名改为行名
row.names(exp1) <- exp1[, 1]
exp1 <- exp1[, -1]

data_merge <- exp1   #第一列基因，第一行样本。TCGA的Count数据，几万个基因

data_merge=apply(data_merge,2,as.numeric)

rownames(data_merge) <- rownames(exp1)

expr1 <- as.matrix(data_merge)

#基因集需要是list为对象
#默认情况下，kcdf="Gaussian"，适用于输入表达式值连续的情况，如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs或log-TPMs。
#当输入表达式值是整数Count计数时，比如那些从RNA-seq实验中得到的值，那么这个参数应该设置为kcdf="Poisson"

library(GSVA)
ssgsea<- gsva(expr1,gcSample, method='ssgsea', kcdf='Poisson',abs.ranking=TRUE)
# 转置结果
ssgsea_transposed <- t(ssgsea)

write.csv(ssgsea_transposed,file = paste("Tumor features_score_t.csv",sep=''))

