


setwd("D:/LUAD/21.0.DrugPredict/1.drug1_GDSC/zong/1.cor_gene")
sample_name <- "TCGA"
library(data.table)
dat <- fread(paste("D:/LUAD/21.0.DrugPredict/1.drug1_GDSC/out_put_",sample_name,"_DEGs_zong.csv",sep=''))
dat <- as.data.frame(dat)

rownames(dat) <- dat[,1]
dat_gene <- dat[,-1]

is.numeric(dat_gene)

# ##保证前13都是up基因，后12down基因
# geneupdown <- c("HSP90AA1","IL1R2","PIK3CA","ATG5","IFNA1","EIF2AK3","PDIA3","HMGB1","CALR","CCL20","IL6","CXCL10","CCL8","CCL19","CD5","MS4A1","IRF4","IGSF6","TNFRSF17","SH2D1A","ICOS","CXCR3","CSF2","LAMP3","CXCL13")
#保证差异基因  前5up，后23down 
geneupdown <- c("CCL20","HSP90AA1","IFNA1","IL1R2","IL6","CASP1","CCL18","CCL19","CCR5","CD4","CD40","CD5","CSF2",
                "CXCL13","CXCR3","CXCR4","FOXP3","ICOS","IGSF6","IRF4","LAMP3","MS4A1","NLRP3","P2RX7","SDC1","SH2D1A",
                "TIGIT","TNFRSF17")
data <- fread(paste('D:/LUAD/21.0.DrugPredict/1.drug1_GDSC/149HH200LL_5up23downDiffExp_TCGA.csv',sep=''), header = T)
data <- as.data.frame(data)
rownames(data) <- data[,1]
dat_im <- data[,-1]    #新代码: 赋值dat_im
dat_im <- t(dat_im) #转置为行名是gene

dat_im <- dat_im[match(geneupdown,rownames(dat_im)),]
dat_im <- t(dat_im)

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene) ))==2))

dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]

colSums(dat_im)
library(psych)
data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

write.csv(data.r,file= paste("data.r_All.csv",sep=''))         #,quote=F 注释掉，否则输出的第一行药物名会以逗号分割
write.csv(data.p,file= paste("data.p_All.csv",sep=''))         #,quote=F
save(data.r,data.p, file = "data.r_p_All.rda")             ##另加代码，保存data.r和data.p  99行时加载

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
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))

#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_",sample_name,"_DEGs_Drug_All.pdf",sep=''),width =80,height = 13)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()



         ###另加代码   加载第60行保存的药物和基因的相关性r和p值data
    load("data.r_p_All.rda")         

            ##GJ另一个脚本中复制的一行代码
##计算前5个up基因中，60%比例药物和5up基因呈负相关；同时后23个down基因中，60%比例药物和23down基因正相关。
#计算前5个up基因data.r<0的比例 
condition1 <- apply(data.r[1:5, ] < 0, 2, function(x) mean(x, na.rm = TRUE) >= 0.6)
#计算14-25个基因data.r>0的比例
condition2 <- apply(data.r[6:28, ] > 0, 2, function(x) mean(x, na.rm = TRUE) >= 0.6)

# 输出符合条件的列  
drug <- names(which(condition1 & condition2))
print(drug)   #80%比例1个药物   #60%比例有3个药物  #70%比例?个药物 

data.r <- data.r[,match(drug,colnames(data.r))]
data.p <- data.p[,match(drug,colnames(data.p))]
write.csv(data.r,file= paste("data.r60.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p60.csv",sep=''),quote=F)
save(data.r,data.p, file = "data.r_p60.rda")

                      ###画图###
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
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.4, 0, length.out=ceiling(paletteLength/2) + 1), 
seq(0.4/paletteLength, 0.8, length.out=floor(paletteLength/2)))

#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_",sample_name,"_DEGs_3Drug_60per.pdf",sep=''),width =2.35,height = 9)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()




#                     ##改师弟的代码
#                   ##药物和基因正相关
#   TLSICDgene_name <- data.r
#                   # 计算每列大于0的比例
#   column_percentages <- colSums(TLSICDgene_name > 0) / nrow(TLSICDgene_name)
#                   # 筛选出大于85%数据大于0的列名
#   selected_columns <- colnames(TLSICDgene_name)[column_percentages > 0.85]
#   drugup <-selected_columns
#   print(drugup)
# writeLines(drugup,"7drugup_85%正相关.txt")
#    data.r0.85 <- data.r[,match(drugup,colnames(data.r))]
#    data.p0.85 <- data.p[,match(drugup,colnames(data.p))]
# write.csv(data.r0.85,file= paste("data.r85%正相关.csv",sep=''),quote=F)
# write.csv(data.p0.85,file= paste("data.p85%正.csv",sep=''),quote=F)
# save(data.r0.85,data.p0.85, file = "data.r_p_0.85Pos.rda")          ##另加代码:保存data.r和data.p 85%比例药物和基因正相关
# 
# 
# load("data.r_p_All.rda") ##另加代码
#                    # 计算每列大于0的比例
#    column_percentages <- colSums(TLSICDgene_name > 0) / nrow(TLSICDgene_name)
#                    # 筛选出大于90%数据大于0的列名
#    selected_columns <- colnames(TLSICDgene_name)[column_percentages > 0.9]
#    drugup <-selected_columns
#    print(drugup)
# writeLines(drugup,"2drugup_90%正相关.txt")
#    data.r0.9 <- data.r[,match(drugup,colnames(data.r))]
#    data.p0.9 <- data.p[,match(drugup,colnames(data.p))]
# write.csv(data.r0.9,file= paste("data.r90%正相关.csv",sep=''),quote=F)
# write.csv(data.p0.9,file= paste("data.p90%正.csv",sep=''),quote=F)
# save(data.r0.9,data.p0.9, file = "data.r_p_0.9Pos.rda")          ##另加代码，保存data.r和data.p 90%比例药物和基因正相关
# 
# 
#    ###再筛选一次 80%比例
# load("data.r_p_All.rda") ##另加代码
# # 计算每列大于0的比例
# column_percentages <- colSums(TLSICDgene_name > 0) / nrow(TLSICDgene_name)
# # 筛选出大于80%数据大于0的列名
# selected_columns <- colnames(TLSICDgene_name)[column_percentages > 0.8]
# drugup <-selected_columns
# print(drugup)
# writeLines(drugup,"12drugup_80%正相关.txt")
# data.r0.8 <- data.r[,match(drugup,colnames(data.r))]
# data.p0.8 <- data.p[,match(drugup,colnames(data.p))]
# write.csv(data.r0.8,file= paste("data.r80%正相关.csv",sep=''),quote=F)
# write.csv(data.p0.8,file= paste("data.p80%正.csv",sep=''),quote=F)
# save(data.r0.8,data.p0.8, file = "data.r_p_0.8Pos.rda")          ##另加代码，保存data.r和data.p 80%比例药物和基因正相关
# 
# 
# 
# 
# load("data.r_p_All.rda") ##另加代码   
#                   ### 药物和基因负相关
# # all_nameFedown <- names(which(table(c(rownames(data.r),dataFe1 ))==2))
# # dataFedown <- data.r[match( all_nameFedown,rownames(data.r)),]
#                  #计算每列小于0的比例
#   TLSICDgene_name <- data.r
#   column_percentages <- colSums(TLSICDgene_name < 0) / nrow(TLSICDgene_name)
#                  #筛选出大于90%数据大于0的列名
#   selected_columns <- colnames(TLSICDgene_name)[column_percentages > 0.90]
#   drugdown <-selected_columns
#   print(drugdown)
# writeLines(drugdown,"4drugdown_90%负相关.txt")
#   data.r0.9 <- data.r[,match(drugdown,colnames(data.r))]
#   data.p0.9 <- data.p[,match(drugdown,colnames(data.p))]
# write.csv(data.r0.9,file= paste("data.r90%负相关.csv",sep=''),quote=F)
# write.csv(data.p0.9,file= paste("data.p90%负.csv",sep=''),quote=F)
# save(data.r0.9,data.p0.9, file = "data.r_p_0.90Neg.rda")          ##另加代码，保存data.r和data.p 90%比例药物和基因负相关
# 
# 
# load("data.r_p_All.rda") ##另加代码   
#                    ### 药物和基因负相关
# # all_nameFedown <- names(which(table(c(rownames(data.r),dataFe1 ))==2))
# # dataFedown <- data.r[match( all_nameFedown,rownames(data.r)),]
#                   #计算每列小于0的比例
# TLSICDgene_name <- data.r
# column_percentages <- colSums(TLSICDgene_name < 0) / nrow(TLSICDgene_name)
#                   #筛选出大于80%数据大于0的列名
# selected_columns <- colnames(TLSICDgene_name)[column_percentages > 0.80]
# drugdown <-selected_columns
# print(drugdown)
# writeLines(drugdown,"37drugdown_80%负相关.txt")
# data.r0.8 <- data.r[,match(drugdown,colnames(data.r))]
# data.p0.8 <- data.p[,match(drugdown,colnames(data.p))]
# write.csv(data.r0.8,file= paste("data.r80%负相关.csv",sep=''),quote=F)
# write.csv(data.p0.8,file= paste("data.p80%负.csv",sep=''),quote=F)
# save(data.r0.8,data.p0.8, file = "data.r_p_0.8Neg.rda")          ##另加代码，保存data.r和data.p 80%比例药物和基因负相关
# 
# 
# load("data.r_p_All.rda") ##另加代码   
#               #计算每列小于0的比例
# TLSICDgene_name <- data.r
# column_percentages <- colSums(TLSICDgene_name < 0) / nrow(TLSICDgene_name)
#           #筛选出大于85%数据大于0的列名
# selected_columns <- colnames(TLSICDgene_name)[column_percentages > 0.85]
# drugdown <-selected_columns
# print(drugdown)
# writeLines(drugdown,"16drugdown_85%负相关.txt")
# data.r0.85 <- data.r[,match(drugdown,colnames(data.r))]
# data.p0.85 <- data.p[,match(drugdown,colnames(data.p))]
# write.csv(data.r0.85,file= paste("data.r85%负相关.csv",sep=''),quote=F)
# write.csv(data.p0.85,file= paste("data.p85%负.csv",sep=''),quote=F)
# save(data.r0.85,data.p0.85, file = "data.r_p_0.85Neg.rda")          ##另加代码，保存data.r和data.p 85%比例药物和基因负相关rug <- c("Daporinad","Vinblastine","Vinblastine","Vinorelbine","Eg5","Dinaciclib","Paclitaxel","CDK9","Vincristine","Podophyllotoxin bromide",
#           "PD0325901","MK-1775","Trametinib","Dihydrorotenone","BMS-754807","AZD5153","Osimertinib","Afatinib","Wee1 Inhibito","Taselisib","AZD5438",
#           "Ulixertinib","Entinostat","YK-4-279","ULK1","AZD5582","VSP34","PAK","OTX015","Afuresertib","Erlotinib","ERK","SCH772984","AZD3759","Sorafenib",
#           "IWP-2","NU7441","JQ1","GSK343","VX-11e","Fulvestran","GSK269962A","Lapatinib","ZM447439","MK-2206","Crizotinib","RO-3306","Gefitinib","AZD6482",
#           "PRT062607","I-BET-762","BMS-345541","VE-822","Elephantin","Ipatasertib","Sinularin","Alpelisib","Tamoxifen","Nilotinib","Palbociclib",
#           "WIKI4","GSK2606414","Oxaliplatin","Linsitinib","TAF1","PF???4708671","Ribociclib","Sapitinib","LGK974","OF-1","VE821","AGI6780",
#           "Wnt-C59","Selumetinib","GSK1904529A","AZD5991","KRAS(G12C)Inhibitor???12","ML323","")

