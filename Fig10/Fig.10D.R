
setwd(paste('D:/LUAD/15.TIDE/11.GeneRank_TIDE',sep=''))
library(data.table)
library(GSVA)
library(randomForestSRC)
library(survival)
library(neuralnet)
library(Boruta)
library(xgboost)
library(Matrix)
library(xgboostExplainer)
library(ggplot2)
library(dplyr)
library(e1071)

#ICD对TLSscore的基因排序，所以读取ICD基因集的表达量；下面31行代码读TLS的得分

#基因对TIDE打分的重要性排序，以下，读取TLSICD基因的表达量。
data_BER <-  fread(paste("D:/LUAD/15.TIDE/11.GeneRank_TIDE/Exp_TCGA_TLSICD.csv",sep=''))
data_BER <- as.data.frame(data_BER)
# n1 <- grep("Tumor",data_BER$Group)
# dat_BER <- data_BER[n1,]
# gene_exp <- dat_BER[,-c(1,2)]
# colnames(gene_exp)[length(colnames(gene_exp))] <- gsub("[\r\n]", "",colnames(gene_exp)[length(colnames(gene_exp))])
# gene_exp$SampleName <- gsub("[\r\n]", "",gene_exp$SampleName)
# drug_exp <- gene_exp
drug_exp <- data_BER

   #读取TLSscore；和19行代码对应
   #读取TIDE打分；和19行代码对应
data_BER <-  fread(paste("D:/LUAD/15.TIDE/11.GeneRank_TIDE/LUAD_TIDE_TCGA.csv",sep=''),header = FALSE) 
data_BER <- as.data.frame(data_BER)
# data_BER <- t(data_BER)  #输入文件TIDE打分本身并未倒置

# 将第一行作为行名 #AI代码
colnames(data_BER) <- data_BER[1,]
# 删除原始数据的第一行（现在已经成为行名）
data_BER <- data_BER[-1, ]

# temp <- data_BER 
temp <- as.data.frame(data_BER)
colnames(temp)[2] <- "risk" #需改数字，temp的第2列
temp$risk <-as.numeric(temp$risk)  #评分


Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

##合并时要加SampleName
temp$SampleName
drug_exp$SampleName
colnames(drug_exp)

train_geneset_data <- merge(drug_exp,temp,by='SampleName')[,-c(1)]


p.obj <- rfsrc(risk ~.,data = train_geneset_data,
               mtry=3,
               nodesize=5,
               ntree=2000,
               tree.err = TRUE,
               importance = TRUE
)
out.rf <- var.select(object=p.obj,conservative = "high",)
RF_order <- data.frame(RF=Fun(out.rf$varselect$vimp),gene=rownames(out.rf$varselect))
full <- predict(p.obj,train_geneset_data ,importance = TRUE)
pdf('train_var_important.pdf',height=length(colnames(train_geneset_data))/3)
plot(full)
dev.off()

var_out<- data.frame(vimp=out.rf$varselect$vimp,group=rownames(out.rf$varselect))
pdf('randomforest_vimp.pdf',width=13,height=5)

library(ggplot2)
library(hrbrthemes)
library(showtext)
showtext_auto()
p <- ggplot(var_out, aes(x = reorder(group, -vimp), y = vimp)) + 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .7,fill='darkred') + aes(fill=vimp)+
  xlab("Gene") + 
  ylab("Vimp")+  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(p)
dev.off()




train_geneset_data0 <- train_geneset_data
#load("D:/R/daima/药物对基因集得分重要性/train_geneset_data.Rda")

train_geneset_data <- train_geneset_data0
colnames(train_geneset_data) <- gsub("[^A-Za-z0-9_]", "_", colnames(train_geneset_data))

#as.numeric(train_geneset_data)
train_geneset_data[] <- lapply(train_geneset_data, as.numeric)
nn <- neuralnet(risk ~ ., data=train_geneset_data, hidden=c(100,100),stepmax=1e6,learningrate=0.1,startweights = "random",linear.output = T, lifesign = "full")
pdf(paste0("ANN_struct.pdf"),width=40,height=40)
p <- plot(nn,col.text = "black",col.entry='black',col.intercept='#50026E',col.out.synapse='#008209',col.hidden.synapse='#0F4FA8',col.entry.synapse='#A40004',col.out='#008209',fontsize=30, col.hidden = "#0F4FA8",radius=0.1, show.weights = FALSE,information=F)
print(p)
dev.off()
weight <- data.frame(sum = rowSums(abs(nn$weights[[1]][[1]]))[-1],row.names = colnames(train_geneset_data)[-ncol(train_geneset_data)])
rownames(weight) <- colnames(train_geneset_data0)[-length(train_geneset_data0)]
ANN_order <- data.frame(ANN=Fun(weight[order(weight$sum),]),gene=rownames(weight)[order(weight$sum)])

library(ggplot2)
pdf("ann_rank.pdf",width=15)
p <- ggplot(weight, aes(x = reorder(rownames(weight), -sum), y = sum)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = 'darkblue') +
  labs(x = "FeatureName", y = "weight") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
#print(p)
dev.off()





boruta_result <- Boruta(train_geneset_data[,c(-ncol(train_geneset_data))], train_geneset_data[,ncol(train_geneset_data)])
Boruta_imp <- attStats(boruta_result)[order(attStats(boruta_result)$meanImp),]
BORUTA_order <- data.frame(BORUTA=Fun(Boruta_imp$meanImp),gene=rownames(Boruta_imp))
pdf('boruta_importance.pdf',width = 20,height=10)
par(oma=c(3,3,3,3)) 
plot(boruta_result,las=2,xlab='')
legend(x = 'topleft', 
       legend = c(paste('P-value:',boruta_result$pValue),sep=''),
       lty = 0,
       bty = 'n')
dev.off()

pdf('boruta_history.pdf',width = 24,height=10)
par(oma=c(3,3,3,3)) 
plot(plotImpHistory(boruta_result),las=2)
legend(x = 'topleft', 
       legend = c(paste('P-value:',boruta_result$pValue),sep=''),
       lty = 0,
       bty = 'n')
dev.off()




dtrain <- xgb.DMatrix(as.matrix(train_geneset_data[,c(-ncol(train_geneset_data))]), label = train_geneset_data[,ncol(train_geneset_data)])

params <- list(
  eta = 0.01,
  max_depth = 3,
  subsample = 0.8,
  colsample_bytree = 0.8
)
xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
feature_importance <- xgb.importance(model = xgb_model)
XGBOOST_order <- data.frame(XGBOOST=Fun(feature_importance$Gain),gene=feature_importance$Feature)

colnames(train_geneset_data)[-length(colnames(train_geneset_data))]

out_gene <- colnames(train_geneset_data)[-length(colnames(train_geneset_data))][which(!(colnames(train_geneset_data)[-length(colnames(train_geneset_data))]%in%XGBOOST_order$gene))]
XGBOOST_order <- rbind(XGBOOST_order,data.frame(gene=out_gene,XGBOOST=rep(0,length(out_gene))))

pdf("XGboost_rank.pdf",width=13)
p <- ggplot(feature_importance, aes(x = reorder(feature_importance$Feature, -feature_importance$Gain), y = feature_importance$Gain)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = 'orange2') +
  labs(x = "FeatureName", y = "Gain") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
#print(p)
dev.off()




data = cbind(risk=train_geneset_data$risk,train_geneset_data[,-ncol(train_geneset_data)])
data$risk = ifelse(data$risk >  median(data$risk), 1, 0)
source("D:/LUAD/15.TIDE/11.GeneRank_TIDE/msvmRFE.R")
nfold = 10
nrows = nrow(data)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))

results = lapply(folds, svmRFE.wrap, data, k=10, halve.above=100)
top.features = WriteFeatures(results, data, save=F)
SVM_order <- data.frame(SVM=Fun(length(top.features$AvgRank)-top.features$AvgRank),gene=top.features$FeatureName)

library(ggplot2)
pdf("svm_rank.pdf",width=13)
p <- ggplot(top.features, aes(x = reorder(FeatureName, AvgRank), y = AvgRank)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = 'lightblue') +
  labs(x = "FeatureName", y = "rank") + theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
#print(p)
dev.off()


save.image(file="project_image.RData")

save(ANN_order,file='ANN_order.Rdata')
save(BORUTA_order,file='BORUTA_order.Rdata')
save(RF_order,file='RF_order.Rdata')
save(XGBOOST_order,file='XGBOOST_order.Rdata')
save(SVM_order,file='SVM_order.Rdata')

load('BORUTA_order.Rdata')
load('XGBOOST_order.Rdata')
load('ANN_order.Rdata')
load('SVM_order.Rdata')
load('RF_order.Rdata')

out <- ANN_order %>% left_join(BORUTA_order,by='gene') %>% left_join(RF_order,by='gene') %>% left_join(XGBOOST_order,by='gene') %>% left_join(SVM_order,by='gene') 
rownames(out) <- out$gene
out <- subset(out,select=-c(gene))
out <- cbind(Mean = rowSums(out)/length(colnames(out)),out)
out <- out[order(-out$Mean),]*3
value <- c()

for(i in 1:length(colnames(out)))
{
  eval(parse(text=paste('value<-append(value,',out[,i],')',sep='')))
}

# color <- c()
# COLORS <- c("#70f3ff","#44cef6","#3eede7","#1685a9","#177cb0","#065279","#003472","#4b5cc4","#a1afc9","#2e4e7e","#3b2e7e","#4a4266","#426666","#425066","#574266","#8d4bbb","#815463","#815476","#4c221b","#003371","#56004f","#801dae","#4c8dae","#b0a4e3","#cca4e3","#edd1d8","#e4c6d0","#ff461f","#ff2d51","#f36838","#ed5736","#ff4777","#f00056","#ffb3a7","#f47983","#db5a6b","#c93756","#f9906f","#f05654","#ff2121","#f20c00","#8c4356","#c83c23","#9d2933","#ff4c00","#ff4e20","#f35336","#dc3023","#ff3300","#cb3a56","#a98175","#b36d61","#ef7a82","#ff0097","#c32136","#be002f","#c91f37","#bf242a","#c3272b","#9d2933","#60281e","#622a1d","#bce672","#c9dd22","#bddd22","#afdd22","#a3d900","#9ed900","#9ed048","#96ce54","#00bc12","#0eb83a","#0eb83a","#0aa344","#16a951","#21a675","#057748","#0c8918","#00e500","#40de5a","#00e079","#00e09e","#3de1ad","#2add9c","#2edfa3","#7fecad","#a4e2c6","#7bcfa6","#1bd1a5","#48c0a3","#549688","#789262","#758a99","#50616d","#424c50","#41555d","#eaff56","#fff143","#faff72","#ffa631","#ffa400","#fa8c35","#ff8c31","#ff8936","#ff7500","#ffb61e","#ffc773","#ffc64b","#f2be45","#f0c239","#e9bb1d","#d9b611","#eacd76","#eedeb0","#d3b17d","#e29c45","#a78e44","#c89b40","#ae7000","#ca6924","#b25d25","#b35c44","#9b4400","#9c5333","#a88462","#896c39","#827100","#6e511e","#7c4b00","#955539","#845a33","#ffffff","#e9e7ef")
# 
# 
# for(j in 1:length(colnames(out)))
# {
#   eval(parse(text=paste('color <- append(color,c(',colnames(out)[j],'="',sample(COLORS,size=1),'"))',sep='')))
# }

#原代码
widelength <- length(rownames(XGBOOST_order))
color <- c("#4b5cc4","#dc3023","#057748","#fff143","#758a99","#177cb0") ##？
# widelength <- length(rownames(XGBOOST_order))
# color <- c("#4b5cc4","#dc3023","#057748","#fff143","#177cb0","#758a99") 


my_data <- data.frame(
  Category = rep(c(colnames(out)[1], colnames(out)[2], colnames(out)[3], colnames(out)[4],colnames(out)[5],colnames(out)[6]), each = length(out[,1])),
  Sample = rep(1:(length(colnames(train_geneset_data))-1), times = length(colnames(out))),
  Value = value,
  Color = color,
  name = rownames(out)
)

my_data$Category = factor(my_data$Category,levels = colnames(out))

color[1] = "#4c221b"



# 定义新罗马字体族
my_font <- "Times"

#   theme(
#     axis.text.x = element_text(angle = 45, size = 14, color = 'black', hjust = 1, family = my_font),  # 设置x轴标签字体为新罗马
#     axis.text.y = element_text(size = 14, color = 'black', family = my_font),  # 设置y轴标签字体为新罗马
#     axis.title = element_text(size = 14, color = 'black', family = my_font),  # 设置轴标题字体为新罗马
#     plot.title = element_text(size = 14, color = 'black', family = my_font),  # 设置图标题字体为新罗马
#     legend.text = element_text(size = 14, color = 'black', family = my_font),  # 设置图例文本字体为新罗马
#     legend.title = element_text(size = 14, color = 'black', family = my_font)  # 设置图例标题字体为新罗马
#   )
# dev.off()

pdf('TLSICDgene_TIDE-algorithm.pdf',width=widelength/3+5,height=10)
ggplot(my_data, aes(x = as.factor(Sample), y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = color) +
  labs(title = "", x = "Gene", y = "Important") +
  scale_x_discrete(labels = my_data$name) +  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,size=14,color='black', hjust = 1, family = my_font),
        axis.text.y = element_text(size=14,color='black', family = my_font))
dev.off()

write.csv(out,file=paste('out.csv',sep=''),row.names = T)


