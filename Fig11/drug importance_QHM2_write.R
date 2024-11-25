
setwd(paste('D:/LUAD/21.0.DrugPredict/1.drug1_GDSC/zong/4.drug_Importance',sep=''))

sample_name <- "TCGA"
#load("project_image.RData")

library(neuralnet)
library(dplyr)
library(Boruta)
library(survival)
library(randomForestSRC)
library(survival)
library(e1071)
library(xgboost)
library(Matrix)
# library(xgboostExplainer)
library(timeROC)
library(survminer)
library(ggplot2)


Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

cancer <- c("TCGA")
sur <- read.csv("D:/LUAD/21.0.DrugPredict/1.drug1_GDSC/zong/4.drug_Importance/sur_TCGA.csv",header = T)
colnames(sur) <- c("SampleName","Time","Status")
sur$Time <-as.numeric(sur$Time)
sur$Status <-as.numeric(sur$Status)
all_sur_data <- sur



#########  EXP                    ##第一个数据库GDSC2
library(data.table)   
dat <- fread(paste("D:/LUAD/21.0.DrugPredict/1.drug1_GDSC/out_put_TCGA_DEGs_zong.csv",sep=''))
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat<- dat[,-1]
                                                          # 第1个数据库GDSC2比例60%，3Drug
load("D:/LUAD/21.0.DrugPredict/1.drug1_GDSC/zong/1.cor_gene/data.r_p60.rda")
drug <- colnames(data.r)
dat <- dat[,match(drug,colnames(dat))]
dat_gene1 <- dat



                                  ##第2个数据库CTRP2
library(data.table)   
dat <- fread(paste("D:/LUAD/21.0.DrugPredict/2.drug2_CTRP/out_put_TCGA_DEGs_zong.csv",sep=''))
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat<- dat[,-1]
                                 # 第2个数据库CTRP2比例80%，20Drug
load("D:/LUAD/21.0.DrugPredict/2.drug2_CTRP/zong/11.cor_gene_/data.r_p80.rda") 
drug <- colnames(data.r)
dat <- dat[,match(drug,colnames(dat))]
dat_gene2 <- dat

# IC50drug <- c("Trametinib","SCH772984","Selumetinib")
# dat_gene2 <- dat_gene2[,match(IC50drug,colnames(dat_gene2))]

all_name <- names(which(table(c(rownames(dat_gene1),rownames(dat_gene2) ))==2))


dat_gene1 <- dat_gene1[match(all_name,rownames(dat_gene1)),]
dat_gene2 <- dat_gene2[match(all_name,rownames(dat_gene2)),]

my_data <- data.frame(SampleName=rownames(dat_gene2),dat_gene1,dat_gene2)
# write.csv(my_data, file = "23DrugExp.csv", row.names = FALSE) #输出23个药物在样本中的表达量文件


gene_exp = my_data %>%
  dplyr::select(SampleName, everything())




#########  T
# gene_exp = read.csv('~/diff/data/T.csv',header=T,check.names=F)
# gene_exp = gene_exp[grep(sample_name,gene_exp$CODE),]
# if(length(rownames(my_data))==0)
# {
#   next
# }
# for(i in unique(gene_exp$T))
# {
#   eval(parse(text=paste0('gene_exp$',i,'=rep(0,',length(rownames(gene_exp)),')')))
#   eval(parse(text=paste0('gene_exp$',i,'[which(gene_exp$T=="',i,'")]','=1' )))
# }
# gene_exp = gene_exp %>% 
#   subset(select = -c(CODE,T)) %>%
#   select(SampleName, everything())




mixed = merge(gene_exp, all_sur_data, by = "SampleName")
mixed = na.omit(mixed)

sample = mixed$SampleName
mixed= subset(mixed, select = -c(SampleName))
rownames(mixed) = sample





##############ANN
set.seed(969659)
data = na.omit(as.data.frame(lapply(mixed,as.numeric)))
data$Time = data$Time
nn = neuralnet(Status+Time~., data=data, hidden=c(5),stepmax=1e8,learningrate=0.1,startweights = "random",linear.output = F, lifesign = "full")

# temp = predict(nn,data)
# temp_median = median(temp[,1])
# temp[,1] = ifelse((temp[,1]) > temp_median, 1, 0)
# rownames(temp) = NULL
# print(length(which(temp[,1]==data$Status)))
# if(length(nn)<11) {nn <- neuralnet(Status+Time ~ ., data=data, hidden=15)}

weight <- data.frame(sum = rowSums(abs(nn$weights[[1]][[1]]))[-1],row.names = colnames(gene_exp)[-1])
ANN_order <- data.frame(ANN=Fun(weight[order(weight$sum),]),gene=rownames(weight)[order(weight$sum)])
# pdf(paste0(name_fold,"/ANN_struct .pdf"),width=20,height=20)
# p =invisible( plot(nn))
#  print(p)
# dev.off()

pdf("ANN_rank.pdf",width=15)
p <- ggplot(weight, aes(x = reorder(rownames(weight), -sum), y = sum)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = 'darkblue') +
  labs(x = "FeatureName", y = "weight") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(p)
dev.off()



#############Boruta
boruta_result <- Boruta(mixed[,-c(ncol(mixed),ncol(mixed)-1)], Surv(time = mixed$Time, event = mixed$Status))
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










##############RF
p.obj <- rfsrc(Surv(Time,Status)~.,data = mixed,
               mtry=3,
               nodesize=5,
               ntree=2000,
               tree.err = TRUE,
               importance = TRUE
)
out.rf <- var.select(object=p.obj,conservative = "high",)
RF_order <- data.frame(RF=Fun(out.rf$varselect$vimp),gene=rownames(out.rf$varselect))
full <- predict(p.obj,mixed ,importance = TRUE)
pdf('RF_var_important.pdf',height=length(colnames(mixed))/5+2)
plot(full)
dev.off()

var_out<- data.frame(vimp=out.rf$varselect$vimp,group=rownames(out.rf$varselect))

pdf('vimp.pdf',width=13,height=5)

library(ggplot2)
library(hrbrthemes)
library(showtext)
showtext_auto()
ggplot(var_out, aes(x = reorder(group, -vimp), y = vimp)) + 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .7,fill='darkred') + aes(fill=vimp)+
  xlab("Gene") + 
  ylab("Vimp")+  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

dev.off()



##############SVM
data = as.data.frame(apply(mixed[,-c(ncol(mixed),ncol(mixed)-1)],2,function(x) as.numeric(as.character(x))))
data = cbind(status=mixed$Status,data)

source("D:/LUAD/21.0.DrugPredict/1.drug1_GDSC/zong/4.drug_Importance/4.drug_Importance/msvmRFE.R")


nfold = 10
nrows = nrow(data)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))

results = lapply(folds, svmRFE.wrap, data, k=10, halve.above=100)
top.features = WriteFeatures(results, data, save=F)
SVM_order <- data.frame(SVM=Fun(length(top.features$AvgRank)-top.features$AvgRank),gene=top.features$FeatureName)
library(ggplot2)

pdf("rank.pdf",width=10)
p <- ggplot(top.features, aes(x = reorder(FeatureName, AvgRank), y = AvgRank)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = 'lightblue') +
  labs(x = "FeatureName", y = "rank") + theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()
# featsweep = lapply(1:5, FeatSweep.wrap, results, data)
# 
# errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
# pdf(paste0(name_fold,"/svm_rfe.pdf"), height = 8, width = 10)
# plot(errors,type=c("b"),col = 4, lty = 2)
# text(x = errors,y = sprintf("%.5f", errors),label  = sprintf("%.5f", errors), pos = 3, offset = 0.5, col = "red")
# dev.off()







##################XGboost
dtrain <- xgb.DMatrix(as.matrix(mixed[,-c(ncol(mixed),ncol(mixed)-1)]), label = mixed$Time,weight = mixed$Status)

params <- list(
  objective = "survival:cox",
  eval_metric = "cox-nloglik",
  eta = 0.01,
  max_depth = 3,
  subsample = 0.8,
  colsample_bytree = 0.8
)
xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
feature_importance <- xgb.importance(model = xgb_model)
XGBOOST_order <- data.frame(XGBOOST=Fun(feature_importance$Gain),gene=feature_importance$Feature)

out_gene <- colnames(mixed)[1:(ncol(mixed)-2)][which(!(colnames(mixed)[1:(ncol(mixed)-2)]%in%XGBOOST_order$gene))]
XGBOOST_order <- rbind(XGBOOST_order,data.frame(gene=out_gene,XGBOOST=rep(0,length(out_gene))))

# predict_full <- predict(xgb_model,dtrain)
# full_time <- mixed$Time/365
# ROC_rt=timeROC(T=full_time,delta=mixed$Status,
#                marker=predict_full,cause=1,
#                weighting='aalen',
#                times=c(1,2,3),ROC=TRUE)
# pdf(paste0(name_fold,'/xgboost_roc_full.pdf'),width=5,height=5)
# plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
# plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
# plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
# legend('bottomright',
#        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
#          paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
#          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
#        col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
# dev.off()

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


save.image(file="project_image.RData")
#load('project_image.RData')


out <- ANN_order %>% left_join(BORUTA_order,by='gene') %>% left_join(RF_order,by='gene') %>% left_join(XGBOOST_order,by='gene') %>% left_join(SVM_order,by='gene') 
rownames(out) <- out$gene
out <- subset(out,select=-c(gene))
out[is.na(out)] <- 0
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


widelength <- length(rownames(XGBOOST_order))
color <- c("#4b5cc4","#dc3023","#057748","#fff143","#758a99","#177cb0")

my_data <- data.frame(
  Category = rep(c(colnames(out)[1], colnames(out)[2], colnames(out)[3], colnames(out)[4],colnames(out)[5],colnames(out)[6]), each = length(out[,1])),
  Sample = rep(1:length(rownames(out)), times = length(colnames(out))),
  Value = value,
  Color = color,
  name = rownames(out)
)

my_data$Category = factor(my_data$Category,levels = colnames(out))

color[1] = "#4c221b"

pdf('algorithm.pdf',width=widelength/3+5,height=10)
ggplot(my_data, aes(x = as.factor(Sample), y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = color) +
  labs(title = "", x = "Drug", y = "Important") +
  scale_x_discrete(labels = my_data$name) +  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,size=14,color='black', hjust = 1,family = "serif"),
        axis.text.y = element_text(size=14,color='black',family = "serif"))
dev.off()


pdf('algorithm2.pdf',width=widelength/3+5,height=70)
ggplot(my_data, aes(x = as.factor(Sample), y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = color) +
  labs(title = "", x = "Drug", y = "Important") +
  scale_x_discrete(labels = my_data$name) +  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,size=14,color='black', hjust = 1),axis.text.y = element_text(size=14,color='black'))
dev.off()


write.csv(out,file=paste('out.csv',sep=''),row.names = T)


