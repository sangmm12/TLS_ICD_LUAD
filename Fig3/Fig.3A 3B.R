
# install.packages("devtools")
# library(devtools)
# devtools::install_github("AppliedDataSciencePartners/xgboostExplainer")

##设置目录
setwd('~/projects/老课题MachineLearning_单独去批次_TLS/TCGA不去批次跑模型')

rm(list=ls())
##载入R包
library(dplyr)
library(tibble)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(neuralnet)
library(xgboost)
library(Matrix)
library(xgboostExplainer)
library(ggplot2)
library(ggrepel)
library(parallel) 
library(lightgbm)
library(kernlab)

Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

##因为有些GSE数据库问题，cox的最大深度以及SPC的fold设置为3
coxphitermax = 1
SPCnfold = 3
seed=123456
#GSE37745死亡样本太少？
final_result <- data.frame(tcga=c(1),tcga11=c(1))#GSE72094=c(1))#,GSE19188=c(1),GSE30219=c(1),GSE31210=c(1),GSE37745=c(1),GSE42127=c(1),GSE50081=c(1),GSE68465=c(1),GSE72094=c(1))
list_data = list()
#c('tcga','GSE12417','GSE37642','GSE146173','GSE106291')
#c('GSE37642','tcga','GSE12417','GSE146173','GSE106291')
#c('GSE106291','tcga','GSE12417','GSE37642','GSE146173')

test_list = c('tcga','tcga11')#,'GSE19188','GSE30219','GSE31210','GSE37745','GSE42127','GSE50081','GSE68465','GSE72094')

colnames(final_result) = test_list

for(sample_name in test_list)
{
  all_sur_data <- read.csv(paste('~/projects/老课题MachineLearning_单独去批次_TLS/TCGA不去批次跑模型/sur_',sample_name,'.csv',sep=''),header=T)
  gene_exp <- t(read.csv(paste('~/projects/老课题MachineLearning_单独去批次_TLS/TCGA不去批次跑模型/exp',sample_name,'.csv',sep=''),check.names = F)[,-1]) %>%
    as.data.frame()
  colnames(gene_exp) = gene_exp[1,]
  gene_exp = gene_exp[-1,]
  gene_exp$SampleName = rownames(gene_exp)
  #if(max(sapply(gene_exp[,-c(1)],max))<10){gene_exp[,-c(1)] = 2^gene_exp[,-c(1)]}
  # temp_sm = gene_exp[,c(1)]
  # gene_exp = log2(gene_exp[,-c(1)]+1000)
  # gene_exp = Fun(gene_exp)
  # gene_exp$SampleName = temp_sm
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','OS.time',colnames(mixed))
  colnames(mixed) <- sub('Status','OS',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  
  mixed[is.na(mixed)] = 0
  mixed = mixed[,which(!colnames(mixed)%in%names(which(colSums(mixed)==0)))]
  mixed$SampleName <- temp_row
  mixed <- dplyr::select(mixed, SampleName,OS, OS.time, everything())
  mixed = mixed[,which(colnames(mixed)!='KIAA0101')]
  #mixed[is.na(mixed)] <- rnorm(length(mixed[is.na(mixed)]))
  eval(parse(text=paste0('list_data$',sample_name,'=mixed')))
}

eval(parse(text=paste0('mixed = list_data$',test_list[1])))

cl <- makeCluster(40)  #最大线程数

clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax"))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(tibble))
clusterEvalQ(cl, library(survival))
clusterEvalQ(cl, library(randomForestSRC))
clusterEvalQ(cl, library(glmnet))
clusterEvalQ(cl, library(plsRcox))
clusterEvalQ(cl, library(superpc))
clusterEvalQ(cl, library(gbm))
clusterEvalQ(cl, library(CoxBoost))
clusterEvalQ(cl, library(survivalsvm))
clusterEvalQ(cl, library(BART))
clusterEvalQ(cl, library(neuralnet))
clusterEvalQ(cl, library(xgboost))
clusterEvalQ(cl, library(Matrix))
clusterEvalQ(cl, library(xgboostExplainer))
clusterEvalQ(cl, library(ggplot2))
clusterEvalQ(cl, library(ggrepel))
clusterEvalQ(cl, library(lightgbm))
clusterEvalQ(cl, library(kernlab))


save(file='listdata.Rda',list_data)
#70% 30%划分
# train_idx <- sample(1:nrow(mixed), 0.7 * nrow(mixed)) # 70% 的数据作为训练集
# mixed_test <- mixed[-train_idx, ]
# mixed <- mixed[train_idx, ]

##所有预测变量
#RS_COXBOOST
#RS_STEPWISE
#RS_RIDGE
#RS_LASSO
#RS_SVM
#RS_GBDT
#RS_SPC
#RS_PLS
#RS_ANN
#RS_XGBOOST

seed <- 123456
predict_list = c("RS_COXBOOST","RS_STEPWISE_f","RS_RIDGE","RS_LASSO","RS_SVM","RS_GBDT","RS_SPC","RS_PLS","RS_ANN_05","RS_XGBOOST_01","RS_LGBM","RS_NUSVC","RS_Enet_01","RS_Enet_02","RS_Enet_03","RS_Enet_04","RS_Enet_05","RS_Enet_06","RS_Enet_07","RS_Enet_08","RS_Enet_09","RS_STEPWISE_b","RS_STEPWISE_bo","RS_XGBOOST_02","RS_XGBOOST_03","RS_XGBOOST_04","RS_XGBOOST_05","RS_XGBOOST_06","RS_XGBOOST_07","RS_XGBOOST_08","RS_XGBOOST_09","RS_XGBOOST_10","RS_ANN_06","RS_ANN_07","RS_ANN_08","RS_ANN_09","RS_ANN_10","RS_ANN_11","RS_ANN_12","RS_ANN_13","RS_ANN_14","RS_ANN_15")
predict_list_name = c("CoxBoost","Stepwise Cox [forward]","Ridge","Lasso","Survival SVM","GBDT","Supervised PCA","plsRcox","ANN [hidden=5]","XGboost [max_depth=1]","LightGBM","Nu-SVC","Enet [a=0.1]","Enet [a=0.2]","Enet [a=0.3]","Enet [a=0.4]","Enet [a=0.5]","Enet [a=0.6]","Enet [a=0.7]","Enet [a=0.8]","Enet [a=0.9]","Stepwise Cox [backward]","Stepwise Cox [both]","XGboost [max_depth=2]","XGboost [max_depth=3]","XGboost [max_depth=4]","XGboost [max_depth=5]","XGboost [max_depth=6]","XGboost [max_depth=7]","XGboost [max_depth=8]","XGboost [max_depth=9]","XGboost [max_depth=10]","ANN [hidden=6]","ANN [hidden=7]","ANN [hidden=8]","ANN [hidden=9]","ANN [hidden=10]","ANN [hidden=11]","ANN [hidden=12]","ANN [hidden=13]","ANN [hidden=14]","ANN [hidden=15]")
comb <- combn(predict_list,2)
temp_remove = c()
for(i in (1:(length(comb[,])/2)))
{
  if(as.character(gsub('[0-9]','',comb[,i][1]))==as.character(gsub('[0-9]','',comb[,i][2])))
  {
    temp_remove = append(i,temp_remove)
  }
}
comb = comb[,-c(temp_remove)]


### CoxBoost ####
predict_comb = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  pen <- optimCoxBoostPenalty(temp_mixed$OS.time,temp_mixed$OS,as.matrix(temp_mixed[,-c(1,2,3)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(temp_mixed$OS.time,temp_mixed$OS,as.matrix(temp_mixed[,-c(1,2,3)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit <- CoxBoost(temp_mixed$OS.time,temp_mixed$OS,as.matrix(temp_mixed[,-c(1,2,3)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  
  as.numeric(predict(fit,newdata=x[,-c(1,2,3)],newtime=x[,3], newstatus=x[,2], type="lp"))
}
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  pen <- optimCoxBoostPenalty(temp_mixed$OS.time,temp_mixed$OS,as.matrix(temp_mixed[,-c(1,2,3)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(temp_mixed$OS.time,temp_mixed$OS,as.matrix(temp_mixed[,-c(1,2,3)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit <- CoxBoost(temp_mixed$OS.time,temp_mixed$OS,as.matrix(temp_mixed[,-c(1,2,3)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  
  rs <- cbind(x[,2:3],RS = as.numeric(predict(fit,newdata=x[,-c(1,2,3)],newtime=x[,3], newstatus=x[,2], type="lp")))
  cc <- summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[[1]]
}
temp_vector = c()
RS_COXBOOST = parLapply(cl,list_data,predict_comb)
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],'CoxBoost')
rm(list=c('fit','predict_comb','temp_vector','temp_comb'))




#### Stepwise Cox ####
for (direction in c("both", "backward", "forward")) {
  clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax","direction"))
  predict_comb = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit <- step(coxph(Surv(OS.time,OS)~.,temp_mixed[,-c(1)],iter.max=coxphitermax),direction = direction,steps = 1)
    as.numeric(predict(fit,type = 'risk',newdata = x[,-c(1)]))
  }
  if(direction == 'forward') {RS_STEPWISE_f = parLapply(cl,list_data,predict_comb)}
  if(direction == 'both') {RS_STEPWISE_bo = parLapply(cl,list_data,predict_comb)}
  if(direction == 'backward') {RS_STEPWISE_b = parLapply(cl,list_data,predict_comb)}
  predict_x = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit <- step(coxph(Surv(OS.time,OS)~.,temp_mixed[,-c(1)],iter.max=coxphitermax),direction = direction,steps = 1)
     
    rs <- cbind(x[,c(2,3)],RS = as.numeric(predict(fit,type = 'risk',newdata = x[,-c(1)])))
    cc <- summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
  }
  temp_vector = c()
  temp_comb = parLapply(cl,list_data,predict_x)
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Stepwise Cox','[',direction,']')))
  rm(list=c('fit','predict_comb','temp_vector','temp_comb'))
}


#### Lasso,Ridge,Enet ####
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax","alpha"))
  predict_comb = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit <- cv.glmnet(as.matrix(temp_mixed[,-c(1,2,3)]), as.matrix(Surv(temp_mixed$OS.time,temp_mixed$OS)),family = "cox",alpha=alpha,nfolds = 10)
    return(as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2,3)]),s=fit$lambda.min)))
  }
  predict_x = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit <- cv.glmnet(as.matrix(temp_mixed[,-c(1,2,3)]), as.matrix(Surv(temp_mixed$OS.time,temp_mixed$OS)),family = "cox",alpha=alpha,nfolds = 10)
    rs <- cbind(x[,2:3],RS = as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2,3)]),s=fit$lambda.min)))
    summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
  }
  temp_vector = c()
  temp_comb = parLapply(cl,list_data,predict_x)
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Enet ','[a=',alpha,']',sep='')))
  if(alpha == 0)
  {
    rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Ridge'))
    RS_RIDGE <- parLapply(cl,list_data,predict_comb)
  }
  if(alpha == 1) 
  {
    rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Lasso'))
    RS_LASSO <- parLapply(cl,list_data,predict_comb)
  }
  if(alpha!=0||alpha!=1)
  {
    eval(parse(text=paste0("RS_Enet_0",alpha*10," = parLapply(cl,list_data,predict_comb)")))
  }
  rm(list=c('fit','predict_comb','temp_vector','temp_comb'))
}


#### Survival SVM ####
predict_comb = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit <- survivalsvm(Surv(OS.time,OS)~., data= temp_mixed[,-c(1)], gamma.mu = 1)
  
  as.numeric(predict(fit, x[,-c(1)])$predicted)
}
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit <- survivalsvm(Surv(OS.time,OS)~., data= temp_mixed[,-c(1)], gamma.mu = 1)
  
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit, x[,-c(1)])$predicted))
  summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
temp_vector = c()
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
RS_SVM = parLapply(cl,list_data,predict_comb)
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Survival SVM'))
rm(list=c('fit','predict_comb','temp_vector','temp_comb'))


#### GBDT ####
predict_comb = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  
  as.numeric(predict(fit,x[,-c(1)],n.trees = best,type = 'link'))
}
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit,x[,-c(1)],n.trees = best,type = 'link')))
  cc <- summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
temp_comb = parLapply(cl,list_data,predict_x)
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
RS_GBDT = parLapply(cl,list_data,predict_comb)
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('GBDT'))
rm(list=c('fit','predict_comb','temp_vector','temp_comb'))


#### Supervised principal components ####
predict_comb = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  data <- list(x=t(temp_mixed[,-c(1,2,3)]),y=temp_mixed$OS.time,censoring.status=temp_mixed$OS,featurenames=colnames(temp_mixed)[-c(1,2,3)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  set.seed(seed)
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                       n.fold = SPCnfold,
                       n.components=3,
                       min.features=5,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  
  as.numeric(superpc.predict(fit,data,list(x=t(x[,-c(1,2,3)]),y=x$OS.time,censoring.status=x$OS,featurenames=colnames(x)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred)
}
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  data <- list(x=t(temp_mixed[,-c(1,2,3)]),y=temp_mixed$OS.time,censoring.status=temp_mixed$OS,featurenames=colnames(temp_mixed)[-c(1,2,3)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  set.seed(seed)
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                       n.fold = SPCnfold,
                       n.components=3,
                       min.features=5,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  
  rs <- cbind(x[,2:3],RS=as.numeric(superpc.predict(fit,data,list(x=t(x[,-c(1,2,3)]),y=x$OS.time,censoring.status=x$OS,featurenames=colnames(x)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred))
  cc <- summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
temp_comb = parLapply(cl,list_data,predict_x)
temp_comb = parLapply(cl,list_data,predict_x)
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
RS_SPC = parLapply(cl,list_data,predict_comb)
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Supervised PCA'))
rm(list=c('fit','predict_comb','temp_vector','temp_comb'))



#### plsRcox ####
pdf('plsRcox.pdf')
predict_comb = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  cv.plsRcox.res=cv.plsRcox(list(x=temp_mixed[,-c(1,2,3)],time=temp_mixed$OS.time,status=temp_mixed$OS),nf=10,verbose = FALSE)
  fit <- plsRcox(temp_mixed[,-c(1,2,3)],time=temp_mixed$OS.time,event=temp_mixed$OS,nt=as.numeric(cv.plsRcox.res[5]))
  
  as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2,3)]))
}
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  cv.plsRcox.res=cv.plsRcox(list(x=temp_mixed[,-c(1,2,3)],time=temp_mixed$OS.time,status=temp_mixed$OS),nf=10,verbose = FALSE)
  fit <- plsRcox(temp_mixed[,-c(1,2,3)],time=temp_mixed$OS.time,event=temp_mixed$OS,nt=as.numeric(cv.plsRcox.res[5]))
  
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2,3)])))
  summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
RS_PLS = parLapply(cl,list_data,predict_comb)
temp_vector = c()
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('plsRcox'))
rm(list=c('fit','predict_comb','temp_vector','temp_comb'))
dev.off()


#ANN
set.seed(seed)
for (max in seq(5,15,1)) 
{
  clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax","max"))
  predict_comb = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit <- neuralnet(OS ~ ., data=temp_mixed[,-c(1,3)], hidden=max,learningrate=0.2,startweights = "random",linear.output = T, lifesign = "full")
    n = 0
    
    while(length(fit)<14)
    {
      fit <- neuralnet(OS ~ ., data=temp_mixed[,-c(1,3)], hidden=max)
      n = n + 1
      if(n >= 3)
      {
        break
      }
    }
    if(n >= 3)
    {
      return(rep(0,length(x[,2])))
    }
    else
    {
      predict(fit,newdata = x[,-c(1,3)])
    }
  }
  predict_x = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit <- neuralnet(OS ~ ., data=temp_mixed[,-c(1,3)], hidden=max,learningrate=0.2,startweights = "random",linear.output = T, lifesign = "full")
    n = 0
    
    while(length(fit)<14)
    {
      fit <- neuralnet(OS ~ ., data=temp_mixed[,-c(1,3)], hidden=max)
      n = n + 1
      if(n >= 3)
      {
        break
      }
    }
    
    if(n >= 3)
    {
      return(0.5)
    }
    else
    {
      rs <- cbind(x[,2:3],RS=predict(fit,newdata = x[,-c(1,3)]))
      return(as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
    }
  }
  
  if(max<10)
  {
    eval(parse(text=paste0("RS_ANN_0",max," = parLapply(cl,list_data,predict_comb)")))
  }
  else
  {
    eval(parse(text=paste0("RS_ANN_",max," = parLapply(cl,list_data,predict_comb)")))
  }

  temp_vector = c()
  temp_comb = parLapply(cl,list_data,predict_x)
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('ANN [hidden=',max,']')))
  rm(list=c('fit','predict_comb','temp_vector','temp_comb'))
}


#XGboost
for (max in seq(1,10,1)) 
{
  clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax","max"))
  predict_comb = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    dtrain <- xgb.DMatrix(as.matrix(temp_mixed[,-c(1,2,3)]), label = temp_mixed$OS.time,weight = temp_mixed$OS)
    params <- list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.01,
      max_depth = max,
      subsample = 0.8,
      colsample_bytree = 0.8
    )
    xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
    
    dtest <- xgb.DMatrix(as.matrix(x[,-c(1,2,3)]), label = x$OS.time,weight = x$OS)
    predict(xgb_model,dtest)
  }
  predict_x = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    dtrain <- xgb.DMatrix(as.matrix(temp_mixed[,-c(1,2,3)]), label = temp_mixed$OS.time,weight = temp_mixed$OS)
    params <- list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.01,
      max_depth = max,
      subsample = 0.8,
      colsample_bytree = 0.8
    )
    xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
    
    dtest <- xgb.DMatrix(as.matrix(x[,-c(1,2,3)]), label = x$OS.time,weight = x$OS)
    rs <- cbind(x[,2:3],RS=predict(xgb_model,dtest))
    summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
  }
  if(max<10)
  {
    eval(parse(text=paste0("RS_XGBOOST_0",max," = parLapply(cl,list_data,predict_comb)")))
  }
  else
  {
    eval(parse(text=paste0("RS_XGBOOST_",max," = parLapply(cl,list_data,predict_comb)")))
  }
  temp_vector = c()
  temp_comb = parLapply(cl,list_data,predict_x)
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('XGboost [max_depth=',max,']')))
  rm(list=c('predict_comb','temp_comb'))
}


#LightGBM
predict_comb = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  train_data = lgb.Dataset(data = as.matrix(temp_mixed[,-c(1,2,3)]), label = temp_mixed$OS)
  params = list(objective = "binary", metric = "binary_logloss")
  fit = lgb.train(params = params, data = train_data, nrounds = 100)
  return(as.numeric(predict(fit,as.matrix(x[,-c(1,2,3)]))))
}
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  train_data = lgb.Dataset(data = as.matrix(temp_mixed[,-c(1,2,3)]), label = temp_mixed$OS)
  params = list(objective = "binary", metric = "binary_logloss")
  fit = lgb.train(params = params, data = train_data, nrounds = 100)
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit,as.matrix(x[,-c(1,2,3)]))))
  summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
RS_LGBM = parLapply(cl,list_data,predict_comb)
temp_vector = c()
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('LightGBM'))
rm(list=c('predict_comb','temp_vector','temp_comb','fit'))


#Nu-Support Vector Classification
predict_comb = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit <- ksvm(temp_mixed$OS ~ ., data = temp_mixed[, -c(1,2,3)], type = "nu-svc", kernel = "rbfdot")
  return(as.numeric(predict(fit, x[, -c(1,2,3)], type = "response")))
}
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit <- ksvm(temp_mixed$OS ~ ., data = temp_mixed[, -c(1,2,3)], type = "nu-svc", kernel = "rbfdot")
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit,x[,-c(1,2,3)], type = "response")))
  summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
RS_NUSVC = parLapply(cl,list_data,predict_comb)
temp_vector = c()
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Nu-SVC'))
rm(list=c('predict_comb','temp_vector','temp_comb','fit'))




for(i in 1:ncol(comb))
{
  temp_vector = c()
  for(data_name in names(list_data))
  {
    eval(parse(text = paste('rs <- cbind(list_data$',data_name,'[,2:3],RS1=',comb[1,i],'$',data_name,',RS2=',comb[2,i],'$',data_name,')',sep='')))
    temp_vector = append(temp_vector,summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1])
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0(predict_list_name[which(predict_list==comb[1,i])],' + ',predict_list_name[which(predict_list==comb[2,i])])))
}


#### Random survival forest ####
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],
                  ntree = 1000,nodesize = 5,
                  splitrule = 'logrank',
                  proximity = T,
                  forest = T,
                  seed = seed)
  
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted))
  summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
temp_vector = c()
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF')))


#ANN + Random survival forest
set.seed(seed)
for (max in seq(5,15,1)) {
  clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax","max"))
  predict_x = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],
                    ntree = 1000,nodesize = 5,
                    splitrule = 'logrank',
                    proximity = T,
                    forest = T,
                    seed = seed)
    
    n = 0
    fit <- neuralnet(OS ~ ., data=temp_mixed[,-c(1,3)], hidden=max)
    while(length(fit)<14)
    {
      fit <- neuralnet(OS ~ ., data=temp_mixed[,-c(1,3)], hidden=max)
      n = n + 1
      if(n >= 3)
      {
        break
      }
    }
    if(n >= 3)
    {
      return(0.5)
    }
    else
    {
      rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=predict(fit,newdata = x[,-c(1,3)]))
      return(summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1])
    }
  }
  temp_vector = c()
  temp_comb = parLapply(cl,list_data,predict_x)
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste('Survival RF + ANN [hidden=',max,']',sep='')))
  rm(list=c('fit'))
}


#GBDT + Random survival forest

predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],
                  ntree = 1000,nodesize = 5,
                  splitrule = 'logrank',
                  proximity = T,
                  forest = T,
                  seed = seed)
  fit_GB <- gbm(formula = Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],distribution = 'coxph',
                n.trees = 10000,
                interaction.depth = 3,
                n.minobsinnode = 10,
                shrinkage = 0.001,
                cv.folds = 10,n.cores = 6)
  best <- which.min(fit_GB$cv.error)
  set.seed(seed)
  fit_GB <- gbm(formula = Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],distribution = 'coxph',
                n.trees = best,
                interaction.depth = 3,
                n.minobsinnode = 10,
                shrinkage = 0.001,
                cv.folds = 10,n.cores = 8)
  
  rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=as.numeric(predict(fit_GB,x[,-c(1)],n.trees = best,type = 'link')))
  summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
}
temp_vector = c()
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Survival RF + GBDT'))
rm(list=c('fit_GB'))





#Stepwise Cox + Random survival forest
for (direction in c("both", "backward", "forward")) {
  clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax","direction"))
  predict_x = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],
                    ntree = 1000,nodesize = 5,
                    splitrule = 'logrank',
                    proximity = T,
                    forest = T,
                    seed = seed)
    fit <- step(coxph(Surv(OS.time,OS)~.,temp_mixed[,-c(1)],iter.max=coxphitermax),direction = direction)
    
    rs = cbind(x[,c(2,3)],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=as.numeric(predict(fit,type = 'risk',newdata = x[,-c(1)])))
    summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
  }
  temp_vector = c()
  temp_comb = parLapply(cl,list_data,predict_x)
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Stepwise Cox',' [',direction,']')))
  rm(list=c('fit'))
}



#Supervised principal components + Random survival forest
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],
                  ntree = 1000,nodesize = 5,
                  splitrule = 'logrank',
                  proximity = T,
                  forest = T,
                  seed = seed)
  
  data <- list(x=t(temp_mixed[,-c(1,2,3)]),y=temp_mixed$OS.time,censoring.status=temp_mixed$OS,featurenames=colnames(temp_mixed)[-c(1,2,3)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  set.seed(seed)
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                       n.fold = SPCnfold,
                       n.components=3,
                       min.features=5,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=as.numeric(superpc.predict(fit,data,list(x=t(x[,-c(1,2,3)]),y=x$OS.time,censoring.status=x$OS,featurenames=colnames(x)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred))
  summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
}
temp_vector = c()
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Supervised PCA')))
rm(list=c('fit'))




#plsRcox + Random survival forest

pdf('pls.pdf')
predict_x = function(x)
{
  set.seed(seed)
  temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
  fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],
                  ntree = 1000,nodesize = 5,
                  splitrule = 'logrank',
                  proximity = T,
                  forest = T,
                  seed = seed)
  
  set.seed(seed)
  cv.plsRcox.res=cv.plsRcox(list(x=temp_mixed[,-c(1,2,3)],time=temp_mixed$OS.time,status=temp_mixed$OS),nt=10,verbose = FALSE)
  fit <- plsRcox(temp_mixed[,-c(1,2,3)],time=temp_mixed$OS.time,event=temp_mixed$OS,nt=as.numeric(cv.plsRcox.res[5]))
  
  rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2,3)])))
  summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
}
temp_vector = c()
temp_comb = parLapply(cl,list_data,predict_x)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + plsRcox')))
rm(list=c('fit'))
dev.off()




#XGboost + Random survival forest
for (max in seq(1,5,1)) {
  clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax","max"))
  predict_x = function(x)
  {
    
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],
                    ntree = 1000,nodesize = 5,
                    splitrule = 'logrank',
                    proximity = T,
                    forest = T,
                    seed = seed)
    
    set.seed(seed)
    dtrain <- xgb.DMatrix(as.matrix(temp_mixed[,-c(1,2,3)]), label = temp_mixed$OS.time,weight = temp_mixed$OS)
    params <- list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.01,
      max_depth = max,
      subsample = 0.8,
      colsample_bytree = 0.8
    )
    xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
    
    dtest <- xgb.DMatrix(as.matrix(x[,-c(1,2,3)]), label = x$OS.time,weight = x$OS)
    rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=predict(xgb_model,dtest))
    summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
  }
  temp_vector = c()
  temp_comb = parLapply(cl,list_data,predict_x)
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + XGboost [max_depth=',max,']')))
}


#Lasso + Random survival forest
for (alpha in seq(0,1,0.1)) {
  clusterExport(cl, c("mixed","seed","SPCnfold","coxphitermax","alpha"))
  predict_x = function(x)
  {
    set.seed(seed)
    temp_mixed = mixed[,which(colnames(mixed)%in%colnames(x))]
    fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = temp_mixed[,-c(1)],
                    ntree = 1000,nodesize = 5,
                    splitrule = 'logrank',
                    proximity = T,
                    forest = T,
                    seed = seed)
    fit <- cv.glmnet(as.matrix(temp_mixed[,-c(1,2,3)]), as.matrix(Surv(temp_mixed$OS.time,temp_mixed$OS)),family = "cox",alpha=alpha,nfolds = 10)
    
    rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2 = as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2,3)]),s=fit$lambda.min)))
    summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
  }
  temp_vector = c()
  temp_comb = parLapply(cl,list_data,predict_x)
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('temp_comb$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Enet',' [a=',alpha,']')))
  if(alpha == 0) {rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Ridge')))}
  if(alpha == 1) {rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Lasso')))}
  rm(list=c('fit'))
}

save(final_result,file='result_final.Rdata')


