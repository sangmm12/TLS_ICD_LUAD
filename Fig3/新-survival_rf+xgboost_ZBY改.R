
rm(list = ls())
library(data.table)
library(dplyr)
library(survival)
library(xgboostExplainer)
library(randomForestSRC)
library(xgboost)

setwd('D:/1.LUAD/6.2.GSE单独去批次MachineLearning/8.TLS_10per_SurvivalRF_XGboost.max_depth1/2.模型生存')

Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

final_list <- data.frame()

# sample_name <- "GSE19188"
for(sample_name in c("TCGA","GSE11969","GSE19188","GSE30219","GSE31210","GSE37745","GSE42127","GSE50081","GSE68465","GSE72094"))
{
  name_fold = paste0('risk/',sample_name)
  if (!dir.exists(name_fold)){
    dir.create(name_fold)
  } else {
    print("Dir already exists!")
  }
  
  all_sur_data <- fread(paste('D:/1.LUAD/6.2.GSE单独去批次MachineLearning/8.TLS_10per_SurvivalRF_XGboost.max_depth1/2.模型生存/sur_',sample_name,'.csv',sep=''),header=T)
  gene_exp <- fread(paste('D:/1.LUAD/6.2.GSE单独去批次MachineLearning/8.TLS_10per_SurvivalRF_XGboost.max_depth1/2.模型生存/Exp_',sample_name,'_TLS.csv',sep=''),header=T)
  # gene_exp = gene_exp[,which(!colnames(gene_exp)=='ADGRE5')]
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','Time',colnames(mixed))
  colnames(mixed) <- sub('Status','Status',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  mixed$SampleName <- temp_row
  mixed <- mixed %>% 
    dplyr::select( SampleName,Status, Time, everything())
  mixed <- na.omit(mixed)
  
  
  if(sample_name=='TCGA')
  {
    mixed_TCGA=mixed
    mixed_TCGA_temp = mixed
  } else
  {
    mixed_TCGA_temp = mixed_TCGA[,colnames(mixed)]
    print(colnames(mixed_TCGA_temp))
    print(colnames(mixed))
    if(file.exists(paste(name_fold,'sample.Rda')))
    {
      load(paste(name_fold,'sample.Rda'))
      print('load....')
    } else {
      mixed_sample = sample(length(rownames(mixed)),length(rownames(mixed))*0.10)
      save(mixed_sample,file=paste(name_fold,'sample.Rda'))
    }
    smixed = mixed[mixed_sample,]     #改比例
    mixed_TCGA_temp = rbind(mixed_TCGA_temp,smixed)
  }
  
  set.seed(123456)
  
  #测试集mixed  训练集mixed_TCGA_temp  
###########################################
  fit_rf <- rfsrc(Surv(Time,Status)~.,data = mixed_TCGA_temp[,-c(1)],
                  ntree = 1000,nodesize = 5,
                  splitrule = 'logrank',
                  proximity = T,
                  forest = T,
                  seed = 123456)
 
  dtest <- xgb.DMatrix(as.matrix(mixed[,-c(1,2,3)]), label = mixed$Time,weight = mixed$Status)
  
  params <- list(
    objective = "survival:cox",
    eval_metric = "cox-nloglik",
    eta = 0.01,
    max_depth = 1,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  xgb_model <- xgb.train(params = params,nrounds = 100,dtest)

  
  x_p = Fun(as.numeric(predict(xgb_model,dtest))+100)
  y_p = Fun(as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted)+100)
#######################################
  
  
  
  rs <- data.frame(Time=mixed$Time,Status=mixed$Status,RS1=x_p,RS2=y_p)
  
  Risk_scroe_temp <- summary(coxph(Surv(Time,Status)~RS1+RS2,rs))$coef[3]*x_p+summary(coxph(Surv(Time,Status)~RS1+RS2,rs))$coef[4]*y_p
  
  threshold <- median(Risk_scroe_temp)
  
  risk <- data.frame(riskscore = Fun(Risk_scroe_temp)*100,group = risk_group <- ifelse((Risk_scroe_temp) > threshold, "High", "Low"),Status=mixed$Status,Time=mixed$Time/365)
  rownames(risk)=mixed$SampleName
  write.csv(risk,file=paste(name_fold,'/risk.csv',sep=''))
  final_list <- rbind(final_list,data.frame(CODE=rep(sample_name,length(risk$group)),SampleName=mixed$SampleName,Risk=risk$group,Status=risk$Status))
  
  rt=risk[order(risk$riskscore),]
  riskClass=rt$group
  lowLength=length(which(riskClass=="Low"))
  highLength=length(which(riskClass=="High"))
  lowMax=max(rt[which(rt$group=="Low"),]$riskscore)
  line=rt$riskscore
  line[line>1000]=1000
  pdf(file=paste(name_fold,'/risk line.pdf',sep=''),width = 6,height = 6)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("lightblue",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
  dev.off()
  color=as.vector(rt$Status)
  color[color==1]="red"
  color[color==0]="lightblue"
  pdf(file=paste(name_fold,'/risk point.pdf',sep=''),width = 6,height = 6)
  plot(rt$Time, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  library(timeROC)
  predict_full <- risk$riskscore
  full_time <- mixed$Time/365
  ROC_rt=timeROC(T=full_time,delta=mixed$Status,
                 marker=predict_full,cause=1,
                 weighting='aalen',
                 times=c(1,2,3),ROC=TRUE)
  pdf(file=paste(name_fold,'/roc.pdf',sep=''),width = 6,height = 6)
  plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
  dev.off()
  
  
  ###原生存曲线代码###  
#   library(survminer)
#   rt=risk[order(risk$riskscore),]
#   diff=survdiff(Surv(Time, Status) ~ group,data = rt)
#   pValue=1-pchisq(diff$chisq,df=1)
#   pValue=signif(pValue,4)
#   pValue=format(pValue, scientific = TRUE)
#   fit <- survfit(Surv(Time, Status) ~ group, data = rt)
#   surPlot=ggsurvplot(fit,
#                      conf.int = TRUE,
#                      data=rt,
#                      pval=paste0("p=",pValue),
#                      pval.size=5,
#                      legend.labs=c("High risk", "Low risk"),
#                      legend.title="Risk",
#                      xlab="Time(years)",
#                      break.time.by = 1,
#                      risk.table.title="",
#                      risk.table=F,
#                      risk.table.height=.25)
#   pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 6,height =6)
#   print(surPlot)
#   dev.off()
# }
# #write.csv(final_list,file='riskRFTcell.csv',row.names=F)

##新颜色  
  library(survminer)
  rt=risk[order(risk$riskscore),]
  rt <- rt[order(rt$group), ] #另加
  diff=survdiff(Surv(Time, Status) ~ group,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(Time, Status) ~ group, data = rt)
  surPlot=ggsurvplot(fit,
                     data=rt,
                     #font.title = paste(connam[i],sep=''),
                     #ggtitle = paste(connam[i],sep=''),
                     #conf.int=TRUE,
                     legend.labs=c(unique(rt$group)),
                     legend = "top",
                     legend.title="Risk",
                     pval=paste0("p=",pValue),
                     pval.size=5,
                     xlab="Time(years)",
                     break.time.by = ceiling((max(rt$Time))/4),
                     risk.table.title="",
                     palette=c("red","green"),
                     risk.table=T,
                     risk.table.height=.25,)
  pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 6,height =6)
  print(surPlot)
  dev.off()
}


  
  