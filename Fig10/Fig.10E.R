
library(rstatix)
setwd("D:/LUAD/16.TIP/1.TIP")
library(data.table)
TIP <-  fread(paste("D:/LUAD/16.TIP/1.TIP/TIP_LUAD.txt",sep=''))
TIP  <- as.data.frame(TIP)
rownames(TIP) <- TIP[,1]
TIP <- TIP[,-1]
dat_immu <- t(TIP)
aa <- substr(rownames(dat_immu), 1, 15)
rownames(dat_immu) <- aa
# 将dat_immu第一列样本名中的"."替换�?"-"
sample_names <- rownames(dat_immu)
sample_names <- gsub("\\.", "-", sample_names)
rownames(dat_immu) <- sample_names


cluster <- read.csv("D:/LUAD/16.TIP/1.TIP/scale_TCGA_HHLL.csv",header = T) ##归一化后risk
cluster <- as.data.frame(cluster)
cluster<-cluster[,c(1,10)] #�?10行group
colnames(cluster)<-c("V1","V2")
rownames(cluster) <- cluster[,1]

all_name <- names(which(table(c(rownames(dat_immu),cluster[,1] ))==2))

dat_cluster <- cluster[match(all_name,cluster[,1]),]
dat_imm <- dat_immu[match(all_name,rownames(dat_immu)),]
dat_im <- as.data.frame(dat_imm)

dat_cluster <- na.omit(dat_cluster)

all_name <- names(which(table(c(rownames(dat_im),dat_cluster[,1] ))==2))


dat_cluster <- dat_cluster[match(all_name,dat_cluster[,1]),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]

dat_cluster <- dat_cluster[c('V2')]

dat <- data.frame()
for(coln in colnames(dat_im))
{
  for(i in all_name)
  {
    dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],dat_im[c(coln)][match(i,rownames(dat_im)),]))
  }
}

dat[,3] = as.numeric(dat[,3])


colnames(dat) <- c("Gene","Group","value")
aa <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")

library(ggpubr)
pdf(paste("cluster_TIP .pdf",sep=''),width=16,height = 6)
p <- ggboxplot(dat, x = "Gene", y = "value",
               color = "Group", palette = c("#FF6666", "#3399FF"), 
               add = "median_q1q3",x.text.angle=60) # palette可以按照期刊选择相应的配色，�?"npg"�?
p <- p+xlab("Step")+ylab("Expression of TIP")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
p
#print(p)
dev.off()

