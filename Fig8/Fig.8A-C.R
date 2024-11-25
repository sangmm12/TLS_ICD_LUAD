
##4.2.2版本跑

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(forcats)
library(dplyr)

library(data.table)
library(limma)
library(edgeR)
library(ggpubr)
library(ggthemes)
library(ggrepel)


rm (list = ls ())

setwd("D:/LUAD/1.volcano_XDJ新/LUAD/TCGA_risk_HHLL/2.GO_new")

library(data.table)

cancer <- "LUAD"

TLS <- c("LAMP3", "CCL2", "CCL3","CCL4", "CCL5", 
         "CCL8", "CCL18", "CCL19", "CCL21","CXCL9", 
         "CXCL10", "CXCL11","CXCR4", "BCL6", "CXCL13", 
         "CD200", "FBLN7", "ICOS", "SGPP2", "SH2D1A", 
         "TIGIT", "PDCD1", "CD4", "CCR5", "CXCR3", "CSF2", 
         "IGSF6", "IL2RA", "CD38", "CD40", "CD5", "MS4A1", "SDC1", 
         "GFI1", "IL1R1", "IL1R2", "IL10", "CCL20", "IRF4", "TRAF6", "STAT5A" ,"TNFRSF17")


ICD<- c("NT5E","ENTPD1","CALR","HMGB1","HSP90AA1","ATG5","BAX",
        "CASP8","PDIA3","EIF2AK3","PIK3CA","CXCR3","IFNA1","IFNB1",
        "IL10","IL6","TNF","CASP1","IL1R1","IL1B","NLRP3",
        "P2RX7","LY96","MYD88","TLR4","CD4","CD8A","CD8B",
        "FOXP3","IFNG","IFNGR1","IL17A","IL17RA","PRF1")

TLSICD <- unique(c(TLS,ICD))


TCGA <- read.table(paste0("D:/LUAD/1.volcano_XDJ新/LUAD/",cancer,"_TCGA.txt"),sep='')
colnames(TCGA) <- TCGA[1,]
TCGA <- TCGA[-1,]
rownames(TCGA) <- TCGA[,1]
TCGA <- TCGA[,-1]

dat_TCGA <- TCGA

# 筛除TCGA的normal
#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(dat_TCGA),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#仅保留肿瘤样本
data = dat_TCGA[,group == 0]



cluster <- read.csv("D:/LUAD/1.volcano_XDJ新/LUAD/TCGA_risk_HHLL/TCGA_risk_HHLL.CSV",header = T)

cluster <- cbind.data.frame(cluster$SampleName,cluster$group)


c1 <- cluster[which(cluster$`cluster$group`=="Low-low"),]
c2 <- cluster[which(cluster$`cluster$group`=="High-high"),]



data_c1 <- as.matrix(data)[,na.omit(match(c1[,1],colnames(data)))]
data_c2 <- as.matrix(data)[,na.omit(match(c2[,1],colnames(data)))]


data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)


rownames(data_merge) <- rownames(data)



###########################

data_merge <- data_merge[complete.cases(data_merge[,2]), ]

d0 <- DGEList(data_merge)

# 注意： calcNormFactors并不会标准化数据，只是计算标准化因子
d0 <- calcNormFactors(d0, method="TMM")
#logCPM<-cpm(d0,log=T,prior.count = 3)

# 过滤低表达
# cutoff <- 22.7097
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# d <- d0[-drop,] 
d <- d0
group_list <- factor(c(rep("c1",dim(data_c1)[2]),
                       rep("c2",dim(data_c2)[2])), levels = c("c1","c2"))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(date)


v <- voom(d, design, plot=F)

fit <- lmFit(v, design)
fit <- eBayes(fit, trend=F)
output <- topTable(fit, coef=2,n=Inf)

res<-subset(output, adj.P.Val<0.05 & abs(logFC) >0.5 )


write.csv(res,file= "DEG_voom-limma.csv")




deg.data <- output

#rownames(deg.data) <- gene.name

deg.data <-as.data.frame(deg.data)

deg.data$log10P <-  (-log10(deg.data$adj.P.Val))
deg.data$Symbol <- rownames(deg.data)

deg.data$Group <- "not-significant"
deg.data$Group[which((deg.data$adj.P.Val<0.05) & (deg.data$logFC > 0.5) )]="up-regulated"
deg.data$Group[which((deg.data$adj.P.Val<0.05) & (deg.data$logFC < -0.5) )]="down-regulated"
table(deg.data$Group)

deg.data$Label=""
deg.data<-deg.data[order(deg.data$adj.P.Val),]



degallgenes <- names(which(table(c(TLSICD,rownames(deg.data)))==2))
deg.top30.genes<- intersect(TLSICD, rownames(res))
deg.data$Label[match(deg.top30.genes,deg.data$Symbol)]<-deg.top30.genes
deg.data <- deg.data[match(degallgenes,rownames(deg.data)),]



# up.genes<-head(deg.data$Symbol[which(deg.data$Group=="up-regulated")],30)
# down.genes<-head(deg.data$Symbol[which(deg.data$Group=="down-regulated")],30)
# deg.top30.genes<-c(as.character(up.genes),as.character(down.genes))
# deg.data$Label[match(deg.top30.genes,deg.data$Symbol)]<-deg.top30.genes
# 
# 
# pdf("volcano1.pdf",width = 10)
# p <- ggscatter(deg.data, x="logFC", y="log10P",
#                color="Group",
#                palette=c("#2f5688","#BBBBBB","#CC0000"),
#                alpha=0.8, size=2.5,
#                label=deg.data$Label,
#                font.label=8,repel=T,
#                xlab="logFC",ylab="-log10(padj)",)+
#   theme_base()+
#   geom_hline(yintercept=1.30,linetype="dashed")+
#   geom_vline(xintercept=c(-1,1),linetype="dashed") +
#   xlim(c(-3, 3)) +
#   labs(x="logFC",y="-log10(padj)")
# p
# dev.off()



#写出up和down

gene_name<-deg.data$Symbol[which(deg.data$Group=="up-regulated")]
write.csv(gene_name,file=paste0(cancer,"-DEGs_up.csv"))


gene_name<-deg.data$Symbol[which(deg.data$Group=="down-regulated")]
write.csv(gene_name,file=paste0(cancer,"-DEGs_down.csv"))


pdf(paste0(cancer,"-HH-LL-volcano.pdf"), width = 10)
p <- ggscatter(deg.data, x = "logFC", y = "log10P",
               color = "Group",
               palette = c("#2f5688", "#BBBBBB", "#CC0000"),
               alpha = 0.8, size = 2.5,
               label = deg.data$Label,
               font.label = 12, repel = TRUE,
               xlab = "logFC", ylab = "-log10(padj)") +
  theme_base() +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  xlim(c(-10, 10)) +
  labs(x = "logFC", y = "-log10(padj)") +
  scale_x_continuous(breaks = c(-10,-8,-6,-4,-2,-1,-0.5,0,0.5,1,2,4,6,8,10), expand = c(0.01, 0))
p
dev.off()




                   ##GO up
###############################################################


gene_name<-deg.data$Symbol[which(deg.data$Group=="up-regulated")]

# write.csv(gene_name,file="DEGs_up.csv")
gene.df <- bitr(gene_name, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)


ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "all",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

write.csv(ego_plot,file=paste0(cancer,"-GO_up_plot.csv"))
  
    ######GO新条形图代码
# 转换富集分析结果为数据框
ego_df <- as.data.frame(ego_plot@result)

# 分别筛选每个类别的Count数最高的前10个GO条目
top10_BP <- ego_df %>% filter(ONTOLOGY == "BP") %>% arrange(desc(Count)) %>% head(10)
top10_CC <- ego_df %>% filter(ONTOLOGY == "CC") %>% arrange(desc(Count)) %>% head(10)
top10_MF <- ego_df %>% filter(ONTOLOGY == "MF") %>% arrange(desc(Count)) %>% head(10)

# 合并每个类别的前10个条目
top10_GO <- rbind(top10_BP, top10_CC, top10_MF)

# 绘制条形图
p <- ggplot(top10_GO, aes(x = reorder(Description, Count), y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +  
  scale_fill_manual(values = c("#25A17C", "#D96622", "#6C73B2")) +
  labs(title = "Top 10 GO Terms by Category",
       x = "GO term",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(),  # 去除次要网格线
        axis.line = element_line(color = "black", size = 0.5))  # 添加坐标轴线

ggsave("GO_up_top10.pdf", plot = p, width = 8, height = 13, dpi = 300) # 保存绘图





  ##原代码
pdf(paste0(cancer,"-GO_up-dotplot.pdf"),height = 9,width=7)
dotplot(ego_plot, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()

#BP#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

# write.csv(ego_plot,file="GO_up_BP_plot.csv")


pdf(paste0(cancer,"-GO_up_BP-dotplot.pdf"),height = 5,width=6)
dotplot(ego_plot)
dev.off()

#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

#write.csv(ego_plot,file="GO_up_CC_dotplot.csv")


pdf(paste0(cancer,"-GO_up_CC-dotplot.pdf"),height = 5,width=6)
dotplot(ego_plot)
dev.off()

#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "MF",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

#write.csv(ego_plot,file="GO_up_MF_dotplot.csv")


pdf(paste0(cancer,"-GO_up_MF-dotplot.pdf"),height = 5,width=6)
dotplot(ego_plot)
dev.off()


                     #####KEGG up
ekk_plot<-enrichKEGG(gene=gene.df$ENTREZID,organism="hsa",
                     keyType = "kegg",pAdjustMethod="none",pvalueCutoff =0.05,qvalueCutoff=1)#"none"


ekk2 <- setReadable(ekk_plot, 'org.Hs.eg.db', 'ENTREZID')
write.csv(ekk2,file=paste0(cancer,"-KEGG_up.csv"))


pdf(paste0(cancer,'-KEGG_up_dotplot.pdf'),height = 5,width=6)
dotplot(ekk2)
dev.off()






            ########GO down###
gene_name<-deg.data$Symbol[which(deg.data$Group=="down-regulated")]
#write.csv(gene_name,file="DEGs_down.csv")

gene.df <- bitr(gene_name, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)


ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "all",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

write.csv(ego_plot,file="GO_down_plot.csv")
   
   ######GO新条形图代码
# 转换富集分析结果为数据框
ego_df <- as.data.frame(ego_plot@result)

# 分别筛选每个类别的Count数最高的前10个GO条目
top10_BP <- ego_df %>% filter(ONTOLOGY == "BP") %>% arrange(desc(Count)) %>% head(10)
top10_CC <- ego_df %>% filter(ONTOLOGY == "CC") %>% arrange(desc(Count)) %>% head(10)
top10_MF <- ego_df %>% filter(ONTOLOGY == "MF") %>% arrange(desc(Count)) %>% head(10)

# 合并每个类别的前10个条目
top10_GO <- rbind(top10_BP, top10_CC, top10_MF)

# 绘制条形图
p <- ggplot(top10_GO, aes(x = reorder(Description, Count), y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +  
  scale_fill_manual(values = c("#25A17C", "#D96622", "#6C73B2")) +
  labs(title = "Top 10 GO Terms by Category",
       x = "GO term",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(),  # 去除次要网格线
        axis.line = element_line(color = "black", size = 0.5))  # 添加坐标轴线

ggsave("GO_down_top10.pdf", plot = p, width = 13, height = 13, dpi = 300) # 保存绘图





    #####原绘图代码
pdf(paste0(cancer,"-GO_down-dotplot.pdf"),height = 10,width=8)
dotplot(ego_plot, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()



#BP#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

write.csv(ego_plot,file="GO_down_BP_plot.csv")

pdf(paste0(cancer,"-GO_down_BP-dotplot.pdf"),height = 5,width=6)
dotplot(ego_plot)
dev.off()

#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"


write.csv(ego_plot,file="GO_down_CC_plot.csv")


pdf(paste0(cancer,"-GO_down_CC-dotplot.pdf"),height = 5,width=6)
dotplot(ego_plot)
dev.off()

#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "MF",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"


write.csv(ego_plot,file="GO_down_MF_plot.csv")


pdf(paste0(cancer,"-GO_down_MF-dotplot.pdf"),height = 5,width=6)
dotplot(ego_plot)
dev.off()



          #####KEGG down###
ekk_plot<-enrichKEGG(gene=gene.df$ENTREZID,organism="hsa",
                     keyType = "kegg",pAdjustMethod="none",pvalueCutoff =0.05,qvalueCutoff=1)#"none"


ekk2 <- setReadable(ekk_plot, 'org.Hs.eg.db', 'ENTREZID')
write.csv(ekk2,file=paste0(cancer,"-KEGG_down.csv"))


pdf(paste0(cancer,'-KEGG_Down_dotplot.pdf'),height = 5,width=6)
dotplot(ekk2)
dev.off()

