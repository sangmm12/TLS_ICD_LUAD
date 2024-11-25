

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


setwd("D:/LUAD/1.volcano_XDJ新/LUAD/TN/GO_KEGG/test")

cancer <- "LUAD"

TLS <- c("LAMP3", "CCL2", "CCL3","CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21","CXCL9", "CXCL10", "CXCL11",
         "CXCR4", "BCL6","CXCL13", "CD200", "FBLN7", "ICOS", "SGPP2", "SH2D1A", "TIGIT", "PDCD1", "CD4", "CCR5",
         "CXCR3", "CSF2", "IGSF6", "IL2RA", "CD38", "CD40","CD5", "MS4A1", "SDC1", "GFI1", "IL1R1", "IL1R2", 
         "IL10", "CCL20", "IRF4", "TRAF6", "STAT5A" ,"TNFRSF17")
ICD<- c("NT5E","ENTPD1","CALR","HMGB1","HSP90AA1","ATG5","BAX","CASP8","PDIA3","EIF2AK3","PIK3CA","CXCR3","IFNA1",
        "IFNB1","IL10","IL6","TNF","CASP1","IL1R1","IL1B","NLRP3","P2RX7","LY96","MYD88","TLR4","CD4","CD8A","CD8B",
        "FOXP3","IFNG","IFNGR1","IL17A","IL17RA","PRF1")

geneall <- unique(c(TLS,ICD))

TCGA <- read.table(paste0("D:/LUAD/1.volcano_XDJ新/LUAD/TN/GO_KEGG/",cancer,"_TCGA.txt"),sep='')
colnames(TCGA) <- TCGA[1,]
TCGA <- TCGA[-1,]
rownames(TCGA) <- TCGA[,1]
TCGA <- TCGA[,-1]

GTEx <- read.table(paste0("D:/LUAD/1.volcano_XDJ新/LUAD/TN/GO_KEGG/",cancer,"_GTEx.txt"),sep='')
colnames(GTEx) <- GTEx[1,]
GTEx <- GTEx[-1,]
rownames(GTEx) <- GTEx[,1]
GTEx <- GTEx[,-1]

all_name <- names(which(table(c(rownames(TCGA),rownames(GTEx) ))==2))  #两个数据库匹配


dat_TCGA <- TCGA[match(all_name,rownames(TCGA)),]
dat_GTEx <- GTEx[match(all_name,rownames(GTEx)),]


# 筛除TCGA的normal
#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(dat_TCGA),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#仅保留正常样本
dat_TCGA_normal = dat_TCGA[,group == 1]     #保留TCGA数据库中59Normal样本
#仅保留肿瘤样本
dat_TCGA= dat_TCGA[,group == 0]             #保留TCGA数据库中515Tumor样本

dat_GTEx <- cbind(dat_GTEx,dat_TCGA_normal) #给GTEx赋新值，把313个GTEx样本和59TCGA正常样本合并
#dat_GTEx意味着所有正常样本



###control为C1
c1 <- colnames(dat_GTEx)
c2 <- colnames(dat_TCGA)

data_c1 <- dat_GTEx
data_c2 <- dat_TCGA

data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)

rownames(data_merge) <- rownames(dat_GTEx)


###########################


data_merge <- data_merge[complete.cases(data_merge[,2]), ]

d0 <- DGEList(data_merge)

# 注意： calcNormFactors并不会标准化数据，只是计算标准化因子
d0 <- calcNormFactors(d0, method="TMM")
#logCPM<-cpm(d0,log=T,prior.count = 3)

# 过滤低表达
# cutoff <- 37.9404
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

write.csv(res,file= paste0(cancer,"-DEG-voom-limma.csv"))


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

degallgenes <- names(which(table(c(geneall,rownames(deg.data)))==2))

deg.top30.genes<- intersect(geneall, rownames(res))

deg.data$Label[match(deg.top30.genes,deg.data$Symbol)]<-deg.top30.genes

deg.data <- deg.data[match(degallgenes,rownames(deg.data)),]

write.csv(deg.data,file= paste0(cancer,"-72gene.csv"))

deg.data_deg <- deg.data[match(deg.top30.genes,rownames(deg.data)),]

write.csv(deg.data_deg,file= paste0(cancer,"-Diffgene.csv"))



#写出up和down

gene_name<-deg.data$Symbol[which(deg.data$Group=="up-regulated")]
write.csv(gene_name,file=paste0(cancer,"-DEGs_up.csv"))

gene_name<-deg.data$Symbol[which(deg.data$Group=="down-regulated")]
write.csv(gene_name,file=paste0(cancer,"-DEGs_down.csv"))


pdf(paste0(cancer,"-volcano.pdf"), width = 10)
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


  ####GO up
gene_name<-deg.data$Symbol[which(deg.data$Group=="up-regulated")]

# write.csv(gene_name,file="DEGs_up.csv")
gene.df <- bitr(gene_name, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)


ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "all",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file=paste0(cancer,"-GO_up_plot.csv"))

  ##2024.5.21新条形图
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

ggsave("GO_up_top10.pdf", plot = p, width = 6, height = 13, dpi = 300) # 保存绘图
#

# # 将富集分析结果转换为数据框
# ego_df <- as.data.frame(ego_plot@result)
# 
# # 分别筛选每个类别的Count数最高的前10个GO条目
# top10_BP <- ego_df %>% filter(ONTOLOGY == "BP") %>% arrange(desc(Count)) %>% head(10)
# top10_CC <- ego_df %>% filter(ONTOLOGY == "CC") %>% arrange(desc(Count)) %>% head(10)
# top10_MF <- ego_df %>% filter(ONTOLOGY == "MF") %>% arrange(desc(Count)) %>% head(10)
# 
# # 合并每个类别的前10个条目
# top10_GO <- rbind(top10_BP, top10_CC, top10_MF)
# 
# # 绘制条形图
# p <- ggplot(top10_GO, aes(x = reorder(Description, Count), y = Count, fill = ONTOLOGY)) +
#   geom_bar(stat = "identity", width = 0.8) +
#   coord_flip() +
#   facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
#   scale_fill_manual(values = c("#25A17C", "#D96622", "#6C73B2")) +
#   labs(title = "Top 10 GO Terms by Category",
#        x = "GO term",
#        y = "Count") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# ggsave("GO_up_top10.pdf", plot = p, width = 6, height = 13, dpi = 300) #保存绘图


 # # 筛选每个类别的前10个GO条目
 # top10_BP <- ego_df %>% filter(ONTOLOGY == "BP") %>% top_n(-10, p.adjust)
 # top10_CC <- ego_df %>% filter(ONTOLOGY == "CC") %>% top_n(-10, p.adjust)
 # top10_MF <- ego_df %>% filter(ONTOLOGY == "MF") %>% top_n(-10, p.adjust)
 # 
 # # 合并每个类别的前10个条目
 # top10_GO <- rbind(top10_BP, top10_CC, top10_MF)
 # 
 # # 创建一个新列用于按-log10(p.adjust)排序
 # top10_GO$log10Padj <- -log10(top10_GO$p.adjust)
 # 
 # # 绘制条形图
 # p <- ggplot(top10_GO, aes(x = reorder(Description, -log10Padj), y = log10Padj, fill = ONTOLOGY)) +
 #   geom_bar(stat = "identity", width = 0.8) +
 #   coord_flip() +
 #   facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
 #   scale_fill_manual(values = c("#25A17C", "#D96622", "#6C73B2")) +
 #   labs(title = "Top 10 GO Terms by Category",
 #        x = "GO Term",
 #        y = "-log10(p.adjust)") +
 #   theme_minimal() +
 #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
 # 
 # ggsave("GO_up_top10.pdf", plot = p, width = 6, height = 13, dpi = 300) #保存绘图

# pdf("GO_up_top10.pdf", width = 14, height = 10) #打印绘图到PDF
# print(p)
# dev.off()



# ###新条形图代码
# WEGO2 <- cbind(ego_plot@result$ID,ego_plot@result$ONTOLOGY,ego_plot@result$Count,
#                ego_plot@result$GeneRatio,ego_plot@result$Description,ego_plot@result$Count)
# colnames(WEGO2) <- c("ID", "Level1", "Number of genes", "Percentage of genes", "Description", "Position")
# head(WEGO2)
# class(WEGO2)
# WEGO2 <- data.frame(WEGO2)
# head(WEGO2)
# class(WEGO2)
# 
# #Percentage of genes分数转化为小数
# WEGO2$Percentage.of.genes <- sapply(WEGO2$Percentage.of.genes, function(x) {
#   parts <- strsplit(x, "/")[[1]]  # 将字符串按斜杠分割
#   numerator <- as.numeric(parts[1])  # 分子
#   denominator <- as.numeric(parts[2])  # 分母
#   numerator / denominator  # 计算小数
# })  
# print(WEGO2) #打印结果
# 
# 
# maxp <- max(WEGO2$Percentage.of.genes)
# class(maxp)
# class(WEGO2$Number.of.genes)
# class(WEGO2$Percentage.of.genes)
# WEGO2$Number.of.genes <- as.numeric(WEGO2$Number.of.genes)
# class(WEGO2$Number.of.genes)
# 
# normalizer <- max(WEGO2$Number.of.genes) / max(WEGO2$Percentage.of.genes)

#  ##绘图
# ggplot(data = WEGO2,
#        aes(x = fct_reorder(Description, Position),
#            y = Percentage.of.genes,
#            fill = Level1)) +
#   scale_fill_manual(values = c("#25A17C","#D96622","#6C73B2")) +
#   geom_col(data = WEGO2,
#            width = 0.5,
#            position = "dodge",
#            aes(x = fct_reorder(Description, Position),
#                y = Number.of.genes/normalizer)) +
#   
#   # 双坐标轴的关键属性sec.axis
#   scale_y_continuous(limits = c(0, maxp+3),
#                      sec.axis = sec_axis(trans = ~.*normalizer,
#                                          name = "Number of genes"),
#                      expand = c(0,0)) +
#   theme_classic() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_text(size = 10),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 8),
#         legend.key.size =  unit(0.2, "in"),
#         legend.position = "none",
#         panel.border = element_blank(),
#         axis.line.y  = element_line(color = "black",
#                                     size = 1,
#                                     lineend = "square"),
#         plot.title = element_text(size = 15,
#                                   hjust = 0.5)) +
#   labs(x = NULL,
#        title = "Gene Ontology Classification",
#        fill = NULL) -> p
# print(p)


# ##分面
# library(ggh4x)
# strip <- strip_themed(background_x = elem_list_rect(color = c("#25A17C","#D96622","#6C73B2")))
# p + facet_wrap2(~Level1,
#                 nrow = 1,
#                 scales = "free_x",
#                 strip.position = "bottom",
#                 strip = strip) +
#   geom_text(aes(y = rep(1, nrow(WEGO2)),
#                 label = Description),
#             size = 2,
#             angle = 90,
#             hjust = 0) +
#   theme(axis.line.x = element_blank(),
#         strip.background = element_rect(size = 2),
#         strip.text = element_text(face = "bold",
#                                   size = rel(1.5))) +
#   annotate(geom = 'segment',
#            y = Inf,
#            yend = Inf,
#            color = 'black',
#            x = -Inf,
#            xend = Inf,
#            linewidth = 2)
# 
# ggsave("go_classify_plot.pdf",width = 14,height = 6,dpi = 300)




#原代码
pdf(paste0(cancer,"-GO_up-dotplot.pdf"),height = 20,width=10)
dotplot(ego_plot, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()

#BP#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

# write.csv(ego_plot,file="GO_up_BP_plot.csv")

pdf(paste0(cancer,"-GO_up_BP-dotplot.pdf"),height =4.5,width=6)
dotplot(ego_plot)
dev.off()

#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

#write.csv(ego_plot,file="GO_up_CC_dotplot.csv")

pdf(paste0(cancer,"-GO_up_CC-dotplot.pdf"),height = 2,width=6)
dotplot(ego_plot)
dev.off()

#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "MF",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

#write.csv(ego_plot,file="GO_up_MF_dotplot.csv")

pdf(paste0(cancer,"-GO_up_MF-dotplot.pdf"),height = 4.5,width=6)
dotplot(ego_plot)
dev.off()



    ###KEGG up
ekk_plot<-enrichKEGG(gene=gene.df$ENTREZID,organism="hsa",
                     keyType = "kegg",pAdjustMethod="none",pvalueCutoff =0.05,qvalueCutoff=1)#"none"

ekk2 <- setReadable(ekk_plot, 'org.Hs.eg.db', 'ENTREZID')
write.csv(ekk2,file=paste0(cancer,"-KEGG_up.csv"))

pdf(paste0(cancer,'-KEGG_up_dotplot.pdf'),width=6.5)
dotplot(ekk2)
dev.off()



     ####GO down
gene_name<-deg.data$Symbol[which(deg.data$Group=="down-regulated")]
#write.csv(gene_name,file="DEGs_down.csv")

gene.df <- bitr(gene_name, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)

ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "all",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

write.csv(ego_plot,file="GO_down_plot.csv")

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
        axis.line = element_line(color = "black", linewidth = 0.5))  # 添加坐标轴线

ggsave("GO_down_top10.pdf", plot = p, width = 13, height = 13, dpi = 300) # 保存绘图
#



#原代码
pdf(paste0(cancer,"-GO_down-dotplot.pdf"),height = 20,width=10)
dotplot(ego_plot, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()


#BP#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"

write.csv(ego_plot,file="GO_down_BP_plot.csv")

pdf(paste0(cancer,"-GO_down_BP-dotplot.pdf"),height = 7,width=8)
dotplot(ego_plot)
dev.off()

#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"


write.csv(ego_plot,file="GO_down_CC_plot.csv")

pdf(paste0(cancer,"-GO_down_CC-dotplot.pdf"),height = 4.5,width=6)
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



EGG down
ekk_plot<-enrichKEGG(gene=gene.df$ENTREZID,organism="hsa",
                     keyType = "kegg",pAdjustMethod="none",pvalueCutoff =0.05,qvalueCutoff=1)#"none"


ekk2 <- setReadable(ekk_plot, 'org.Hs.eg.db', 'ENTREZID')
write.csv(ekk2,file=paste0(cancer,"-KEGG_down.csv"))


pdf(paste0(cancer,'-KEGG_Down_dotplot.pdf'),width=6.5)
dotplot(ekk2)
dev.off()

