

# install.packages("ComplexHeatmap")
# install.packages(c("circlize", "grid", "gridExtra", "heatmap3"))
# library(ComplexHeatmap)

rm(list=ls())
setwd("D:/1.LUAD/4.GSCA/3.TLSICD")

             ###SNV

OS <- read.csv("SNV_OS.csv")
DFI <- read.csv("SNV_DFI.csv")     #原DFS改为DFI
DSS <- read.csv("SNV_DSS.csv")
PFS <- read.csv("SNV_PFS.csv")


# 将HR列中大于1的值赋值为2
OS$HR[OS$HR > 1] <- 2
OS$HR[OS$HR < 1] <- 1
DFI$HR[DFI$HR > 1] <- 2
DFI$HR[DFI$HR < 1] <- 1
DSS$HR[DSS$HR > 1] <- 2
DSS$HR[DSS$HR < 1] <- 1
PFS$HR[PFS$HR > 1] <- 2
PFS$HR[PFS$HR < 1] <- 1

OS$HR[OS$pvalue > 0.05] <- 0
DFI$HR[DFI$pvalue > 0.05] <- 0
DSS$HR[DSS$pvalue > 0.05] <- 0
PFS$HR[PFS$pvalue > 0.05] <- 0

# OS$CancerCode <- gsub("\\(.*\\)", "", OS$CancerCode)
# DSS$CancerCode <- gsub("\\(.*\\)", "", DSS$CancerCode)
# DFI$CancerCode <- gsub("\\(.*\\)", "", DFI$CancerCode)
# PFS$CancerCode <- gsub("\\(.*\\)", "", PFS$CancerCode)


rownames(OS) <- OS$GeneSymbol
rownames(DSS) <- DSS$GeneSymbol
rownames(DFI) <- DFI$GeneSymbol
rownames(PFS) <- PFS$GeneSymbol


# 合并数据框的"HR"列，根据"Name"列的值进行合并
merged_df <- merge(merge(merge(OS, DSS, by = "GeneSymbol", all = TRUE), DFI, by = "GeneSymbol", all = TRUE), PFS, by = "GeneSymbol", all = TRUE)
# 
# # 删除第一列中 ID 为 "name1" 和 "name2" 的行数据
# merged_df<- subset(merged_df, !(GeneSymbol %in% c("TCGA-COADREAD", "TCGA-GBMLGG","TCGA-KIPAN","TCGA-SKCM-M",
#                                                     "TCGA-SKCM-P","TCGA-STES")))  #这6个癌症不属于泛癌

# 选择合并后的数据框中的"HR"列
final_merged_HR <- merged_df[, grepl("HR", names(merged_df))]

rownames(final_merged_HR) <- merged_df$GeneSymbol

colnames(final_merged_HR) <- c("OS","DSS","DFI","PFS")

final_sheet_temp <- final_merged_HR

final_sheet_temp <- data.frame(final_sheet_temp)

final_sheet_temp$gene <- rownames(final_sheet_temp)

colnames(final_sheet_temp)[1:4] = paste0('SNV ',colnames(final_sheet_temp)[1:4])#RNA改为SNV





                 ###CNV
OS <- read.csv("CNV_OS.csv")
DFI <- read.csv("CNV_DFI.csv")
DSS <- read.csv("CNV_DSS.csv")
PFS <- read.csv("CNV_PFS.csv")

# 将HR列中大于1的值赋值为2
OS$HR[OS$HR > 1] <- 2
OS$HR[OS$HR < 1] <- 1
DFI$HR[DFI$HR > 1] <- 2
DFI$HR[DFI$HR < 1] <- 1
DSS$HR[DSS$HR > 1] <- 2
DSS$HR[DSS$HR < 1] <- 1
PFS$HR[PFS$HR > 1] <- 2
PFS$HR[PFS$HR < 1] <- 1

OS$HR[OS$pvalue > 0.05] <- 0
DFI$HR[DFI$pvalue > 0.05] <- 0
DSS$HR[DSS$pvalue > 0.05] <- 0
PFS$HR[PFS$pvalue > 0.05] <- 0

# OS$CancerCode <- gsub("\\(.*\\)", "", OS$CancerCode)
# DSS$CancerCode <- gsub("\\(.*\\)", "", DSS$CancerCode)
# DFI$CancerCode <- gsub("\\(.*\\)", "", DFI$CancerCode)
# PFS$CancerCode <- gsub("\\(.*\\)", "", PFS$CancerCode)

rownames(OS) <- OS$GeneSymbol
rownames(DSS) <- DSS$GeneSymbol
rownames(DFI) <- DFI$GeneSymbol
rownames(PFS) <- PFS$GeneSymbol

# 合并数据框的"HR"列，根据"Name"列的值进行合并
merged_df <- merge(merge(merge(OS, DSS, by = "GeneSymbol", all = TRUE), DFI, by = "GeneSymbol", all = TRUE), PFS, by = "GeneSymbol", all = TRUE)

# 选择合并后的数据框中的"HR"列
final_merged_HR <- merged_df[, grepl("HR", names(merged_df))]


rownames(final_merged_HR) <- merged_df$GeneSymbol

colnames(final_merged_HR) <- c("OS","DSS","DFI","PFS")

final_merged_HR$gene <- rownames(final_merged_HR)

colnames(final_merged_HR)[1:4] = paste0('CNV ',colnames(final_merged_HR)[1:4])


# 加载dplyr包
library(dplyr)

# 使用full_join()根据"gene"列合并数据框
final_sheet_temp <- full_join(final_sheet_temp, final_merged_HR, by = "gene")





             ### METHY
OS <- read.csv("Methy_OS.csv")
DFI <- read.csv("Methy_DFI.csv")
DSS <- read.csv("Methy_DSS.csv")
PFS <- read.csv("Methy_PFS.csv")

# 将HR列中大于1的值赋值为2
OS$HR[OS$HR > 1] <- 2
OS$HR[OS$HR < 1] <- 1
DFI$HR[DFI$HR > 1] <- 2
DFI$HR[DFI$HR < 1] <- 1
DSS$HR[DSS$HR > 1] <- 2
DSS$HR[DSS$HR < 1] <- 1
PFS$HR[PFS$HR > 1] <- 2
PFS$HR[PFS$HR < 1] <- 1

OS$HR[OS$pvalue > 0.05] <- 0
DFI$HR[DFI$pvalue > 0.05] <- 0
DSS$HR[DSS$pvalue > 0.05] <- 0
PFS$HR[PFS$pvalue > 0.05] <- 0

# OS$CancerCode <- gsub("\\(.*\\)", "", OS$CancerCode)
# DSS$CancerCode <- gsub("\\(.*\\)", "", DSS$CancerCode)
# DFI$CancerCode <- gsub("\\(.*\\)", "", DFI$CancerCode)
# PFS$CancerCode <- gsub("\\(.*\\)", "", PFS$CancerCode)

rownames(OS) <- OS$GeneSymbol
rownames(DSS) <- DSS$GeneSymbol
rownames(DFI) <- DFI$GeneSymbol
rownames(PFS) <- PFS$GeneSymbol

# 合并数据框的"HR"列，根据"Name"列的值进行合并
merged_df <- merge(merge(merge(OS, DSS, by = "GeneSymbol", all = TRUE), DFI, by = "GeneSymbol", all = TRUE), PFS, by = "GeneSymbol", all = TRUE)

# 选择合并后的数据框中的"HR"列
final_merged_HR <- merged_df[, grepl("HR", names(merged_df))]


rownames(final_merged_HR) <- merged_df$GeneSymbol

colnames(final_merged_HR) <- c("OS","DSS","DFI","PFS")

final_merged_HR$gene <- rownames(final_merged_HR)

colnames(final_merged_HR)[1:4] = paste0('Methylaton ',colnames(final_merged_HR)[1:4])

# 加载dplyr包
library(dplyr)

# 使用full_join()根据"gene"列合并数据框
final_sheet_temp <- full_join(final_sheet_temp, final_merged_HR, by = "gene")




                            ###RNA
OS <- read.csv("RNAsanger_OS.csv")
DFI <- read.csv("RNAsanger_DFI.csv")
DSS <- read.csv("RNAsanger_DSS.csv")
PFI <- read.csv("RNAsanger_PFI.csv")

# 将HR列中大于1的值赋值为2
OS$HR[OS$HR > 1] <- 2
OS$HR[OS$HR < 1] <- 1
DFI$HR[DFI$HR > 1] <- 2
DFI$HR[DFI$HR < 1] <- 1
DSS$HR[DSS$HR > 1] <- 2
DSS$HR[DSS$HR < 1] <- 1
PFI$HR[PFI$HR > 1] <- 2
PFI$HR[PFI$HR < 1] <- 1

OS$HR[OS$pvalue > 0.05] <- 0
DFI$HR[DFI$pvalue > 0.05] <- 0
DSS$HR[DSS$pvalue > 0.05] <- 0
PFI$HR[PFI$pvalue > 0.05] <- 0

rownames(OS) <- OS$GeneSymbol
rownames(DSS) <- DSS$GeneSymbol
rownames(DFI) <- DFI$GeneSymbol
rownames(PFI) <- PFI$GeneSymbol


# 合并数据框的"HR"列，根据"Name"列的值进行合并
merged_df <- merge(merge(merge(OS, DSS, by = "GeneSymbol", all = TRUE), DFI, by = "GeneSymbol", all = TRUE), PFI, by = "GeneSymbol", all = TRUE)

# 选择合并后的数据框中的"HR"列
final_merged_HR <- merged_df[, grepl("HR", names(merged_df))]

rownames(final_merged_HR) <- merged_df$GeneSymbol

colnames(final_merged_HR) <- c("OS","DSS","DFI","PFI")

final_merged_HR$gene <- rownames(final_merged_HR)

colnames(final_merged_HR)[1:4] = paste0('RNA ',colnames(final_merged_HR)[1:4])


# 加载dplyr包
library(dplyr)

# 使用full_join()根据"gene"列合并数据框
final_sheet_temp <- full_join(final_sheet_temp, final_merged_HR, by = "gene")

# # 删除第一列中 ID 为 "name1" 和 "name2" 的行数据
# final_sheet_temp<- subset(final_sheet_temp, !(gene %in% c("TCGA-COADREAD", "TCGA-GBMLGG","TCGA-KIPAN","TCGA-SKCM-M",
# "TCGA-SKCM-P","TCGA-STES")))

# # 去除 "TCGA-" 前缀
# final_sheet_temp$gene <- gsub('TCGA-', '', final_sheet_temp$gene)

rownames(final_sheet_temp) <- final_sheet_temp$gene

final_sheet_temp <- final_sheet_temp[,-5]

final_sheet <- final_sheet_temp

final_sheet <- as.matrix(final_sheet)
final_sheet_temp <- as.matrix(final_sheet_temp)




library(ComplexHeatmap)
colors = structure(c("white","#00A087B2","#BB362F"), names = c(0,1,2))
row_boxplot <- rowAnnotation(
  . = anno_barplot(rowSums(final_sheet_temp)[order(-rowSums(final_sheet_temp))]/2,
                   bar_width = 0.8,
                   axis = F,
                   border=F,
                   gp = gpar(fill = 1:10),
                   add_numbers=T,
                   numbers_gp = gpar(fontsize = 12),
                   numbers_offset = unit(4, "mm"),
                   width = unit(20, "mm"))
)




pdf('complexheatmap_bind4.pdf',width=24,height=40)
Heatmap(final_sheet,
        na_col='grey',
        col=colors,
        cluster_columns = F,
        cluster_rows = F,
        height=nrow(final_sheet)*unit(8, "mm"),
        width=ncol(final_sheet)*unit(10, "mm"),
        right_annotation = row_boxplot,
        rect_gp = gpar(col = "black", lwd = 1),
        column_split = c(rep("SNV",4),rep("CNV",4),rep("Methylation",4),rep("RNA",4)),
        row_names_side = 'left',
        column_gap = unit(5, "mm"),
        column_names_gp = gpar(col = c("#94ffca", "#c9dbff", "#e5ffaa","#b4f1ff")),
        column_title_gp = gpar(fill = c("#94ffca", "#c9dbff", "#e5ffaa","#b4f1ff")))
dev.off()





