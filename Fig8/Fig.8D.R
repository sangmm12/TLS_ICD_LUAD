
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
#BiocManager::install("GSVA",type="binary")
#BiocManager::install("msigdbr")
library(GSVA)
library(GSEABase)
library(limma)
library(Seurat)
library(msigdbr)
setwd("D:/1.LUAD/14.GSVA/1.GSVA/11.13")

#all_gene_sets = msigdbr(species = "Mus musculus")
human <- msigdbr(species = "Homo sapiens")
#hallmark
human_GO_bp = msigdbr(species = "human",
                      category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)


human_GO_bp_Set = human_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)


s.sets <- human_GO_bp_Set


###########

library(data.table)                      ##TCGA所有基因几万个的Count表达（整数）
data <- read.table("D:/1.LUAD/14.GSVA/1.GSVA/11.13/LUAD_Count_Exp.txt",row.names = 1)
colnames(data) <- data[1,]
data <- data[-1,]

cluster <- read.csv("D:/1.LUAD/14.GSVA/1.GSVA/11.13/riskTCGA.csv",header = T) #改为header=T,有表头
c1 <- cluster[which(cluster$group=="High-high"),]
c2 <- cluster[which(cluster$group=="Low-low"),]

data_c1 <- as.matrix(data)[,na.omit(match(c1[,1],colnames(data)))]
data_c2 <- as.matrix(data)[,na.omit(match(c2[,1],colnames(data)))]


data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)


rownames(data_merge) <- rownames(data)

##########

expr1 <- as.matrix(data_merge)


expr <- expr1



es.matrix <- gsva(
  expr,
  s.sets,
  min.sz = 10,
  max.sz = Inf,
  tau = 1,
  method = "gsva",
  abs.ranking = FALSE,
  verbose = TRUE,
  parallel.sz = 1
)

saveRDS(es.matrix,file="hallmark.rds")



n1 <- 1:dim(data_c1)[2]
#grep("male",seurat_obj@meta.data$sex)

n2 <- dim(data_c1)[2]+1:dim(data_c2)[2]

#grep("female",seurat_obj@meta.data$sex)

es.matrix.1 <-
  as.data.frame(es.matrix[, n1],
                row.names = row.names(es.matrix))
es.matrix.2 <-
  as.data.frame(es.matrix[, n2],
                row.names = row.names(es.matrix))


es.matrix.f <- cbind(es.matrix.1, es.matrix.2)

grouP <-
  c(rep("case", dim(es.matrix.1)[2]),
    rep("control", dim(es.matrix.2)[2]))

grouP <- as.factor( grouP)
design <- model.matrix(~ grouP + 0)



row.names(design)<-c(colnames(es.matrix.1), colnames(es.matrix.2))

comparE <-
  makeContrasts(grouPcase - grouPcontrol, levels = design)

fit <- lmFit(es.matrix, design)
fit2 <- contrasts.fit(fit, comparE)
fit3 <- eBayes(fit2)


diff <- topTable(fit3, coef = 1, number = dim(es.matrix)[1])

t_results <-
  # as.data.frame(diff$t, row.names = rownames(es.matrix))
  as.data.frame(diff$t, row.names = rownames(diff))  #2024.11.13改

head(t_results)
colnames(t_results) <- c("t_value")




saveRDS(t_results,file="t_results_hallmark.rds")
write.csv(t_results,file="t_results_hallmark.csv")

#t_results <- readRDS("t_results_c5.go.bp_Elongating.rds")

library(ggplot2)

library(pheatmap)
rownames(t_results) <- gsub("HALLMARK_", "", rownames(t_results))
rownames(t_results) <- gsub("_", " ", rownames(t_results))
focus.cluster <- "t_value"
sub_t_results <- as.data.frame(t_results[, focus.cluster],
                               row.names = rownames(t_results))
sub_t_results$hallmark <- rownames(sub_t_results)
colnames(sub_t_results) <- c("t", "hallmark")

sub_t_results$hallmark = with(sub_t_results, reorder(hallmark, t))
sub_t_results$fill <- ""
sub_t_results[sub_t_results$t >= 2.58,]$fill <-
  "up"
sub_t_results[sub_t_results$t <= -2.58,]$fill <-
  "down"
sub_t_results[abs(sub_t_results$t) < 2.58,]$fill <-
  "no"
sub_t_results$color <- ""
sub_t_results[abs(sub_t_results$t) < 2.58,]$color <-
  "n"
sub_t_results[abs(sub_t_results$t) >= 2.58,]$color <-
  "y"


sub_t_results <- sub_t_results[c(1:50),]

p <-
  ggplot(sub_t_results, aes(x = hallmark, y = t, fill = fill)) +
  geom_bar(stat = "identity", width = .8) +
  scale_fill_manual(
    values = c(
      "down" = "#36648b",
      "up" = "#e94644",
      "no" = "#cccccc"
    ),
    guide = F
  ) + ylim(-5,11.5)+     #改最小最大值
  geom_hline(
    yintercept = c(-2.58, 2.58),
    color = "white",
    linetype = "dotted",
    size = 0.5
  ) +
  coord_flip() +
  xlab("") +
  geom_text(
    data = subset(sub_t_results, t < 0),
    aes(
      x = hallmark,
      y = 0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 1.8,
    hjust = "inward"  #互换
  ) +geom_text(
    data = subset(sub_t_results, t > 0),
    aes(
      x = hallmark,
      y = -0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 1.8,
    hjust = "outward"  #互换
  ) +
  scale_colour_manual(values = c("y" = "black", "n" = "#cccccc"),
                      guide = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),
  )
#p
ggsave(
  filename = "hallmark.pdf",
  plot = p,
  height = 7.5,
  width =6
)

