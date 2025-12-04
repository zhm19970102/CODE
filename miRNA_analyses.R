```r
# Microarray data visualization
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

# Load packages
library(GEOquery)
library(dplyr)
library(tidyverse)
library(data.table)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ggvenn)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(ggraph)

# Load dataset 1: GSE28858
gset1 <- getGEO(filename = "GSE28858_series_matrix.txt.gz", 
                GSEMatrix = TRUE, 
                getGPL = TRUE)

pdata1 = pData(gset1)
exprSet1 = exprs(gset1)
anno_data1 = gset1@featureData@data

rownames(exprSet1) = anno_data1$ILMN_Gene
pdata1 = pdata1 %>% filter(.data[["age:ch1"]] <= 50)
exprSet1 = exprSet1[, colnames(exprSet1) %in% rownames(pdata1)]
anno_exprSet1 = as.data.frame(exprSet1)
pdata1$Group = c(rep('PCAD', 10), rep('non-PCAD', 9))
pdata1 = pdata1[, c(1, 37)]
rownames(pdata1) = pdata1[, 1]

anno_exprSet1 = na.omit(anno_exprSet1)
group_list1 = factor(pdata1$Group)
group_list1 <- relevel(group_list1, ref = "non-PCAD")

exp_anno1 = normalizeBetweenArrays(anno_exprSet1)

design1 = model.matrix(~ group_list1)
colnames(design1) <- levels(group_list1)
rownames(design1) <- colnames(exp_anno1)

fit1 = lmFit(exp_anno1, design1)
fit1 = eBayes(fit1)
allDiff1 = topTable(fit1, coef = 2, adjust = 'fdr', number = Inf)
write.table(allDiff1, file = "allDiff1.txt", sep = "\t", col.names = NA)

# PCA analysis for GSE28858
exprSet_clean1 <- na.omit(exp_anno1)
exprSet_filtered1 <- exprSet_clean1[rowSums(exprSet_clean1) > 0, ]
pca_result1 <- prcomp(t(exprSet_filtered1), scale. = TRUE, center = TRUE)
pca_summary1 <- summary(pca_result1)
pca_data1 <- as.data.frame(pca_result1$x)
pca_data1$Sample <- colnames(exprSet_filtered1)
pca_data1$Group <- c(rep('PCAD', 10), rep('non-PCAD', 9))

ggplot(pca_data1, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  labs(x = paste("PC1 (", round(pca_summary1$importance[2, 1]*100, 1), "%)"),
       y = paste("PC2 (", round(pca_summary1$importance[2, 2]*100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("PCAD" = "red", "non-PCAD" = "blue"))

# Volcano plot for GSE28858
data1 <- allDiff1
data1$significant = "Non_significant"
data1$significant[data1$logFC >= 0.26 & data1$P.Value < 0.05] = "Up"
data1$significant[data1$logFC <= -0.26 & data1$P.Value < 0.05] = "Down"

tiff("volcano_plot_GSE28858.tiff", 
     width = 10, 
     height = 6, 
     units = "in", 
     res = 300,
     compression = "lzw")

ggplot(data1, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = significant), size = 1.5, alpha = 0.3) + 
  theme_minimal() +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept = 1.3, linetype = 4, size = 0.8) +
  geom_vline(xintercept = c(-0.26, 0.26), linetype = 4, size = 0.8) +
  theme(title = element_text(size = 16), 
        text = element_text(size = 16),
        plot.caption = element_text(size = 18, color = "gray40", hjust = 1)) +
  labs(x = "log2(fold_change)", 
       y = "-log10(p_value)",
       caption = "Data source: GSE28858")
dev.off()

# Heatmap for GSE28858
annotation_col1 = data.frame(
  Data = c(rep("GSE28858", 19)),
  Group = c(rep("PCAD", 10), rep("non-PCAD", 9))
)
rownames(annotation_col1) = colnames(exp_anno1)
sorted_index <- order(allDiff1$P.Value, decreasing = FALSE)
sorted_rownames <- rownames(allDiff1)[sorted_index][1:20]
exp_anno1_1 <- exp_anno1[sorted_rownames, ]

tiff("heatmap_GSE28858.tiff", 
     width = 10, 
     height = 6, 
     units = "in", 
     res = 300,
     compression = "lzw")

pheatmap(exp_anno1_1,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col1,
         fontsize = 16,
         show_colnames = FALSE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# Load dataset 2: GSE59421
gset2 <- getGEO(filename = "GSE59421_series_matrix.txt.gz", 
                GSEMatrix = TRUE, 
                getGPL = TRUE)

pdata2 = pData(gset2)
exprSet2 = exprs(gset2)
pdata2 = pdata2 %>% filter(.data[["age:ch1"]] <= 50)
exprSet2 = exprSet2[, colnames(exprSet2) %in% rownames(pdata2)]
pdata2 = fread('pdata_GSE59421_2.CSV')
rownames(pdata2) = pdata2$Sample
exprSet2 = exprSet2[, colnames(exprSet2) %in% rownames(pdata2)]
anno_exprSet2 = na.omit(exprSet2)

group_list2 = factor(pdata2$Group)
group_list2 <- relevel(group_list2, ref = "non-PCAD")
exp_anno2 = normalizeBetweenArrays(anno_exprSet2)

design2 = model.matrix(~ group_list2)
colnames(design2) <- levels(group_list2)
rownames(design2) <- colnames(exp_anno2)

fit2 = lmFit(exp_anno2, design2)
fit2 = eBayes(fit2)
allDiff2 = topTable(fit2, coef = 2, adjust = 'fdr', number = Inf)
write.table(allDiff2, file = "allDiff2.txt", sep = "\t", col.names = NA)

# PCA analysis for GSE59421
exprSet_clean2 <- na.omit(exp_anno2)
exprSet_filtered2 <- exprSet_clean2[rowSums(exprSet_clean2) > 0, ]
gene_variance <- apply(exprSet_filtered2, 1, var, na.rm = TRUE)
nonzero_variance_genes <- gene_variance > 0
exprSet_filtered2 <- exprSet_filtered2[nonzero_variance_genes, ]

pca_result2 <- prcomp(t(exprSet_filtered2), scale. = TRUE, center = TRUE)
pca_summary2 <- summary(pca_result2)
pca_data2 <- as.data.frame(pca_result2$x)
pca_data2$Sample <- colnames(exprSet_filtered2)
pca_data2$Group <- pdata2$Group

ggplot(pca_data2, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  labs(x = paste("PC1 (", round(pca_summary2$importance[2, 1]*100, 1), "%)"),
       y = paste("PC2 (", round(pca_summary2$importance[2, 2]*100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("PCAD" = "red", "non-PCAD" = "blue"))

# Volcano plot for GSE59421
data2 <- allDiff2
data2$significant = "Non_significant"
data2$significant[data2$logFC >= 0.26 & data2$P.Value < 0.05] = "Up"
data2$significant[data2$logFC <= -0.26 & data2$P.Value < 0.05] = "Down"

tiff("volcano_plot_GSE59421.tiff", 
     width = 10, 
     height = 6, 
     units = "in", 
     res = 300,
     compression = "lzw")

ggplot(data2, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = significant), size = 1.5, alpha = 0.3) + 
  theme_minimal() +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept = 1.3, linetype = 4, size = 0.5) +
  geom_vline(xintercept = c(-0.26, 0.26), linetype = 4, size = 0.5) +
  theme(title = element_text(size = 16), 
        text = element_text(size = 16),
        plot.caption = element_text(size = 12, color = "gray40", hjust = 1)) +
  labs(x = "log2(fold_change)", 
       y = "-log10(p_value)",
       caption = "Data source: GSE59421")
dev.off()

# Heatmap for GSE59421
annotation_col2 = data.frame(
  Data = c(rep("GSE59421", 31)),
  Group = pdata2$Group
)
rownames(annotation_col2) = colnames(exp_anno2)
sorted_index2 <- order(allDiff2$P.Value, decreasing = FALSE)
sorted_rownames2 <- rownames(allDiff2)[sorted_index2][1:20]
exp_anno2_2 <- exp_anno2[sorted_rownames2, ]
exp_anno2_2 <- exp_anno2_2[, c(1, 5, 6, 7, 10, 11, 13, 14, 15, 
                               18, 22, 23, 24, 26, 28, 30, 31, 2,
                               3, 4, 8, 9, 12, 16, 17, 19, 20, 21, 25, 27, 29)]

tiff("heatmap_GSE59421.tiff", 
     width = 10, 
     height = 6, 
     units = "in", 
     res = 300,
     compression = "lzw")

pheatmap(exp_anno2_2,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col2,
         fontsize = 16,
         show_colnames = FALSE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# Get miRNA lists from both datasets
data1$significant = "Non_significant"
data1$significant[data1$P.Value < 0.05] = "Significant"
GSE28858_list = data1 %>% 
  filter(significant != 'Non_significant') %>% row.names()
GSE28858_list <- gsub("\\.$", "", GSE28858_list)
GSE28858_list <- gsub("\\.", "-", GSE28858_list)

data2$significant = "Non_significant"
data2$significant[data2$P.Value < 0.05] = "Significant"
GSE59421_list = data2 %>% 
  filter(significant != 'Non_significant') %>% row.names()

inter_miRNAs <- intersect(GSE28858_list, GSE59421_list)
write.table(inter_miRNAs, file = 'inter_miRNAs.tsv', sep = '\t')

# Venn diagram
venn_list <- list(
  GSE28858 = GSE28858_list,
  GSE59421 = GSE59421_list
)

p = ggvenn(
  venn_list,
  fill_color = c("#00AFBB", "#E7B800"),
  stroke_size = 0.8,
  set_name_size = 9,
  text_size = 9
)

ggsave("venn_plot.tiff", p, width = 10, height = 6, dpi = 300)

# Target gene prediction
inter_list = fread('converted_mirnas.txt', header = FALSE)
inter_list = as.data.frame(inter_list)
inter_list = inter_list$V1

library(multiMiR)
target_validated <- get_multimir(
  org = "hsa",
  mirna = inter_list,
  table = "validated",
  summary = TRUE
)

target_df <- as.data.frame(target_validated@data)
target_high_conf <- subset(target_df, support_type == "Functional MTI")
target_high_conf <- unique(target_high_conf[, c("mature_mirna_id", "target_symbol")])
target_genes <- unique(target_high_conf$target_symbol)

write.table(target_genes, "validated_target_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Gene ID conversion
gene_id <- bitr(
  target_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# GO enrichment analysis
ego <- enrichGO(
  gene = gene_id$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE,
)

write.table(as.data.frame(ego), "GO_enrichment_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# KEGG enrichment analysis
ekegg <- enrichKEGG(
  gene = gene_id$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

write.table(as.data.frame(ekegg), "KEGG_enrichment_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot enrichment results
dotplot(ego, showCategory = 12, title = "GO Biological Process (validated targets)") + 
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )
ggsave("GO_dotplot.tiff", width = 10, height = 8, dpi = 300)

dotplot(ekegg, showCategory = 12, title = "KEGG Pathway (validated targets)") + 
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )
ggsave("KEGG_dotplot.tiff", width = 10, height = 8, dpi = 300)

# miRNA-gene interaction network
edges <- target_high_conf
colnames(edges) <- c("from", "to")
fwrite(edges, file = 'target_gene_net.txt', sep = '\t')

net <- graph_from_data_frame(edges, directed = TRUE)

ggraph(net, layout = "fr") + 
  geom_edge_link(alpha = 0.5) + 
  geom_node_point(size = 5, color = "skyblue") +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()
ggsave("miRNAs_gene_interaction.tiff", width = 10, height = 8, dpi = 300)
```