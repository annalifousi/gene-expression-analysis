# Install required packages
install.packages("readxl")
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ComplexHeatmap", "ggplot2", "msigdbr", "pheatmap", "matrixStats", "fgsea"))

# Load libraries
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(matrixStats)

# Load data
tpm <- read_excel("GSE281695_TPM_values_all_samples.xlsx")
normalized_counts <- read_excel("GSE281695_normalized_counts_all_samples.xlsx")
raw_counts <- read_excel("GSE281695_raw_counts_all_samples.xlsx")
rlog_counts <- read_excel("GSE281695_rlog_transformed_counts_all_samples.xlsx")

# Rename first column to "Gene" if it contains gene IDs
colnames(raw_counts)[1] <- "Gene"
raw_counts <- column_to_rownames(raw_counts, "Gene")

# Remove constant rows before PCA analysis
constant_cols <- apply(raw_counts, 1, var) == 0
raw_counts_pca <- raw_counts[!constant_cols, ]

# PCA plot
pca_count <- prcomp(t(raw_counts_pca), scale. = TRUE)
df <- data.frame(PC1=pca_count$x[,1], PC2=pca_count$x[,2])
df$Group <- colnames(raw_counts)
ggplot(df, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point() + 
  labs(x = "PC1", y = "PC2", title = "PCA Plot") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Create metadata table for DESeq2 analysis
sample_metadata <- data.frame(
  sample_id = colnames(raw_counts),
  condition = c("IR_24h", "NonIR_24h", "IR_7d", "NonIR_7d", "IR_24h", "NonIR_7d")
)
rownames(sample_metadata) <- sample_metadata$sample_id

# DESeq2 dataset and analysis
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = sample_metadata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Volcano plot of DESeq2 results
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene") %>%
  mutate(significant = padj < 0.03)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) + 
  xlim(c(-3, 3)) + 
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot") +
  theme_minimal()

# DESeq2 comparison function
perform_DESeq2 <- function(dds, condition1, condition2) {
  res <- results(dds, contrast=c("condition", condition1, condition2))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("Gene") %>%
    mutate(significant = padj < 0.05, condition = ifelse(log2FoldChange > 0, condition1, condition2))
  res_df <- na.omit(res_df)
  ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant, shape=condition)) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values = c("turquoise", "red")) +
    scale_shape_manual(values = c(16, 17)) +
    labs(title = paste("Volcano Plot:", condition1, "vs", condition2),
         x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
    theme_minimal() +
    theme(legend.position = "top")
}

# Plot volcano plots for each condition pair
plot1 <- perform_DESeq2(dds, "IR_24h", "NonIR_24h")
plot2 <- perform_DESeq2(dds, "IR_7d", "NonIR_7d")
plot3 <- perform_DESeq2(dds, "IR_24h", "IR_7d")
plot4 <- perform_DESeq2(dds, "NonIR_24h", "NonIR_7d")

# Calculate variance-stabilized data
vsd <- vst(dds)
row_variances <- rowVars(assay(vsd))

# Generate heatmap of top 100 most variable genes
top_genes_indices <- order(row_variances, decreasing = TRUE)[1:100]
mat <- assay(vsd)[top_genes_indices, ]
annotation_df <- sample_metadata[colnames(mat), , drop = FALSE]
pheatmap(mat, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", clustering_method = "complete", 
         annotation_col = annotation_df, show_rownames = FALSE, show_colnames = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(50), main = "Top 100 Most Variable Genes")

# Load gene sets for GSEA
mouse_gene_sets <- msigdbr(species = "Mus musculus", category = "C5")
gene_sets_list <- split(mouse_gene_sets$ensembl_gene, mouse_gene_sets$gs_name)

# Rank genes for GSEA analysis
genes_ranked <- res$log2FoldChange / res$pvalue
names(genes_ranked) <- rownames(res)
genes_ranked <- genes_ranked[!is.na(genes_ranked)]

# Perform GSEA
fgsea_results <- fgsea(pathways = gene_sets_list, stats = genes_ranked, minSize = 15, maxSize = 500)

# Plot top 10 enriched pathways
topPathways <- fgsea_results %>%
  dplyr::arrange(padj) %>%
  dplyr::slice(1:10)
ggplot(topPathways, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
  geom_col(show.legend = FALSE) +
  scale_fill_gradient(low = "blue", high = "lightblue") +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", title = "Top 10 Enriched Pathways") +
  theme_minimal()


