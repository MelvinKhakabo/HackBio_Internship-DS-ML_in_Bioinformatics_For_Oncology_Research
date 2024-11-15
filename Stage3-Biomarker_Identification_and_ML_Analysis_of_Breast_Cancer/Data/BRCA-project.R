install.packages("gplots")
BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
library(ggplot2)
library(reshape2)
install.packages("ggrepel")
library(ggrepel)
# load the required libraries
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(curatedTCGAData)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
library(gplots)

# read the data 
raw.data <- read.csv("/home/balqees/Downloads/rcounts.tsv", sep = '\t')
rownames(raw.data) <- raw.data$gene_id
raw.data <- raw.data[,-1]
data.subset <- raw.data[,sample_info$ids]
sample_info <- read.csv("/home/balqees/Desktop/Hackbio-Internship/Stage3-BRCA/brca_sample_info.csv")

rownames(sample_info) <- sample_info$X
sample_info <- sample_info[,-1]

# making the rownames and column names identical
all(rownames(sample_info) %in% colnames(data.subset))
all(rownames(sample_info) == colnames(data.subset))
## DESeq2 Analysis

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(data.subset),
                              colData = sample_info,
                              design = ~ tissue_type)

# Pre-filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set factor levels
dds$tissue_type <- relevel(dds$tissue_type, ref = "normal")

# Run DESeq
dds <- DESeq(dds)
res <- results(dds)

# Explore Results ----------------
summary(res)
res0.01 <- results(dds, alpha = 0.05)

# Filter for log2 fold change > 1 or < -1
res_filtered <- subset(res0.01, log2FoldChange > 3 | log2FoldChange < -3)
summary(res_filtered)

# volcano visualization
# Convert DESeq2 results to a data frame
res_df <- as.data.frame(res_filtered)

# Add a column to identify significantly differentially expressed genes
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                             "Significant", "Not Significant")

# Create the volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value")

# Optional: Add labels for top genes
top_genes <- res_df[order(res_df$padj)[1:10], ]
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes))) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value")

# Save the plot as a PNG file
ggsave("enhanced_volcano.png", plot = p, width = 10, height = 8, dpi = 300)

#==============================================================================
#-----------01# BoxPLot for the 20 upregulated genes in tumor vs normal ----------
#------------------------------------------------------------------------------
#-=============================================================================
library(DESeq2)
library(tidyverse)
library(ggplot2)

# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for upregulated genes (log2FoldChange > 1 and adjusted p-value < 0.05)
res_upregulated <- res_filtered[which(res_filtered$log2FoldChange > 3 & res_filtered$padj < 0.05),]

# Order by log2FoldChange
res_upregulated <- res_upregulated[order(res_upregulated$log2FoldChange, decreasing = TRUE),]

# Select top 18 upregulated genes
top_upregulated <- head(res_upregulated, 20)

# Extract normalized counts for these genes
normalized_counts <- counts(dds, normalized=TRUE)[rownames(top_upregulated),]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(sample_info, by = c("sample" = "ids"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(res_filtered) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = tissue_type, y = log2(count + 1), fill = tissue_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("normal" = "lightgray", "tumor" = "purple")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = tissue_type[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)


# Save the plot
ggsave("top_20_upregulated_genes.png", width = 15, height = 10, dpi = 300)

write.csv(res_upregulated, "top_20_upregulated_genes.csv")
#==============================================================================
#----- #02 BoxPLot for the top 20  downregulated genes in tumor vs normal -----
#------------------------------------------------------------------------------
#-=============================================================================
library(DESeq2)
library(tidyverse)
library(ggplot2)

# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for upregulated genes (log2FoldChange > 1 and adjusted p-value < 0.05)
res_downregulated <- res_filtered[which(res_filtered$log2FoldChange < -3 & res_filtered$padj < 0.05),]

# Order by log2FoldChange
res_downregulated <- res_downregulated[order(res_downregulated$log2FoldChange, decreasing = FALSE),]

# Select top 18 upregulated genes
top_downregulated <- head(res_downregulated, 20)

# Extract normalized counts for these genes
normalized_counts <- counts(dds, normalized=TRUE)[rownames(top_downregulated),]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(sample_info, by = c("sample" = "ids"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(res_filtered) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = tissue_type, y = log2(count + 1), fill = tissue_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("normal" = "lightpink", "tumor" = "purple")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = tissue_type[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)


# Save the plot
ggsave("top_20_downregulated_genes.png", width = 15, height = 10, dpi = 300)
write.csv(res_downregulated,"Downregulated_20genes.csv")

library(clusterProfiler)
library(org.Hs.eg.db)
# Convert gene symbols to Entrez IDs
gene_symbols <- rownames(res_filtered)
gene_list <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Filter out any NA values that may result from conversion
gene_list <- gene_list[!is.na(gene_list$ENTREZID), "ENTREZID"]

ekegg <- enrichKEGG(gene         = gene_list,
                    organism     = 'hsa', # Human
                    pvalueCutoff = 0.05)
barplot(ekegg, showCategory = 15)


ekegg_df <- as.data.frame(ekegg@result)

# Assuming 'Count' and 'GeneRatio' are columns in your data frame
ekegg_df$RichFactor <- ekegg_df$Count / as.numeric(sapply(strsplit(ekegg_df$GeneRatio, "/"), `[`, 2))

# Create a lollipop plot using the calculated RichFactor
ggplot(ekegg_df[1:20, ], aes(x = RichFactor, y = fct_reorder(Description, RichFactor))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10") +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  xlab("Rich Factor") +
  ylab(NULL) +
  ggtitle("KEGG Enrichment Lollipop Plot")
ggsave("lollipop_KEGG20_pathways.png", width = 15, height = 10, dpi = 300)







