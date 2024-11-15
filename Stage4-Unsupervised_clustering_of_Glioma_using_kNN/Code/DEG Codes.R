library(TCGAbiolinks)
library(DESeq2)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(curatedTCGAData)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
library(gplots)

# Query and download data from TCGA
query <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)

GDCdownload(query)
data <- GDCprepare(query)

table(data@colData@listData[["definition"]])
keep_samples <- data@colData@listData[["definition"]] != "Recurrent Solid Tumor"
data_516 <- data[, keep_samples]
table(data_516@colData@listData[["definition"]])
# Extract gene expression matrix and clinical data
expr_matrix <- assay(data_516)
coldata <- as.data.frame(data_516@colData@listData[["paper_IDH.status"]])
rownames(coldata) <- data_516@colData@listData$barcode
names(coldata)[1] <-"idh_status"

######## 5 convert Ensemble Ids to gene symbols ########

# Install and load org.Hs.eg.db package
library(org.Hs.eg.db)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(expr_matrix)

# Remove version numbers for compatibility
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove NA values and keep only mapped genes
mapped_genes <- !is.na(gene_symbols)
expr_matrix <- expr_matrix[mapped_genes, ]
gene_symbols <- gene_symbols[mapped_genes]

unique_gene_symbols <- make.unique(gene_symbols, sep = "_")
length(unique_gene_symbols)

# Assign gene symbols to prostate_matrix
rownames(expr_matrix) <- unique_gene_symbols
coldata <- na.omit(coldata)

expr_matrix <- expr_matrix[,rownames(coldata)]

# making the rownames and column names identical
all(rownames(coldata) %in% colnames(expr_matrix))
all(rownames(coldata) == colnames(expr_matrix))

expr_matrix <- round(expr_matrix)


# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = expr_matrix,
                              colData = coldata,
                              design = ~ idh_status)
# Filter out low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set the reference level for the factor of interest
dds$idh_status <- relevel(dds$idh_status, ref = "WT")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Explore results
summary(res)

# Filter significant genes
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 2)

summary(sig_genes)
# Visualize results
plotMA(res, ylim=c(-5,5))










