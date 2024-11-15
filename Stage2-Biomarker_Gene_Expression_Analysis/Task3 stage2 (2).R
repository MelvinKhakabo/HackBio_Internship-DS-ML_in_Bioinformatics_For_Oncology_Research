# loading the data
glioblastoma_data <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv", row.names= 1) 
# This argument uses the first column of the dataset (which contains gene names) as row names. 
# This ensures that only the remaining columns (which contain the expression values) are treated as numeric data.

View(glioblastoma_data)


library(gplots)                              

# generating heatmaps and clustering samples to look at the pattern of gene expression and figure out how to group them
heatmap.2(as.matrix(glioblastoma_data), 
          trace = 'none', 
          scale='row', # for each row, the values are centered (mean = 0) and scaled (standard deviation = 1).
          dendrogram = 'col' , 
          Rowv = FALSE , 
          Colv = TRUE # Cluster samples 
          )

# using diverging color palette 
heatmap.2(as.matrix(glioblastoma_data),
          main = "Diverging Heatmap of Gene Expression" , 
          trace = 'none', 
          scale='row', 
          dendrogram = 'col' , 
          Rowv = FALSE , 
          Colv = TRUE , 
          col = colorRampPalette(c("blue", "white", "red"))(100)) 

# using sequential color palette 
heatmap.2(as.matrix(glioblastoma_data), 
           main = "Sequential heatmap of Gene Expression" , 
           trace = 'none', 
           scale='row', 
           dendrogram = 'col' , 
           Rowv = FALSE , 
           Colv = TRUE ,
          col = colorRampPalette(c("white", "darkred"))(100))

# variants of the heatmap: 1- clustered rows
heatmap.2(as.matrix(glioblastoma_data), 
          main = "Clustered rows (genes)" , 
          trace = 'none', 
          scale='row', 
          dendrogram = 'row' , 
          Rowv = TRUE , 
          Colv = FALSE , 
          col = colorRampPalette(c("blue", "white", "red"))(100)) 
# 2- clustered columns
heatmap.2(as.matrix(glioblastoma_data), 
          main = "Clustered columns (samples)" , 
          trace = 'none', 
          scale='row', 
          dendrogram = 'col' , 
          Rowv = FALSE , 
          Colv = TRUE , 
          col = colorRampPalette(c("blue", "white", "red"))(100)) 

# clustering both column and rows
heatmap.2(as.matrix(glioblastoma_data), 
          main = "Clustered rows & columns" , 
          trace = 'none', 
          scale='row', 
          dendrogram = 'both' , 
          Rowv = TRUE , 
          Colv = TRUE , 
          col = colorRampPalette(c("blue", "white", "red"))(100)) 


# showing samples (columns) names 
colnames(glioblastoma_data)

# manually grouping the samples with the help of the heatmap
# the first 5 groups are simliar and the last 5 groups are similar
samples1 <- c(1,2,3,4,5)
samples2 <- c(6,7,8,9,10)

samples1_data <- glioblastoma_data[, samples1]
samples2_data <- glioblastoma_data[, samples2]

head(samples1_data)
head(samples2_data)

# calculating means 
samples1_mean <- rowMeans(samples1_data)
samples2_mean <- rowMeans(samples2_data)

head(samples1_mean)
head(samples2_mean)

# calculating fold change 
fold_change <- samples2_mean / samples1_mean
fold_change

# calculating log2 fold change
log2_fold_change <- log2(fold_change)
log2_fold_change

# calculating p-value
p_values <- apply(glioblastoma_data, 1, function(gene_expression_row) {
  t.test(gene_expression_row[1:5], gene_expression_row[6:10])$p.value })
# t.test (gene_expression_row[]) only accept vectors

# visualize the log2 fold change and negative log of p-values to see distribution 
plot(log2_fold_change, -log10(p_values))


# sub-setting up-regulated and down-regulated genes

View(glioblastoma_data)
upregulated_genes <- subset(glioblastoma_data, log2_fold_change > 1 & p_values < 2)
downregulated_genes <- subset(glioblastoma_data, log2_fold_change < -1 & p_values < 2)

head(upregulated_genes)
head(downregulated_genes)

# Number of upregulated genes
nrow(upregulated_genes)

# Number of downregulated genes
nrow(downregulated_genes)

# Save upregulated and downregulated genes as CSV for further functional enrichment analysis
write.csv(upregulated_genes, "upregulated_genes.csv", row.names = TRUE)
write.csv(downregulated_genes, "downregulated_genes.csv", row.names = TRUE)

library(ggplot2)

# upregulated pathways plotting
upregulated_enrichments <- read.csv("~/Desktop/enrichment/up.csv")

View(upregulated_enrichments)

# Calculate -log10 of the FDR for significance
upregulated_enrichments$log_FDR <- -log10(upregulated_enrichments$Enrichment.FDR)

# Select the top 5 pathways based on FDR 
top_pathways_up <- upregulated_enrichments[order(upregulated_enrichments$Enrichment.FDR), ][1:5, ]

View(top_pathways_up)

# Create the lollipop plot
ggplot(top_pathways_up, aes(x = reorder(Pathway, nGenes), y = nGenes)) +
  geom_segment(aes(xend = Pathway, yend = 0), color = "gray", linewidth = 1) +
  geom_point(aes(size = log_FDR), color = "red") + 
  coord_flip() +
  labs(x = "Pathway", y = "Number of Genes", 
       title = "Top 5 Pathways of Upregulated Genes") +
  scale_size_continuous(name = "-log10(FDR)", breaks = c(1, 1.5, 2, 2.5, 3) ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

# downregulated pathways plotting
downregulated_enrichment <- read.csv("~/Desktop/enrichment/down.csv")

View(downregulated_enrichment)

# Calculate -log10 of the FDR for significance
downregulated_enrichment$log_FDR <- -log10(downregulated_enrichment$Enrichment.FDR)

# Select the top 5 pathways based on FDR 
top_pathways_down <- downregulated_enrichment[order(downregulated_enrichment$Enrichment.FDR), ][1:5, ]

View(top_pathways_down)

# Create the lollipop plot
ggplot(top_pathways_down, aes(x = reorder(Pathway, nGenes), y = nGenes)) +
  geom_segment(aes(xend = Pathway, yend = 0), color = "gray", linewidth = 1) +
  geom_point(aes(size = log_FDR), color = "red") +  
  coord_flip() +
  labs(x = "Pathway", y = "Number of Genes", 
       title = "Top 5 Pathways of Downregulated genes") +
  scale_size_continuous(name = "-log10(FDR)", breaks = c(1, 1.5, 2, 2.5, 3)) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

