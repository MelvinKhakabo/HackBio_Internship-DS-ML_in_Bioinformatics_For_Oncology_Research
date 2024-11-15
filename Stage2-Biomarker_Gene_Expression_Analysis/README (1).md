# stage-2-hackbio
# Gene Expression Analysis with Heatmaps and Pathway Enrichment

This repository contains R code for analyzing gene expression data, including generating heatmaps and identifying upregulated and downregulated genes. 

## Requirements

- R (version 4.0 or later)
- Required R packages: gplots, ggplot2, RColorBrewer

## Installation

Install the necessary R packages using the following commands:

install.packages(c("gplots", "RColorBrewer", "ggplot2"))

## Usage

1. *Load the Data*
   - The script loads gene expression data from a provided URL. Ensure you have an active internet connection or download the data manually if needed.

2. *Generate Heatmaps*
   - The script generates various heatmaps to visualize gene expression patterns with different color palettes and clustering options.

3. *Analyze Gene Expression*
   - The script compares gene expression between two groups of samples. It calculates fold changes, log2 fold changes, and p-values for each gene to determine differential expression.

4. *Identify Differentially Expressed Genes*
   - Upregulated and downregulated genes are identified based on fold change and p-values. These genes are saved to CSV files for further analysis.

5. *Visualize Pathway Enrichment*
   - The script creates lollipop plots to visualize the top pathways associated with upregulated and downregulated genes using pre-processed enrichment data.

## Files

- glioblastoma.csv: Gene expression data used for the analysis.
- upregulated_genes.csv: CSV file containing identified upregulated genes.
- downregulated_genes.csv: CSV file containing identified downregulated genes.
- up.csv: Enrichment data for upregulated pathways.
- down.csv: Enrichment data for downregulated pathways.

## Instructions

1. place the script in your working directory.
2. Run the R script.
