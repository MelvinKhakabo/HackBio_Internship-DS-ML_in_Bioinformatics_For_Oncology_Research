**Authors:** Logy Khaled (@Logy), Alaa Hewela (@Alaa253), Rahma mamdouh Mohammed (@Rahmamam2000), Melvin Khakabo (@MELmostly), Tawfek Ahmed Tawfek (@Tawfekahmed25)

## **GithubRepo:**

https://github.com/Tawfek25/Stage2-Biomarker-HackBio

## **Goal**

This analysis aimed to visualize differential gene expression in a glioblastoma dataset using heatmaps and interpreting functional enrichment analysis results to identify the top biological pathways associated with this cancer.

.  
## **Visualization of the data**

At first, we made a heatmap in R using the function (heatmap.2) of the package (gplots) that had clustered columns, since the columns represent the samples, this was done to identify the pattern of gene expression and know how to group the samples.

## **Heatmap with diverging color palette**   
![image](https://github.com/user-attachments/assets/a24b610b-63a2-4283-a35c-99aeb1600c87))

**Diverging color palettes** highlight extremes of expression using two colors, allowing us to identify samples with similar gene expression patterns for grouping. This aids in calculating fold change and p-values. The heatmap indicates which genes are upregulated or downregulated in each sample, with red representing highly expressed genes.

## **Heatmap using sequential color palette**  
![image](https://github.com/user-attachments/assets/0d6e36ad-5bfb-4bbb-836a-202e16bef2aa)

**Sequential color palette** visualizes continuous value gradients, ideal for representing data from low to high without a clear midpoint. This aids in understanding the magnitude of changes in gene expression.

**Using the diverging color palette is more significant in visualizing this dataset.**

## **Variants of the heatmap:**

### **With clustering of rows (genes)**

![image](https://github.com/user-attachments/assets/f789f09f-2da4-41bb-bbc3-6011618af00f)

### **With clustering of columns (samples)**

![image](https://github.com/user-attachments/assets/34643c8f-ad90-470a-808d-50fc5efe1a1e)

### **With clustering of rows and columns**

![image](https://github.com/user-attachments/assets/84ea2ddb-4557-4773-851a-bb2b5de8abc4)

## **Functional enrichment analysis**

After grouping the samples, we calculated the fold change and p-values, then plotted them to determine the cutoffs.

![image](https://github.com/user-attachments/assets/6d624ff7-6472-4a24-9422-f6c36adc5f33)

Since there is a clear distinction between positive and negative, we set the log2 fold change cutoff to 1 and the p-value cutoff to 2 to include all genes for better results.  
**Upregulated genes** have **a positive log2 fold change** and **a significant p-value**.

**Downregulated genes** have **a negative log2 fold change** and **a significant p-value**.  
After subsetting the upregulated and downregulated genes, we ran functional enrichment analysis using ShinyGo to identify the top pathways involved.

### **The results were then plotted as lollipop plot**     
![image](https://github.com/user-attachments/assets/877cb6a5-5a8e-4508-896e-df31d062f97a)

## **The top 5 enriched pathways of upregulated genes includes:**   

Cytokine-cytokine receptor interaction, IL-17 signaling pathway, Viral protein interaction with cytokine and cytokine receptor, JAK-STAT signaling pathway, Malaria

Cytokine signaling is crucial for immune regulation and tumor progression. The IL-17 pathway promotes inflammation and immune cell recruitment, contributing to tumor growth and poor prognosis. The viral protein interaction pathway suggests potential viral influences that may alter immune responses. Dysregulation of these pathways supports a tumor microenvironment that facilitates immune evasion and progression in glioblastoma, highlighting the complex interplay between immune signaling and tumor biology.
