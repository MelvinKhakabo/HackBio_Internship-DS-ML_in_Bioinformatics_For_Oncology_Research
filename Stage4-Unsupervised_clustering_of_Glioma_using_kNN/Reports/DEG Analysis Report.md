Contributors: Ajisegiri Sewedo B, Balqees Mansour, Salaam Ridwan, Adekoya Adepeju, Melvin Khakabo, Elizabeth Abodunrin, Otuekong Emmanuel, Samuel Aladegbaiye

### Differential Expression Analysis of Newer Lower-Grade Gliomas Dataset Using Bioinformatics Approach 

**Introduction**

Lower-grade Gliomas (LGG) are brain tumors that occur in children and young adults, often with benign or low-grade behavior (Greuter _et al.,_ 2021). They account for about 10% of all brain tumors and can be associated with genetic mutations, such as IDH1, IDH2, or CHK2 (Lin _et al.,_ 2020). Isocitrate dehydrogenase (IDH) mutations, in particular, have been associated with a favorable prognosis of LGG and exhibit distinct profiles (Greuter _et al.,_ 2021).&#x20;

Previous studies on LGG have characterized the TCGA dataset for IDH status classes and identification of relevant pathways (Ceccarelli _et al.,_ 2016).  However, given the heterogeneous nature of tumors and improved technologies, it becomes crucial to analyze the recent LGG dataset for any disparities. This study uses differential expression analysis to assess the newer TCGA-LGG dataset and compare it with previous findings.

**Method**

This study analyzed the TCGA-LGG dataset retrieved from the TCGA database. The GDCprepare function was used to preprocess the data. A DESeqDataSet object was created from the count matrix and sample information.  A model design based on IDH type (Wild-type vs. Mutant) was specified. Genes with low counts across samples were filtered out to improve statistical power. The reference level for tissue type was set to "WT" to facilitate comparisons against this baseline before the execution of the DESeq function.

**Results**

A heatmap was used to visualize the cluster of IDH statuses, including WT and Mutant, combined (Fig. 1). Meanwhile, the top 20 upregulated genes and the top 20 downregulated genes are shown using the boxplots (Fig. 2-3)

![LGG_gene_expression_heatmap](https://github.com/user-attachments/assets/dc697e3c-3816-4773-a099-2b63b417fd41)
**Figure 1: A heatmap showing the gene expression profile of IDH statuses.**

![top_20_upregulated_genes_LGG](https://github.com/user-attachments/assets/7178161d-107f-418c-bac5-bb112cdc51f8)
**Figure 2: A boxplot showing the top 20 upregulated genes from LGG dataset.**

![top_20_downregulated_genes_LGG](https://github.com/user-attachments/assets/8b7d46c0-d66f-433d-964d-6a4abff4c2b2)
**Figure 3: A boxplot showing the top 20 downregulated genes from the LGG dataset.**

**Discussion**

Based on DEG analysis, our findings align with the outcome of Ceccarelli _et al._ (2016), which identified WT and mutant IDH classes. Also, this study identified ADAMS20, ANKRD62P1, and DEFB119 as the top three upregulated genes, while AGR3, APCDD1L, and C6ORF15 are the top three downregulated genes. These differ from genes such as mutated TP53, amplified EGFR, and mutated BRAF in the previous study.

**Conclusions**

This study showed that molecular profiling and characterization of LGG can benefit from an extensive differential expression analysis of newer datasets. This can help identify previously unknown biomarkers.

**References**

Ceccarelli, M., Barthel, Floris P., Malta, Tathiane M., Sabedot, Thais S., Salama, Sofie R., Murray, Bradley A., Morozova, O., Newton, Y., Radenbaugh, A., Pagnotta, Stefano M., Anjum, S., Wang, J., Manyam, G., Zoppoli, P., Ling, S., Rao, Arjun A., Grifford, M., Cherniack, Andrew D., Zhang, H., & Poisson, L. (2016). Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. _Cell_, _164_(3), 550–563. https\://doi.org/10.1016/j.cell.2015.12.028

Greuter, L., Guzman, R., & Soleman, J. (2021). Pediatric and Adult Low-Grade Gliomas: Where Do the Differences Lie? _Children_, _8_(11), 1075. https\://doi.org/10.3390/children8111075

Lin, L., Cai, J., Tan, Z., Meng, X., Li, R., Li, Y., & Jiang, C. (2020). Mutant IDH1 Enhances Temozolomide Sensitivity via Regulation of the ATM/CHK2 Pathway in Glioma. _Cancer Research and Treatment_. https\://doi.org/10.4143/crt.2020.506

 

 
