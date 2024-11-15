Overview:

In this task, my team and I reproduceD the unsupervised clustering of gliomas (n = 516 LGG) based on methylation levels to classify the two IDH statuses(mutant and wild type) as reported in a specified paper.
 
 I (the data scientist) employed the K-Nearest Neighbors (KNN) method for the classification. The goal being to demonstrate our ability to perform unsupervised clustering using methylation data and to evaluate the clustering performance. 
 
 The biomarker hunters developed expression biomarkers that can distinguish any 2 IDH classes.

Steps
Read the paper: https://www.cell.com/cell/fulltext/S0092-8674(15)01692-X 
Carefully read the paper that reports the clustering of gliomas using methylation levels to classify the six IDH statuses.
Understaood the methodology and the specific findings related to unsupervised clustering and the use of methylation data.
Obtained the datasets from the TCGA (more up to date) using TCGABiolinks in R.
Ensured the identitiy of the IDH status for each sample.
Data was Preprocessed: Cleaned and preprocessed according to our pipeline and needs.
For ML engineers, feature selection wass important: For ML, perform kNN and random forest
For Biomarker-Devs, normalization was important: performed DEA

Generated final report containing:
Introduction to gliomas, IDH status, and the significance of methylation levels.
Description of the dataset and preprocessing steps.
Methodology for implementing and applying the KNN algorithm.
Results of the kNN, including visualizations and evaluation metrics as presented in the paper.
Comparison with the findings reported in the target paper. (Do you think we need more clusters based on newer datasets)
Conclusion and insights gained from the analysis.
