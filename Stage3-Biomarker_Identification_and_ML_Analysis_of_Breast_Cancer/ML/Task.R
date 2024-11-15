# Load required packages
install.packages("caret")
install.packages("DALEX")
install.packages("pROC")
install.packages("ggplot2")
install.packages("lattice")

library(caret)
library(DALEX)
library(pROC)
library(ggplot2)
library(lattice)

# Set a seed for reproducibility
set.seed(34567)

# Load the BRCA data
brca_data <- read.csv("count_data_brca.csv", row.names = 1)
meta <- read.csv("brca_sample_info.csv", row.names = 1)

# Ensure that the sample IDs match between count data and metadata
brca_data <- brca_data[, rownames(meta)]

# Check the distribution of tissue types (tumor vs normal)
table(meta$tissue_type)

# Transpose the data so that genes are columns and samples are rows
all.trans <- data.frame(t(brca_data))

# Select the top 1000 most variable genes based on standard deviation
SDs <- apply(all.trans, 2, sd)
topPreds <- order(SDs, decreasing = TRUE)[1:1000]
all.trans <- all.trans[, topPreds]
dim(all.trans)

# Merge the gene expression data with metadata (including tissue_type labels)
all.trans <- merge(all.trans, meta, by = "row.names")
dim(all.trans)
head(all.trans[, 1:5])

# Set row names to sample IDs and remove the redundant column
rownames(all.trans) <- all.trans$Row.names
all.trans <- all.trans[,-1]
head(all.trans[, 1:5])

# Preprocessing steps
# 1. Remove near-zero variance predictors
all.zero <- preProcess(all.trans, method = 'nzv', uniqueCut = 15)
all.trans <- predict(all.zero, all.trans)

# 2. Center the data
all.center <- preProcess(all.trans, method = 'center')
all.trans <- predict(all.center, all.trans)

# 3. Remove highly correlated features
all.corr <- preProcess(all.trans, method = 'corr', cutoff = 0.5)
all.trans <- predict(all.corr, all.trans)

dim(all.trans)

# Splitting into training and testing sets (70:30 split)
intrain <- createDataPartition(y = all.trans$tissue_type, p = 0.7, list = FALSE)
length(intrain)

# Separate training and test sets
train.brca <- all.trans[intrain, ]
test.brca <- all.trans[-intrain, ]

# Check dimensions of training and test sets
dim(train.brca)
dim(test.brca)

#predict
trainPred <- predict(knn.brca, newdata = train.brca)
testPred <- predict(knn.brca, newdata = test.brca)

#interpretation
#confusion matrix
confusionMatrix(trainPred, train.brca$tissue_type)
confusionMatrix(testPred, train.brca$tissue_type)
