# Load required packages
install.packages("caret")
install.packages("DALEX")
install.packages("pROC")
install.packages("ggplot2")
install.packages("lattice")
install.packages("ranger")

library(caret)
library(DALEX)
library(pROC)
library(ggplot2)
library(lattice)
library(ranger)

# Set a seed for reproducibility
set.seed(34567)

# Load the BRCA count data
brca_data <- read.csv("count_data_brca.csv", row.names = 1)
meta <- read.csv("brca_sample_info.csv", row.names = 1) 

# Ensure that the sample IDs match between count data and metadata
brca_data <- brca_data[, rownames(meta)]

#samples under each subset
table(meta$tissue_type)

#preview normalized data
#boxplot(brca_data[,1:100], las=2)
#par(oma = c(10,0,0,0)) #converts data to log10 or easy visibility
#boxplot(log10(brca_data[1:100] +1), ylim = c(0,10), las=2)

#preprocessing1
# Transpose the data so that genes are columns and samples are rows
all.trans <- data.frame(t(brca_data))


# Select the top 1000 most variable genes based on standard deviation
SDs <- apply(all.trans, 2, sd)
topPreds <- order(SDs, decreasing = TRUE)[1:1000]
all.trans <- all.trans[, topPreds]

#merge the gene expression data with metadata
all.trans <- merge(all.trans, meta, by = "row.names")
dim(all.trans)
head(all.trans[, 1:5])

# PreProcessing steps
# 1. Remove near-zero variance predictors
all.zero <- preProcess(all.trans, method = 'nzv', uniqueCut = 15)
all.trans <- predict(all.zero, all.trans)

# 2. Center the data
all.center <- preProcess(all.trans, method = 'center')
all.trans <- predict(all.center, all.trans)

# 3. Remove highly correlated features
all.corr <- preProcess(all.trans, method = 'corr', cutoff = 0.5)
all.trans <- predict(all.corr, all.trans)


# Split the data into training and testing sets
intrain <- createDataPartition(y = all.trans$tissue_type, p = 0.8) [[1]]

# Create the training and testing datasets
train.brca <- all.trans[intrain, ]
test.brca <- all.trans[-intrain, ]
#check length
dim(train.brca)
dim(test.brca)

# Random forest model training
rf.ctrl <- trainControl(method = 'cv', number = 10)  # Use 10-fold cross-validation(avoid bias)

rf.brca <- train(tissue_type ~ .,               # Model formula using your target variable
                 data = train.brca,                # Training data
                 method = 'ranger',                # Random forest method from the ranger package
                 trControl = rf.ctrl,              # Corrected 'trControl' argument
                 importance = 'permutation',       # Permutation-based importance
                 tuneGrid = data.frame(mtry = 100,
                                       min.node.size = 1,
                                       splitrule = 'gini'))

# Error rate
rf.brca$finalModel$prediction.error
plot(varImp(rf.brca), top = 10)

