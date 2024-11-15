
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

# Load the BRCA count data
brca_data <- read.csv("count_data_brca.csv", row.names = 1,)
meta <- read.csv("brca_sample_info.csv", row.names = 1)


# Transpose the data so that genes are columns and samples are rows
all.trans <- data.frame(t(brca_data))


# Select the top 1000 most variable genes based on standard deviation
SDs <- apply(all.trans, 2, sd)
topPreds <- order(SDs, decreasing = TRUE)[1:1000]
all.trans <- all.trans[, topPreds]

#merge the gene expression data with metadata
all.trans <- merge(all.trans, meta, by = "row.names")
row.names(all.trans) = as.character(all.trans$Row.names)
all.trans <- all.trans[, 2:ncol(all.trans)]
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
#all.corr <- preProcess(all.trans, method = 'corr', cutoff = 0.1)
#all.trans <- predict(all.corr, all.trans)
dim(all.trans)



#TO WORK ON
# Splitting into training and testing sets (70:30 split)
intrain <- createDataPartition(y = all.trans$tissue_type, p = 0.7) [[1]]
length(intrain)


# Separate training and test sets
train.brca <- all.trans[intrain, ]
test.brca <- all.trans[-intrain, ]

# Check dimensions of training and test sets
dim(train.brca)
dim(test.brca)

#train
#control group
ctrl.brca <- trainControl(method = 'cv', number = 5)


#train
knn.brca <- train(tissue_type~.,
                  data = train.brca,
                  method = 'knn',
                  trControl = ctrl.brca,
                  tuneGrid = data.frame (k=1:20))

#best k
knn.brca$bestTune

#check factor levels:
levels(train.brca)  # replace with your actual training data column name  
levels(test.brca)# replace with your actual test data column name
# Make sure the factor levels in test data match those in the training data
test.brca$ids <- factor(test.brca$ids, levels = levels(train.brca$ids))
# Check which levels in test.brca$ids are not present in train.brca$ids
unseen_levels <- setdiff(levels(test.brca$ids), levels(train.brca$ids))
# Option 1: Remove samples with unseen levels in the test data
test.brca <- test.brca[!test.brca$ids %in% unseen_levels, ]
# Option 2: Adjust levels so that unseen levels are ignored
test.brca$ids <- factor(test.brca$ids, levels = levels(train.brca$ids))
# Exclude rows with NA values in predictions (after ensuring levels match)
testPred <- predict(knn.brca, newdata = test.brca)
testPred <- na.omit(testPred)


#predict
trainPred <- predict(knn.brca, newdata = train.brca)
testPred <- predict(knn.brca, newdata = test.brca)


#interpretation
#confusion matrix
#convert to factors
train.brca$tissue_type <- factor(train.brca$tissue_type)
test.brca$tissue_type <- factor(test.brca$tissue_type)
levels(testPred)  
levels(train.brca$tissue_type)
#check length
length(testPred)  
length(test.brca$tissue_type)
# Check for NA values in the predictions
sum(is.na(testPred))  # Count the number of NA predictions
# Convert to factor if necessary
test.brca$tissue_type <- factor(test.brca$tissue_type, levels = levels(train.brca$tissue_type))

# Ensure test set has only the same columns as the training set
test.brca <- test.brca[, colnames(train.brca)]

testPred <- predict(knn.brca, newdata = test.brca)


#confusion matrix
library(caret)  
confusionMatrix(trainPred, train.brca$tissue_type)
confusionMatrix(testPred, test.brca$tissue_type)



#determine variable importance
explainer.brca <- explain(knn.brca,
                          data = train.brca,
                          label = 'knn',
                          y = as.numeric(train.brca$tissue_type),
importance.brca <- feature_importance(explainer.brca, n-sample = 40, type = 'difference')

head(importance.brca$variable)
tail(importance.brca$variable)
plot(importance.brca)





