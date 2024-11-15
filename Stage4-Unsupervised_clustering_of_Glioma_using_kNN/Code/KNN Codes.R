# Load required packages
install.packages("caret")
install.packages("DALEX")
install.packages("pROC")
install.packages("ggplot2")
install.packages("lattice")
install.packages("iml")

library(caret)
library(DALEX)
library(pROC)
library(ggplot2)
library(lattice)
library(iml)

# Set a seed for reproducibility
set.seed(34567)

# Load the LGG count data

lgg_data <- read.csv("glioma_rawcounts.csv", row.names = 1,)
meta <- read.csv("glioma_metadata.csv", row.names = 1)


dim(lgg_data)
rownames(lgg_data) <- lgg_data$X

lgg_data$X <- NULL

colnames(lgg_data) <- gsub("\\.", "-", colnames(lgg_data))
lgg_data <- log10(lgg_data + 1)


# Transpose the data so that genes are columns and samples are rows
all.trans <- data.frame(t(lgg_data))


# Select the top 1000 most variable genes based on standard deviation
SDs <- apply(all.trans, 2, sd)
topPreds <- order(SDs, decreasing = TRUE)[1:1000]
topPreds <- as.character(topPreds)


rownames(meta) <- meta$barcode

meta$barcode <- NULL

#merge data
all.trans <- merge(lgg_data, meta, by = "row.names")

rownames(all.trans) <- merge.data$Row.names

all.trans$Row.names <- NULL




# PreProcessing steps
library(caret)
# 1. Remove near-zero variance predictors

all.zero <- preProcess(all.trans, method = 'nzv', uniqueCut = 15)
all.trans <- predict(all.zero, all.trans)

# 2. Center the data
all.center <- preProcess(all.trans, method = 'center')
all.trans <- predict(all.center, all.trans)

# 3. Remove highly correlated features
all.corr <- preProcess(all.trans, method = 'corr', cutoff = 0.1)
all.trans <- predict(all.corr, all.trans)
dim(all.trans)



#TO WORK ON
# Splitting into training and testing sets (70:30 split)
intrain <- createDataPartition(y = all.trans$IDH.status, p = 0.7) [[1]]
length(intrain)


# Separate training and test sets
train.lgg <- all.trans[intrain, ]
test.lgg <- all.trans[-intrain, ]

# Check dimensions of training and test sets
dim(train.lgg)
dim(test.lgg)

#train
#control group
ctrl.lgg <- trainControl(method = 'cv', number = 5)


#train
knn.lgg <- train(IDH.status~.,
                  data = train.lgg,
                  method = 'knn',
                  trControl = ctrl.lgg,
                  tuneGrid = data.frame (k=1:20))

#best k
knn.lgg$bestTune


#predict
trainPred <- predict(knn.lgg, newdata = train.lgg)
testPred <- predict(knn.lgg, newdata = test.lgg)


#interpretation
#confusion matrix
#convert to factors
train.lgg$IDH.status <- factor(train.lgg$IDH.status)
test.lgg$IDH.status <- factor(test.lgg$IDH.status)
levels(testPred)  
levels(train.lgg$IDH.status)
#check length
length(testPred)  
length(test.lgg$IDH.status)
# Check for NA values in the predictions
sum(is.na(testPred))  # Count the number of NA predictions
# Convert to factor if necessary
test.lgg$IDH.status <- factor(test.lgg$IDH.status, levels = levels(train.lgg$IDH.status))

# Ensure test set has only the same columns as the training set
test.lgg <- test.lgg[, colnames(train.lgg)]

testPred <- predict(knn.lgg, newdata = test.lgg)


#confusion matrix
library(caret)  
confusionMatrix(trainPred, train.lgg$IDH.status)
confusionMatrix(testPred, test.lgg$IDH.status)



#determine variable importance
library(iml)
# First, create the explainer object  
explainer.lgg <- explain(knn.lgg,  
                         data = train.lgg,  
                         label = 'knn',  
                         y = as.numeric(train.lgg$IDH.status))  

# Now, calculate feature importance separately  
importance.lgg <- feature_importance(explainer.lgg, n_sample = 100)  

# Print the head and tail of the important variables  
head(importance.lgg$variable)  
tail(importance.lgg$variable)  

# Plot the feature importance  
plot(importance.lgg)
                    
                          
