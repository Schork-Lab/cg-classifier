### This training script trains 5 separate
### ML algorthims, runs a grid search and 
### then tests performance on a new set

library("caret")
library("randomForest")
library("kernlab")
library("C50")
library("pROC")
library("doMC")
library("ROCR")
library("plyr")
library("dplyr")

### LOAD UP SOME DATA

setwd("~/Projects/variant_filter/explore")
load("sample_data.Rdata")

### TRAINING MODELS

registerDoMC(10)

## control object for training all models
ctrl = trainControl(method="repeatedcv",
                    repeats=1,
                    number=3,
                    summaryFunction=twoClassSummary,
                    classProbs=T)

## Random Forest
rf = train(Y ~ .,
           data=d,
           trControl=ctrl,
           method="rf",
           metric="ROC")
print(rf)

## Stocastic gradient boosting
gb = train(Y~.,
           data=d,
           trControl=ctrl,
           metric="ROC",
           method="gbm",
           tuneGrid=expand.grid(shrinkage=c(1e-3, 1e-2, 1e-1),
                                n.trees=c(100,1000),
                                interaction.depth=c(1,3)))
print(gb)

## Support vector machine with RBF kernel
svm.rbf = train(Y ~ .,
            data=d,
            trControl=ctrl,
            method="svmRadialCost",
            metric="ROC",
            tuneGrid=expand.grid(C=c( 0.25, 0.5, 1)))
print(svm.rbf)

## Linear support vector machine
svm.linear = train(Y ~ .,
                data=d,
                trControl=ctrl,
                method="svmLinear",
                metric="ROC",
                tuneGrid=expand.grid(C=c(0.25, 0.5, 1)))
print(svm.linear)

## C5.0 boosting
c5 = train(Y ~ .,
           data=d,
           trControl=ctrl,
           method="C5.0",
           metric="ROC")
print(c5)

### EVALUTATE ON TESTING SET

## create list of models
models = list(rf, gb, svm.rbf, svm.linear, c5)
names(models) = c("Random Forest",
                  "Stocastic Gradient Boosting",
                  "SVM with RBF kernel",
                  "SVM with Linear Kernel",
                  "C5.0 boosting")

## make predictions for new dataset
predictions = lapply(models, 
                     function(x) predict(x, newdata=d.testing, type="prob"))

## calculate AUCs
aucs = lapply(predictions, 
              function(x) auc(roc(d.testing$Y, x[,1])))

## ...other stuff here