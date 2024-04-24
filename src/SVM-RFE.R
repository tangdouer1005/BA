a<-read.csv("k.csv",row.names = 1)
library(tidyverse)
library(glmnet)
source('msvmRFE.R')
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
train<-a
input <- train
#k-fold crossValidation
svmRFE(input, k = 5, halve.above = 100)
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100)
top.features = WriteFeatures(results, input, save=F) 
head(top.features)
write.csv(top.features,"feature_svm.csv")

featsweep = lapply(1:7, FeatSweep.wrap, results, input) #7 variables
featsweep

#Draw
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

dev.new(width=4, height=4, bg='white')
pdf("B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) #View error rate
dev.off()

#The location of the red circle in the figure represents the lowest error rate point
which.min(errors) 
top<-top.features[1:which.min(errors), "FeatureName"]
write.csv(top,"top.csv")