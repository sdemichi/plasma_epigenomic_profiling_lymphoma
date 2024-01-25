# train/test split performed previously

# libraries
library(tidyverse)
library(sva)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(MetBrewer)
library(caret)
library(glmnet)
library(e1071)
library(edgeR)
library(umap)
library(DMwR)
library(pROC)
library(gbm)
library(klaR)

# import train data. Feature selection (with differential analysis) performed prior to this step. Imput features passed significance threshold from differential analysis.
setwd("/path/to/training/data/")
load("train_data_differential_analysis_featureset.RData")

# specify tune grid for hyperparameter tuning; tuning hyperparameters depend on model used.
tune_grid_rf <- expand.grid(
  mtry = seq(2, 10, by = 2) # for rf
)

# for glmnet
#tune_grid_glmnet <- expand.grid(
#  alpha = seq(0, 1, by = 0.1),    # Elastic Net mixing parameter (0: Ridge, 1: Lasso)
#  lambda = seq(0.001, 1, by = 0.1) # Regularization strength
#)

# Perform Recursive Feature Elimination (RFE) with cross-validation.
ctrl <- trainControl(
  method = "repeatedcv",  # use cross-validation for evaluation
  repeats = 5, # number of model iterations
  number = 10,  # number of folds for cross-validation
  # search = "grid", # grid search if tune grid is not specified
  sampling = "smote", # synthetic minority oversampling; sample features from similar samples in the minority set. Use this to construct new samples, to address class imbalance.
  classProbs = TRUE,
  verbose = FALSE
)

# transpose feature matrix
test <- t(de_features_matrix)

# sink("04_08_2023_NB_noPCA_10kb_without_hyperparameter_tuning_accuracy_error_log.txt") # if the model fails to generate, use this to troubleshoot.
model <- train(test, 
               as.factor(y_train), 
               method = "rf", # rf being recognized
               #method = "glmnet", # elastic net regression
               #method = "glm", # logistic regression; turn off tuneGrid.
               #method = "nb", # naive bayes; turn off tuneGrid.
               trControl = ctrl, 
               tuneGrid = tune_grid_rf, # for rf
               #tuneGrid = tune_grid_glmnet, #for glmnet
               metric = "Accuracy",
               #metric = "ROC",
               #metric = "Kappa", # for unbalanced sets; since we're using SMOTE, accuracy should be okay
               ntree = 1000 # for rf
)
# sink()
# print(model) # see results

# predictors(model) # extract predictors for further analysis (aka feature validation)

# save workspace image for further analysis
setwd("/path/to/training/data/")
save.image("caret_rf_model_after_differential_analysis_lymphoma_non_lymphoma_controls.RData")
