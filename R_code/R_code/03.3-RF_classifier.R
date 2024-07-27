library(dplyr)
library(magrittr)
library(randomForest)
library(caret)
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))
# Read host genes VST values
vsd_mat <- read.csv(cv_vst_csv, row.names=1)
# Output variable to predict (needs to convert numeric into factors)
y <- as.factor(classifier_meta$LRTI_adjudication)
# Note it's No. Evidence -- R variable
levels(y) <- make.names(levels(y))
# Table of features selected by lasso on each train-test split
lasso_features <- read.csv(paste0(lasso_results_prefix, "coefs.csv"),
                           stringsAsFactors = F)
# Table of features selected by key_lasso on each train-test split
key_lasso_features <- read.csv(paste0(key_lasso_results_prefix, "coefs.csv"),
                               stringsAsFactors = F)

# Number of trees to use for random forest
n_rf_trees <- 10000

# RF classifier trained on key_lasso genes ----------------------------------------------------------------

nested_cv_rf_key <- function(test_fold) {
  # Extract the 7 features from the lasso on key proteins
  keep <- key_lasso_features %>%
    dplyr::filter(gene != '(Intercept)') %>%
    .$gene
  
  # Subset the data to the selected genes
  X <- t(vsd_mat)[, keep]
  
  # Define the outer training and test sets
  train_indices <- classifier_meta$fold != test_fold
  test_indices <- classifier_meta$fold == test_fold
  
  X_train <- X[train_indices,]
  y_train <- y[train_indices]
  X_test <- X[test_indices,]
  y_test <- y[test_indices]
  
  # Define the control for inner cross-validation
  inner_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  # The number of variables to randomly sample as candidates at each split
  tune_grid <- expand.grid(mtry = c(1, 2, 3, 4)) 
  
  # Train the model with inner cross-validation for hyperparameter tuning
  set.seed(702074539 + test_fold) # Different seed for each outer fold for variability
  rf_tuned <- train(x = X_train, y = y_train,
                    method = "rf",
                    ntree = n_rf_trees,
                    tuneGrid = tune_grid,
                    trControl = inner_control,
                    metric = "ROC")
  
  # Select the best model from inner CV
  best_mtry <- rf_tuned$bestTune$mtry
  
  # Train the final model on the entire outer training set using the best hyperparameters
  rf_final <- randomForest(X_train, y_train, mtry = best_mtry, ntree = n_rf_trees)
  
  # Predict on the outer test set
  pred_prob <- predict(rf_final, newdata = X_test, type = 'prob')[,"Definite"]
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}

# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_rf_key (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(key_lassoRF_results_prefix, "preds.csv"), row.names=F)

# RF classifier trained on lasso genes ---------------------------------------------------------------

nested_cv_rf_14 <- function(test_fold) {
  # Extract the lasso features from the original paper
  keep <- lasso_features %>%
    dplyr::filter(gene != '(Intercept)') %>%
    .$gene
  
  # Subset the data to the selected genes
  X <- t(vsd_mat)[, keep]
  
  # Define the outer training and test sets
  train_indices <- classifier_meta$fold != test_fold
  test_indices <- classifier_meta$fold == test_fold
  
  X_train <- X[train_indices,]
  y_train <- y[train_indices]
  X_test <- X[test_indices,]
  y_test <- y[test_indices]
  
  # Define the control for inner cross-validation
  inner_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  # The number of variables to randomly sample as candidates at each split
  tune_grid <- expand.grid(mtry = c(1, 2, 3, 4)) 
  
  # Train the model with inner cross-validation for hyperparameter tuning
  set.seed(702074539 + test_fold) # Different seed for each outer fold for variability
  rf_tuned <- train(x = X_train, y = y_train,
                    method = "rf",
                    ntree = n_rf_trees,
                    tuneGrid = tune_grid,
                    trControl = inner_control,
                    metric = "ROC")
  
  # Select the best model from inner CV
  best_mtry <- rf_tuned$bestTune$mtry
  
  # Train the final model on the entire outer training set using the best hyperparameters
  rf_final <- randomForest(X_train, y_train, mtry = best_mtry, ntree = n_rf_trees)
  
  # Predict on the outer test set
  pred_prob <- predict(rf_final, newdata = X_test, type = 'prob')[,"Definite"]
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}

# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_rf_14 (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(lassoRF_results_prefix, "preds.csv"), row.names=F)


# RF classifier trained on all genes ----------------------------------------------------------------------

nested_cv_rf_all <- function(test_fold) {
  
  # Subset the data to the selected genes
  X <- t(vsd_mat)[, keep]
  
  # Define the outer training and test sets
  train_indices <- classifier_meta$fold != test_fold
  test_indices <- classifier_meta$fold == test_fold
  
  X_train <- X[train_indices,]
  y_train <- y[train_indices]
  X_test <- X[test_indices,]
  y_test <- y[test_indices]
  
  # Define the control for inner cross-validation
  inner_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  # The number of variables to randomly sample as candidates at each split
  tune_grid <- expand.grid(mtry = c(1, 2, 3, 4)) 
  
  # Train the model with inner cross-validation for hyperparameter tuning
  set.seed(702074539 + test_fold) # Different seed for each outer fold for variability
  rf_tuned <- train(x = X_train, y = y_train,
                    method = "rf",
                    ntree = n_rf_trees,
                    tuneGrid = tune_grid,
                    trControl = inner_control,
                    metric = "ROC")
  
  # Select the best model from inner CV
  best_mtry <- rf_tuned$bestTune$mtry
  
  # Train the final model on the entire outer training set using the best hyperparameters
  rf_final <- randomForest(X_train, y_train, mtry = best_mtry, ntree = n_rf_trees)
  
  # Predict on the outer test set
  pred_prob <- predict(rf_final, newdata = X_test, type = 'prob')[,"Definite"]
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}

# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_rf_all (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(RF_results_prefix, "preds.csv"), row.names=F)