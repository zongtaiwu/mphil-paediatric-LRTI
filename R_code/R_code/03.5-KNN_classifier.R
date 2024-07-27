library(dplyr)
library(magrittr)
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


# KNN classifier trained on key_lasso genes ---------------------------------------------------------------
nested_cv_knn_key <- function(test_fold) {
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
  
  # Define the grid of hyperparameters to search
  tune_grid <- expand.grid(k = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 30))
  
  # Train the model with inner cross-validation for hyperparameter tuning
  set.seed(352607443 + test_fold) # Different seed for each outer fold for variability
  knn_tuned <- train(x = X_train, y = y_train,
                     method = "knn",
                     tuneGrid = tune_grid,
                     trControl = inner_control,
                     metric = "ROC")
  
  # Print the best parameters found
  best_k <- knn_tuned$bestTune$k
  print(paste("Best k:", best_k))
  
  # Train the final model on the entire outer training set using the best hyperparameters
  final_model <- knn3(x = X_train, y = y_train, k = best_k)
  
  # Predict on the outer test set
  pred_prob <- predict(final_model, newdata = X_test, type = "prob")[, "Definite"]
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}

# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_knn_key (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(key_lassoKNN_results_prefix, "preds.csv"), row.names=F)



# KNN classifier trained on lasso genes --------------------------------------------------

nested_cv_knn_14 <- function(test_fold) {
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
  
  # Define the grid of hyperparameters to search
  tune_grid <- expand.grid(k = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 30))
  
  # Train the model with inner cross-validation for hyperparameter tuning
  set.seed(352607443 + test_fold) # Different seed for each outer fold for variability
  knn_tuned <- train(x = X_train, y = y_train,
                     method = "knn",
                     tuneGrid = tune_grid,
                     trControl = inner_control,
                     metric = "ROC")
  
  # Print the best parameters found
  best_k <- knn_tuned$bestTune$k
  print(paste("Best k:", best_k))
  
  # Train the final model on the entire outer training set using the best hyperparameters
  final_model <- knn3(x = X_train, y = y_train, k = best_k)
  
  # Predict on the outer test set
  pred_prob <- predict(final_model, newdata = X_test, type = "prob")[, "Definite"]
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}

# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_knn_14 (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(lassoKNN_results_prefix, "preds.csv"), row.names=F)

# KNN classifier trained on all genes---------------------------------------------------------------------

nested_cv_knn_all <- function(test_fold) {
  
  # all
  X <- t(vsd_mat)
  
  # Define the outer training and test sets
  train_indices <- classifier_meta$fold != test_fold
  test_indices <- classifier_meta$fold == test_fold
  
  X_train <- X[train_indices,]
  y_train <- y[train_indices]
  X_test <- X[test_indices,]
  y_test <- y[test_indices]
  
  # Define the control for inner cross-validation
  inner_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  # Define the grid of hyperparameters to search
  tune_grid <- expand.grid(k = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 30))
  
  # Train the model with inner cross-validation for hyperparameter tuning
  set.seed(352607443 + test_fold) # Different seed for each outer fold for variability
  knn_tuned <- train(x = X_train, y = y_train,
                     method = "knn",
                     tuneGrid = tune_grid,
                     trControl = inner_control,
                     metric = "ROC")
  
  # Print the best parameters found
  best_k <- knn_tuned$bestTune$k
  print(paste("Best k:", best_k))
  
  # Train the final model on the entire outer training set using the best hyperparameters
  final_model <- knn3(x = X_train, y = y_train, k = best_k)
  
  # Predict on the outer test set
  pred_prob <- predict(final_model, newdata = X_test, type = "prob")[, "Definite"]
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}

# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_knn_all (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(KNN_results_prefix, "preds.csv"), row.names=F)