library(dplyr)
library(magrittr)
library(xgboost)
library(tibble)
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))
# Read host genes VST values
vsd_mat <- read.csv(cv_vst_csv, row.names=1)
# Output variable to predict (needs to convert numeric into factors)
y <- as.factor(classifier_meta$LRTI_adjudication)
# Need to change the order and minus one to match 0/1
y <- factor(y, levels = c("No Evidence", "Definite"))
y_numeric <- as.numeric(y) - 1
# Table of features selected by lasso on each train-test split
lasso_features <- read.csv(paste0(lasso_results_prefix, "coefs.csv"),
                           stringsAsFactors = F)
# Table of features selected by key_lasso on each train-test split
key_lasso_features <- read.csv(paste0(key_lasso_results_prefix, "coefs.csv"),
                               stringsAsFactors = F)


# XGBoost classifier trained on key_lasso genes -----------------------------------------------------------

nested_cv_xgb_key <- function(test_fold) {
  # Extract the lasso features from the original paper
  keep <- key_lasso_features %>%
    dplyr::filter(gene != '(Intercept)') %>%
    .$gene
  
  # Subset the data to the selected genes
  X <- t(vsd_mat)[, keep]
  
  # Define the outer training and test sets
  train_indices <- classifier_meta$fold != test_fold
  test_indices <- classifier_meta$fold == test_fold
  
  # Prepare data to be loaded into model training
  # Dataframe containing training set X and label y
  train_df <- cbind(as.data.frame(X[train_indices, ]), y=y_numeric[train_indices])
  # Convert X training set to matrix
  train_matrix <- as.matrix(train_df[, -ncol(train_df)])
  # Convert X testing set to matrix
  test_matrix <-  as.matrix(X[test_indices, ])
  
  # Set the seed for reproducibility
  set.seed(571731518 + test_fold) # Different seed for each outer fold for variability
  
  # Default booster parameters for XGBoost
  params_booster <- list(booster = 'gbtree', eta = 0.3, gamma = 0, 
                         max.depth = 8, subsample = 1, colsample_bytree = 1, min_child_weight = 1, 
                         objective = "binary:logistic", eval_metric = "error")
  
  # Train the best XGBoost model with cross-validation
  bst.cv <- xgb.cv(data = train_matrix, 
                   label = train_df$y, 
                   params = params_booster,
                   nrounds = 300, 
                   nfold = 5,
                   showsd = T,
                   print_every_n = 10,
                   early_stopping_round = 20
  )
  
  xgb_cv_res <- data.frame(TRAINING_ERROR = bst.cv$evaluation_log$train_error_mean, 
                           VALIDATION_ERROR = bst.cv$evaluation_log$test_error_mean, # Don't confuse this with the test data set. 
                           ITERATION = bst.cv$evaluation_log$iter) %>%
    mutate(MIN = VALIDATION_ERROR == min(VALIDATION_ERROR))
  # Obtain the best nrounds
  best_nrounds <- xgb_cv_res %>%
    filter(MIN) %>%
    pull(ITERATION)
  
  # Train the model using the best model
  xgb_model <- xgboost(data = train_matrix, 
                       label = train_df$y, 
                       nrounds = best_nrounds, 
                       params = params_booster,
                       set.seed(571731518 + test_fold))
  
  # Predictions on the testing set
  predictions <- predict(xgb_model, test_matrix)
  # Assign sample names because caret predict doesn't preserve row name
  names(predictions) <- rownames(X)[test_indices]
  
  # Return XGBoost and predictions on the test set
  list(test_fold=test_fold, mod=xgb_model, pred=predictions)
}

# Run on all the train-test splits
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_xgb_key(i))

# Output the probabilities
cv_results %>%
  lapply(function(x) x$pred) %>%
  do.call(what=c) %>%
  {data.frame(pred=.)} %>%
  tibble::rownames_to_column("sample_name") %>%
  left_join(., classifier_meta %>% 
              dplyr::select(sample_name, LRTI_adjudication), by="sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after=sample_name) %>% 
  write.csv(paste0(key_lassoXGBoost_results_prefix, "preds.csv"), row.names=F)

# XGBoost classifier trained on lasso genes ----------------------------------------------

nested_cv_xgb_14 <- function(test_fold) {
  # Extract the lasso features from the original paper
  keep <- lasso_features %>%
    dplyr::filter(gene != '(Intercept)') %>%
    .$gene
  
  # Subset the data to the selected genes
  X <- t(vsd_mat)[, keep]
  
  # Define the outer training and test sets
  train_indices <- classifier_meta$fold != test_fold
  test_indices <- classifier_meta$fold == test_fold
  
  # Prepare data to be loaded into model training
  # Dataframe containing training set X and label y
  train_df <- cbind(as.data.frame(X[train_indices, ]), y=y_numeric[train_indices])
  # Convert X training set to matrix
  train_matrix <- as.matrix(train_df[, -ncol(train_df)])
  # Convert X testing set to matrix
  test_matrix <-  as.matrix(X[test_indices, ])
  
  # Set the seed for reproducibility
  set.seed(571731518 + test_fold) # Different seed for each outer fold for variability
  
  # Default booster parameters for XGBoost
  params_booster <- list(booster = 'gbtree', eta = 0.3, gamma = 0, 
                         max.depth = 8, subsample = 1, colsample_bytree = 1, min_child_weight = 1, 
                         objective = "binary:logistic", eval_metric = "error")
  
  # Train the best XGBoost model with cross-validation
  bst.cv <- xgb.cv(data = train_matrix, 
                   label = train_df$y, 
                   params = params_booster,
                   nrounds = 300, 
                   nfold = 5,
                   showsd = T,
                   print_every_n = 10,
                   early_stopping_round = 20
  )
  
  xgb_cv_res <- data.frame(TRAINING_ERROR = bst.cv$evaluation_log$train_error_mean, 
                           VALIDATION_ERROR = bst.cv$evaluation_log$test_error_mean, # Don't confuse this with the test data set. 
                           ITERATION = bst.cv$evaluation_log$iter) %>%
    mutate(MIN = VALIDATION_ERROR == min(VALIDATION_ERROR))
  # Obtain the best nrounds
  best_nrounds <- xgb_cv_res %>%
    filter(MIN) %>%
    pull(ITERATION)
  
  # Train the model using the best model
  xgb_model <- xgboost(data = train_matrix, 
                       label = train_df$y, 
                       nrounds = best_nrounds, 
                       params = params_booster,
                       set.seed(571731518 + test_fold))
  
  # Predictions on the testing set
  predictions <- predict(xgb_model, test_matrix)
  # Assign sample names because caret predict doesn't preserve row name
  names(predictions) <- rownames(X)[test_indices]

  # Return XGBoost and predictions on the test set
  list(test_fold=test_fold, mod=xgb_model, pred=predictions)
}

# Run on all the train-test splits
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_xgb_14(i))

# Output the probabilities
cv_results %>%
  lapply(function(x) x$pred) %>%
  do.call(what=c) %>%
  {data.frame(pred=.)} %>%
  tibble::rownames_to_column("sample_name") %>%
  left_join(., classifier_meta %>% 
              dplyr::select(sample_name, LRTI_adjudication), by="sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after=sample_name) %>% 
  write.csv(paste0(lassoXGBoost_results_prefix, "preds.csv"), row.names=F)

# XGBoost classifier trained on all genes -----------------------------------------------------------------

nested_cv_xgb_all <- function(test_fold) {
  # all
  X <- t(vsd_mat)
  
  # Define the outer training and test sets
  train_indices <- classifier_meta$fold != test_fold
  test_indices <- classifier_meta$fold == test_fold
  
  # Prepare data to be loaded into model training
  # Dataframe containing training set X and label y
  train_df <- cbind(as.data.frame(X[train_indices, ]), y=y_numeric[train_indices])
  # Convert X training set to matrix
  train_matrix <- as.matrix(train_df[, -ncol(train_df)])
  # Convert X testing set to matrix
  test_matrix <-  as.matrix(X[test_indices, ])
  
  # Set the seed for reproducibility
  set.seed(571731518 + test_fold) # Different seed for each outer fold for variability
  
  # Default booster parameters for XGBoost
  params_booster <- list(booster = 'gbtree', eta = 0.3, gamma = 0, 
                         max.depth = 8, subsample = 1, colsample_bytree = 1, min_child_weight = 1, 
                         objective = "binary:logistic", eval_metric = "error")
  
  # Train the best XGBoost model with cross-validation
  bst.cv <- xgb.cv(data = train_matrix, 
                   label = train_df$y, 
                   params = params_booster,
                   nrounds = 300, 
                   nfold = 5,
                   showsd = T,
                   print_every_n = 10,
                   early_stopping_round = 20
  )
  
  xgb_cv_res <- data.frame(TRAINING_ERROR = bst.cv$evaluation_log$train_error_mean, 
                           VALIDATION_ERROR = bst.cv$evaluation_log$test_error_mean, # Don't confuse this with the test data set. 
                           ITERATION = bst.cv$evaluation_log$iter) %>%
    mutate(MIN = VALIDATION_ERROR == min(VALIDATION_ERROR))
  # Obtain the best nrounds
  best_nrounds <- xgb_cv_res %>%
    filter(MIN) %>%
    pull(ITERATION)
  
  # Train the model using the best model
  xgb_model <- xgboost(data = train_matrix, 
                       label = train_df$y, 
                       nrounds = best_nrounds, 
                       params = params_booster,
                       set.seed(571731518 + test_fold))
  
  # Predictions on the testing set
  predictions <- predict(xgb_model, test_matrix)
  # Assign sample names because caret predict doesn't preserve row name
  names(predictions) <- rownames(X)[test_indices]
  
  # Return XGBoost and predictions on the test set
  list(test_fold=test_fold, mod=xgb_model, pred=predictions)
}

# Run on all the train-test splits
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) nested_cv_xgb_all(i))

# Output the probabilities
cv_results %>%
  lapply(function(x) x$pred) %>%
  do.call(what=c) %>%
  {data.frame(pred=.)} %>%
  tibble::rownames_to_column("sample_name") %>%
  left_join(., classifier_meta %>% 
              dplyr::select(sample_name, LRTI_adjudication), by="sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after=sample_name) %>% 
  write.csv(paste0(XGBoost_results_prefix, "preds.csv"), row.names=F)