library(dplyr)
library(magrittr)
library(tibble)
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
# Table of lasso-selected features on key protein genes
key_lasso_features <- read.csv(paste0(key_lasso_results_prefix, "coefs.csv"),
                               stringsAsFactors = F)
# Table of lasso-selected features on all genes
lasso_features <- read.csv(paste0(lasso_results_prefix, "coefs.csv"),
                           stringsAsFactors = F)


# LR classifier trained on key_lasso genes ------------------------------------------------------

lasso_cv_LR_key <- function(test_fold) {
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
  
  # Ensure the levels are correctly set so "Definite" is 1 and "No Evidence" is 0
  y_train <- ifelse(y_train == "Definite", 1, 0)
  
  # Train the logistic regression model
  logistic_model <- glm(y_train ~ ., data = data.frame(X_train), family = binomial(link = "logit"))
  
  # Predict on the outer test set
  pred_prob <- predict(logistic_model, newdata = data.frame(X_test), type = "response")
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}


# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) lasso_cv_LR_key (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(key_LR_results_prefix, "preds.csv"), row.names=F)


# LR classifier trained on lasso genes ----------------------------------------------------------

lasso_cv_LR_14 <- function(test_fold) {
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
  
  # Ensure the levels are correctly set so "Definite" is 1 and "No Evidence" is 0
  y_train <- ifelse(y_train == "Definite", 1, 0)
  
  # Train the logistic regression model
  logistic_model <- glm(y_train ~ ., data = data.frame(X_train), family = binomial(link = "logit"))
  
  # Predict on the outer test set
  pred_prob <- predict(logistic_model, newdata = data.frame(X_test), type = "response")
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}


# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) lasso_cv_LR_14 (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(lasso_LR_results_prefix, "preds.csv"), row.names=F)


# LR classifier trained on all genes ---------------------------------------------------------

lasso_cv_LR_all <- function(test_fold) {
  # all
  X <- t(vsd_mat)
  
  # Define the outer training and test sets
  train_indices <- classifier_meta$fold != test_fold
  test_indices <- classifier_meta$fold == test_fold
  
  X_train <- X[train_indices,]
  y_train <- y[train_indices]
  X_test <- X[test_indices,]
  y_test <- y[test_indices]
  
  # Ensure the levels are correctly set so "Definite" is 1 and "No Evidence" is 0
  y_train <- ifelse(y_train == "Definite", 1, 0)
  
  # Train the logistic regression model
  logistic_model <- glm(y_train ~ ., data = data.frame(X_train), family = binomial(link = "logit"))
  
  # Predict on the outer test set
  pred_prob <- predict(logistic_model, newdata = data.frame(X_test), type = "response")
  
  # Return the predictions and actual labels for evaluation
  data.frame(sample_name = rownames(X_test), 
             LRTI_adjudication = y_test, 
             pred = pred_prob)
}


# Run the nested cross-validation for all outer folds
cv_results <- lapply(1:max(classifier_meta$fold),
                     function(i) lasso_cv_LR_all (i))
cv_results <- do.call(rbind, cv_results)
# Convert back to No Evidence
cv_results <- cv_results %>%
  mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                    "No.Evidence" = "No Evidence"))
# Output the probabilities
cv_results %>%
  write.csv(paste0(LR_results_prefix, "preds.csv"), row.names=F)