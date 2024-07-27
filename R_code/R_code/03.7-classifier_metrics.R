library(dplyr)
library(magrittr)
library(pROC)
library(caret)
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))

# Function to output the per-fold metrics from the pred probability,
# and to summarise across all folds
# Set threshold = 0.5
output_metrics_prob <- function(preds_df, output_suffix) {
  
  # Calculate metrics for each fold
  metrics_df <- preds_df %>%
    dplyr::inner_join(classifier_meta %>% dplyr::select(sample_name, fold), by = "sample_name") %>%
    dplyr::group_by(fold) %>%
    dplyr::reframe(
      accuracy = {
        threshold_pred <- ifelse(pred >= 0.5, 'Definite', 'No Evidence')
        confusion <- caret::confusionMatrix(as.factor(threshold_pred), as.factor(LRTI_adjudication), positive = 'Definite')
        confusion$overall["Accuracy"]
      },
      precision = {
        threshold_pred <- ifelse(pred >= 0.5, 'Definite', 'No Evidence')
        confusion <- caret::confusionMatrix(as.factor(threshold_pred), as.factor(LRTI_adjudication), positive = 'Definite')
        confusion$byClass["Pos Pred Value"]
      },
      recall = {
        threshold_pred <- ifelse(pred >= 0.5, 'Definite', 'No Evidence')
        confusion <- caret::confusionMatrix(as.factor(threshold_pred), as.factor(LRTI_adjudication), positive = 'Definite')
        confusion$byClass["Sensitivity"]
      },
      f1 = {
        threshold_pred <- ifelse(pred >= 0.5, 'Definite', 'No Evidence')
        confusion <- caret::confusionMatrix(as.factor(threshold_pred), as.factor(LRTI_adjudication), positive = 'Definite')
        precision <- confusion$byClass["Pos Pred Value"]
        recall <- confusion$byClass["Sensitivity"]
        if (is.na(precision) | is.na(recall) | precision + recall == 0) {
          0  # Handle cases where precision or recall is NA or where precision + recall is zero
        } else {
          2 * (precision * recall) / (precision + recall)
        }
      },
      roc = {
        roc_curve <- pROC::roc(as.factor(LRTI_adjudication), pred, levels = c("No Evidence", "Definite"), direction = "<")
        as.numeric(pROC::auc(roc_curve))
      }
    ) %>%
    dplyr::ungroup()
  
  # Write individual fold metrics to CSV
  write.csv(metrics_df, paste0(cv_metrics_prefix, output_suffix, "_fold.csv"), row.names = FALSE)
  
  # Summarize metrics across folds
  summary_df <- metrics_df %>%
    dplyr::reframe(
      Metric = c("Accuracy", "Precision", "Recall", "F1 Score", "ROC AUC"),
      Mean = c(mean(accuracy, na.rm = TRUE), mean(precision, na.rm = TRUE), mean(recall, na.rm = TRUE), mean(f1, na.rm = TRUE), mean(roc, na.rm = TRUE)),
      SD = c(sd(accuracy, na.rm = TRUE), sd(precision, na.rm = TRUE), sd(recall, na.rm = TRUE), sd(f1, na.rm = TRUE), sd(roc, na.rm = TRUE)),
      Median = c(median(accuracy, na.rm = TRUE), median(precision, na.rm = TRUE), median(recall, na.rm = TRUE), median(f1, na.rm = TRUE), median(roc, na.rm = TRUE)),
      IQR = c(IQR(accuracy, na.rm = TRUE), IQR(precision, na.rm = TRUE), IQR(recall, na.rm = TRUE), IQR(f1, na.rm = TRUE), IQR(roc, na.rm = TRUE))
    )
  
  # Write summary metrics to CSV
  write.csv(summary_df, paste0(cv_metrics_prefix, output_suffix, "_summary.csv"), row.names = FALSE)
}


# F1 score -------------------------------------------------------------------------
# Output the pred probability to 
# LR
read.csv(paste0(key_LR_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "key_lasso_LR")

read.csv(paste0(lasso_LR_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "lasso_LR")

read.csv(paste0(LR_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "LR")


# RF
read.csv(paste0(key_lassoRF_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "key_lasso_RF")

read.csv(paste0(lassoRF_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "lasso_RF")

read.csv(paste0(RF_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "RF")


# SVM
read.csv(paste0(key_lassoSVM_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "key_lasso_SVM")

read.csv(paste0(lassoSVM_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "lasso_SVM")

read.csv(paste0(SVM_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "SVM")

# KNN
read.csv(paste0(key_lassoKNN_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "key_lassoKNN")

read.csv(paste0(lassoKNN_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "lasso_KNN")

read.csv(paste0(KNN_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "KNN")

# XGBoost
read.csv(paste0(key_lassoXGBoost_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "key_lassoXGBoost")

read.csv(paste0(lassoXGBoost_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "lassoXGBoost")

read.csv(paste0(XGBoost_results_prefix, "preds.csv"), stringsAsFactors = F) ->
  host_preds_df
output_metrics_prob(host_preds_df, "XGBoost")