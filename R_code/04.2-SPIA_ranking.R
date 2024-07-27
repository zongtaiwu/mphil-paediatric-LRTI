library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))

# Consensus results from DEG and key SPIA results
# IMPORT THE TABLE AFTER VENN DIAGRAM AND EXAMINATION
spia_results <- read.csv(consensus_spia_csv) # 41 pathways

# Read host genes VST values and filter the ones in genes of interest
vsd_mat_spia <- read.csv(cv_vst_csv, row.names=1)

# Prepare the 41 pathways with genes & Entrez IDs in each pathway
spia_rank <- deseq_key_spia_selected %>%
  filter(KEGG_pathway %in% spia_results$KEGG_pathway) %>%
  select(KEGG_pathway, Entrez_id) %>% # Remove the pathways with only one gene
  filter(grepl(" ", Entrez_id)) # 36 pathways


# LR and ranking by ROC-AUC----------------------------------------------------------------------
# Empty df
spia_ranking_df <- data.frame(Name = character(), Fold1 = numeric(), Fold2 = numeric(), 
                              Fold3 = numeric(), Fold4 = numeric(), Fold5 = numeric())

# Prepare y
y <- as.factor(classifier_meta$LRTI_adjudication)
# Note it's No. Evidence -- R variable
levels(y) <- make.names(levels(y))

# Function that calculates F1 score in each pathway
metrics_spia <- function(preds_df) {
  
  # Calculate metrics for each fold
  metrics_df <- preds_df %>%
    dplyr::inner_join(classifier_meta %>% dplyr::select(sample_name, fold), by = "sample_name") %>%
    dplyr::group_by(fold) %>%
    dplyr::reframe(
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
      }
    ) %>%
    dplyr::ungroup()
  
  return(metrics_df)
}

# Iterate over 36 pathways
for (i in 1:36){
  # Group of genes in each pathway
  numbers <- unlist(strsplit(spia_rank$Entrez_id[i], " "))
  list <- as.character(numbers)
  # Extract genes within each spia_rank pathway
  # res_key has gene id matched with entrez id
  keep <- row.names(res_key[res_key$entrez %in% list, ])
  
  # Conduct a 5 fold cross validation to calcuate 
  LR_spia <- function(test_fold) {
    
    X <- t(vsd_mat_spia)[ ,keep]
    print(paste("Dimensions of X:", dim(X)))
    
    # Define the outer training and test sets
    train_indices <- classifier_meta$fold != test_fold
    test_indices <- classifier_meta$fold == test_fold
    
    X_train <- X[train_indices, ]
    y_train <- y[train_indices]
    X_test <- X[test_indices, ]
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
  
  # Output the prob into cv_results
  cv_results <- lapply(1:max(classifier_meta$fold),
                       function(i) LR_spia (i))
  cv_results <- do.call(rbind, cv_results)
  # Convert back to No Evidence
  cv_results <- cv_results %>%
    mutate(LRTI_adjudication = recode(LRTI_adjudication, 
                                      "No.Evidence" = "No Evidence"))
  
  # Calculate F1 scores
  f1_spia <- metrics_spia(cv_results)
  # List of F1 scores in folds
  f1_list <- list(
    Name = spia_rank$KEGG_pathway[i],
    `Fold1` = f1_spia$f1[1],
    `Fold2` = f1_spia$f1[2],
    `Fold3` = f1_spia$f1[3],
    `Fold4` = f1_spia$f1[4],
    `Fold5` = f1_spia$f1[5]
  )
 
  # Add into the summary dataframe
  new_row <-  as.data.frame(f1_list, stringsAsFactors = FALSE)
  spia_ranking_df <- rbind(spia_ranking_df, new_row)
}


# Rank and plot -----------------------------------------------------------
# Convert to the long format
long_spia_ranking_df <- spia_ranking_df %>%
  pivot_longer(
    cols = starts_with("Fold"),
    names_to = "Fold",
    values_to = "Value"
  )

# Rank and Plot the ranking by mean value
mean_values <- long_spia_ranking_df %>%
  group_by(Name) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE)) %>%
  arrange(MeanValue)

long_spia_ranking_df$Name <- factor(long_spia_ranking_df$Name, levels = mean_values$Name)

# Pathway and their activated/inhibited status
status <- deg_spia %>%
  select(Name, Status)

long_spia_ranking_df <- merge(long_spia_ranking_df, status, by = "Name")

# plot
ggplot(long_spia_ranking_df, aes(x = Name, y = Value, fill = Status)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Activated" = "coral1", "Inhibited" = "cornflowerblue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Pathway Name", y = "F1 score", title = "LR classifier using key genes in each pathway") +
  coord_flip() # Flip coordinates to make the plot horizontal



### TOP 10 pathways
# Prepare plot data
long_spia_ranking_df <- spia_ranking_df %>%
  pivot_longer(
    cols = starts_with("Fold"),
    names_to = "Fold",
    values_to = "Value"
  )

top_10_pathways <- long_spia_ranking_df %>%
  group_by(Name) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE)) %>%
  arrange(desc(MeanValue)) %>%
  slice_head(n = 10)

# Filter the long_spia_ranking_df to include only the top 10 pathways
filtered_spia_ranking_df <- long_spia_ranking_df %>%
  filter(Name %in% top_10_pathways$Name)

# Rank and Plot the ranking by mean value
mean_values <- filtered_spia_ranking_df %>%
  group_by(Name) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE)) %>%
  arrange(MeanValue)

filtered_spia_ranking_df$Name <- factor(filtered_spia_ranking_df$Name, levels = mean_values$Name)

filtered_spia_ranking_df <- merge(filtered_spia_ranking_df, status, by = "Name")

# plot
ggplot(filtered_spia_ranking_df, aes(x = Name, y = Value, fill = Status)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Activated" = "coral1", "Inhibited" = "cornflowerblue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Pathway Name", y = "F1 score", title = "Top 10 LR classifier using key genes in each pathway") +
  coord_flip() # Flip coordinates to make the plot horizontal

