library(dplyr)
library(magrittr)
library(stringr)
source("constants.R")

# Samples with NA values removed
columnData_virus <- classifier_meta %>%
  filter (sample_name %in% colnames(virus_counts)) 

# Normalised viral counts
norm_virus <- counts(dds_virus, normalized = TRUE)
selected_virus <- norm_virus["Human orthopneumovirus", ]
# Normalised bacterial/fungal counts
norm_bf_species <- counts(dds_bf_species,  normalized = TRUE)

# Samples containing both viral and bacteria/fungal species
aligned_samples <- intersect(names(selected_virus), colnames(norm_bf_species))

virus_data <- selected_virus[aligned_samples, drop = FALSE] 
bacteria_data <- norm_bf_species[, aligned_samples]

# Function to compute kendall's tau correlational test
compute_correlations <- function(virus, bacteria, method) {
  result <- apply(bacteria, 1, function(bacteria_row) {
    # Check for zero standard deviation in bacteria row
    if (sd(bacteria_row) == 0) {
      return(c(Correlation = NA, P.value = NA))
    } else {
      test <- cor.test(as.numeric(virus), as.numeric(bacteria_row), method = method)
      return(c(Correlation = test$estimate, P.value = test$p.value))
    }
  })
  
  # Convert the result to a data frame and properly name the columns
  correlation_df <- as.data.frame(t(result))
  correlation_df$Bacteria <- rownames(bacteria)
  correlation_df <- correlation_df %>%
    rename(Correlation = V1, P.value = V2) %>%  # Rename the columns correctly
    select(Bacteria, Correlation, P.value)
  
  return(correlation_df)
}

kendall_results <- compute_correlations(virus_data, bacteria_data, "kendall") %>%
  arrange(P.value) %>%
  filter(Correlation > 0)

# Output the table
write.csv(kendall_results, RSV_corr_csv,
          row.names = FALSE)
