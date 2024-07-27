library(dplyr)
library(magrittr)
library(stringr)
source("constants.R")

# Samples with NA values removed
columnData_virus <- classifier_meta %>%
  filter (sample_name %in% colnames(virus_counts)) 
columnData_species <- classifier_meta %>%
  filter (sample_name %in% colnames(bf_counts_species))

# CHANGE
norm_microbe <- counts(dds_bf_species, normalized = TRUE) 
#norm_microbe <- counts(dds_bf_genus, normalized = TRUE) 
#norm_microbe <- counts(dds_virus, normalized = TRUE) 
col = columnData_species
#col = columnData_viru

norm_df <- as.data.frame(norm_microbe, stringsAsFactors = FALSE)

# Gene counts
selected_df <- norm_df %>%
  filter(rownames(.) == "Haemophilus influenzae") # CHANGE
selected_df = t(selected_df)
selected_df <- pivot_longer(
  data = as.data.frame(t(selected_df)),
  cols = everything(),
  names_to = "Sample",
  values_to = "Gene"
)

# Group
group <- as.factor(col$LRTI_adjudication)

# Combine
plot_data <- data.frame(
  Gene = selected_df$Gene,
  Group = group
)

# t test
t.test(Gene ~ group, data = plot_data, alternative = "two.sided", var.equal = FALSE)

# Plot
ggplot(plot_data, aes(x = Group, y = Gene)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "lightblue") +  # Remove default outliers to avoid overplotting
  geom_jitter(color = "black", size = 1, width = 0.2, height = 0) +  # Add jitter points
  theme_minimal() +
  labs(
    title = "Expression of Haemophilus influenzae",  # Corrected typo in the title
    x = "Group",
    y = "Normalized nt counts"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title, bold, and adjust size
    axis.title = element_text(face = "bold"),  # Optionally bold axis titles for better visibility
    plot.title.position = "plot"  # Center title in the entire plot area
  )
