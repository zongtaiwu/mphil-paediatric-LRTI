library(data.table)
library(readxl)
library(dplyr)
library(ggplot2)
library(rstatix)
microbe_raw <- read_excel("/Users/zongtaiwu/Documents/RASCALS analysis/RASCALSCriticallyIll-DiagnosticTaqManResu_DATA_LABELS_2022-02-14_1619.xlsx")

microbe <- microbe_raw %>%
  filter(number_hits > 0)
demo <- read_excel("/Users/zongtaiwu/Documents/RASCALS analysis/Demographics for cytokines.xlsx")
# Change to unclassified if Nil
demo$`Consultant consensus`[demo$`Consultant consensus` %in% c("Nil", "Unclassified")] <- "Non-bacterial"
demo$`Consultant consensus`[demo$`Consultant consensus` == "bacterial"] <- "Bacterial"
demo$`Consultant consensus`[demo$`Consultant consensus` == "non-bacterial"] <- "Non-bacterial"

demo$consultant <- demo$`Consultant consensus`

# First analysis: RSV and severity----------------------------------------------------------
# Figure 1a. Between clinical diagnosis
ggplot(data = demo, aes(x = as.factor(consultant), y = ventfree28)) +
  geom_boxplot(fill = "lightblue1", outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.7, size = 1.5) +
  ylim(0, 28) +
  labs(y = "Vent-free days",
       x = "Group",
       title = "") +
  theme_minimal()

demo$consultant <- as.factor(demo$consultant) # Not normal, but equal variance
levene_test(ventfree28 ~ consultant, data = demo)
kruskal_result <- kruskal.test(ventfree28 ~ consultant, data = demo)
print(kruskal_result)
dunn_test(ventfree28 ~ consultant, data = demo)

# Figure 1b. Between RSV status, focus on Non-bacterial
demo_non_bac <- demo%>% 
  filter(consultant != "Bacterial")
demo_bac <- demo %>% 
  filter(consultant == "Bacterial") %>%
  filter(`RASCALS Study ID` %in% microbe$`RASCALS Study ID`)


rsv <- microbe[apply(microbe, 1, function(row) any(grepl("RSV", row))), ] 
# remove 86, 107, 111 because the count doesn't meet the threshold
rsv <- rsv[!rsv$`RASCALS Study ID` %in% c("C086", "C107", "C111"), ]


rsv_id <- rsv$`RASCALS Study ID`

demo_non_bac$group <- ifelse(demo_non_bac $`RASCALS Study ID` %in% rsv_id, "RSV", "Non-RSV")

demo_non_bac$group <- factor(demo_non_bac$group, levels = c("RSV", "Non-RSV"))
ggplot(data = demo_non_bac, aes(x = as.factor(group), y = ventfree28,
                                fill = as.factor(group))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.7, size = 1.5) +
  ylim(0, 28) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("RSV(+)", "RSV(-)")) +
  labs(y = "Vent-free days",
       x = "Group",
       title = "") +
  theme_minimal()

wilcox.test(ventfree28 ~ group, data = demo_non_bac, alternative = "two.sided")


#### 
# First Conclusion: in non bacterial patients (by clinicians), RSV signal ~ more MV 

# Second analysis: severity between RSV in LRTI groups-------------------------------------------------------------
demo$group <- ifelse(demo$`RASCALS Study ID` %in% rsv_id, "RSV", "Non-RSV")
demo$diag <- ifelse(demo$`Primary respiratory problem type` %in% 
                            c("Bronchiolitis", "Lower respiratory tract infection - other",
                              "Pneumonia", "Pneumonitis"),
                          "LRTI", "Non-LRTI")

# Table: RSV vs. LRTI
table_1 <- table(demo$group,
                 demo$diag)
# all the RSV carriers have LRTI!
chisq.test(table_1)


# Figure 2. Boxplot of severity between RSV in LRTI groups
demo_LRTI <- demo %>% filter(diag == "LRTI") 
demo_LRTI$group <- factor(demo_LRTI$group, levels = c("RSV", "Non-RSV"))
demo_nonLRTI <- demo %>% filter(diag != "LRTI") 
ggplot(data = demo_LRTI, aes(x = as.factor(group), y = ventfree28, fill = as.factor(group))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.7, size = 1.5) +
  scale_fill_brewer(palette = "Set2") +
  ylim(0, 28) +
  labs(y = "Vent-free days",
       x = "Group",
       title = "") +
  theme_minimal()

wilcox.test(ventfree28 ~ group, data = demo_LRTI, alternative = "two.sided", var.eqaul = F, exact = FALSE)


# Third analysis: Cytokine and endpoint for Non-RSV, classify LRTI ---------------------------------------------------
cytokine <- read_excel("/Users/zongtaiwu/Documents/RASCALS analysis/Cleaned combined cytokine assays correction.xlsx")
cytokine <- cytokine %>%
  select("id", "sample", "IL-1b", "IL-2", "IL-4", "IL-6", "IL-10", "IL-13", "IP-10", "MCP-1(MCAF)",
         "MIG", "MIP-1b", "RANTES", "TNF-a", "TRAIL", "VEGF", "IL-8", "IFN-g")
cytokine_blood <- cytokine %>%
  filter(sample == "plasma")
cytokine_bal <- cytokine %>%
  filter(sample == "bal")
# ID and diag of those without RSV
col_no_rsv <- demo %>% filter(!`RASCALS Study ID` %in% rsv_id) %>% select(`RASCALS Study ID`, group, diag)
col_all <- demo %>% select(`RASCALS Study ID`, group, diag)


# Merge information about diagnosis to cytokine counts
all <- cytokine %>%
  right_join(demo, by = c("id" = "RASCALS Study ID")) %>%
  relocate(group, diag, .after = id) %>%
  filter(!is.na(sample))


value <- cytokine_bal$`IL-8` #SET
ratio <- cytokine_blood$`IL-4`/cytokine_blood$TRAIL #SET

norm_value <- qnorm((rank(value, na.last="keep")-0.5)/sum(!is.na(value)))
norm_ratio <- qnorm((rank(ratio ,na.last="keep")-0.5)/sum(!is.na(ratio))) #SET

cytokine_bal$norm_value <- norm_value #Change here
cyt_no_rsv <- cytokine_bal %>% # Change here
  right_join(col_no_rsv, by = c("id" = "RASCALS Study ID")) %>%
  relocate(group, diag, .after = id) %>%
  relocate(norm_value, .after = id) %>%
  filter(!is.na(sample))

# In non-rsv individuals, LRTI vs. Non-LRTI
wilcox.test(cyt_no_rsv$norm_value ~ cyt_no_rsv$diag, alternative = "two.sided", var.eqaul = F)

# Figure 3.1
data <- data.frame(norm_value = cyt_no_rsv$norm_value,
                   diag = cyt_no_rsv$diag)
ggplot(data = data, aes(x = as.factor(diag), y = norm_value, fill = as.factor(diag))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.7, size = 1.5) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "IL-4:TRAIL Ratio",
       x = "Diagnosis",
       title = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  theme(legend.position = "right")

# In all individuals, RSV-LRTI, non-RSV-LRTI vs. non-RSV-Non-LRTI
value <- cytokine_blood$`TRAIL` #SET
ratio <- cytokine_blood$`IL-4`/cytokine_blood$TRAIL #SET
norm_value <- qnorm((rank(value, na.last="keep")-0.5)/sum(!is.na(value)))
norm_value <- qnorm((rank(ratio ,na.last="keep")-0.5)/sum(!is.na(ratio))) #SET
cytokine_blood$norm_value <- norm_value #Change here
cyt_all_rsv <- cytokine_blood %>% # Change here
  right_join(col_all, by = c("id" = "RASCALS Study ID")) %>%
  relocate(group, diag, .after = id) %>%
  relocate(norm_value, .after = id) %>%
  mutate(
    group_diag = case_when(
      group == "Non-RSV" & diag == "Non-LRTI" ~ "RSV(-) Non-LRTI",
      group == "Non-RSV" & diag == "LRTI"     ~ "RSV(-) LRTI",
      group == "RSV"     & diag == "LRTI"     ~ "RSV(+) LRTI",
      TRUE                                    ~ NA_character_
    )) %>% 
  relocate(group_diag, .after = id) %>%
  filter(!is.na(sample))

data2 <- data.frame(norm_value = cyt_all_rsv$norm_value,
                    diag = cyt_all_rsv$group_diag)
data2$diag <- factor(data2$diag, levels = c("RSV(+) LRTI", "RSV(-) LRTI", "RSV(-) Non-LRTI"))

# Figure 3.2
library(RColorBrewer)
set1_2 <- brewer.pal(3, "Set1")[1:2]
# make a named vector matching your factor levels
custom_cols <- c(
  "RSV(+) LRTI"     = "orange2",
  "RSV(-) LRTI"     = set1_2[1],
  "RSV(-) Non-LRTI" = set1_2[2]
)

ggplot(data = data2, aes(x = as.factor(diag), y = norm_value, fill = as.factor(diag))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.7, size = 1.5) +
  scale_fill_manual(values = custom_cols) +
  labs(y = "IL-4:TRAIL Level",
       x = "Diagnosis") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+
  theme(legend.position = "none")

dunn_test(norm_value ~ diag, data = data2)


# Fourth analysis: Evaluate model---------------------------------------------------
library(pROC)
library(caret)

value <- cytokine_bal$`IL-8` #SET
ratio <- cytokine_blood$`IL-4`/cytokine_blood$TRAIL #SET

norm_value <- qnorm((rank(value, na.last="keep")-0.5)/sum(!is.na(value)))
norm_value <- qnorm((rank(ratio ,na.last="keep")-0.5)/sum(!is.na(ratio))) #SET

cytokine_blood$norm_value <- norm_value
cyt_no_rsv <- cytokine_blood %>% # Change here
  right_join(col_no_rsv, by = c("id" = "RASCALS Study ID")) %>%
  relocate(group, diag, .after = id) %>%
  relocate(norm_value, .after = id) %>%
  filter(!is.na(sample))
m_data <- data.frame(id = cyt_no_rsv$id,
                     norm_value = cyt_no_rsv$norm_value,
                   diag = cyt_no_rsv$diag)
m_data$diag <- factor(m_data$diag, levels = c("Non-LRTI", "LRTI"))
# 1. Model
model <- glm(diag ~ norm_value, data = m_data, family = "binomial")
# 2. Predict probabilities
m_data$predicted_prob <- predict(model, type = "response")
# 3. Create ROC object
roc_obj <- roc(response = m_data$diag, predictor = m_data$predicted_prob, levels = c("Non-LRTI", "LRTI"), direction = "<")

# 4. Find optimal threshold closest to (0,1)
optimal_threshold_closest <- coords(roc_obj, "best", ret = "threshold", best.method = "closest.topleft")
optimal_coords_closest <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "closest.topleft")

print(optimal_coords_closest)
cat(paste0("Optimal Threshold (Closest to Top-Left): ", round(optimal_threshold_closest, 3), "\n"))
cat(paste0("Sensitivity", round(optimal_coords_closest["sensitivity"], 3), "\n"))
cat(paste0("Specificity", round(optimal_coords_closest["specificity"], 3), "\n\n"))

# 5. Apply the chosen optimal threshold to make predictions and get a confusion matrix
chosen_threshold <- as.numeric(optimal_threshold_closest[1]) # Take the first if multiple are returned
m_data$predicted_class_optimized <- ifelse(m_data$predicted_prob > chosen_threshold, "LRTI", "Non-LRTI")
m_data$predicted_class_optimized <- factor(m_data$predicted_class_optimized, levels = levels(m_data$diag))
conf_matrix_optimized <- confusionMatrix(data = m_data$predicted_class_optimized, reference = m_data$diag, positive = "LRTI")
print(conf_matrix_optimized)


col_diag <- col_no_rsv %>%
  rename(id = `RASCALS Study ID`) %>%
  select(id, true_diag = diag)

bal_m_data_1 <- bal_m_data %>% 
  arrange(id)
bal_m_data_1 <- bal_m_data_1[-c(13, 24, 65), ]
rownames(bal_m_data_1) <- NULL
bal_pred <- bal_m_data_1 %>%
  select(id, bal_prediction = predicted_class_optimized)
bal_model_data <- na.omit(bal_pred %>% full_join(col_diag, by = "id"))

ratio_m_data_1 <- ratio_m_data %>% 
  arrange(id)
ratio_m_data_1 <- ratio_m_data_1[-14, ]
ratio_pred <- ratio_m_data_1 %>%
  select(id, blood_prediction = predicted_class_optimized)
blood_model_data <- na.omit(ratio_pred %>% full_join(col_diag, by = "id"))

# Find those not included in the cytokine model
a <- demo%>%filter(diag == "Non-LRTI")%>%select(`RASCALS Study ID`)
setdiff(a$`RASCALS Study ID`, blood_model_data$id)

b <- demo%>%filter(diag == "LRTI")%>%filter(group == "Non-RSV")%>%select(`RASCALS Study ID`)
setdiff(b$`RASCALS Study ID`, blood_model_data$id)


# 6. Visualize the chosen threshold on the ROC curve
plot(roc_obj, col = "#1c61b6", lwd = 2)

# Plot the point corresponding to the closest-to-top-left threshold
# Extract specificity and sensitivity for plotting
plot_specificity <- as.numeric(optimal_coords_closest["specificity"])
plot_sensitivity <- as.numeric(optimal_coords_closest["sensitivity"])
# Use (1 - specificity) for the x-coordinate (False Positive Rate)
points(x = plot_specificity, y = plot_sensitivity,
       col = "darkgreen", pch = 16, cex = 1.5)
# Add text label for the point, also using 1-specificity for x
text(x = plot_specificity, y = plot_sensitivity,
     labels = paste0("Thres: ", round(chosen_threshold, 3), "\nSens: ", round(plot_sensitivity, 3), "\nSpec: ", round(plot_specificity, 3)),
     pos = 2, col = "darkgreen")