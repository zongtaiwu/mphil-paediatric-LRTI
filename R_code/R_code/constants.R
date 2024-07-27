## Constants

### CHANGE THIS TO POINT TO THE LOCAL REPO PATH ON YOUR SYSTEM
project_root <- "/Users/zongt1/VAP/project-main"

# "raw data" folder
raw_data_dir <- file.path(project_root, "data", "raw")

# "processed data" folder
processed_data_dir <- file.path(project_root, "data", "processed")

# Raw sample metadata
metadata_csv <- file.path(raw_data_dir, "sample_metadata.csv")

# Classifier metadata (with cross-validation folds)
classifier_meta_csv <- file.path(processed_data_dir, "classifier_metadata.csv")

# Raw host gene counts
host_counts_csv <- file.path(raw_data_dir, "host_gene_counts.csv")

# "results" folder
results_dir <- file.path(project_root, "results")

# Host differential analysis results
host_DESeq_csv <- file.path(results_dir, "/DE_analysis/DE_host_DESeq.csv")
host_edgeR_csv <- file.path(results_dir, "/DE_analysis/DE_host_edgeR.csv")

# Host gene VST values for cross-validation
cv_vst_csv <- file.path(processed_data_dir, "cv_vst.csv")

# Key genes for LASSO 
# IMPORT FROM PROTEIN NETWORK ANALYSIS KEY PROTEIN LIST
key_genes_csv <- file.path(processed_data_dir, "LRTI_key_gene_id.csv") 

# Consensus differential genes 
DEG_genes_csv <- file.path(processed_data_dir, "input_gene_id.csv")

# Filename prefix for logistic-lasso results
lasso_results_prefix <- file.path(results_dir, "lasso_")
key_lasso_results_prefix <- file.path(results_dir, "key_lasso_")

# Filename prefix for classifier prediction results
key_LR_results_prefix <- file.path(results_dir, "Pred/key_lasso_LR_")
lasso_LR_results_prefix <- file.path(results_dir, "Pred/lasso_LR_")
LR_results_prefix <- file.path(results_dir, "Pred/LR_")

key_lassoRF_results_prefix <- file.path(results_dir, "Pred/key_lasso_RF_")
lassoRF_results_prefix <- file.path(results_dir, "Pred/lasso_RF_")
RF_results_prefix <- file.path(results_dir, "Pred/RF_")

key_lassoSVM_results_prefix <- file.path(results_dir, "Pred/key_lasso_SVM_")
lassoSVM_results_prefix <- file.path(results_dir, "Pred/lasso_SVM_")
SVM_results_prefix <- file.path(results_dir, "Pred/SVM_")

key_lassoKNN_results_prefix <- file.path(results_dir, "Pred/key_lasso_KNN_")
lassoKNN_results_prefix <- file.path(results_dir, "Pred/lasso_KNN_")
KNN_results_prefix <- file.path(results_dir, "Pred/KNN_")

key_lassoXGBoost_results_prefix <- file.path(results_dir, "Pred/key_lasso_XGBoost_")
lassoXGBoost_results_prefix <- file.path(results_dir, "Pred/lasso_XGBoost_")
XGBoost_results_prefix <- file.path(results_dir, "Pred/XGBoost_")


# Filename prefix for classifier metrics
cv_metrics_prefix <- file.path(results_dir, "Metrics/metrics_")

# SPIA
SPIA_key_DESeq_csv <- file.path(results_dir, "SPIA/deseq_key_spia.csv")
SPIA_deg_DESeq_csv <- file.path(results_dir, "SPIA/deseq_deg_spia.csv")
SPIA_key_edgeR_csv <- file.path(results_dir, "SPIA/edgeR_key_spia.csv")
SPIA_deg_edgeR_csv <- file.path(results_dir, "SPIA/edgeR_deg_spia.csv")

# Consensus results from DEG and key SPIA results
### IMPORT THE TABLE AFTER VENN DIAGRAM AND EXAMINATION
consensus_spia_csv <- file.path(results_dir, "SPIA/consensus_spia.csv")

# Path to raw microbe reports
microbe_reports_csv <- file.path(raw_data_dir, "microbe_reports.csv")

# Path to microbe reports with background filtering stats
microbe_reports_bgfilter_csv <- file.path(processed_data_dir, "microbe_reports_bgfilter.csv")

# Microbe DE analysis results
virus_DESeq_csv <- file.path(results_dir, "/DE_analysis/DE_virus_DESeq.csv")
virus_edgeR_csv <- file.path(results_dir, "/DE_analysis/DE_virus_edgeR.csv")
bf_genus_DESeq_csv <- file.path(results_dir, "/DE_analysis/DE_bf_genus_DESeq.csv")
bf_genus_edgeR_csv <- file.path(results_dir, "/DE_analysis/DE_bf_genus_edgeR.csv")
bf_species_DESeq_csv <- file.path(results_dir, "/DE_analysis/DE_bf_species_DESeq.csv")
bf_species_edgeR_csv <- file.path(results_dir, "/DE_analysis/DE_bf_species_edgeR.csv")

# Correlation test result
RSV_corr_csv <- file.path(results_dir, "RSV_correlation.csv")
