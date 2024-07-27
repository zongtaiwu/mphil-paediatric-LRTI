library(dplyr)
library(magrittr)
library(ggplot2)
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))
# Read host genes variance stabilizing transform values
vsd_mat <- read.csv(cv_vst_csv, row.names=1)

# Key protein host gene counts
key_genes <- read.csv(key_genes_csv) # 105 genes
key_vsd_mat <- vsd_mat %>%
  filter(row.names(vsd_mat) %in% key_genes$Gene.stable.ID) # 72 genes

# Vector mapping ensembl gene IDs to gene symbols
read.csv(host_counts_csv, stringsAsFactors = F) %>%
  {`names<-`(.$gene_symbol, .$X)} ->
  ensg2gene


# Define lasso functions --------------------------------------------------
# Function to extract coefficients from fitted lasso model
lasso_coef_df <- function(mod) {
  # Use the most regularized value of the tuning parameter that is
  # within 1 standard error of the optimum
  coef(mod, s='lambda.1se', gamma=c("gamma.1se"))[,1] %>%
    # Filter to nonzero coefficients
    .[. != 0] %>%
    {data.frame(coef=.)} %>%
    # Convert rowname to column, then grab the gene symbol for the
    # ensembl gene name
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(gene_symbol=ensg2gene[gene])
}

# Function to run the model on the chosen gene set and output results to the out_prefix path.
lasso_main <- function(out_prefix, ...) {
  # Fit lasso logistic regression on the whole dataset
  mod <- glmnet::cv.glmnet(X, y, family='binomial', foldid=classifier_meta$fold, ...)
  
  # Get the nonzero coefficients
  lasso_coef_df(mod) -> coefs
  
  # Save coefficients to file
  coefs %>%
    write.csv(paste0(out_prefix, "coefs.csv"), row.names=F)
  
  list(mod=mod, coefs=coefs)
}

# Lasso on key proteins ---------------------------------------------------
# Regression variables: only key genes
X <- t(key_vsd_mat)
# Labels
y <- classifier_meta$LRTI_adjudication == 'Definite' 

# Run the lasso
lasso_main(key_lasso_results_prefix)

# Lasso on all genes ------------------------------------------------------
# Regression variables: all genes
X <- t(vsd_mat)
# Labels
y <- classifier_meta$LRTI_adjudication == 'Definite' 

# Run the lasso
lasso_main(lasso_results_prefix)