library(dplyr)
library(magrittr)
library(ggplot2)
library(gtools)
library(tidyr)
library(stringr)
library(taxize)
library(taxonomizr)
library(ggsignif)

library(DESeq2)
library(edgeR)
library("EnhancedVolcano")
source("constants.R")

# Filtering ---------------------------------------------------------------
# Read background filtered microbe reports
microbe_reports_bgf <- read.csv(microbe_reports_bgfilter_csv, stringsAsFactors = F)

# Virus filtering
virus_results <- microbe_reports_bgf %>%
  filter(category == "viruses",
         p_adj < 0.05) %>%
  left_join(., y = metadata_df %>% select(sample_name, LRTI_adjudication), by = "sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after = sample_name) %>%
  select(sample_name, LRTI_adjudication, name, tax_id, genus_tax_id, nt_count, nt_rpm)

### BF filtering
# BF pathogens
# Limit to bacterial/fungal hits that passed background filtering
bf_reports <- microbe_reports_bgf %>% 
  filter(category == "bacteria" | category == "eukaryota") %>% 
  filter(p_adj < 0.05)

# Loop through all the samples and group by genus
samplenames <- metadata_df %>% pull(sample_name)
# Create an empty dataframe for results
bf_results_genus <- data.frame()
bf_results_species <- data.frame()

# Grouped by genus and reported as genus name
for(sn in samplenames){
  
  bf_reports %>% 
    dplyr::filter(sample_name == sn,
                  genus_tax_id > 0) -> # require a genus taxid (removes "uncultured bacteria")
    report
  
  report %>% 
    group_by(genus_tax_id) %>% 
    mutate(genus_sum_nt_rpm = sum(nt_rpm)) %>% # calculate genus level nt_rpm
    mutate(genus_sum_nt = sum(nt_count)) %>% #calculate genus level nt_count
    mutate(genus = word(name, 1)) %>% # take the genus name
    ungroup() %>%
    distinct(genus_tax_id, .keep_all = TRUE)->
    genus_report
  
  bf_results_genus <- rbind(bf_results_genus, genus_report)
}

# Grouped by genus and reported as the most dominant species
for(sn in samplenames){
  
  bf_reports %>% 
    dplyr::filter(sample_name == sn,
                  genus_tax_id > 0) -> # require a genus taxid (removes "uncultured bacteria")
    report
  
  report %>% 
    group_by(genus_tax_id) %>% 
    mutate(genus_sum_nt_rpm = sum(nt_rpm)) %>% # calculate genus level nt_rpm
    mutate(genus_sum_nt = sum(nt_count)) %>% #calculate genus level nt_count
    ungroup() ->
    report
  
  report %>%
    group_by(genus_tax_id) %>%          # group by genus
    top_n(1, nt_rpm) %>%                # select the top species within the genus by nt_rpm
    ungroup() ->
    top_species_in_genus_report
  
  bf_results_species <- rbind(bf_results_species, top_species_in_genus_report)
}


# Write out RBM results
bf_results_genus <- bf_results_genus %>%
  left_join(., y = metadata_df %>% select(sample_name, LRTI_adjudication, nonhost_reads), by = "sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after = sample_name) %>%
  select(sample_name, LRTI_adjudication, genus, genus_tax_id,nt_count, nt_rpm, genus_sum_nt_rpm, genus_sum_nt, 
         nonhost_reads)

bf_results_species <- bf_results_species %>%
  left_join(., y = metadata_df %>% select(sample_name, LRTI_adjudication, nonhost_reads), by = "sample_name") %>% 
  dplyr::relocate(LRTI_adjudication, .after = sample_name) %>%
  select(sample_name, LRTI_adjudication, name, tax_id, genus_tax_id,nt_count, nt_rpm, genus_sum_nt_rpm, genus_sum_nt, 
         nonhost_reads)

# Prepare count data: Change count dataset into pivot wide for DE analysis ------------------------------------------------------------------

# VIRUS counts
virus_counts <- virus_results %>%
  filter(LRTI_adjudication == "Definite" | LRTI_adjudication == "No Evidence") %>%
  select(sample_name, name, nt_count) %>%
  mutate(sample_num = as.numeric(str_replace(sample_name, "P", ""))) %>%
  arrange(sample_num) %>%
  select(-sample_num) %>%
  pivot_wider(names_from = sample_name, values_from = nt_count, values_fill = list(nt_count = 0))
virus_counts <- as.data.frame(virus_counts)
row.names(virus_counts) <- virus_counts$name
virus_counts <- virus_counts[, -1]
dim(virus_counts)
# 51 x 125

# Counts for each GENUS, include tax_id because there might be multiple ids for a bracket genus
bf_counts_genus <- bf_results_genus %>%
  filter(LRTI_adjudication == "Definite" | LRTI_adjudication == "No Evidence") %>%
  select(sample_name, genus_tax_id, nt_count) %>%
  mutate(sample_num = as.numeric(str_replace(sample_name, "P", ""))) %>%
  arrange(sample_num) %>%
  select(-sample_num) %>%
  pivot_wider(names_from = sample_name, values_from = nt_count, values_fill = list(nt_count = 0))
bf_counts_genus <- as.data.frame(bf_counts_genus)
row.names(bf_counts_genus) <- bf_counts_genus$genus_tax_id # Use genus_tax_id, need to convert to name later
bf_counts_genus <- bf_counts_genus[, -1]
# 895 x 167
# Filter low count genes
keep <- rowSums(bf_counts_genus > 10) > 2
bf_counts_genus <- bf_counts_genus[keep, ]
bf_counts_genus <- bf_counts_genus[, colSums(bf_counts_genus) > 0]
dim(bf_counts_genus)
# 237 x 167

# Counts for SPECIES that top each genus
bf_counts_species <- bf_results_species %>%
  filter(LRTI_adjudication == "Definite" | LRTI_adjudication == "No Evidence") %>%
  select(sample_name, name, nt_count) %>%
  mutate(sample_num = as.numeric(str_replace(sample_name, "P", ""))) %>%
  arrange(sample_num) %>%
  select(-sample_num) %>%
  pivot_wider(names_from = sample_name, values_from = nt_count, values_fill = list(nt_count = 0))
bf_counts_species <- as.data.frame(bf_counts_species)
row.names(bf_counts_species) <- bf_counts_species$name
bf_counts_species <- bf_counts_species[, -1]
# 1641 x 167
# Filter out low count genes
keep <- rowSums(bf_counts_species > 10) > 2
bf_counts_species <- bf_counts_species[keep, ]
bf_counts_species <- bf_counts_species[, colSums(bf_counts_species) > 0]
dim(bf_counts_species)
# 260 x 166


# Read classifier metadata table, limiting to Definite/No Evidence samples that are used for CV
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))


# DE analysis for virus ---------------------------------------------------

# Remove NA samples
columnData_virus <- classifier_meta %>%
  filter (sample_name %in% colnames(virus_counts)) 

### DESeq
dds_virus <- DESeq2::DESeqDataSetFromMatrix(countData = virus_counts,
                                      colData = columnData_virus,
                                      design = ~ LRTI_adjudication)
dds_virus <- estimateSizeFactors(dds_virus, type = "poscounts") # Address disperse data
# Apply variance stabilizing transform (for t test later)
vst_virus <- assay(
  DESeq2::varianceStabilizingTransformation(dds_virus, fitType = "local")
)
# Log fold change
dds_virus <- DESeq(dds_virus) # Fit a model
res_virus <- results(dds_virus, contrast=c("LRTI_adjudication", "Definite", "No Evidence"), pAdjustMethod = "BH")
DESeq_res_virus <- res_virus[order(res_virus$pvalue),]
DESeq_res_virus <- as.data.frame(DESeq_res_virus) 
write.csv(DESeq_res_virus, virus_DESeq_csv)

# Volcano plot
EnhancedVolcano(toptable = DESeq_res_virus,
                x = "log2FoldChange",
                y = "pvalue",
                lab = rownames(DESeq_res_virus),
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4.0,
                FCcutoff = 2, 
                colAlpha = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                #col=c('black', 'black', 'black', 'red3'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,
                title = "Definite versus No Evidence (DESeq2 results)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
)


### EdgeR
group = factor(columnData$LRTI_adjudication)
dge_virus <- DGEList(counts=virus_counts, 
               group=group)
dge_virus <- calcNormFactors(dge_virus, method="TMMwsp") # Normalize, but we don't filter low counts
# Normalisation and Variance stabilising transform (for t test later)
v_virus <- voom(dge_virus, plot = FALSE)$E
# Log fold change analysis
design <- model.matrix(~0 + group, data=dge_virus$samples)
colnames(design) <- levels(group)
dge_virus <- estimateDisp(dge_virus, design)
# Fit a model
fit <- glmFit(dge_virus, design)
lrt <- glmLRT(fit, contrast=c(1,-1)) 
edgeR_res_virus <- topTags(lrt, n=Inf) 
edgeR_res_virus <- edgeR_res_virus[order(edgeR_res_virus$table$PValue),]
edgeR_res_virus <- as.data.frame(edgeR_res_virus)
write.csv(edgeR_res_virus, virus_edgeR_csv)
# t test
t_test_plots(v_virus, columnData_virus)


EnhancedVolcano(toptable = edgeR_res_virus,
                x = "logFC",
                y = "PValue",
                lab = rownames(edgeR_res_virus),
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4.0,
                FCcutoff = 2, 
                colAlpha = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                #col=c('black', 'black', 'black', 'red3'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,
                title = "Definite versus No Evidence (edgeR results)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
)


# DE analysis for BF (genus name) ------------------------------------------------------
# Prepare the SQLite database
sqlFile <- "nameNode.sqlite"
prepareDatabase(sqlFile)

# Genus level
### DESeq
dds_bf_genus <- DESeq2::DESeqDataSetFromMatrix(countData = bf_counts_genus,
                                            colData = classifier_meta,
                                            design = ~ LRTI_adjudication)
dds_bf_genus <- estimateSizeFactors(dds_bf_genus, type = "poscounts") # Address disperse data
# Apply variance stabilizing transform (for t test later)
vst_bf_genus <- assay(
  DESeq2::varianceStabilizingTransformation(dds_bf_genus, fitType = "local")
)
# Log fold change
dds_bf_genus <- DESeq(dds_bf_genus) # Fit a model
res_bf_genus <- results(dds_bf_genus, contrast=c("LRTI_adjudication", "Definite", "No Evidence"), pAdjustMethod = "BH")
DESeq_res_bf_genus <- res_bf_genus[order(res_bf_genus$pvalue),]
# Turn the taxid into genus
DESeq_res_bf_genus <- as.data.frame(DESeq_res_bf_genus)
DESeq_res_bf_genus <- DESeq_res_bf_genus %>%
  mutate(genus_tax_id = as.numeric(rownames(DESeq_res_bf_genus))) 
genus_tax_ids <- DESeq_res_bf_genus$genus_tax_id
taxonomy <- getTaxonomy(genus_tax_ids, sqlFile)
genus_names <- taxonomy[, "genus"]
DESeq_res_bf_genus$genus <- genus_names
write.csv(DESeq_res_bf_genus, bf_genus_DESeq_csv)

# Volcano plot
# List of species to label
genus_to_label <- c("Haemophilus", "Pasteurella")
DESeq_res_bf_genus$label <- ifelse(DESeq_res_bf_genus$genus %in% genus_to_label, 
                                   DESeq_res_bf_genus$genus, 
                                   NA)

EnhancedVolcano(toptable = DESeq_res_bf_genus,
                x = "log2FoldChange",
                y = "pvalue",
                lab = DESeq_res_bf_genus$genus,
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4.0,
                FCcutoff = 2, 
                colAlpha = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                #col=c('black', 'black', 'black', 'red3'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,
                title = "Definite versus No Evidence (DESeq2 results)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
)

### EdgeR
group=factor(classifier_meta$LRTI_adjudication)
dge_bf_genus <- DGEList(counts=bf_counts_genus, 
                        group=group)
dge_bf_genus <- calcNormFactors(dge_bf_genus, method="TMMwsp")# Normalize,handling zero inflated data
# Normalisation and Variance stabilising transform (for t test later)
v_genus <- voom(dge_bf_genus, plot = FALSE)$E
# Comparing Definite with No evidence
design <- model.matrix(~0 + group, data=dge_bf_genus$samples)
colnames(design) <- levels(group)
dge_bf_genus <- estimateDisp(dge_bf_genus, design)
# Fit a model
fit <- glmFit(dge_bf_genus, design)
lrt <- glmLRT(fit, contrast=c(1,-1)) 
edgeR_res_bf_genus <- topTags(lrt, n=Inf) 
edgeR_res_bf_genus <- edgeR_res_bf_genus[order(edgeR_res_bf_genus$table$PValue),]
# Turn the taxid into genus
edgeR_res_bf_genus <- as.data.frame(edgeR_res_bf_genus)
edgeR_res_bf_genus <- edgeR_res_bf_genus %>%
  mutate(genus_tax_id = as.numeric(rownames(edgeR_res_bf_genus))) 
genus_tax_ids <- edgeR_res_bf_genus$genus_tax_id
taxonomy <- getTaxonomy(genus_tax_ids, sqlFile)
genus_names <- taxonomy[, "genus"]
edgeR_res_bf_genus$genus <- genus_names
write.csv(edgeR_res_bf_genus, bf_genus_edgeR_csv)
# t test
t_test_plots(v_genus, classifier_meta, 0.3)

# Volcano plot
EnhancedVolcano(toptable = edgeR_res_bf_genus,
                x = "logFC",
                y = "PValue",
                lab = edgeR_res_bf_genus$genus,
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4.0,
                FCcutoff = 2, 
                colAlpha = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                #col=c('black', 'black', 'black', 'red3'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,
                title = "Definite versus No Evidence (edgeR results)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
)

# DE analysis for BF (species name) -------------------------------------------------------------------------

# Remove NA samples
columnData_species <- classifier_meta %>%
  filter (sample_name %in% colnames(bf_counts_species)) 

# DESeq
dds_bf_species <- DESeq2::DESeqDataSetFromMatrix(countData = bf_counts_species,
                                               colData = columnData_species,
                                               design = ~ LRTI_adjudication)
dds_bf_species <- estimateSizeFactors(dds_bf_species, type = "poscounts") # Address disperse data
# Apply variance stabilizing transform (for t test later)
vst_bf_species <- assay(
  DESeq2::varianceStabilizingTransformation(dds_bf_species, fitType = "local")
)
# Log fold change
dds_bf_species <- DESeq(dds_bf_species) # Fit a model
res_bf_species <- results(dds_bf_species, contrast=c("LRTI_adjudication", "Definite", "No Evidence"), pAdjustMethod = "BH")
DESeq_res_bf_species <- res_bf_species[order(res_bf_species$pvalue),]
DESeq_res_bf_species <- data.frame(DESeq_res_bf_species)
write.csv(DESeq_res_bf_species, bf_species_DESeq_csv)

#Volcano plot
EnhancedVolcano(toptable = DESeq_res_bf_species,
                x = "log2FoldChange",
                y = "pvalue",
                lab = rownames(DESeq_res_bf_species),
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4.0,
                FCcutoff = 2, 
                colAlpha = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                #col=c('black', 'black', 'black', 'red3'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,
                title = "Definite versus No Evidence (DESeq2 results)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
)

# edgeR
columnData_species <- classifier_meta %>%
  filter (sample_name %in% colnames(bf_counts_species))

group=factor(columnData_species$LRTI_adjudication)
dge_bf_species <- DGEList(counts=bf_counts_species, 
               group=group)
dge_bf_species <- calcNormFactors(dge_bf_species, method="TMMwsp") # Normalize,handling zero inflated data
# Normalisation and Variance stabilising transform (for t test later)
v_species <- voom(dge_bf_species, plot = FALSE)$E
# Comparing Definite with No evidence
design <- model.matrix(~0 + group, data=dge_bf_species$samples)
colnames(design) <- levels(group)
dge_bf_species <- estimateDisp(dge_bf_species, design)
# Fit a model
fit <- glmFit(dge_bf_species, design)
lrt <- glmLRT(fit, contrast=c(1,-1)) 
edgeR_res_bf_species <- topTags(lrt, n=Inf) 
edgeR_res_bf_species <- edgeR_res_bf_species[order(edgeR_res_bf_species$table$PValue),]
edgeR_res_bf_species <- as.data.frame(edgeR_res_bf_species)
write.csv(edgeR_res_bf_species, bf_species_edgeR_csv)

# Volcano plot
EnhancedVolcano(toptable = edgeR_res_bf_species,
                x = "logFC",
                y = "PValue",
                lab = rownames(edgeR_res_bf_species),
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4.0,
                FCcutoff = 2, 
                colAlpha = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                #col=c('black', 'black', 'black', 'red3'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,
                title = "Definite versus No Evidence (edgeR results)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
)
