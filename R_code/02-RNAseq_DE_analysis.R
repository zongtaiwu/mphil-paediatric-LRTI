library(dplyr)
library(magrittr)
library(DESeq2)
library(edgeR)
library(ggplot2)
library("EnhancedVolcano")
source("constants.R")

# Read classifier metadata table, limiting to Definite/No Evidence samples
classifier_meta <- read.csv(classifier_meta_csv, stringsAsFactors = F) %>% 
  filter(LRTI_adjudication %in% c("Definite", "No Evidence"))

# Read host counts, omitting the first column (gene symbol)
host_counts <- read.csv(host_counts_csv, stringsAsFactors = F, row.names = 1)[,-1] #dim: 36816, 261
# Filter host counts to the samples with known LRTI status
host_counts[, classifier_meta$sample_name] ->
  cv_host_counts #dim: 36816, 167

# Filter for genes with >10 counts in >20% of samples 
keep <- rowSums(cv_host_counts > 10) > 0.2*ncol(cv_host_counts)
cv_host_counts <- cv_host_counts[keep, ] #dim: 13323, 167

# Vector mapping ensembl gene IDs to gene symbols (For volcano plot and later feature selection)
read.csv(host_counts_csv, stringsAsFactors = F) %>%
  {`names<-`(.$gene_symbol, .$X)} ->
  ensg2gene
ensg2gene_unqiue <- ensg2gene[!duplicated(ensg2gene)]

# DE analysis using DESeq -----------------------------------------------
dds_host <- DESeq2::DESeqDataSetFromMatrix(countData = cv_host_counts,
                                           colData = classifier_meta,
                                           design = ~ LRTI_adjudication)
dds_host <- DESeq(dds_host)
res_host_Deseq <- results(dds_host, contrast=c("LRTI_adjudication", "Definite", "No Evidence"), alpha=0.05)
res_host_Deseq <- res_host_Deseq[order(res_host_Deseq$pvalue),]

# Significant DE genes
sig_genes_deseq <- subset(res_host_Deseq, padj < 0.05 & pvalue < 0.05) 
sig_genes_deseq <- sig_genes_deseq %>%
  as.data.frame() %>%
  filter(log2FoldChange > 2 | log2FoldChange < -2)

# Output the table
write.csv(sig_genes_deseq, host_DESeq_csv)

# DE analysis using edgeR -------------------------------------------------
group <- factor(classifier_meta$LRTI_adjudication) 
dge <- DGEList(counts=cv_host_counts, group=group)
keep <- filterByExpr(dge) # Filter out lowly expressed genes
dge <- dge[keep, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM") # Normalize
design <- model.matrix(~0 + group, data=dge$samples) # Comparing Definite with No evidence
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design)

# Fit a model
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, contrast=c(1,-1)) 
res_host_edgeR <- topTags(lrt, n=Inf)

# Significantly DE genes
sig_genes_edgeR <- res_host_edgeR$table[res_host_edgeR$table$PValue < 0.05 & res_host_edgeR$table$FDR < 0.05, ] # 4864
sig_genes_edgeR <- sig_genes_edgeR %>%
  filter(logFC > 2 | logFC < -2) # 144 in total

# Output the table
write.csv(sig_genes_edgeR, host_edgeR_csv)

# Apply variance stabilizing transform using DESeq-------------------------------------------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cv_host_counts,
                                      colData = classifier_meta,
                                      design = ~1)
vsd <- DESeq2::varianceStabilizingTransformation(dds)
# Round transformed values to 2 decimal digits
vsd_mat <- SummarizedExperiment::assay(vsd) %>% 
  round(., digits=2)

# Save VST output
write.csv(vsd_mat, cv_vst_csv)

# Visualisation (volcano plots) ------------------------------------------------------------------------

# DESeq
# Convert gene ID to gene names
symbols <- ensg2gene[rownames(res_host_Deseq)]
volcano_host_Deseq <- res_host_Deseq
rownames(volcano_host_Deseq) <- symbols

EnhancedVolcano(toptable = volcano_host_Deseq,
                x = "log2FoldChange",
                y = "pvalue",
                lab = rownames(volcano_host_Deseq),
                pCutoff = 0.05,
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
# Map unique gene symbols to row names
symbols <- ensg2gene_unique[rownames(res_host_edgeR)]
volcano_host_edgeR <- as.data.frame(res_host_edgeR)
# Apply unique symbols to row names, ensuring lengths match
valid_indices <- !is.na(symbols)
volcano_host_edgeR <- volcano_host_edgeR[valid_indices, ]
symbols <- symbols[valid_indices]
# Now set the row names of the filtered data frame
rownames(volcano_host_edgeR) <- symbols

EnhancedVolcano(toptable = volcano_host_edgeR,
                x = "logFC",
                y = "PValue",
                lab = rownames(volcano_host_edgeR),
                pCutoff = 0.05,
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

