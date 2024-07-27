library(SPIA)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(dplyr)
library(knitr)
library(tidyverse)
library(stringr)
source("constants.R")


# Gene list from network input and key genes
deg_genes <- read.csv(DEG_genes_csv, header = FALSE) # length = 99
colnames(deg_genes) <- c("gene_id")
key_genes <- read.csv(key_genes_csv, header = TRUE) # length = 105
colnames(key_genes) <- c("gene_id", "gene")

# CHOOSE LOGFC VALUES FROM DESEQ2 OR EDGER DATASET
res_Spia <- as.data.frame(res_host_Deseq)
#res_Spia <- as.data.frame(res_host_edgeR)

# key DESeq---------------------------------------------------------------------
res_key <- res_Spia[row.names(res_Spia) %in% key_genes$gene_id, ] # set key_genes
background <- res_Spia

#Map Ensembl to EntrezID
ens.str <- substr(rownames(res_key), 1, 15)
res_key$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res_key$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# Set background
ens.str1 <- substr(rownames(background), 1, 15)
background$entrez <- mapIds(org.Hs.eg.db,
                      keys=ens.str1,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

# DE and ALL input for SPIA
res_key<-res_key[!is.na(res_key$entrez),]
res_key<-res_key[!duplicated(res_key$entrez),]
res_key <- res_key[res_key$padj < .05, ]
DE_VAP = res_key$log2FoldChange
names(DE_VAP) <- as.vector(res_key$entrez)
ALL_VAP <- background$entrez
ALL_VAP=unname(unlist(ALL_VAP))

ALL_T <- background$log2FoldChange
names(ALL_T) <- as.vector(background$entrez)

# Run SPIA
key_spia=spia(de=DE_VAP,all=ALL_VAP,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine = "norminv")

# Extract Entrez-id in each pathway
geneids <- str_extract_all(key_spia$KEGGLINK, "\\d+")

geneids <- sapply(geneids, function(x) paste(x[-1], collapse = " "))

key_spia_selected <- key_spia %>%
  dplyr::select(Name, pNDE, pPERT, pG, pGFdr, pGFWER, Status)%>%
  #filter(pG < 0.06) %>%
  rename(KEGG_pathway = Name)

# Obtain the key proteins in the pathways, added to results
numbers_list <- str_extract_all(key_spia$KEGGLINK, "\\d+")
collapsed_vectors <- sapply(numbers_list, function(x) paste(x[-1], collapse = " "))
deseq_key_spia_selected <- cbind(key_spia_selected, Entrez_id=collapsed_vectors)

write.csv(deseq_key_spia_selected, SPIA_key_DESeq_csv)

# DEG DESeq-------------------------------------------------------------------
res_deg <- res_Spia[row.names(res_Spia) %in% deg_genes$gene_id, ] # set deg genes
background <- res_Spia

#Map Ensembl to EntrezID
ens.str <- substr(rownames(res_deg), 1, 15)
res_deg$symbol <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
res_deg$entrez <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")

# Set background
ens.str1 <- substr(rownames(background), 1, 15)
background$entrez <- mapIds(org.Hs.eg.db,
                            keys=ens.str1,
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")

# DE and ALL input for SPIA
res_deg<-res_deg[!is.na(res_deg$entrez),]
res_deg<-res_deg[!duplicated(res_deg$entrez),]
res_deg <- res_deg[res_deg$padj < .05, ] #differentially expressed genes between VAP and normal
DE_VAP = res_deg$log2FoldChange
names(DE_VAP) <- as.vector(res_deg$entrez)
ALL_VAP <- background$entrez
ALL_VAP=unname(unlist(ALL_VAP))

ALL_T <- background$log2FoldChange
names(ALL_T) <- as.vector(background$entrez)

# Run SPIA
deg_spia=spia(de=DE_VAP,all=ALL_VAP,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine = "norminv")

deg_spia_selected <- deg_spia %>%
  select(Name, pNDE, pPERT, pG, pGFdr, pGFWER, Status)%>%
  #filter(pG < 0.06) %>%
  rename(KEGG_pathway = Name)

# Obtain the key proteins in the pathways, added to results
numbers_list <- str_extract_all(deg_spia$KEGGLINK, "\\d+")
collapsed_vectors <- sapply(numbers_list, function(x) paste(x[-1], collapse = " "))
deseq_deg_spia_selected <- cbind(deg_spia_selected, Entrez_id=collapsed_vectors)

write.csv(deseq_deg_spia_selected, SPIA_deg_DESeq_csv)


# key edgeR---------------------------------------------------------------------
res_key <- res_Spia[row.names(res_Spia) %in% key_genes$gene_id, ] # set key_genes
background <- res_Spia


#Map Ensembl to EntrezID
ens.str <- substr(rownames(res_key), 1, 15)
res_key$symbol <- mapIds(org.Hs.eg.db,
                         keys=ens.str,
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
res_key$entrez <- mapIds(org.Hs.eg.db,
                         keys=ens.str,
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")

# Set background
ens.str1 <- substr(rownames(background), 1, 15)
background$entrez <- mapIds(org.Hs.eg.db,
                            keys=ens.str1,
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")

# DE and ALL input for SPIA
res_key<-res_key[!is.na(res_key$entrez),]
res_key<-res_key[!duplicated(res_key$entrez),]
res_key <- res_key[res_key$FDR < .05, ]
DE_VAP = res_key$logFC
names(DE_VAP) <- as.vector(res_key$entrez)
ALL_VAP <- background$entrez
ALL_VAP=unname(unlist(ALL_VAP))

ALL_T <- background$logFC
names(ALL_T) <- as.vector(background$entrez)

# Run SPIA
key_spia=spia(de=DE_VAP,all=ALL_VAP,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine = "norminv")

# Extract Entrez-id in each pathway
geneids <- str_extract_all(key_spia$KEGGLINK, "\\d+")

geneids <- sapply(geneids, function(x) paste(x[-1], collapse = " "))

key_spia_selected <- key_spia %>%
  dplyr::select(Name, pNDE, pPERT, pG, pGFdr, pGFWER, Status)%>%
  #filter(pG < 0.06) %>%
  rename(KEGG_pathway = Name)

# Obtain the key proteins in the pathways, added to results
numbers_list <- str_extract_all(key_spia$KEGGLINK, "\\d+")
collapsed_vectors <- sapply(numbers_list, function(x) paste(x[-1], collapse = " "))
edgeR_key_spia_selected <- cbind(key_spia_selected, Entrez_id=collapsed_vectors)

write.csv(edgeR_key_spia_selected, SPIA_key_edgeR_csv)

# DEG edgeR-------------------------------------------------------------------

res_deg <- res_Spia[row.names(res_Spia) %in% deg_genes$gene_id, ] # set deg genes
background <- res_Spia

#Map Ensembl to EntrezID
ens.str <- substr(rownames(res_deg), 1, 15)
res_deg$symbol <- mapIds(org.Hs.eg.db,
                         keys=ens.str,
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
res_deg$entrez <- mapIds(org.Hs.eg.db,
                         keys=ens.str,
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")

# Set background
ens.str1 <- substr(rownames(background), 1, 15)
background$entrez <- mapIds(org.Hs.eg.db,
                            keys=ens.str1,
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")

# DE and ALL input for SPIA
res_deg<-res_deg[!is.na(res_deg$entrez),]
res_deg<-res_deg[!duplicated(res_deg$entrez),]
res_deg <- res_deg[res_deg$FDR < .05, ] #differentially expressed genes between VAP and normal
DE_VAP = res_deg$logFC
names(DE_VAP) <- as.vector(res_deg$entrez)
ALL_VAP <- background$entrez
ALL_VAP=unname(unlist(ALL_VAP))

ALL_T <- background$logFC
names(ALL_T) <- as.vector(background$entrez)

# Run SPIA
deg_spia=spia(de=DE_VAP,all=ALL_VAP,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine = "norminv")
#res_spia$Name=substr(res_spia$Name,1,10)
#input_spia[1:20 ,-12]
#subset(input_spia, ID == "05322")

deg_spia_selected <- deg_spia %>%
  select(Name, pNDE, pPERT, pG, pGFdr, pGFWER, Status)%>%
  #filter(pG < 0.06) %>%
  rename(KEGG_pathway = Name)


# Obtain the key proteins in the pathways, added to results
numbers_list <- str_extract_all(deg_spia$KEGGLINK, "\\d+")
collapsed_vectors <- sapply(numbers_list, function(x) paste(x[-1], collapse = " "))
edgeR_deg_spia_selected <- cbind(deg_spia_selected, Entrez_id=collapsed_vectors)

write.csv(edgeR_deg_spia_selected, SPIA_deg_edgeR_csv)