####################################################
# Analysis setup

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path

dirpath <- dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# Libraries
library("WebGestaltR")
library("tidyverse")

#################
# User input
#################

# read DEG results
data <- read.table(
  "./results/TCGA_LUSC_DEG_Tumor__Normal.csv",
  header = TRUE,
  row.names = 1,
  sep = ","
)

# select only the edge R results
data <- data %>%
  select(gene_id, matches("QLT"))

colnames(data) <-
  c("gene_symbol", "logFC", "lfcSE", "stat", "P.Value", "adj.P.Val")

data <-
  data %>% dplyr::select(., c("gene_symbol", "logFC", "P.Value", "adj.P.Val"))

# Sort based on the adj p value
data <- data %>% dplyr::arrange(., adj.P.Val)

# Keep unique ENS id
data <- data %>% dplyr::filter(!is.na(gene_symbol))
data <- data %>% dplyr::distinct(., gene_symbol, .keep_all = T)

significant_genes <- data %>%
  as.data.frame() %>%
  filter(adj.P.Val <= 0.01, abs(logFC) >= 2) %>%
  dplyr::select(c("gene_symbol")) %>% column_to_rownames(var = "gene_symbol") %>%
  rownames()

head(significant_genes)
length(significant_genes)

### over-representation test

library(clusterProfiler)

#convert gene symbol to Entrez ID for 

significant_genes_map<- clusterProfiler::bitr(geneID = significant_genes,
                                              fromType="SYMBOL", toType="ENTREZID",
                                              OrgDb="org.Hs.eg.db")

head(significant_genes_map)

## background genes are genes that are detected in the RNAseq experiment 
background_genes<- data %>%
  as.data.frame() %>%
  dplyr::select(c("gene_symbol")) %>% column_to_rownames(var = "gene_symbol") %>%
  rownames()

background_genes_map<- bitr(geneID = background_genes, 
                            fromType="SYMBOL", 
                            toType="ENTREZID",
                            OrgDb="org.Hs.eg.db")
head(background_genes_map)

# Go term enrichment

ego <- enrichGO(gene          = significant_genes_map$ENTREZID,
                universe      = background_genes_map$ENTREZID,
                OrgDb         = "org.Hs.eg.db",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)


# visualization of GO

library(enrichplot)
barplot(ego, showCategory=20) 
dotplot(ego)

# MsigdbR
# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures

library(msigdbr)

m_df <- msigdbr(species = "Homo sapiens")
head(m_df)

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)


table(m_t2g$gs_name)
head(m_t2g)

em <- enricher(significant_genes_map$ENTREZID, TERM2GENE=m_t2g, 
               universe = background_genes_map$ENTREZID )
head(em)