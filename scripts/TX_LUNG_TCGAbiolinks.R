# Install TCGAbiolinks package
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
# BiocManager::install("EDASeq")

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path

dirpath <- dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# Load necessary libraries
library(TCGAbiolinks)  # For accessing and downloading TCGA data
library(SummarizedExperiment)  # For handling TCGA data structures
library(dplyr)  # For data manipulation

# Explore TCGA
# List availabple projects at TCGA
GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])

# TCGA_READdetails
TCGAbiolinks:::getProjectSummary("TCGA-LUSC")

# Define the cancer type. Use "LUSC" for Lung Squamous Cell Carcinoma
# LUAD: Lung adenocarcinoma
# LUSC: Lung squamous cell carcinoma
cancer_type <- "TCGA-LUSC"

# Query for RNA-seq data (HTSeq - Counts)
rna_query <- GDCquery(
  project = cancer_type,
  # Cancer project (LUSC)
  data.category = "Transcriptome Profiling",
  # RNA-seq data category
  data.type = "Gene Expression Quantification",
  # RNA-seq type
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts"# Workflow type (raw counts)
)

# Download the RNA-seq data using api
GDCdownload(rna_query, method = "api", directory = "./data/")

# Prepare the RNA-seq data as a SummarizedExperiment object
rna_data <- GDCprepare(rna_query, directory = "./data/")

# Query and download clinical data
clinical_query <- GDCquery(
  project = cancer_type,
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "bcr xml"
)
GDCdownload(clinical_query, method = 'api', directory = "./data/")
clinical_data <- GDCprepare_clinic(clinical_query,
                                   clinical.info = "patient",
                                   directory = "./data/")

# Explore RNA-seq data
# View the list of assays available in the data (e.g., stranded, unstranded counts)
assayNames(rna_data)
str(rna_data@assays)

# If you find both unstranded and stranded counts, you can specify the appropriate one
rna_matrix <- assay(rna_data, "stranded_second")

dim(rna_matrix)
rna_matrix[1:4, 1:4]

# Explore clinical data
dim(clinical_data)
head(clinical_data)
summary(clinical_data)

clinical_data[1:4, 1:4]

dim(rna_matrix)
dim(clinical_data)

# Data Filtering
# Extract gene annotations
gene_info <- rowRanges(rna_data)

# Filter for protein-coding genes
if ("gene_type" %in% colnames(mcols(gene_info))) {
  protein_coding_genes <- mcols(gene_info)$gene_type == "protein_coding"
} else {
  stop("Column 'gene_type' not found in gene annotations.")
}

# Subset to retain only protein-coding genes
rna_data_filtered <- rna_data[protein_coding_genes, ]

# Extract gene annotations (including HGNC symbols)
gene_annotations <- rowRanges(rna_data_filtered)

# Extract HGNC symbols
hgnc_symbols <- mcols(gene_annotations)$gene_name

# Check valid gene names or annotations
# Filter out rows where the HGNC symbol is missing (if you want to remove genes without HGNC symbols)
valid_genes <- which(hgnc_symbols != "" & !is.na(hgnc_symbols))

# Final filtering step
rna_data_filtered <- rna_data_filtered[valid_genes, ]
filtered_hgnc_symbols <- hgnc_symbols[valid_genes]

# Replace rownames of filtered counts matrix with filtered HGNC symbols
rownames(rna_data_filtered) <- filtered_hgnc_symbols

# Confirm the filtering
dim(rna_data_filtered)
head(rowRanges(rna_data_filtered))

# Filtered raw counts
filtered_counts <- assay(rna_data_filtered, "stranded_second")
dim(filtered_counts)
filtered_counts[1:4, 1:4]

# Filter lowly expressed genes (at least 10 counts in 20% of samples)
min_counts <- 10
min_samples <- floor(ncol(filtered_counts) * 0.20)
filtered_counts <- filtered_counts[rowSums(filtered_counts >= min_counts) >= min_samples, ]

dim(filtered_counts)

# Data Normalization
# Load DESeq2 package
library(DESeq2)

# Create a DESeqDataSet from filtered counts matrix
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = colData(rna_data_filtered),
  design = ~ 1
)

# Estimate size factors (normalization factors)
dds <- estimateSizeFactors(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Check the first few rows of normalized counts
head(normalized_counts)

# Perform VST transformation
vst_transformed_data <- vst(dds, blind = FALSE)

# Extract VST-transformed counts matrix
vst_counts <- assay(vst_transformed_data)

# Check the first few rows of VST-transformed counts
head(vst_counts)

# Exploratory data analysis
# Perform PCA on the rlog-transformed data
plotPCA(vst_transformed_data, intgroup = "sample_type")

# Perform Hierarchical Clustering
library("nclust")

scaled_counts <- t(scale(t(vst_counts)))

hist <- coldmap(scaled_counts, method = "ward")

meta_data <- vst_transformed_data@colData

table(meta_data$sample_type)
table(meta_data$tissue_type)

coldmap(
  scaled_counts,
  clust = hist,
  saturation = TRUE,
  ctag = make_tag(
    meta_data,
    varnames = c("sample_type"),
    cols = c("violet")
  ),
  ctag.space = 1.5,
  rmarg = 1,
  rlab = list(
    list(
      c("CD3[DEG]")
    )
  )
)

# Save R object
# Create a list with all important variables
important_variables <- list(
  rna_data_filtered = rna_data_filtered,
  clinical_data = clinical_data,
  filtered_counts = filtered_counts,
  dds = dds,
  normalized_counts = normalized_counts,
  vst_counts = vst_counts,
  scaled_counts = scaled_counts,
  meta_data = meta_data,
  hist = hist
)

# Save the list as a single RDS file
saveRDS(important_variables, file = "./RDS/TCGA_LUSC.rds")

rm(list = ls())
gc()
.rs.restartR()