# script to run survival analysis using TCGA data
# setwd("~/Desktop/demo/survivalAnalysis")
# source: https://github.com/kpatel427/YouTubeTutorials/blob/main/KM_survival.R

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path

dirpath <- dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# Libraries
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

# getting clinical data for TCGA-lusc cohort -------------------
clinical_lusc <- GDCquery_clinic("TCGA-LUSC")

any(colnames(clinical_lusc) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_lusc) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_lusc[,c(9,39,45)]


# looking at some variables associated with survival 
table(clinical_lusc$vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# change certain values the way they are encoded
clinical_lusc$deceased <- ifelse(clinical_lusc$vital_status == "Alive", FALSE, TRUE)
table(clinical_lusc$deceased)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_lusc$overall_survival <- ifelse(clinical_lusc$vital_status == "Alive",
                                         clinical_lusc$days_to_last_follow_up,
                                         clinical_lusc$days_to_death)

# get gene expression data -----------

# build a query to get gene expression data for entire cohort
query_lusc_all = GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

output_lusc <- getResults(query_lusc_all)

# get 20 primary tissue sample barcodes
tumor <- output_lusc$cases

# OR
tumor <- output_lusc[output_lusc$sample_type == "Primary Tumor", "cases"]
tumor

# # get gene expression data from 20 primary tumors 
query_lusc <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor)

# download data
GDCdownload(query_lusc)

# get counts
tcga_lusc_data <- GDCprepare(query_lusc, summarizedExperiment = TRUE)
lusc_matrix <- assay(tcga_lusc_data, "unstranded")
lusc_matrix[1:10,1:10]

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_lusc_data))
coldata <- as.data.frame(colData(tcga_lusc_data))

# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = lusc_matrix,
                              colData = coldata,
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# vst 
vsd <- vst(dds, blind=FALSE)
lusc_matrix_vst <- assay(vsd)
lusc_matrix_vst[1:10,1:10]

# Get data for TP53 gene and add gene metadata information to it -------------
lusc_goi <- lusc_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "CDKN2A")


# get median value
median_value <- median(lusc_goi$counts)

# denote which cases have higher or lower expression than median count
lusc_goi$strata <- ifelse(lusc_goi$counts >= median_value, "HIGH", "LOW")
table(lusc_goi$strata )

# Add clinical information to lusc_goi
lusc_goi$case_id <- gsub('-01.*', '', lusc_goi$case_id)
lusc_goi <- merge(lusc_goi, clinical_lusc, by.x = 'case_id', by.y = 'submitter_id')

# fitting survival curve -----------
# Fits a Kaplan-Meier survival model based on high and low expression of the gene
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = lusc_goi)
fit
ggsurvplot(fit,
           data = lusc_goi,
           pval = T,
           risk.table = T)

# Performs a log-rank test to statistically compare the survival differences between the high and low expression groups.
fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = lusc_goi)