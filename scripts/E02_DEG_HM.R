# HeatMap of DEGs

# Set options for future globals and memory limit
options(future.globals.maxSize = 10000 * 1024 ^ 2)  # Set maximum size for future globals
library(BiocParallel)

# Register parallel backend
register(MulticoreParam(12))  # Use 12 cores

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path

dirpath <- dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# tidyverse core packages
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(nclust)

# Restore the object
TCGA_LUSC <- readRDS("../00_RAW_DATA/RDS/TCGA_LUSC.rds")

# filtered counts
vst_counts <- TCGA_LUSC$vst_counts # If it's a dense matrix

# Create colData and rowData
meta_data <- TCGA_LUSC$meta_data

# update the sample type
# Recode sample types
meta_data$sample_type <- recode(meta_data$sample_type,
                                "Primary Tumor" = "Tumor",
                                "Solid Tissue Normal" = "Normal")

# remove TCGA_LUSC
rm(TCGA_LUSC)

##--
# DEGS
DEG <- read.table(
  "./results/TCGA_LUSC_DEG_Tumor__Normal.csv",
  header = TRUE,
  row.names = 1,
  sep = ","
)

# Filter rows where any "FDR" or "padj" column has value < 0.05
filtered_DEG <- DEG %>%
  filter(if_any(contains(c("FDR", "padj")), ~ . < 0.01))

# Filter rows where any "logFC" or "log2FoldChange" column has absolute value > 0.58
filtered_DEG <- filtered_DEG %>%
  filter(if_any(contains(c(
    "logFC", "log2FoldChange"
  )), ~ abs(.) > 0.58))

head(filtered_DEG)

# DEG counts
# Subset the matrix based on matching row names
deg_counts <- vst_counts[rownames(vst_counts) %in% filtered_DEG$gene_id, ]

# Clustering
scaled_counts <- t(scale(t(deg_counts)))

hist <- coldmap(scaled_counts, method = "ward")

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
      c(
        "CD3[DEG]",
        "TTF1",
        "PDCD1",
        "CD274",
        "KRT7",
        "TP63",
        "TP73",
        "EGFR",
        "ALK",
        "ROS1",
        "BRAF",
        "KRAS",
        "HER2",
        "MET",
        "RET",
        "NTRK",
        "MLH1",
        "MSH2",
        "MSH6",
        "PMS2"
      )
    ),
    list(
      c("FLRT3", "PPP2R2C", "MMP3", "MMP12", "CAPN8", "FILIP1", "SPP1")
    ),
    list(c("KLK6", "MUC22", "CSN1S1")),
    list(c("DPPA", "TTTY16", "TRIM58", "HKDC1")),
    list(c("ITLN1", "RS1", "ANGPT4", "CD300LG", "KRT31", "S100A7", "KRT14"))
  )
)
