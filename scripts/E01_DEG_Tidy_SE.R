# Source: https://stemangiola.github.io/biocasia2020_tidytranscriptomics/articles/tidytranscriptomics.html

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

# tidyverse-friendly packages
library(plotly)
library(ggrepel)
library(GGally)
library(tidybulk)
library(tidySummarizedExperiment) # we'll load this below to show what it can do
library(org.Hs.eg.db)
library(data.table)

# Use colourblind-friendly colours
friendly_cols <- dittoSeq::dittoColors()

# Set theme
my_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(
          t = 10,
          r = 10,
          b = 10,
          l = 10
        )),
        axis.title.y = element_text(margin = margin(
          t = 10,
          r = 10,
          b = 10,
          l = 10
        )),
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          vjust = 1
        )
      )
  )

# Restore the object
TCGA_LUSC <- readRDS("../00_RAW_DATA/RDS/TCGA_LUSC.rds")

# filtered counts
filtered_counts <- TCGA_LUSC$filtered_counts # If it's a dense matrix

# Create colData and rowData
colData <- TCGA_LUSC$meta_data

# update the sample type
# Recode sample types
colData$sample_type <- recode(colData$sample_type,
                              "Primary Tumor" = "Tumor",
                              "Solid Tissue Normal" = "Normal")

rowData <- DataFrame(gene_id = rownames(filtered_counts))  # Adjust based on your actual row names

# Create the SummarizedExperiment object
TCGA_LUSC_tidy <- SummarizedExperiment(
  assays = list(counts = filtered_counts),
  colData = colData,
  rowData = rowData
)

rm("filtered_counts", "colData", "rowData", "TCGA_LUSC")

# Filter variable transcripts
se.norm.variable = TCGA_LUSC_tidy |> keep_variable()

# TSNE
se.norm.tSNE =
  TCGA_LUSC_tidy |>
  identify_abundant() |>
  reduce_dimensions(method = "tSNE",
                    perplexity = 5,
                    pca_scale = TRUE)

se.norm.tSNE |>
  pivot_sample() |>
  dplyr::select(contains("tSNE"), everything())

se.norm.tSNE |>
  pivot_sample() |>
  ggplot(aes(x = `tSNE1`, y = `tSNE2`, color = sample_type)) + geom_point() + my_theme

# Identify abundant transcripts
se.aggr <- TCGA_LUSC_tidy |> identify_abundant(factor_of_interest = sample_type)

# DE analysis
se.de =
  se.aggr |>
  test_differential_abundance(~ 0 + sample_type, action = "get")
se.de

# DE with contrast
se.de =
  se.aggr |>
  identify_abundant(factor_of_interest = sample_type) |>
  test_differential_abundance(
    ~ 0 + sample_type,
    contrasts = c("sample_typeNormal - sample_typeTumor"),
    action = "get"
  )

# Perform differential abundance analysis using various methods
de_all <- se.aggr |>
  identify_abundant(factor_of_interest = sample_type) |>
  
  # EdgeR QLT
  test_differential_abundance(
    ~ 0 + sample_type,
    contrasts = c("sample_typeNormal - sample_typeTumor"),
    method = "edgeR_quasi_likelihood",
    prefix = "edgerQLT_"
  ) %>%
  
  # EdgeR LRT
  test_differential_abundance(
    ~ 0 + sample_type,
    contrasts = c("sample_typeNormal - sample_typeTumor"),
    method = "edgeR_likelihood_ratio",
    prefix = "edgerLR_"
  ) %>%
  
  # # Limma-voom
  # test_differential_abundance(
  #   ~ 0 + sample_type,
  #   #contrasts = c("sample_typeNormal - sample_typeTumor"),
  #   method = "limma_voom",
  #   prefix = "voom_"
  # ) %>%
  
  # DESeq2
  test_differential_abundance(
    ~ 0 + sample_type,
    contrasts = list(c("sample_type", "Normal", "Tumor")), # Ensure correct format
    method = "deseq2",
    prefix = "deseq2_"
  )

# Comparison of methods

de_all %>%
  pivot_transcript() %>% colnames()

de_all %>%
  pivot_transcript() %>%
  dplyr::select(
    edgerQLT_PValue___sample_typeNormal...sample_typeTumor,
    edgerLR_PValue___sample_typeNormal...sample_typeTumor,
    #voom_P.Value___sample_typeNormal...sample_typeTumor,
    deseq2_pvalue___sample_type.Normal.Tumor,
    .feature
  ) %>%
  ggpairs(1:3)

deg_results <- de_all %>%
  pivot_transcript()

deg_results <- as.data.frame(deg_results)
deg_results[1:4, 1:5]

write.csv(deg_results, file = "./results/TCGA_LUSC_DEG_Tumor__Normal.csv")
