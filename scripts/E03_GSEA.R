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
  c("gene_symbol",
    "logFC",
    "lfcSE",
    "stat",
    "P.Value",
    "adj.P.Val")

data <-
  data %>% dplyr::select(., c("gene_symbol", "logFC", "P.Value", "adj.P.Val"))

# Sort based on the adj p value
data <- data %>% dplyr::arrange(., adj.P.Val)

# Keep unique ENS id
data <- data %>% dplyr::filter(!is.na(gene_symbol))
data <- data %>% dplyr::distinct(., gene_symbol, .keep_all = T)

# select columns of interest
er <-
  data %>% dplyr::select(., c("gene_symbol", "logFC", "P.Value"))

# Pre-ranking of all genes
er$fcsign <- sign(er$logFC)
er$logP = -log10(er$P.Value)
er$metric = er$logP / er$fcsign

er <- er[complete.cases(er), ]

er <- dplyr::filter(er, metric != Inf)

ranks <- er[, c("gene_symbol", "metric")]

ranks <- ranks %>% filter(., !is.na(metric))
ranks <- ranks %>% filter(., !is.na(gene_symbol))

# Write rank file
write.table(
  ranks,
  file = paste("./", "gsea", "/", "deg_rank", ".rnk", sep = ""),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)

# Directory to store GSEA results
outputDirectory <- paste("./", "gsea", sep = "")
projectName <- "TCGA_LUSC_GSEA"

# Rank file
rankFile <-
  paste("./", "gsea", "/", "deg_rank", ".rnk", sep = "") #give path

print(projectName)

# list all gene sets available
# listGeneSet()

# Run GSEA Analysis
###################
suppressWarnings(
  enrichResult <-
    WebGestaltR(
      enrichMethod = "GSEA",
      organism = "hsapiens",
      enrichDatabase = c("pathway_KEGG", "pathway_Panther", "pathway_Reactome"),
      interestGeneFile = rankFile,
      interestGeneType = "genesymbol",
      sigMethod = "top",
      topThr = 100,
      minNum = 5,
      nThreads = 10,
      outputDirectory = outputDirectory,
      projectName = projectName,
      gseaPlotFormat = "png"
    )
)

