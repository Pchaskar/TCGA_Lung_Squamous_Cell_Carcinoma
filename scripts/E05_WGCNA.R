#WGCNA

# Source: https://github.com/kpatel427/YouTubeTutorials/blob/main/WGCNA.R
# https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0

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

# WGCNA packages
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()          # allow multi-threading (optional)

# Restore the object
TCGA_LUSC <- readRDS("../00_RAW_DATA/RDS/TCGA_LUSC.rds")

# filtered counts
filtered_counts <- TCGA_LUSC$filtered_counts # If it's a dense matrix

# Create colData and rowData
colData <- TCGA_LUSC$meta_data

# remove the TCGA ongect
rm(TCGA_LUSC)

# Create a DESeqDataSet from filtered counts matrix
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = colData,
                              design = ~ 1)# not spcifying model

## remove all genes with counts < 15 in more than 75% of samples (0.75*553=414.75)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 415, ]
nrow(dds75) # 13284 genes

# perform variance stabilization
dds_norm <- vst(dds75)

# get normalized counts
norm.counts <- assay(dds_norm) %>%
  t()
dim(norm.counts)
norm.counts[1:4, 1:4]

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(
  norm.counts,
  powerVector = power,
  networkType = "signed",
  verbose = 5
)

sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(
  norm.counts,
  maxBlockSize = 14000,
  TOMType = "signed",
  power = soft_power,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  randomSeed = 1234,
  verbose = 3
)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  cbind(bwnet$unmergedColors, bwnet$colors),
  c("unmerged", "merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05
)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

# create traits file - binarize categorical variables
traits <- as.data.frame(colData) %>%
  dplyr::mutate(sample_type_bin = ifelse(grepl('Tumor', sample_type), 1, 0)) %>%
  select(85)

# binarize categorical variables

#colData$severity <- factor(colData$severity,
#                           levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))

#severity.out <- binarizeCategoricalColumns(
#  colData$severity,
#  includePairwise = FALSE,
#  includeLevelVsAll = TRUE,
#  minCount = 1
#)


#traits <- cbind(traits, severity.out)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>%
  column_to_rownames(var = 'Row.names')


CorLevelPlot(
  heatmap.data,
  x = names(heatmap.data)[8],
  y = names(heatmap.data)[1:7],
  col = c("blue1", "skyblue", "white", "pink", "red")
)



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>%
  filter(`bwnet$colors` == 'blue') %>%
  rownames()

blueME <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'blue') %>%
  rownames()

turquoiseME <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'turquoise') %>%
  rownames()

yellowME <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'yellow') %>%
  rownames()

# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile.
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

dim(module.membership.measure.pvals)

module.membership.measure.pvals[1:7, 1:10]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$sample_type_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


hub_genes <- gene.signf.corr.pvals %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(25) %>% rownames()


# Using the gene significance you can identify genes that have a high significance for trait of interest
# Using the module membership measures you can identify genes with high module membership in interesting modules.


# Perform Hierarchical Clustering
library("nclust")

scaled_counts <- t(scale(norm.counts))

hist <- coldmap(scaled_counts, method = "ward")

meta_data <- as.data.frame(colData)

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
      blueME[1:25]
    ),
    list(
      turquoiseME[1:25]
    ),
    list(
      yellowME[1:25]
    )
  )
)

# Examine Expression Profiles
module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = labels2colors(bwnet$colors)
)

write_delim(module_df,
            file = "./results/gene_modules_WGCNA.txt",
            delim = "\t")

# pick out a few modules of interest here
modules_of_interest = c("blue","yellow", "turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized <- t(norm.counts)

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

p1 <- submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

# Generate and Export Networks
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = soft_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "./results/edgelist_WGCNA.tsv",
            delim = "\t")

# Perform Hierarchical Clustering
library("nclust")

scaled_counts <- t(scale(t(subexpr)))
dim(scaled_counts)

hist <- coldmap(scaled_counts, method = "ward")

meta_data <- as.data.frame(colData)

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
      blueME[1:25]
    ),
    list(
      turquoiseME[1:25]
    ),
    list(
      yellowME[1:25]
    )
  )
)
