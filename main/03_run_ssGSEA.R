#
# Script: 03_run_ssGSEA.R
#
# Description:
# Performs single-sample Gene Set Enrichment Analysis (ssGSEA) using the
# GSVA package on a gene expression matrix. It utilizes custom KEGG gene
# sets, normalizes the scores, calculates differential pathway activity
# between groups (based on TCGA IDs), and generates a heatmap.
#
# Project:      MM SRGs Prognostic Model
#setwd("C:/MMcode/main/")
# Note:
# Assumes input expression matrix has genes (rows) and samples (columns).
# Assumes the input KEGG file is a tab-separated table with 'ID' and
# a '/' delimited 'geneID' column.
#
# ---

# ---
# Section 1: Load Libraries
# ---
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(reshape2)
library(ggpubr)

# ---
# Section 2: Configuration (Define Paths)
# ---
# All paths are relative to the project root (the 'C:\MMcode\' folder).
# This script (03_run_ssGSEA.R) is in the 'main/' subfolder.

# Input files
EXPRESSION_FILE <- "../data/input/expression_matrix_ssgsea.txt" # Expression data
KEGG_GENESET_FILE <- "../data/input/KEGG_gene_sets.txt" # KEGG gene sets list

# Output files
OUTPUT_SCORES_FILE <- "../results/tables/ssgsea_scores_raw.txt"
OUTPUT_HEATMAP_FILE <- "../results/figures/ssgsea_heatmap.pdf"

# ---
# Section 3: Load Data and Prepare Gene Sets
# ---

# Load expression data
data <- read.delim(EXPRESSION_FILE, row.names = 1, header = TRUE)

# Convert to matrix and ensure numeric type
dimnames <- list(rownames(data), colnames(data))
data <- matrix(as.numeric(as.matrix(data)), nrow = nrow(data), dimnames = dimnames)

# Load custom KEGG gene set list
kegg_results <- read.table(KEGG_GENESET_FILE, header = TRUE, sep = "\t")

# Create a list to store gene sets
gene_sets_list <- list()

# Iterate through the KEGG results table to build the gene set list
for (i in 1:nrow(kegg_results)) {
  kegg_id <- kegg_results$ID[i]
  # Split the gene IDs by the '/' delimiter
  gene_list <- unlist(strsplit(kegg_results$geneID[i], "/"))
  
  # Create a unique name (ID and row number)
  unique_name <- paste(kegg_id, i, sep = "_")
  
  # Store the gene list under the unique name
  gene_sets_list[[unique_name]] <- gene_list
}

# Convert the list to a GeneSetCollection object for GSVA
geneSets <- GeneSetCollection(lapply(names(gene_sets_list), function(set_name) {
  GeneSet(geneIds = gene_sets_list[[set_name]], setName = set_name)
}))

# ---
# Section 4: Run ssGSEA Calculation
# ---

# Create ssgseaParam parameter object (modern GSVA API)
ssgsea_params <- ssgseaParam(expr = data, geneSets = geneSets)

# Compute ssGSEA scores
gsvaResult <- gsva(ssgsea_params, verbose = TRUE)

# Normalize scores to the [0, 1] range
normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# Apply normalization row-wise (pathway scores)
gsvaResult <- apply(gsvaResult, 1, normalize)
# Transpose back (pathways as rows, samples as columns)
gsvaResult <- t(gsvaResult)

# Optional: Calculate and add mean values as the last row
mean_values <- colMeans(gsvaResult, na.rm = TRUE)
gsvaResult <- rbind(gsvaResult, Average = mean_values)


# ---
# Section 5: Save Raw Scores
# ---

# Prepare output table including sample IDs in the first row
gsvaOut <- rbind(id = colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file = OUTPUT_SCORES_FILE, sep = "\t", quote = FALSE, col.names = FALSE)
cat(paste("Raw ssGSEA scores saved to:", OUTPUT_SCORES_FILE, "\n"))


