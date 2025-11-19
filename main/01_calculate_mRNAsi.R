#
# Script: 01_calculate_mRNAsi.R
#
# Description:
# This script calculates the mRNA stemness index (mRNAsi) based on the
# OCLR model weights (Malta et al., 2018). It takes a gene expression
# matrix as input, computes the Spearman correlation against the
# model weights, and outputs a normalized [0, 1] mRNAsi score per sample.
#
# Project:      MM SRGs Prognostic Model
#
# Note:
# Assumes input expression matrix has genes (HUGO Symbols) as rows
# and samples as columns.
#
# ---

# ---
# Section 1: Load Libraries
# ---
# 'org.Hs.eg.db' was not used in the final script logic.
library(tidyverse)
library(dplyr)

# ---
# Section 2: Configuration (Define Paths)
# ---
# All paths are relative to the project root (the 'C:\MMcode\' folder).
# This script (01_calculate_mRNAsi.R) is in the 'main/' subfolder.
# Therefore, '..' goes UP one level to the root directory.

# Input files
WEIGHT_FILE           <- "../data/input/Stemness_index.rda"
INPUT_EXPRESSION_FILE <- "../data/input/TCGA_MMRF_expression.txt" # Renamed from original

# Output files
OUTPUT_SCORES_FILE    <- "../results/tables/mRNAsi_scores.txt"
OUTPUT_COMBINED_FILE  <- "../results/tables/expression_with_mRNAsi.txt"

# ---
# Section 3: Load Data and Model
# ---

# Load the OCLR model weights
# This loads an object named 'mRNAsi'
load(WEIGHT_FILE)

# Load the input gene expression data
tpm_data <- read.delim(INPUT_EXPRESSION_FILE, row.names = 1, header = TRUE)

# ---
# Section 4: Process Data and Calculate Scores
# ---

# Find genes common to both the expression data and the mRNAsi model
common_genes <- intersect(rownames(tpm_data), mRNAsi$HUGO)

# Filter the expression data to only these common genes
filtered_data <- tpm_data[common_genes, ]

# Extract weights from the model, ensuring they are in the
# same order as the genes in 'filtered_data'
weights_df <- mRNAsi[match(common_genes, mRNAsi$HUGO), ]
weights <- as.numeric(weights_df$Weight)

# Sanity check
if (length(weights) != nrow(filtered_data)) {
  stop("Error: Weight vector length does not match gene count in filtered data.")
}

# Calculate mRNAsi scores (Spearman correlation) for each sample
mRNAsi_scores <- apply(filtered_data, 2, function(sample_expression) {
  cor(sample_expression, weights, method = "spearman", use = "complete.obs")
})

# Normalize scores to the [0, 1] range
mRNAsi_scores_normalized <- (mRNAsi_scores - min(mRNAsi_scores, na.rm = TRUE)) /
  (max(mRNAsi_scores, na.rm = TRUE) - min(mRNAsi_scores, na.rm = TRUE))

# ---
# Section 5: Save Results
# ---

# 1. Save a file with only the sample names and their scores
mRNAsi_results_df <- data.frame(Sample = names(mRNAsi_scores_normalized),
                                mRNAsi = mRNAsi_scores_normalized,
                                row.names = names(mRNAsi_scores_normalized))

write.table(mRNAsi_results_df,
            file = OUTPUT_SCORES_FILE,
            sep = "\t",
            row.names = FALSE, # Set to FALSE to avoid writing row names
            quote = FALSE)

# 2. Save the expression matrix with the mRNAsi scores added as the last row.
# Note: This file combines normalized expression data with scores.
combined_data <- rbind(filtered_data, mRNAsi = mRNAsi_scores_normalized)

write.table(combined_data,
            file = OUTPUT_COMBINED_FILE,
            sep = "\t",
            row.names = TRUE,
            quote = FALSE)

cat(paste("Success: mRNAsi scores calculated and saved to:\n",
          OUTPUT_SCORES_FILE, "\n",
          OUTPUT_COMBINED_FILE, "\n"))