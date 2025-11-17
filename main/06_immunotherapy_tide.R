#
# Script: 06_immunotherapy_tide.R
#
# Description:
# Performs differential analysis of TIDE (Tumor Immune Dysfunction and Exclusion)
# prediction scores between the high and low risk groups, generating a
# violin plot with statistical comparison.
#
# Key Steps:
# 1. Prepare and standardize expression data for TIDE submission (row-wise centering).
# 2. Process external TIDE score results, clean TCGA barcodes, and deduplicate.
# 3. Merge TIDE scores with the risk group data.
# 4. Generate a violin plot comparing TIDE scores between risk groups.
#
# Project:      MM SRGs Prognostic Model
#
# ---

# ---
# Section 1: Load Libraries and Configuration
# ---
library(limma) # Needed for avereps()
library(ggpubr)
library(ggplot2) # ggpubr is built on ggplot2

# Input Files
EXPRESSION_FILE <- "../data/input/expression_TCGA_full.txt"
TIDE_SCORE_FILE <- "../data/input/tide_scores_raw.txt"
RISK_GROUP_FILE <- "../data/input/risk_group_tide.txt"

# Output Directories
FIG_DIR <- "../results/figures/tide/"
TABLE_DIR <- "../results/tables/"

# Create directories if they don't exist
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---
# Section 2: Prepare Expression Data for TIDE Submission
# ---
# TIDE requires row-wise centered (gene-wise mean = 0) expression data (Gene x Sample).

# Load expression data (Assuming Gene x Sample: Rows=Genes, Columns=Samples)
exp <- read.table(EXPRESSION_FILE, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Ensure data is numeric
exp <- data.matrix(exp)

# Perform row-wise centering (standardization for TIDE input)
# Expr[i, j] = Expr[i, j] - mean(Expr[i, ])
Expr_centered <- t(apply(exp, 1, function(x) x - mean(x, na.rm = TRUE)))

# Save the normalized file for potential external TIDE resubmission
write.table(data.frame(ID = rownames(Expr_centered), Expr_centered, check.names = FALSE),
            file.path(TABLE_DIR, 'tcga_normalize_for_tide.txt'), 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("Normalized expression data for TIDE saved to:", file.path(TABLE_DIR, 'tcga_normalize_for_tide.txt'), "\n"))


# ---
# Section 3: Process TIDE Scores and Merge with Risk Group
# ---

# Load TIDE results (Assuming TIDE.txt contains scores with full TCGA barcodes as row names)
tide <- read.table(TIDE_SCORE_FILE, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 1. Filter: Retain only Primary Solid Tumor samples (TCGA barcode position 14-15 < 10)
# This assumes TIDE results include scores for both tumor and normal samples.
tide <- tide[substr(rownames(tide), 14, 15) < "10", , drop = FALSE]

# 2. Simplify Sample ID to primary barcode (12 characters, e.g., TCGA-XX-XXXX)
rownames(tide) <- substr(rownames(tide), 1, 12)

# 3. Deduplicate Samples: Use avereps to compute the average TIDE score if one patient has multiple samples
# This is necessary since the filtering step above might leave multiple tumor samples per patient (e.g., primary vs recurrent).
tide <- avereps(tide)

# Extract only the TIDE score column and ensure it is a dataframe
tide <- as.data.frame(tide[, "TIDE", drop = FALSE])
colnames(tide) <- "TIDE_Score"


# Load risk group data (Assuming one column named 'risk')
risk <- read.table(RISK_GROUP_FILE, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)


# Merge data
sameSample <- intersect(rownames(tide), rownames(risk))
tide <- tide[sameSample, , drop = FALSE]
risk <- risk[sameSample, "risk", drop = FALSE]
data_merged <- cbind(tide, risk)
data_merged$risk <- as.character(data_merged$risk) # Ensure risk is treated as character/factor

# Save merged data
write.csv(data_merged, file.path(TABLE_DIR, "tide_risk_merged_data.csv"), row.names = TRUE)


# ---
# Section 4: Plot TIDE Scores
# ---

# Standardize risk labels for clear plotting
data_merged$risk <- factor(data_merged$risk, 
                           levels = c("low", "high"),
                           labels = c("Low-risk", "High-risk"))

# Set up comparisons
group_levels <- levels(data_merged$risk)
my_comparisons <- list(c(group_levels[1], group_levels[2]))


# Generate Violin Plot
tide_plot <- ggviolin(data_merged, x = "risk", y = "TIDE_Score", fill = "risk",
                      xlab = "", ylab = "TIDE Score",
                      palette = c("DodgerBlue1", "Firebrick2"), # Low-risk (Blue), High-risk (Red)
                      legend.title = "Risk Group",
                      add = "boxplot", add.params = list(fill = "white")) +
  # Add p-value calculated using Wilcoxon test
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Save the plot
pdf(file = file.path(FIG_DIR, "TIDE_score_violin_plot.pdf"), width = 5, height = 4.5)
print(tide_plot)
dev.off()
cat(paste("TIDE score violin plot saved to:", FIG_DIR, "\n"))