#
# Script: 05_drug_sensitivity_oncopredict.R
#
# Description:
# Predicts drug sensitivity (IC50) for patient samples using the
# oncoPredict package and the GDSC2 reference dataset. It then performs
# differential analysis of predicted IC50 between high/low risk groups
# and generates boxplots for significantly different drugs.
#
# Project:      MM SRGs Prognostic Model
#
# ---

# ---
# Section 1: Load Libraries and Configuration
# ---
library(oncoPredict)
library(parallel)
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)

# Input Files
EXPRESSION_FILE <- "../data/input/expression_TCGA_full.txt"
RISK_GROUP_FILE <- "../data/input/risk_group_data.txt"

# Reference Files (GDSC2 from oncoPredict)
GDSC2_EXPR_REF <- "../data/ref/GDSC2_Expr.rds"
GDSC2_RES_REF <- "../data/ref/GDSC2_Res.rds"

# Output Directories
FIG_DIR <- "../results/figures/drug_sensitivity/"
TABLE_DIR <- "../results/tables/"

# Create directories if they don't exist
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# Set random seed for reproducibility
set.seed(999)

# ---
# Section 2: Load Data and Reference Files
# ---

# Load patient expression data (Features x Samples)
data <- read.table(EXPRESSION_FILE, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data <- t(data)

# Convert to numeric matrix
dimnames <- list(rownames(data), colnames(data))
data <- matrix(as.numeric(as.matrix(data)), nrow = nrow(data), dimnames = dimnames)
data <- t(data) # Transpose back to Gene (Row) x Sample (Column) for oncoPredict

# Load reference files
GDSC2_Expr <- readRDS(file = GDSC2_EXPR_REF)
GDSC2_Res <- readRDS(file = GDSC2_RES_REF)
GDSC2_Res <- exp(GDSC2_Res) # Back-transform from log(IC50) to IC50

# ---
# Section 3: Run Drug Sensitivity Prediction (oncoPredict)
# ---

cat("--- Running calcPhenotype (Drug Sensitivity Prediction) ---\n")
# Note: This step might take a long time (e.g., ~30 minutes)
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = data,
              batchCorrect = 'eb', # Empirical Bayes batch correction
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              printOutput = TRUE,
              removeLowVaringGenesFrom = 'rawData')
cat("--- Prediction Complete ---\n")

# Read output
senstivity <- read.csv("calcPhenotype_Output/DrugPredictions.csv", header = TRUE, sep = ",", check.names = FALSE, row.names = 1)

# Clean drug names (remove trailing numbers/batch identifiers if present)
colnames(senstivity) <- gsub("(.*)\\_(\\d+)$", "\\1", colnames(senstivity))

# Clean sample names (keep only the first 12 characters, e.g., TCGA primary ID)
rownames(senstivity) <- substr(rownames(senstivity), 1, 12)

# Log-transform IC50 (standard practice for plotting)
senstivity[is.na(senstivity)] = 0
senstivity <- log2(senstivity + 1)

# Load risk group data
risk <- read.table(RISK_GROUP_FILE, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
risk$risk2 <- factor(risk$risk2, levels = c("low", "high"))

# ---
# Section 4: Merge Data and Differential Analysis
# ---

# Merge sensitivity and risk data
sameSample <- intersect(rownames(risk), rownames(senstivity))
risk <- risk[sameSample, "risk2", drop = FALSE]
senstivity <- senstivity[sameSample, , drop = FALSE]

rt <- cbind(risk, senstivity)
rt <- as.data.frame(rt)

# Ensure all sensitivity columns are numeric
for (i in colnames(rt)[2:ncol(rt)]) {
  rt[[i]] <- as.numeric(rt[[i]])
}

# Find significantly differential drugs (p < 0.001)
sig_drugs <- c()
P_VALUE_CUTOFF <- 0.001

for (i in colnames(rt)[2:ncol(rt)]) {
  # Optionally filter drugs with very low variance (Standard Deviation < 0.05)
  if (sd(rt[, i], na.rm = TRUE) < 0.05) { next }
  
  wilcoxTest <- wilcox.test(rt[, i] ~ rt[, "risk2"])
  pvalue <- wilcoxTest$p.value
  
  if (!is.na(pvalue) && pvalue < P_VALUE_CUTOFF) {
    sig_drugs <- c(sig_drugs, i)
  }
}

# Subset the dataframe to include only the risk group and significant drugs
rt_filtered <- rt[, c("risk2", sig_drugs)]

# Data transformation for ggplot2 (long format)
rt_melted <- melt(rt_filtered, id.vars = "risk2")
colnames(rt_melted) <- c("Risk_Group", "Drug_Name", "Predicted_IC50_log2")

# Set up comparisons for ggpubr (though stat_compare_means handles this for two groups)
my_comparisons <- list(c("low", "high"))

# ---
# Section 5: Plot Boxplots for Significant Drugs
# ---

# Generate Boxplot
boxplot_plot <- ggboxplot(rt_melted, x = "Drug_Name", y = "Predicted_IC50_log2", fill = "Risk_Group",
                          xlab = "",
                          ylab = expression("Predicted IC"[50]~"(Log"[2]~" Value)"),
                          legend.title = "Risk Group",
                          width = 0.8,
                          palette = c("DodgerBlue1", "Firebrick2")) +
  rotate_x_text(45) +
  stat_compare_means(aes(group = Risk_Group),
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")), label = "p.signif") +
  theme(axis.text = element_text(face = "bold", colour = "#441718", size = 10),
        axis.title = element_text(face = "bold", colour = "#441718", size = 12),
        plot.title = element_text(face = "bold", colour = "#441718", size = 14),
        legend.text = element_text(face = "bold", size = 10),
        panel.border = element_rect(fill = NA, color = "#35A79D", size = 1.5, linetype = "solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size = 0.5, linetype = "dotdash"),
        legend.title = element_text(face = "bold", size = 10),
        axis.line = element_blank()
  )

# Save the plot
pdf(file = file.path(FIG_DIR, "Drug_Senstivity_Boxplots.pdf"), width = 20, height = 8)
print(boxplot_plot)
dev.off()
cat(paste("Boxplots of significant drugs saved to:", FIG_DIR, "\n"))


# ---
# Section 6: Generate and Save Statistical Summary Tables
# ---

# Save merged input/sensitivity data
write.csv(rt, file.path(TABLE_DIR, "01_Drug_Sensitivity_FullData.csv"), row.names = TRUE)

# Generate statistical results table
drug_stats <- data.frame(
  Drug = colnames(senstivity),
  p_value = sapply(colnames(senstivity), function(x) {
    wilcox.test(senstivity[, x] ~ risk$risk2)$p.value
  }),
  # Calculate IC50 change: (Median_High - Median_Low) / Median_Low * 100
  IC50_change_percent = sapply(colnames(senstivity), function(x) {
    (median(senstivity[risk$risk2 == "high", x]) - median(senstivity[risk$risk2 == "low", x])) /
      median(senstivity[risk$risk2 == "low", x]) * 100
  })
)

write.csv(drug_stats, file.path(TABLE_DIR, "02_Drug_Sensitivity_Statistics.csv"), row.names = FALSE)
cat(paste("Statistical summary saved to:", TABLE_DIR, "\n"))

# ---
# Section 7: Focused Drug Analysis
# ---

# Define target drugs for focused search (supports partial/fuzzy matching with regex)
target_drugs <- c(
  "PD173074", "RO[_.-]?3306", "Tozasertib", "IAP[_.-]?5620", "AZD5582",
  "Ribociclib", "RAK4[_.-]?4710", "AZD5991", "PAK[_.-]?5339",
  "TAF1[_.-]?5496Dihydrorotenone", "Gallibiscoquinazole", "ULK1[_.-]?498",
  "AZD4547", "GDC0810", "MIRA[_.-]?1", "AZD6482"
)

# Perform fuzzy matching
matched_drugs <- grep(paste(target_drugs, collapse = "|"),
                      drug_stats$Drug, value = TRUE, ignore.case = TRUE)

# Extract and format target drug results
target_results <- drug_stats[drug_stats$Drug %in% matched_drugs, ]

if (nrow(target_results) > 0) {
  target_results$IC50_change_percent <- round(target_results$IC50_change_percent, 1)
  target_results$p_value <- format.pval(target_results$p_value, digits = 3)
  
  cat("\n----- TARGETED DRUG SENSITIVITY CHANGES -----\n")
  print(target_results[order(target_results$IC50_change_percent), ], row.names = FALSE)
} else {
  cat("\n----- WARNING: No target drugs matched. Example of actual drug names (top 20) -----\n")
  print(head(drug_stats$Drug, 20))
}

# Filter significantly increased sensitivity (High Risk is MORE sensitive, IC50 is LOWER)
sig_drugs_more_sensitive <- subset(drug_stats, p_value < 0.05 & IC50_change_percent < 0)

if (nrow(sig_drugs_more_sensitive) > 0) {
  cat("\n----- DRUGS WITH SIGNIFICANTLY INCREASED SENSITIVITY IN HIGH-RISK GROUP -----\n")
  sig_drugs_more_sensitive$IC50_change_percent <- round(sig_drugs_more_sensitive$IC50_change_percent, 1)
  print(sig_drugs_more_sensitive[order(sig_drugs_more_sensitive$IC50_change_percent), ], row.names = FALSE)
} else {
  cat("\n----- No drugs found with significantly increased sensitivity in the High-Risk Group. -----\n")
}