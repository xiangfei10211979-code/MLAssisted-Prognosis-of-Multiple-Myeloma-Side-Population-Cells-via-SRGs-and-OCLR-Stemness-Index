#
# Script: 08_wgcna_module_preservationM2_AugmentOnly.R
#
# Purpose: Generate an augmented reference expression matrix (SP dataset) using
#          a Principal Component Analysis (PCA)-based resampling method.
#          This augmented matrix (ref_aug) preserves the co-expression structure
#          and is used as the reference network (M2) for subsequent WGCNA
#          module preservation analysis.
#
# Notes: This script only contains the data loading, gene matching, and
#        PCA-based augmentation steps. The WGCNA::modulePreservation
#        analysis code has been explicitly removed as requested.
#
# ---

# ---
# Section 1: Configuration and Libraries
# ---
rm(list = ls()); options(stringsAsFactors = FALSE)

# Input Files (relative to C:/MMcode/main/)
MODULE_BLUE_FILE <- "../data/input/module_blue.txt"
MODULE_GREEN_FILE <- "../data/input/module_green.txt"
SP_REFERENCE_FILE <- "../data/input/SP_expression.txt"

# Test Cohorts (Expression Matrices) - Used only for preparing sub-matrices
DATASETS_FILES <- list(
  TCGA     = "../data/input/TCGA_expression.txt",
  GSE24080 = "../data/input/GSE24080_expression.txt",
  GSE57317 = "../data/input/GSE57317_expression.txt",
  GSE19784 = "../data/input/GSE19784_expression.txt"
)

# Output Directory (relative to C:/MMcode/main/)
OUT_DIR <- "../results/wgcna_preservation/"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Parameters for Augmentation
TARGET_N_REF <- 30        # Target reference sample size (including original)
PRES_RANDOM_SEED <- 123   # Random seed for PCA-based augmentation

# Load Required Packages (WGCNA and non-WGCNA dependencies)
required_pkgs <- c("reshape2", "pheatmap", "ggplot2", "MASS")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cran.rstudio.com/")
  suppressMessages(library(p))
}

# WGCNA check (Only needed for thread management, but modulePreservation is removed)
WGCNA_AVAILABLE <- FALSE
if (requireNamespace("WGCNA", quietly = TRUE)) {
  library(WGCNA); enableWGCNAThreads(); WGCNA_AVAILABLE <- TRUE
}

set.seed(PRES_RANDOM_SEED)

# ---
# Section 2: Helper Functions
# ---

clean_genes <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- gsub("\\s+", "", x); x <- gsub("\\..*$", "", x)
}

match_genes_to_matrix <- function(module_genes, mat_colnames) {
  mg0 <- unique(clean_genes(module_genes)); rn0 <- unique(clean_genes(mat_colnames))
  common <- intersect(mg0, rn0)
  if (length(common) > 0) return(list(matched = common, missing = setdiff(mg0, common)))
  
  common2 <- intersect(toupper(mg0), toupper(rn0))
  if (length(common2) > 0) { matched <- rn0[toupper(rn0) %in% common2]; return(list(matched = matched, missing = setdiff(mg0, matched))) }
  
  mg2 <- gsub("[-_]", "", toupper(mg0)); rn2 <- gsub("[-_]", "", toupper(rn0)); common3 <- intersect(mg2, rn2)
  if (length(common3) > 0) { matched <- rn0[rn2 %in% common3]; return(list(matched = matched, missing = setdiff(mg0, matched))) }
  
  list(matched = character(0), missing = mg0)
}

read_expr_matrix_auto <- function(path) {
  if (!file.exists(path)) { warning("File not found:", path); return(NULL) }
  dat <- tryCatch(read.table(path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE),
                  error = function(e) { warning("read.table failed:", e$message); return(NULL) })
  if (is.null(dat)) return(NULL)
  mat <- as.matrix(dat)
  # WGCNA standard: Samples x Genes (rows x columns)
  if (nrow(mat) >= ncol(mat)) mat <- t(mat)
  colnames(mat) <- clean_genes(colnames(mat)); rownames(mat) <- clean_genes(rownames(mat)); mat
}

# ---
# Section 3: Data Loading and Reference Matrix Preparation
# ---

# 1. Load Module Genes (union of blue and green)
if (!file.exists(MODULE_BLUE_FILE) || !file.exists(MODULE_GREEN_FILE)) stop("One or more module files not found.")
blue_genes <- clean_genes(read.table(MODULE_BLUE_FILE, header = FALSE, stringsAsFactors = FALSE)[, 1])
green_genes <- clean_genes(read.table(MODULE_GREEN_FILE, header = FALSE, stringsAsFactors = FALSE)[, 1])
combined_genes <- sort(unique(c(blue_genes, green_genes)))
cat(sprintf("[INFO] Blue genes: %d, Green genes: %d, Union genes: %d\n", 
            length(blue_genes), length(green_genes), length(combined_genes)))

# 2. Load SP Reference Matrix (Reference)
sp_mat_full <- read_expr_matrix_auto(SP_REFERENCE_FILE)
if (is.null(sp_mat_full)) stop("Failed to read SP reference file.")
cat(paste("[INFO] SP raw dim (Samples x Genes):", paste(dim(sp_mat_full), collapse = " x "), "\n"))

# Match union genes to SP
mm_sp <- match_genes_to_matrix(combined_genes, colnames(sp_mat_full))
if (length(mm_sp$matched) == 0) stop("No genes from the combined modules matched the SP reference data.")

# Build Reference Matrix (Samples x Genes, with only common genes)
ref_mat <- sp_mat_full[, mm_sp$matched, drop = FALSE]

# Remove zero variance columns (essential for PCA)
var_ref <- apply(ref_mat, 2, var, na.rm = TRUE)
if (any(var_ref == 0)) {
  zgs <- names(var_ref[var_ref == 0]); 
  cat(paste("[WARN] SP reference: removed", length(zgs), "zero variance genes.\n"))
  ref_mat <- ref_mat[, var_ref > 0, drop = FALSE]
}
cat(paste("[INFO] Reference matrix dim (Genes for preservation):", paste(dim(ref_mat), collapse = " x "), "\n"))

# ---
# Section 4: PCA-based Augmentation (M2)
# ---

# Assumption: ref_mat is a log2-transformed expression matrix (Samples x Genes)

n_add <- max(0, TARGET_N_REF - nrow(ref_mat))

if (n_add <= 0) {
  ref_aug <- ref_mat
  cat("[INFO] Original SP samples >= Target. Skipping augmentation.\n")
} else {
  set.seed(PRES_RANDOM_SEED)
  cat(paste("[INFO] PCA-based augmentation: Target samples to add =", n_add, "\n"))
  
  # 1) PCA decomposition on the original data (Samples x Genes)
  #    - center = TRUE: Centers the data (subtracts the mean of each gene).
  #    - scale = FALSE: Does not scale by SD, assuming log2 data is already comparable.
  pca <- prcomp(ref_mat, center = TRUE, scale. = FALSE)
  
  # Determine the number of Principal Components (PCs) to keep.
  # Max PCs = min(N_samples - 1, N_genes). We take the minimum effective rank.
  n_pc <- min(nrow(ref_mat) - 1, ncol(ref_mat))
  if (n_pc < 1) stop("Insufficient samples for PCA-based augmentation.")
  cat(paste("[INFO] PCA: Number of PCs retained =", n_pc, "\n"))
  
  # Extract original sample scores (n_orig x n_pc)
  scores_orig <- pca$x[, seq_len(n_pc), drop = FALSE]
  pc_means <- colMeans(scores_orig)
  pc_sds   <- apply(scores_orig, 2, sd)
  pc_sds[pc_sds == 0] <- 1e-6 # Regularization for stability
  
  # 2) Independent Normal Sampling in PC Space
  #    Assumes PC scores are independently and normally distributed.
  scores_new <- matrix(NA, nrow = n_add, ncol = n_pc)
  for (j in seq_len(n_pc)) {
    # Sample new PC scores based on the mean and SD of the original PC scores
    scores_new[, j] <- rnorm(n_add, mean = pc_means[j], sd = pc_sds[j])
  }
  colnames(scores_new) <- colnames(scores_orig)
  
  # 3) Reconstruction to Gene Space
  #    Rotation matrix: (N_genes x N_PCs)
  rotation_sub <- pca$rotation[, seq_len(n_pc), drop = FALSE]
  
  # New samples' expression (centered): X_centered ≈ scores %*% t(rotation)
  new_centered <- scores_new %*% t(rotation_sub)
  
  # Add back the mean center vector (original gene means) to get full expression
  if (!is.null(pca$center)) {
    new_full <- sweep(new_centered, 2, pca$center, "+") # Add means back
  } else {
    new_full <- new_centered
  }
  
  # Combine original and synthetic samples
  ref_aug <- rbind(ref_mat, as.matrix(new_full))
  rownames(ref_aug) <- c(rownames(ref_mat), paste0("SP_sim_pca_", seq_len(n_add)))
  cat(paste("[INFO] ref_aug dim (after PCA-based aug):", paste(dim(ref_aug), collapse = " x "), "\n"))
  
  # 4) Save augmented matrix and stats
  write.table(t(ref_aug), file = file.path(OUT_DIR, "SP_ref_union_augmented_matrix_M2.txt"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.csv(data.frame(gene = colnames(ref_aug),
                       mean = colMeans(ref_aug),
                       sd = apply(ref_aug, 2, sd)),
            file = file.path(OUT_DIR, "SP_ref_union_aug_stats_M2.csv"),
            row.names = FALSE)
  
  cat("[INFO] PCA-based augmentation (M2) complete.\n")
}

# ---
# Section 5: Prepare Test Data and Module Map for Preservation Input
# ---

# 1. Process Test Datasets and Save Sub-Matrices (using ref_aug's gene set)
summary_list <- list()
for (ds in names(DATASETS_FILES)) {
  p <- DATASETS_FILES[[ds]]; cat(paste("Processing test dataset:", ds, "\n"))
  mat <- read_expr_matrix_auto(p)
  
  if (is.null(mat)) { summary_list[[ds]] <- list(status = "file_missing"); next }
  
  # Match genes against the augmented reference genes
  mm <- match_genes_to_matrix(colnames(ref_aug), colnames(mat))
  matched <- mm$matched; missing <- mm$missing
  
  if (length(matched) == 0) { summary_list[[ds]] <- list(status = "no_matched_genes", matched_genes = 0); next }
  
  # Create sub-matrix (Samples x Genes) and ensure gene order matches ref_aug
  mat_sub <- mat[, colnames(mat) %in% matched, drop = FALSE]
  mat_sub <- mat_sub[, colnames(ref_aug)[colnames(ref_aug) %in% matched], drop = FALSE]
  
  varcols <- apply(mat_sub, 2, var, na.rm = TRUE)
  mat_sub <- mat_sub[, varcols > 0, drop = FALSE]
  
  # Save the sub-matrix (Genes x Samples format)
  write.table(t(mat_sub), file = file.path(OUT_DIR, paste0(ds, "_union_matrix_M2.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  summary_list[[ds]] <- list(
    status = "ok", 
    original_samples = nrow(mat), original_genes = ncol(mat), 
    matched_genes = length(matched), final_genes = ncol(mat_sub), 
    mat_sub = mat_sub
  )
}

# Save overall summary table
summary_df <- do.call(rbind, lapply(names(summary_list), function(ds) {
  s <- summary_list[[ds]]
  data.frame(
    dataset = ds, status = s$status,
    original_samples = if (is.null(s$original_samples)) NA else s$original_samples,
    original_genes   = if (is.null(s$original_genes)) NA else s$original_genes,
    matched_genes    = if (is.null(s$matched_genes)) NA else s$matched_genes,
    final_genes      = if (is.null(s$final_genes)) NA else s$final_genes,
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_df, file = file.path(OUT_DIR, "datasets_summary_union_M2.csv"), row.names = FALSE)


# 2. Build Module Color Map (for WGCNA input)
module_label_map <- sapply(colnames(ref_aug), function(g) {
  if (g %in% blue_genes) return("blue")
  if (g %in% green_genes) return("green")
  return("grey")
})

# Save mapping
write.table(data.frame(gene = colnames(ref_aug), module = module_label_map), 
            file = file.path(OUT_DIR, "union_gene_module_map_M2.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n[DONE] PCA-based augmentation (M2) and data preparation complete. Output files saved in:", OUT_DIR, "\n")