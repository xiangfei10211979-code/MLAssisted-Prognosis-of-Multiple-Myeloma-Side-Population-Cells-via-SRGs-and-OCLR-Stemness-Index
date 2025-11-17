#
# Script: 08_wgcna_module_preservation.R
#
# Description:
# Performs WGCNA module preservation analysis using the SP dataset (augmented)
# as the reference network and validating preservation in multiple GEO/TCGA
# cohorts (Test Networks). The analysis focuses on the combined set of genes
# from the 'blue' and 'green' modules.
#
# Key Steps:
# 1. Load module genes and expression data.
# 2. Augment the SP reference matrix with noise-added samples for robustness.
# 3. Extract common genes and prepare sub-matrices for all test cohorts.
# 4. Run WGCNA::modulePreservation for each test set against the augmented SP reference.
# 5. Extract and report Z.summary and medianRank statistics.
#
# Project:      MM SRGs Prognostic Model
#
# ---

# ---
# Section 1: Configuration and Libraries
# ---
rm(list = ls()); options(stringsAsFactors = FALSE)

# Input Files
MODULE_BLUE_FILE <- "../data/input/module_blue.txt"
MODULE_GREEN_FILE <- "../data/input/module_green.txt"
SP_REFERENCE_FILE <- "../data/input/SP_expression.txt"

# Test Cohorts (Expression Matrices)
DATASETS_FILES <- list(
  TCGA = "../data/input/TCGA_expression.txt",
  GSE24080 = "../data/input/GSE24080_expression.txt",
  GSE57317 = "../data/input/GSE57317_expression.txt",
  GSE19784 = "../data/input/GSE19784_expression.txt"
)

# Output Directory
OUT_DIR <- "../results/wgcna_preservation/"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Parameters (Adjustable)
TARGET_N_REF <- 30         # Target reference sample size (including original)
NOISE_FRACTION <- 0.08      # Noise intensity for augmentation: Noise SD = NOISE_FRACTION * Gene_SD
N_PERMUTATIONS <- 200     # Number of permutations for Z-score calculation (Use >= 1000 for final run)
PRES_RANDOM_SEED <- 123   # Random seed for preservation

# Load Required Packages
required_pkgs <- c("reshape2", "pheatmap", "ggplot2", "MASS")
for (p in required_pkgs) if (!requireNamespace(p, quietly = TRUE)) stop(paste("Package", p, "is required."))

WGCNA_AVAILABLE <- FALSE
if (requireNamespace("WGCNA", quietly = TRUE)) {
  library(WGCNA); enableWGCNAThreads(); WGCNA_AVAILABLE <- TRUE
} else {
  message("[WARN] WGCNA not installed/available. Module preservation will be skipped.")
}
suppressMessages({ library(reshape2); library(pheatmap); library(ggplot2); library(MASS) })

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
  
  # Attempt case-insensitive or stripped matching if direct match fails (conservative approach)
  common2 <- intersect(toupper(mg0), toupper(rn0))
  if (length(common2) > 0) { matched <- rn0[toupper(rn0) %in% common2]; return(list(matched = matched, missing = setdiff(mg0, matched))) }
  
  list(matched = character(0), missing = mg0)
}

read_expr_matrix_auto <- function(path) {
  if (!file.exists(path)) { warning("File not found:", path); return(NULL) }
  dat <- tryCatch(read.table(path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE),
                  error = function(e) { warning("read.table failed:", e$message); return(NULL) })
  if (is.null(dat)) return(NULL)
  mat <- as.matrix(dat)
  # WGCNA expects Samples x Genes (rows x columns)
  if (nrow(mat) >= ncol(mat)) mat <- t(mat)
  colnames(mat) <- clean_genes(colnames(mat)); rownames(mat) <- clean_genes(rownames(mat)); mat
}

# ---
# Section 3: Data Loading, Combining Modules, and Reference Augmentation
# ---

# 1. Load Module Genes (union of blue and green)
blue_genes <- clean_genes(read.table(MODULE_BLUE_FILE, header = FALSE, stringsAsFactors = FALSE)[, 1])
# NOTE: Original script sets green_genes to character(0) later, so we use blue only for now
green_genes_raw <- clean_genes(read.table(MODULE_GREEN_FILE, header = FALSE, stringsAsFactors = FALSE)[, 1])
combined_genes <- unique(c(blue_genes, green_genes_raw))

cat(sprintf("[INFO] Blue genes: %d, Green genes: %d, Union genes: %d\n", 
            length(blue_genes), length(green_genes_raw), length(combined_genes)))

# 2. Load SP Reference Matrix
sp_mat_full <- read_expr_matrix_auto(SP_REFERENCE_FILE)
if (is.null(sp_mat_full)) stop("Failed to read SP reference file.")
cat(paste("[INFO] SP raw dim (Samples x Genes):", paste(dim(sp_mat_full), collapse = " x "), "\n"))

# Match union genes to SP (Reference)
mm_sp <- match_genes_to_matrix(combined_genes, colnames(sp_mat_full))
if (length(mm_sp$matched) == 0) stop("No genes from the combined modules matched the SP reference data.")

# Build Reference Matrix (only with common genes)
ref_mat <- sp_mat_full[, mm_sp$matched, drop = FALSE]

# Remove zero variance columns
var_ref <- apply(ref_mat, 2, var, na.rm = TRUE)
if (any(var_ref == 0)) {
  zgs <- names(var_ref[var_ref == 0]);
  cat(paste("[WARN] SP reference: removed", length(zgs), "zero variance genes.\n"))
  ref_mat <- ref_mat[, var_ref > 0, drop = FALSE]
}
cat(paste("[INFO] Reference matrix dim (Genes for preservation):", paste(dim(ref_mat), collapse = " x "), "\n"))

# 3. Augmentation -> ref_aug (Adding simulated samples)
n_orig <- nrow(ref_mat)
if (n_orig >= TARGET_N_REF) {
  ref_aug <- ref_mat
  cat("[INFO] Original SP samples >= Target. Skipping augmentation.\n")
} else {
  n_add <- TARGET_N_REF - n_orig
  gene_sds <- apply(ref_mat, 2, sd, na.rm = TRUE);
  gene_sds[is.na(gene_sds) | gene_sds == 0] <- median(gene_sds[gene_sds > 0], na.rm = TRUE)
  if (is.na(gene_sds[1])) stop("Gene standard deviations are all zero or NA. Cannot augment.")
  
  ref_aug <- ref_mat
  for (i in seq_len(n_add)) {
    idx <- sample(seq_len(n_orig), 1)
    samp <- ref_mat[idx, ]
    noise <- rnorm(length(samp), mean = 0, sd = NOISE_FRACTION * gene_sds)
    new_sample <- as.numeric(samp) + noise
    ref_aug <- rbind(ref_aug, new_sample)
  }
  rownames(ref_aug) <- c(rownames(ref_mat), paste0("SP_sim_", seq_len(n_add)))
  cat(paste("[INFO] Augmented reference dim:", paste(dim(ref_aug), collapse = " x "), "\n"))
}

# Save augmented matrix and stats
write.table(t(ref_aug), file = file.path(OUT_DIR, "SP_ref_union_augmented_matrix.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# ---
# Section 4: Process Test Datasets and Save Sub-Matrices
# ---

summary_list <- list()
ref_genes <- colnames(ref_aug)
for (ds in names(DATASETS_FILES)) {
  p <- DATASETS_FILES[[ds]]; cat(paste("Processing", ds, "...\n"))
  mat <- read_expr_matrix_auto(p)
  
  if (is.null(mat)) {
    summary_list[[ds]] <- list(status = "file_missing"); next
  }
  
  # Match genes against the augmented reference genes
  mm <- match_genes_to_matrix(ref_genes, colnames(mat))
  matched <- mm$matched; missing <- mm$missing
  
  cat(sprintf("[INFO] %s matched genes: %d, missing: %d\n", ds, length(matched), length(missing)))
  
  if (length(matched) == 0) {
    summary_list[[ds]] <- list(status = "no_matched_genes", matched_genes = 0); next
  }
  
  # Create sub-matrix (Samples x Genes) and ensure gene order matches ref_genes
  mat_sub <- mat[, colnames(mat) %in% matched, drop = FALSE]
  mat_sub <- mat_sub[, ref_genes[ref_genes %in% matched], drop = FALSE]
  
  # Remove zero variance columns in the sub-matrix
  varcols <- apply(mat_sub, 2, var, na.rm = TRUE)
  if (any(varcols == 0)) {
    zg <- names(varcols[varcols == 0]);
    cat(paste("[WARN]", ds, ": removed", length(zg), "zero variance genes.\n"))
    mat_sub <- mat_sub[, varcols > 0, drop = FALSE]
  }
  
  # Save the sub-matrix (Genes x Samples format for general use)
  write.table(t(mat_sub), file = file.path(OUT_DIR, paste0(ds, "_union_matrix.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
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
    original_genes = if (is.null(s$original_genes)) NA else s$original_genes,
    matched_genes = if (is.null(s$matched_genes)) NA else s$matched_genes,
    final_genes = if (is.null(s$final_genes)) NA else s$final_genes,
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_df, file = file.path(OUT_DIR, "datasets_summary_union.csv"), row.names = FALSE)

# ---
# Section 5: Module Color Mapping and Preservation Analysis
# ---

# Define the green_genes to be empty as per the original script's logic
green_genes <- character(0)

# Build module label map for the common genes
module_label_map <- sapply(colnames(ref_aug), function(g) {
  if (g %in% blue_genes) return("blue")
  if (g %in% green_genes_raw) return("green") # If green was intended to be separate
  return("grey")
})

write.table(data.frame(gene = colnames(ref_aug), module = module_label_map),
            file = file.path(OUT_DIR, "union_gene_module_map.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

report_file <- file.path(OUT_DIR, "modulePreservation_union_report.txt")
cat("modulePreservation (Union Modules) Report\n", file = report_file)
cat(paste("Reference (Dataset 1): SP Augmented (N =", nrow(ref_aug), ")\n"), file = report_file, append = TRUE)


if (!WGCNA_AVAILABLE) {
  cat("[WARN] WGCNA not available, skipping modulePreservation.\n", file = report_file, append = TRUE)
} else {
  # Prepare Reference Data (must be matrix with valid R names)
  ref_mat2 <- as.matrix(ref_aug)
  colnames(ref_mat2) <- make.names(colnames(ref_mat2))
  
  for (ds in names(summary_list)) {
    s <- summary_list[[ds]]
    if (s$status != "ok") {
      cat(paste("[INFO]", ds, "skipped due to status:", s$status, "\n"), file = report_file, append = TRUE); next
    }
    
    # Prepare Test Data (must be matrix with valid R names)
    test_sub <- as.matrix(s$mat_sub); colnames(test_sub) <- make.names(colnames(test_sub))
    
    # Identify Common Genes (after renaming to R-compatible names)
    common_genes <- intersect(colnames(ref_mat2), colnames(test_sub))
    
    # Filter for genes present in both, and with non-zero variance (already done above, but re-checked)
    ref_sub <- ref_mat2[, common_genes, drop = FALSE]
    test_sub2 <- test_sub[, common_genes, drop = FALSE]
    
    # Filter again for non-zero variance in both final matrices
    var_check <- apply(ref_sub, 2, var, na.rm = TRUE) > 0 & apply(test_sub2, 2, var, na.rm = TRUE) > 0
    ref_sub <- ref_sub[, var_check, drop = FALSE]
    test_sub2 <- test_sub2[, var_check, drop = FALSE]
    
    common_genes2 <- colnames(ref_sub)
    
    cat(paste("[INFO]", ds, "final common genes for preservation:", length(common_genes2), "\n"), file = report_file, append = TRUE)
    
    if (length(common_genes2) < 30) {
      cat(paste("[WARN]", ds, "common genes after filtering < 30, skipping preservation.\n"), file = report_file, append = TRUE); next
    }
    
    # Construct module color vectors for common genes
    module_colors <- sapply(common_genes2, function(g) {
      g_orig <- sub("\\.", "-", g) # Attempt to revert R naming for lookup
      if (g_orig %in% blue_genes) return("blue")
      if (g_orig %in% green_genes_raw) return("green")
      return("grey")
    })
    
    multiExpr <- list(Reference = list(data = as.data.frame(ref_sub)),
                      Test = list(data = as.data.frame(test_sub2)))
    colorList <- list(Reference = module_colors, Test = module_colors)
    
    mp_file <- file.path(OUT_DIR, paste0(ds, "_modulePreservation_union_raw.rds"))
    mp_err_file <- file.path(OUT_DIR, paste0(ds, "_modulePreservation_union_error.txt"))
    
    mp <- NULL
    # Standard call attempt
    mp <- tryCatch({
      modulePreservation(multiExpr, colorList, referenceNetworks = 1, 
                         nPermutations = N_PERMUTATIONS, randomSeed = PRES_RANDOM_SEED, verbose = 0)
    }, error = function(e) {
      writeLines(e$message, con = mp_err_file)
      cat(paste("[ERROR] modulePreservation standard call failed for", ds, ":", e$message, "\n"), file = report_file, append = TRUE)
      return(NULL)
    })
    
    if (is.null(mp)) {
      # Alternative: use calculateQ=FALSE
      cat(paste("[INFO] Attempting alternate call with calculateQ=FALSE for", ds, "\n"), file = report_file, append = TRUE)
      mp <- tryCatch({
        modulePreservation(multiExpr, colorList, referenceNetworks = 1, 
                           nPermutations = N_PERMUTATIONS, randomSeed = PRES_RANDOM_SEED, 
                           calculateQ = FALSE, verbose = 0)
      }, error = function(e) {
        writeLines(e$message, con = mp_err_file, append = TRUE)
        cat(paste("[ERROR] Alternate call also failed:", e$message, "\n"), file = report_file, append = TRUE)
        return(NULL)
      })
    }
    
    # Save results
    if (!is.null(mp)) {
      saveRDS(mp, mp_file)
      cat(paste("[INFO] modulePreservation saved to:", mp_file, "\n"), file = report_file, append = TRUE)
      
      # Extract Z.summary if available
      if (!is.null(mp$preservation) && is.data.frame(mp$preservation$Z)) {
        presZ <- mp$preservation$Z
        idxZ <- which(presZ$reference.set == 1 & presZ$test.set == 2)
        if (length(idxZ) > 0) {
          z_df <- presZ[idxZ, c("module", "Z.summary"), drop = FALSE]
          write.table(z_df, file = file.path(OUT_DIR, paste0(ds, "_union_Zsummary.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
          cat(paste("[INFO] Z.summary saved for", ds, "\n"), file = report_file, append = TRUE)
        }
      }
    }
  } # end ds loop
} # end WGCNA check


# ---
# Section 6: Report Generation from RDS Files (Extraction and Plotting)
# ---

cat("\n\n--- Section 6: Report Generation from Saved RDS Files ---\n", file = report_file, append = TRUE)

# Function to extract preservation statistics and generate reports/plots
extract_and_report_mp <- function(mp, out_prefix, out_dir) {
  if (is.null(mp)) { warning("Input mp is NULL"); return(NULL) }
  
  # 1) Extract Z summary table
  presZ <- NULL
  if (!is.null(mp$preservation) && !is.null(mp$preservation$Z)) {
    presZ <- as.data.frame(mp$preservation$Z)
    presZ <- presZ[presZ$reference.set == 1 & presZ$test.set == 2, , drop = FALSE]
  }
  
  # 2) Extract medianRank
  medRank <- NULL
  if (!is.null(mp$preservation) && !is.null(mp$preservation$medianRank)) {
    medRank <- as.data.frame(mp$preservation$medianRank)
    medRank <- medRank[medRank$reference.set == 1 & medRank$test.set == 2, , drop = FALSE]
  }
  
  # 3) Build final result table
  res_df <- presZ
  if (!is.null(res_df) && !is.null(medRank)) {
    modcol <- intersect(c("module", "moduleLabel", "moduleName"), colnames(res_df))[1]
    medcol <- intersect(c("module", "moduleLabel", "moduleName"), colnames(medRank))[1]
    
    if (!is.null(modcol) && !is.null(medcol)) {
      # Use dplyr::left_join for cleaner merge if available, otherwise base merge
      if (requireNamespace("dplyr", quietly = TRUE)) {
        res_df <- dplyr::left_join(res_df, medRank[, c(medcol, "medianRank.pres"), drop = FALSE], 
                                   by = setNames(medcol, modcol))
      } else {
        res_df <- merge(res_df, medRank[, c(medcol, "medianRank.pres"), drop = FALSE], 
                        by.x = modcol, by.y = medcol, all.x = TRUE)
      }
    }
  }
  
  # 4) Save table
  table_file <- file.path(out_dir, paste0(out_prefix, "_preservation_table.csv"))
  if (!is.null(res_df)) {
    write.csv(res_df, file = table_file, row.names = FALSE)
  }
  
  # 5) Plot Z.summary barplot
  plot_file <- file.path(out_dir, paste0(out_prefix, "_Zsummary_barplot.pdf"))
  if (!is.null(res_df) && any(grepl("Z.summary", colnames(res_df))) && nrow(res_df) > 0) {
    zcol <- grep("Z.summary", colnames(res_df), value = TRUE)[1]
    modcol_name <- intersect(c("module", "moduleLabel", "moduleName"), colnames(res_df))[1]
    
    if (is.null(modcol_name)) modcol_name <- "module"
    
    # Order for plotting
    ord <- order(res_df[[zcol]], decreasing = TRUE)
    
    pdf(plot_file, width = 8, height = 4)
    par(mar = c(6, 5, 3, 1))
    barplot(res_df[[zcol]][ord],
            names.arg = res_df[[modcol_name]][ord],
            las = 2,
            ylim = c(min(0, min(res_df[[zcol]], na.rm = TRUE) - 1), max(res_df[[zcol]], na.rm = TRUE) + 2),
            main = paste0("Module Preservation Z.summary: ", out_prefix),
            ylab = "Z.summary", col = c("blue", "green", "grey")[match(res_df[[modcol_name]][ord], c("blue", "green", "grey"))])
    abline(h = 10, col = "red", lty = 2) # Strong preservation threshold
    abline(h = 2, col = "orange", lty = 2) # Moderate preservation threshold
    dev.off()
  }
  
  # 6) Save interpretation text
  interp_file <- file.path(out_dir, paste0(out_prefix, "_interpretation.txt"))
  interp_lines <- c(
    paste0("Module Preservation Report: ", out_prefix),
    paste0("Date: ", Sys.time()),
    "",
    "Interpretation thresholds (WGCNA/Langfelder et al.):",
    "  Z.summary > 10 : Strongly preserved",
    "  2 < Z.summary <= 10 : Moderately preserved",
    "  Z.summary <= 2 : No evidence of preservation",
    "",
    "Note: medianRank.pres is a composite measure (lower is better, more conserved).",
    paste("Raw table saved to:", table_file)
  )
  writeLines(interp_lines, con = interp_file)
  
  return(list(table_file = table_file, plot_file = plot_file, interp_file = interp_file))
}

# Batch process saved RDS files
mp_list <- list()
DATASETS <- names(DATASETS_FILES)
for (ds in DATASETS) {
  mp_file <- file.path(OUT_DIR, paste0(ds, "_modulePreservation_union_raw.rds"))
  if (file.exists(mp_file)) {
    mp_list[[ds]] <- readRDS(mp_file)
  } else {
    warning(paste("RDS file not found for", ds, " - skipping report generation."))
  }
}

# Generate reports for successfully loaded results
for (ds in names(mp_list)) {
  mp <- mp_list[[ds]]
  prefix <- paste0(ds, "_vs_SP_union")
  res <- extract_and_report_mp(mp, out_prefix = prefix, out_dir = OUT_DIR)
  cat(paste("[REPORT] Generated report for", ds, ". Table saved to:", res$table_file, "\n"), file = report_file, append = TRUE)
}

cat(paste("\n[DONE] Union-module preservation analysis completed. All outputs are in:", OUT_DIR, "\n"))