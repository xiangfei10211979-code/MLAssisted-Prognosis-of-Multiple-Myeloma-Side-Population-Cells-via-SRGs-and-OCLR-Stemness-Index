#
# Script: 07_prognostic_validation.R
#
# Description:
# Performs comprehensive prognostic validation for multiple markers using
# single-factor Cox regression, time-dependent AUC (timeROC), non-proportional
# hazards modeling, bootstrap CI estimation, and survival subgroup analysis.
# This script focuses on short-term prognostic evaluation (90, 180, 365 days).
#
# Project:      MM SRGs Prognostic Model Validation (GSE24080)
#
# ---

# ---
# Section 1: Configuration and Libraries
# ---
library(survival)
library(dplyr)
library(timeROC)
library(boot)
library(grDevices) # For plotting colors

# Input/Output Configuration
INPUT_FILE <- "../data/input/GSE24080_prognostic_data.txt"
REPORT_OUTPUT <- "../results/reports/prognosis/GSE24080_singlefactors_report.txt"
TIME_ROC_OUTPUT <- "../results/tables/GSE24080_timeROC_shortwindow.txt"
FIG_PREFIX <- "../results/figures/prognosis/GSE24080_plot_"

# Analysis Parameters
VARS_TO_ANALYZE <- c("riskex1", "Average", "hsa05222_6", "mRNAsi")
TIMES_EVAL <- c(90, 180, 365) # Short-term windows (days)
R_BOOTSTRAP <- 500 # Number of bootstrap replicates

# Create output directories
dir.create(dirname(REPORT_OUTPUT), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(TIME_ROC_OUTPUT), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(FIG_PREFIX), showWarnings = FALSE, recursive = TRUE)

# Set random seed for reproducible bootstrap/timeROC
set.seed(123)

# ---
# Section 2: Data Loading and Preparation
# ---

data <- read.delim(INPUT_FILE, header = TRUE, stringsAsFactors = FALSE)

# Check required columns
if (!all(c("time", "state") %in% colnames(data))) {
  stop("Input data must contain 'time' and 'state' columns.")
}
miss_vars <- VARS_TO_ANALYZE[!VARS_TO_ANALYZE %in% colnames(data)]
if (length(miss_vars) > 0) {
  stop(paste("Missing analysis variables:", paste(miss_vars, collapse = ", ")))
}

# Generate Z-score versions for numerical stability in time-dependent models/timeROC
for (v in VARS_TO_ANALYZE) {
  zname <- paste0(v, "_z")
  data[[zname]] <- as.numeric(scale(data[[v]]))
}

# Create survival object (1 = Event, 0 = Censored)
surv_obj <- Surv(time = data$time, event = data$state)


# ---
# Section 3: Core Analysis Function (Univariate Cox, C-index, timeROC, Bootstrap)
# ---

cox_one <- function(var, data, surv_obj, times_eval, R_boot) {
  # --- 1. Cox Regression ---
  fmla <- as.formula(paste("surv_obj ~", var))
  fit <- tryCatch(coxph(fmla, data = data), error = function(e) NULL)
  
  if (is.null(fit)) {
    warning(paste("Cox regression failed for variable:", var))
    return(list(var = var, coef = NA, p = NA, C_index = NA, AUCs = rep(NA, length(times_eval)), boot_hr = c(NA, NA)))
  }
  
  s <- summary(fit)
  coefv <- s$coefficients[1, "coef"]
  sev <- s$coefficients[1, "se(coef)"]
  z <- s$coefficients[1, "z"]
  p <- s$coefficients[1, "Pr(>|z|)"]
  
  # Calculate 95% CI for HR
  lower_coef <- coefv - 1.96 * sev
  upper_coef <- coefv + 1.96 * sev
  
  # HR (per 1 unit)
  HR_1 <- exp(coefv); HR_1_ci <- exp(c(lower_coef, upper_coef))
  # HR (per 0.1 unit)
  HR_0.1 <- exp(coefv * 0.1); HR_0.1_ci <- exp(c(lower_coef, upper_coef) * 0.1)
  
  # C-index
  conc_obj <- survConcordance(fmla, data = data)
  c_index <- if (!is.null(conc_obj$concordance)) conc_obj$concordance else NA
  
  # --- 2. Time-Dependent AUC (using Z-score version) ---
  aucs <- rep(NA, length(times_eval))
  zvar <- paste0(var, "_z")
  valid_idx <- !is.na(data[[zvar]]) & !is.na(data$time) & !is.na(data$state)
  
  for (i in seq_along(times_eval)) {
    t <- times_eval[i]
    n_event <- sum(data$time <= t & data$state == 1, na.rm = TRUE)
    
    # Check minimum event count and sample size
    if (n_event >= 5 && sum(valid_idx) >= 30) {
      roc_obj <- tryCatch({
        timeROC(T = data$time[valid_idx],
                delta = data$state[valid_idx],
                marker = data[[zvar]][valid_idx],
                cause = 1,
                times = t,
                iid = TRUE)
      }, error = function(e) NULL)
      if (!is.null(roc_obj)) aucs[i] <- roc_obj$AUC[1]
    }
  }
  
  # --- 3. Bootstrap CI for Cox Coefficient ---
  boot_stat <- function(d, inds) {
    dd <- d[inds, , drop = FALSE]
    fitb <- tryCatch(coxph(as.formula(paste("Surv(time, state) ~", var)), data = dd), error = function(e) NA)
    if (is.na(fitb[1])) return(NA)
    return(coef(fitb))
  }
  
  boot_res <- tryCatch(boot(data = data, statistic = boot_stat, R = R_boot), error = function(e) NULL)
  boot_ci_coef <- c(NA, NA)
  
  if (!is.null(boot_res) && !all(is.na(boot_res$t))) {
    tvals <- boot_res$t; tvals <- tvals[!is.na(tvals)]
    # Use percentile method for CI
    if (length(tvals) >= 20) boot_ci_coef <- quantile(tvals, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  
  boot_hr <- exp(boot_ci_coef); boot_hr_0.1 <- exp(boot_ci_coef * 0.1)
  
  return(list(var = var, coef = coefv, SE = sev, z = z, p = p,
              HR_1 = HR_1, HR_1_ci = HR_1_ci,
              HR_0.1 = HR_0.1, HR_0.1_ci = HR_0.1_ci,
              C_index = c_index, AUCs = aucs,
              boot_hr = boot_hr, boot_hr_0.1 = boot_hr_0.1))
}

# Execute Core Analysis for all variables
res_list <- lapply(VARS_TO_ANALYZE, function(v) cox_one(v, data, surv_obj, TIMES_EVAL, R_BOOTSTRAP))


# ---
# Section 4: Export Time-Dependent AUC Table
# ---

timeROC_results <- data.frame(Variable = character(), Time = integer(), Events = integer(), AtRisk = integer(), AUC = numeric(), stringsAsFactors = FALSE)
for (v in VARS_TO_ANALYZE) {
  zvar <- paste0(v, "_z")
  valid_idx <- !is.na(data[[zvar]]) & !is.na(data$time) & !is.na(data$state)
  
  for (t in TIMES_EVAL) {
    n_risk <- sum(data$time >= t & valid_idx, na.rm = TRUE)
    n_event <- sum(data$time <= t & data$state == 1 & valid_idx, na.rm = TRUE)
    auc_val <- NA
    
    if (n_event >= 5) {
      roc_obj <- tryCatch({
        timeROC(T = data$time[valid_idx],
                delta = data$state[valid_idx],
                marker = data[[zvar]][valid_idx],
                cause = 1,
                times = t,
                iid = TRUE)
      }, error = function(e) NULL)
      if (!is.null(roc_obj)) auc_val <- roc_obj$AUC[1]
    }
    timeROC_results <- rbind(timeROC_results, data.frame(Variable = v, Time = t, Events = n_event, AtRisk = n_risk, AUC = auc_val))
  }
}
write.table(timeROC_results, file = TIME_ROC_OUTPUT, sep = "\t", quote = FALSE, row.names = FALSE)


# ---
# Section 5: Generate Report (Sink Output)
# ---

sink(REPORT_OUTPUT)
cat("========= GSE24080 Single-Factor Cox & TimeROC Report (Short Windows) =========\n\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Bootstrap Replicates (R): %d\n", R_BOOTSTRAP))
cat(sprintf("Time Windows (days): %s\n\n", paste(TIMES_EVAL, collapse = ", ")))
cat("Note: Raw Cox uses original variable scale; TimeROC uses Z-score scale for stability.\n\n")

for (r in res_list) {
  vname <- r$var
  cat("---- Variable:", vname, "----\n")
  cat(sprintf("Coef = %.4f, SE = %.4f, Z = %.3f, P = %.4g\n", r$coef, r$SE, r$z, r$p))
  cat(sprintf("HR (per 1 unit) = %.4f (95%% CI: %.4f - %.4f)\n", r$HR_1, r$HR_1_ci[1], r$HR_1_ci[2]))
  cat(sprintf("HR (per 0.1 unit) = %.4f (95%% CI: %.4f - %.4f)\n", r$HR_0.1, r$HR_0.1_ci[1], r$HR_0.1_ci[2]))
  cat(sprintf("C-index = %s\n", ifelse(is.na(r$C_index), "NA", format(round(r$C_index, 4), nsmall = 4))))
  
  # TimeROC AUC
  if (all(is.na(r$AUCs))) {
    cat("Time-dependent AUCs: NA (Insufficent events/samples or calculation failure)\n")
  } else {
    for (i in seq_along(TIMES_EVAL)) {
      cat(sprintf("  Time %d days: AUC = %s\n", TIMES_EVAL[i], ifelse(is.na(r$AUCs[i]), "NA", format(round(r$AUCs[i], 4), nsmall = 4))))
    }
  }
  
  # Bootstrap CI
  if (all(is.na(r$boot_hr))) {
    cat("Bootstrap HR CI: NA (Bootstrap failed or insufficient data)\n")
  } else {
    cat(sprintf("Bootstrap HR CI (per 1 unit): %.4f - %.4f\n", r$boot_hr[1], r$boot_hr[2]))
    cat(sprintf("Bootstrap HR CI (per 0.1 unit): %.4f - %.4f\n", r$boot_hr_0.1[1], r$boot_hr_0.1[2]))
  }
  cat("\n")
}
cat(sprintf("========= timeROC table saved to: %s =========\n", TIME_ROC_OUTPUT))
sink()


# ---
# Section 6: Subgroup Survival Analysis (KM by Clinical Variables)
# ---
cat("\n--- Subgroup Analysis (Reported to Console) ---\n")
candidate_clin <- c("ISS", "iss", "Stage", "stage", "age", "Age", "gender", "Gender")
present_clin <- intersect(candidate_clin, colnames(data))
MAIN_VAR <- "riskex1" # Variable used for Cox in subgroups

if (length(present_clin) == 0) {
  cat("No common clinical variables (ISS/age/gender/stage) found. Skipping subgroup analysis.\n")
} else {
  for (clin in present_clin) {
    cat(sprintf("\nProcessing Subgroup: %s\n", clin))
    colv <- data[[clin]]
    
    # Define groups (median split for numeric, direct factor for categorical)
    if (is.numeric(colv)) {
      grp <- factor(ifelse(colv > median(colv, na.rm = TRUE), "high", "low"), levels = c("low", "high"))
      clin_title <- paste(clin, "(Median Split)")
    } else {
      grp <- as.factor(colv)
      clin_title <- clin
    }
    
    idx <- !is.na(grp) & !is.na(data$time) & !is.na(data$state) & !is.na(data[[MAIN_VAR]])
    
    if (sum(idx) < 30) { cat("  Skipping: Insufficient valid samples.\n"); next }
    
    dat_sub <- data[idx, , drop = FALSE]; grp_f <- grp[idx]
    levs <- levels(grp_f)
    
    for (lv in levs) {
      subd <- dat_sub[grp_f == lv, , drop = FALSE]
      if (nrow(subd) < 10) { next }
      
      # Subgroup Cox analysis
      fit <- coxph(Surv(time, state) ~ riskex1, data = subd)
      s <- summary(fit)
      coefv <- s$coefficients[1, "coef"]; pval <- s$coefficients[1, "Pr(>|z|)"]
      HR <- exp(coefv); CI <- exp(coefv + c(-1, 1) * 1.96 * s$coefficients[1, "se(coef)"])
      
      cat(sprintf("  Subgroup %s=%s: n=%d, coef=%.4f, HR=%.4f (95%%CI %.4f-%.4f), p=%.4g\n",
                  clin, lv, nrow(subd), coefv, HR, CI[1], CI[2], pval))
    }
    
    # Save Subgroup KM Plot
    try({
      fit_km <- survfit(Surv(time, state) ~ grp_f, data = dat_sub)
      png(filename = paste0(FIG_PREFIX, "KM_by_", clin, ".png"), width = 800, height = 600)
      plot(fit_km, xlab = "Time (days)", ylab = "Survival", main = paste("KM by", clin_title), lwd = 2)
      legend("topright", legend = levels(grp_f), lty = 1:length(levels(grp_f)), col = 1:length(levels(grp_f)), bty = "n")
      dev.off()
    }, silent = TRUE)
  }
}

# ---
# Section 7: KM (Median Split) and Cumulative Event Plot for Main Variable (riskex1)
# ---

idx_all <- !is.na(data[[MAIN_VAR]]) & !is.na(data$time) & !is.na(data$state)
dat0 <- data[idx_all, , drop = FALSE]

# Create risk group based on median split of MAIN_VAR
median_val <- median(dat0[[MAIN_VAR]], na.rm = TRUE)
dat0$risk_group <- factor(ifelse(dat0[[MAIN_VAR]] > median_val, "high", "low"), levels = c("low", "high"))
fit_km_rg <- survfit(Surv(time, state) ~ risk_group, data = dat0)

# 1. KM Plot (Median Split)
png(filename = paste0(FIG_PREFIX, "KM_", MAIN_VAR, "_median.png"), width = 700, height = 500)
plot(fit_km_rg, xlab = "Time (days)", ylab = "Survival", main = paste("KM by", MAIN_VAR, "(median split)"), lwd = 2, col = c("blue", "red"))
legend("topright", legend = levels(dat0$risk_group), col = c("blue", "red"), lwd = 2, bty = "n")
dev.off()
cat(sprintf("\nKM Plot (median split) saved: %s\n", paste0(FIG_PREFIX, "KM_", MAIN_VAR, "_median.png")))

# Log-rank test and Cox summary
lr <- survdiff(Surv(time, state) ~ risk_group, data = dat0)
p_lr <- 1 - pchisq(lr$chisq, df = 1)
fit_cox_rg <- coxph(Surv(time, state) ~ risk_group, data = dat0)
srg <- summary(fit_cox_rg)
coefv <- srg$coefficients[1, "coef"]; HR <- exp(coefv); CI <- exp(coefv + c(-1, 1) * 1.96 * srg$coefficients[1, "se(coef)"])
cat(sprintf("KM (median split) Summary: n_low=%d, n_high=%d, Log-rank p=%.4g\n", sum(dat0$risk_group == "low"), sum(dat0$risk_group == "high"), p_lr))
cat(sprintf("Cox (risk_group high vs low): HR=%.4f (95%%CI %.4f-%.4f), p=%.4g\n", HR, CI[1], CI[2], srg$coefficients[1, "Pr(>|z|)"]))


# 2. Cumulative Event Plot (1 - Survival)
png(filename = paste0(FIG_PREFIX, "CumulativeEvents_", MAIN_VAR, "_median.png"), width = 700, height = 500)
sums <- summary(fit_km_rg)
strata_sizes <- as.numeric(table(sums$strata))
start_idx <- 1; cols <- c("blue", "red")
plot(NA, xlim = range(sums$time, na.rm = TRUE), ylim = c(0, 1), xlab = "Time (days)", ylab = "Cumulative event rate",
     main = "Cumulative Events by Risk Group (Median Split)", type = "n")

for (i in seq_along(strata_sizes)) {
  npoints <- strata_sizes[i]
  times_i <- sums$time[start_idx:(start_idx + npoints - 1)]
  surv_i <- sums$surv[start_idx:(start_idx + npoints - 1)]
  lines(times_i, 1 - surv_i, col = cols[i], type = "s", lwd = 2)
  start_idx <- start_idx + npoints
}
legend("topleft", legend = levels(dat0$risk_group), col = cols, lwd = 2, bty = "n")
dev.off()
cat(sprintf("Cumulative Events Plot saved: %s\n", paste0(FIG_PREFIX, "CumulativeEvents_", MAIN_VAR, "_median.png")))


# ---
# Section 8: Proportional Hazards Check and Time-Dependent Modeling
# ---

fit_overall <- coxph(Surv(time, state) ~ riskex1, data = data)
cat("\n--- Proportional Hazards Check (Variable: riskex1) ---\n")
zph <- cox.zph(fit_overall)
cat("cox.zph Test Results:\n"); print(zph); cat("\n")

if (any(zph$table[, 3, drop = FALSE] < 0.05, na.rm = TRUE)) {
  cat("Non-proportional hazards detected (p<0.05). Fitting time-dependent model (tt()).\n")
  
  fit_tt <- tryCatch({
    coxph(Surv(time, state) ~ riskex1 + tt(riskex1), data = data, tt = function(x, t, ...) x * log(t + 1))
  }, error = function(e) NULL)
  
  if (!is.null(fit_tt)) {
    cat("Time-Dependent Cox (tt) Model Summary:\n"); print(summary(fit_tt)); cat("\n")
    coefs <- coef(fit_tt)
    
    if (length(coefs) >= 2) {
      beta0 <- coefs[1]; beta_t <- coefs[2]
      tt_seq <- seq(30, max(data$time, na.rm = TRUE), length.out = 200)
      HR_t <- exp(beta0 + beta_t * log(tt_seq + 1))
      
      # Plot HR(t)
      png(filename = paste0(FIG_PREFIX, "HR_over_time_tt.png"), width = 800, height = 500)
      plot(tt_seq, HR_t, type = "l", lwd = 2, xlab = "Time (days)", ylab = "HR(t)", main = "Riskex1 HR over time (tt model)")
      abline(h = 1, lty = 2)
      dev.off()
      cat(sprintf("HR(t) Plot (tt model) saved: %s\n", paste0(FIG_PREFIX, "HR_over_time_tt.png")))
    }
  }
} else {
  cat("cox.zph did not detect significant non-proportional hazards (p>=0.05). Proportional hazards assumption holds.\n")
}

# survSplit Piecewise Model (using z-score for stability)
try({
  cuts <- c(90, 180, 365)
  if (max(data$time, na.rm = TRUE) > min(cuts)) {
    ds <- survSplit(Surv(time, state) ~ ., data = data, cut = cuts, episode = "period", id = "idtemp")
    ds$logt <- log(ds$time + 1)
    
    # Use z-score version of the variable
    zvar_split <- paste0(MAIN_VAR, "_z")
    fit_split <- coxph(as.formula(paste("Surv(tstart, time, state) ~", zvar_split, "+ I(", zvar_split, "* logt)")), data = ds)
    
    cat("SurvSplit Piecewise Model Summary (using z-score):\n"); print(summary(fit_split)); cat("\n")
    coef_s <- coef(fit_split)
    
    if (!is.null(coef_s) && any(grepl(zvar_split, names(coef_s)))) {
      b0_name <- names(coef_s)[grep(paste0("^", zvar_split, "$"), names(coef_s))]
      bt_name <- names(coef_s)[grep(paste0("I\\(", zvar_split), names(coef_s))]
      
      b0 <- coef_s[b0_name]; bt <- coef_s[bt_name]
      ttseq <- seq(30, max(data$time, na.rm = TRUE), length.out = 200)
      HR_ts <- exp(b0 + bt * log(ttseq + 1))
      
      # Plot HR(t)
      png(filename = paste0(FIG_PREFIX, "HR_over_time_split.png"), width = 800, height = 500)
      plot(ttseq, HR_ts, type = "l", lwd = 2, xlab = "Time (days)", ylab = "HR(t)", main = paste(MAIN_VAR, "HR over time (survSplit model)"))
      abline(h = 1, lty = 2)
      dev.off()
      cat(sprintf("HR(t) Piecewise Plot saved: %s\n", paste0(FIG_PREFIX, "HR_over_time_split.png")))
    }
  }
}, silent = TRUE)

cat("\nAnalysis complete. Detailed report saved to:", REPORT_OUTPUT, "\n")