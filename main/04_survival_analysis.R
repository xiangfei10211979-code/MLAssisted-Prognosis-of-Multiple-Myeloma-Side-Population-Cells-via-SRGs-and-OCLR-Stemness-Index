#
# Script: 04_survival_analysis.R
#
# Description:
# Performs comprehensive prognostic analysis including Kaplan-Meier (KM)
# survival curves, Cox regression (Forest Plot), and non-parametric
# hazard function estimation for the derived risk stratification model.
#
# Key Steps:
# 1. Load survival data (time, state, risk group).
# 2. Run KM analysis and plot the survival curve with the risk table.
# 3. Run Univariate Cox regression and generate a Forest Plot.
# 4. Estimate and plot the non-parametric hazard function.
#
# Project:      MM SRGs Prognostic Model
#
# ---

# ---
# Section 1: Load Libraries and Configuration
# ---
library(survival)
library(survminer)
library(ggplot2)
library(muhaz) # Required for the Hazard Function Plot

# Input Files
INPUT_DATA_FILE <- "../data/input/survival_data_risk_group.txt" # Requires time, state, risk columns

# Output Directories
FIG_DIR <- "../results/figures/survival/"

# Create directory if it doesn't exist
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Define plotting constants
RISK_LEVELS <- c("Low", "High") # Expected levels for the 'risk' factor
PLOT_PALETTE <- c("Low" = "blue", "High" = "red")

# ---
# Section 2: Load and Prepare Data
# ---

data <- read.table(INPUT_DATA_FILE, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Ensure 'risk' column is a factor with defined levels
data$risk <- factor(data$risk, levels = sort(unique(data$risk)))

# Map 'risk' levels to descriptive names for plotting (e.g., if levels are 0, 1)
# NOTE: If your input 'risk' column contains "Low" and "High" text directly,
# you can skip this remapping and only update RISK_LEVELS/PLOT_PALETTE above.
if (length(levels(data$risk)) == 2) {
  levels(data$risk) <- RISK_LEVELS
} else {
  stop("Error: 'risk' factor does not have exactly two levels.")
}

# Create survival object: state=0 (censored), state=1 (event/death)
surv_object <- Surv(time = data$time, event = ifelse(data$state == 0, 1, 0))

# ---
# Section 3: Kaplan-Meier Survival Analysis
# ---

fit <- survfit(surv_object ~ risk, data = data)

# Plot KM Curve with Risk Table
survival_plot_full <- ggsurvplot(fit,
                                 data = data,
                                 surv.median.line = "hv",
                                 pval = TRUE,
                                 conf.int = TRUE,
                                 risk.table = TRUE,
                                 tables.height = 0.2,
                                 risk.table.y.text = TRUE,
                                 risk.table.col = "strata",
                                 risk.table.fontsize = 4,
                                 xlab = "Time (days)",
                                 ylab = "Survival Probability",
                                 title = "Kaplan-Meier Survival Curve by Risk Group",
                                 legend.labs = RISK_LEVELS,
                                 legend.title = "Risk Group",
                                 palette = PLOT_PALETTE,
                                 ggtheme = theme_bw(),
                                 break.time.by = max(data$time, na.rm = TRUE) / 5
)

# Save the full plot (curve + risk table)
ggsave(file.path(FIG_DIR, "01_KM_curve_full.pdf"),
       plot = survival_plot_full$plot,
       width = 8, height = 6)

# Save the risk table only (for separate insertion in manuscripts)
ggsave(file.path(FIG_DIR, "02_risk_table_only.pdf"),
       plot = survival_plot_full$table,
       width = 8, height = 3)

cat(paste("Kaplan-Meier plots saved to:", FIG_DIR, "\n"))

# ---
# Section 4: Univariate Cox Regression (Forest Plot)
# ---

cox_model <- coxph(surv_object ~ risk, data = data)
cox_summary <- summary(cox_model)

# Extract statistics
hr <- cox_summary$conf.int[, "exp(coef)"]
hr_lower <- cox_summary$conf.int[, "lower .95"]
hr_upper <- cox_summary$conf.int[, "upper .95"]
p_value <- cox_summary$coefficients[, "Pr(>|z|)"]

forest_data <- data.frame(
  Feature = "Risk (High vs Low)",
  HR = hr,
  lower = hr_lower,
  upper = hr_upper,
  p.value = p_value
)

# Generate Forest Plot
forest_plot <- ggplot(forest_data, aes(y = Feature, x = HR, xmin = lower, xmax = upper)) +
  geom_errorbarh(height = 0.2, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_point(size = 3, color = PLOT_PALETTE["High"]) +
  labs(title = "Forest Plot for Risk Group (Univariate Cox Regression)",
       y = "",
       x = "Hazard Ratio (95% CI)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(FIG_DIR, "03_forest_plot_risk.pdf"), plot = forest_plot, width = 6, height = 2)
cat(paste("Forest plot saved to:", FIG_DIR, "\n"))

# ---
# Section 5: Hazard Function Plot (muhaz)
# ---
# Filter data by risk group
data_low_risk <- subset(data, risk == RISK_LEVELS[1])
data_high_risk <- subset(data, risk == RISK_LEVELS[2])

# Estimate hazard function for each group (using product-limit estimator)
hazard_estimate_low <- est.haz(Surv(time, ifelse(state == 0, 1, 0)), data = data_low_risk)
hazard_estimate_high <- est.haz(Surv(time, ifelse(state == 0, 1, 0)), data = data_high_risk)

# Convert estimates to data frames for ggplot
hazard_df_low <- data.frame(time = hazard_estimate_low$time, hazard = hazard_estimate_low$hazard, Group = RISK_LEVELS[1])
hazard_df_high <- data.frame(time = hazard_estimate_high$time, hazard = hazard_estimate_high$hazard, Group = RISK_LEVELS[2])
hazard_df <- rbind(hazard_df_low, hazard_df_high)

# Generate Hazard Plot
hazard_plot <- ggplot(hazard_df, aes(x = time, y = hazard, color = Group)) +
  geom_line(linewidth = 1) +
  labs(title = "Estimated Hazard Function by Risk Group",
       x = "Time (days)",
       y = "Hazard Rate") +
  scale_color_manual(values = PLOT_PALETTE) +
  theme_bw()

ggsave(file.path(FIG_DIR, "04_hazard_function_plot.pdf"), plot = hazard_plot, width = 8, height = 6)
cat(paste("Hazard function plot saved to:", FIG_DIR, "\n"))