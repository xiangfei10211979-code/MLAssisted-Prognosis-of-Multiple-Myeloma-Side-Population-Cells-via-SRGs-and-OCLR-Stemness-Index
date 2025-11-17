# ---README ---

Project purpose: reproducible pipeline for mRNAsi, WGCNA, ssGSEA, pathway-based prognostic optimization and validation using R and MATLAB. Follow the "Gemini" standards below before running analyses.

Top-level instructions (must be applied across all scripts):
- Remove absolute local paths; use relative paths only. Define every path in a single Configuration block at the top of each script and never hard-code elsewhere.
- Remove all non-English comments and any obsolete/commented-out code blocks. Use a one-line header comment at each script top: purpose, inputs, outputs, author (optional).
- Load all libraries at script top. List external data/file dependencies clearly in the Configuration.
- Save all generated intermediate and final outputs under `results/` subfolders. Do not write outputs to parent or system folders.

Recommended Configuration snippet (copy into the top of every R / MATLAB driver script and edit once per project clone):
CONFIGURATION (place exactly once at script top)
ROOT_DIR = "../"                       # project root relative to script
DATA_DIR = file.path(ROOT_DIR, "data")
INPUT_DIR = file.path(DATA_DIR, "input")
PROCESSED_DIR = file.path(DATA_DIR, "processed")
REF_DIR = file.path(DATA_DIR, "ref")
MAIN_DIR = file.path(ROOT_DIR, "main")
RESULTS_DIR = file.path(ROOT_DIR, "results")
FIGURES_DIR = file.path(RESULTS_DIR, "figures")
TABLES_DIR = file.path(RESULTS_DIR, "tables")
# Example: in R use file.path(CONFIG$INPUT_DIR, "TCGA_MMRF_expression.txt")

Directory layout (project root = C:/MMcode when local; use relative paths in scripts):
README.md
data/
  input/
    Stemness_index.rda
    TCGA_MMRF_expression.txt
    GSE109651_SP_cells.txt
    pathway_gene_lists.gmt
    GSE19784_expression.txt
    GSE24080_expression.txt
    GSE24080_prognostic_data.txt
    GSE57317_expression.txt
    KEGG_gene_sets.txt
    KEGG_gene_sets2.txt
    module_blue.txt
    module_green.txt
    risk_group_data.txt
    risk_group_tide.txt
    score.csv
    SP_expression.txt
    survival_data_risk_group.txt
    TCGA_expression.txt
    tide_scores_raw.txt
    wgcna_clinical_traits.txt
    wgcna_expression_input.txt
  processed/
    TCGA_MMRF_combat.txt
    hubGenesblue.txt
    module_all.txt
    module_blue.txt
    preservation_blue.csv
  ref/
    GDSC2_Expr.rds
    GDSC2_Res.rds
main/
  01_calculate_mRNAsi.R
  02_run_WGCNA.R
  03_run_ssGSEA.R
  04_survival_analysis.R
  05_drug_sensitivity_oncopredict.R
  06_immunotherapy_tide.R
  07_prognostic_validation.R
  08_01wgcna_module_preservationM1.R
  08_02wgcna_module_preservationM2.R
  08_03wgcna_module_preservationM3.R
  03_run_OptW.m
  Run_OptW.m
  Cox_Opt.m
  HRConfidenceConstraint.m
results/
  figures/
    KM_plot_GSE24080.pdf
  tables/
    mRNAsi_scores.txt
    OptW_weights.csv

Mapping of existing absolute paths (how to convert to relative):
- "C:\MMcode\data\input\TCGA_MMRF_expression.txt" -> file.path(INPUT_DIR, "TCGA_MMRF_expression.txt")
- "C:\MMcode\main\03_run_ssGSEA.R" -> file.path(MAIN_DIR, "03_run_ssGSEA.R")
Apply the same mapping for every absolute path you shared.

Execution order (minimal reproducible run):
1) Place raw inputs into data/input/ and reference files into data/ref/.
2) From project root run stepwise:
   - R: Rscript --vanilla main/01_calculate_mRNAsi.R
   - R: Rscript --vanilla main/02_run_WGCNA.R
   - R: Rscript --vanilla main/03_run_ssGSEA.R
   - MATLAB (local session or batch): open MATLAB, setpwd to project root and run "Run_OptW.m" or call: matlab -batch "cd('main'); Run_OptW"
   - R: Rscript --vanilla main/04_survival_analysis.R
   - R: Rscript --vanilla main/05_drug_sensitivity_oncopredict.R
   - R: Rscript --vanilla main/06_immunotherapy_tide.R
   - R: Rscript --vanilla main/07_prognostic_validation.R
3) Verify outputs saved under results/figures and results/tables after each step.

Dependencies and environment:
- R (recommended >= 4.3.0; use actual version in manuscript) and MATLAB R2023b (Optimization Toolbox & Statistics and Machine Learning Toolbox required).
- Key R packages to confirm/ install at script start: WGCNA, survival, survminer, timeROC, boot, GSVA, GSEABase, oncoPredict, limma, ggpubr, MASS, muhaz, pheatmap, reshape2. Use a small helper in R to check and install missing packages:
  if (!requireNamespace("pkgname", quietly=TRUE)) install.packages("pkgname") or BiocManager::install("pkgname")
- Put package checks in a single "00_setup_environment.R" if desired. Save packageVersion() outputs to results/tables/R_packages_versions.txt for reproducibility.

Coding standards checklist (must be enforced before commit):
- All scripts contain a one-paragraph header (purpose, inputs, outputs).
- All file reads/writes use relative paths from CONFIG.
- No Chinese text or inline absolute paths.
- All library() calls at the top.
- No long commented blocks or dead code. Keep functions small and documented with short comments.
- Intermediate files written only to data/processed/ and results/ subfolders.

Reproducibility tips:
- Commit scripts and a small metadata file listing R and MATLAB versions and package versions.
- Save a copy of sessionInfo() output into results/tables/sessionInfo_R.txt.
- Save MATLAB version using ver command into results/tables/sessionInfo_MATLAB.txt.
- If using Git, add a .gitignore that excludes large raw datasets but includes small reference files and scripts.

Citations and methods text:
- Keep the short Methods paragraph (package names + citation keys) in your manuscript Methods and include the bib entries in the repo bib file. Do not generate citations in scripts; place citation keys in manuscript only.

Quick troubleshooting:
- If a table compile error occurs in LaTeX when using \resizebox, ensure \usepackage{graphicx} is loaded and that \caption is outside the resizebox block.
- If any script fails due to missing file, confirm path = file.path(CONFIG$INPUT_DIR, "<filename>") and that file exists.

End of README. Follow the checklist above and run scripts in the order listed. If any step fails, report the exact error and the script name for targeted fixes.
