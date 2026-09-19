#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages(
  library(BayesPrism)
)

ROOT <- "/home/mazekai/IBD_EcoTyper"

S1_RDS <- file.path(
  ROOT,
  "03_reference/combined_reference/output/14_bayesprism_broad_run",
  "bayesprism_broad14_result.rds"
)

S2_RDS <- file.path(
  ROOT,
  "03_reference/combined_reference/output/20_bayesprism_epithelial_full_run",
  "bayesprism_epithelial7_result.rds"
)

OUT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/22_hierarchical_fraction_merge"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("Loading BayesPrism results...\n")

bp1 <- readRDS(S1_RDS)
bp2 <- readRDS(S2_RDS)

s1_final <- bp1@posterior.theta_f@theta
s1_initial <- bp1@posterior.initial.cellType@theta

s2_final <- bp2@posterior.theta_f@theta
s2_initial <- bp2@posterior.initial.cellType@theta

cat(
  "Stage1 final:",
  nrow(s1_final), "x", ncol(s1_final), "\n"
)

cat(
  "Stage2 final:",
  nrow(s2_final), "x", ncol(s2_final), "\n"
)

stopifnot(
  ncol(s1_final) == 14,
  "Epithelial" %in% colnames(s1_final)
)

epi7 <- c(
  "Colonic_absorptive",
  "Enteroendocrine",
  "Goblet",
  "Ileal_absorptive",
  "Paneth",
  "Stem_TA_progenitor",
  "Tuft"
)

stopifnot(
  setequal(colnames(s2_final), epi7)
)

s2_final <- s2_final[, epi7, drop = FALSE]
s2_initial <- s2_initial[, epi7, drop = FALSE]

# ------------------------------------------------------------------
# Sample alignment
# ------------------------------------------------------------------

stopifnot(
  setequal(rownames(s1_final), rownames(s2_final))
)

s2_final <- s2_final[
  rownames(s1_final),
  ,
  drop = FALSE
]

s2_initial <- s2_initial[
  rownames(s1_final),
  ,
  drop = FALSE
]

cat(
  "Sample alignment:",
  identical(rownames(s1_final), rownames(s2_final)),
  "\n"
)

# ------------------------------------------------------------------
# Closure before merge
# ------------------------------------------------------------------

s1_err <- max(abs(rowSums(s1_final) - 1))
s2_err <- max(abs(rowSums(s2_final) - 1))

cat(
  "Stage1 final max closure error:",
  format(s1_err, scientific = TRUE),
  "\n"
)

cat(
  "Stage2 final max closure error:",
  format(s2_err, scientific = TRUE),
  "\n"
)

stopifnot(
  s1_err < 1e-8,
  s2_err < 1e-8
)

# ------------------------------------------------------------------
# Hierarchical reconstruction
# ------------------------------------------------------------------

stage1_epi <- s1_final[, "Epithelial"]

epi_whole <- sweep(
  s2_final,
  MARGIN = 1,
  STATS = stage1_epi,
  FUN = "*"
)

epi_parent_err <- max(
  abs(
    rowSums(epi_whole) -
      stage1_epi
  )
)

cat(
  "Reconstructed epithelial parent-mass error:",
  format(epi_parent_err, scientific = TRUE),
  "\n"
)

non_epi13 <- setdiff(
  colnames(s1_final),
  "Epithelial"
)

full20 <- cbind(
  s1_final[, non_epi13, drop = FALSE],
  epi_whole
)

full20_err <- max(
  abs(
    rowSums(full20) - 1
  )
)

cat(
  "Full20 max closure error:",
  format(full20_err, scientific = TRUE),
  "\n"
)

stopifnot(
  epi_parent_err < 1e-8,
  full20_err < 1e-8
)

# ------------------------------------------------------------------
# Primary panel
# ------------------------------------------------------------------

primary_epi4 <- c(
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine"
)

excluded_epi3 <- c(
  "Paneth",
  "Stem_TA_progenitor",
  "Tuft"
)

primary17 <- cbind(
  s1_final[, non_epi13, drop = FALSE],
  epi_whole[, primary_epi4, drop = FALSE]
)

excluded3 <- epi_whole[
  ,
  excluded_epi3,
  drop = FALSE
]

mass_identity_err <- max(
  abs(
    rowSums(primary17) +
      rowSums(excluded3) -
      1
  )
)

stopifnot(
  mass_identity_err < 1e-8
)

# ------------------------------------------------------------------
# Initial vs final sensitivity only
# ------------------------------------------------------------------

rho_epi <- cor(
  s1_initial[, "Epithelial"],
  s1_final[, "Epithelial"],
  method = "spearman"
)

mae_epi <- mean(
  abs(
    s1_initial[, "Epithelial"] -
      s1_final[, "Epithelial"]
  )
)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------

saveRDS(
  full20,
  file.path(
    OUT,
    "whole_biopsy_fraction_full20.rds"
  )
)

saveRDS(
  primary17,
  file.path(
    OUT,
    "whole_biopsy_fraction_primary17.rds"
  )
)

saveRDS(
  excluded3,
  file.path(
    OUT,
    "whole_biopsy_fraction_excluded3.rds"
  )
)

write.table(
  data.frame(
    sample_id = rownames(full20),
    full20,
    check.names = FALSE
  ),
  gzfile(
    file.path(
      OUT,
      "whole_biopsy_fraction_full20.tsv.gz"
    ),
    "wt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sample_qc <- data.frame(
  sample_id = rownames(full20),
  stage1_epithelial_final = stage1_epi,
  reconstructed_epi_sum = rowSums(epi_whole),
  full20_sum = rowSums(full20),
  primary17_sum = rowSums(primary17),
  excluded3_sum = rowSums(excluded3),
  stringsAsFactors = FALSE
)

write.table(
  sample_qc,
  file.path(
    OUT,
    "fraction_merge_sample_qc.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_qc <- data.frame(
  metric = c(
    "n_samples",
    "stage1_closure_max_error",
    "stage2_closure_max_error",
    "epithelial_parent_max_error",
    "full20_closure_max_error",
    "primary17_plus_excluded3_max_error",
    "stage1_initial_final_epi_spearman",
    "stage1_initial_final_epi_MAE"
  ),
  value = c(
    nrow(full20),
    s1_err,
    s2_err,
    epi_parent_err,
    full20_err,
    mass_identity_err,
    rho_epi,
    mae_epi
  )
)

write.table(
  summary_qc,
  file.path(
    OUT,
    "fraction_merge_qc_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================\n")
cat("FRACTION MERGE QC\n")
cat("============================================\n")
print(summary_qc, row.names = FALSE)

cat("\nPrimary17 row-sum summary:\n")
print(summary(rowSums(primary17)))

cat("\nExcluded3 mass summary:\n")
print(summary(rowSums(excluded3)))

cat("\nFINAL STATUS: PASS\n")
cat("Output:", OUT, "\n")

