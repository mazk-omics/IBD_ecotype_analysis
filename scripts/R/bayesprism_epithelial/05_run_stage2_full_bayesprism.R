#!/usr/bin/env Rscript

# ==============================================================================
# IBD EcoTyper
# Stage 2 epithelial-only BayesPrism
#
# S2-C: FULL FORMAL RUN
#
# Runs in one call:
#   1. Initial Gibbs sampling
#   2. Reference update
#   3. Final Gibbs sampling
#
# Input:
#   Formal Stage 2 epithelial prism from S2-B2b
#
# Output:
#   - complete BayesPrism result
#   - initial epithelial fractions
#   - final epithelial fractions
#   - 7 epithelial subtype Z expression matrices
#   - updated reference
#   - read-mass conservation audit
#
# Stage 2 design:
#   2490 samples
#   13668 genes
#   7 epithelial components
#   cell.type = cell.state
# ==============================================================================


# ==============================================================================
# 0. Settings
# ==============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(BayesPrism)
})


# ==============================================================================
# 1. Project paths
# ==============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

PRISM_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "19_bayesprism_epithelial_prism"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "20_bayesprism_epithelial_full_run"
)

EXPRESSION_DIR <- file.path(
  OUTPUT_DIR,
  "celltype_expression"
)

SCRIPT_DIR <- file.path(
  PROJECT_ROOT,
  "scripts",
  "R",
  "bayesprism_epithelial"
)

dir.create(
  SCRIPT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  EXPRESSION_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Input
# ==============================================================================

PRISM_RDS <- file.path(
  PRISM_DIR,
  "bayesprism_epithelial7_prism_object.rds"
)


# ==============================================================================
# 3. Outputs
# ==============================================================================

OUT_RESULT <- file.path(
  OUTPUT_DIR,
  "bayesprism_epithelial7_result.rds"
)

OUT_REFERENCE_UPDATE <- file.path(
  OUTPUT_DIR,
  "bayesprism_epithelial7_reference_update.rds"
)

OUT_INITIAL_FRACTION <- file.path(
  OUTPUT_DIR,
  "bayesprism_epithelial7_fraction_initial.tsv"
)

OUT_FINAL_FRACTION <- file.path(
  OUTPUT_DIR,
  "bayesprism_epithelial7_fraction_final.tsv"
)

OUT_CONSERVATION_SAMPLE <- file.path(
  OUTPUT_DIR,
  "stage2_read_mass_conservation_by_sample.tsv"
)

OUT_COMPONENT_SUMMARY <- file.path(
  OUTPUT_DIR,
  "stage2_epithelial_component_summary.tsv"
)

OUT_AUDIT <- file.path(
  OUTPUT_DIR,
  "stage2_full_run_audit.tsv"
)

OUT_CONTROL <- file.path(
  OUTPUT_DIR,
  "stage2_full_run_control.tsv"
)

OUT_MANIFEST <- file.path(
  OUTPUT_DIR,
  "stage2_full_run_manifest.tsv"
)

OUT_LOG <- file.path(
  OUTPUT_DIR,
  "stage2_full_run.log"
)

OUT_SESSIONINFO <- file.path(
  OUTPUT_DIR,
  "sessionInfo.txt"
)


# ==============================================================================
# 4. Logging
# ==============================================================================

log_con <- file(
  OUT_LOG,
  open = "wt"
)

sink(
  log_con,
  type = "output",
  split = TRUE
)

sink(
  log_con,
  type = "message"
)

on.exit({

  while (sink.number(type = "message") > 0) {
    sink(type = "message")
  }

  while (sink.number(type = "output") > 0) {
    sink(type = "output")
  }

  close(log_con)

}, add = TRUE)


cat("============================================================\n")
cat("Stage 2 epithelial BayesPrism\n")
cat("S2-C: FULL formal BayesPrism run\n")
cat("============================================================\n\n")

cat("START:", format(Sys.time()), "\n")

cat(
  "BayesPrism version:",
  as.character(packageVersion("BayesPrism")),
  "\n\n"
)


# ==============================================================================
# 5. Helpers
# ==============================================================================

assert_true <- function(condition, message) {

  if (!isTRUE(condition)) {
    stop(
      message,
      call. = FALSE
    )
  }
}


write_tsv <- function(x, path) {

  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


audit_row <- function(
  item,
  observed,
  expected = NA_character_,
  status
) {

  data.frame(
    item = as.character(item),
    observed = as.character(observed),
    expected = as.character(expected),
    status = as.character(status),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 6. Load formal Stage 2 prism
# ==============================================================================

assert_true(
  file.exists(PRISM_RDS),
  paste(
    "Missing formal Stage 2 prism:",
    PRISM_RDS
  )
)

cat("Loading formal Stage 2 prism...\n")

prism <- readRDS(
  PRISM_RDS
)

assert_true(
  inherits(prism, "prism"),
  "Input RDS is not a BayesPrism prism object."
)


# ==============================================================================
# 7. Pre-run validation
# ==============================================================================

assert_true(
  nrow(prism@mixture) == 2490L,
  paste0(
    "Unexpected sample count: ",
    nrow(prism@mixture)
  )
)

assert_true(
  ncol(prism@mixture) == 13668L,
  paste0(
    "Unexpected gene count: ",
    ncol(prism@mixture)
  )
)

assert_true(
  nrow(prism@phi_cellType@phi) == 7L,
  paste0(
    "Unexpected cell-type count: ",
    nrow(prism@phi_cellType@phi)
  )
)

assert_true(
  nrow(prism@phi_cellState@phi) == 7L,
  paste0(
    "Unexpected cell-state count: ",
    nrow(prism@phi_cellState@phi)
  )
)

assert_true(
  identical(
    colnames(prism@mixture),
    colnames(prism@phi_cellType@phi)
  ),
  "Mixture and cell-type reference gene order differ."
)

assert_true(
  identical(
    colnames(prism@mixture),
    colnames(prism@phi_cellState@phi)
  ),
  "Mixture and cell-state reference gene order differ."
)


EXPECTED_COMPONENTS <- c(
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine",
  "Tuft",
  "Paneth",
  "Stem_TA_progenitor"
)

assert_true(
  setequal(
    rownames(prism@phi_cellType@phi),
    EXPECTED_COMPONENTS
  ),
  "Unexpected epithelial component identities."
)

cat(
  "[PASS] Prism validated:",
  nrow(prism@mixture),
  "samples x",
  ncol(prism@mixture),
  "genes;",
  nrow(prism@phi_cellType@phi),
  "epithelial components\n\n"
)


# ==============================================================================
# 8. Formal computational parameters
# ==============================================================================

N_CORES <- 25L


GIBBS_CONTROL <- list(

  chain.length = 1000,

  burn.in = 500,

  thinning = 2,

  n.cores = N_CORES,

  seed = 123,

  alpha = 1
)


OPT_CONTROL <- list(

  maxit = 100000,

  maximize = FALSE,

  trace = 0,

  eps = 1e-7,

  dowarn = TRUE,

  tol = 0,

  maxNA = 500,

  n.cores = N_CORES,

  optimizer = "MAP",

  sigma = 2
)


cat("Formal run configuration:\n")

cat("  n.cores       :", N_CORES, "\n")
cat("  chain.length  :", GIBBS_CONTROL$chain.length, "\n")
cat("  burn.in       :", GIBBS_CONTROL$burn.in, "\n")
cat("  thinning      :", GIBBS_CONTROL$thinning, "\n")
cat("  seed          :", GIBBS_CONTROL$seed, "\n")
cat("  alpha         :", GIBBS_CONTROL$alpha, "\n")
cat("  optimizer     :", OPT_CONTROL$optimizer, "\n")
cat("  sigma         :", OPT_CONTROL$sigma, "\n")
cat("  update.gibbs  : TRUE\n\n")


# ==============================================================================
# 9. Save computational provenance before run
# ==============================================================================

control_table <- data.frame(

  parameter = c(
    "BayesPrism_version",
    "n.cores",
    "chain.length",
    "burn.in",
    "thinning",
    "seed",
    "alpha",
    "update.gibbs",
    "optimizer",
    "sigma",
    "maxit",
    "eps",
    "samples",
    "genes",
    "components"
  ),

  value = c(
    as.character(packageVersion("BayesPrism")),
    N_CORES,
    GIBBS_CONTROL$chain.length,
    GIBBS_CONTROL$burn.in,
    GIBBS_CONTROL$thinning,
    GIBBS_CONTROL$seed,
    GIBBS_CONTROL$alpha,
    TRUE,
    OPT_CONTROL$optimizer,
    OPT_CONTROL$sigma,
    OPT_CONTROL$maxit,
    OPT_CONTROL$eps,
    nrow(prism@mixture),
    ncol(prism@mixture),
    nrow(prism@phi_cellType@phi)
  ),

  stringsAsFactors = FALSE
)

write_tsv(
  control_table,
  OUT_CONTROL
)


# ==============================================================================
# 10. FULL BayesPrism run
# ==============================================================================

cat("============================================================\n")
cat("Starting FULL Stage 2 BayesPrism run\n")
cat("Initial Gibbs -> Reference update -> Final Gibbs\n")
cat("============================================================\n\n")

start_time <- Sys.time()


bp <- run.prism(

  prism = prism,

  n.cores = N_CORES,

  update.gibbs = TRUE,

  gibbs.control = GIBBS_CONTROL,

  opt.control = OPT_CONTROL
)


end_time <- Sys.time()


elapsed_seconds <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "secs"
  )
)


cat("\n")
cat("============================================================\n")
cat("FULL BayesPrism run completed\n")
cat("============================================================\n")

cat(
  "Elapsed seconds:",
  elapsed_seconds,
  "\n"
)

cat(
  "Elapsed hours:",
  elapsed_seconds / 3600,
  "\n\n"
)


# ==============================================================================
# 11. Validate final BayesPrism object
# ==============================================================================

assert_true(
  inherits(bp, "BayesPrism"),
  "run.prism() did not return a BayesPrism object."
)

assert_true(
  isTRUE(
    bp@control_param$update.gibbs
  ),
  "Returned result does not record update.gibbs = TRUE."
)

assert_true(
  !is.null(bp@reference.update),
  "Updated reference is missing."
)

assert_true(
  !is.null(bp@posterior.theta_f),
  "Final posterior theta is missing."
)


# ==============================================================================
# 12. Extract posterior results
# ==============================================================================

initial_state <- bp@posterior.initial.cellState
initial_type <- bp@posterior.initial.cellType

Z_state <- initial_state@Z
Z_type <- initial_type@Z

theta_state <- initial_state@theta
theta_initial <- initial_type@theta

theta_final <- bp@posterior.theta_f@theta


# ==============================================================================
# 13. Structural validation
# ==============================================================================

assert_true(
  identical(
    dim(Z_type),
    c(
      2490L,
      13668L,
      7L
    )
  ),
  paste(
    "Unexpected initial cell-type Z dimensions:",
    paste(
      dim(Z_type),
      collapse = " x "
    )
  )
)

assert_true(
  identical(
    dim(Z_state),
    c(
      2490L,
      13668L,
      7L
    )
  ),
  paste(
    "Unexpected initial cell-state Z dimensions:",
    paste(
      dim(Z_state),
      collapse = " x "
    )
  )
)

assert_true(
  nrow(theta_initial) == 2490L &&
    ncol(theta_initial) == 7L,
  "Initial theta is not 2490 x 7."
)

assert_true(
  nrow(theta_final) == 2490L &&
    ncol(theta_final) == 7L,
  "Final theta is not 2490 x 7."
)

assert_true(
  identical(
    rownames(theta_initial),
    rownames(theta_final)
  ),
  "Initial and final theta sample order differs."
)

assert_true(
  identical(
    colnames(theta_initial),
    colnames(theta_final)
  ),
  "Initial and final theta component order differs."
)

assert_true(
  setequal(
    colnames(theta_initial),
    EXPECTED_COMPONENTS
  ),
  "Initial theta epithelial identities differ from expectations."
)

assert_true(
  setequal(
    colnames(theta_final),
    EXPECTED_COMPONENTS
  ),
  "Final theta epithelial identities differ from expectations."
)


# ==============================================================================
# 14. Theta sanity checks
# ==============================================================================

initial_rowsum_error <- max(
  abs(
    rowSums(theta_initial) - 1
  )
)

final_rowsum_error <- max(
  abs(
    rowSums(theta_final) - 1
  )
)


assert_true(
  initial_rowsum_error < 1e-8,
  paste0(
    "Initial theta row-sum error = ",
    initial_rowsum_error
  )
)

assert_true(
  final_rowsum_error < 1e-8,
  paste0(
    "Final theta row-sum error = ",
    final_rowsum_error
  )
)

assert_true(
  all(is.finite(theta_initial)),
  "Initial theta contains non-finite values."
)

assert_true(
  all(is.finite(theta_final)),
  "Final theta contains non-finite values."
)

assert_true(
  all(theta_initial >= 0),
  "Initial theta contains negative values."
)

assert_true(
  all(theta_final >= 0),
  "Final theta contains negative values."
)


# ==============================================================================
# 15. Initial Z read-mass conservation
# ==============================================================================

cat("Auditing Stage 2 posterior read-mass conservation...\n")


Z_sum <- apply(
  Z_type,
  c(1, 2),
  sum
)


mixture_input <- bp@prism@mixture


assert_true(
  identical(
    dim(Z_sum),
    dim(mixture_input)
  ),
  "Summed Stage 2 Z dimensions differ from pseudo-mixture."
)


entry_abs_diff <- abs(
  Z_sum - mixture_input
)


max_entry_abs_diff <- max(
  entry_abs_diff
)

mean_entry_abs_diff <- mean(
  entry_abs_diff
)


input_mass <- rowSums(
  mixture_input
)

allocated_mass <- rowSums(
  Z_sum
)


mass_difference <- allocated_mass -
  input_mass


relative_error <- abs(
  mass_difference
) / input_mass


median_relative_error <- median(
  relative_error
)

p95_relative_error <- as.numeric(
  quantile(
    relative_error,
    probs = 0.95,
    names = FALSE
  )
)

max_relative_error <- max(
  relative_error
)


conservation_table <- data.frame(

  sample_id = rownames(
    mixture_input
  ),

  input_epithelial_mass =
    input_mass,

  allocated_stage2_mass =
    allocated_mass,

  mass_difference =
    mass_difference,

  relative_error =
    relative_error,

  stringsAsFactors = FALSE
)


write_tsv(
  conservation_table,
  OUT_CONSERVATION_SAMPLE
)


cat("\nRead-mass conservation:\n")

cat(
  "  max entry absolute difference :",
  max_entry_abs_diff,
  "\n"
)

cat(
  "  mean entry absolute difference:",
  mean_entry_abs_diff,
  "\n"
)

cat(
  "  median sample relative error  :",
  median_relative_error,
  "\n"
)

cat(
  "  P95 sample relative error     :",
  p95_relative_error,
  "\n"
)

cat(
  "  max sample relative error     :",
  max_relative_error,
  "\n\n"
)


# ==============================================================================
# 16. Save fractions
# ==============================================================================

initial_fraction_out <- data.frame(

  sample_id =
    rownames(theta_initial),

  theta_initial,

  check.names = FALSE,

  stringsAsFactors = FALSE
)


final_fraction_out <- data.frame(

  sample_id =
    rownames(theta_final),

  theta_final,

  check.names = FALSE,

  stringsAsFactors = FALSE
)


write_tsv(
  initial_fraction_out,
  OUT_INITIAL_FRACTION
)

write_tsv(
  final_fraction_out,
  OUT_FINAL_FRACTION
)


cat("[PASS] Saved initial and final epithelial fractions\n")


# ==============================================================================
# 17. Save seven epithelial subtype Z matrices
# ==============================================================================

component_names <- dimnames(
  Z_type
)[[3]]


assert_true(
  length(component_names) == 7L,
  "Z_type does not contain exactly seven named components."
)


cat("\nSaving epithelial subtype Z matrices...\n")


expression_files <- character(
  length(component_names)
)


for (i in seq_along(component_names)) {

  component <- component_names[i]

  component_matrix <- Z_type[
    ,
    ,
    component,
    drop = TRUE
  ]

  assert_true(
    nrow(component_matrix) == 2490L,
    paste(
      component,
      "Z matrix sample count mismatch."
    )
  )

  assert_true(
    ncol(component_matrix) == 13668L,
    paste(
      component,
      "Z matrix gene count mismatch."
    )
  )

  out_file <- file.path(
    EXPRESSION_DIR,
    paste0(
      "bayesprism_expression_",
      component,
      ".rds"
    )
  )

  saveRDS(
    component_matrix,
    out_file,
    compress = FALSE
  )

  expression_files[i] <- out_file

  cat(
    "  [SAVED]",
    component,
    "\n"
  )

  rm(component_matrix)
  gc(verbose = FALSE)
}


# ==============================================================================
# 18. Epithelial component summary
# ==============================================================================

component_summary <- data.frame(

  component =
    colnames(theta_final),

  initial_mean =
    colMeans(theta_initial)[
      colnames(theta_final)
    ],

  initial_median =
    apply(
      theta_initial[
        ,
        colnames(theta_final),
        drop = FALSE
      ],
      2,
      median
    ),

  final_mean =
    colMeans(theta_final),

  final_median =
    apply(
      theta_final,
      2,
      median
    ),

  stringsAsFactors = FALSE
)


write_tsv(
  component_summary,
  OUT_COMPONENT_SUMMARY
)


cat("\nEpithelial component fractions:\n")

print(
  component_summary,
  row.names = FALSE
)


# ==============================================================================
# 19. Save updated reference
# ==============================================================================

saveRDS(
  bp@reference.update,
  OUT_REFERENCE_UPDATE,
  compress = TRUE
)

cat(
  "\n[PASS] Saved updated Stage 2 reference\n"
)


# ==============================================================================
# 20. Save complete BayesPrism result
#
# Do this after the run has successfully completed.
# compress = FALSE reduces post-run serialization overhead.
# ==============================================================================

cat("\nSaving complete Stage 2 BayesPrism result...\n")

saveRDS(
  bp,
  OUT_RESULT,
  compress = FALSE
)

cat(
  "[PASS] Saved complete Stage 2 BayesPrism result:\n",
  OUT_RESULT,
  "\n\n"
)


# ==============================================================================
# 21. Formal audit table
# ==============================================================================

audit <- do.call(
  rbind,
  list(

    audit_row(
      "BayesPrism_class",
      inherits(bp, "BayesPrism"),
      TRUE,
      "PASS"
    ),

    audit_row(
      "update_gibbs",
      bp@control_param$update.gibbs,
      TRUE,
      ifelse(
        isTRUE(
          bp@control_param$update.gibbs
        ),
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "samples",
      nrow(theta_final),
      2490,
      ifelse(
        nrow(theta_final) == 2490L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "genes",
      dim(Z_type)[2],
      13668,
      ifelse(
        dim(Z_type)[2] == 13668L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "epithelial_components",
      ncol(theta_final),
      7,
      ifelse(
        ncol(theta_final) == 7L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "initial_theta_rowsum_max_error",
      initial_rowsum_error,
      "<1e-8",
      ifelse(
        initial_rowsum_error < 1e-8,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "final_theta_rowsum_max_error",
      final_rowsum_error,
      "<1e-8",
      ifelse(
        final_rowsum_error < 1e-8,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "initial_theta_all_finite",
      all(is.finite(theta_initial)),
      TRUE,
      ifelse(
        all(is.finite(theta_initial)),
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "final_theta_all_finite",
      all(is.finite(theta_final)),
      TRUE,
      ifelse(
        all(is.finite(theta_final)),
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "initial_theta_nonnegative",
      all(theta_initial >= 0),
      TRUE,
      ifelse(
        all(theta_initial >= 0),
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "final_theta_nonnegative",
      all(theta_final >= 0),
      TRUE,
      ifelse(
        all(theta_final >= 0),
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "median_Z_mass_relative_error",
      median_relative_error,
      "informational",
      "INFO"
    ),

    audit_row(
      "p95_Z_mass_relative_error",
      p95_relative_error,
      "informational",
      "INFO"
    ),

    audit_row(
      "max_Z_mass_relative_error",
      max_relative_error,
      "informational",
      "INFO"
    ),

    audit_row(
      "elapsed_seconds",
      elapsed_seconds,
      "informational",
      "INFO"
    )
  )
)


write_tsv(
  audit,
  OUT_AUDIT
)


cat("\nFormal Stage 2 full-run audit:\n")

print(
  audit,
  row.names = FALSE
)


assert_true(
  !any(audit$status == "FAIL"),
  "Stage 2 full BayesPrism audit contains FAIL."
)


# ==============================================================================
# 22. Manifest
# ==============================================================================

manifest <- data.frame(

  role = c(
    "SOURCE_STAGE2_PRISM",
    "FULL_BAYESPRISM_RESULT",
    "UPDATED_REFERENCE",
    "INITIAL_FRACTIONS",
    "FINAL_FRACTIONS",
    "READ_MASS_CONSERVATION",
    "COMPONENT_SUMMARY",
    "FULL_RUN_AUDIT",
    "CONTROL",
    "LOG",
    "SESSION_INFO"
  ),

  path = c(
    PRISM_RDS,
    OUT_RESULT,
    OUT_REFERENCE_UPDATE,
    OUT_INITIAL_FRACTION,
    OUT_FINAL_FRACTION,
    OUT_CONSERVATION_SAMPLE,
    OUT_COMPONENT_SUMMARY,
    OUT_AUDIT,
    OUT_CONTROL,
    OUT_LOG,
    OUT_SESSIONINFO
  ),

  stringsAsFactors = FALSE
)


write_tsv(
  manifest,
  OUT_MANIFEST
)


# ==============================================================================
# 23. sessionInfo
# ==============================================================================

capture.output(
  sessionInfo(),
  file = OUT_SESSIONINFO
)


# ==============================================================================
# 24. Final summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("S2-C FULL RUN SUMMARY\n")
cat("============================================================\n")

cat(
  "Samples                         :",
  nrow(theta_final),
  "\n"
)

cat(
  "Genes                           :",
  dim(Z_type)[2],
  "\n"
)

cat(
  "Epithelial components           :",
  ncol(theta_final),
  "\n"
)

cat(
  "Initial Gibbs                   : COMPLETED\n"
)

cat(
  "Reference update                : COMPLETED\n"
)

cat(
  "Final Gibbs                     : COMPLETED\n"
)

cat(
  "Initial theta row-sum max error :",
  initial_rowsum_error,
  "\n"
)

cat(
  "Final theta row-sum max error   :",
  final_rowsum_error,
  "\n"
)

cat(
  "Median Z mass relative error    :",
  median_relative_error,
  "\n"
)

cat(
  "P95 Z mass relative error       :",
  p95_relative_error,
  "\n"
)

cat(
  "Max Z mass relative error       :",
  max_relative_error,
  "\n"
)

cat(
  "Elapsed seconds                 :",
  elapsed_seconds,
  "\n"
)

cat(
  "Elapsed hours                   :",
  elapsed_seconds / 3600,
  "\n"
)

cat("\n")

cat(
  "IMPORTANT:\n",
  "Stage 2 BayesPrism computational run is complete.\n",
  "Biological QC has NOT yet been declared PASS.\n",
  "Do NOT start EcoTyper yet.\n"
)

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\n")
cat("============================================================\n")
cat("S2-C FULL BAYESPRISM RUN COMPLETED\n")
cat("END:", format(Sys.time()), "\n")
cat("============================================================\n")