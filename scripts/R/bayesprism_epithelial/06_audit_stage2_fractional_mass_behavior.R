#!/usr/bin/env Rscript

# ==============================================================================
# IBD EcoTyper
# Stage 2 epithelial BayesPrism
#
# QC-01:
# Audit how fractional Stage 1 epithelial posterior mass was handled
# internally by Stage 2 BayesPrism Gibbs sampling.
#
# NO model fitting.
# NO modification of existing results.
# ==============================================================================


options(
  stringsAsFactors = FALSE,
  warn = 1
)


# ==============================================================================
# 1. Paths
# ==============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "17_bayesprism_epithelial_input"
)

RUN_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "20_bayesprism_epithelial_full_run"
)

EXPRESSION_DIR <- file.path(
  RUN_DIR,
  "celltype_expression"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "21_bayesprism_epithelial_qc"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


MIXTURE_RDS <- file.path(
  INPUT_DIR,
  "bayesprism_stage2_epithelial_pseudomixture.rds"
)

CONSERVATION_TSV <- file.path(
  RUN_DIR,
  "stage2_read_mass_conservation_by_sample.tsv"
)


# ==============================================================================
# 2. Outputs
# ==============================================================================

OUT_SUMMARY <- file.path(
  OUTPUT_DIR,
  "stage2_fractional_mass_audit_summary.tsv"
)

OUT_SAMPLE <- file.path(
  OUTPUT_DIR,
  "stage2_fractional_mass_audit_by_sample.tsv"
)

OUT_WORST <- file.path(
  OUTPUT_DIR,
  "stage2_fractional_mass_worst_samples.tsv"
)

OUT_LOG <- file.path(
  OUTPUT_DIR,
  "stage2_fractional_mass_audit.log"
)

OUT_SESSION <- file.path(
  OUTPUT_DIR,
  "fractional_mass_sessionInfo.txt"
)


# ==============================================================================
# 3. Logging
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
cat("Stage 2 fractional pseudo-mixture forensic QC\n")
cat("============================================================\n\n")

cat("START:", format(Sys.time()), "\n\n")


# ==============================================================================
# 4. Helpers
# ==============================================================================

assert_true <- function(condition, message) {

  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
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


# ==============================================================================
# 5. Load original Stage 2 pseudo-mixture
# ==============================================================================

assert_true(
  file.exists(MIXTURE_RDS),
  paste(
    "Missing pseudo-mixture:",
    MIXTURE_RDS
  )
)

cat("Loading Stage 2 pseudo-mixture...\n")

mixture <- readRDS(
  MIXTURE_RDS
)

assert_true(
  nrow(mixture) == 2490L,
  "Unexpected pseudo-mixture sample count."
)

assert_true(
  ncol(mixture) == 13668L,
  "Unexpected pseudo-mixture gene count."
)

cat(
  "[PASS]",
  nrow(mixture),
  "samples x",
  ncol(mixture),
  "genes\n\n"
)


# ==============================================================================
# 6. Sum all 7 posterior subtype Z matrices
# ==============================================================================

COMPONENTS <- c(
  "Colonic_absorptive",
  "Enteroendocrine",
  "Goblet",
  "Ileal_absorptive",
  "Paneth",
  "Stem_TA_progenitor",
  "Tuft"
)


Z_sum <- matrix(
  0,
  nrow = nrow(mixture),
  ncol = ncol(mixture),
  dimnames = dimnames(mixture)
)


cat("Accumulating 7 Stage 2 Z matrices...\n")


for (component in COMPONENTS) {

  z_file <- file.path(
    EXPRESSION_DIR,
    paste0(
      "bayesprism_expression_",
      component,
      ".rds"
    )
  )

  assert_true(
    file.exists(z_file),
    paste(
      "Missing Z matrix:",
      z_file
    )
  )

  cat(
    "  loading:",
    component,
    "\n"
  )

  z <- readRDS(
    z_file
  )

  assert_true(
    identical(
      dim(z),
      dim(mixture)
    ),
    paste(
      component,
      "dimensions differ from pseudo-mixture."
    )
  )

  assert_true(
    identical(
      rownames(z),
      rownames(mixture)
    ),
    paste(
      component,
      "sample order mismatch."
    )
  )

  assert_true(
    identical(
      colnames(z),
      colnames(mixture)
    ),
    paste(
      component,
      "gene order mismatch."
    )
  )

  Z_sum <- Z_sum + z

  rm(z)

  gc(verbose = FALSE)
}


cat("\n[PASS] Seven posterior Z matrices accumulated\n\n")


# ==============================================================================
# 7. Construct candidate integerizations
# ==============================================================================

mixture_floor <- floor(
  mixture
)

mixture_round <- round(
  mixture
)

mixture_ceil <- ceiling(
  mixture
)


# ==============================================================================
# 8. Entry-wise concordance
# ==============================================================================

diff_original <- Z_sum - mixture

diff_floor <- Z_sum - mixture_floor

diff_round <- Z_sum - mixture_round

diff_ceil <- Z_sum - mixture_ceil


max_abs_original <- max(
  abs(diff_original)
)

max_abs_floor <- max(
  abs(diff_floor)
)

max_abs_round <- max(
  abs(diff_round)
)

max_abs_ceil <- max(
  abs(diff_ceil)
)


mean_abs_original <- mean(
  abs(diff_original)
)

mean_abs_floor <- mean(
  abs(diff_floor)
)

mean_abs_round <- mean(
  abs(diff_round)
)

mean_abs_ceil <- mean(
  abs(diff_ceil)
)


TOL <- 1e-10


prop_equal_original <- mean(
  abs(diff_original) < TOL
)

prop_equal_floor <- mean(
  abs(diff_floor) < TOL
)

prop_equal_round <- mean(
  abs(diff_round) < TOL
)

prop_equal_ceil <- mean(
  abs(diff_ceil) < TOL
)


# ==============================================================================
# 9. Fractional-part audit
# ==============================================================================

fractional_part <- mixture -
  floor(mixture)


fractional_n <- sum(
  fractional_part > TOL
)

mean_fractional_part <- mean(
  fractional_part[
    fractional_part > TOL
  ]
)

max_fractional_part <- max(
  fractional_part
)


# If Gibbs uses floor/truncation:
#
# mixture - sum(Z)
#
# should approximately equal the fractional part of mixture.

residual_expected_floor <- mixture -
  Z_sum

floor_residual_error <- residual_expected_floor -
  fractional_part


max_floor_residual_error <- max(
  abs(
    floor_residual_error
  )
)

mean_floor_residual_error <- mean(
  abs(
    floor_residual_error
  )
)


# ==============================================================================
# 10. Sample-wise mass audit
# ==============================================================================

original_mass <- rowSums(
  mixture
)

floor_mass <- rowSums(
  mixture_floor
)

round_mass <- rowSums(
  mixture_round
)

allocated_mass <- rowSums(
  Z_sum
)


lost_mass <- original_mass -
  allocated_mass


relative_mass_loss <- lost_mass /
  original_mass


floor_relative_mass_loss <- (
  original_mass -
    floor_mass
) / original_mass


sample_table <- data.frame(

  sample_id =
    rownames(mixture),

  original_mass =
    original_mass,

  floor_mass =
    floor_mass,

  round_mass =
    round_mass,

  posterior_allocated_mass =
    allocated_mass,

  mass_lost =
    lost_mass,

  relative_mass_loss =
    relative_mass_loss,

  expected_floor_relative_loss =
    floor_relative_mass_loss,

  allocated_minus_floor =
    allocated_mass -
      floor_mass,

  stringsAsFactors = FALSE
)


sample_table$abs_allocated_minus_floor <- abs(
  sample_table$allocated_minus_floor
)


# ==============================================================================
# 11. Add previous conservation results when available
# ==============================================================================

if (file.exists(CONSERVATION_TSV)) {

  old_qc <- read.delim(
    CONSERVATION_TSV,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (
    "sample_id" %in% colnames(old_qc) &&
    all(
      sample_table$sample_id %in%
        old_qc$sample_id
    )
  ) {

    idx <- match(
      sample_table$sample_id,
      old_qc$sample_id
    )

    extra_cols <- setdiff(
      colnames(old_qc),
      "sample_id"
    )

    for (nm in extra_cols) {

      new_col_name <- paste0(
        "previous_",
        nm
      )

      sample_table[[new_col_name]] <-
        old_qc[[nm]][idx]
    }

    cat(
      "[PASS] Previous conservation QC merged by sample_id\n\n"
    )

  } else {

    cat(
      "[WARN] Previous conservation table could not be safely merged.\n\n"
    )
  }

} else {

  cat(
    "[INFO] Previous conservation table not found; skipping merge.\n\n"
  )
}


# ==============================================================================
# 12. Worst samples
# ==============================================================================

sample_table <- sample_table[
  order(
    -sample_table$relative_mass_loss
  ),
  ,
  drop = FALSE
]


worst_samples <- head(
  sample_table,
  50
)


write_tsv(
  sample_table,
  OUT_SAMPLE
)

write_tsv(
  worst_samples,
  OUT_WORST
)


# ==============================================================================
# 13. Relationship between loss and epithelial mass
# ==============================================================================

spearman_mass_loss <- suppressWarnings(
  cor(
    original_mass,
    relative_mass_loss,
    method = "spearman"
  )
)


mass_quantiles <- quantile(
  original_mass,
  probs = c(
    0,
    0.001,
    0.01,
    0.05,
    0.5,
    0.95,
    0.99,
    1
  ),
  names = FALSE
)


# ==============================================================================
# 14. Formal summary
# ==============================================================================

summary_table <- data.frame(

  metric = c(

    "samples",
    "genes",
    "fractional_entries",

    "mean_fractional_part_nonzero",
    "max_fractional_part",

    "max_abs_difference_original",
    "mean_abs_difference_original",
    "proportion_exact_original",

    "max_abs_difference_floor",
    "mean_abs_difference_floor",
    "proportion_exact_floor",

    "max_abs_difference_round",
    "mean_abs_difference_round",
    "proportion_exact_round",

    "max_abs_difference_ceil",
    "mean_abs_difference_ceil",
    "proportion_exact_ceil",

    "max_floor_residual_error",
    "mean_floor_residual_error",

    "median_sample_relative_mass_loss",
    "p95_sample_relative_mass_loss",
    "max_sample_relative_mass_loss",

    "spearman_original_mass_vs_relative_loss",

    "original_mass_min",
    "original_mass_p0.1",
    "original_mass_p1",
    "original_mass_p5",
    "original_mass_median",
    "original_mass_p95",
    "original_mass_p99",
    "original_mass_max"
  ),

  value = c(

    nrow(mixture),
    ncol(mixture),
    fractional_n,

    mean_fractional_part,
    max_fractional_part,

    max_abs_original,
    mean_abs_original,
    prop_equal_original,

    max_abs_floor,
    mean_abs_floor,
    prop_equal_floor,

    max_abs_round,
    mean_abs_round,
    prop_equal_round,

    max_abs_ceil,
    mean_abs_ceil,
    prop_equal_ceil,

    max_floor_residual_error,
    mean_floor_residual_error,

    median(
      relative_mass_loss
    ),

    as.numeric(
      quantile(
        relative_mass_loss,
        0.95,
        names = FALSE
      )
    ),

    max(
      relative_mass_loss
    ),

    spearman_mass_loss,

    mass_quantiles
  ),

  stringsAsFactors = FALSE
)


write_tsv(
  summary_table,
  OUT_SUMMARY
)


# ==============================================================================
# 15. Console report
# ==============================================================================

cat("============================================================\n")
cat("ENTRY-WISE CONCORDANCE\n")
cat("============================================================\n\n")

cat(
  "Against ORIGINAL pseudo-mixture:\n",
  "  max abs diff  =", max_abs_original, "\n",
  "  mean abs diff =", mean_abs_original, "\n",
  "  exact prop    =", prop_equal_original, "\n\n"
)

cat(
  "Against FLOOR(pseudo-mixture):\n",
  "  max abs diff  =", max_abs_floor, "\n",
  "  mean abs diff =", mean_abs_floor, "\n",
  "  exact prop    =", prop_equal_floor, "\n\n"
)

cat(
  "Against ROUND(pseudo-mixture):\n",
  "  max abs diff  =", max_abs_round, "\n",
  "  mean abs diff =", mean_abs_round, "\n",
  "  exact prop    =", prop_equal_round, "\n\n"
)

cat(
  "Against CEILING(pseudo-mixture):\n",
  "  max abs diff  =", max_abs_ceil, "\n",
  "  mean abs diff =", mean_abs_ceil, "\n",
  "  exact prop    =", prop_equal_ceil, "\n\n"
)


cat("============================================================\n")
cat("FRACTIONAL MASS\n")
cat("============================================================\n\n")

cat(
  "Fractional entries                  :",
  fractional_n,
  "\n"
)

cat(
  "Mean non-zero fractional part       :",
  mean_fractional_part,
  "\n"
)

cat(
  "Max fractional part                 :",
  max_fractional_part,
  "\n"
)

cat(
  "Max residual error vs fractional part:",
  max_floor_residual_error,
  "\n\n"
)


cat("============================================================\n")
cat("SAMPLE-WISE MASS LOSS\n")
cat("============================================================\n\n")

cat(
  "Median relative loss :",
  median(relative_mass_loss),
  "\n"
)

cat(
  "P95 relative loss    :",
  quantile(
    relative_mass_loss,
    0.95,
    names = FALSE
  ),
  "\n"
)

cat(
  "Max relative loss    :",
  max(relative_mass_loss),
  "\n"
)

cat(
  "Spearman input mass vs relative loss:",
  spearman_mass_loss,
  "\n\n"
)


cat("Worst 20 samples:\n\n")

print(
  head(
    worst_samples[
      ,
      c(
        "sample_id",
        "original_mass",
        "floor_mass",
        "posterior_allocated_mass",
        "mass_lost",
        "relative_mass_loss"
      )
    ],
    20
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Automatic interpretation flag
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("AUTOMATIC DIAGNOSTIC\n")
cat("============================================================\n\n")


if (
  max_abs_floor < 1e-8 &&
  prop_equal_floor > 0.999999
) {

  cat(
    "[STRONG EVIDENCE] Sum of Stage 2 subtype Z equals ",
    "FLOOR(pseudo-mixture).\n"
  )

  cat(
    "Fractional posterior mass was effectively discarded ",
    "at the multinomial sampling step.\n"
  )

} else {

  cat(
    "[NOT RESOLVED] Z sum does not exactly match simple floor/truncation.\n"
  )

  cat(
    "Further investigation of BayesPrism numerical handling is required.\n"
  )
}


# ==============================================================================
# 17. sessionInfo
# ==============================================================================

capture.output(
  sessionInfo(),
  file = OUT_SESSION
)


cat("\n")
cat("Output directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\nEND:", format(Sys.time()), "\n")