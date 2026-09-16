#!/usr/bin/env Rscript

# ==============================================================================
# IBD EcoTyper
# Stage 2 epithelial-only BayesPrism
#
# S2-B2a:
# Audit the gene filtering that new.prism() WOULD apply to the
# Stage 1 Epithelial-Z pseudo-mixture.
#
# IMPORTANT:
#   - does NOT modify the mixture
#   - does NOT construct a prism object
#   - does NOT run Gibbs sampling
# ==============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(BayesPrism)
})


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

PREPROCESS_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "18_bayesprism_epithelial_preprocessing"
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
  PREPROCESS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


MIXTURE_RDS <- file.path(
  INPUT_DIR,
  "bayesprism_stage2_epithelial_pseudomixture.rds"
)

COMMON_GENES_TXT <- file.path(
  PREPROCESS_DIR,
  "stage2_reference_mixture_common_genes.txt"
)


# ==============================================================================
# 2. Outputs
# ==============================================================================

OUT_ALL_GENES <- file.path(
  PREPROCESS_DIR,
  "stage2_newprism_default_outlier_gene_audit.tsv"
)

OUT_OUTLIERS <- file.path(
  PREPROCESS_DIR,
  "stage2_newprism_default_outliers.tsv"
)

OUT_SUMMARY <- file.path(
  PREPROCESS_DIR,
  "stage2_newprism_default_outlier_summary.tsv"
)

OUT_RETAINED <- file.path(
  PREPROCESS_DIR,
  "stage2_newprism_expected_retained_genes.txt"
)

OUT_LOG <- file.path(
  PREPROCESS_DIR,
  "stage2_newprism_outlier_audit.log"
)

OUT_SESSIONINFO <- file.path(
  PREPROCESS_DIR,
  "stage2_newprism_outlier_sessionInfo.txt"
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
cat("Stage 2 epithelial BayesPrism\n")
cat("S2-B2a: new.prism default outlier audit\n")
cat("============================================================\n\n")

cat("START:", format(Sys.time()), "\n")
cat(
  "BayesPrism version:",
  as.character(packageVersion("BayesPrism")),
  "\n\n"
)


# ==============================================================================
# 4. Helpers
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


# ==============================================================================
# 5. Validate inputs
# ==============================================================================

assert_true(
  file.exists(MIXTURE_RDS),
  paste(
    "Missing Stage 2 pseudo-mixture:",
    MIXTURE_RDS
  )
)

assert_true(
  file.exists(COMMON_GENES_TXT),
  paste(
    "Missing S2-B1 common gene list:",
    COMMON_GENES_TXT
  )
)


# ==============================================================================
# 6. Load frozen epithelial pseudo-mixture
# ==============================================================================

cat("Loading Stage 1 Epithelial-Z pseudo-mixture...\n")

mixture <- readRDS(
  MIXTURE_RDS
)

assert_true(
  is.matrix(mixture) ||
    inherits(mixture, "Matrix"),
  "Pseudo-mixture is not matrix-like."
)

assert_true(
  nrow(mixture) == 2490L,
  paste0(
    "Unexpected sample count: ",
    nrow(mixture)
  )
)

assert_true(
  ncol(mixture) == 13668L,
  paste0(
    "Unexpected gene count: ",
    ncol(mixture)
  )
)

assert_true(
  !is.null(rownames(mixture)),
  "Pseudo-mixture has no sample IDs."
)

assert_true(
  !is.null(colnames(mixture)),
  "Pseudo-mixture has no gene names."
)

assert_true(
  anyDuplicated(colnames(mixture)) == 0L,
  "Duplicated pseudo-mixture genes."
)


cat(
  "[PASS] Pseudo-mixture:",
  nrow(mixture),
  "samples x",
  ncol(mixture),
  "genes\n\n"
)


# ==============================================================================
# 7. Confirm all mixture genes are covered by cleaned epithelial reference
# ==============================================================================

common_genes <- readLines(
  COMMON_GENES_TXT
)

common_genes <- common_genes[
  nzchar(common_genes)
]

assert_true(
  length(common_genes) == 13668L,
  paste0(
    "Unexpected S2-B1 common gene count: ",
    length(common_genes)
  )
)

assert_true(
  setequal(
    colnames(mixture),
    common_genes
  ),
  paste0(
    "Pseudo-mixture gene universe differs from ",
    "S2-B1 reference/mixture common gene universe."
  )
)

cat(
  "[PASS] All 13,668 mixture genes are represented ",
  "in the cleaned epithelial reference\n\n"
)


# ==============================================================================
# 8. Reproduce BayesPrism new.prism() default outlier logic
# ==============================================================================

OUTLIER_CUT <- 0.01
OUTLIER_FRACTION <- 0.10

cat("Simulating BayesPrism new.prism() default outlier filter:\n")
cat("  outlier.cut      =", OUTLIER_CUT, "\n")
cat("  outlier.fraction =", OUTLIER_FRACTION, "\n\n")


# ------------------------------------------------------------------------------
# BayesPrism source logic:
#
# mixture.norm <- mixture / rowSums(mixture)
#
# outlier.idx <-
#   colSums(mixture.norm > outlier.cut) /
#   nrow(mixture.norm) >
#   outlier.fraction
#
# We reproduce this exactly, but calculate in gene blocks
# to avoid unnecessarily retaining another ~2490 x 13668 dense matrix.
# ------------------------------------------------------------------------------

mixture_totals <- rowSums(
  mixture
)

assert_true(
  all(is.finite(mixture_totals)),
  "Non-finite pseudo-mixture row sums detected."
)

assert_true(
  all(mixture_totals > 0),
  "One or more samples have non-positive epithelial total mass."
)


n_samples <- nrow(
  mixture
)

n_genes <- ncol(
  mixture
)

genes <- colnames(
  mixture
)


# ==============================================================================
# 9. Per-gene statistics
# ==============================================================================

n_gt_cut <- integer(
  n_genes
)

mean_fraction <- numeric(
  n_genes
)

max_fraction <- numeric(
  n_genes
)


BLOCK_SIZE <- 500L

starts <- seq.int(
  1L,
  n_genes,
  by = BLOCK_SIZE
)

cat(
  "Auditing",
  n_genes,
  "genes in",
  length(starts),
  "blocks...\n"
)


for (b in seq_along(starts)) {

  s <- starts[b]

  e <- min(
    s + BLOCK_SIZE - 1L,
    n_genes
  )

  idx <- s:e

  block <- mixture[
    ,
    idx,
    drop = FALSE
  ]

  block_fraction <- sweep(
    block,
    1,
    mixture_totals,
    "/"
  )

  n_gt_cut[idx] <- colSums(
    block_fraction > OUTLIER_CUT
  )

  mean_fraction[idx] <- colMeans(
    block_fraction
  )

  max_fraction[idx] <- apply(
    block_fraction,
    2,
    max
  )

  rm(
    block,
    block_fraction
  )

  if (
    b %% 5L == 0L ||
    b == length(starts)
  ) {

    cat(
      "  completed block",
      b,
      "/",
      length(starts),
      "\n"
    )

    gc(
      verbose = FALSE
    )
  }
}


# ==============================================================================
# 10. Apply exact default criterion
# ==============================================================================

fraction_samples_gt_cut <-
  n_gt_cut /
  n_samples


default_outlier <-
  fraction_samples_gt_cut >
  OUTLIER_FRACTION


audit_table <- data.frame(

  gene = genes,

  n_samples_gt_1pct =
    n_gt_cut,

  fraction_samples_gt_1pct =
    fraction_samples_gt_cut,

  mean_epithelial_fraction =
    mean_fraction,

  max_epithelial_fraction =
    max_fraction,

  default_newprism_outlier =
    default_outlier,

  stringsAsFactors = FALSE
)


# ==============================================================================
# 11. Sort and save
# ==============================================================================

audit_table <- audit_table[
  order(
    -audit_table$fraction_samples_gt_1pct,
    -audit_table$mean_epithelial_fraction,
    audit_table$gene
  ),
  ,
  drop = FALSE
]

write_tsv(
  audit_table,
  OUT_ALL_GENES
)


outlier_table <- audit_table[
  audit_table$default_newprism_outlier,
  ,
  drop = FALSE
]

write_tsv(
  outlier_table,
  OUT_OUTLIERS
)


retained_genes <- genes[
  !default_outlier
]

writeLines(
  retained_genes,
  OUT_RETAINED
)


# ==============================================================================
# 12. Summary
# ==============================================================================

summary_table <- data.frame(

  metric = c(
    "BayesPrism_version",
    "samples",
    "input_genes",
    "outlier_cut",
    "outlier_fraction",
    "default_outlier_genes",
    "expected_retained_genes"
  ),

  value = c(
    as.character(
      packageVersion("BayesPrism")
    ),
    n_samples,
    n_genes,
    OUTLIER_CUT,
    OUTLIER_FRACTION,
    sum(default_outlier),
    length(retained_genes)
  ),

  stringsAsFactors = FALSE
)

write_tsv(
  summary_table,
  OUT_SUMMARY
)


# ==============================================================================
# 13. sessionInfo
# ==============================================================================

capture.output(
  sessionInfo(),
  file = OUT_SESSIONINFO
)


# ==============================================================================
# 14. Console summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("S2-B2a SUMMARY\n")
cat("============================================================\n")

cat(
  "Input Stage1 Epithelial-Z genes :",
  n_genes,
  "\n"
)

cat(
  "Default new.prism outlier genes :",
  sum(default_outlier),
  "\n"
)

cat(
  "Expected retained genes         :",
  length(retained_genes),
  "\n"
)

cat("\n")


if (sum(default_outlier) == 0L) {

  cat(
    "[INFO] Default new.prism filtering would remove ZERO genes.\n"
  )

} else {

  cat(
    "Genes that default new.prism would remove:\n\n"
  )

  print(
    outlier_table,
    row.names = FALSE
  )
}


cat("\nIMPORTANT:\n")
cat(
  "No genes have actually been removed from the frozen pseudo-mixture.\n"
)
cat(
  "new.prism() has NOT been run.\n"
)
cat(
  "run.prism() has NOT been run.\n"
)

cat("\nOutput directory:\n")
cat(PREPROCESS_DIR, "\n")

cat("\n")
cat("============================================================\n")
cat("S2-B2a OUTLIER AUDIT COMPLETED\n")
cat("END:", format(Sys.time()), "\n")
cat("============================================================\n")