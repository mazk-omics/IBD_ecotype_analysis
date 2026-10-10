#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

# Summarize the fraction TSV files produced by:
#   17_export_cell_type_fractions_tsv.R
#
# Usage:
#   Rscript 18_summarize_cell_type_fractions.R [PROJECT_ROOT]
#
# Main output:
#   bayesprism_hierarchical_full20_fraction_mean_median.tsv
#
# Each output row is one cell type. Fractions are reported on both the
# 0-to-1 scale and the 0-to-100 percent scale.

args <- commandArgs(trailingOnly = TRUE)

default_root <- Sys.getenv(
  "IBD_ECOTYPER_ROOT",
  unset = "/home/mazekai/IBD_EcoTyper"
)

if (length(args) > 1L) {
  stop("Usage: Rscript 18_summarize_cell_type_fractions.R [PROJECT_ROOT]")
}

ROOT <- if (length(args) == 1L) args[[1L]] else default_root
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = FALSE)

OUT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/22_hierarchical_fraction_merge"
)

input_files <- c(
  stage1_broad14 = file.path(
    OUT,
    "bayesprism_stage1_broad14_final_fraction.tsv"
  ),
  stage2_epithelial7_conditional = file.path(
    OUT,
    "bayesprism_stage2_epithelial7_conditional_final_fraction.tsv"
  ),
  stage2_epithelial7_whole_biopsy = file.path(
    OUT,
    "bayesprism_stage2_epithelial7_whole_biopsy_fraction.tsv"
  ),
  hierarchical_full20 = file.path(
    OUT,
    "bayesprism_hierarchical_full20_final_fraction.tsv"
  )
)

missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing fraction TSV file(s). Run 17_export_cell_type_fractions_tsv.R ",
    "first:\n  ",
    paste(unname(missing_inputs), collapse = "\n  ")
  )
}

read_fraction_tsv <- function(path, scope) {
  x <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(x) == 0L || ncol(x) < 2L) {
    stop(scope, ": input table is empty or has no fraction columns: ", path)
  }
  if (!identical(names(x)[1L], "sample_id")) {
    stop(scope, ": first column must be named sample_id: ", path)
  }
  if (anyNA(x$sample_id) || any(!nzchar(x$sample_id)) ||
      anyDuplicated(x$sample_id)) {
    stop(scope, ": sample_id values must be unique and non-empty")
  }

  fraction_names <- names(x)[-1L]
  if (anyDuplicated(fraction_names) || any(!nzchar(fraction_names))) {
    stop(scope, ": fraction column names must be unique and non-empty")
  }

  is_numeric <- vapply(x[-1L], is.numeric, logical(1L))
  if (!all(is_numeric)) {
    stop(
      scope,
      ": non-numeric fraction column(s): ",
      paste(fraction_names[!is_numeric], collapse = ", ")
    )
  }

  values <- as.matrix(x[-1L])
  storage.mode(values) <- "double"
  rownames(values) <- x$sample_id

  if (any(!is.finite(values))) {
    stop(scope, ": fraction matrix contains NA, NaN, or infinite values")
  }
  if (any(values < -1e-12) || any(values > 1 + 1e-12)) {
    stop(scope, ": fraction values fall outside the expected [0, 1] range")
  }

  values
}

summarize_fraction_matrix <- function(x, scope, denominator) {
  means <- colMeans(x)
  medians <- apply(x, 2L, median)

  data.frame(
    fraction_scope = scope,
    denominator = denominator,
    cell_type = colnames(x),
    n_samples = nrow(x),
    mean_fraction = unname(means),
    median_fraction = unname(medians),
    mean_percent = unname(means * 100),
    median_percent = unname(medians * 100),
    rank_by_mean = rank(-means, ties.method = "min"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_tsv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  tmp <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(tmp), add = TRUE)

  write.table(
    x,
    file = tmp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  if (!file.copy(tmp, path, overwrite = TRUE)) {
    stop("Could not write output: ", path)
  }

  invisible(path)
}

denominators <- c(
  stage1_broad14 = "whole_biopsy",
  stage2_epithelial7_conditional = "epithelial_compartment",
  stage2_epithelial7_whole_biopsy = "whole_biopsy",
  hierarchical_full20 = "whole_biopsy"
)

fraction_matrices <- Map(
  function(path, scope) read_fraction_tsv(path, scope),
  path = unname(input_files),
  scope = names(input_files)
)
names(fraction_matrices) <- names(input_files)

# Closure expectations differ by scope:
# - Stage1, Stage2 conditional, and full20 each sum to 1.
# - Stage2 whole-biopsy subtypes sum to the Stage1 Epithelial fraction,
#   so they must not be independently renormalized to 1.
closed_scopes <- c(
  "stage1_broad14",
  "stage2_epithelial7_conditional",
  "hierarchical_full20"
)

for (scope in closed_scopes) {
  closure_error <- max(abs(rowSums(fraction_matrices[[scope]]) - 1))
  if (closure_error >= 1e-8) {
    stop(
      scope,
      ": rows do not sum to 1; maximum absolute closure error = ",
      format(closure_error, scientific = TRUE)
    )
  }
}

all_summary <- do.call(
  rbind,
  lapply(names(fraction_matrices), function(scope) {
    summarize_fraction_matrix(
      fraction_matrices[[scope]],
      scope = scope,
      denominator = denominators[[scope]]
    )
  })
)
rownames(all_summary) <- NULL

full20_summary <- all_summary[
  all_summary$fraction_scope == "hierarchical_full20",
  c(
    "cell_type",
    "n_samples",
    "mean_fraction",
    "median_fraction",
    "mean_percent",
    "median_percent",
    "rank_by_mean"
  ),
  drop = FALSE
]

# Display the most abundant cell types first while retaining an explicit rank.
full20_summary <- full20_summary[
  order(full20_summary$rank_by_mean, full20_summary$cell_type),
  ,
  drop = FALSE
]
rownames(full20_summary) <- NULL

full20_output <- file.path(
  OUT,
  "bayesprism_hierarchical_full20_fraction_mean_median.tsv"
)
all_output <- file.path(
  OUT,
  "bayesprism_all_fraction_tables_mean_median.tsv"
)

write_tsv_atomic(full20_summary, full20_output)
write_tsv_atomic(all_summary, all_output)

cat("PASS: fraction summary completed.\n")
cat("Samples summarized: ", nrow(fraction_matrices$hierarchical_full20), "\n", sep = "")
cat("Cell types in final full20 table: ", nrow(full20_summary), "\n", sep = "")
cat("\nFinal full20 mean/median table:\n")
print(full20_summary, row.names = FALSE)
cat("\nOutputs:\n")
cat("  ", full20_output, "\n", sep = "")
cat("  ", all_output, "\n", sep = "")
