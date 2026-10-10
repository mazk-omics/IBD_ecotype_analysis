#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

suppressPackageStartupMessages({
  library(BayesPrism)
})

# Export BayesPrism Stage1/Stage2 final cell-type fractions as plain TSV files.
#
# Usage:
#   Rscript 17_export_cell_type_fractions_tsv.R [PROJECT_ROOT]
#
# PROJECT_ROOT defaults to IBD_ECOTYPER_ROOT, then to the server project path.
# Stage2 theta is conditional on the epithelial compartment. Its whole-biopsy
# fraction is therefore:
#
#   Stage1 final Epithelial fraction * Stage2 final conditional fraction

args <- commandArgs(trailingOnly = TRUE)

default_root <- Sys.getenv(
  "IBD_ECOTYPER_ROOT",
  unset = "/home/mazekai/IBD_EcoTyper"
)

if (length(args) > 1L) {
  stop("Usage: Rscript 17_export_cell_type_fractions_tsv.R [PROJECT_ROOT]")
}

ROOT <- if (length(args) == 1L) args[[1L]] else default_root
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = FALSE)

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

for (path in c(S1_RDS, S2_RDS)) {
  if (!file.exists(path)) {
    stop("Required BayesPrism result does not exist: ", path)
  }
}

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

assert_fraction_matrix <- function(x, label, expected_ncol = NULL) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(label, " must be a numeric matrix")
  }
  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop(label, " is empty")
  }
  if (is.null(rownames(x)) || anyNA(rownames(x)) ||
      any(!nzchar(rownames(x))) || anyDuplicated(rownames(x))) {
    stop(label, " must have unique, non-empty sample row names")
  }
  if (is.null(colnames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(colnames(x))) || anyDuplicated(colnames(x))) {
    stop(label, " must have unique, non-empty cell-type column names")
  }
  if (!is.null(expected_ncol) && ncol(x) != expected_ncol) {
    stop(label, " has ", ncol(x), " columns; expected ", expected_ncol)
  }
  if (any(!is.finite(x))) {
    stop(label, " contains NA, NaN, or infinite values")
  }
  if (any(x < -1e-12)) {
    stop(label, " contains negative fractions below numerical tolerance")
  }

  closure_error <- max(abs(rowSums(x) - 1))
  if (closure_error >= 1e-8) {
    stop(
      label,
      " rows do not sum to 1; maximum absolute closure error = ",
      format(closure_error, scientific = TRUE)
    )
  }

  invisible(closure_error)
}

write_fraction_tsv <- function(x, path, scale = 1) {
  out <- data.frame(
    sample_id = rownames(x),
    x * scale,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  tmp <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(tmp), add = TRUE)

  write.table(
    out,
    file = tmp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace existing output: ", path)
  }
  if (!file.rename(tmp, path)) {
    stop("Could not move temporary output into place: ", path)
  }

  invisible(path)
}

cat("Project root: ", ROOT, "\n", sep = "")
cat("Loading Stage1: ", S1_RDS, "\n", sep = "")
bp1 <- readRDS(S1_RDS)
cat("Loading Stage2: ", S2_RDS, "\n", sep = "")
bp2 <- readRDS(S2_RDS)

# BayesPrism 2.2.3 updated/final cell-type fractions.
s1_final <- bp1@posterior.theta_f@theta
s2_final <- bp2@posterior.theta_f@theta

s1_error <- assert_fraction_matrix(
  s1_final,
  "Stage1 final theta",
  expected_ncol = 14L
)
s2_error <- assert_fraction_matrix(
  s2_final,
  "Stage2 final theta",
  expected_ncol = 7L
)

if (!("Epithelial" %in% colnames(s1_final))) {
  stop("Stage1 final theta is missing the Epithelial column")
}

epi7 <- c(
  "Colonic_absorptive",
  "Enteroendocrine",
  "Goblet",
  "Ileal_absorptive",
  "Paneth",
  "Stem_TA_progenitor",
  "Tuft"
)

if (!setequal(colnames(s2_final), epi7)) {
  stop(
    "Unexpected Stage2 cell types. Observed: ",
    paste(colnames(s2_final), collapse = ", ")
  )
}
s2_final <- s2_final[, epi7, drop = FALSE]

if (!setequal(rownames(s1_final), rownames(s2_final))) {
  only_s1 <- setdiff(rownames(s1_final), rownames(s2_final))
  only_s2 <- setdiff(rownames(s2_final), rownames(s1_final))
  stop(
    "Stage1/Stage2 sample sets differ. Stage1-only = ", length(only_s1),
    "; Stage2-only = ", length(only_s2)
  )
}
s2_final <- s2_final[rownames(s1_final), , drop = FALSE]

stage1_epi <- s1_final[, "Epithelial"]
s2_whole <- sweep(
  s2_final,
  MARGIN = 1L,
  STATS = stage1_epi,
  FUN = "*"
)

parent_error <- max(abs(rowSums(s2_whole) - stage1_epi))
if (parent_error >= 1e-8) {
  stop(
    "Stage2 whole-biopsy fractions do not reconstruct the Stage1 ",
    "Epithelial parent; maximum absolute error = ",
    format(parent_error, scientific = TRUE)
  )
}

non_epi13 <- setdiff(colnames(s1_final), "Epithelial")
full20 <- cbind(
  s1_final[, non_epi13, drop = FALSE],
  s2_whole
)

full20_error <- assert_fraction_matrix(
  full20,
  "Hierarchical full20 fraction",
  expected_ncol = 20L
)

outputs <- c(
  stage1_fraction = file.path(
    OUT,
    "bayesprism_stage1_broad14_final_fraction.tsv"
  ),
  stage2_conditional_fraction = file.path(
    OUT,
    "bayesprism_stage2_epithelial7_conditional_final_fraction.tsv"
  ),
  stage2_whole_biopsy_fraction = file.path(
    OUT,
    "bayesprism_stage2_epithelial7_whole_biopsy_fraction.tsv"
  ),
  hierarchical_full20_fraction = file.path(
    OUT,
    "bayesprism_hierarchical_full20_final_fraction.tsv"
  ),
  hierarchical_full20_percent = file.path(
    OUT,
    "bayesprism_hierarchical_full20_final_percent.tsv"
  )
)

write_fraction_tsv(s1_final, outputs[["stage1_fraction"]])
write_fraction_tsv(s2_final, outputs[["stage2_conditional_fraction"]])
write_fraction_tsv(s2_whole, outputs[["stage2_whole_biopsy_fraction"]])
write_fraction_tsv(full20, outputs[["hierarchical_full20_fraction"]])
write_fraction_tsv(
  full20,
  outputs[["hierarchical_full20_percent"]],
  scale = 100
)

qc <- data.frame(
  metric = c(
    "n_samples",
    "stage1_n_cell_types",
    "stage2_n_epithelial_subtypes",
    "full_n_cell_types",
    "stage1_closure_max_error",
    "stage2_conditional_closure_max_error",
    "stage2_parent_reconstruction_max_error",
    "full20_closure_max_error"
  ),
  value = c(
    nrow(full20),
    ncol(s1_final),
    ncol(s2_final),
    ncol(full20),
    s1_error,
    s2_error,
    parent_error,
    full20_error
  ),
  stringsAsFactors = FALSE
)

qc_path <- file.path(OUT, "bayesprism_cell_type_fraction_export_qc.tsv")
write.table(
  qc,
  file = qc_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nPASS: cell-type fraction TSV export completed.\n")
print(qc, row.names = FALSE)
cat("\nOutputs:\n")
cat(paste0("  ", unname(outputs), collapse = "\n"), "\n", sep = "")
cat("  ", qc_path, "\n", sep = "")
