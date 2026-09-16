#!/usr/bin/env Rscript

# ==============================================================================
# IBD EcoTyper Project
# Stage 2 epithelial-only BayesPrism
#
# S2-A: Build and audit formal Stage 2 inputs
#
# This script ONLY:
#   1. Creates Stage 2 directories
#   2. Extracts the 7 epithelial identities from the CLOSED unified reference
#   3. Loads Stage 1 Broad Epithelial posterior Z as pseudo-mixture
#   4. Audits dimensions / identities / sample order / genes / numeric properties
#   5. Saves frozen Stage 2 input objects and audit files
#
# This script DOES NOT:
#   - reselect cells
#   - downsample cells
#   - rebuild the reference corpus
#   - run cleanup.genes()
#   - run new.prism()
#   - run run.prism()
#   - normalize / log-transform / round the pseudo-mixture
#
# Project:
#   /home/mazekai/IBD_EcoTyper
# ==============================================================================


# ==============================================================================
# 0. Basic settings
# ==============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(Matrix)
})


# ==============================================================================
# 1. Project paths
# ==============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

SCRIPT_DIR <- file.path(
  PROJECT_ROOT,
  "scripts",
  "R",
  "bayesprism_epithelial"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "17_bayesprism_epithelial_input"
)


# ------------------------------------------------------------------------------
# Existing CLOSED inputs
# ------------------------------------------------------------------------------

REFERENCE_RDS <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "11_bayesprism_full_input",
  "bayesprism_full_reference_count_matrix.rds"
)

EPITHELIAL_Z_RDS <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "14_bayesprism_broad_run",
  "celltype_expression",
  "bayesprism_expression_Epithelial.rds"
)

BROAD14_FRACTION_INITIAL <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "14_bayesprism_broad_run",
  "bayesprism_broad14_fraction_initial.tsv"
)


# ==============================================================================
# 2. Create directories
# ==============================================================================

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

cat("============================================================\n")
cat("Stage 2 epithelial BayesPrism input build\n")
cat("============================================================\n")
cat("Project root :", PROJECT_ROOT, "\n")
cat("Script dir   :", SCRIPT_DIR, "\n")
cat("Output dir   :", OUTPUT_DIR, "\n")
cat("\n")


# ==============================================================================
# 3. Output files
# ==============================================================================

OUT_REFERENCE_RDS <- file.path(
  OUTPUT_DIR,
  "bayesprism_epithelial7_reference_raw.rds"
)

OUT_MIXTURE_RDS <- file.path(
  OUTPUT_DIR,
  "bayesprism_stage2_epithelial_pseudomixture.rds"
)

OUT_REFERENCE_COUNTS <- file.path(
  OUTPUT_DIR,
  "stage2_reference_cell_counts.tsv"
)

OUT_REFERENCE_AUDIT <- file.path(
  OUTPUT_DIR,
  "stage2_reference_audit.tsv"
)

OUT_MIXTURE_AUDIT <- file.path(
  OUTPUT_DIR,
  "stage2_mixture_audit.tsv"
)

OUT_GENE_OVERLAP <- file.path(
  OUTPUT_DIR,
  "stage2_precleanup_gene_overlap.tsv"
)

OUT_GENE_SUMMARY <- file.path(
  OUTPUT_DIR,
  "stage2_precleanup_gene_overlap_summary.tsv"
)

OUT_MANIFEST <- file.path(
  OUTPUT_DIR,
  "stage2_input_manifest.tsv"
)

OUT_SESSIONINFO <- file.path(
  OUTPUT_DIR,
  "sessionInfo.txt"
)

OUT_LOG <- file.path(
  OUTPUT_DIR,
  "stage2_input_build.log"
)


# ==============================================================================
# 4. Logging
# ==============================================================================

log_con <- file(OUT_LOG, open = "wt")

sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

on.exit({
  while (sink.number(type = "message") > 0) {
    sink(type = "message")
  }

  while (sink.number(type = "output") > 0) {
    sink(type = "output")
  }

  close(log_con)
}, add = TRUE)

cat("\n============================================================\n")
cat("START:", format(Sys.time()), "\n")
cat("============================================================\n\n")


# ==============================================================================
# 5. Helper functions
# ==============================================================================

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required input file does not exist:\n",
      path,
      call. = FALSE
    )
  }
}


assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}


audit_row <- function(item, observed, expected = NA_character_, status) {
  data.frame(
    item = as.character(item),
    observed = as.character(observed),
    expected = as.character(expected),
    status = as.character(status),
    stringsAsFactors = FALSE
  )
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


# Safely inspect numerical values without unnecessarily densifying sparse matrices.
numeric_matrix_stats <- function(x, chunk_size = 100L) {

  nr <- nrow(x)

  negative_n <- 0
  nonfinite_n <- 0
  na_n <- 0
  fractional_n <- 0

  min_value <- Inf
  max_value <- -Inf

  starts <- seq.int(1L, nr, by = chunk_size)

  for (s in starts) {

    e <- min(s + chunk_size - 1L, nr)

    block <- x[s:e, , drop = FALSE]

    if (inherits(block, "sparseMatrix")) {
      vals <- block@x
    } else {
      vals <- as.vector(block)
    }

    if (length(vals) == 0L) {
      next
    }

    na_here <- is.na(vals)

    na_n <- na_n + sum(na_here)

    finite_here <- is.finite(vals)

    nonfinite_n <- nonfinite_n + sum(!finite_here)

    usable <- vals[finite_here]

    if (length(usable) > 0L) {

      negative_n <- negative_n + sum(usable < 0)

      fractional_n <- fractional_n +
        sum(abs(usable - round(usable)) > 1e-8)

      min_value <- min(min_value, usable)
      max_value <- max(max_value, usable)
    }
  }

  if (!is.finite(min_value)) min_value <- NA_real_
  if (!is.finite(max_value)) max_value <- NA_real_

  list(
    na_n = na_n,
    nonfinite_n = nonfinite_n,
    negative_n = negative_n,
    fractional_n = fractional_n,
    min_value = min_value,
    max_value = max_value
  )
}


# Stage 1 Epithelial Z is expected to be a matrix.
# This helper only provides conservative support if it was saved inside a
# one-matrix list.
unwrap_matrix <- function(obj, object_name) {

  if (
    is.matrix(obj) ||
    inherits(obj, "Matrix")
  ) {
    return(obj)
  }

  if (is.data.frame(obj)) {
    return(as.matrix(obj))
  }

  if (is.list(obj)) {

    matrix_like <- vapply(
      obj,
      function(z) {
        is.matrix(z) ||
          inherits(z, "Matrix") ||
          is.data.frame(z)
      },
      logical(1)
    )

    if (sum(matrix_like) == 1L) {

      nm <- names(obj)[which(matrix_like)]

      cat(
        object_name,
        "was stored as a list; using matrix-like element:",
        nm,
        "\n"
      )

      z <- obj[[which(matrix_like)]]

      if (is.data.frame(z)) {
        z <- as.matrix(z)
      }

      return(z)
    }
  }

  stop(
    object_name,
    " is not a recognized matrix-like object.",
    call. = FALSE
  )
}


# Attempt to recover sample order from the Stage 1 fraction table.
read_stage1_sample_ids <- function(path, expected_n) {

  if (!file.exists(path)) {
    return(NULL)
  }

  # Common write.table(..., row.names = TRUE) format
  x1 <- try(
    read.delim(
      path,
      header = TRUE,
      check.names = FALSE,
      row.names = 1
    ),
    silent = TRUE
  )

  if (!inherits(x1, "try-error")) {

    ids <- rownames(x1)

    if (
      length(ids) == expected_n &&
      !anyNA(ids) &&
      !anyDuplicated(ids)
    ) {
      return(ids)
    }
  }

  # Fallback: search an explicit sample-ID-like column
  x2 <- try(
    read.delim(
      path,
      header = TRUE,
      check.names = FALSE
    ),
    silent = TRUE
  )

  if (inherits(x2, "try-error")) {
    return(NULL)
  }

  preferred_names <- c(
    "sample",
    "Sample",
    "sample_id",
    "sampleID",
    "SampleID",
    "sample.id",
    "Sample.ID"
  )

  hit <- intersect(preferred_names, colnames(x2))

  if (length(hit) >= 1L) {

    ids <- as.character(x2[[hit[1]]])

    if (
      length(ids) == expected_n &&
      !anyNA(ids) &&
      !anyDuplicated(ids)
    ) {
      return(ids)
    }
  }

  # Conservative final fallback: first unique non-numeric column
  for (j in seq_len(ncol(x2))) {

    v <- x2[[j]]

    if (
      !is.numeric(v) &&
      length(v) == expected_n &&
      !anyNA(v) &&
      !anyDuplicated(v)
    ) {
      return(as.character(v))
    }
  }

  NULL
}


# ==============================================================================
# 6. Validate required source files
# ==============================================================================

cat("Checking source files...\n")

stop_if_missing(REFERENCE_RDS)
stop_if_missing(EPITHELIAL_Z_RDS)

cat("[PASS] Unified reference RDS exists\n")
cat("[PASS] Stage 1 Epithelial Z RDS exists\n")

if (file.exists(BROAD14_FRACTION_INITIAL)) {
  cat("[PASS] Broad14 initial fraction table exists\n")
} else {
  cat(
    "[WARN] Broad14 initial fraction table not found; ",
    "sample-order cross-check will be unavailable.\n"
  )
}

cat("\n")


# ==============================================================================
# 7. Load CLOSED unified reference
# ==============================================================================

cat("Loading CLOSED unified reference...\n")

reference_obj <- readRDS(REFERENCE_RDS)

assert_true(
  is.list(reference_obj),
  "Unified reference RDS is not a list."
)

required_fields <- c(
  "reference",
  "input.type",
  "cell.type.labels",
  "cell.state.labels",
  "metadata"
)

missing_fields <- setdiff(
  required_fields,
  names(reference_obj)
)

assert_true(
  length(missing_fields) == 0L,
  paste(
    "Unified reference is missing required fields:",
    paste(missing_fields, collapse = ", ")
  )
)

reference_raw <- reference_obj$reference

cat(
  "Reference class:",
  paste(class(reference_raw), collapse = ", "),
  "\n"
)

cat(
  "Reference dimensions:",
  nrow(reference_raw),
  "cells x",
  ncol(reference_raw),
  "genes\n"
)


# ------------------------------------------------------------------------------
# Locked source-reference expectations
# ------------------------------------------------------------------------------

assert_true(
  nrow(reference_raw) == 323144L,
  paste0(
    "Unexpected unified reference cell count: ",
    nrow(reference_raw),
    " != 323144"
  )
)

assert_true(
  ncol(reference_raw) == 14431L,
  paste0(
    "Unexpected unified reference gene count: ",
    ncol(reference_raw),
    " != 14431"
  )
)

assert_true(
  length(reference_obj$cell.type.labels) == nrow(reference_raw),
  "cell.type.labels length does not match reference rows."
)

assert_true(
  length(reference_obj$cell.state.labels) == nrow(reference_raw),
  "cell.state.labels length does not match reference rows."
)

assert_true(
  !is.null(rownames(reference_raw)),
  "Reference matrix has no cell rownames."
)

assert_true(
  !is.null(colnames(reference_raw)),
  "Reference matrix has no gene colnames."
)

assert_true(
  anyDuplicated(rownames(reference_raw)) == 0L,
  "Reference matrix contains duplicated cell IDs."
)

assert_true(
  anyDuplicated(colnames(reference_raw)) == 0L,
  "Reference matrix contains duplicated gene names."
)

cat("[PASS] Unified reference dimensions and identifiers match locked input\n\n")


# ==============================================================================
# 8. Define Stage 2 epithelial identities
# ==============================================================================

EPITHELIAL_IDENTITIES <- c(
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine",
  "Tuft",
  "Paneth",
  "Stem_TA_progenitor"
)

EXPECTED_COUNTS <- c(
  Colonic_absorptive = 19654L,
  Ileal_absorptive = 11343L,
  Goblet = 23023L,
  Enteroendocrine = 3751L,
  Tuft = 4523L,
  Paneth = 4062L,
  Stem_TA_progenitor = 27030L
)

source_type_labels <- as.character(
  reference_obj$cell.type.labels
)

missing_epi_identities <- setdiff(
  EPITHELIAL_IDENTITIES,
  unique(source_type_labels)
)

assert_true(
  length(missing_epi_identities) == 0L,
  paste(
    "Missing epithelial identities in unified reference:",
    paste(missing_epi_identities, collapse = ", ")
  )
)


# ==============================================================================
# 9. Extract Stage 2 epithelial reference
# ==============================================================================

cat("Extracting 7 epithelial identities...\n")

epi_idx <- which(
  source_type_labels %in% EPITHELIAL_IDENTITIES
)

assert_true(
  length(epi_idx) == 93386L,
  paste0(
    "Unexpected epithelial cell count: ",
    length(epi_idx),
    " != 93386"
  )
)

epithelial_reference <- reference_raw[
  epi_idx,
  ,
  drop = FALSE
]

epithelial_labels <- source_type_labels[
  epi_idx
]

names(epithelial_labels) <- rownames(
  epithelial_reference
)


# ------------------------------------------------------------------------------
# Stage 2 BayesPrism hierarchy is intentionally 1:1
# cell.type = cell.state
# ------------------------------------------------------------------------------

epithelial_cell_type_labels <- epithelial_labels
epithelial_cell_state_labels <- epithelial_labels


# ==============================================================================
# 10. Handle metadata safely
# ==============================================================================

source_metadata <- reference_obj$metadata

cat("Inspecting source metadata structure...\n")

cat(
  "Source metadata class:",
  paste(class(source_metadata), collapse = ", "),
  "\n"
)

if (is.list(source_metadata)) {
  cat(
    "Source metadata names:",
    paste(names(source_metadata), collapse = ", "),
    "\n"
  )
}

if (!is.null(dim(source_metadata))) {
  cat(
    "Source metadata dimensions:",
    paste(dim(source_metadata), collapse = " x "),
    "\n"
  )
}

# ------------------------------------------------------------------------------
# IMPORTANT
#
# The unified reference RDS "metadata" field is not assumed to be a per-cell
# metadata table.
#
# Stage 2 BayesPrism requires deterministic alignment of:
#   reference rows
#   cell.type.labels
#   cell.state.labels
#
# Therefore:
#
#   1. If source metadata can be aligned unambiguously, retain it.
#   2. Otherwise, construct a minimal deterministic per-cell metadata table.
#
# We DO NOT guess or force alignment.
# ------------------------------------------------------------------------------

metadata_alignment_mode <- NA_character_


# ------------------------------------------------------------------------------
# Case 1:
# metadata is row-wise and has exactly one row per unified reference cell
# ------------------------------------------------------------------------------

if (
  (is.data.frame(source_metadata) || is.matrix(source_metadata)) &&
  !is.null(nrow(source_metadata)) &&
  nrow(source_metadata) == nrow(reference_raw)
) {

  epithelial_metadata <- source_metadata[
    epi_idx,
    ,
    drop = FALSE
  ]

  rownames(epithelial_metadata) <- rownames(
    epithelial_reference
  )

  metadata_alignment_mode <- "subset_by_reference_row_index"

  cat(
    "[PASS] Source metadata aligned by reference row index\n"
  )


# ------------------------------------------------------------------------------
# Case 2:
# metadata has cell IDs as rownames and epithelial cells can be matched exactly
# ------------------------------------------------------------------------------

} else if (
  (is.data.frame(source_metadata) || is.matrix(source_metadata)) &&
  !is.null(rownames(source_metadata)) &&
  all(
    rownames(epithelial_reference) %in%
      rownames(source_metadata)
  )
) {

  epithelial_metadata <- source_metadata[
    rownames(epithelial_reference),
    ,
    drop = FALSE
  ]

  metadata_alignment_mode <- "subset_by_cell_id"

  cat(
    "[PASS] Source metadata aligned by cell ID\n"
  )


# ------------------------------------------------------------------------------
# Case 3:
# source metadata is reference-level/provenance metadata rather than
# a per-cell table.
#
# Build minimal deterministic cell metadata.
# ------------------------------------------------------------------------------

} else {

  epithelial_metadata <- data.frame(

    cell_id = rownames(
      epithelial_reference
    ),

    original_identity =
      as.character(epithelial_labels),

    stringsAsFactors = FALSE,

    row.names = rownames(
      epithelial_reference
    )
  )

  metadata_alignment_mode <-
    "minimal_cell_metadata_reconstructed"

  cat(
    "[INFO] Source metadata is not a deterministically alignable ",
    "per-cell table.\n"
  )

  cat(
    "[INFO] Constructed minimal Stage 2 cell metadata from ",
    "reference rownames + locked epithelial labels.\n"
  )
}


# ------------------------------------------------------------------------------
# Validate Stage 2 metadata
# ------------------------------------------------------------------------------

assert_true(
  nrow(epithelial_metadata) == nrow(epithelial_reference),
  paste0(
    "Stage 2 metadata row count mismatch: ",
    nrow(epithelial_metadata),
    " != ",
    nrow(epithelial_reference)
  )
)

assert_true(
  identical(
    rownames(epithelial_metadata),
    rownames(epithelial_reference)
  ),
  "Stage 2 metadata rownames are not aligned with epithelial reference."
)

assert_true(
  anyDuplicated(
    rownames(epithelial_metadata)
  ) == 0L,
  "Stage 2 metadata contains duplicated cell IDs."
)

cat(
  "Metadata alignment mode:",
  metadata_alignment_mode,
  "\n"
)

cat(
  "Stage 2 metadata dimensions:",
  nrow(epithelial_metadata),
  "rows x",
  ncol(epithelial_metadata),
  "columns\n"
)

cat(
  "[PASS] Stage 2 metadata deterministically aligned with epithelial reference\n\n"
)


# ==============================================================================
# 11. Reference cell-count audit
# ==============================================================================

observed_counts <- table(
  factor(
    epithelial_labels,
    levels = EPITHELIAL_IDENTITIES
  )
)

observed_counts <- as.integer(observed_counts)
names(observed_counts) <- EPITHELIAL_IDENTITIES

cell_count_table <- data.frame(
  identity = EPITHELIAL_IDENTITIES,
  observed_cells = observed_counts[EPITHELIAL_IDENTITIES],
  expected_cells = EXPECTED_COUNTS[EPITHELIAL_IDENTITIES],
  difference =
    observed_counts[EPITHELIAL_IDENTITIES] -
    EXPECTED_COUNTS[EPITHELIAL_IDENTITIES],
  stringsAsFactors = FALSE
)

cell_count_table$status <- ifelse(
  cell_count_table$difference == 0,
  "PASS",
  "FAIL"
)

write_tsv(
  cell_count_table,
  OUT_REFERENCE_COUNTS
)

print(cell_count_table)

assert_true(
  all(cell_count_table$status == "PASS"),
  "One or more epithelial identity cell counts do not match locked expectations."
)

assert_true(
  sum(observed_counts) == 93386L,
  "Epithelial reference total is not 93,386 cells."
)

cat("\n[PASS] All 7 epithelial cell counts match locked values\n\n")


# ==============================================================================
# 12. Reference numerical audit
# ==============================================================================

cat("Auditing epithelial reference numerical properties...\n")

reference_stats <- numeric_matrix_stats(
  epithelial_reference
)

reference_audit <- do.call(
  rbind,
  list(

    audit_row(
      "n_cells",
      nrow(epithelial_reference),
      93386,
      ifelse(
        nrow(epithelial_reference) == 93386L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "n_genes_before_cleanup",
      ncol(epithelial_reference),
      14431,
      ifelse(
        ncol(epithelial_reference) == 14431L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "n_identities",
      length(unique(epithelial_labels)),
      7,
      ifelse(
        length(unique(epithelial_labels)) == 7L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "duplicated_cell_ids",
      anyDuplicated(rownames(epithelial_reference)),
      0,
      ifelse(
        anyDuplicated(rownames(epithelial_reference)) == 0L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "duplicated_gene_names",
      anyDuplicated(colnames(epithelial_reference)),
      0,
      ifelse(
        anyDuplicated(colnames(epithelial_reference)) == 0L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "NA_entries",
      reference_stats$na_n,
      0,
      ifelse(
        reference_stats$na_n == 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "nonfinite_entries",
      reference_stats$nonfinite_n,
      0,
      ifelse(
        reference_stats$nonfinite_n == 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "negative_entries",
      reference_stats$negative_n,
      0,
      ifelse(
        reference_stats$negative_n == 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "minimum_value",
      reference_stats$min_value,
      ">=0",
      ifelse(
        reference_stats$min_value >= 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "maximum_value",
      reference_stats$max_value,
      "informational",
      "INFO"
    ),

    audit_row(
      "fractional_entries",
      reference_stats$fractional_n,
      0,
      ifelse(
        reference_stats$fractional_n == 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "metadata_rows",
      nrow(epithelial_metadata),
      93386,
      ifelse(
        nrow(epithelial_metadata) == 93386L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "metadata_alignment_mode",
      metadata_alignment_mode,
      "deterministic",
      "PASS"
    )
  )
)

write_tsv(
  reference_audit,
  OUT_REFERENCE_AUDIT
)

print(reference_audit)

assert_true(
  !any(reference_audit$status == "FAIL"),
  "Stage 2 epithelial reference audit contains FAIL."
)

cat("\n[PASS] Stage 2 epithelial raw-count reference audit\n\n")


# ==============================================================================
# 13. Build formal Stage 2 reference object
# ==============================================================================

stage2_reference_object <- list(

  reference = epithelial_reference,

  input.type = reference_obj$input.type,

  cell.type.labels = epithelial_cell_type_labels,

  cell.state.labels = epithelial_cell_state_labels,

  metadata = epithelial_metadata,

  provenance = list(

    project = "IBD_EcoTyper",

    stage = "Stage2_epithelial_BayesPrism",

    source_reference = REFERENCE_RDS,

    source_reference_cells = 323144L,

    source_reference_genes = 14431L,

    selected_epithelial_cells = 93386L,

    selected_epithelial_identities = EPITHELIAL_IDENTITIES,

    construction =
      "Direct subset of CLOSED unified raw-count reference; no reselection or downsampling",

    bayesprism_hierarchy =
      "cell.type = cell.state; seven epithelial components 1:1",

    source_metadata_class =
      paste(class(source_metadata), collapse = ","),

    metadata_alignment_mode =
      metadata_alignment_mode,

    created_at = as.character(Sys.time())
  )
)

saveRDS(
  stage2_reference_object,
  OUT_REFERENCE_RDS,
  compress = TRUE
)

cat(
  "[PASS] Saved Stage 2 epithelial reference:\n",
  OUT_REFERENCE_RDS,
  "\n\n"
)


# ==============================================================================
# 14. Release full unified reference from memory
# ==============================================================================

rm(
  reference_raw,
  reference_obj,
  source_metadata,
  source_type_labels
)

gc(verbose = FALSE)


# ==============================================================================
# 15. Load Stage 1 Broad Epithelial Z pseudo-mixture
# ==============================================================================

cat("Loading Stage 1 Broad Epithelial Z...\n")

epithelial_z_obj <- readRDS(
  EPITHELIAL_Z_RDS
)

epithelial_pseudomixture <- unwrap_matrix(
  epithelial_z_obj,
  "Stage 1 Epithelial Z"
)

rm(epithelial_z_obj)
gc(verbose = FALSE)

cat(
  "Pseudo-mixture class:",
  paste(class(epithelial_pseudomixture), collapse = ", "),
  "\n"
)

cat(
  "Pseudo-mixture dimensions:",
  nrow(epithelial_pseudomixture),
  "samples x",
  ncol(epithelial_pseudomixture),
  "genes\n"
)


# ==============================================================================
# 16. Structural audit of Stage 1 Epithelial Z
# ==============================================================================

assert_true(
  nrow(epithelial_pseudomixture) == 2490L,
  paste0(
    "Unexpected Stage 1 Epithelial Z sample count: ",
    nrow(epithelial_pseudomixture),
    " != 2490"
  )
)

assert_true(
  ncol(epithelial_pseudomixture) == 13668L,
  paste0(
    "Unexpected Stage 1 Epithelial Z gene count: ",
    ncol(epithelial_pseudomixture),
    " != 13668"
  )
)

assert_true(
  !is.null(rownames(epithelial_pseudomixture)),
  "Epithelial pseudo-mixture has no sample rownames."
)

assert_true(
  !is.null(colnames(epithelial_pseudomixture)),
  "Epithelial pseudo-mixture has no gene colnames."
)

assert_true(
  anyDuplicated(rownames(epithelial_pseudomixture)) == 0L,
  "Duplicated sample IDs found in Epithelial pseudo-mixture."
)

assert_true(
  anyDuplicated(colnames(epithelial_pseudomixture)) == 0L,
  "Duplicated gene names found in Epithelial pseudo-mixture."
)


# ==============================================================================
# 17. Sample-order audit against Stage 1 Broad14 fraction table
# ==============================================================================

stage1_sample_ids <- read_stage1_sample_ids(
  BROAD14_FRACTION_INITIAL,
  expected_n = 2490L
)

sample_order_status <- "NOT_CHECKED"

if (is.null(stage1_sample_ids)) {

  cat(
    "[WARN] Could not independently recover sample IDs from ",
    "Broad14 initial fraction table.\n"
  )

  sample_order_status <- "WARN"

} else {

  assert_true(
    setequal(
      rownames(epithelial_pseudomixture),
      stage1_sample_ids
    ),
    paste0(
      "Pseudo-mixture and Broad14 fraction table ",
      "do not contain the same sample IDs."
    )
  )

  assert_true(
    identical(
      rownames(epithelial_pseudomixture),
      stage1_sample_ids
    ),
    paste0(
      "Pseudo-mixture sample IDs are identical as a set, ",
      "but the row order differs from Stage 1 Broad14."
    )
  )

  sample_order_status <- "PASS"

  cat(
    "[PASS] Pseudo-mixture sample IDs and order ",
    "match Stage 1 Broad14\n"
  )
}

cat("\n")


# ==============================================================================
# 18. Numerical audit of epithelial pseudo-mixture
# ==============================================================================

cat("Auditing Stage 1 Epithelial Z numerical properties...\n")

mixture_stats <- numeric_matrix_stats(
  epithelial_pseudomixture,
  chunk_size = 50L
)

mixture_row_sums <- Matrix::rowSums(
  epithelial_pseudomixture
)

zero_or_negative_samples <- sum(
  mixture_row_sums <= 0 |
    !is.finite(mixture_row_sums)
)

mixture_audit <- do.call(
  rbind,
  list(

    audit_row(
      "n_samples",
      nrow(epithelial_pseudomixture),
      2490,
      ifelse(
        nrow(epithelial_pseudomixture) == 2490L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "n_genes",
      ncol(epithelial_pseudomixture),
      13668,
      ifelse(
        ncol(epithelial_pseudomixture) == 13668L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "duplicated_sample_ids",
      anyDuplicated(rownames(epithelial_pseudomixture)),
      0,
      ifelse(
        anyDuplicated(rownames(epithelial_pseudomixture)) == 0L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "duplicated_gene_names",
      anyDuplicated(colnames(epithelial_pseudomixture)),
      0,
      ifelse(
        anyDuplicated(colnames(epithelial_pseudomixture)) == 0L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "NA_entries",
      mixture_stats$na_n,
      0,
      ifelse(
        mixture_stats$na_n == 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "nonfinite_entries",
      mixture_stats$nonfinite_n,
      0,
      ifelse(
        mixture_stats$nonfinite_n == 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "negative_entries",
      mixture_stats$negative_n,
      0,
      ifelse(
        mixture_stats$negative_n == 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "fractional_entries",
      mixture_stats$fractional_n,
      ">0 expected",
      ifelse(
        mixture_stats$fractional_n > 0,
        "PASS",
        "WARN"
      )
    ),

    audit_row(
      "samples_with_nonpositive_total_mass",
      zero_or_negative_samples,
      0,
      ifelse(
        zero_or_negative_samples == 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "minimum_value",
      mixture_stats$min_value,
      ">=0",
      ifelse(
        mixture_stats$min_value >= 0,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "maximum_value",
      mixture_stats$max_value,
      "informational",
      "INFO"
    ),

    audit_row(
      "sample_order_vs_stage1",
      sample_order_status,
      "PASS",
      sample_order_status
    )
  )
)

write_tsv(
  mixture_audit,
  OUT_MIXTURE_AUDIT
)

print(mixture_audit)

assert_true(
  !any(mixture_audit$status == "FAIL"),
  "Stage 2 epithelial pseudo-mixture audit contains FAIL."
)

cat(
  "\n[PASS] Stage 1 Epithelial Z structural/numerical audit\n\n"
)


# ==============================================================================
# 19. Gene-space overlap audit BEFORE Stage 2 cleanup
# ==============================================================================

cat("Auditing pre-cleanup gene-space overlap...\n")

reference_genes <- colnames(
  epithelial_reference
)

mixture_genes <- colnames(
  epithelial_pseudomixture
)

gene_union <- union(
  reference_genes,
  mixture_genes
)

gene_overlap_table <- data.frame(
  gene = gene_union,
  in_epithelial_reference =
    gene_union %in% reference_genes,
  in_stage1_epithelial_Z =
    gene_union %in% mixture_genes,
  stringsAsFactors = FALSE
)

gene_overlap_table$category <- ifelse(
  gene_overlap_table$in_epithelial_reference &
    gene_overlap_table$in_stage1_epithelial_Z,
  "overlap",
  ifelse(
    gene_overlap_table$in_epithelial_reference,
    "reference_only",
    "stage1_Z_only"
  )
)

write_tsv(
  gene_overlap_table,
  OUT_GENE_OVERLAP
)

n_ref_genes <- length(reference_genes)
n_mix_genes <- length(mixture_genes)
n_overlap <- length(
  intersect(reference_genes, mixture_genes)
)
n_ref_only <- length(
  setdiff(reference_genes, mixture_genes)
)
n_mix_only <- length(
  setdiff(mixture_genes, reference_genes)
)

gene_summary <- data.frame(
  metric = c(
    "epithelial_reference_genes_before_cleanup",
    "stage1_epithelial_Z_genes",
    "overlap_genes",
    "reference_only_genes",
    "stage1_Z_only_genes"
  ),
  value = c(
    n_ref_genes,
    n_mix_genes,
    n_overlap,
    n_ref_only,
    n_mix_only
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  gene_summary,
  OUT_GENE_SUMMARY
)

print(gene_summary)

assert_true(
  n_ref_genes == 14431L,
  "Unexpected pre-cleanup epithelial reference gene count."
)

assert_true(
  n_mix_genes == 13668L,
  "Unexpected Stage 1 Epithelial Z gene count."
)

assert_true(
  n_mix_only == 0L,
  paste0(
    n_mix_only,
    " Stage 1 Epithelial Z genes are absent from ",
    "the original epithelial raw-count reference."
  )
)

cat(
  "\nPre-cleanup overlap:",
  n_overlap,
  "genes\n"
)

cat(
  "Reference-only genes:",
  n_ref_only,
  "\n"
)

cat(
  "Stage1-Z-only genes:",
  n_mix_only,
  "\n\n"
)

cat(
  "[PASS] Stage 1 Epithelial Z gene universe is contained ",
  "within the raw epithelial reference universe\n\n"
)


# ==============================================================================
# 20. Save frozen Stage 2 pseudo-mixture
# ==============================================================================

saveRDS(
  epithelial_pseudomixture,
  OUT_MIXTURE_RDS,
  compress = TRUE
)

cat(
  "[PASS] Saved Stage 2 epithelial pseudo-mixture:\n",
  OUT_MIXTURE_RDS,
  "\n\n"
)


# ==============================================================================
# 21. Input/output manifest
# ==============================================================================

manifest <- data.frame(

  role = c(
    "SOURCE_REFERENCE",
    "SOURCE_STAGE1_EPITHELIAL_Z",
    "SOURCE_STAGE1_FRACTION_TABLE",
    "STAGE2_REFERENCE",
    "STAGE2_PSEUDOMIXTURE",
    "REFERENCE_CELL_COUNTS",
    "REFERENCE_AUDIT",
    "MIXTURE_AUDIT",
    "GENE_OVERLAP",
    "GENE_OVERLAP_SUMMARY",
    "LOG",
    "SESSION_INFO"
  ),

  path = c(
    REFERENCE_RDS,
    EPITHELIAL_Z_RDS,
    BROAD14_FRACTION_INITIAL,
    OUT_REFERENCE_RDS,
    OUT_MIXTURE_RDS,
    OUT_REFERENCE_COUNTS,
    OUT_REFERENCE_AUDIT,
    OUT_MIXTURE_AUDIT,
    OUT_GENE_OVERLAP,
    OUT_GENE_SUMMARY,
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
# 22. sessionInfo()
# ==============================================================================

capture.output(
  sessionInfo(),
  file = OUT_SESSIONINFO
)


# ==============================================================================
# 23. Final summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("S2-A SUMMARY\n")
cat("============================================================\n")

cat(
  "Stage 2 raw epithelial reference :",
  nrow(epithelial_reference),
  "cells x",
  ncol(epithelial_reference),
  "genes\n"
)

cat(
  "Stage 2 pseudo-mixture           :",
  nrow(epithelial_pseudomixture),
  "samples x",
  ncol(epithelial_pseudomixture),
  "genes\n"
)

cat(
  "Epithelial identities            :",
  length(unique(epithelial_labels)),
  "\n"
)

cat(
  "Pre-cleanup gene overlap         :",
  n_overlap,
  "\n"
)

cat(
  "Reference-only genes             :",
  n_ref_only,
  "\n"
)

cat(
  "Stage1-Z-only genes              :",
  n_mix_only,
  "\n"
)

cat(
  "Fractional pseudo-mixture values :",
  mixture_stats$fractional_n,
  "\n"
)

cat(
  "Sample-order audit               :",
  sample_order_status,
  "\n"
)

cat("\n")

cat(
  "IMPORTANT:\n",
  "No Stage 2 gene cleanup or BayesPrism model has been run.\n",
  "Do NOT proceed to run.prism() before reviewing these audit outputs.\n"
)

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\n============================================================\n")
cat("S2-A INPUT CONSTRUCTION COMPLETED\n")
cat("END:", format(Sys.time()), "\n")
cat("============================================================\n")