#!/usr/bin/env Rscript

## ============================================================================
## 04_extract_GSE282122_raw_counts.R
##
## Purpose:
##   Extract the final selected GSE282122 cells from the authoritative raw-count
##   H5AD source, preserving the exact order defined by
##   selected_cell_manifest.parquet.
##
## Input:
##   1. selected_cell_manifest.parquet
##   2. TAURUS_raw_counts_annotated_final.h5ad
##
## Output:
##   - GSE282122_selected_raw_counts.h5
##   - GSE282122_selected_cell_manifest.parquet
##   - GSE282122_gene_metadata.parquet
##   - GSE282122_raw_count_extraction_audit.tsv
##   - GSE282122_raw_count_extraction_sessionInfo.txt
##
## Key safeguards:
##   - Exact source_cell_key -> raw source cell mapping required
##   - No silent cell loss
##   - /X must pass raw integer-count audit
##   - Manifest order is preserved exactly
##   - Output matrix is written as sparse HDF5
##   - Random post-write concordance test required
##
## Date:
##   2026-09-13
## ============================================================================


## ============================================================================
## 0. Configuration
## ============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

manifest_path <- paste0(
  "/home/mazekai/IBD_EcoTyper/",
  "03_reference/combined_reference/output/",
  "05_final_cell_selection/selected_cell_manifest.parquet"
)

source_h5ad <- paste0(
  "/home/mazekai/enteric_glia/data_scRNA/GSE282122/",
  "TAURUS_raw_counts_annotated_final.h5ad"
)

output_dir <- paste0(
  "/home/mazekai/IBD_EcoTyper/",
  "03_reference/combined_reference/output/",
  "06_raw_count_extraction/GSE282122"
)

expected_n_selected <- 103538L

## Set environment variable IBD_REF_OVERWRITE=1 only when intentional.
overwrite <- identical(
  Sys.getenv("IBD_REF_OVERWRITE"),
  "1"
)


## ============================================================================
## 1. Package checks
## ============================================================================

required_pkgs <- c(
  "arrow",
  "zellkonverter",
  "SummarizedExperiment",
  "SingleCellExperiment",
  "HDF5Array",
  "rhdf5"
)

pkg_ok <- vapply(
  required_pkgs,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

if (!all(pkg_ok)) {
  stop(
    "Missing required R packages: ",
    paste(required_pkgs[!pkg_ok], collapse = ", ")
  )
}

cat("\n========================================\n")
cat("GSE282122 RAW-COUNT EXTRACTION\n")
cat("========================================\n")


## ============================================================================
## 2. Input checks
## ============================================================================

stopifnot(file.exists(manifest_path))
stopifnot(file.exists(source_h5ad))

cat("\nManifest:\n", manifest_path, "\n", sep = "")
cat("\nRaw source:\n", source_h5ad, "\n", sep = "")


## ============================================================================
## 3. Load and validate final cell manifest
## ============================================================================

cat("\n[1] Loading selected-cell manifest...\n")

manifest <- arrow::read_parquet(manifest_path)

required_manifest_columns <- c(
  "global_cell_id",
  "dataset",
  "cell_id",
  "source_cell_key",
  "source_bundle",
  "donor_id",
  "sample_id",
  "reference_identity",
  "final_reference_action"
)

missing_manifest_columns <- setdiff(
  required_manifest_columns,
  colnames(manifest)
)

if (length(missing_manifest_columns) > 0L) {
  stop(
    "Manifest is missing required columns: ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}

## Global manifest key integrity
stopifnot(
  !base::anyDuplicated(manifest$global_cell_id),
  !base::anyDuplicated(manifest$cell_id),
  !base::anyDuplicated(manifest$source_cell_key),
  all(manifest$cell_id == manifest$source_cell_key),
  all(
    manifest$global_cell_id ==
      paste0(manifest$dataset, "::", manifest$source_cell_key)
  )
)

manifest_gse <- manifest[
  manifest$dataset == "GSE282122",
  ,
  drop = FALSE
]

stopifnot(
  nrow(manifest_gse) == expected_n_selected,
  !base::anyDuplicated(manifest_gse$source_cell_key),
  !base::anyDuplicated(manifest_gse$global_cell_id)
)

cat(
  "Selected GSE282122 cells:",
  nrow(manifest_gse),
  "\n"
)

cat("Manifest integrity: PASS\n")


## ============================================================================
## 4. Audit H5AD storage architecture
## ============================================================================

cat("\n[2] Auditing H5AD structure...\n")

h5_structure <- rhdf5::h5ls(
  source_h5ad,
  recursive = TRUE
)

has_x <- any(
  h5_structure$group == "/" &
    h5_structure$name == "X"
)

has_raw <- any(
  h5_structure$group == "/" &
    h5_structure$name == "raw"
)

has_layers <- any(
  h5_structure$group == "/" &
    h5_structure$name == "layers"
)

if (!has_x) {
  stop("H5AD does not contain /X")
}

cat("Contains /X     :", has_x, "\n")
cat("Contains /raw   :", has_raw, "\n")
cat("Contains /layers:", has_layers, "\n")


## ============================================================================
## 5. Load H5AD as HDF5-backed SingleCellExperiment
## ============================================================================

cat("\n[3] Reading H5AD as HDF5-backed SCE...\n")

sce <- zellkonverter::readH5AD(
  source_h5ad,
  reader = "R",
  use_hdf5 = TRUE,
  skip_assays = FALSE,
  raw = FALSE,
  layers = FALSE,
  verbose = TRUE
)

cat(
  "Source dimensions:",
  paste(dim(sce), collapse = " x "),
  "\n"
)

cat(
  "Assays:",
  paste(
    SummarizedExperiment::assayNames(sce),
    collapse = ", "
  ),
  "\n"
)

if (!"X" %in% SummarizedExperiment::assayNames(sce)) {
  stop("Expected assay 'X' was not found.")
}

stopifnot(
  nrow(sce) == 33075L,
  ncol(sce) == 987743L,
  !base::anyDuplicated(colnames(sce))
)

x <- SummarizedExperiment::assay(
  sce,
  "X"
)


## ============================================================================
## 6. Exact selected-cell mapping audit
## ============================================================================

cat("\n[4] Auditing selected-cell mapping...\n")

manifest_keys <- manifest_gse$source_cell_key
source_keys   <- colnames(sce)

idx <- base::match(
  manifest_keys,
  source_keys
)

n_matched <- sum(!is.na(idx))
n_missing <- sum(is.na(idx))

stopifnot(
  length(manifest_keys) == expected_n_selected,
  n_matched == expected_n_selected,
  n_missing == 0L,
  !base::anyDuplicated(idx),
  identical(
    source_keys[idx],
    manifest_keys
  )
)

cat("Manifest selected cells :", length(manifest_keys), "\n")
cat("Raw source cells        :", length(source_keys), "\n")
cat("Matched                 :", n_matched, "\n")
cat("Missing                 :", n_missing, "\n")
cat("Match rate              : 100.000000%\n")
cat("Exact 1:1 mapping       : PASS\n")


## ============================================================================
## 7. Direct /X/data raw-count audit
## ============================================================================

cat("\n[5] Directly auditing /X/data...\n")

x_data_sample <- rhdf5::h5read(
  source_h5ad,
  "/X/data",
  index = list(seq_len(100000L))
)

direct_negative <- sum(x_data_sample < 0)

direct_integer_fraction <- mean(
  abs(x_data_sample - round(x_data_sample)) < 1e-8
)

cat("Direct values sampled :", length(x_data_sample), "\n")
cat("Minimum               :", min(x_data_sample), "\n")
cat("Maximum               :", max(x_data_sample), "\n")
cat("Negative values       :", direct_negative, "\n")
cat(
  "Integer fraction      :",
  sprintf("%.10f", direct_integer_fraction),
  "\n"
)

stopifnot(
  direct_negative == 0L,
  abs(direct_integer_fraction - 1) < 1e-12
)

cat("Direct /X/data audit  : PASS\n")


## ============================================================================
## 8. DelayedMatrix raw-count audit
## ============================================================================

cat("\n[6] Auditing sampled cells through DelayedMatrix...\n")

gene_idx <- unique(
  round(
    seq(
      1,
      nrow(x),
      length.out = 500L
    )
  )
)

selected_positions <- unique(
  round(
    seq(
      1,
      length(idx),
      length.out = 200L
    )
  )
)

cell_idx <- idx[selected_positions]

test_block <- as.matrix(
  x[
    gene_idx,
    cell_idx,
    drop = FALSE
  ]
)

vals <- as.numeric(test_block)
vals <- vals[is.finite(vals)]

nz <- vals[vals != 0]

matrix_negative <- sum(vals < 0)

matrix_integer_fraction <- mean(
  abs(nz - round(nz)) < 1e-8
)

cat("Sampled entries       :", length(vals), "\n")
cat("Non-zero entries      :", length(nz), "\n")
cat("Minimum               :", min(vals), "\n")
cat("Maximum               :", max(vals), "\n")
cat("Negative values       :", matrix_negative, "\n")
cat(
  "Integer fraction      :",
  sprintf("%.10f", matrix_integer_fraction),
  "\n"
)

stopifnot(
  matrix_negative == 0L,
  length(nz) > 0L,
  abs(matrix_integer_fraction - 1) < 1e-12
)

cat("DelayedMatrix audit   : PASS\n")


## ============================================================================
## 9. Gene identifier audit
## ============================================================================

cat("\n[7] Auditing gene identifiers...\n")

gene_ids <- rownames(sce)

if (is.null(gene_ids)) {
  stop("Source object has no gene rownames.")
}

n_gene_na <- sum(is.na(gene_ids))
n_gene_blank <- sum(gene_ids == "")
n_gene_duplicate <- sum(duplicated(gene_ids))

cat("Genes                  :", length(gene_ids), "\n")
cat("Unique genes           :", length(unique(gene_ids)), "\n")
cat("NA gene IDs            :", n_gene_na, "\n")
cat("Blank gene IDs         :", n_gene_blank, "\n")
cat("Duplicated gene IDs    :", n_gene_duplicate, "\n")

stopifnot(
  length(gene_ids) == nrow(x),
  n_gene_na == 0L,
  n_gene_blank == 0L,
  n_gene_duplicate == 0L
)

cat("Gene identifier audit : PASS\n")


## ============================================================================
## 10. Build selected raw-count matrix
## ============================================================================

cat("\n[8] Constructing selected delayed matrix...\n")

counts_selected <- x[
  ,
  idx,
  drop = FALSE
]

rownames(counts_selected) <- gene_ids
colnames(counts_selected) <- manifest_gse$global_cell_id

stopifnot(
  nrow(counts_selected) == 33075L,
  ncol(counts_selected) == expected_n_selected,
  identical(
    colnames(counts_selected),
    manifest_gse$global_cell_id
  ),
  !base::anyDuplicated(colnames(counts_selected))
)

cat(
  "Selected dimensions:",
  paste(dim(counts_selected), collapse = " x "),
  "\n"
)

cat("Selected matrix construction: PASS\n")


## ============================================================================
## 11. Prepare output
## ============================================================================

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(dir.exists(output_dir))

counts_file <- file.path(
  output_dir,
  "GSE282122_selected_raw_counts.h5"
)

manifest_out <- file.path(
  output_dir,
  "GSE282122_selected_cell_manifest.parquet"
)

gene_metadata_out <- file.path(
  output_dir,
  "GSE282122_gene_metadata.parquet"
)

audit_out <- file.path(
  output_dir,
  "GSE282122_raw_count_extraction_audit.tsv"
)

session_out <- file.path(
  output_dir,
  "GSE282122_raw_count_extraction_sessionInfo.txt"
)

if (file.exists(counts_file)) {

  if (!overwrite) {
    stop(
      "Output already exists: ",
      counts_file,
      "\nRefusing to overwrite. ",
      "Set IBD_REF_OVERWRITE=1 only for an intentional rerun."
    )
  }

  cat(
    "IBD_REF_OVERWRITE=1: removing existing count matrix.\n"
  )

  unlink(counts_file)
}


## ============================================================================
## 12. Save metadata
## ============================================================================

cat("\n[9] Writing metadata...\n")

arrow::write_parquet(
  manifest_gse,
  manifest_out
)

gene_metadata <- as.data.frame(
  SummarizedExperiment::rowData(sce)
)

gene_metadata$source_gene_id <- gene_ids

gene_metadata <- gene_metadata[
  c(
    "source_gene_id",
    setdiff(
      colnames(gene_metadata),
      "source_gene_id"
    )
  )
]

arrow::write_parquet(
  gene_metadata,
  gene_metadata_out
)

cat("Metadata written: PASS\n")


## ============================================================================
## 13. Write sparse raw-count matrix
## ============================================================================

cat("\n[10] Writing sparse HDF5 raw-count matrix...\n")

counts_h5 <- HDF5Array::writeTENxMatrix(
  counts_selected,
  filepath = counts_file,
  group = "raw_counts",
  level = 6,
  verbose = TRUE
)

stopifnot(
  identical(
    dim(counts_h5),
    c(33075L, expected_n_selected)
  )
)

cat("Sparse HDF5 write: PASS\n")


## ============================================================================
## 14. Post-write concordance audit
## ============================================================================

cat("\n[11] Running post-write concordance audit...\n")

rownames_identical <- identical(
  rownames(counts_h5),
  rownames(counts_selected)
)

colnames_identical <- identical(
  colnames(counts_h5),
  manifest_gse$global_cell_id
)

stopifnot(
  rownames_identical,
  colnames_identical
)

set.seed(20260913L)

test_genes <- sample(
  seq_len(nrow(counts_selected)),
  100L
)

test_cells <- sample(
  seq_len(ncol(counts_selected)),
  100L
)

before <- as.matrix(
  counts_selected[
    test_genes,
    test_cells,
    drop = FALSE
  ]
)

after <- as.matrix(
  counts_h5[
    test_genes,
    test_cells,
    drop = FALSE
  ]
)

value_concordance <- identical(
  as.numeric(before),
  as.numeric(after)
)

stopifnot(value_concordance)

cat("Row names identical    :", rownames_identical, "\n")
cat("Column names identical :", colnames_identical, "\n")
cat("Random concordance     :", value_concordance, "\n")
cat("Post-write audit       : PASS\n")


## ============================================================================
## 15. Write audit summary and session information
## ============================================================================

audit_summary <- data.frame(
  dataset = "GSE282122",
  source_file = source_h5ad,
  source_n_genes = nrow(sce),
  source_n_cells = ncol(sce),
  manifest_n_selected = nrow(manifest_gse),
  n_matched = n_matched,
  n_missing = n_missing,
  match_rate = n_matched / nrow(manifest_gse),
  x_direct_negative_values = direct_negative,
  x_direct_integer_fraction = direct_integer_fraction,
  x_matrix_negative_values = matrix_negative,
  x_matrix_integer_fraction = matrix_integer_fraction,
  n_gene_duplicates = n_gene_duplicate,
  output_n_genes = nrow(counts_h5),
  output_n_cells = ncol(counts_h5),
  rownames_identical = rownames_identical,
  colnames_identical = colnames_identical,
  random_value_concordance = value_concordance,
  stringsAsFactors = FALSE
)

utils::write.table(
  audit_summary,
  audit_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = session_out
)


## ============================================================================
## 16. Final summary
## ============================================================================

cat("\n========================================\n")
cat("GSE282122 EXTRACTION COMPLETE\n")
cat("========================================\n")

cat("Selected cells :", ncol(counts_h5), "\n")
cat("Genes          :", nrow(counts_h5), "\n")
cat("Mapping        : PASS\n")
cat("Raw-count audit: PASS\n")
cat("Gene audit     : PASS\n")
cat("Sparse write   : PASS\n")
cat("Concordance    : PASS\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nRaw-count matrix:\n")
cat(counts_file, "\n")

cat("\nALL CHECKS PASSED\n")
