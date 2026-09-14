#!/usr/bin/env Rscript

## ============================================================================
## 05_extract_SCP1884_raw_counts.R
##
## Purpose:
##   Extract raw counts for the final selected SCP1884 reference cells from
##   six source bundles:
##
##     CO_EPI, CO_IMM, CO_STR
##     TI_EPI, TI_IMM, TI_STR
##
##   Colon and terminal-ileum source matrices use different gene universes:
##
##     Colon: 28,663 genes
##     TI:    28,923 genes
##
##   The canonical SCP1884 reference gene universe is therefore defined as
##   the ordered intersection of Ensembl IDs:
##
##     27,830 common genes
##
##   Ordering is inherited from the Colon feature file.
##
##   Each bundle is:
##     1. audited independently
##     2. matched exactly to selected_cell_manifest.parquet
##     3. harmonized to the canonical gene universe
##     4. written as sparse HDF5
##     5. checked for post-write concordance
##
##   The six extracted bundles are finally concatenated and reordered to
##   exactly match the authoritative SCP1884 manifest cell order.
##
## Inputs:
##   - selected_cell_manifest.parquet
##   - six SCP1884 raw.mtx matrices
##   - six barcode files
##   - six feature files
##
## Outputs:
##   - SCP1884_selected_raw_counts.h5
##   - SCP1884_selected_cell_manifest.parquet
##   - SCP1884_gene_metadata.parquet
##   - SCP1884_bundle_extraction_audit.tsv
##   - SCP1884_gene_harmonization_audit.tsv
##   - SCP1884_raw_count_extraction_sessionInfo.txt
##   - six bundle-level sparse HDF5 matrices
##
## Safety:
##   Existing outputs are NOT overwritten by default.
##   To intentionally rerun from scratch:
##
##     IBD_REF_OVERWRITE=1 Rscript --vanilla \
##       scripts/R/reference/05_extract_SCP1884_raw_counts.R
##
## Date:
##   2026-09-14
## ============================================================================


## ============================================================================
## 0. Global configuration
## ============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

overwrite <- identical(
  Sys.getenv("IBD_REF_OVERWRITE"),
  "1"
)

project_root <- "/home/mazekai/IBD_EcoTyper"

expression_root <- paste0(
  "/home/mazekai/enteric_glia/",
  "data_scRNA/SCP1884/expression"
)

manifest_path <- file.path(
  project_root,
  "03_reference/combined_reference/output/",
  "05_final_cell_selection/",
  "selected_cell_manifest.parquet"
)

output_dir <- file.path(
  project_root,
  "03_reference/combined_reference/output/",
  "06_raw_count_extraction/SCP1884"
)

bundle_output_dir <- file.path(
  output_dir,
  "bundles"
)

expected_total_selected <- 155384L
expected_common_genes <- 27830L


## ============================================================================
## 1. Package checks
## ============================================================================

required_pkgs <- c(
  "Matrix",
  "arrow",
  "DelayedArray",
  "HDF5Array"
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
    paste(
      required_pkgs[!pkg_ok],
      collapse = ", "
    )
  )
}


## ============================================================================
## 2. Helper functions
## ============================================================================

read_features <- function(path) {

  x <- utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c(
      "ensembl_id",
      "gene_symbol"
    )
  )

  x
}


read_barcodes <- function(path) {

  utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )[[1]]
}


## Full audit of sparse non-zero values without creating several
## multi-GB temporary vectors at once.
audit_sparse_values <- function(
    values,
    chunk_size = 5000000L
) {

  n <- length(values)

  if (n == 0L) {
    stop("Sparse matrix contains no non-zero values.")
  }

  n_negative <- 0
  n_noninteger <- 0

  min_value <- Inf
  max_value <- -Inf

  starts <- seq.int(
    from = 1L,
    to = n,
    by = chunk_size
  )

  for (start in starts) {

    end <- min(
      start + chunk_size - 1L,
      n
    )

    v <- values[start:end]

    n_negative <- n_negative +
      sum(v < 0)

    n_noninteger <- n_noninteger +
      sum(
        abs(v - round(v)) > 1e-8
      )

    min_value <- min(
      min_value,
      min(v)
    )

    max_value <- max(
      max_value,
      max(v)
    )
  }

  list(
    n_values = n,
    n_negative = n_negative,
    n_noninteger = n_noninteger,
    integer_fraction =
      1 - (n_noninteger / n),
    min_value = min_value,
    max_value = max_value
  )
}


## ============================================================================
## 3. Input checks
## ============================================================================

stopifnot(
  dir.exists(expression_root),
  file.exists(manifest_path)
)

cat("\n")
cat("========================================\n")
cat("SCP1884 RAW-COUNT EXTRACTION\n")
cat("========================================\n")

cat("\nExpression root:\n")
cat(expression_root, "\n")

cat("\nManifest:\n")
cat(manifest_path, "\n")


## ============================================================================
## 4. Bundle configuration
## ============================================================================

bundle_config <- data.frame(

  source_bundle = c(
    "CO_EPI",
    "CO_IMM",
    "CO_STR",
    "TI_EPI",
    "TI_IMM",
    "TI_STR"
  ),

  subdir = c(
    "Colon/Epithelial",
    "Colon/Immune",
    "Colon/Stromal",
    "TerminalIleum/Epithelial",
    "TerminalIleum/Immune",
    "TerminalIleum/Stromal"
  ),

  expected_selected = c(
    19501L,
    35699L,
    11500L,
    27312L,
    48271L,
    13101L
  ),

  expected_source_genes = c(
    28663L,
    28663L,
    28663L,
    28923L,
    28923L,
    28923L
  ),

  expected_source_cells = c(
    97788L,
    152509L,
    39433L,
    154136L,
    201072L,
    75695L
  ),

  stringsAsFactors = FALSE
)

stopifnot(
  sum(bundle_config$expected_selected) ==
    expected_total_selected
)


## ============================================================================
## 5. Manifest audit
## ============================================================================

cat("\n[1] Loading authoritative manifest...\n")

manifest <- as.data.frame(
  arrow::read_parquet(
    manifest_path
  )
)

required_manifest_columns <- c(
  "global_cell_id",
  "dataset",
  "cell_id",
  "source_cell_key",
  "source_bundle",
  "donor_id",
  "sample_id",
  "reference_identity"
)

missing_manifest_columns <- setdiff(
  required_manifest_columns,
  colnames(manifest)
)

if (length(missing_manifest_columns) > 0L) {
  stop(
    "Manifest missing columns: ",
    paste(
      missing_manifest_columns,
      collapse = ", "
    )
  )
}

manifest_scp1884 <- manifest[
  manifest$dataset == "SCP1884",
  ,
  drop = FALSE
]

stopifnot(
  nrow(manifest_scp1884) ==
    expected_total_selected,

  !base::anyDuplicated(
    manifest_scp1884$global_cell_id
  ),

  !base::anyDuplicated(
    manifest_scp1884$source_cell_key
  ),

  all(
    manifest_scp1884$global_cell_id ==
      paste0(
        "SCP1884::",
        manifest_scp1884$source_cell_key
      )
  )
)

manifest_bundle_counts <- table(
  manifest_scp1884$source_bundle
)

observed_bundle_counts <- as.integer(
  manifest_bundle_counts[
    bundle_config$source_bundle
  ]
)

stopifnot(
  identical(
    observed_bundle_counts,
    bundle_config$expected_selected
  )
)

cat(
  "Selected SCP1884 cells:",
  nrow(manifest_scp1884),
  "\n"
)

cat("Manifest integrity: PASS\n")


## ============================================================================
## 6. Source file audit
## ============================================================================

cat("\n[2] Auditing source files...\n")

for (i in seq_len(nrow(bundle_config))) {

  bundle <- bundle_config$source_bundle[i]
  subdir <- bundle_config$subdir[i]

  raw_file <- file.path(
    expression_root,
    subdir,
    "raw.mtx"
  )

  barcode_file <- file.path(
    expression_root,
    subdir,
    paste0(
      bundle,
      ".barcodes.tsv"
    )
  )

  feature_file <- file.path(
    expression_root,
    subdir,
    paste0(
      bundle,
      ".features.tsv"
    )
  )

  stopifnot(
    file.exists(raw_file),
    file.exists(barcode_file),
    file.exists(feature_file)
  )
}

cat("All source files present: PASS\n")


## ============================================================================
## 7. Establish canonical gene universe
## ============================================================================

cat("\n[3] Establishing canonical SCP1884 gene universe...\n")

colon_features <- read_features(
  file.path(
    expression_root,
    "Colon/Stromal/",
    "CO_STR.features.tsv"
  )
)

ti_features <- read_features(
  file.path(
    expression_root,
    "TerminalIleum/Stromal/",
    "TI_STR.features.tsv"
  )
)

stopifnot(
  nrow(colon_features) == 28663L,
  nrow(ti_features) == 28923L,

  !base::anyDuplicated(
    colon_features$ensembl_id
  ),

  !base::anyDuplicated(
    ti_features$ensembl_id
  ),

  !base::anyDuplicated(
    colon_features$gene_symbol
  ),

  !base::anyDuplicated(
    ti_features$gene_symbol
  )
)

common_gene_order <-
  colon_features$ensembl_id[
    colon_features$ensembl_id %in%
      ti_features$ensembl_id
  ]

colon_only <- setdiff(
  colon_features$ensembl_id,
  ti_features$ensembl_id
)

ti_only <- setdiff(
  ti_features$ensembl_id,
  colon_features$ensembl_id
)

stopifnot(
  length(common_gene_order) ==
    expected_common_genes,

  length(colon_only) == 833L,

  length(ti_only) == 1093L,

  !base::anyDuplicated(
    common_gene_order
  )
)

colon_common <- colon_features[
  base::match(
    common_gene_order,
    colon_features$ensembl_id
  ),
  ,
  drop = FALSE
]

ti_common <- ti_features[
  base::match(
    common_gene_order,
    ti_features$ensembl_id
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colon_common$ensembl_id,
    ti_common$ensembl_id
  ),

  identical(
    colon_common$gene_symbol,
    ti_common$gene_symbol
  )
)

common_gene_metadata <- colon_common

stopifnot(
  identical(
    common_gene_metadata$ensembl_id,
    common_gene_order
  )
)

cat("Colon genes      : 28663\n")
cat("TI genes         : 28923\n")
cat("Common genes     : 27830\n")
cat("Colon-only genes : 833\n")
cat("TI-only genes    : 1093\n")
cat("Symbol conflicts : 0\n")
cat("Canonical gene universe: PASS\n")


## ============================================================================
## 8. Output preflight
## ============================================================================

dir.create(
  bundle_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

bundle_output_files <- file.path(
  bundle_output_dir,
  paste0(
    "SCP1884_",
    bundle_config$source_bundle,
    "_selected_raw_counts.h5"
  )
)

final_counts_file <- file.path(
  output_dir,
  "SCP1884_selected_raw_counts.h5"
)

manifest_output <- file.path(
  output_dir,
  "SCP1884_selected_cell_manifest.parquet"
)

gene_metadata_output <- file.path(
  output_dir,
  "SCP1884_gene_metadata.parquet"
)

bundle_audit_output <- file.path(
  output_dir,
  "SCP1884_bundle_extraction_audit.tsv"
)

gene_audit_output <- file.path(
  output_dir,
  "SCP1884_gene_harmonization_audit.tsv"
)

session_output <- file.path(
  output_dir,
  "SCP1884_raw_count_extraction_sessionInfo.txt"
)

target_outputs <- c(
  bundle_output_files,
  final_counts_file,
  manifest_output,
  gene_metadata_output,
  bundle_audit_output,
  gene_audit_output,
  session_output
)

existing_outputs <- target_outputs[
  file.exists(target_outputs)
]

if (length(existing_outputs) > 0L) {

  if (!overwrite) {

    stop(
      paste0(
        "Existing SCP1884 outputs detected.\n",
        "Refusing to overwrite.\n\n",
        paste(
          existing_outputs,
          collapse = "\n"
        ),
        "\n\nSet IBD_REF_OVERWRITE=1 only for an intentional clean rerun."
      )
    )

  } else {

    cat(
      "\nIBD_REF_OVERWRITE=1 detected.\n"
    )

    cat(
      "Removing existing SCP1884 outputs...\n"
    )

    unlink(
      existing_outputs,
      recursive = TRUE
    )
  }
}

cat("Output preflight: PASS\n")


## ============================================================================
## 9. Bundle extraction function
## ============================================================================

extract_scp1884_bundle <- function(
    bundle,
    subdir,
    expected_selected,
    expected_source_genes,
    expected_source_cells
) {

  cat("\n")
  cat("============================================================\n")
  cat("PROCESSING:", bundle, "\n")
  cat("============================================================\n")

  raw_file <- file.path(
    expression_root,
    subdir,
    "raw.mtx"
  )

  barcode_file <- file.path(
    expression_root,
    subdir,
    paste0(
      bundle,
      ".barcodes.tsv"
    )
  )

  feature_file <- file.path(
    expression_root,
    subdir,
    paste0(
      bundle,
      ".features.tsv"
    )
  )


  ## --------------------------------------------------------------------------
  ## Feature and barcode audit
  ## --------------------------------------------------------------------------

  features <- read_features(
    feature_file
  )

  barcodes <- read_barcodes(
    barcode_file
  )

  stopifnot(
    nrow(features) ==
      expected_source_genes,

    length(barcodes) ==
      expected_source_cells,

    !base::anyDuplicated(
      features$ensembl_id
    ),

    !base::anyDuplicated(
      barcodes
    )
  )

  ## Verify that all Colon bundles share the exact same feature universe,
  ## and likewise for all TI bundles.

  reference_features <- if (
    startsWith(
      bundle,
      "CO_"
    )
  ) {
    colon_features
  } else {
    ti_features
  }

  stopifnot(
    identical(
      features$ensembl_id,
      reference_features$ensembl_id
    ),

    identical(
      features$gene_symbol,
      reference_features$gene_symbol
    )
  )

  cat(
    "Features/barcodes:",
    nrow(features),
    "genes /",
    length(barcodes),
    "cells\n"
  )


  ## --------------------------------------------------------------------------
  ## Manifest mapping before reading matrix
  ## --------------------------------------------------------------------------

  manifest_sub <- manifest_scp1884[
    manifest_scp1884$source_bundle ==
      bundle,
    ,
    drop = FALSE
  ]

  stopifnot(
    nrow(manifest_sub) ==
      expected_selected,

    !base::anyDuplicated(
      manifest_sub$source_cell_key
    )
  )

  cell_idx <- base::match(
    manifest_sub$source_cell_key,
    barcodes
  )

  stopifnot(
    !anyNA(cell_idx),
    !base::anyDuplicated(cell_idx),

    identical(
      barcodes[cell_idx],
      manifest_sub$source_cell_key
    )
  )

  gene_idx <- base::match(
    common_gene_order,
    features$ensembl_id
  )

  stopifnot(
    !anyNA(gene_idx),
    !base::anyDuplicated(gene_idx)
  )

  cat(
    "Selected-cell mapping:",
    length(cell_idx),
    "/",
    expected_selected,
    "matched\n"
  )

  cat("Exact cell mapping: PASS\n")


  ## --------------------------------------------------------------------------
  ## Load raw MatrixMarket source
  ## --------------------------------------------------------------------------

  cat("Reading raw.mtx...\n")

  raw <- Matrix::readMM(
    raw_file
  )

  stopifnot(
    nrow(raw) ==
      expected_source_genes,

    ncol(raw) ==
      expected_source_cells
  )

  cat(
    "Source matrix:",
    paste(
      dim(raw),
      collapse = " x "
    ),
    "\n"
  )

  cat(
    "Source NNZ:",
    length(raw@x),
    "\n"
  )


  ## --------------------------------------------------------------------------
  ## Full raw-count audit
  ## --------------------------------------------------------------------------

  count_audit <- audit_sparse_values(
    raw@x
  )

  stopifnot(
    count_audit$n_negative == 0,
    count_audit$n_noninteger == 0,
    abs(
      count_audit$integer_fraction - 1
    ) < 1e-12
  )

  cat(
    "Raw-count audit:",
    "negative =",
    count_audit$n_negative,
    "| non-integer =",
    count_audit$n_noninteger,
    "| min =",
    count_audit$min_value,
    "| max =",
    count_audit$max_value,
    "\n"
  )

  cat("Raw-count audit: PASS\n")


  ## --------------------------------------------------------------------------
  ## Directly subset canonical genes AND selected cells
  ##
  ## Important:
  ##   Do not convert the entire source matrix to CsparseMatrix first.
  ##   The largest source bundles contain >160 million non-zero entries.
  ##   Subsetting first reduces peak memory usage.
  ## --------------------------------------------------------------------------

  cat(
    "Extracting canonical genes and selected cells...\n"
  )

  selected <- raw[
    gene_idx,
    cell_idx,
    drop = FALSE
  ]

  rm(raw)
  invisible(gc())

  selected <- methods::as(
    selected,
    "CsparseMatrix"
  )

  rownames(selected) <-
    common_gene_order

  colnames(selected) <-
    manifest_sub$global_cell_id

  stopifnot(
    nrow(selected) ==
      expected_common_genes,

    ncol(selected) ==
      expected_selected,

    identical(
      rownames(selected),
      common_gene_order
    ),

    identical(
      colnames(selected),
      manifest_sub$global_cell_id
    )
  )

  selected_nnz <- length(
    selected@x
  )

  cat(
    "Selected matrix:",
    paste(
      dim(selected),
      collapse = " x "
    ),
    "\n"
  )

  cat(
    "Selected NNZ:",
    selected_nnz,
    "\n"
  )

  cat("Selected extraction: PASS\n")


  ## --------------------------------------------------------------------------
  ## Bundle-level sparse HDF5 output
  ## --------------------------------------------------------------------------

  output_file <- file.path(
    bundle_output_dir,
    paste0(
      "SCP1884_",
      bundle,
      "_selected_raw_counts.h5"
    )
  )

  if (file.exists(output_file)) {
    stop(
      "Unexpected existing bundle output: ",
      output_file
    )
  }

  cat(
    "Writing bundle sparse HDF5...\n"
  )

  written <- HDF5Array::writeTENxMatrix(
    selected,
    filepath = output_file,
    group = "raw_counts",
    level = 6,
    verbose = TRUE
  )

  stopifnot(
    all(
      dim(written) ==
        c(
          expected_common_genes,
          expected_selected
        )
    ),

    identical(
      rownames(written),
      common_gene_order
    ),

    identical(
      colnames(written),
      manifest_sub$global_cell_id
    )
  )


  ## --------------------------------------------------------------------------
  ## Post-write concordance
  ## --------------------------------------------------------------------------

  bundle_position <- base::match(
    bundle,
    bundle_config$source_bundle
  )

  set.seed(
    20260914L +
      bundle_position
  )

  test_genes <- sample(
    seq_len(
      nrow(selected)
    ),
    100L
  )

  test_cells <- sample(
    seq_len(
      ncol(selected)
    ),
    min(
      100L,
      ncol(selected)
    )
  )

  before <- as.matrix(
    selected[
      test_genes,
      test_cells,
      drop = FALSE
    ]
  )

  after <- as.matrix(
    written[
      test_genes,
      test_cells,
      drop = FALSE
    ]
  )

  concordance <- identical(
    as.numeric(before),
    as.numeric(after)
  )

  stopifnot(concordance)

  cat("Post-write concordance: PASS\n")


  ## --------------------------------------------------------------------------
  ## Bundle audit record
  ## --------------------------------------------------------------------------

  audit <- data.frame(
    source_bundle = bundle,

    source_genes =
      expected_source_genes,

    source_cells =
      expected_source_cells,

    source_nnz =
      count_audit$n_values,

    common_genes =
      expected_common_genes,

    selected_cells =
      expected_selected,

    selected_nnz =
      selected_nnz,

    missing_selected_cells =
      sum(is.na(cell_idx)),

    negative_values =
      count_audit$n_negative,

    noninteger_values =
      count_audit$n_noninteger,

    integer_fraction =
      count_audit$integer_fraction,

    min_raw_count =
      count_audit$min_value,

    max_raw_count =
      count_audit$max_value,

    output_genes =
      nrow(written),

    output_cells =
      ncol(written),

    concordance =
      concordance,

    output_file =
      output_file,

    stringsAsFactors = FALSE
  )

  rm(
    selected,
    written,
    before,
    after,
    features,
    barcodes,
    manifest_sub
  )

  invisible(gc())

  cat(bundle, ": COMPLETE\n")

  audit
}


## ============================================================================
## 10. Sequential extraction of six bundles
## ============================================================================

cat("\n[4] Starting sequential bundle extraction...\n")

## Smaller bundles first.
run_order <- c(
  "CO_STR",
  "TI_STR",
  "CO_EPI",
  "CO_IMM",
  "TI_EPI",
  "TI_IMM"
)

audit_list <- vector(
  "list",
  length(run_order)
)

names(audit_list) <- run_order

for (bundle in run_order) {

  config <- bundle_config[
    bundle_config$source_bundle ==
      bundle,
    ,
    drop = FALSE
  ]

  stopifnot(
    nrow(config) == 1L
  )

  audit_list[[bundle]] <-
    extract_scp1884_bundle(
      bundle =
        config$source_bundle,

      subdir =
        config$subdir,

      expected_selected =
        config$expected_selected,

      expected_source_genes =
        config$expected_source_genes,

      expected_source_cells =
        config$expected_source_cells
    )
}


## ============================================================================
## 11. All-bundle audit
## ============================================================================

bundle_audit <- do.call(
  rbind,
  audit_list
)

rownames(bundle_audit) <- NULL

stopifnot(
  nrow(bundle_audit) == 6L,

  sum(
    bundle_audit$selected_cells
  ) ==
    expected_total_selected,

  all(
    bundle_audit$missing_selected_cells ==
      0L
  ),

  all(
    bundle_audit$negative_values ==
      0L
  ),

  all(
    bundle_audit$noninteger_values ==
      0L
  ),

  all(
    abs(
      bundle_audit$integer_fraction -
        1
    ) < 1e-12
  ),

  all(
    bundle_audit$concordance
  )
)

cat("\n")
cat("========================================\n")
cat("ALL SCP1884 BUNDLES: PASS\n")
cat("========================================\n")

print(bundle_audit)


## ============================================================================
## 12. Save bundle audit before final combination
## ============================================================================

utils::write.table(
  bundle_audit,
  bundle_audit_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ============================================================================
## 13. Re-open bundle matrices as disk-backed TENxMatrix objects
## ============================================================================

cat("\n[5] Re-opening bundle-level HDF5 matrices...\n")

bundle_mats <- lapply(
  bundle_config$source_bundle,
  function(bundle) {

    f <- file.path(
      bundle_output_dir,
      paste0(
        "SCP1884_",
        bundle,
        "_selected_raw_counts.h5"
      )
    )

    stopifnot(
      file.exists(f)
    )

    HDF5Array::TENxMatrix(
      f,
      group = "raw_counts"
    )
  }
)

names(bundle_mats) <-
  bundle_config$source_bundle

stopifnot(
  all(
    vapply(
      bundle_mats,
      function(x) {

        nrow(x) ==
          expected_common_genes &&
          identical(
            rownames(x),
            common_gene_order
          )
      },
      logical(1)
    )
  )
)

cat("All bundle gene orders identical: PASS\n")


## ============================================================================
## 14. Combine bundles
## ============================================================================

cat("\n[6] Combining six bundle matrices...\n")

combined_by_bundle <- do.call(
  DelayedArray::acbind,
  unname(bundle_mats)
)

stopifnot(
  nrow(combined_by_bundle) ==
    expected_common_genes,

  ncol(combined_by_bundle) ==
    expected_total_selected,

  !base::anyDuplicated(
    colnames(combined_by_bundle)
  )
)

cat(
  "Combined-by-bundle dimensions:",
  paste(
    dim(combined_by_bundle),
    collapse = " x "
  ),
  "\n"
)


## ============================================================================
## 15. Restore authoritative manifest order
## ============================================================================

cat("\n[7] Restoring authoritative manifest cell order...\n")

final_cell_idx <- base::match(
  manifest_scp1884$global_cell_id,
  colnames(combined_by_bundle)
)

stopifnot(
  !anyNA(final_cell_idx),
  !base::anyDuplicated(final_cell_idx)
)

scp1884_final <- combined_by_bundle[
  ,
  final_cell_idx,
  drop = FALSE
]

stopifnot(
  nrow(scp1884_final) ==
    expected_common_genes,

  ncol(scp1884_final) ==
    expected_total_selected,

  identical(
    rownames(scp1884_final),
    common_gene_order
  ),

  identical(
    colnames(scp1884_final),
    manifest_scp1884$global_cell_id
  )
)

cat(
  "Final matrix dimensions:",
  paste(
    dim(scp1884_final),
    collapse = " x "
  ),
  "\n"
)

cat("Final manifest ordering: PASS\n")


## ============================================================================
## 16. Final sparse HDF5 output
## ============================================================================

cat("\n[8] Writing final SCP1884 sparse HDF5 matrix...\n")

scp1884_h5 <- HDF5Array::writeTENxMatrix(
  scp1884_final,
  filepath = final_counts_file,
  group = "raw_counts",
  level = 6,
  verbose = TRUE
)

stopifnot(
  all(
    dim(scp1884_h5) ==
      c(
        expected_common_genes,
        expected_total_selected
      )
  ),

  identical(
    rownames(scp1884_h5),
    common_gene_order
  ),

  identical(
    colnames(scp1884_h5),
    manifest_scp1884$global_cell_id
  )
)

cat("Final sparse HDF5 write: PASS\n")


## ============================================================================
## 17. Final post-write concordance
## ============================================================================

cat("\n[9] Final post-write concordance audit...\n")

set.seed(20260914L)

test_genes <- sample(
  seq_len(
    nrow(scp1884_final)
  ),
  200L
)

test_cells <- sample(
  seq_len(
    ncol(scp1884_final)
  ),
  200L
)

before <- as.matrix(
  scp1884_final[
    test_genes,
    test_cells,
    drop = FALSE
  ]
)

after <- as.matrix(
  scp1884_h5[
    test_genes,
    test_cells,
    drop = FALSE
  ]
)

final_concordance <- identical(
  as.numeric(before),
  as.numeric(after)
)

stopifnot(
  final_concordance
)

cat("Final random concordance: PASS\n")


## ============================================================================
## 18. Save metadata
## ============================================================================

cat("\n[10] Writing final metadata...\n")

arrow::write_parquet(
  manifest_scp1884,
  manifest_output
)

arrow::write_parquet(
  common_gene_metadata,
  gene_metadata_output
)

gene_harmonization_audit <- data.frame(
  colon_total_genes =
    nrow(colon_features),

  ti_total_genes =
    nrow(ti_features),

  common_genes =
    length(common_gene_order),

  colon_only_genes =
    length(colon_only),

  ti_only_genes =
    length(ti_only),

  shared_symbol_disagreements =
    sum(
      colon_common$gene_symbol !=
        ti_common$gene_symbol
    ),

  canonical_order_source =
    "Colon feature order",

  canonical_gene_key =
    "Ensembl ID",

  stringsAsFactors = FALSE
)

utils::write.table(
  gene_harmonization_audit,
  gene_audit_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = session_output
)

cat("Metadata output: PASS\n")


## ============================================================================
## 19. Final audit
## ============================================================================

stopifnot(
  file.exists(
    final_counts_file
  ),

  file.exists(
    manifest_output
  ),

  file.exists(
    gene_metadata_output
  ),

  file.exists(
    bundle_audit_output
  ),

  file.exists(
    gene_audit_output
  ),

  file.exists(
    session_output
  ),

  final_concordance,

  identical(
    colnames(scp1884_h5),
    manifest_scp1884$global_cell_id
  ),

  identical(
    rownames(scp1884_h5),
    common_gene_order
  )
)


## ============================================================================
## 20. Final summary
## ============================================================================

cat("\n")
cat("========================================\n")
cat("SCP1884 EXTRACTION COMPLETE\n")
cat("========================================\n")

cat(
  "Source bundles    : 6\n"
)

cat(
  "Canonical genes   :",
  nrow(scp1884_h5),
  "\n"
)

cat(
  "Selected cells    :",
  ncol(scp1884_h5),
  "\n"
)

cat("Bundle mapping     : PASS\n")
cat("Raw-count audit    : PASS\n")
cat("Gene harmonization : PASS\n")
cat("Bundle writes      : PASS\n")
cat("Manifest ordering  : PASS\n")
cat("Final sparse write : PASS\n")
cat("Final concordance  : PASS\n")

cat("\nFinal output:\n")
cat(
  final_counts_file,
  "\n"
)

cat("\n")
cat("ALL CHECKS PASSED\n")
