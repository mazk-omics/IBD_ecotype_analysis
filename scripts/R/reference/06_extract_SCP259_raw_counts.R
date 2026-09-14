#!/usr/bin/env Rscript

## ============================================================================
## 06_extract_SCP259_raw_counts.R
##
## Extract selected raw counts from SCP259 source matrices.
##
## Source compartments:
##   Epithelial <- Epi
##   Stromal    <- Fib
##   Immune     <- Imm
##
## Source gene identifiers:
##   gene symbols only
##
## Canonical SCP259 gene universe:
##   intersection of Epi / Fib / Imm gene symbols
##   = 18,172 genes
##
## Final expected matrix:
##   18,172 genes x 64,222 selected cells
##
## Existing outputs are not overwritten unless:
##
##   IBD_REF_OVERWRITE=1 Rscript --vanilla \
##     scripts/R/reference/06_extract_SCP259_raw_counts.R
## ============================================================================


options(
  stringsAsFactors = FALSE,
  warn = 1
)

overwrite <- identical(
  Sys.getenv("IBD_REF_OVERWRITE"),
  "1"
)


## ============================================================================
## 1. Paths
## ============================================================================

project_root <- "/home/mazekai/IBD_EcoTyper"

expression_root <- file.path(
  project_root,
  "01_raw_data/scRNA/SCP259/raw/SCP259/expression"
)

manifest_path <- file.path(
  project_root,
  "03_reference/combined_reference/output",
  "05_final_cell_selection",
  "selected_cell_manifest.parquet"
)

output_dir <- file.path(
  project_root,
  "03_reference/combined_reference/output",
  "06_raw_count_extraction",
  "SCP259"
)

bundle_output_dir <- file.path(
  output_dir,
  "bundles"
)

expected_total_selected <- 64222L
expected_common_genes <- 18172L


## ============================================================================
## 2. Packages
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
    "Missing required packages: ",
    paste(
      required_pkgs[!pkg_ok],
      collapse = ", "
    )
  )
}


## ============================================================================
## 3. Helpers
## ============================================================================

read_genes <- function(path) {

  x <- utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(x) != 1L) {
    stop(
      "Expected one-column gene file: ",
      path
    )
  }

  colnames(x) <- "gene_symbol"

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


audit_sparse_values <- function(
    values,
    chunk_size = 5000000L
) {

  n <- length(values)

  if (n == 0L) {
    stop("Sparse matrix contains no non-zero values.")
  }

  n_negative <- 0L
  n_noninteger <- 0L

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
      1 - n_noninteger / n,
    min_value = min_value,
    max_value = max_value
  )
}


## ============================================================================
## 4. Bundle configuration
##
## IMPORTANT:
##   manifest_bundle != source file prefix
##
##   Epithelial -> Epi
##   Stromal    -> Fib
##   Immune     -> Imm
## ============================================================================

bundle_config <- data.frame(

  manifest_bundle = c(
    "Epithelial",
    "Stromal",
    "Immune"
  ),

  file_prefix = c(
    "Epi",
    "Fib",
    "Imm"
  ),

  source_hash = c(
    "5cdc540d328cee7a2efc2348",
    "5cdc540d328cee7a2efc2349",
    "5cdc540d328cee7a2efc234a"
  ),

  expected_selected = c(
    18098L,
    11652L,
    34472L
  ),

  expected_source_genes = c(
    20028L,
    19076L,
    20529L
  ),

  expected_source_cells = c(
    123006L,
    31872L,
    210614L
  ),

  stringsAsFactors = FALSE
)

stopifnot(
  sum(bundle_config$expected_selected) ==
    expected_total_selected
)

bundle_config$source_dir <- file.path(
  expression_root,
  bundle_config$source_hash
)

bundle_config$matrix_file <- file.path(
  bundle_config$source_dir,
  paste0(
    "gene_sorted-",
    bundle_config$file_prefix,
    ".matrix.mtx"
  )
)

bundle_config$barcode_file <- file.path(
  bundle_config$source_dir,
  paste0(
    bundle_config$file_prefix,
    ".barcodes2.tsv"
  )
)

bundle_config$gene_file <- file.path(
  bundle_config$source_dir,
  paste0(
    bundle_config$file_prefix,
    ".genes.tsv"
  )
)


## ============================================================================
## 5. Source file checks
## ============================================================================

stopifnot(
  file.exists(manifest_path),
  all(file.exists(bundle_config$matrix_file)),
  all(file.exists(bundle_config$barcode_file)),
  all(file.exists(bundle_config$gene_file))
)

cat("\n")
cat("========================================\n")
cat("SCP259 RAW-COUNT EXTRACTION\n")
cat("========================================\n")

cat("All required source files: PASS\n")


## ============================================================================
## 6. Manifest audit
## ============================================================================

manifest <- as.data.frame(
  arrow::read_parquet(
    manifest_path
  )
)

manifest_scp259 <- manifest[
  manifest$dataset == "SCP259",
  ,
  drop = FALSE
]

stopifnot(
  nrow(manifest_scp259) ==
    expected_total_selected,

  !base::anyDuplicated(
    manifest_scp259$source_cell_key
  ),

  !base::anyDuplicated(
    manifest_scp259$global_cell_id
  ),

  all(
    manifest_scp259$global_cell_id ==
      paste0(
        "SCP259::",
        manifest_scp259$source_cell_key
      )
  )
)

observed_bundle_counts <- as.integer(
  table(
    factor(
      manifest_scp259$source_bundle,
      levels = bundle_config$manifest_bundle
    )
  )
)

stopifnot(
  identical(
    observed_bundle_counts,
    bundle_config$expected_selected
  )
)

cat(
  "Selected SCP259 cells:",
  nrow(manifest_scp259),
  "\n"
)

cat("Manifest audit: PASS\n")


## ============================================================================
## 7. Load source gene universes
## ============================================================================

epi_genes <- read_genes(
  bundle_config$gene_file[
    bundle_config$file_prefix == "Epi"
  ]
)

fib_genes <- read_genes(
  bundle_config$gene_file[
    bundle_config$file_prefix == "Fib"
  ]
)

imm_genes <- read_genes(
  bundle_config$gene_file[
    bundle_config$file_prefix == "Imm"
  ]
)

stopifnot(
  nrow(epi_genes) == 20028L,
  nrow(fib_genes) == 19076L,
  nrow(imm_genes) == 20529L,

  !anyNA(epi_genes$gene_symbol),
  !anyNA(fib_genes$gene_symbol),
  !anyNA(imm_genes$gene_symbol),

  !any(epi_genes$gene_symbol == ""),
  !any(fib_genes$gene_symbol == ""),
  !any(imm_genes$gene_symbol == ""),

  !base::anyDuplicated(
    epi_genes$gene_symbol
  ),

  !base::anyDuplicated(
    fib_genes$gene_symbol
  ),

  !base::anyDuplicated(
    imm_genes$gene_symbol
  )
)


## ============================================================================
## 8. Canonical SCP259 gene universe
##
## Preserve Epi gene order.
## ============================================================================

common_gene_order <-
  epi_genes$gene_symbol[
    epi_genes$gene_symbol %in%
      fib_genes$gene_symbol &
      epi_genes$gene_symbol %in%
      imm_genes$gene_symbol
  ]

stopifnot(
  length(common_gene_order) ==
    expected_common_genes,

  !base::anyDuplicated(
    common_gene_order
  )
)

epi_gene_idx <- base::match(
  common_gene_order,
  epi_genes$gene_symbol
)

fib_gene_idx <- base::match(
  common_gene_order,
  fib_genes$gene_symbol
)

imm_gene_idx <- base::match(
  common_gene_order,
  imm_genes$gene_symbol
)

stopifnot(
  !anyNA(epi_gene_idx),
  !anyNA(fib_gene_idx),
  !anyNA(imm_gene_idx),

  !base::anyDuplicated(epi_gene_idx),
  !base::anyDuplicated(fib_gene_idx),
  !base::anyDuplicated(imm_gene_idx),

  identical(
    epi_genes$gene_symbol[
      epi_gene_idx
    ],
    common_gene_order
  ),

  identical(
    fib_genes$gene_symbol[
      fib_gene_idx
    ],
    common_gene_order
  ),

  identical(
    imm_genes$gene_symbol[
      imm_gene_idx
    ],
    common_gene_order
  )
)

common_gene_metadata <- data.frame(
  gene_symbol = common_gene_order,
  stringsAsFactors = FALSE
)

cat("\nGene universe:\n")
cat("Epi genes    :", nrow(epi_genes), "\n")
cat("Fib genes    :", nrow(fib_genes), "\n")
cat("Imm genes    :", nrow(imm_genes), "\n")
cat("Common genes :", length(common_gene_order), "\n")

cat("Canonical gene universe: PASS\n")


## ============================================================================
## 9. Output preflight
## ============================================================================

dir.create(
  bundle_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

bundle_output_files <- file.path(
  bundle_output_dir,
  paste0(
    "SCP259_",
    bundle_config$manifest_bundle,
    "_selected_raw_counts.h5"
  )
)

final_counts_file <- file.path(
  output_dir,
  "SCP259_selected_raw_counts.h5"
)

manifest_output <- file.path(
  output_dir,
  "SCP259_selected_cell_manifest.parquet"
)

gene_metadata_output <- file.path(
  output_dir,
  "SCP259_gene_metadata.parquet"
)

bundle_audit_output <- file.path(
  output_dir,
  "SCP259_bundle_extraction_audit.tsv"
)

gene_audit_output <- file.path(
  output_dir,
  "SCP259_gene_harmonization_audit.tsv"
)

session_output <- file.path(
  output_dir,
  "SCP259_raw_count_extraction_sessionInfo.txt"
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
        "Existing SCP259 outputs detected.\n",
        "Refusing to overwrite:\n",
        paste(
          existing_outputs,
          collapse = "\n"
        )
      )
    )

  } else {

    cat(
      "IBD_REF_OVERWRITE=1 detected.\n"
    )

    unlink(
      existing_outputs,
      recursive = TRUE
    )
  }
}

cat("Output preflight: PASS\n")


## ============================================================================
## 10. Extraction function
## ============================================================================

extract_scp259_bundle <- function(
    manifest_bundle,
    file_prefix,
    matrix_file,
    barcode_file,
    gene_file,
    expected_selected,
    expected_source_genes,
    expected_source_cells
) {

  cat("\n")
  cat("============================================================\n")
  cat(
    "PROCESSING:",
    manifest_bundle,
    "<-",
    file_prefix,
    "\n"
  )
  cat("============================================================\n")


  ## --------------------------------------------------------------------------
  ## Genes / barcodes
  ## --------------------------------------------------------------------------

  genes <- read_genes(
    gene_file
  )

  barcodes <- read_barcodes(
    barcode_file
  )

  stopifnot(
    nrow(genes) ==
      expected_source_genes,

    length(barcodes) ==
      expected_source_cells,

    !base::anyDuplicated(
      genes$gene_symbol
    ),

    !base::anyDuplicated(
      barcodes
    )
  )

  expected_gene_vector <- switch(
    file_prefix,
    Epi = epi_genes$gene_symbol,
    Fib = fib_genes$gene_symbol,
    Imm = imm_genes$gene_symbol,
    stop(
      "Unexpected file_prefix: ",
      file_prefix
    )
  )

  stopifnot(
    identical(
      genes$gene_symbol,
      expected_gene_vector
    )
  )

  cat(
    "Source:",
    nrow(genes),
    "genes /",
    length(barcodes),
    "cells\n"
  )


  ## --------------------------------------------------------------------------
  ## Manifest mapping
  ## --------------------------------------------------------------------------

  manifest_sub <- manifest_scp259[
    manifest_scp259$source_bundle ==
      manifest_bundle,
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

    !base::anyDuplicated(
      cell_idx
    ),

    identical(
      barcodes[cell_idx],
      manifest_sub$source_cell_key
    )
  )

  gene_idx <- base::match(
    common_gene_order,
    genes$gene_symbol
  )

  stopifnot(
    !anyNA(gene_idx),

    !base::anyDuplicated(
      gene_idx
    ),

    identical(
      genes$gene_symbol[
        gene_idx
      ],
      common_gene_order
    )
  )

  cat(
    "Selected-cell mapping:",
    length(cell_idx),
    "/",
    expected_selected,
    "\n"
  )

  cat("Exact mapping: PASS\n")


  ## --------------------------------------------------------------------------
  ## Read source raw counts
  ## --------------------------------------------------------------------------

  cat("Reading raw MatrixMarket matrix...\n")

  raw <- Matrix::readMM(
    matrix_file
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
    "Raw counts:",
    "negative =",
    count_audit$n_negative,
    "| noninteger =",
    count_audit$n_noninteger,
    "| min =",
    count_audit$min_value,
    "| max =",
    count_audit$max_value,
    "\n"
  )

  cat("Full raw-count audit: PASS\n")


  ## --------------------------------------------------------------------------
  ## Subset first, then convert sparse representation
  ## --------------------------------------------------------------------------

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

  cat("Selected extraction: PASS\n")


  ## --------------------------------------------------------------------------
  ## Write bundle HDF5
  ## --------------------------------------------------------------------------

  output_file <- file.path(
    bundle_output_dir,
    paste0(
      "SCP259_",
      manifest_bundle,
      "_selected_raw_counts.h5"
    )
  )

  cat("Writing bundle HDF5...\n")

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

  set.seed(
    20260914L +
      base::match(
        manifest_bundle,
        bundle_config$manifest_bundle
      )
  )

  test_genes <- sample(
    seq_len(nrow(selected)),
    100L
  )

  test_cells <- sample(
    seq_len(ncol(selected)),
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
  ## Audit
  ## --------------------------------------------------------------------------

  audit <- data.frame(

    manifest_bundle =
      manifest_bundle,

    source_prefix =
      file_prefix,

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
    genes,
    barcodes,
    manifest_sub
  )

  invisible(gc())

  cat(
    manifest_bundle,
    ": COMPLETE\n"
  )

  audit
}


## ============================================================================
## 11. Sequential extraction
##
## Fib is smallest, so process it first.
## ============================================================================

cat("\nStarting SCP259 bundle extraction...\n")

run_order <- c(
  "Stromal",
  "Epithelial",
  "Immune"
)

audit_list <- vector(
  "list",
  length(run_order)
)

names(audit_list) <- run_order

for (bundle in run_order) {

  config <- bundle_config[
    bundle_config$manifest_bundle ==
      bundle,
    ,
    drop = FALSE
  ]

  stopifnot(
    nrow(config) == 1L
  )

  audit_list[[bundle]] <-
    extract_scp259_bundle(

      manifest_bundle =
        config$manifest_bundle,

      file_prefix =
        config$file_prefix,

      matrix_file =
        config$matrix_file,

      barcode_file =
        config$barcode_file,

      gene_file =
        config$gene_file,

      expected_selected =
        config$expected_selected,

      expected_source_genes =
        config$expected_source_genes,

      expected_source_cells =
        config$expected_source_cells
    )
}


## ============================================================================
## 12. Bundle audit
## ============================================================================

bundle_audit <- do.call(
  rbind,
  audit_list
)

rownames(bundle_audit) <- NULL

stopifnot(
  nrow(bundle_audit) == 3L,

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
    bundle_audit$concordance
  )
)

cat("\n")
cat("========================================\n")
cat("ALL SCP259 BUNDLES: PASS\n")
cat("========================================\n")

print(bundle_audit)

utils::write.table(
  bundle_audit,
  bundle_audit_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ============================================================================
## 13. Re-open bundle HDF5 matrices
## ============================================================================

cat("\nRe-opening bundle HDF5 matrices...\n")

bundle_order <- bundle_config$manifest_bundle

bundle_mats <- lapply(
  bundle_order,
  function(bundle) {

    f <- file.path(
      bundle_output_dir,
      paste0(
        "SCP259_",
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

names(bundle_mats) <- bundle_order

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
## 14. Disk-backed combination
## ============================================================================

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
  "Combined matrix:",
  paste(
    dim(combined_by_bundle),
    collapse = " x "
  ),
  "\n"
)

cat("Bundle combination: PASS\n")


## ============================================================================
## 15. Restore authoritative manifest order
## ============================================================================

final_cell_idx <- base::match(
  manifest_scp259$global_cell_id,
  colnames(combined_by_bundle)
)

stopifnot(
  !anyNA(final_cell_idx),
  !base::anyDuplicated(final_cell_idx)
)

scp259_final <- combined_by_bundle[
  ,
  final_cell_idx,
  drop = FALSE
]

stopifnot(
  nrow(scp259_final) ==
    expected_common_genes,

  ncol(scp259_final) ==
    expected_total_selected,

  identical(
    rownames(scp259_final),
    common_gene_order
  ),

  identical(
    colnames(scp259_final),
    manifest_scp259$global_cell_id
  )
)

cat("Authoritative manifest order: PASS\n")


## ============================================================================
## 16. Final HDF5
## ============================================================================

cat("Writing final SCP259 HDF5...\n")

scp259_h5 <- HDF5Array::writeTENxMatrix(
  scp259_final,
  filepath = final_counts_file,
  group = "raw_counts",
  level = 6,
  verbose = TRUE
)

stopifnot(
  all(
    dim(scp259_h5) ==
      c(
        expected_common_genes,
        expected_total_selected
      )
  ),

  identical(
    rownames(scp259_h5),
    common_gene_order
  ),

  identical(
    colnames(scp259_h5),
    manifest_scp259$global_cell_id
  )
)

cat("Final sparse write: PASS\n")


## ============================================================================
## 17. Independent final concordance
## ============================================================================

set.seed(20260914L)

test_genes <- sample(
  seq_len(
    nrow(scp259_final)
  ),
  200L
)

test_cells <- sample(
  seq_len(
    ncol(scp259_final)
  ),
  200L
)

before <- as.matrix(
  scp259_final[
    test_genes,
    test_cells,
    drop = FALSE
  ]
)

after <- as.matrix(
  scp259_h5[
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
## 18. Metadata
## ============================================================================

arrow::write_parquet(
  manifest_scp259,
  manifest_output
)

arrow::write_parquet(
  common_gene_metadata,
  gene_metadata_output
)

gene_harmonization_audit <- data.frame(

  Epi_total_genes =
    nrow(epi_genes),

  Fib_total_genes =
    nrow(fib_genes),

  Imm_total_genes =
    nrow(imm_genes),

  common_genes =
    length(common_gene_order),

  Epi_not_in_common =
    nrow(epi_genes) -
      length(common_gene_order),

  Fib_not_in_common =
    nrow(fib_genes) -
      length(common_gene_order),

  Imm_not_in_common =
    nrow(imm_genes) -
      length(common_gene_order),

  canonical_gene_key =
    "gene_symbol",

  canonical_order_source =
    "Epi gene order",

  harmonization_strategy =
    "Three-way intersection of Epi, Fib and Imm gene symbols",

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

  nrow(scp259_h5) ==
    expected_common_genes,

  ncol(scp259_h5) ==
    expected_total_selected,

  identical(
    rownames(scp259_h5),
    common_gene_order
  ),

  identical(
    colnames(scp259_h5),
    manifest_scp259$global_cell_id
  ),

  final_concordance
)


## ============================================================================
## 20. Summary
## ============================================================================

cat("\n")
cat("========================================\n")
cat("SCP259 EXTRACTION COMPLETE\n")
cat("========================================\n")

cat("Source bundles      : 3\n")
cat(
  "Canonical genes     :",
  nrow(scp259_h5),
  "\n"
)
cat(
  "Selected cells      :",
  ncol(scp259_h5),
  "\n"
)

cat("Raw-count integrity : PASS\n")
cat("Cell mapping        : PASS\n")
cat("Gene harmonization  : PASS\n")
cat("Bundle writes       : PASS\n")
cat("Bundle combination  : PASS\n")
cat("Manifest ordering   : PASS\n")
cat("Final HDF5          : PASS\n")
cat("Final concordance   : PASS\n")

cat("\nFinal output:\n")
cat(
  final_counts_file,
  "\n"
)

cat("\nALL CHECKS PASSED\n")
