#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

required_pkgs <- c(
  "arrow",
  "HDF5Array",
  "DelayedArray"
)

ok <- vapply(
  required_pkgs,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

if (!all(ok)) {
  stop(
    "Missing packages: ",
    paste(required_pkgs[!ok], collapse = ", ")
  )
}


## ============================================================================
## Paths
## ============================================================================

project_root <- "/home/mazekai/IBD_EcoTyper"

raw_root <- file.path(
  project_root,
  "03_reference/combined_reference/output",
  "06_raw_count_extraction"
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
  "07_unified_reference"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================================
## Input paths
## ============================================================================

gse_h5 <- file.path(
  raw_root,
  "GSE282122",
  "GSE282122_selected_raw_counts.h5"
)

scp1884_h5 <- file.path(
  raw_root,
  "SCP1884",
  "SCP1884_selected_raw_counts.h5"
)

scp259_h5 <- file.path(
  raw_root,
  "SCP259",
  "SCP259_selected_raw_counts.h5"
)


gse_gene_file <- file.path(
  raw_root,
  "GSE282122",
  "GSE282122_gene_metadata.parquet"
)

scp1884_gene_file <- file.path(
  raw_root,
  "SCP1884",
  "SCP1884_gene_metadata.parquet"
)

scp259_gene_file <- file.path(
  raw_root,
  "SCP259",
  "SCP259_gene_metadata.parquet"
)


stopifnot(
  file.exists(gse_h5),
  file.exists(scp1884_h5),
  file.exists(scp259_h5),
  file.exists(gse_gene_file),
  file.exists(scp1884_gene_file),
  file.exists(scp259_gene_file),
  file.exists(manifest_path)
)


## ============================================================================
## Load final matrices
## ============================================================================

cat("\n========================================\n")
cat("UNIFIED REFERENCE BUILD\n")
cat("========================================\n")

cat("\n[1] Opening final dataset matrices...\n")

gse <- HDF5Array::TENxMatrix(
  gse_h5,
  group = "raw_counts"
)

scp1884 <- HDF5Array::TENxMatrix(
  scp1884_h5,
  group = "raw_counts"
)

scp259 <- HDF5Array::TENxMatrix(
  scp259_h5,
  group = "raw_counts"
)

cat(
  "GSE282122:",
  paste(dim(gse), collapse = " x "),
  "\n"
)

cat(
  "SCP1884  :",
  paste(dim(scp1884), collapse = " x "),
  "\n"
)

cat(
  "SCP259   :",
  paste(dim(scp259), collapse = " x "),
  "\n"
)

stopifnot(
  ncol(gse) == 103538L,
  ncol(scp1884) == 155384L,
  ncol(scp259) == 64222L
)

cat("Dataset dimensions: PASS\n")


## ============================================================================
## Load gene metadata
## ============================================================================

cat("\n[2] Loading gene metadata...\n")

gse_gene <- as.data.frame(
  arrow::read_parquet(
    gse_gene_file
  )
)

scp1884_gene <- as.data.frame(
  arrow::read_parquet(
    scp1884_gene_file
  )
)

scp259_gene <- as.data.frame(
  arrow::read_parquet(
    scp259_gene_file
  )
)


## ============================================================================
## GSE282122: matrix rows -> gene symbol
## ============================================================================

stopifnot(
  "source_gene_id" %in% colnames(gse_gene),
  "gene_symbol" %in% colnames(gse_gene)
)

gse_meta_idx <- base::match(
  rownames(gse),
  gse_gene$source_gene_id
)

stopifnot(
  !anyNA(gse_meta_idx),
  !anyDuplicated(gse_meta_idx)
)

gse_symbol <- gse_gene$gene_symbol[
  gse_meta_idx
]

stopifnot(
  !anyNA(gse_symbol),
  !any(gse_symbol == "")
)

cat(
  "GSE282122 duplicated symbols:",
  sum(duplicated(gse_symbol)),
  "\n"
)


## ============================================================================
## SCP1884: matrix Ensembl rows -> gene symbol
## ============================================================================

stopifnot(
  "ensembl_id" %in% colnames(scp1884_gene),
  "gene_symbol" %in% colnames(scp1884_gene)
)

scp1884_meta_idx <- base::match(
  rownames(scp1884),
  scp1884_gene$ensembl_id
)

stopifnot(
  !anyNA(scp1884_meta_idx),
  !anyDuplicated(scp1884_meta_idx)
)

scp1884_symbol <-
  scp1884_gene$gene_symbol[
    scp1884_meta_idx
  ]

stopifnot(
  !anyNA(scp1884_symbol),
  !any(scp1884_symbol == "")
)

cat(
  "SCP1884 duplicated symbols:",
  sum(duplicated(scp1884_symbol)),
  "\n"
)


## ============================================================================
## SCP259: matrix rows already gene symbols
## ============================================================================

stopifnot(
  "gene_symbol" %in% colnames(scp259_gene)
)

scp259_meta_idx <- base::match(
  rownames(scp259),
  scp259_gene$gene_symbol
)

stopifnot(
  !anyNA(scp259_meta_idx),
  !anyDuplicated(scp259_meta_idx)
)

scp259_symbol <-
  scp259_gene$gene_symbol[
    scp259_meta_idx
  ]

stopifnot(
  identical(
    rownames(scp259),
    scp259_symbol
  )
)

cat(
  "SCP259 duplicated symbols:",
  sum(duplicated(scp259_symbol)),
  "\n"
)


## ============================================================================
## Gene-symbol uniqueness is required
## ============================================================================

stopifnot(
  !anyDuplicated(gse_symbol),
  !anyDuplicated(scp1884_symbol),
  !anyDuplicated(scp259_symbol)
)

cat("Cross-dataset gene-symbol uniqueness: PASS\n")


## ============================================================================
## Three-way gene intersection
##
## Preserve GSE282122 row order as canonical order.
## ============================================================================

cat("\n[3] Calculating cross-dataset intersection...\n")

common_gene_order <- gse_symbol[
  gse_symbol %in% scp1884_symbol &
  gse_symbol %in% scp259_symbol
]

stopifnot(
  !anyNA(common_gene_order),
  !anyDuplicated(common_gene_order)
)

cat(
  "GSE282122 symbols :",
  length(gse_symbol),
  "\n"
)

cat(
  "SCP1884 symbols   :",
  length(scp1884_symbol),
  "\n"
)

cat(
  "SCP259 symbols    :",
  length(scp259_symbol),
  "\n"
)

cat(
  "Three-way common  :",
  length(common_gene_order),
  "\n"
)

cat(
  "GSE not common    :",
  length(gse_symbol) -
    length(common_gene_order),
  "\n"
)

cat(
  "SCP1884 not common:",
  length(scp1884_symbol) -
    length(common_gene_order),
  "\n"
)

cat(
  "SCP259 not common :",
  length(scp259_symbol) -
    length(common_gene_order),
  "\n"
)


## ============================================================================
## Match common genes in all matrices
## ============================================================================

gse_idx <- base::match(
  common_gene_order,
  gse_symbol
)

scp1884_idx <- base::match(
  common_gene_order,
  scp1884_symbol
)

scp259_idx <- base::match(
  common_gene_order,
  scp259_symbol
)

stopifnot(
  !anyNA(gse_idx),
  !anyNA(scp1884_idx),
  !anyNA(scp259_idx),

  !anyDuplicated(gse_idx),
  !anyDuplicated(scp1884_idx),
  !anyDuplicated(scp259_idx),

  identical(
    gse_symbol[gse_idx],
    common_gene_order
  ),

  identical(
    scp1884_symbol[scp1884_idx],
    common_gene_order
  ),

  identical(
    scp259_symbol[scp259_idx],
    common_gene_order
  )
)

cat("Three-dataset gene matching: PASS\n")


## ============================================================================
## Subset and rename rows to canonical gene symbol
## ============================================================================

cat("\n[4] Subsetting matrices to common gene space...\n")

gse_common <- gse[
  gse_idx,
  ,
  drop = FALSE
]

scp1884_common <- scp1884[
  scp1884_idx,
  ,
  drop = FALSE
]

scp259_common <- scp259[
  scp259_idx,
  ,
  drop = FALSE
]

rownames(gse_common) <-
  common_gene_order

rownames(scp1884_common) <-
  common_gene_order

rownames(scp259_common) <-
  common_gene_order

stopifnot(
  identical(
    rownames(gse_common),
    common_gene_order
  ),

  identical(
    rownames(scp1884_common),
    common_gene_order
  ),

  identical(
    rownames(scp259_common),
    common_gene_order
  )
)

cat("Dataset gene harmonization: PASS\n")


## ============================================================================
## Combine disk-backed matrices
## ============================================================================

cat("\n[5] Combining three references...\n")

unified_by_dataset <- DelayedArray::acbind(
  gse_common,
  scp1884_common,
  scp259_common
)

stopifnot(
  ncol(unified_by_dataset) == 323144L,
  nrow(unified_by_dataset) ==
    length(common_gene_order),
  !anyDuplicated(
    colnames(unified_by_dataset)
  )
)

cat(
  "Combined dimensions:",
  paste(
    dim(unified_by_dataset),
    collapse = " x "
  ),
  "\n"
)

cat("Three-reference combination: PASS\n")


## ============================================================================
## Restore authoritative global manifest order
## ============================================================================

cat("\n[6] Restoring authoritative cell order...\n")

manifest <- as.data.frame(
  arrow::read_parquet(
    manifest_path
  )
)

stopifnot(
  nrow(manifest) == 323144L,
  !anyDuplicated(
    manifest$global_cell_id
  ),
  sum(
    is.na(
      manifest$reference_identity
    ) |
      manifest$reference_identity == ""
  ) == 0L
)

final_cell_idx <- base::match(
  manifest$global_cell_id,
  colnames(unified_by_dataset)
)

stopifnot(
  !anyNA(final_cell_idx),
  !anyDuplicated(final_cell_idx)
)

unified_final <- unified_by_dataset[
  ,
  final_cell_idx,
  drop = FALSE
]

stopifnot(
  ncol(unified_final) == 323144L,

  identical(
    colnames(unified_final),
    manifest$global_cell_id
  ),

  identical(
    rownames(unified_final),
    common_gene_order
  )
)

cat("Global manifest ordering: PASS\n")


## ============================================================================
## Output
## ============================================================================

final_h5 <- file.path(
  output_dir,
  "unified_reference_raw_counts.h5"
)

if (file.exists(final_h5)) {
  stop(
    "Output already exists: ",
    final_h5
  )
}

cat("\n[7] Writing unified reference HDF5...\n")

unified_h5 <- HDF5Array::writeTENxMatrix(
  unified_final,
  filepath = final_h5,
  group = "raw_counts",
  level = 6,
  verbose = TRUE
)

stopifnot(
  identical(
    dim(unified_h5),
    dim(unified_final)
  ),

  identical(
    rownames(unified_h5),
    common_gene_order
  ),

  identical(
    colnames(unified_h5),
    manifest$global_cell_id
  )
)

cat("Unified HDF5 structure: PASS\n")


## ============================================================================
## Random concordance
## ============================================================================

set.seed(20260914L)

test_genes <- sample(
  seq_len(nrow(unified_final)),
  min(200L, nrow(unified_final))
)

test_cells <- sample(
  seq_len(ncol(unified_final)),
  200L
)

before <- as.matrix(
  unified_final[
    test_genes,
    test_cells,
    drop = FALSE
  ]
)

after <- as.matrix(
  unified_h5[
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

cat("Final random concordance: PASS\n")


## ============================================================================
## Metadata / audits
## ============================================================================

arrow::write_parquet(
  manifest,
  file.path(
    output_dir,
    "unified_reference_cell_manifest.parquet"
  )
)

gene_metadata <- data.frame(
  gene_symbol =
    common_gene_order,
  stringsAsFactors = FALSE
)

arrow::write_parquet(
  gene_metadata,
  file.path(
    output_dir,
    "unified_reference_gene_metadata.parquet"
  )
)

gene_audit <- data.frame(
  GSE282122_genes =
    length(gse_symbol),

  SCP1884_genes =
    length(scp1884_symbol),

  SCP259_genes =
    length(scp259_symbol),

  common_genes =
    length(common_gene_order),

  canonical_gene_key =
    "gene_symbol",

  canonical_order =
    "GSE282122 final matrix order",

  stringsAsFactors = FALSE
)

utils::write.table(
  gene_audit,
  file.path(
    output_dir,
    "unified_reference_gene_harmonization_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

identity_counts <- as.data.frame(
  table(
    manifest$reference_identity,
    manifest$dataset
  )
)

colnames(identity_counts) <- c(
  "reference_identity",
  "dataset",
  "n_cells"
)

utils::write.table(
  identity_counts,
  file.path(
    output_dir,
    "unified_reference_identity_counts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "unified_reference_sessionInfo.txt"
  )
)


## ============================================================================
## Summary
## ============================================================================

cat("\n")
cat("========================================\n")
cat("UNIFIED REFERENCE COMPLETE\n")
cat("========================================\n")

cat(
  "Genes      :",
  nrow(unified_h5),
  "\n"
)

cat(
  "Cells      :",
  ncol(unified_h5),
  "\n"
)

cat(
  "Identities :",
  length(
    unique(
      manifest$reference_identity
    )
  ),
  "\n"
)

cat("Gene harmonization : PASS\n")
cat("Cell harmonization : PASS\n")
cat("Manifest ordering  : PASS\n")
cat("HDF5 write         : PASS\n")
cat("Concordance        : PASS\n")

cat("\nALL CHECKS PASSED\n")
