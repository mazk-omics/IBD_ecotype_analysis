#!/usr/bin/env Rscript

## ============================================================================
## 08_build_candidate_marker_pseudobulk.R
##
## Build raw-count pseudobulk matrices for exhaustive marker auditing.
##
## Primary grouping:
##
##   dataset x biological_sample x coarse_identity
##
## coarse_identity = candidate cell type harmonized from author annotations.
##
## IMPORTANT
## ---------
## 1. Use the FULL candidate universe, not the sampled 323,144-cell reference.
## 2. Do NOT filter by primary_reference_action.
## 3. Do NOT group by primary_reference_action.
## 4. Preserve author annotations, state flags and historical actions as QC.
## 5. Do NOT normalize or calculate marker genes in this script.
##
## Candidate universe:
##
##   GSE282122 : 505,591 cells
##   SCP1884   : 720,633 cells
##   SCP259    : 365,492 cells
##
##   Total     : 1,591,716 cells
##
## ============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)


## ============================================================================
## 0. Configuration
## ============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

OUTDIR <- file.path(
  PROJECT_ROOT,
  "03_reference/combined_reference/output/08_marker_pseudobulk"
)

OVERWRITE <- identical(
  Sys.getenv("IBD_MARKER_PB_OVERWRITE"),
  "1"
)


## --------------------------------------------------------------------------
## Candidate metadata
## --------------------------------------------------------------------------

GSE_META <- file.path(
  PROJECT_ROOT,
  "03_reference/GSE282122/audit_v2/output/04_identity_coverage_v1/mapped_candidate_cells_v1.parquet"
)

SCP1884_META <- file.path(
  PROJECT_ROOT,
  "03_reference/SCP1884/audit_v1/output/02_candidate_taxonomy/scp1884_candidate_cells_provisional_taxonomy.parquet"
)

SCP259_META <- file.path(
  PROJECT_ROOT,
  "03_reference/SCP259/audit_v1/output/01_candidate_taxonomy_coverage/scp259_candidate_cells_provisional_taxonomy.parquet"
)


## --------------------------------------------------------------------------
## Raw-count sources
## --------------------------------------------------------------------------

GSE_H5AD <- paste0(
  "/home/mazekai/enteric_glia/data_scRNA/GSE282122/",
  "TAURUS_raw_counts_annotated_final.h5ad"
)

GSE_GENE_META <- file.path(
  PROJECT_ROOT,
  "03_reference/combined_reference/output/06_raw_count_extraction/GSE282122/GSE282122_gene_metadata.parquet"
)

SCP1884_EXPR <- paste0(
  "/home/mazekai/enteric_glia/",
  "data_scRNA/SCP1884/expression"
)

SCP259_EXPR <- file.path(
  PROJECT_ROOT,
  "01_raw_data/scRNA/SCP259/raw/SCP259/expression"
)


EXPECTED_CANDIDATES <- c(
  GSE282122 = 505591L,
  SCP1884 = 720633L,
  SCP259 = 365492L
)

EXPECTED_TOTAL <- sum(
  EXPECTED_CANDIDATES
)


## Sparse-processing parameters
GSE_CELL_BLOCK <- 2000L
MTX_NNZ_BLOCK <- 5000000L


## ============================================================================
## 1. Package checks
## ============================================================================

required_pkgs <- c(
  "arrow",
  "Matrix",
  "rhdf5",
  "zellkonverter",
  "SummarizedExperiment",
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
    ),
    call. = FALSE
  )
}


## ============================================================================
## 2. Helpers
## ============================================================================

fail <- function(...) {

  stop(
    ...,
    call. = FALSE
  )
}


clean_character <- function(x) {

  x <- as.character(x)

  x[
    is.na(x) |
      x == ""
  ] <- "<NA>"

  x
}


pick_col <- function(
    df,
    candidates,
    required = TRUE
) {

  hit <- candidates[
    candidates %in%
      colnames(df)
  ]

  if (length(hit) > 0L) {

    return(
      hit[[1L]]
    )
  }

  if (required) {

    fail(
      "None of these columns exists: ",
      paste(
        candidates,
        collapse = ", "
      )
    )
  }

  NA_character_
}


collapse_unique <- function(x) {

  paste(
    sort(
      unique(
        clean_character(x)
      )
    ),
    collapse = " | "
  )
}


collapse_counts <- function(x) {

  tab <- sort(
    table(
      clean_character(x)
    ),
    decreasing = TRUE
  )

  paste(
    paste0(
      names(tab),
      "=",
      as.integer(tab)
    ),
    collapse = " | "
  )
}


read_one_col <- function(path) {

  utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )[[1L]]
}


read_scp1884_features <- function(path) {

  utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c(
      "ensembl_id",
      "gene_symbol"
    )
  )
}


read_scp259_genes <- function(path) {

  x <- utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(x) != 1L) {

    fail(
      "Expected one-column gene file: ",
      path
    )
  }

  colnames(x) <- "gene_symbol"

  x
}


## ============================================================================
## 3. Standardize candidate metadata
## ============================================================================

standardize_gse <- function(df) {

  donor_col <- pick_col(
    df,
    c(
      "Patient",
      "donor_id"
    )
  )

  site_col <- pick_col(
    df,
    c(
      "Site",
      "Ileum_vs_Colon"
    )
  )

  disease_col <- pick_col(
    df,
    c(
      "Disease",
      "disease"
    )
  )

  inflammation_col <- pick_col(
    df,
    c(
      "Inflammation",
      "inflammation"
    )
  )

  path_cols <- c(
    "bucket",
    "sub_bucket",
    "major",
    "minor",
    "final_analysis"
  )

  path_cols <- path_cols[
    path_cols %in%
      colnames(df)
  ]

  if (length(path_cols) == 0L) {

    fail(
      "GSE282122 has no author-annotation path columns."
    )
  }

  annotation_parts <- lapply(
    df[path_cols],
    clean_character
  )

  source_annotation <- do.call(
    paste,
    c(
      annotation_parts,
      list(
        sep = " > "
      )
    )
  )

  data.frame(
    dataset = "GSE282122",
    cell_id = clean_character(
      df$cell_id
    ),
    donor_id = clean_character(
      df[[donor_col]]
    ),
    sample_id = clean_character(
      df$sample_id
    ),
    site = clean_character(
      df[[site_col]]
    ),
    disease = clean_character(
      df[[disease_col]]
    ),
    inflammation = clean_character(
      df[[inflammation_col]]
    ),
    coarse_identity = clean_character(
      df$coarse_identity
    ),
    primary_reference_action = clean_character(
      df$primary_reference_action
    ),
    state_flag = clean_character(
      df$state_flag
    ),
    source_annotation = source_annotation,
    source_bundle = "GSE282122_H5AD",
    stringsAsFactors = FALSE
  )
}


standardize_scp1884 <- function(df) {

  donor_col <- pick_col(
    df,
    c(
      "donor_id",
      "PubID"
    )
  )

  site_col <- pick_col(
    df,
    c(
      "site_canonical",
      "Site"
    )
  )

  bundle_col <- pick_col(
    df,
    c(
      "expression_bundle",
      "source_bundle"
    )
  )

  source_annotation <- paste(
    clean_character(
      df$anno_overall
    ),
    clean_character(
      df$anno2
    ),
    sep = " > "
  )

  data.frame(
    dataset = "SCP1884",
    cell_id = clean_character(
      df$cell_id
    ),
    donor_id = clean_character(
      df[[donor_col]]
    ),
    sample_id = clean_character(
      df$biological_sample_id
    ),
    site = clean_character(
      df[[site_col]]
    ),
    disease = clean_character(
      df$disease
    ),
    inflammation = clean_character(
      df$inflammation
    ),
    coarse_identity = clean_character(
      df$provisional_coarse_identity
    ),
    primary_reference_action = clean_character(
      df$primary_reference_action
    ),
    state_flag = clean_character(
      df$state_flag
    ),
    source_annotation = source_annotation,
    source_bundle = clean_character(
      df[[bundle_col]]
    ),
    stringsAsFactors = FALSE
  )
}


standardize_scp259 <- function(df) {

  bundle_col <- pick_col(
    df,
    c(
      "expression_bundle",
      "source_bundle"
    )
  )

  annotation_col <- pick_col(
    df,
    c(
      "author_annotation",
      "Cluster"
    )
  )

  cell_col <- pick_col(
    df,
    c(
      "cell_id",
      "NAME"
    )
  )

  data.frame(
    dataset = "SCP259",
    cell_id = clean_character(
      df[[cell_col]]
    ),
    donor_id = clean_character(
      df$donor_id
    ),
    sample_id = clean_character(
      df$sample_id
    ),
    site = clean_character(
      df$site
    ),
    disease = clean_character(
      df$disease
    ),
    inflammation = clean_character(
      df$inflammation
    ),
    coarse_identity = clean_character(
      df$coarse_identity
    ),
    primary_reference_action = clean_character(
      df$primary_reference_action
    ),
    state_flag = clean_character(
      df$state_flag
    ),
    source_annotation = clean_character(
      df[[annotation_col]]
    ),
    source_bundle = clean_character(
      df[[bundle_col]]
    ),
    stringsAsFactors = FALSE
  )
}


## ============================================================================
## 4. Build dataset x sample x coarse_identity design
## ============================================================================

build_pb_design <- function(meta) {

  required <- c(
    "dataset",
    "cell_id",
    "donor_id",
    "sample_id",
    "site",
    "disease",
    "inflammation",
    "coarse_identity",
    "primary_reference_action",
    "state_flag",
    "source_annotation",
    "source_bundle"
  )

  missing <- setdiff(
    required,
    colnames(meta)
  )

  if (length(missing) > 0L) {

    fail(
      "Metadata missing required columns: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  if (anyDuplicated(
    meta$cell_id
  )) {

    fail(
      "cell_id is not unique within dataset ",
      unique(meta$dataset)
    )
  }

  if (
    anyNA(meta$sample_id) ||
    any(meta$sample_id == "")
  ) {

    fail(
      "sample_id contains missing/blank values in ",
      unique(meta$dataset)
    )
  }

  if (
    anyNA(meta$coarse_identity) ||
    any(meta$coarse_identity == "")
  ) {

    fail(
      "coarse_identity contains missing/blank values in ",
      unique(meta$dataset)
    )
  }

  ## Primary grouping:
  ## dataset x biological sample x coarse identity

  meta$pb_id <- paste(
    meta$dataset,
    meta$sample_id,
    meta$coarse_identity,
    sep = "::"
  )

  group_levels <- unique(
    meta$pb_id
  )

  meta$pb_index <- match(
    meta$pb_id,
    group_levels
  )

  n_groups <- length(
    group_levels
  )

  if (anyNA(
    meta$pb_index
  )) {

    fail(
      "Failed to construct pb_index for ",
      unique(meta$dataset)
    )
  }

  first_idx <- match(
    group_levels,
    meta$pb_id
  )

  pb_meta <- meta[
    first_idx,
    c(
      "dataset",
      "donor_id",
      "sample_id",
      "site",
      "disease",
      "inflammation",
      "coarse_identity"
    ),
    drop = FALSE
  ]

  pb_meta$pb_id <- group_levels

  pb_meta$n_cells <- tabulate(
    meta$pb_index,
    nbins = n_groups
  )

  ## A biological sample must have internally consistent
  ## donor/site/disease/inflammation metadata.

  fixed_fields <- c(
    "donor_id",
    "sample_id",
    "site",
    "disease",
    "inflammation"
  )

  for (field in fixed_fields) {

    n_unique <- tapply(
      meta[[field]],
      meta$pb_index,
      function(z) {
        length(
          unique(z)
        )
      }
    )

    if (any(
      n_unique != 1L
    )) {

      fail(
        "Non-unique ",
        field,
        " inside pseudobulk groups for ",
        unique(meta$dataset)
      )
    }
  }

  idx_names <- as.character(
    seq_len(
      n_groups
    )
  )

  ## Preserve author annotation provenance.

  source_annotations <- tapply(
    meta$source_annotation,
    meta$pb_index,
    collapse_unique
  )

  source_annotation_composition <- tapply(
    meta$source_annotation,
    meta$pb_index,
    collapse_counts
  )

  n_source_annotations <- tapply(
    meta$source_annotation,
    meta$pb_index,
    function(z) {
      length(
        unique(
          clean_character(z)
        )
      )
    }
  )

  ## Preserve state information.

  state_flags <- tapply(
    meta$state_flag,
    meta$pb_index,
    collapse_unique
  )

  state_flag_composition <- tapply(
    meta$state_flag,
    meta$pb_index,
    collapse_counts
  )

  ## Preserve source bundle information.

  source_bundles <- tapply(
    meta$source_bundle,
    meta$pb_index,
    collapse_unique
  )

  ## Historical action is retained ONLY as provenance.
  ## It is NOT used for grouping or filtering.

  actions <- tapply(
    meta$primary_reference_action,
    meta$pb_index,
    collapse_unique
  )

  action_composition <- tapply(
    meta$primary_reference_action,
    meta$pb_index,
    collapse_counts
  )

  pb_meta$source_annotations <- unname(
    source_annotations[
      idx_names
    ]
  )

  pb_meta$source_annotation_composition <- unname(
    source_annotation_composition[
      idx_names
    ]
  )

  pb_meta$n_source_annotations <- as.integer(
    unname(
      n_source_annotations[
        idx_names
      ]
    )
  )

  pb_meta$state_flags <- unname(
    state_flags[
      idx_names
    ]
  )

  pb_meta$state_flag_composition <- unname(
    state_flag_composition[
      idx_names
    ]
  )

  pb_meta$source_bundles <- unname(
    source_bundles[
      idx_names
    ]
  )

  pb_meta$primary_reference_actions <- unname(
    actions[
      idx_names
    ]
  )

  pb_meta$action_composition <- unname(
    action_composition[
      idx_names
    ]
  )

  pb_meta <- pb_meta[
    c(
      "pb_id",
      "dataset",
      "donor_id",
      "sample_id",
      "site",
      "disease",
      "inflammation",
      "coarse_identity",
      "n_cells",
      "n_source_annotations",
      "source_annotations",
      "source_annotation_composition",
      "state_flags",
      "state_flag_composition",
      "source_bundles",
      "primary_reference_actions",
      "action_composition"
    )
  ]

  stopifnot(
    nrow(pb_meta) == n_groups,
    identical(
      pb_meta$pb_id,
      group_levels
    ),
    sum(
      pb_meta$n_cells
    ) == nrow(meta),
    !anyDuplicated(
      pb_meta$pb_id
    )
  )

  cat(
    "  ",
    format(
      nrow(meta),
      big.mark = ","
    ),
    " cells -> ",
    format(
      nrow(pb_meta),
      big.mark = ","
    ),
    " sample x coarse-identity pseudobulks\n",
    sep = ""
  )

  cat(
    "  Coarse identities: ",
    length(
      unique(
        meta$coarse_identity
      )
    ),
    "\n",
    sep = ""
  )

  list(
    cell_meta = meta,
    pb_meta = pb_meta
  )
}


## ============================================================================
## 5. Aggregate one MatrixMarket bundle
## ============================================================================

aggregate_triplet_bundle <- function(
    raw_file,
    barcode_file,
    source_gene_to_output,
    meta,
    n_output_genes,
    n_output_groups,
    require_all_source_cells = TRUE,
    nnz_block = MTX_NNZ_BLOCK
) {

  barcodes <- read_one_col(
    barcode_file
  )

  if (anyDuplicated(
    barcodes
  )) {

    fail(
      "Duplicated source barcodes: ",
      barcode_file
    )
  }

  cell_idx <- match(
    meta$cell_id,
    barcodes
  )

  if (anyNA(
    cell_idx
  )) {

    fail(
      "Candidate cells missing from source bundle: ",
      sum(
        is.na(
          cell_idx
        )
      ),
      "\n",
      raw_file
    )
  }

  if (anyDuplicated(
    cell_idx
  )) {

    fail(
      "Candidate-to-source mapping is not 1:1: ",
      raw_file
    )
  }

  if (require_all_source_cells) {

    if (
      nrow(meta) !=
        length(barcodes)
    ) {

      fail(
        "Candidate/source cell-count mismatch for ",
        raw_file,
        ": metadata=",
        nrow(meta),
        ", source=",
        length(barcodes)
      )
    }

    if (!identical(
      sort(
        cell_idx
      ),
      seq_along(
        barcodes
      )
    )) {

      fail(
        "Candidate metadata does not close exactly to source cells: ",
        raw_file
      )
    }
  }

  source_col_to_group <- integer(
    length(
      barcodes
    )
  )

  source_col_to_group[
    cell_idx
  ] <- meta$pb_index

  cat(
    "    Reading MatrixMarket source...\n"
  )

  raw <- Matrix::readMM(
    raw_file
  )

  if (!methods::is(
    raw,
    "TsparseMatrix"
  )) {

    raw <- methods::as(
      raw,
      "dgTMatrix"
    )
  }

  if (
    nrow(raw) !=
      length(
        source_gene_to_output
      )
  ) {

    fail(
      "Source gene dimension mismatch: ",
      raw_file
    )
  }

  if (
    ncol(raw) !=
      length(
        barcodes
      )
  ) {

    fail(
      "Source cell dimension mismatch: ",
      raw_file
    )
  }

  n_nnz <- length(
    raw@x
  )

  cat(
    "    Source matrix: ",
    format(
      nrow(raw),
      big.mark = ","
    ),
    " genes x ",
    format(
      ncol(raw),
      big.mark = ","
    ),
    " cells; NNZ=",
    format(
      n_nnz,
      big.mark = ","
    ),
    "\n",
    sep = ""
  )

  out <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    dims = c(
      n_output_genes,
      n_output_groups
    ),
    giveCsparse = TRUE
  )

  input_mass <- 0

  starts <- seq.int(
    1L,
    n_nnz,
    by = nnz_block
  )

  for (block_no in seq_along(
    starts
  )) {

    from <- starts[[block_no]
    ]

    to <- min(
      from + nnz_block - 1L,
      n_nnz
    )

    src_gene <- raw@i[
      from:to
    ] + 1L

    src_cell <- raw@j[
      from:to
    ] + 1L

    block_x <- raw@x[
      from:to
    ]

    dst_gene <- source_gene_to_output[
      src_gene
    ]

    dst_group <- source_col_to_group[
      src_cell
    ]

    keep <- (
      dst_gene > 0L &
        dst_group > 0L
    )

    if (any(
      keep
    )) {

      selected_x <- block_x[
        keep
      ]

      input_mass <- input_mass +
        sum(
          selected_x
        )

      partial <- Matrix::sparseMatrix(
        i = dst_gene[
          keep
        ],
        j = dst_group[
          keep
        ],
        x = selected_x,
        dims = c(
          n_output_genes,
          n_output_groups
        ),
        giveCsparse = TRUE
      )

      out <- out +
        partial

      rm(
        partial,
        selected_x
      )
    }

    rm(
      src_gene,
      src_cell,
      block_x,
      dst_gene,
      dst_group,
      keep
    )

    if (
      block_no %% 10L == 0L ||
      block_no == length(starts)
    ) {

      cat(
        "      NNZ block ",
        block_no,
        "/",
        length(starts),
        "\n",
        sep = ""
      )
    }

    invisible(
      gc(FALSE)
    )
  }

  rm(
    raw
  )

  invisible(
    gc()
  )

  output_mass <- sum(
    out
  )

  if (
    abs(
      output_mass -
        input_mass
    ) >
      1e-6
  ) {

    fail(
      "Raw-count mass closure failed for ",
      raw_file,
      ": input=",
      input_mass,
      ", output=",
      output_mass
    )
  }

  list(
    counts = out,
    input_mass = input_mass
  )
}


## ============================================================================
## 6. Save one dataset
## ============================================================================

save_dataset <- function(
    dataset,
    counts,
    pb_meta,
    gene_symbols,
    input_mass
) {

  dataset_dir <- file.path(
    OUTDIR,
    dataset
  )

  dir.create(
    dataset_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  counts_file <- file.path(
    dataset_dir,
    paste0(
      dataset,
      "_candidate_marker_pseudobulk_raw_counts.h5"
    )
  )

  meta_file <- file.path(
    dataset_dir,
    paste0(
      dataset,
      "_candidate_marker_pseudobulk_metadata.parquet"
    )
  )

  meta_tsv <- file.path(
    dataset_dir,
    paste0(
      dataset,
      "_candidate_marker_pseudobulk_metadata.tsv"
    )
  )

  gene_file <- file.path(
    dataset_dir,
    paste0(
      dataset,
      "_candidate_marker_pseudobulk_genes.tsv"
    )
  )

  audit_file <- file.path(
    dataset_dir,
    paste0(
      dataset,
      "_candidate_marker_pseudobulk_audit.tsv"
    )
  )

  targets <- c(
    counts_file,
    meta_file,
    meta_tsv,
    gene_file,
    audit_file
  )

  existing <- targets[
    file.exists(
      targets
    )
  ]

  if (length(
    existing
  ) > 0L) {

    if (!OVERWRITE) {

      fail(
        "Existing outputs detected for ",
        dataset,
        ". Set IBD_MARKER_PB_OVERWRITE=1 only for an intentional rerun."
      )
    }

    unlink(
      existing,
      recursive = TRUE
    )
  }

  stopifnot(
    nrow(counts) == length(gene_symbols),
    ncol(counts) == nrow(pb_meta),
    !anyDuplicated(gene_symbols),
    !anyDuplicated(pb_meta$pb_id)
  )

  rownames(
    counts
  ) <- gene_symbols

  colnames(
    counts
  ) <- pb_meta$pb_id

  if (methods::is(
    counts,
    "sparseMatrix"
  )) {

    values <- counts@x

    if (
      any(
        !is.finite(
          values
        )
      ) ||
      any(
        values < 0
      ) ||
      any(
        abs(
          values -
            round(
              values
            )
        ) >
          1e-8
      )
    ) {

      fail(
        dataset,
        ": pseudobulk matrix failed raw-count audit."
      )
    }
  }

  output_mass <- sum(
    counts
  )

  if (
    abs(
      output_mass -
        input_mass
    ) >
      1e-6
  ) {

    fail(
      dataset,
      ": dataset-level mass closure failed: input=",
      input_mass,
      ", output=",
      output_mass
    )
  }

  library_size <- as.numeric(
    Matrix::colSums(
      counts
    )
  )

  if (any(
    library_size <= 0
  )) {

    fail(
      dataset,
      ": at least one pseudobulk has zero library size."
    )
  }

  pb_meta$raw_library_size <-
    library_size

  cat(
    "  Writing ",
    dataset,
    " pseudobulk HDF5...\n",
    sep = ""
  )

  written <- HDF5Array::writeTENxMatrix(
    counts,
    filepath = counts_file,
    group = "raw_counts",
    level = 6,
    verbose = TRUE
  )

  stopifnot(
    identical(
      dim(written),
      dim(counts)
    ),
    identical(
      rownames(written),
      gene_symbols
    ),
    identical(
      colnames(written),
      pb_meta$pb_id
    )
  )

  set.seed(
    20261003L
  )

  test_genes <- sample(
    seq_len(
      nrow(counts)
    ),
    min(
      100L,
      nrow(counts)
    )
  )

  test_pb <- sample(
    seq_len(
      ncol(counts)
    ),
    min(
      100L,
      ncol(counts)
    )
  )

  before <- as.matrix(
    counts[
      test_genes,
      test_pb,
      drop = FALSE
    ]
  )

  after <- as.matrix(
    written[
      test_genes,
      test_pb,
      drop = FALSE
    ]
  )

  concordance <- identical(
    as.numeric(before),
    as.numeric(after)
  )

  if (!concordance) {

    fail(
      dataset,
      ": post-write concordance failed."
    )
  }

  arrow::write_parquet(
    pb_meta,
    meta_file
  )

  utils::write.table(
    pb_meta,
    meta_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  utils::write.table(
    data.frame(
      gene_symbol = gene_symbols,
      stringsAsFactors = FALSE
    ),
    gene_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  audit <- data.frame(
    dataset = dataset,
    n_genes = nrow(counts),
    n_pseudobulks = ncol(counts),
    n_candidate_cells = sum(
      pb_meta$n_cells
    ),
    n_samples = length(
      unique(
        pb_meta$sample_id
      )
    ),
    n_donors = length(
      unique(
        pb_meta$donor_id
      )
    ),
    n_coarse_identities = length(
      unique(
        pb_meta$coarse_identity
      )
    ),
    input_raw_count_mass = input_mass,
    output_raw_count_mass = output_mass,
    mass_closure = abs(
      input_mass -
        output_mass
    ) <= 1e-6,
    min_cells_per_pseudobulk = min(
      pb_meta$n_cells
    ),
    median_cells_per_pseudobulk = stats::median(
      pb_meta$n_cells
    ),
    max_cells_per_pseudobulk = max(
      pb_meta$n_cells
    ),
    min_library_size = min(
      library_size
    ),
    median_library_size = stats::median(
      library_size
    ),
    max_library_size = max(
      library_size
    ),
    postwrite_random_concordance = concordance,
    stringsAsFactors = FALSE
  )

  utils::write.table(
    audit,
    audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  audit
}


## ============================================================================
## 7. Input preflight
## ============================================================================

required_files <- c(
  GSE_META,
  SCP1884_META,
  SCP259_META,
  GSE_H5AD,
  GSE_GENE_META
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(
  missing_files
) > 0L) {

  fail(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

if (!dir.exists(
  SCP1884_EXPR
)) {

  fail(
    "Missing SCP1884 expression root: ",
    SCP1884_EXPR
  )
}

if (!dir.exists(
  SCP259_EXPR
)) {

  fail(
    "Missing SCP259 expression root: ",
    SCP259_EXPR
  )
}

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n")
cat("====================================================================\n")
cat("IBD SIGNATURE V2: CANDIDATE MARKER PSEUDOBULK BUILD\n")
cat("====================================================================\n")
cat("Grouping: dataset x biological_sample x coarse_identity\n")
cat("Historical Action: metadata only; no filtering and no grouping\n")


## ============================================================================
## 8. Load and standardize candidate universes
## ============================================================================

cat(
  "\n[1] Loading candidate metadata...\n"
)

gse_raw_meta <- as.data.frame(
  arrow::read_parquet(
    GSE_META
  )
)

scp1884_raw_meta <- as.data.frame(
  arrow::read_parquet(
    SCP1884_META
  )
)

scp259_raw_meta <- as.data.frame(
  arrow::read_parquet(
    SCP259_META
  )
)

stopifnot(
  nrow(
    gse_raw_meta
  ) == EXPECTED_CANDIDATES[[
    "GSE282122"
  ]],
  nrow(
    scp1884_raw_meta
  ) == EXPECTED_CANDIDATES[[
    "SCP1884"
  ]],
  nrow(
    scp259_raw_meta
  ) == EXPECTED_CANDIDATES[[
    "SCP259"
  ]]
)

gse_meta <- standardize_gse(
  gse_raw_meta
)

scp1884_meta <- standardize_scp1884(
  scp1884_raw_meta
)

scp259_meta <- standardize_scp259(
  scp259_raw_meta
)

rm(
  gse_raw_meta,
  scp1884_raw_meta,
  scp259_raw_meta
)

invisible(
  gc()
)

stopifnot(
  !anyDuplicated(
    gse_meta$cell_id
  ),
  !anyDuplicated(
    scp1884_meta$cell_id
  ),
  !anyDuplicated(
    scp259_meta$cell_id
  ),
  nrow(gse_meta) +
    nrow(scp1884_meta) +
    nrow(scp259_meta) ==
    EXPECTED_TOTAL
)

cat(
  "  GSE282122: ",
  format(
    nrow(gse_meta),
    big.mark = ","
  ),
  " cells\n",
  sep = ""
)

cat(
  "  SCP1884  : ",
  format(
    nrow(scp1884_meta),
    big.mark = ","
  ),
  " cells\n",
  sep = ""
)

cat(
  "  SCP259   : ",
  format(
    nrow(scp259_meta),
    big.mark = ","
  ),
  " cells\n",
  sep = ""
)

cat(
  "  TOTAL    : ",
  format(
    EXPECTED_TOTAL,
    big.mark = ","
  ),
  " cells\n",
  sep = ""
)


## ============================================================================
## 9. Build pseudobulk designs
## ============================================================================

cat(
  "\n[2] Building pseudobulk designs...\n"
)

cat(
  "\nGSE282122\n"
)

gse_design <- build_pb_design(
  gse_meta
)

gse_meta <- gse_design$cell_meta


cat(
  "\nSCP1884\n"
)

scp1884_design <- build_pb_design(
  scp1884_meta
)

scp1884_meta <- scp1884_design$cell_meta


cat(
  "\nSCP259\n"
)

scp259_design <- build_pb_design(
  scp259_meta
)

scp259_meta <- scp259_design$cell_meta


summary_rows <- list()


## ============================================================================
## 10. GSE282122 aggregation
## ============================================================================

cat("\n")
cat("====================================================================\n")
cat("[3] GSE282122\n")
cat("====================================================================\n")

sce <- zellkonverter::readH5AD(
  GSE_H5AD,
  reader = "R",
  use_hdf5 = TRUE,
  skip_assays = FALSE,
  raw = FALSE,
  layers = FALSE,
  verbose = TRUE
)

stopifnot(
  nrow(sce) == 33075L,
  ncol(sce) == 987743L,
  !anyDuplicated(
    colnames(sce)
  )
)

if (!"X" %in%
    SummarizedExperiment::assayNames(
      sce
    )) {

  fail(
    "GSE282122 expected assay X was not found."
  )
}


## --------------------------------------------------------------------------
## Gene symbols
## --------------------------------------------------------------------------

gse_gene <- as.data.frame(
  arrow::read_parquet(
    GSE_GENE_META
  )
)

stopifnot(
  "source_gene_id" %in%
    colnames(gse_gene),
  "gene_symbol" %in%
    colnames(gse_gene)
)

gene_match <- match(
  rownames(sce),
  gse_gene$source_gene_id
)

if (anyNA(
  gene_match
)) {

  fail(
    "GSE282122 source gene IDs failed gene-metadata matching."
  )
}

gse_symbols <- as.character(
  gse_gene$gene_symbol[
    gene_match
  ]
)

if (
  anyNA(
    gse_symbols
  ) ||
  any(
    gse_symbols == ""
  ) ||
  anyDuplicated(
    gse_symbols
  )
) {

  fail(
    "GSE282122 gene symbols are missing/blank/duplicated."
  )
}


## --------------------------------------------------------------------------
## Candidate cell mapping
## --------------------------------------------------------------------------

source_cells <- colnames(
  sce
)

candidate_idx <- match(
  gse_meta$cell_id,
  source_cells
)

if (
  anyNA(
    candidate_idx
  ) ||
  anyDuplicated(
    candidate_idx
  )
) {

  fail(
    "GSE282122 candidate cells failed exact source-cell matching."
  )
}

source_to_group <- integer(
  length(
    source_cells
  )
)

source_to_group[
  candidate_idx
] <- gse_meta$pb_index


## --------------------------------------------------------------------------
## Audit direct sparse H5AD X structure
## --------------------------------------------------------------------------

h5_structure <- rhdf5::h5ls(
  GSE_H5AD,
  recursive = TRUE
)

required_x_nodes <- c(
  "/X/data",
  "/X/indices",
  "/X/indptr"
)

observed_x_nodes <- paste0(
  h5_structure$group,
  "/",
  h5_structure$name
)

missing_x_nodes <- setdiff(
  required_x_nodes,
  observed_x_nodes
)

if (length(
  missing_x_nodes
) > 0L) {

  fail(
    "GSE282122 H5AD is missing required sparse-X nodes: ",
    paste(
      missing_x_nodes,
      collapse = ", "
    )
  )
}

indptr <- as.numeric(
  rhdf5::h5read(
    GSE_H5AD,
    "/X/indptr"
  )
)

## A cell-major CSR matrix has one indptr entry per source cell + 1.
if (
  length(indptr) !=
    ncol(sce) + 1L
) {

  fail(
    "GSE282122 /X/indptr length mismatch: observed ",
    length(indptr),
    ", expected ",
    ncol(sce) + 1L
  )
}


## --------------------------------------------------------------------------
## Aggregate GSE282122 CSR directly
## --------------------------------------------------------------------------

n_gse_pb <- nrow(
  gse_design$pb_meta
)

gse_pb <- Matrix::sparseMatrix(
  i = integer(0),
  j = integer(0),
  dims = c(
    nrow(sce),
    n_gse_pb
  ),
  giveCsparse = TRUE
)

gse_input_mass <- 0

starts <- seq.int(
  1L,
  ncol(sce),
  by = GSE_CELL_BLOCK
)

for (block_no in seq_along(
  starts
)) {

  first_cell <- starts[[
    block_no
  ]]

  last_cell <- min(
    first_cell +
      GSE_CELL_BLOCK -
      1L,
    ncol(sce)
  )

  block_ptr <- indptr[
    first_cell:
      (last_cell + 1L)
  ]

  nnz_per_cell <- diff(
    block_ptr
  )

  start0 <- block_ptr[[
    1L
  ]]

  end0 <- block_ptr[[
    length(
      block_ptr
    )
  ]]

  if (end0 > start0) {

    h5_index <- seq.int(
      start0 + 1,
      end0
    )

    gene0 <- as.integer(
      rhdf5::h5read(
        GSE_H5AD,
        "/X/indices",
        index = list(
          h5_index
        )
      )
    )

    block_x <- as.numeric(
      rhdf5::h5read(
        GSE_H5AD,
        "/X/data",
        index = list(
          h5_index
        )
      )
    )

    entry_group <- rep.int(
      source_to_group[
        first_cell:
          last_cell
      ],
      times = as.integer(
        nnz_per_cell
      )
    )

    stopifnot(
      length(gene0) ==
        length(block_x),
      length(block_x) ==
        length(entry_group)
    )

    keep <- entry_group > 0L

    if (any(
      keep
    )) {

      selected_x <- block_x[
        keep
      ]

      gse_input_mass <- gse_input_mass +
        sum(
          selected_x
        )

      partial <- Matrix::sparseMatrix(
        i = gene0[
          keep
        ] + 1L,
        j = entry_group[
          keep
        ],
        x = selected_x,
        dims = c(
          nrow(sce),
          n_gse_pb
        ),
        giveCsparse = TRUE
      )

      gse_pb <- gse_pb +
        partial

      rm(
        partial,
        selected_x
      )
    }

    rm(
      gene0,
      block_x,
      entry_group,
      keep
    )
  }

  if (
    block_no %% 25L == 0L ||
    block_no == length(starts)
  ) {

    cat(
      "  H5AD cell block ",
      block_no,
      "/",
      length(starts),
      "\n",
      sep = ""
    )
  }

  invisible(
    gc(FALSE)
  )
}

if (
  abs(
    sum(gse_pb) -
      gse_input_mass
  ) >
    1e-6
) {

  fail(
    "GSE282122 raw-count mass closure failed."
  )
}

summary_rows[[
  "GSE282122"
]] <- save_dataset(
  dataset = "GSE282122",
  counts = gse_pb,
  pb_meta = gse_design$pb_meta,
  gene_symbols = gse_symbols,
  input_mass = gse_input_mass
)

rm(
  sce,
  gse_pb,
  indptr,
  source_to_group,
  source_cells,
  gse_gene
)

invisible(
  gc()
)


## ============================================================================
## 11. SCP1884 aggregation
## ============================================================================

cat("\n")
cat("====================================================================\n")
cat("[4] SCP1884\n")
cat("====================================================================\n")

scp1884_config <- data.frame(
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
  stringsAsFactors = FALSE
)

colon_features <- read_scp1884_features(
  file.path(
    SCP1884_EXPR,
    "Colon/Stromal",
    "CO_STR.features.tsv"
  )
)

ti_features <- read_scp1884_features(
  file.path(
    SCP1884_EXPR,
    "TerminalIleum/Stromal",
    "TI_STR.features.tsv"
  )
)

common_ens <- colon_features$ensembl_id[
  colon_features$ensembl_id %in%
    ti_features$ensembl_id
]

if (
  length(
    common_ens
  ) !=
    27830L
) {

  fail(
    "SCP1884 expected 27,830 common genes; observed ",
    length(common_ens)
  )
}

colon_idx <- match(
  common_ens,
  colon_features$ensembl_id
)

ti_idx <- match(
  common_ens,
  ti_features$ensembl_id
)

if (
  anyNA(
    colon_idx
  ) ||
  anyNA(
    ti_idx
  )
) {

  fail(
    "SCP1884 common-gene matching failed."
  )
}

if (!identical(
  colon_features$gene_symbol[
    colon_idx
  ],
  ti_features$gene_symbol[
    ti_idx
  ]
)) {

  fail(
    "SCP1884 Colon/TI gene-symbol alignment failed."
  )
}

scp1884_symbols <- as.character(
  colon_features$gene_symbol[
    colon_idx
  ]
)

if (
  anyNA(
    scp1884_symbols
  ) ||
  any(
    scp1884_symbols == ""
  ) ||
  anyDuplicated(
    scp1884_symbols
  )
) {

  fail(
    "SCP1884 canonical gene symbols are invalid."
  )
}

n_scp1884_genes <- length(
  scp1884_symbols
)

n_scp1884_pb <- nrow(
  scp1884_design$pb_meta
)

scp1884_pb <- Matrix::sparseMatrix(
  i = integer(0),
  j = integer(0),
  dims = c(
    n_scp1884_genes,
    n_scp1884_pb
  ),
  giveCsparse = TRUE
)

scp1884_input_mass <- 0

for (ii in seq_len(
  nrow(
    scp1884_config
  )
)) {

  bundle <- scp1884_config$source_bundle[
    ii
  ]

  subdir <- scp1884_config$subdir[
    ii
  ]

  cat(
    "\n  Bundle: ",
    bundle,
    "\n",
    sep = ""
  )

  raw_file <- file.path(
    SCP1884_EXPR,
    subdir,
    "raw.mtx"
  )

  barcode_file <- file.path(
    SCP1884_EXPR,
    subdir,
    paste0(
      bundle,
      ".barcodes.tsv"
    )
  )

  feature_file <- file.path(
    SCP1884_EXPR,
    subdir,
    paste0(
      bundle,
      ".features.tsv"
    )
  )

  if (!all(
    file.exists(
      c(
        raw_file,
        barcode_file,
        feature_file
      )
    )
  )) {

    fail(
      "Missing SCP1884 source files for bundle ",
      bundle
    )
  }

  features <- read_scp1884_features(
    feature_file
  )

  gene_idx <- match(
    common_ens,
    features$ensembl_id
  )

  if (anyNA(
    gene_idx
  )) {

    fail(
      "SCP1884 canonical gene matching failed in ",
      bundle
    )
  }

  if (!identical(
    as.character(
      features$gene_symbol[
        gene_idx
      ]
    ),
    scp1884_symbols
  )) {

    fail(
      "SCP1884 gene-symbol order mismatch in ",
      bundle
    )
  }

  source_gene_to_output <- integer(
    nrow(
      features
    )
  )

  source_gene_to_output[
    gene_idx
  ] <- seq_along(
    common_ens
  )

  meta_sub <- scp1884_meta[
    scp1884_meta$source_bundle ==
      bundle,
    ,
    drop = FALSE
  ]

  if (
    nrow(
      meta_sub
    ) ==
      0L
  ) {

    fail(
      "No SCP1884 candidate metadata for bundle ",
      bundle
    )
  }

  agg <- aggregate_triplet_bundle(
    raw_file = raw_file,
    barcode_file = barcode_file,
    source_gene_to_output = source_gene_to_output,
    meta = meta_sub,
    n_output_genes = n_scp1884_genes,
    n_output_groups = n_scp1884_pb,
    require_all_source_cells = TRUE
  )

  scp1884_pb <- scp1884_pb +
    agg$counts

  scp1884_input_mass <- scp1884_input_mass +
    agg$input_mass

  rm(
    agg,
    features,
    gene_idx,
    source_gene_to_output,
    meta_sub
  )

  invisible(
    gc()
  )
}

if (
  abs(
    sum(scp1884_pb) -
      scp1884_input_mass
  ) >
    1e-6
) {

  fail(
    "SCP1884 dataset-level raw-count mass closure failed."
  )
}

summary_rows[[
  "SCP1884"
]] <- save_dataset(
  dataset = "SCP1884",
  counts = scp1884_pb,
  pb_meta = scp1884_design$pb_meta,
  gene_symbols = scp1884_symbols,
  input_mass = scp1884_input_mass
)

rm(
  scp1884_pb,
  colon_features,
  ti_features
)

invisible(
  gc()
)


## ============================================================================
## 12. SCP259 aggregation
## ============================================================================

cat("\n")
cat("====================================================================\n")
cat("[5] SCP259\n")
cat("====================================================================\n")

scp259_config <- data.frame(
  source_bundle = c(
    "Epithelial",
    "Stromal",
    "Immune"
  ),
  prefix = c(
    "Epi",
    "Fib",
    "Imm"
  ),
  source_hash = c(
    "5cdc540d328cee7a2efc2348",
    "5cdc540d328cee7a2efc2349",
    "5cdc540d328cee7a2efc234a"
  ),
  stringsAsFactors = FALSE
)

scp259_config$source_dir <- file.path(
  SCP259_EXPR,
  scp259_config$source_hash
)

scp259_config$raw_file <- file.path(
  scp259_config$source_dir,
  paste0(
    "gene_sorted-",
    scp259_config$prefix,
    ".matrix.mtx"
  )
)

scp259_config$barcode_file <- file.path(
  scp259_config$source_dir,
  paste0(
    scp259_config$prefix,
    ".barcodes2.tsv"
  )
)

scp259_config$gene_file <- file.path(
  scp259_config$source_dir,
  paste0(
    scp259_config$prefix,
    ".genes.tsv"
  )
)

all_scp259_sources <- c(
  scp259_config$raw_file,
  scp259_config$barcode_file,
  scp259_config$gene_file
)

if (!all(
  file.exists(
    all_scp259_sources
  )
)) {

  fail(
    "Missing SCP259 source files:\n",
    paste(
      all_scp259_sources[
        !file.exists(
          all_scp259_sources
        )
      ],
      collapse = "\n"
    )
  )
}

epi_gene <- read_scp259_genes(
  scp259_config$gene_file[
    scp259_config$prefix ==
      "Epi"
  ]
)

fib_gene <- read_scp259_genes(
  scp259_config$gene_file[
    scp259_config$prefix ==
      "Fib"
  ]
)

imm_gene <- read_scp259_genes(
  scp259_config$gene_file[
    scp259_config$prefix ==
      "Imm"
  ]
)

scp259_symbols <- epi_gene$gene_symbol[
  epi_gene$gene_symbol %in%
    fib_gene$gene_symbol &
    epi_gene$gene_symbol %in%
    imm_gene$gene_symbol
]

if (
  length(
    scp259_symbols
  ) !=
    18172L
) {

  fail(
    "SCP259 expected 18,172 common genes; observed ",
    length(
      scp259_symbols
    )
  )
}

scp259_symbols <- as.character(
  scp259_symbols
)

if (
  anyNA(
    scp259_symbols
  ) ||
  any(
    scp259_symbols == ""
  ) ||
  anyDuplicated(
    scp259_symbols
  )
) {

  fail(
    "SCP259 canonical gene symbols are invalid."
  )
}

n_scp259_genes <- length(
  scp259_symbols
)

n_scp259_pb <- nrow(
  scp259_design$pb_meta
)

scp259_pb <- Matrix::sparseMatrix(
  i = integer(0),
  j = integer(0),
  dims = c(
    n_scp259_genes,
    n_scp259_pb
  ),
  giveCsparse = TRUE
)

scp259_input_mass <- 0

for (ii in seq_len(
  nrow(
    scp259_config
  )
)) {

  bundle <- scp259_config$source_bundle[
    ii
  ]

  cat(
    "\n  Bundle: ",
    bundle,
    "\n",
    sep = ""
  )

  genes <- read_scp259_genes(
    scp259_config$gene_file[
      ii
    ]
  )

  gene_idx <- match(
    scp259_symbols,
    genes$gene_symbol
  )

  if (anyNA(
    gene_idx
  )) {

    fail(
      "SCP259 canonical gene matching failed in ",
      bundle
    )
  }

  if (!identical(
    as.character(
      genes$gene_symbol[
        gene_idx
      ]
    ),
    scp259_symbols
  )) {

    fail(
      "SCP259 canonical gene order mismatch in ",
      bundle
    )
  }

  source_gene_to_output <- integer(
    nrow(
      genes
    )
  )

  source_gene_to_output[
    gene_idx
  ] <- seq_along(
    scp259_symbols
  )

  meta_sub <- scp259_meta[
    scp259_meta$source_bundle ==
      bundle,
    ,
    drop = FALSE
  ]

  if (
    nrow(
      meta_sub
    ) ==
      0L
  ) {

    fail(
      "No SCP259 candidate metadata for bundle ",
      bundle
    )
  }

  agg <- aggregate_triplet_bundle(
    raw_file = scp259_config$raw_file[
      ii
    ],
    barcode_file = scp259_config$barcode_file[
      ii
    ],
    source_gene_to_output = source_gene_to_output,
    meta = meta_sub,
    n_output_genes = n_scp259_genes,
    n_output_groups = n_scp259_pb,
    require_all_source_cells = TRUE
  )

  scp259_pb <- scp259_pb +
    agg$counts

  scp259_input_mass <- scp259_input_mass +
    agg$input_mass

  rm(
    agg,
    genes,
    gene_idx,
    source_gene_to_output,
    meta_sub
  )

  invisible(
    gc()
  )
}

if (
  abs(
    sum(scp259_pb) -
      scp259_input_mass
  ) >
    1e-6
) {

  fail(
    "SCP259 dataset-level raw-count mass closure failed."
  )
}

summary_rows[[
  "SCP259"
]] <- save_dataset(
  dataset = "SCP259",
  counts = scp259_pb,
  pb_meta = scp259_design$pb_meta,
  gene_symbols = scp259_symbols,
  input_mass = scp259_input_mass
)

rm(
  scp259_pb,
  epi_gene,
  fib_gene,
  imm_gene
)

invisible(
  gc()
)


## ============================================================================
## 13. Cross-dataset metadata and coverage
## ============================================================================

cat("\n")
cat("====================================================================\n")
cat("[6] Cross-dataset audit outputs\n")
cat("====================================================================\n")

all_pb_meta <- rbind(
  gse_design$pb_meta,
  scp1884_design$pb_meta,
  scp259_design$pb_meta
)

if (
  sum(
    all_pb_meta$n_cells
  ) !=
    EXPECTED_TOTAL
) {

  fail(
    "Combined pseudobulk metadata does not close to full candidate universe."
  )
}

if (anyDuplicated(
  all_pb_meta$pb_id
)) {

  fail(
    "Combined pseudobulk IDs are not globally unique."
  )
}

utils::write.table(
  all_pb_meta,
  file.path(
    OUTDIR,
    "all_datasets_pseudobulk_metadata.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

arrow::write_parquet(
  all_pb_meta,
  file.path(
    OUTDIR,
    "all_datasets_pseudobulk_metadata.parquet"
  )
)


## --------------------------------------------------------------------------
## Dataset x coarse-identity coverage
## --------------------------------------------------------------------------

coverage_keys <- unique(
  all_pb_meta[
    c(
      "dataset",
      "coarse_identity"
    )
  ]
)

coverage_rows <- vector(
  "list",
  nrow(
    coverage_keys
  )
)

for (ii in seq_len(
  nrow(
    coverage_keys
  )
)) {

  ds <- coverage_keys$dataset[
    ii
  ]

  ct <- coverage_keys$coarse_identity[
    ii
  ]

  sub <- all_pb_meta[
    all_pb_meta$dataset ==
      ds &
      all_pb_meta$coarse_identity ==
        ct,
    ,
    drop = FALSE
  ]

  coverage_rows[[
    ii
  ]] <- data.frame(
    dataset = ds,
    coarse_identity = ct,
    n_cells = sum(
      sub$n_cells
    ),
    n_pseudobulks = nrow(
      sub
    ),
    n_samples = length(
      unique(
        sub$sample_id
      )
    ),
    n_donors = length(
      unique(
        sub$donor_id
      )
    ),
    n_sites = length(
      unique(
        sub$site
      )
    ),
    sites = collapse_unique(
      sub$site
    ),
    diseases = collapse_unique(
      sub$disease
    ),
    inflammation_states = collapse_unique(
      sub$inflammation
    ),
    min_cells_per_pseudobulk = min(
      sub$n_cells
    ),
    median_cells_per_pseudobulk = stats::median(
      sub$n_cells
    ),
    max_cells_per_pseudobulk = max(
      sub$n_cells
    ),
    stringsAsFactors = FALSE
  )
}

coverage_out <- do.call(
  rbind,
  coverage_rows
)

coverage_out <- coverage_out[
  order(
    coverage_out$dataset,
    -coverage_out$n_cells,
    coverage_out$coarse_identity
  ),
  ,
  drop = FALSE
]

utils::write.table(
  coverage_out,
  file.path(
    OUTDIR,
    "pseudobulk_coarse_identity_coverage.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## --------------------------------------------------------------------------
## Dataset summary
## --------------------------------------------------------------------------

summary_table <- do.call(
  rbind,
  summary_rows
)

utils::write.table(
  summary_table,
  file.path(
    OUTDIR,
    "pseudobulk_build_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## --------------------------------------------------------------------------
## Final integrity
## --------------------------------------------------------------------------

stopifnot(
  sum(
    summary_table$n_candidate_cells
  ) ==
    EXPECTED_TOTAL,
  sum(
    all_pb_meta$n_cells
  ) ==
    EXPECTED_TOTAL,
  all(
    summary_table$mass_closure
  ),
  all(
    summary_table$postwrite_random_concordance
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    "pseudobulk_build_sessionInfo.txt"
  )
)


## ============================================================================
## 14. Final report
## ============================================================================

cat("\n")
cat("====================================================================\n")
cat("CANDIDATE MARKER PSEUDOBULK BUILD COMPLETE\n")
cat("====================================================================\n")

print(
  summary_table,
  row.names = FALSE
)

cat(
  "\nCandidate cells represented: ",
  format(
    sum(
      summary_table$n_candidate_cells
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "Total pseudobulks: ",
  format(
    nrow(
      all_pb_meta
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "Grouping: Dataset x BiologicalSample x CoarseIdentity\n"
)

cat(
  "Historical Action: metadata only; not used for filtering/grouping\n"
)

cat(
  "Output directory: ",
  OUTDIR,
  "\n",
  sep = ""
)

cat(
  "FINAL STATUS: PASS\n"
)
