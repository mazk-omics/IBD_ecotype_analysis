#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# ==============================================================================
# IBD EcoTyper
# Stage 2 epithelial BayesPrism
# QC-04: Tuft reference source-balance audit
#
# Questions:
# 1. Which source datasets contribute Tuft cells?
# 2. Is Tuft reference dominated by one dataset?
# 3. Are dataset-specific Tuft pseudo-bulk profiles concordant?
# 4. Are strong Tuft-specific genes consistently expressed across datasets?
#
# NO BayesPrism fitting.
# NO modification of reference.
# ==============================================================================


PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

REF_RDS <- file.path(
  PROJECT_ROOT,
  "03_reference/combined_reference/output/11_bayesprism_full_input",
  "bayesprism_full_reference_count_matrix.rds"
)

MANIFEST_PARQUET <- file.path(
  PROJECT_ROOT,
  "03_reference/combined_reference/output/05_final_cell_selection",
  "selected_cell_manifest.parquet"
)

QC_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference/combined_reference/output/21_bayesprism_epithelial_qc"
)

STRONG_TUFT_TSV <- file.path(
  QC_DIR,
  "stage2_tuft_strong_specific_genes.tsv"
)

dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


OUT_COUNTS <- file.path(
  QC_DIR,
  "tuft_reference_cell_counts_by_dataset.tsv"
)

OUT_PROFILE <- file.path(
  QC_DIR,
  "tuft_reference_profile_similarity_by_dataset.tsv"
)

OUT_MARKERS <- file.path(
  QC_DIR,
  "tuft_reference_marker_expression_by_dataset.tsv"
)

OUT_DOMINANCE <- file.path(
  QC_DIR,
  "tuft_specific_gene_dataset_dominance.tsv"
)

OUT_LOG <- file.path(
  QC_DIR,
  "tuft_reference_source_balance_audit.log"
)

OUT_SESSION <- file.path(
  QC_DIR,
  "tuft_reference_source_balance_sessionInfo.txt"
)


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

  while (sink.number(type = "message") > 0)
    sink(type = "message")

  while (sink.number(type = "output") > 0)
    sink(type = "output")

  close(log_con)

}, add = TRUE)


assert_true <- function(x, msg) {

  if (!isTRUE(x))
    stop(msg, call. = FALSE)
}


write_tsv <- function(x, path) {

  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


cat("============================================================\n")
cat("Tuft reference source-balance audit\n")
cat("============================================================\n\n")

cat("START:", format(Sys.time()), "\n\n")


# ==============================================================================
# Packages
# ==============================================================================

suppressPackageStartupMessages(
  library(Matrix)
)

suppressPackageStartupMessages(
  library(arrow)
)


# ==============================================================================
# Load formal reference
# ==============================================================================

assert_true(
  file.exists(REF_RDS),
  paste("Reference missing:", REF_RDS)
)

cat("Loading unified reference...\n")

ref_obj <- readRDS(
  REF_RDS
)


assert_true(
  is.list(ref_obj),
  "Reference RDS is not expected list."
)

assert_true(
  "reference" %in% names(ref_obj),
  "reference slot absent."
)

X <- ref_obj$reference

labels <- ref_obj$cell.type.labels


assert_true(
  inherits(X, "Matrix") ||
    is.matrix(X),
  "Reference matrix has unexpected class."
)

assert_true(
  nrow(X) == length(labels),
  "Reference labels do not match cell count."
)


if (is.null(names(labels))) {

  names(labels) <- rownames(X)
}


cat(
  "[PASS] Reference:",
  nrow(X),
  "cells x",
  ncol(X),
  "genes\n\n"
)


# ==============================================================================
# Load selected-cell manifest
# ==============================================================================

assert_true(
  file.exists(MANIFEST_PARQUET),
  paste("Manifest missing:", MANIFEST_PARQUET)
)

manifest <- as.data.frame(
  arrow::read_parquet(
    MANIFEST_PARQUET
  )
)


cat(
  "[PASS] Manifest:",
  nrow(manifest),
  "rows x",
  ncol(manifest),
  "columns\n"
)

cat(
  "Manifest columns:\n",
  paste(colnames(manifest), collapse = ", "),
  "\n\n"
)


# ==============================================================================
# Identify manifest alignment key
# ==============================================================================

dataset_col <- "dataset"

assert_true(
  dataset_col %in% colnames(manifest),
  "dataset column missing from manifest."
)


# Candidate cell-ID columns, ordered by preference.
# For a unified multi-dataset reference, global_cell_id should normally
# be the correct primary key.

cell_candidates <- c(
  "global_cell_id",
  "source_cell_key",
  "cell_id"
)


cell_candidates <- cell_candidates[
  cell_candidates %in% colnames(manifest)
]


assert_true(
  length(cell_candidates) > 0,
  "No candidate cell-ID columns found in manifest."
)


reference_ids <- rownames(X)


assert_true(
  !is.null(reference_ids),
  "Reference matrix has no rownames."
)


cat("Testing manifest alignment keys...\n\n")


alignment_test <- do.call(
  rbind,
  lapply(
    cell_candidates,
    function(nm) {

      matched <- reference_ids %in%
        manifest[[nm]]

      data.frame(
        candidate = nm,
        matched_reference_cells = sum(matched),
        total_reference_cells = length(reference_ids),
        match_fraction = mean(matched),
        stringsAsFactors = FALSE
      )
    }
  )
)


print(
  alignment_test,
  row.names = FALSE
)


best_row <- which.max(
  alignment_test$matched_reference_cells
)


cell_col <- alignment_test$candidate[
  best_row
]


best_match_n <- alignment_test$matched_reference_cells[
  best_row
]


cat(
  "\nSelected alignment column:",
  cell_col,
  "\n"
)

cat(
  "Matched cells:",
  best_match_n,
  "/",
  length(reference_ids),
  "\n\n"
)


assert_true(
  best_match_n == length(reference_ids),
  paste(
    "No manifest column provides complete alignment.",
    "Best candidate:",
    cell_col,
    "matched",
    best_match_n,
    "of",
    length(reference_ids),
    "reference cells."
  )
)


# ==============================================================================
# Align manifest to formal reference
# ==============================================================================

idx <- match(
  reference_ids,
  manifest[[cell_col]]
)


assert_true(
  !anyNA(idx),
  paste(
    sum(is.na(idx)),
    "reference cells absent from manifest after best-key alignment."
  )
)


# Check uniqueness of the selected key.
# This matters because match() silently uses the first duplicated entry.

selected_ids <- manifest[[cell_col]][idx]


assert_true(
  !anyDuplicated(selected_ids),
  paste(
    "Selected manifest key is not unique:",
    cell_col
  )
)


dataset <- as.character(
  manifest[[dataset_col]][idx]
)


assert_true(
  !anyNA(dataset),
  "Missing dataset labels after alignment."
)


cat(
  "[PASS] Manifest aligned using:",
  cell_col,
  "\n"
)

cat(
  "[PASS] Aligned cells:",
  length(idx),
  "\n\n"
)


# ==============================================================================
# Align manifest to formal reference
# ==============================================================================

idx <- match(
  rownames(X),
  manifest[[cell_col]]
)


assert_true(
  !anyNA(idx),
  paste(
    sum(is.na(idx)),
    "reference cells absent from manifest."
  )
)


dataset <- as.character(
  manifest[[dataset_col]][idx]
)


assert_true(
  !anyNA(dataset),
  "Missing dataset labels after alignment."
)


# ==============================================================================
# Epithelial reference counts by dataset
# ==============================================================================

epi_types <- c(
  "Colonic_absorptive",
  "Enteroendocrine",
  "Goblet",
  "Ileal_absorptive",
  "Paneth",
  "Stem_TA_progenitor",
  "Tuft"
)


epi_idx <- labels %in% epi_types


count_table <- as.data.frame.matrix(
  table(
    dataset[epi_idx],
    labels[epi_idx]
  )
)


count_table$dataset <- rownames(
  count_table
)

rownames(count_table) <- NULL


for (nm in setdiff(
  epi_types,
  colnames(count_table)
)) {

  count_table[[nm]] <- 0
}


count_table$total_epithelial <- rowSums(
  count_table[, epi_types, drop = FALSE]
)

count_table$Tuft_fraction_within_reference_epithelium <-
  count_table$Tuft /
  count_table$total_epithelial


total_tuft <- sum(
  labels == "Tuft"
)

count_table$Tuft_fraction_of_all_reference_Tuft <-
  count_table$Tuft /
  total_tuft


count_table <- count_table[
  order(
    -count_table$Tuft
  ),
  ,
  drop = FALSE
]


write_tsv(
  count_table,
  OUT_COUNTS
)


cat("============================================================\n")
cat("EPITHELIAL REFERENCE COUNTS BY DATASET\n")
cat("============================================================\n\n")

print(
  count_table[
    ,
    c(
      "dataset",
      epi_types,
      "total_epithelial",
      "Tuft_fraction_within_reference_epithelium",
      "Tuft_fraction_of_all_reference_Tuft"
    )
  ],
  row.names = FALSE
)


# ==============================================================================
# Tuft cells only
# ==============================================================================

tuft_idx <- which(
  labels == "Tuft"
)


assert_true(
  length(tuft_idx) > 0,
  "No Tuft cells found."
)


cat(
  "\nTotal Tuft cells:",
  length(tuft_idx),
  "\n\n"
)


tuft_datasets <- sort(
  unique(
    dataset[tuft_idx]
  )
)


# ==============================================================================
# Global Tuft pseudobulk
# ==============================================================================

global_counts <- Matrix::colSums(
  X[
    tuft_idx,
    ,
    drop = FALSE
  ]
)


global_profile <- global_counts /
  sum(global_counts)


# ==============================================================================
# Canonical markers
# ==============================================================================

canonical <- c(
  "POU2F3",
  "PLCB2",
  "GFI1B",
  "AVIL",
  "PTGS1",
  "SH2D6",
  "DCLK1"
)


canonical_present <- intersect(
  canonical,
  colnames(X)
)


# ==============================================================================
# Strong Tuft-specific genes from QC-03
# ==============================================================================

strong_genes <- character()


if (file.exists(STRONG_TUFT_TSV)) {

  strong_tbl <- read.delim(
    STRONG_TUFT_TSV,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  strong_genes <- intersect(
    strong_tbl$gene,
    colnames(X)
  )
}


cat(
  "Strong Tuft-specific genes available:",
  length(strong_genes),
  "\n\n"
)


# ==============================================================================
# Dataset-specific Tuft pseudobulks
# ==============================================================================

profile_rows <- list()
marker_rows <- list()

dataset_gene_counts <- matrix(
  0,
  nrow = length(tuft_datasets),
  ncol = length(strong_genes),
  dimnames = list(
    tuft_datasets,
    strong_genes
  )
)


for (ds in tuft_datasets) {

  current_idx <- tuft_idx[
    dataset[tuft_idx] == ds
  ]


  counts <- Matrix::colSums(
    X[
      current_idx,
      ,
      drop = FALSE
    ]
  )


  total_counts <- sum(
    counts
  )


  profile <- counts /
    total_counts


  pearson <- suppressWarnings(
    cor(
      global_profile,
      profile,
      method = "pearson"
    )
  )


  spearman <- suppressWarnings(
    cor(
      global_profile,
      profile,
      method = "spearman"
    )
  )


  cosine <- sum(
    global_profile * profile
  ) /
    sqrt(
      sum(global_profile^2) *
        sum(profile^2)
    )


  strong_mass_fraction <- NA_real_

  if (length(strong_genes) > 0) {

    strong_mass_fraction <- sum(
      profile[
        strong_genes
      ]
    )

    dataset_gene_counts[
      ds,
      strong_genes
    ] <- counts[
      strong_genes
    ]
  }


  profile_rows[[ds]] <- data.frame(

    dataset =
      ds,

    Tuft_cells =
      length(current_idx),

    total_raw_counts =
      total_counts,

    pearson_vs_global_Tuft =
      pearson,

    spearman_vs_global_Tuft =
      spearman,

    cosine_vs_global_Tuft =
      cosine,

    strong_Tuft_gene_mass_fraction =
      strong_mass_fraction,

    stringsAsFactors = FALSE
  )


  if (length(canonical_present) > 0) {

    marker_rows[[ds]] <- data.frame(

      dataset =
        ds,

      gene =
        canonical_present,

      normalized_CPM =
        profile[
          canonical_present
        ] * 1e6,

      raw_counts =
        counts[
          canonical_present
        ],

      stringsAsFactors = FALSE
    )
  }
}


profile_table <- do.call(
  rbind,
  profile_rows
)

rownames(profile_table) <- NULL


profile_table <- profile_table[
  order(
    -profile_table$Tuft_cells
  ),
  ,
  drop = FALSE
]


write_tsv(
  profile_table,
  OUT_PROFILE
)


marker_table <- do.call(
  rbind,
  marker_rows
)

rownames(marker_table) <- NULL


write_tsv(
  marker_table,
  OUT_MARKERS
)


# ==============================================================================
# Dataset dominance of strong Tuft genes
# ==============================================================================

if (length(strong_genes) > 0) {

  gene_totals <- colSums(
    dataset_gene_counts
  )


  dominant_dataset <- apply(
    dataset_gene_counts,
    2,
    function(x) {

      if (sum(x) == 0)
        return(NA_character_)

      names(x)[
        which.max(x)
      ]
    }
  )


  max_dataset_fraction <- apply(
    dataset_gene_counts,
    2,
    function(x) {

      if (sum(x) == 0)
        return(NA_real_)

      max(x) /
        sum(x)
    }
  )


  dominance_table <- data.frame(

    gene =
      strong_genes,

    total_Tuft_counts =
      gene_totals,

    dominant_dataset =
      dominant_dataset,

    max_dataset_fraction =
      max_dataset_fraction,

    stringsAsFactors = FALSE
  )


  dominance_table <- dominance_table[
    order(
      -dominance_table$max_dataset_fraction
    ),
    ,
    drop = FALSE
  ]


  write_tsv(
    dominance_table,
    OUT_DOMINANCE
  )

} else {

  dominance_table <- NULL
}


# ==============================================================================
# Console summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("TUFT PSEUDOBULK SIMILARITY BY DATASET\n")
cat("============================================================\n\n")

print(
  profile_table,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("CANONICAL TUFT MARKERS BY DATASET\n")
cat("============================================================\n\n")

print(
  marker_table,
  row.names = FALSE
)


if (!is.null(dominance_table)) {

  cat("\n")
  cat("============================================================\n")
  cat("STRONG TUFT-GENE DATASET DOMINANCE\n")
  cat("============================================================\n\n")

  cat(
    "Median max-dataset contribution:",
    median(
      dominance_table$max_dataset_fraction,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "P95 max-dataset contribution:",
    quantile(
      dominance_table$max_dataset_fraction,
      0.95,
      na.rm = TRUE
    ),
    "\n\n"
  )


  cat(
    "Most dataset-dominated strong genes:\n\n"
  )

  print(
    head(
      dominance_table,
      30
    ),
    row.names = FALSE
  )
}


capture.output(
  sessionInfo(),
  file = OUT_SESSION
)


cat("\nOutput directory:\n")
cat(QC_DIR, "\n")

cat(
  "\nEND:",
  format(Sys.time()),
  "\n"
)