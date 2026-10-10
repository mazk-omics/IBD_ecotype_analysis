#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Paths
# ============================================================

ROOT <- "/home/mazekai/IBD_EcoTyper"

INPUT_ROOT <- file.path(
  ROOT,
  "05_EcoTyper/input/MSCCR_CORE12"
)

SOURCE_MATRIX_DIR <- file.path(
  INPUT_ROOT,
  "expression_matrices"
)

SOURCE_ANNOTATION <- file.path(
  INPUT_ROOT,
  "annotation/MSCCR_CORE12_annotation.tsv"
)

OUT_IBD <- file.path(
  INPUT_ROOT,
  "region_ibd"
)

OUT_ACTIVE <- file.path(
  INPUT_ROOT,
  "region_activeibd"
)

# ============================================================
# Locked CORE12
# ============================================================

CORE12 <- c(
  "B",
  "CD4_T",
  "CD8_T",
  "Endothelial",
  "Enteroendocrine",
  "Fibroblast",
  "Glial",
  "Goblet",
  "ILC",
  "Mast",
  "Monocyte_Macrophage",
  "Pericyte"
)

EXPECTED_ALL <- 2486L
EXPECTED_IBD <- 2025L
EXPECTED_ACTIVE <- 708L

# ============================================================
# Helpers
# ============================================================

fail <- function(...) {
  stop(
    paste0(...),
    call. = FALSE
  )
}

section <- function(x) {
  cat(
    "\n============================================================\n",
    x,
    "\n============================================================\n",
    sep = ""
  )
}

# ============================================================
# 1. Read and audit baseline annotation
# ============================================================

section("1. BASELINE ANNOTATION")

if (!file.exists(SOURCE_ANNOTATION)) {
  fail(
    "Missing annotation: ",
    SOURCE_ANNOTATION
  )
}

ann <- fread(
  SOURCE_ANNOTATION,
  data.table = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "ID",
  "Disease",
  "Inflammation",
  "Region3"
)

missing_cols <- setdiff(
  required_cols,
  colnames(ann)
)

if (length(missing_cols) > 0) {
  fail(
    "Annotation is missing required columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

if (nrow(ann) != EXPECTED_ALL) {
  fail(
    "Expected ",
    EXPECTED_ALL,
    " annotation rows, found ",
    nrow(ann)
  )
}

if (anyDuplicated(ann$ID)) {
  fail("Duplicated sample IDs detected in annotation.")
}

if (anyNA(ann$ID) || any(ann$ID == "")) {
  fail("Missing/empty sample IDs detected.")
}

if (anyNA(ann$Disease)) {
  fail("Missing Disease values detected.")
}

if (anyNA(ann$Inflammation)) {
  fail("Missing Inflammation values detected.")
}

if (anyNA(ann$Region3)) {
  fail("Missing Region3 values detected.")
}

cat("Baseline samples:", nrow(ann), "\n\n")

cat("Disease:\n")
print(table(ann$Disease))

cat("\nInflammation:\n")
print(table(ann$Inflammation))

cat("\nRegion3:\n")
print(table(ann$Region3))

# ============================================================
# 2. Define C and D cohorts
# ============================================================

section("2. DEFINE SENSITIVITY COHORTS")

# C:
# CD + UC, regardless of inflammation status
ann_ibd <- ann[
  ann$Disease %in% c("CD", "UC"),
  ,
  drop = FALSE
]

# D:
# CD + UC AND Inflamed
ann_active <- ann[
  ann$Disease %in% c("CD", "UC") &
    ann$Inflammation == "Inflamed",
  ,
  drop = FALSE
]

if (nrow(ann_ibd) != EXPECTED_IBD) {
  fail(
    "IBD-only cohort expected ",
    EXPECTED_IBD,
    " samples, found ",
    nrow(ann_ibd)
  )
}

if (nrow(ann_active) != EXPECTED_ACTIVE) {
  fail(
    "Active-IBD cohort expected ",
    EXPECTED_ACTIVE,
    " samples, found ",
    nrow(ann_active)
  )
}

if (!all(ann_active$ID %in% ann_ibd$ID)) {
  fail(
    "Active-IBD sample set is not a subset of IBD-only cohort."
  )
}

cat("\nC = IBD-only\n")
cat("Samples:", nrow(ann_ibd), "\n")
print(
  table(
    ann_ibd$Disease,
    ann_ibd$Region3
  )
)

cat("\nD = Active IBD\n")
cat("Samples:", nrow(ann_active), "\n")
print(
  table(
    ann_active$Disease,
    ann_active$Region3
  )
)

# Hard expectations from audited metadata
expected_ibd_region <- c(
  Colon = 799L,
  Ileum = 550L,
  Rectum = 676L
)

expected_active_region <- c(
  Colon = 290L,
  Ileum = 161L,
  Rectum = 257L
)

actual_ibd_region <- table(
  factor(
    ann_ibd$Region3,
    levels = names(expected_ibd_region)
  )
)

actual_active_region <- table(
  factor(
    ann_active$Region3,
    levels = names(expected_active_region)
  )
)

if (!identical(
  as.integer(actual_ibd_region),
  as.integer(expected_ibd_region)
)) {
  fail(
    "IBD-only Region3 counts do not match expected values."
  )
}

if (!identical(
  as.integer(actual_active_region),
  as.integer(expected_active_region)
)) {
  fail(
    "Active-IBD Region3 counts do not match expected values."
  )
}

cat("\n[PASS] Cohort definitions match expected counts.\n")

# ============================================================
# 3. Create output directories
# ============================================================

section("3. CREATE OUTPUT DIRECTORIES")

dirs <- c(
  file.path(
    OUT_IBD,
    "expression_matrices"
  ),
  file.path(
    OUT_IBD,
    "annotation"
  ),
  file.path(
    OUT_ACTIVE,
    "expression_matrices"
  ),
  file.path(
    OUT_ACTIVE,
    "annotation"
  )
)

for (d in dirs) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cat(
  "IBD output:",
  OUT_IBD,
  "\n"
)

cat(
  "Active-IBD output:",
  OUT_ACTIVE,
  "\n"
)

# ============================================================
# 4. Save cohort annotation + manifests
# ============================================================

section("4. WRITE ANNOTATION AND MANIFESTS")

ibd_annotation_file <- file.path(
  OUT_IBD,
  "annotation/MSCCR_CORE12_IBD_annotation.tsv"
)

active_annotation_file <- file.path(
  OUT_ACTIVE,
  "annotation/MSCCR_CORE12_ACTIVEIBD_annotation.tsv"
)

fwrite(
  ann_ibd,
  ibd_annotation_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  ann_active,
  active_annotation_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

writeLines(
  ann_ibd$ID,
  file.path(
    OUT_IBD,
    "sample_manifest.txt"
  )
)

writeLines(
  ann_active$ID,
  file.path(
    OUT_ACTIVE,
    "sample_manifest.txt"
  )
)

cat("[PASS] Annotation/manifests written.\n")

# ============================================================
# 5. Audit source expression matrices
# ============================================================

section("5. SOURCE MATRIX AUDIT")

expected_files <- file.path(
  SOURCE_MATRIX_DIR,
  paste0(
    CORE12,
    ".txt"
  )
)

missing_files <- expected_files[
  !file.exists(expected_files)
]

if (length(missing_files) > 0) {
  fail(
    "Missing CORE12 source matrices:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

qc_rows <- vector(
  "list",
  length(CORE12)
)

names(qc_rows) <- CORE12

# ============================================================
# 6. Subset each matrix
# ============================================================

section("6. BUILD C/D EXPRESSION MATRICES")

for (i in seq_along(CORE12)) {

  ct <- CORE12[[i]]

  cat(
    "\n[",
    i,
    "/",
    length(CORE12),
    "] ",
    ct,
    "\n",
    sep = ""
  )

  src <- file.path(
    SOURCE_MATRIX_DIR,
    paste0(
      ct,
      ".txt"
    )
  )

  cat("  reading source matrix...\n")

  x <- fread(
    src,
    data.table = TRUE,
    check.names = FALSE,
    showProgress = TRUE
  )

  if (ncol(x) < 2L) {
    fail(
      ct,
      ": matrix has fewer than 2 columns."
    )
  }

  gene_col <- colnames(x)[1L]
  source_samples <- colnames(x)[-1L]

  if (length(source_samples) != EXPECTED_ALL) {
    fail(
      ct,
      ": expected ",
      EXPECTED_ALL,
      " source samples, found ",
      length(source_samples)
    )
  }

  if (anyDuplicated(source_samples)) {
    fail(
      ct,
      ": duplicated sample columns detected."
    )
  }

  if (!setequal(
    source_samples,
    ann$ID
  )) {
    missing_in_matrix <- setdiff(
      ann$ID,
      source_samples
    )

    extra_in_matrix <- setdiff(
      source_samples,
      ann$ID
    )

    fail(
      ct,
      ": source sample set does not match annotation.\n",
      "Missing in matrix: ",
      paste(
        head(missing_in_matrix, 10),
        collapse = ", "
      ),
      "\nExtra in matrix: ",
      paste(
        head(extra_in_matrix, 10),
        collapse = ", "
      )
    )
  }

  # ----------------------------------------------------------
  # Reorder using annotation-derived canonical sample order.
  # ----------------------------------------------------------

  ibd_cols <- c(
    gene_col,
    ann_ibd$ID
  )

  active_cols <- c(
    gene_col,
    ann_active$ID
  )

  if (!all(
    ibd_cols %in% colnames(x)
  )) {
    fail(
      ct,
      ": IBD subset columns are incomplete."
    )
  }

  if (!all(
    active_cols %in% colnames(x)
  )) {
    fail(
      ct,
      ": Active-IBD subset columns are incomplete."
    )
  }

  x_ibd <- x[
    ,
    ..ibd_cols
  ]

  x_active <- x[
    ,
    ..active_cols
  ]

  # ----------------------------------------------------------
  # Hard dimensional/order checks before writing
  # ----------------------------------------------------------

  if (ncol(x_ibd) != EXPECTED_IBD + 1L) {
    fail(
      ct,
      ": incorrect IBD-only column count."
    )
  }

  if (ncol(x_active) != EXPECTED_ACTIVE + 1L) {
    fail(
      ct,
      ": incorrect Active-IBD column count."
    )
  }

  if (!identical(
    colnames(x_ibd)[-1L],
    ann_ibd$ID
  )) {
    fail(
      ct,
      ": IBD-only sample order mismatch."
    )
  }

  if (!identical(
    colnames(x_active)[-1L],
    ann_active$ID
  )) {
    fail(
      ct,
      ": Active-IBD sample order mismatch."
    )
  }

  out_ibd <- file.path(
    OUT_IBD,
    "expression_matrices",
    paste0(
      ct,
      ".txt"
    )
  )

  out_active <- file.path(
    OUT_ACTIVE,
    "expression_matrices",
    paste0(
      ct,
      ".txt"
    )
  )

  cat(
    "  writing IBD-only: ",
    EXPECTED_IBD,
    " samples\n",
    sep = ""
  )

  fwrite(
    x_ibd,
    out_ibd,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  cat(
    "  writing Active-IBD: ",
    EXPECTED_ACTIVE,
    " samples\n",
    sep = ""
  )

  fwrite(
    x_active,
    out_active,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  # ----------------------------------------------------------
  # Re-read header only / lightweight output QC
  # ----------------------------------------------------------

  h_ibd <- fread(
    out_ibd,
    nrows = 0,
    check.names = FALSE
  )

  h_active <- fread(
    out_active,
    nrows = 0,
    check.names = FALSE
  )

  if (!identical(
    colnames(h_ibd)[-1L],
    ann_ibd$ID
  )) {
    fail(
      ct,
      ": written IBD matrix header mismatch."
    )
  }

  if (!identical(
    colnames(h_active)[-1L],
    ann_active$ID
  )) {
    fail(
      ct,
      ": written Active-IBD matrix header mismatch."
    )
  }

  qc_rows[[i]] <- data.frame(
    cell_type = ct,
    gene_column = gene_col,
    n_genes = nrow(x),
    n_source_samples = length(source_samples),
    n_ibd_samples = ncol(x_ibd) - 1L,
    n_active_samples = ncol(x_active) - 1L,
    ibd_header_order_ok = TRUE,
    active_header_order_ok = TRUE,
    stringsAsFactors = FALSE
  )

  cat("  [PASS]\n")

  rm(
    x,
    x_ibd,
    x_active,
    h_ibd,
    h_active
  )

  gc()
}

# ============================================================
# 7. Final QC
# ============================================================

section("7. FINAL QC")

qc <- do.call(
  rbind,
  qc_rows
)

fwrite(
  qc,
  file.path(
    OUT_IBD,
    "input_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

fwrite(
  qc,
  file.path(
    OUT_ACTIVE,
    "input_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

print(
  qc,
  row.names = FALSE
)

all_pass <- (
  nrow(qc) == 12L &&
    all(qc$n_source_samples == EXPECTED_ALL) &&
    all(qc$n_ibd_samples == EXPECTED_IBD) &&
    all(qc$n_active_samples == EXPECTED_ACTIVE) &&
    all(qc$ibd_header_order_ok) &&
    all(qc$active_header_order_ok)
)

cat(
  "\nIBD-only annotation:",
  nrow(ann_ibd),
  "samples\n"
)

cat(
  "Active-IBD annotation:",
  nrow(ann_active),
  "samples\n"
)

cat(
  "IBD-only Region3:\n"
)

print(
  table(
    ann_ibd$Region3
  )
)

cat(
  "\nActive-IBD Region3:\n"
)

print(
  table(
    ann_active$Region3
  )
)

if (!all_pass) {
  cat(
    "\nFINAL STATUS: REVIEW_REQUIRED\n"
  )

  quit(
    save = "no",
    status = 1
  )
}

cat(
  "\nFINAL STATUS: PASS_SENSITIVITY_INPUTS_READY\n"
)
