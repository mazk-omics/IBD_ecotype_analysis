suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE57945
# Stage 1: build authoritative 322-sample expression matrix
#
# Input:
#   raw/per_sample/GSM*_CCFA_Risk_*.txt.gz
#
# Important:
#   - NO gene remapping
#   - NO symbol correction
#   - NO filtering
#   - NO transformation
#   - NO normalization
#
# The purpose is to reconstruct and freeze the current GEO
# sample-level processed expression dataset.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

input_dir <- file.path(
  base,
  "raw",
  "per_sample"
)

out_dir <- file.path(
  base,
  "prepared",
  "01_authoritative_matrix"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE57945 authoritative 322-sample matrix construction\n")
cat("============================================================\n\n")


# ============================================================
# 1. Locate sample files
# ============================================================

files <- list.files(
  input_dir,
  pattern = "^GSM[0-9]+_CCFA_Risk_[0-9]+\\.txt\\.gz$",
  full.names = TRUE
)

if (length(files) != 322) {
  stop(
    "Expected 322 individual sample files; found ",
    length(files)
  )
}


# ============================================================
# 2. Parse filenames
# ============================================================

bn <- basename(files)

gsm <- sub(
  "^(GSM[0-9]+)_.*$",
  "\\1",
  bn
)

risk_id <- sub(
  "^GSM[0-9]+_(CCFA_Risk_[0-9]+)\\.txt\\.gz$",
  "\\1",
  bn
)

risk_number <- as.integer(
  sub("^CCFA_Risk_", "", risk_id)
)

manifest <- data.table(
  file = files,
  filename = bn,
  GSM = gsm,
  CCFA_Risk_ID = risk_id,
  RiskNumber = risk_number
)

setorder(
  manifest,
  RiskNumber
)


# ============================================================
# 3. Filename-level checks
# ============================================================

stopifnot(
  nrow(manifest) == 322,
  uniqueN(manifest$GSM) == 322,
  uniqueN(manifest$CCFA_Risk_ID) == 322,
  identical(
    manifest$RiskNumber,
    1:322
  )
)

expected_ids <- sprintf(
  "CCFA_Risk_%03d",
  1:322
)

if (!identical(
  manifest$CCFA_Risk_ID,
  expected_ids
)) {
  stop(
    "CCFA_Risk IDs are not exactly 001-322 in order."
  )
}

cat("[PASS] 322 sample files identified\n")
cat("[PASS] GSM unique: 322\n")
cat("[PASS] CCFA_Risk_ID unique: 322\n")
cat("[PASS] RiskNumber exactly 1-322\n\n")


# ============================================================
# 4. Read reference sample
# ============================================================

cat("[1] Reading reference sample:\n")
cat("    ", manifest$filename[1], "\n\n")

ref <- fread(
  manifest$file[1],
  showProgress = FALSE
)

if (ncol(ref) != 3) {
  stop(
    "Reference file does not contain exactly 3 columns."
  )
}

n_gene <- nrow(ref)

if (n_gene != 36372) {
  stop(
    "Expected 36,372 gene rows; found ",
    n_gene
  )
}

gene_id <- as.character(
  ref[[1]]
)

gene_symbol_raw <- as.character(
  ref[[2]]
)

if (anyNA(gene_id)) {
  stop("Gene ID contains NA.")
}

if (anyDuplicated(gene_id)) {
  stop("Gene ID is not unique.")
}

cat("[PASS] Reference rows:", n_gene, "\n")
cat("[PASS] Gene IDs unique:", uniqueN(gene_id), "\n\n")


# ============================================================
# 5. Allocate expression matrix
# ============================================================

expr <- matrix(
  NA_real_,
  nrow = n_gene,
  ncol = 322,
  dimnames = list(
    gene_id,
    manifest$CCFA_Risk_ID
  )
)


# ============================================================
# 6. Read all 322 files with strict assertions
# ============================================================

qc_list <- vector(
  "list",
  322
)

cat("[2] Reading and validating all individual samples...\n\n")

for (i in seq_len(nrow(manifest))) {

  f <- manifest$file[i]

  dt <- fread(
    f,
    showProgress = FALSE
  )

  if (nrow(dt) != n_gene) {
    stop(
      "Row-count mismatch in ",
      basename(f),
      ": ",
      nrow(dt)
    )
  }

  if (ncol(dt) != 3) {
    stop(
      "Column-count mismatch in ",
      basename(f)
    )
  }

  this_gene_id <- as.character(
    dt[[1]]
  )

  this_symbol <- as.character(
    dt[[2]]
  )

  this_sample_col <- names(dt)[3]

  # Gene ID order must match exactly
  if (!identical(
    this_gene_id,
    gene_id
  )) {
    stop(
      "Gene ID/order mismatch in ",
      basename(f)
    )
  }

  # Raw symbol annotation must match the other individual files
  if (!identical(
    this_symbol,
    gene_symbol_raw
  )) {
    stop(
      "Raw Gene Symbol annotation mismatch in ",
      basename(f)
    )
  }

  # Third-column sample name must match filename
  expected_sample <- manifest$CCFA_Risk_ID[i]

  if (!identical(
    this_sample_col,
    expected_sample
  )) {
    stop(
      "Sample-column mismatch in ",
      basename(f),
      ". Expected ",
      expected_sample,
      "; observed ",
      this_sample_col
    )
  }

  x <- suppressWarnings(
    as.numeric(dt[[3]])
  )

  if (anyNA(x)) {
    stop(
      "NA/non-numeric expression value in ",
      basename(f)
    )
  }

  if (any(!is.finite(x))) {
    stop(
      "Non-finite expression value in ",
      basename(f)
    )
  }

  if (any(x < 0)) {
    stop(
      "Negative expression value in ",
      basename(f)
    )
  }

  expr[, i] <- x

  positive <- x[x > 0]

  qc_list[[i]] <- data.table(
    GSM = manifest$GSM[i],
    CCFA_Risk_ID = expected_sample,
    RiskNumber = manifest$RiskNumber[i],

    TotalExpression = sum(x),
    DetectedGenes = sum(x > 0),
    ZeroFraction = mean(x == 0),

    MedianExpression = median(x),
    MedianNonzeroExpression =
      if (length(positive) > 0)
        median(positive)
      else
        NA_real_,

    MaxExpression = max(x)
  )

  if (
    i %% 25 == 0 ||
    i == 1 ||
    i == 322
  ) {
    cat(
      sprintf(
        "  [%3d/322] %s\n",
        i,
        expected_sample
      )
    )
  }
}

sample_qc <- rbindlist(
  qc_list
)

cat("\n[PASS] All 322 individual files validated\n\n")


# ============================================================
# 7. Global expression checks
# ============================================================

if (anyNA(expr)) {
  stop("Final expression matrix contains NA.")
}

if (any(!is.finite(expr))) {
  stop("Final expression matrix contains non-finite values.")
}

if (any(expr < 0)) {
  stop("Final expression matrix contains negative values.")
}

if (!identical(
  rownames(expr),
  gene_id
)) {
  stop("Expression matrix row names changed unexpectedly.")
}

if (!identical(
  colnames(expr),
  expected_ids
)) {
  stop("Expression matrix sample order changed unexpectedly.")
}


global_zero_fraction <- mean(
  expr == 0
)

global_min <- min(expr)
global_max <- max(expr)


# ============================================================
# 8. Raw gene annotation table
# ============================================================

gene_annotation <- data.table(
  GeneID_raw = gene_id,
  GeneSymbol_raw = gene_symbol_raw
)

gene_annotation[
  ,
  RawSymbolMissing :=
    is.na(GeneSymbol_raw) |
    GeneSymbol_raw == ""
]

gene_annotation[
  ,
  RawSymbolDuplicated :=
    duplicated(GeneSymbol_raw) |
    duplicated(
      GeneSymbol_raw,
      fromLast = TRUE
    )
]

n_unique_symbols <- uniqueN(
  gene_symbol_raw
)

n_duplicate_symbol_rows <- sum(
  gene_annotation$RawSymbolDuplicated
)

n_duplicate_symbol_names <- gene_annotation[
  RawSymbolDuplicated == TRUE,
  uniqueN(GeneSymbol_raw)
]


# ============================================================
# 9. Object provenance
# ============================================================

provenance <- list(
  accession = "GSE57945",

  source = paste0(
    "Current GEO individual processed files ",
    "from GSE57945_RAW.tar"
  ),

  n_samples = 322L,
  n_features = n_gene,

  sample_range =
    "CCFA_Risk_001-CCFA_Risk_322",

  expression_status =
    "GEO-provided processed expression values; unmodified",

  gene_id_status =
    "Raw Gene ID preserved exactly as GEO",

  gene_symbol_status = paste0(
    "Raw Gene Symbol preserved exactly as GEO individual files; ",
    "not considered canonical because historical symbol corruption ",
    "has been detected and will be corrected in a later stage."
  ),

  transformation =
    "None",

  filtering =
    "None",

  normalization =
    "None",

  build_time =
    as.character(Sys.time()),

  R_version =
    R.version.string
)


# ============================================================
# 10. Save R object
# ============================================================

obj <- list(
  expression = expr,
  gene_annotation = gene_annotation,
  sample_manifest = manifest[
    ,
    .(
      GSM,
      CCFA_Risk_ID,
      RiskNumber,
      filename
    )
  ],
  sample_qc = sample_qc,
  provenance = provenance
)

rds_file <- file.path(
  out_dir,
  "GSE57945_authoritative_expression_322.rds"
)

saveRDS(
  obj,
  rds_file,
  compress = "gzip"
)


# ============================================================
# 11. Save expression TSV
# ============================================================

cat("[3] Writing compressed expression matrix...\n")

expr_dt <- data.table(
  GeneID_raw = gene_id,
  GeneSymbol_raw = gene_symbol_raw
)

expr_dt <- cbind(
  expr_dt,
  as.data.table(expr)
)

tsv_file <- file.path(
  out_dir,
  "GSE57945_authoritative_expression_322.tsv.gz"
)

fwrite(
  expr_dt,
  tsv_file,
  sep = "\t",
  compress = "gzip"
)

rm(expr_dt)
gc()


# ============================================================
# 12. Save annotation / manifest / QC
# ============================================================

fwrite(
  gene_annotation,
  file.path(
    out_dir,
    "GSE57945_gene_annotation_raw.tsv"
  ),
  sep = "\t"
)

fwrite(
  manifest[
    ,
    .(
      GSM,
      CCFA_Risk_ID,
      RiskNumber,
      filename
    )
  ],
  file.path(
    out_dir,
    "GSE57945_sample_manifest.tsv"
  ),
  sep = "\t"
)

fwrite(
  sample_qc,
  file.path(
    out_dir,
    "GSE57945_sample_expression_QC.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 13. Summary
# ============================================================

summary_dt <- data.table(
  Metric = c(
    "Samples",
    "Genes",
    "Unique_Gene_IDs",
    "Unique_raw_symbols",
    "Duplicated_raw_symbol_names",
    "Rows_with_duplicated_raw_symbol",
    "Missing_raw_symbol_rows",
    "NA_expression_values",
    "Negative_expression_values",
    "Global_zero_fraction",
    "Global_min_expression",
    "Global_max_expression"
  ),

  Value = as.character(c(
    ncol(expr),
    nrow(expr),
    uniqueN(gene_id),
    n_unique_symbols,
    n_duplicate_symbol_names,
    n_duplicate_symbol_rows,
    sum(gene_annotation$RawSymbolMissing),
    sum(is.na(expr)),
    sum(expr < 0),
    global_zero_fraction,
    global_min,
    global_max
  ))
)

fwrite(
  summary_dt,
  file.path(
    out_dir,
    "GSE57945_authoritative_matrix_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 14. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL MATRIX SUMMARY\n")
cat("============================================================\n")

cat(
  "Dimensions                :",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)

cat(
  "Unique Gene IDs           :",
  uniqueN(gene_id),
  "\n"
)

cat(
  "Unique raw Gene Symbols   :",
  n_unique_symbols,
  "\n"
)

cat(
  "Duplicated symbol names   :",
  n_duplicate_symbol_names,
  "\n"
)

cat(
  "Rows with duplicate symbol:",
  n_duplicate_symbol_rows,
  "\n"
)

cat(
  "Missing raw symbols       :",
  sum(gene_annotation$RawSymbolMissing),
  "\n"
)

cat(
  "NA expression             :",
  sum(is.na(expr)),
  "\n"
)

cat(
  "Negative expression       :",
  sum(expr < 0),
  "\n"
)

cat(
  "Global zero fraction      :",
  sprintf("%.6f", global_zero_fraction),
  "\n"
)

cat(
  "Expression range          :",
  global_min,
  "to",
  global_max,
  "\n"
)

cat(
  "Detected genes/sample     :",
  min(sample_qc$DetectedGenes),
  "to",
  max(sample_qc$DetectedGenes),
  "\n"
)

cat(
  "Median detected genes     :",
  median(sample_qc$DetectedGenes),
  "\n"
)

cat("\nOutputs:\n")
cat("  ", rds_file, "\n")
cat("  ", tsv_file, "\n")

cat("\n")
cat("FINAL STATUS:\n")
cat("PASS_AUTHORITATIVE_322_MATRIX_BUILT\n")
cat("============================================================\n")
