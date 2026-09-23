suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 1B
#
# Full processed-expression audit +
# authoritative 13,769 x 245 matrix reconstruction
#
# IMPORTANT:
#   - Negative values are allowed because they exist in the
#     deposited GEO processed files.
#   - No expression transformation is performed.
#   - No normalization is performed.
#   - No gene harmonization is performed.
#   - No sample is removed.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

sample_dir <- file.path(
  base,
  "raw/per_sample"
)

out_dir <- file.path(
  base,
  "prepared/01_authoritative_matrix"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE93624 authoritative 245-sample matrix build + scale audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Locate processed files
# ============================================================

files <- list.files(
  sample_dir,
  pattern = "^GSM[0-9]+_Sample[0-9]+\\.txt\\.gz$",
  full.names = TRUE
)

stopifnot(
  length(files) == 245
)

manifest <- data.table(
  File = files,
  Filename = basename(files)
)

manifest[
  ,
  GSM := sub(
    "_.*$",
    "",
    Filename
  )
]

manifest[
  ,
  SampleTitle := sub(
    "^GSM[0-9]+_(Sample[0-9]+)\\.txt\\.gz$",
    "\\1",
    Filename
  )
]

manifest[
  ,
  SampleNumber := as.integer(
    sub(
      "^Sample",
      "",
      SampleTitle
    )
  )
]

stopifnot(
  !anyNA(manifest$SampleNumber),
  uniqueN(manifest$GSM) == 245,
  uniqueN(manifest$SampleTitle) == 245,
  uniqueN(manifest$SampleNumber) == 245
)

setorder(
  manifest,
  SampleNumber
)

stopifnot(
  identical(
    manifest$SampleNumber,
    1:245
  )
)

cat(
  "[PASS] Found Sample1-Sample245 with 245 unique GSM IDs\n"
)


# ============================================================
# 2. Safe reader
# ============================================================

read_sample <- function(f) {

  x <- fread(
    f,
    header = TRUE,
    sep = "\t",
    na.strings = c(
      "NA",
      "NaN",
      "Inf",
      "-Inf"
    ),
    showProgress = FALSE
  )

  if (ncol(x) != 2) {
    stop(
      basename(f),
      " has ",
      ncol(x),
      " columns; expected 2."
    )
  }

  x
}


# ============================================================
# 3. Reference feature schema
# ============================================================

ref <- read_sample(
  manifest$File[1]
)

stopifnot(
  nrow(ref) == 13769,
  names(ref)[1] == "gene"
)

ref_genes <- as.character(
  ref[[1]]
)

stopifnot(
  length(ref_genes) == 13769,
  !anyNA(ref_genes),
  all(ref_genes != "")
)

cat(
  "Reference feature rows:",
  length(ref_genes),
  "\n"
)

cat(
  "Unique raw gene symbols:",
  uniqueN(ref_genes),
  "\n"
)

dup_mask <-
  duplicated(ref_genes) |
  duplicated(
    ref_genes,
    fromLast = TRUE
  )

cat(
  "Duplicated raw gene-symbol rows:",
  sum(dup_mask),
  "\n\n"
)


# ============================================================
# 4. Allocate authoritative matrix
# ============================================================

expr <- matrix(
  NA_real_,
  nrow = 13769,
  ncol = 245,
  dimnames = list(
    ref_genes,
    manifest$SampleTitle
  )
)

audit_list <- vector(
  "list",
  245
)


# ============================================================
# 5. Audit all 245 files
# ============================================================

for (i in seq_len(nrow(manifest))) {

  f <- manifest$File[i]

  x <- read_sample(f)

  expected_sample <-
    manifest$SampleTitle[i]

  # ----------------------------------------------------------
  # Schema
  # ----------------------------------------------------------

  if (nrow(x) != 13769) {
    stop(
      basename(f),
      " has ",
      nrow(x),
      " feature rows; expected 13,769."
    )
  }

  if (names(x)[1] != "gene") {
    stop(
      basename(f),
      ": first column is not 'gene'."
    )
  }

  if (names(x)[2] != expected_sample) {
    stop(
      basename(f),
      ": expression column = ",
      names(x)[2],
      "; expected ",
      expected_sample
    )
  }

  genes <- as.character(
    x[[1]]
  )

  same_gene_set <- setequal(
    genes,
    ref_genes
  )

  same_gene_order <- identical(
    genes,
    ref_genes
  )

  if (!same_gene_set) {
    stop(
      basename(f),
      " has a different gene universe."
    )
  }

  if (!same_gene_order) {
    stop(
      basename(f),
      " has a different gene order."
    )
  }


  # ----------------------------------------------------------
  # Numeric values
  # ----------------------------------------------------------

  values <- suppressWarnings(
    as.numeric(
      x[[2]]
    )
  )

  if (length(values) != 13769) {
    stop(
      basename(f),
      ": numeric vector length mismatch."
    )
  }

  if (anyNA(values)) {
    stop(
      basename(f),
      " contains NA or non-numeric values."
    )
  }

  if (any(!is.finite(values))) {
    stop(
      basename(f),
      " contains non-finite values."
    )
  }

  # IMPORTANT:
  # Negative values are legitimate deposited values
  # and are NOT grounds for stopping.
  expr[, i] <- values


  # ----------------------------------------------------------
  # Numerical diagnostics
  # ----------------------------------------------------------

  fractional_mask <-
    abs(
      values -
      round(values)
    ) > 1e-10

  audit_list[[i]] <- data.table(
    SampleNumber =
      manifest$SampleNumber[i],

    SampleTitle =
      expected_sample,

    GSM =
      manifest$GSM[i],

    Filename =
      manifest$Filename[i],

    NFeatures =
      length(values),

    SameGeneSet =
      same_gene_set,

    SameGeneOrder =
      same_gene_order,

    NA_N =
      sum(is.na(values)),

    NonFinite_N =
      sum(!is.finite(values)),

    Negative_N =
      sum(values < 0),

    Zero_N =
      sum(values == 0),

    Positive_N =
      sum(values > 0),

    NegativeFraction =
      mean(values < 0),

    ZeroFraction =
      mean(values == 0),

    Fractional_N =
      sum(fractional_mask),

    FractionalFraction =
      mean(fractional_mask),

    Min =
      min(values),

    Q01 =
      as.numeric(
        quantile(
          values,
          0.01,
          names = FALSE
        )
      ),

    Q25 =
      as.numeric(
        quantile(
          values,
          0.25,
          names = FALSE
        )
      ),

    Median =
      median(values),

    Mean =
      mean(values),

    Q75 =
      as.numeric(
        quantile(
          values,
          0.75,
          names = FALSE
        )
      ),

    Q99 =
      as.numeric(
        quantile(
          values,
          0.99,
          names = FALSE
        )
      ),

    Max =
      max(values)
  )
}

sample_audit <- rbindlist(
  audit_list
)

cat("[PASS] All 245 files read successfully\n")
cat("[PASS] All files share identical gene universe\n")
cat("[PASS] All files share identical gene order\n")
cat("[PASS] No NA or non-finite expression values\n")
cat("[INFO] Negative deposited expression values are preserved\n\n")


# ============================================================
# 6. Matrix integrity
# ============================================================

stopifnot(
  nrow(expr) == 13769,
  ncol(expr) == 245,
  identical(
    rownames(expr),
    ref_genes
  ),
  identical(
    colnames(expr),
    manifest$SampleTitle
  ),
  !anyNA(expr),
  all(is.finite(expr))
)

cat("============================================================\n")
cat("AUTHORITATIVE MATRIX\n")
cat("============================================================\n")

cat(
  "Dimensions:",
  nrow(expr),
  "x",
  ncol(expr),
  "\n"
)

cat(
  "Unique raw gene symbols:",
  uniqueN(rownames(expr)),
  "\n"
)

cat(
  "Expression range:",
  min(expr),
  "to",
  max(expr),
  "\n"
)

cat(
  "Negative values:",
  sum(expr < 0),
  "/",
  length(expr),
  sprintf(
    "(%.4f%%)",
    100 * mean(expr < 0)
  ),
  "\n"
)

cat(
  "Zero values:",
  sum(expr == 0),
  "/",
  length(expr),
  sprintf(
    "(%.4f%%)",
    100 * mean(expr == 0)
  ),
  "\n"
)

cat(
  "Positive values:",
  sum(expr > 0),
  "/",
  length(expr),
  sprintf(
    "(%.4f%%)",
    100 * mean(expr > 0)
  ),
  "\n"
)


# ============================================================
# 7. Global value distribution
# ============================================================

all_values <- as.numeric(
  expr
)

global_quantiles <- data.table(
  Quantile = c(
    "Min",
    "0.1%",
    "1%",
    "5%",
    "25%",
    "Median",
    "75%",
    "95%",
    "99%",
    "99.9%",
    "Max"
  ),

  Value = c(
    min(all_values),

    as.numeric(
      quantile(
        all_values,
        0.001,
        names = FALSE
      )
    ),

    as.numeric(
      quantile(
        all_values,
        0.01,
        names = FALSE
      )
    ),

    as.numeric(
      quantile(
        all_values,
        0.05,
        names = FALSE
      )
    ),

    as.numeric(
      quantile(
        all_values,
        0.25,
        names = FALSE
      )
    ),

    median(all_values),

    as.numeric(
      quantile(
        all_values,
        0.75,
        names = FALSE
      )
    ),

    as.numeric(
      quantile(
        all_values,
        0.95,
        names = FALSE
      )
    ),

    as.numeric(
      quantile(
        all_values,
        0.99,
        names = FALSE
      )
    ),

    as.numeric(
      quantile(
        all_values,
        0.999,
        names = FALSE
      )
    ),

    max(all_values)
  )
)

cat("\n============================================================\n")
cat("GLOBAL VALUE DISTRIBUTION\n")
cat("============================================================\n")

print(
  global_quantiles
)


# ============================================================
# 8. Decimal / integer diagnostic
# ============================================================

fractional_global <- sum(
  abs(
    all_values -
    round(all_values)
  ) > 1e-10
)

cat("\n============================================================\n")
cat("NUMERICAL FORMAT\n")
cat("============================================================\n")

cat(
  "Total values:",
  length(all_values),
  "\n"
)

cat(
  "Fractional values:",
  fractional_global,
  "\n"
)

cat(
  "Fractional fraction:",
  sprintf(
    "%.6f",
    fractional_global /
      length(all_values)
  ),
  "\n"
)


# ============================================================
# 9. Gene-wise distribution forensic
#
# Paper reports that expression estimates used for the TRS
# analysis were inverse-rank transformed gene-by-gene into
# approximately N(0,1).
#
# If the GEO deposited matrix were that transformed matrix:
#   row means should be approximately 0
#   row SDs should be approximately 1
#
# We test this empirically but do not transform anything.
# ============================================================

cat("\n============================================================\n")
cat("GENE-WISE SCALE FORENSIC\n")
cat("============================================================\n")

gene_mean <- rowMeans(
  expr
)

gene_sd <- apply(
  expr,
  1,
  sd
)

gene_min <- apply(
  expr,
  1,
  min
)

gene_max <- apply(
  expr,
  1,
  max
)

gene_stats <- data.table(
  GeneSymbol_raw =
    rownames(expr),

  MeanAcrossSamples =
    gene_mean,

  SDAcrossSamples =
    gene_sd,

  MinAcrossSamples =
    gene_min,

  MaxAcrossSamples =
    gene_max
)

cat(
  "Median gene mean:",
  median(gene_mean),
  "\n"
)

cat(
  "Median gene SD:",
  median(gene_sd),
  "\n"
)

cat(
  "Genes with |mean| < 0.05:",
  sum(
    abs(gene_mean) < 0.05
  ),
  "/",
  length(gene_mean),
  "\n"
)

cat(
  "Genes with SD between 0.95 and 1.05:",
  sum(
    gene_sd >= 0.95 &
    gene_sd <= 1.05
  ),
  "/",
  length(gene_sd),
  "\n"
)

cat(
  "Genes satisfying BOTH:",
  sum(
    abs(gene_mean) < 0.05 &
    gene_sd >= 0.95 &
    gene_sd <= 1.05
  ),
  "/",
  length(gene_sd),
  "\n"
)

# Very large magnitudes are also informative because
# rank-normalized values across only 245 subjects should not
# ordinarily reach values near 5.
cat(
  "Global |expression| > 4 values:",
  sum(
    abs(expr) > 4
  ),
  "\n"
)

cat(
  "Global |expression| > 5 values:",
  sum(
    abs(expr) > 5
  ),
  "\n"
)


# ============================================================
# 10. Raw gene-symbol audit
# ============================================================

feature_annotation <- data.table(
  FeatureOrder =
    seq_along(ref_genes),

  GeneSymbol_raw =
    ref_genes,

  SymbolDuplicated =
    dup_mask
)

dup_symbols <- feature_annotation[
  SymbolDuplicated == TRUE
]

cat("\n============================================================\n")
cat("RAW GENE-SYMBOL DUPLICATION\n")
cat("============================================================\n")

cat(
  "Unique raw symbols:",
  uniqueN(ref_genes),
  "\n"
)

cat(
  "Rows in duplicated-symbol groups:",
  nrow(dup_symbols),
  "\n"
)


# ============================================================
# 11. Sample-level summary
# ============================================================

cat("\n============================================================\n")
cat("SAMPLE-LEVEL DISTRIBUTION SUMMARY\n")
cat("============================================================\n")

sample_summary <- data.table(
  Metric = c(
    "NegativeFraction",
    "ZeroFraction",
    "SampleMin",
    "SampleMedian",
    "SampleMean",
    "SampleMax"
  ),

  Min = c(
    min(sample_audit$NegativeFraction),
    min(sample_audit$ZeroFraction),
    min(sample_audit$Min),
    min(sample_audit$Median),
    min(sample_audit$Mean),
    min(sample_audit$Max)
  ),

  Q1 = c(
    quantile(sample_audit$NegativeFraction, 0.25),
    quantile(sample_audit$ZeroFraction, 0.25),
    quantile(sample_audit$Min, 0.25),
    quantile(sample_audit$Median, 0.25),
    quantile(sample_audit$Mean, 0.25),
    quantile(sample_audit$Max, 0.25)
  ),

  Median = c(
    median(sample_audit$NegativeFraction),
    median(sample_audit$ZeroFraction),
    median(sample_audit$Min),
    median(sample_audit$Median),
    median(sample_audit$Mean),
    median(sample_audit$Max)
  ),

  Q3 = c(
    quantile(sample_audit$NegativeFraction, 0.75),
    quantile(sample_audit$ZeroFraction, 0.75),
    quantile(sample_audit$Min, 0.75),
    quantile(sample_audit$Median, 0.75),
    quantile(sample_audit$Mean, 0.75),
    quantile(sample_audit$Max, 0.75)
  ),

  Max = c(
    max(sample_audit$NegativeFraction),
    max(sample_audit$ZeroFraction),
    max(sample_audit$Min),
    max(sample_audit$Median),
    max(sample_audit$Mean),
    max(sample_audit$Max)
  )
)

print(
  sample_summary
)


# ============================================================
# 12. Conservative scale classification
# ============================================================

inverse_rank_like_fraction <- mean(
  abs(gene_mean) < 0.05 &
  gene_sd >= 0.95 &
  gene_sd <= 1.05
)

if (
  inverse_rank_like_fraction >= 0.90
) {

  scale_class <-
    "CONSISTENT_WITH_GENEWISE_INVERSE_NORMAL_TRANSFORM"

} else if (
  any(expr < 0) &&
  fractional_global /
    length(all_values) > 0.90
) {

  scale_class <-
    paste0(
      "CONTINUOUS_NORMALIZED_EXPRESSION_WITH_NEGATIVE_VALUES;",
      "_NOT_RAW_COUNT_SCALE;",
      "_EXACT_TRANSFORMATION_UNRESOLVED"
    )

} else {

  scale_class <-
    "PROCESSED_EXPRESSION_SCALE_UNRESOLVED"
}

cat("\n============================================================\n")
cat("CONSERVATIVE SCALE CLASSIFICATION\n")
cat("============================================================\n")

cat(
  "Classification:",
  scale_class,
  "\n"
)

cat(
  "NOTE: this classification does NOT infer a specific",
  " undocumented log transformation.\n"
)


# ============================================================
# 13. Save authoritative object
# ============================================================

authoritative <- list(
  expression = expr,

  feature_annotation =
    feature_annotation,

  sample_manifest =
    manifest[
      ,
      .(
        SampleNumber,
        SampleTitle,
        GSM,
        Filename
      )
    ],

  sample_expression_audit =
    sample_audit,

  gene_expression_stats =
    gene_stats,

  provenance = list(
    accession = "GSE93624",

    source =
      "GEO sample-level supplementary processed expression files",

    source_gene_identifier =
      "gene symbol",

    n_features =
      13769L,

    n_samples =
      245L,

    contains_negative_values =
      any(expr < 0),

    scale_classification =
      scale_class,

    geo_documentation =
      paste(
        "GEO describes supplementary files as normalized",
        "gene-level read counts after edgeR TMM normalization;",
        "deposited numerical values include negative values,",
        "so they are not treated as raw/nonnegative count-scale data."
      )
  )
)

saveRDS(
  authoritative,
  file.path(
    out_dir,
    "GSE93624_authoritative_expression_245.rds"
  ),
  compress = FALSE
)


# ============================================================
# 14. Save matrix TSV
# ============================================================

export_dt <- data.table(
  gene =
    rownames(expr)
)

export_dt <- cbind(
  export_dt,
  as.data.table(
    expr
  )
)

fwrite(
  export_dt,
  file.path(
    out_dir,
    "GSE93624_authoritative_expression_13769x245.tsv.gz"
  ),
  sep = "\t"
)


# ============================================================
# 15. Save audit tables
# ============================================================

fwrite(
  feature_annotation,
  file.path(
    out_dir,
    "GSE93624_feature_annotation_raw.tsv"
  ),
  sep = "\t"
)

fwrite(
  manifest[
    ,
    .(
      SampleNumber,
      SampleTitle,
      GSM,
      Filename
    )
  ],
  file.path(
    out_dir,
    "GSE93624_sample_file_manifest.tsv"
  ),
  sep = "\t"
)

fwrite(
  sample_audit,
  file.path(
    out_dir,
    "GSE93624_per_sample_expression_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_stats,
  file.path(
    out_dir,
    "GSE93624_gene_expression_scale_stats.tsv"
  ),
  sep = "\t"
)

fwrite(
  global_quantiles,
  file.path(
    out_dir,
    "GSE93624_global_expression_quantiles.tsv"
  ),
  sep = "\t"
)

fwrite(
  sample_summary,
  file.path(
    out_dir,
    "GSE93624_sample_expression_distribution_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 16. Final hard checks
# ============================================================

stopifnot(
  all(
    sample_audit$SameGeneSet
  ),

  all(
    sample_audit$SameGeneOrder
  ),

  sum(
    sample_audit$NA_N
  ) == 0,

  sum(
    sample_audit$NonFinite_N
  ) == 0,

  nrow(expr) == 13769,

  ncol(expr) == 245,

  uniqueN(
    rownames(expr)
  ) == 13769
)


cat("\nIMPORTANT:\n")
cat(
  "- Negative values are part of the deposited GEO processed expression data.\n",
  "- No values were clipped, shifted, logged, exponentiated, or otherwise transformed.\n",
  "- No normalization was performed by this script.\n",
  "- No gene symbols were renamed or collapsed.\n",
  "- No samples were excluded.\n",
  "- The exact numerical transformation remains unresolved unless supported by source documentation.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_AUTHORITATIVE_245_MATRIX_BUILT\n")
cat("============================================================\n")
