suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 1: full processed-file audit +
# authoritative expression matrix reconstruction
#
# No normalization.
# No gene annotation modification.
# No sample exclusion.
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
cat("GSE93624 authoritative 245-sample matrix build\n")
cat("============================================================\n\n")


# ============================================================
# 1. Locate files
# ============================================================

files <- list.files(
  sample_dir,
  pattern = "^GSM[0-9]+_Sample[0-9]+\\.txt\\.gz$",
  full.names = TRUE
)

if (length(files) != 245) {
  stop(
    "Expected 245 processed sample files; found ",
    length(files)
  )
}


# ------------------------------------------------------------
# Parse filename metadata
# ------------------------------------------------------------

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

cat("[PASS] Found Sample1-Sample245 with 245 unique GSM IDs\n")


# ============================================================
# 2. Read first file as schema reference
# ============================================================

read_sample <- function(f) {

  x <- fread(
    f,
    header = TRUE,
    sep = "\t",
    na.strings = c(
      "NA",
      "NaN"
    ),
    showProgress = FALSE
  )

  if (ncol(x) != 2) {
    stop(
      "Unexpected number of columns in ",
      basename(f),
      ": ",
      ncol(x)
    )
  }

  x
}


ref <- read_sample(
  manifest$File[1]
)

if (nrow(ref) != 13769) {
  stop(
    "Reference file does not contain 13,769 genes."
  )
}

if (names(ref)[1] != "gene") {
  stop(
    "First column is not named 'gene'."
  )
}

ref_genes <- as.character(
  ref[[1]]
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

cat(
  "Duplicated raw gene-symbol rows:",
  sum(duplicated(ref_genes) |
      duplicated(ref_genes, fromLast = TRUE)),
  "\n\n"
)


# ============================================================
# 3. Allocate matrix + audit all 245 files
# ============================================================

expr <- matrix(
  NA_real_,
  nrow = length(ref_genes),
  ncol = 245
)

rownames(expr) <- ref_genes
colnames(expr) <- manifest$SampleTitle


audit_list <- vector(
  "list",
  length(files)
)


for (i in seq_len(nrow(manifest))) {

  f <- manifest$File[i]

  x <- read_sample(f)

  expected_sample <- manifest$SampleTitle[i]

  # ----------------------------------------------------------
  # Schema checks
  # ----------------------------------------------------------

  if (nrow(x) != 13769) {
    stop(
      basename(f),
      " has ",
      nrow(x),
      " rows instead of 13,769."
    )
  }

  if (names(x)[1] != "gene") {
    stop(
      basename(f),
      " first column is not 'gene'."
    )
  }

  if (names(x)[2] != expected_sample) {
    stop(
      basename(f),
      " expression column name is ",
      names(x)[2],
      " but expected ",
      expected_sample
    )
  }

  genes <- as.character(
    x[[1]]
  )

  # ----------------------------------------------------------
  # Exact feature-universe/order check
  # ----------------------------------------------------------

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
      " has the same genes but a different order."
    )
  }


  # ----------------------------------------------------------
  # Numeric-expression checks
  # ----------------------------------------------------------

  values <- suppressWarnings(
    as.numeric(
      x[[2]]
    )
  )

  if (length(values) != 13769) {
    stop(
      basename(f),
      " expression length mismatch."
    )
  }

  if (anyNA(values)) {
    stop(
      basename(f),
      " contains NA/non-numeric expression values."
    )
  }

  if (any(!is.finite(values))) {
    stop(
      basename(f),
      " contains non-finite expression values."
    )
  }

  if (any(values < 0)) {
    stop(
      basename(f),
      " contains negative expression values."
    )
  }


  expr[, i] <- values


  # ----------------------------------------------------------
  # Integer-vs-fractional diagnostic
  # ----------------------------------------------------------

  fractional_n <- sum(
    abs(
      values -
      round(values)
    ) > 1e-10
  )

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

    Negative_N =
      sum(values < 0),

    Zero_N =
      sum(values == 0),

    ZeroFraction =
      mean(values == 0),

    Min =
      min(values),

    Median =
      median(values),

    Mean =
      mean(values),

    Max =
      max(values),

    FractionalValue_N =
      fractional_n,

    FractionalValueFraction =
      fractional_n /
      length(values)
  )
}


sample_audit <- rbindlist(
  audit_list
)

cat("[PASS] All 245 files read successfully\n")
cat("[PASS] All 245 files contain identical gene universe\n")
cat("[PASS] All 245 files contain identical gene order\n")
cat("[PASS] No NA, non-finite or negative expression values\n\n")


# ============================================================
# 4. Matrix integrity
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
  all(is.finite(expr)),
  all(expr >= 0)
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
  "Global zero fraction:",
  sprintf(
    "%.6f",
    mean(expr == 0)
  ),
  "\n"
)


# ============================================================
# 5. Feature annotation table
# ============================================================

feature_annotation <- data.table(
  FeatureOrder =
    seq_along(ref_genes),

  GeneSymbol_raw =
    ref_genes
)

feature_annotation[
  ,
  SymbolDuplicated :=
    duplicated(GeneSymbol_raw) |
    duplicated(
      GeneSymbol_raw,
      fromLast = TRUE
    )
]


# ============================================================
# 6. Duplicate-symbol audit
# ============================================================

dup_symbols <- feature_annotation[
  SymbolDuplicated == TRUE
]

cat("\n============================================================\n")
cat("RAW GENE-SYMBOL DUPLICATION\n")
cat("============================================================\n")

cat(
  "Rows in duplicated-symbol groups:",
  nrow(dup_symbols),
  "\n"
)

cat(
  "Unique duplicated symbols:",
  uniqueN(
    dup_symbols$GeneSymbol_raw
  ),
  "\n"
)

if (nrow(dup_symbols) > 0) {

  print(
    dup_symbols[
      ,
      .N,
      by = GeneSymbol_raw
    ][
      order(
        -N,
        GeneSymbol_raw
      )
    ]
  )
}


# ============================================================
# 7. Global numerical diagnostic
# ============================================================

cat("\n============================================================\n")
cat("NUMERICAL FORMAT DIAGNOSTIC\n")
cat("============================================================\n")

all_values <- as.numeric(
  expr
)

fractional_global <- sum(
  abs(
    all_values -
    round(all_values)
  ) > 1e-10
)

cat(
  "Total expression values:",
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

cat(
  "Interpretation: diagnostic only; no transformation applied.\n"
)


# ============================================================
# 8. Sample-level descriptive summary
# ============================================================

cat("\n============================================================\n")
cat("SAMPLE-LEVEL SUMMARY\n")
cat("============================================================\n")

print(
  sample_audit[
    ,
    .(
      Metric = c(
        "DetectedGenes",
        "ZeroFraction",
        "SampleMedian",
        "SampleMean",
        "SampleMax"
      ),

      Min = c(
        min(13769 - Zero_N),
        min(ZeroFraction),
        min(Median),
        min(Mean),
        min(Max)
      ),

      Q1 = c(
        quantile(13769 - Zero_N, 0.25),
        quantile(ZeroFraction, 0.25),
        quantile(Median, 0.25),
        quantile(Mean, 0.25),
        quantile(Max, 0.25)
      ),

      Median = c(
        median(13769 - Zero_N),
        median(ZeroFraction),
        median(Median),
        median(Mean),
        median(Max)
      ),

      Q3 = c(
        quantile(13769 - Zero_N, 0.75),
        quantile(ZeroFraction, 0.75),
        quantile(Median, 0.75),
        quantile(Mean, 0.75),
        quantile(Max, 0.75)
      ),

      Max = c(
        max(13769 - Zero_N),
        max(ZeroFraction),
        max(Median),
        max(Mean),
        max(Max)
      )
    )
  ]
)


# ============================================================
# 9. Save authoritative object
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
  provenance = list(
    accession = "GSE93624",
    expression_type =
      "GEO sample-level processed expression",
    source_gene_identifier =
      "gene symbol",
    geo_processing_description =
      paste(
        "HTSeq gene counts;",
        "genes retained using RPKM/read-count criteria;",
        "original raw counts for retained genes",
        "normalized using edgeR TMM"
      ),
    n_features = 13769L,
    n_samples = 245L
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
# 10. Save expression TSV
# ============================================================

export_dt <- data.table(
  gene = rownames(expr)
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
    "GSE93624_authoritative_expression_245.tsv.gz"
  ),
  sep = "\t"
)


# ============================================================
# 11. Save audit resources
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
  dup_symbols,
  file.path(
    out_dir,
    "GSE93624_duplicated_raw_gene_symbols.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 12. Final hard checks
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
    sample_audit$Negative_N
  ) == 0
)


cat("\nIMPORTANT:\n")
cat(
  "- The 245 GEO sample-level processed files are treated as the source-of-truth expression files.\n",
  "- No normalization was performed by this script.\n",
  "- No log transformation was performed.\n",
  "- No gene symbols were renamed or collapsed.\n",
  "- No samples were excluded.\n",
  "- Gene-symbol harmonization will be handled in a separate stage.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_AUTHORITATIVE_245_MATRIX_BUILT\n")
cat("============================================================\n")
