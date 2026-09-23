suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 3
#
# Build canonical expression matrix from FINAL frozen
# gene decision table.
#
# RULES:
#   - subset retained source rows only
#   - preserve deposited expression values EXACTLY
#   - no normalization
#   - no log / exp / shift
#   - no imputation
#   - no averaging / summing
#   - no duplicate collapse
#   - only row identity changes:
#         SourceSymbol -> FinalSymbol
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

expr_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE93624_authoritative_expression_245.rds"
)

decision_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE93624_canonical_gene_decision_table_FINAL.tsv"
)

out_dir <- file.path(
  base,
  "prepared/03_canonical_expression"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE93624 canonical expression matrix build\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load authoritative expression
# ============================================================

obj <- readRDS(
  expr_file
)

expr <- obj$expression

stopifnot(
  is.matrix(expr),
  nrow(expr) == 13769,
  ncol(expr) == 245,
  uniqueN(rownames(expr)) == 13769,
  !anyNA(expr),
  all(is.finite(expr))
)

cat(
  "[PASS] Authoritative expression:",
  nrow(expr),
  "x",
  ncol(expr),
  "\n"
)


# ============================================================
# 2. Load FINAL decision table
# ============================================================

decision <- fread(
  decision_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  nrow(decision) == 13769,
  !anyDuplicated(
    decision$SourceSymbol
  )
)

cat(
  "[PASS] FINAL decision table:",
  nrow(decision),
  "features\n"
)


# ============================================================
# 3. Robust logical parsing
# ============================================================

parse_bool <- function(x) {

  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA,
    length(x)
  )

  out[
    x %chin% c(
      "TRUE",
      "T",
      "1"
    )
  ] <- TRUE

  out[
    x %chin% c(
      "FALSE",
      "F",
      "0"
    )
  ] <- FALSE

  out
}


decision[
  ,
  IncludeCanonical_bool :=
    parse_bool(
      IncludeCanonical
    )
]

stopifnot(
  !anyNA(
    decision$IncludeCanonical_bool
  )
)


# ============================================================
# 4. Authoritative feature-order concordance
# ============================================================

decision[
  ,
  FeatureOrder_int :=
    as.integer(
      FeatureOrder
    )
]

stopifnot(
  !anyNA(
    decision$FeatureOrder_int
  )
)

setorder(
  decision,
  FeatureOrder_int
)

stopifnot(
  identical(
    decision$FeatureOrder_int,
    1:13769
  ),
  identical(
    decision$SourceSymbol,
    rownames(expr)
  )
)

cat(
  "[PASS] FINAL decision table is exactly aligned to",
  " authoritative feature order\n"
)


# ============================================================
# 5. Select included source features
# ============================================================

included <- decision[
  IncludeCanonical_bool == TRUE
]

excluded <- decision[
  IncludeCanonical_bool == FALSE
]

stopifnot(
  nrow(included) == 13323,
  nrow(excluded) == 446
)

stopifnot(
  all(
    !is.na(
      included$FinalGeneID
    )
  ),
  all(
    !is.na(
      included$FinalSymbol
    )
  ),
  uniqueN(
    included$FinalGeneID
  ) == 13323,
  uniqueN(
    included$FinalSymbol
  ) == 13323
)

cat(
  "[PASS] Included canonical features:",
  nrow(included),
  "\n"
)

cat(
  "[PASS] Excluded raw features:",
  nrow(excluded),
  "\n"
)


# ============================================================
# 6. Extract source rows
#
# IMPORTANT:
# This is the only expression operation in this script.
# ============================================================

source_symbols <-
  included$SourceSymbol

stopifnot(
  all(
    source_symbols %chin%
      rownames(expr)
  )
)

canonical_source <- expr[
  source_symbols,
  ,
  drop = FALSE
]

stopifnot(
  nrow(canonical_source) == 13323,
  ncol(canonical_source) == 245,
  identical(
    rownames(canonical_source),
    source_symbols
  ),
  identical(
    colnames(canonical_source),
    colnames(expr)
  )
)


# ============================================================
# 7. Exact source-value concordance BEFORE renaming
# ============================================================

source_reference <- expr[
  match(
    source_symbols,
    rownames(expr)
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    canonical_source,
    source_reference
  )
)

cat(
  "[PASS] Extracted expression values are bitwise-identical",
  " to authoritative source rows\n"
)


# ============================================================
# 8. Rename rows to FINAL current symbols
#
# Expression values remain untouched.
# ============================================================

canonical_expr <-
  canonical_source

rownames(
  canonical_expr
) <- included$FinalSymbol


stopifnot(
  nrow(canonical_expr) == 13323,
  ncol(canonical_expr) == 245,
  uniqueN(
    rownames(canonical_expr)
  ) == 13323,
  identical(
    rownames(canonical_expr),
    included$FinalSymbol
  ),
  identical(
    colnames(canonical_expr),
    colnames(expr)
  )
)


# ============================================================
# 9. Exact value concordance AFTER renaming
#
# Ignore dimnames; values must remain exactly identical.
# ============================================================

stopifnot(
  identical(
    unname(canonical_expr),
    unname(canonical_source)
  )
)

cat(
  "[PASS] Row renaming changed no expression value\n"
)


# ============================================================
# 10. Feature manifest
# ============================================================

feature_manifest <- included[
  ,
  .(
    CanonicalRow =
      seq_len(.N),

    OriginalFeatureOrder =
      FeatureOrder_int,

    SourceSymbol,

    FinalGeneID,

    FinalSymbol,

    EvidenceClass,

    EvidenceDetail,

    OriginalCandidateGeneID,

    OriginalCandidateSymbol,

    FinalResolvedGeneID,

    FinalResolvedSymbol,

    FinalGeneType,

    FinalDescription,

    DynamicGeneIDRedirect,

    HistoricalGeneIDRedirect,

    Decision,
    DecisionReason
  )
]


stopifnot(
  identical(
    feature_manifest$FinalSymbol,
    rownames(canonical_expr)
  ),
  uniqueN(
    feature_manifest$FinalGeneID
  ) == nrow(feature_manifest),
  uniqueN(
    feature_manifest$FinalSymbol
  ) == nrow(feature_manifest)
)


# ============================================================
# 11. Sample manifest
# ============================================================

if (
  !is.null(
    obj$sample_manifest
  )
) {

  sample_manifest <-
    obj$sample_manifest

} else {

  sample_manifest <- data.table(
    SampleTitle =
      colnames(expr)
  )
}


stopifnot(
  nrow(sample_manifest) == 245
)


# ============================================================
# 12. Expression distribution audit
# ============================================================

matrix_audit <- data.table(
  Metric = c(
    "Rows",
    "Columns",
    "Min",
    "Max",
    "Negative_N",
    "NegativeFraction",
    "Zero_N",
    "ZeroFraction",
    "Positive_N",
    "PositiveFraction"
  ),

  Value = c(
    nrow(canonical_expr),
    ncol(canonical_expr),

    min(canonical_expr),
    max(canonical_expr),

    sum(
      canonical_expr < 0
    ),

    mean(
      canonical_expr < 0
    ),

    sum(
      canonical_expr == 0
    ),

    mean(
      canonical_expr == 0
    ),

    sum(
      canonical_expr > 0
    ),

    mean(
      canonical_expr > 0
    )
  )
)


cat("\n============================================================\n")
cat("CANONICAL MATRIX AUDIT\n")
cat("============================================================\n")

print(
  matrix_audit
)


# ============================================================
# 13. Per-sample audit
# ============================================================

sample_audit <- data.table(
  SampleTitle =
    colnames(canonical_expr),

  Negative_N =
    colSums(
      canonical_expr < 0
    ),

  Zero_N =
    colSums(
      canonical_expr == 0
    ),

  Positive_N =
    colSums(
      canonical_expr > 0
    ),

  NegativeFraction =
    colMeans(
      canonical_expr < 0
    ),

  ZeroFraction =
    colMeans(
      canonical_expr == 0
    ),

  Min =
    apply(
      canonical_expr,
      2,
      min
    ),

  Median =
    apply(
      canonical_expr,
      2,
      median
    ),

  Mean =
    colMeans(
      canonical_expr
    ),

  Max =
    apply(
      canonical_expr,
      2,
      max
    )
)


# ============================================================
# 14. Source-to-canonical transformation audit
# ============================================================

transformation_audit <- data.table(
  Check = c(
    "Raw_authoritative_rows",
    "Canonical_rows",
    "Excluded_rows",
    "Samples",
    "Canonical_unique_FinalGeneID",
    "Canonical_unique_FinalSymbol",
    "Expression_values_changed",
    "Expression_normalization_performed",
    "Expression_aggregation_performed",
    "Expression_imputation_performed"
  ),

  Value = c(
    "13769",
    as.character(
      nrow(canonical_expr)
    ),
    as.character(
      nrow(excluded)
    ),
    as.character(
      ncol(canonical_expr)
    ),
    as.character(
      uniqueN(
        feature_manifest$FinalGeneID
      )
    ),
    as.character(
      uniqueN(
        feature_manifest$FinalSymbol
      )
    ),
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE"
  )
)


# ============================================================
# 15. Build canonical object
# ============================================================

canonical_obj <- list(

  expression =
    canonical_expr,

  feature_manifest =
    feature_manifest,

  excluded_feature_manifest =
    excluded,

  sample_manifest =
    sample_manifest,

  provenance = list(

    accession =
      "GSE93624",

    source_expression =
      basename(expr_file),

    source_decision_table =
      basename(decision_file),

    raw_dimensions =
      c(
        genes = 13769L,
        samples = 245L
      ),

    canonical_dimensions =
      c(
        genes = 13323L,
        samples = 245L
      ),

    expression_operation =
      paste0(
        "Row subset only followed by rowname replacement ",
        "SourceSymbol -> FinalSymbol"
      ),

    expression_values_modified =
      FALSE,

    normalization_performed =
      FALSE,

    aggregation_performed =
      FALSE,

    imputation_performed =
      FALSE,

    scale =
      paste0(
        "GEO deposited continuous normalized processed ",
        "expression; exact post-TMM transformation unresolved"
      )
  )
)


# ============================================================
# 16. Save RDS
# ============================================================

rds_file <- file.path(
  out_dir,
  "GSE93624_canonical_expression_13323x245.rds"
)

saveRDS(
  canonical_obj,
  rds_file,
  compress = FALSE
)


# ============================================================
# 17. Save expression TSV
# ============================================================

export_expr <- data.table(
  GeneSymbol =
    rownames(canonical_expr),

  GeneID =
    feature_manifest$FinalGeneID
)

export_expr <- cbind(
  export_expr,
  as.data.table(
    canonical_expr
  )
)


tsv_file <- file.path(
  out_dir,
  "GSE93624_canonical_expression_13323x245.tsv.gz"
)

fwrite(
  export_expr,
  tsv_file,
  sep = "\t",
  na = ""
)


# ============================================================
# 18. Save manifests / audits
# ============================================================

fwrite(
  feature_manifest,
  file.path(
    out_dir,
    "GSE93624_canonical_feature_manifest.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  excluded,
  file.path(
    out_dir,
    "GSE93624_excluded_feature_manifest.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  sample_manifest,
  file.path(
    out_dir,
    "GSE93624_sample_manifest_expression_stage.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  matrix_audit,
  file.path(
    out_dir,
    "GSE93624_canonical_matrix_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  sample_audit,
  file.path(
    out_dir,
    "GSE93624_canonical_per_sample_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  transformation_audit,
  file.path(
    out_dir,
    "GSE93624_expression_transformation_audit.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 19. Reload RDS and verify disk roundtrip
# ============================================================

reload <- readRDS(
  rds_file
)

stopifnot(
  identical(
    reload$expression,
    canonical_expr
  ),
  identical(
    reload$feature_manifest$FinalSymbol,
    rownames(canonical_expr)
  )
)

cat(
  "\n[PASS] Saved RDS roundtrip is exact\n"
)


# ============================================================
# 20. MD5 provenance
# ============================================================

md5 <- tools::md5sum(
  c(
    rds_file,
    tsv_file,
    decision_file,
    expr_file
  )
)

md5_table <- data.table(
  File =
    names(md5),

  MD5 =
    unname(md5)
)

fwrite(
  md5_table,
  file.path(
    out_dir,
    "GSE93624_canonical_expression_md5.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 21. FINAL hard guards
# ============================================================

cat("\n============================================================\n")
cat("FINAL INTEGRITY CHECKS\n")
cat("============================================================\n")

stopifnot(
  nrow(canonical_expr) == 13323,
  ncol(canonical_expr) == 245,

  uniqueN(
    rownames(canonical_expr)
  ) == 13323,

  uniqueN(
    feature_manifest$FinalGeneID
  ) == 13323,

  identical(
    unname(canonical_expr),
    unname(
      expr[
        feature_manifest$SourceSymbol,
        ,
        drop = FALSE
      ]
    )
  ),

  identical(
    colnames(canonical_expr),
    colnames(expr)
  )
)

cat("[PASS] Canonical dimensions = 13,323 x 245\n")
cat("[PASS] Canonical FinalGeneID values are unique\n")
cat("[PASS] Canonical FinalSymbol values are unique\n")
cat("[PASS] Sample order unchanged\n")
cat("[PASS] Expression values unchanged exactly\n")
cat("[PASS] No normalization performed\n")
cat("[PASS] No aggregation performed\n")
cat("[PASS] No imputation performed\n")

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_CANONICAL_EXPRESSION_MATRIX_BUILT\n")
cat("============================================================\n")
