suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE57945
# Build final canonical gene-level expression matrix
#
# Input:
#   authoritative 36,372 x 322 matrix
#   frozen canonical decision table
#
# Output:
#   34,368 x 322 matrix
#
# No normalization.
# No expression aggregation.
# No imputation.
# No locus collapse.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

expr_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE57945_authoritative_expression_322.rds"
)

decision_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_canonical_gene_decision_table.tsv"
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
cat("GSE57945 canonical expression matrix build\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Load authoritative matrix
# ------------------------------------------------------------

obj <- readRDS(expr_file)

expr <- obj$expression

cat(
  "Authoritative matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)

stopifnot(
  nrow(expr) == 36372,
  ncol(expr) == 322,
  !anyDuplicated(rownames(expr)),
  !anyDuplicated(colnames(expr))
)


# ------------------------------------------------------------
# 2. Load frozen canonical decision table
# ------------------------------------------------------------

d <- fread(
  decision_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  nrow(d) == 36372
)

included <- d[
  CanonicalAction %chin% c(
    "KEEP_DIRECT_CURRENT",
    "KEEP_SAFE_REDIRECT"
  )
]

cat(
  "Included canonical features:",
  nrow(included),
  "\n"
)

stopifnot(
  nrow(included) == 34368,
  uniqueN(included$FinalGeneID) == 34368,
  uniqueN(included$FinalSymbol) == 34368
)


# ------------------------------------------------------------
# 3. Validate source rows
# ------------------------------------------------------------

stopifnot(
  all(included$GeneID_raw %chin% rownames(expr))
)

source_index <- match(
  included$GeneID_raw,
  rownames(expr)
)

stopifnot(
  !anyNA(source_index)
)


# ------------------------------------------------------------
# 4. Subset expression
# ------------------------------------------------------------

canonical_expr <- expr[
  source_index,
  ,
  drop = FALSE
]

# source order exactly follows decision table
stopifnot(
  identical(
    rownames(canonical_expr),
    included$GeneID_raw
  )
)


# ------------------------------------------------------------
# 5. Replace rownames with current canonical symbols
# ------------------------------------------------------------

rownames(canonical_expr) <- included$FinalSymbol

stopifnot(
  nrow(canonical_expr) == 34368,
  ncol(canonical_expr) == 322,
  !anyDuplicated(rownames(canonical_expr)),
  !anyDuplicated(colnames(canonical_expr))
)


# ------------------------------------------------------------
# 6. Expression integrity
# ------------------------------------------------------------

stopifnot(
  all(is.finite(canonical_expr)),
  !anyNA(canonical_expr),
  all(canonical_expr >= 0)
)

cat(
  "Expression range:",
  min(canonical_expr),
  "to",
  max(canonical_expr),
  "\n"
)

cat(
  "Global zero fraction:",
  sprintf(
    "%.6f",
    mean(canonical_expr == 0)
  ),
  "\n"
)


# ------------------------------------------------------------
# 7. Verify expression values were not altered
# ------------------------------------------------------------

original_subset <- expr[
  source_index,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    as.numeric(canonical_expr),
    as.numeric(original_subset)
  )
)

cat(
  "[PASS] Expression values identical to authoritative source rows\n"
)


# ------------------------------------------------------------
# 8. Mapping manifest
# ------------------------------------------------------------

manifest <- included[
  ,
  .(
    SourceGeneID = GeneID_raw,
    SourceRawSymbol = GeneSymbol_raw,
    FinalGeneID,
    FinalSymbol,
    CanonicalAction,
    CanonicalReason
  )
]

stopifnot(
  identical(
    manifest$FinalSymbol,
    rownames(canonical_expr)
  )
)


# ------------------------------------------------------------
# 9. Save RDS
# ------------------------------------------------------------

saveRDS(
  canonical_expr,
  file.path(
    out_dir,
    "GSE57945_canonical_expression_34368x322.rds"
  ),
  compress = FALSE
)


# ------------------------------------------------------------
# 10. Save TSV
# ------------------------------------------------------------

export_dt <- data.table(
  GeneSymbol = rownames(canonical_expr),
  as.data.table(
    canonical_expr,
    keep.rownames = FALSE
  )
)

fwrite(
  export_dt,
  file.path(
    out_dir,
    "GSE57945_canonical_expression_34368x322.tsv.gz"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 11. Save manifest
# ------------------------------------------------------------

fwrite(
  manifest,
  file.path(
    out_dir,
    "GSE57945_canonical_gene_manifest.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 12. Save excluded-feature manifest
# ------------------------------------------------------------

excluded <- d[
  !CanonicalAction %chin% c(
    "KEEP_DIRECT_CURRENT",
    "KEEP_SAFE_REDIRECT"
  )
]

fwrite(
  excluded,
  file.path(
    out_dir,
    "GSE57945_excluded_feature_manifest.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 13. Summary
# ------------------------------------------------------------

summary <- data.table(
  Metric = c(
    "Original_features",
    "Canonical_features",
    "Samples",
    "Unique_final_symbols",
    "Unique_final_geneids",
    "Excluded_features",
    "Zero_fraction"
  ),
  Value = c(
    nrow(expr),
    nrow(canonical_expr),
    ncol(canonical_expr),
    uniqueN(rownames(canonical_expr)),
    uniqueN(manifest$FinalGeneID),
    nrow(excluded),
    mean(canonical_expr == 0)
  )
)

fwrite(
  summary,
  file.path(
    out_dir,
    "GSE57945_canonical_expression_summary.tsv"
  ),
  sep = "\t"
)


cat("\n============================================================\n")
cat("FINAL MATRIX\n")
cat("============================================================\n")

cat(
  "Dimensions:",
  nrow(canonical_expr),
  "x",
  ncol(canonical_expr),
  "\n"
)

cat(
  "Unique symbols:",
  uniqueN(rownames(canonical_expr)),
  "\n"
)

cat(
  "Unique GeneIDs:",
  uniqueN(manifest$FinalGeneID),
  "\n"
)

cat(
  "Excluded original features:",
  nrow(excluded),
  "\n"
)

cat("\nCanonical action composition:\n")

print(
  manifest[
    ,
    .N,
    by = CanonicalAction
  ]
)


cat("\nIMPORTANT:\n")
cat(
  "- Expression values were copied without modification.\n",
  "- No normalization was performed.\n",
  "- No expression rows were summed, averaged, or collapsed.\n",
  "- Row names are current canonical gene symbols.\n",
  "- FinalGeneID mapping is preserved in the manifest.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_CANONICAL_EXPRESSION_MATRIX_BUILT\n")
cat("============================================================\n")
