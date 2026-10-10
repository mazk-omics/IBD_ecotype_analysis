suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE72819 Stage 3
# Build canonical expression matrix
#
# Frozen inputs:
#   Stage 1 authoritative RPKM:
#       26,468 genes x 73 samples
#
#   Stage 2D FINAL gene decision table:
#       24,656 included
#        1,812 excluded
#
# Allowed:
#   - subset FINAL Included rows
#   - rename retained rows to FINAL NCBI symbols
#
# Forbidden:
#   - normalization
#   - log transformation
#   - batch correction
#   - aggregation
#   - imputation
#   - expression recalculation
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

source_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE72819_authoritative_RPKM_26468x73.rds"
)

decision_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE72819_canonical_gene_decision_table_FINAL.tsv"
)

out_dir <- file.path(
  root,
  "prepared/03_canonical_expression"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE72819 Stage 3: build canonical expression matrix\n")
cat("============================================================\n\n")


# ============================================================
# Helper
# ============================================================

check <- function(condition, message) {

  if (!isTRUE(condition)) {
    stop(
      paste0("[FAIL] ", message),
      call. = FALSE
    )
  }

  cat("[PASS] ", message, "\n", sep = "")
}


same_text <- function(x, y) {

  if (length(x) != length(y)) {
    return(FALSE)
  }

  x <- as.character(x)
  y <- as.character(y)

  all(
    !is.na(x) &
    !is.na(y) &
    x == y
  )
}


# ============================================================
# 1. Load authoritative RPKM
# ============================================================

source <- readRDS(source_file)

check(
  is.matrix(source),
  "Authoritative expression object is a matrix"
)

check(
  nrow(source) == 26468,
  "Authoritative expression contains 26,468 genes"
)

check(
  ncol(source) == 73,
  "Authoritative expression contains 73 samples"
)

check(
  !is.null(rownames(source)),
  "Authoritative matrix has gene identifiers"
)

check(
  !is.null(colnames(source)),
  "Authoritative matrix has sample identifiers"
)

sample_ids <- as.character(
  colnames(source)
)

source_names <- as.character(
  rownames(source)
)

check(
  length(unique(sample_ids)) == 73,
  "All 73 source sample identifiers are unique"
)

cat(
  "\n[INFO] Loaded authoritative RPKM: ",
  nrow(source),
  " x ",
  ncol(source),
  "\n",
  sep = ""
)


# ============================================================
# 2. Load FINAL gene decision table
# ============================================================

decision <- fread(
  decision_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("")
)

required_cols <- c(
  "SourceRow",
  "SourceName",
  "SourceGeneID",
  "ResolveStatus",
  "GeneIDChain",
  "FinalGeneID",
  "FinalSymbol",
  "FinalDecision",
  "Included"
)

check(
  all(required_cols %in% names(decision)),
  "FINAL decision table contains all required columns"
)

check(
  nrow(decision) == 26468,
  "FINAL decision table contains 26,468 source rows"
)

decision[, SourceRow := as.integer(SourceRow)]
decision[, SourceName := as.character(SourceName)]
decision[, SourceGeneID := as.character(SourceGeneID)]
decision[, ResolveStatus := as.character(ResolveStatus)]
decision[, GeneIDChain := as.character(GeneIDChain)]
decision[, FinalGeneID := as.character(FinalGeneID)]
decision[, FinalSymbol := as.character(FinalSymbol)]
decision[, FinalDecision := as.character(FinalDecision)]


if (is.logical(decision$Included)) {

  decision[, IncludedLogical := Included]

} else {

  decision[
    ,
    IncludedLogical :=
      toupper(trimws(as.character(Included))) == "TRUE"
  ]
}


check(
  all(
    decision$SourceRow ==
      seq_len(nrow(decision))
  ),
  "SourceRow exactly follows original source row order"
)

check(
  sum(decision$IncludedLogical) == 24656,
  "FINAL included gene count = 24,656"
)

check(
  sum(!decision$IncludedLogical) == 1812,
  "FINAL excluded gene count = 1,812"
)


# ============================================================
# 3. Validate source identity alignment
# ============================================================

check(
  same_text(
    decision$SourceName,
    source_names
  ),
  "FINAL decision table matches authoritative source row names"
)

source_geneids <- sub(
  "^GeneID:",
  "",
  source_names
)

check(
  same_text(
    decision$SourceGeneID,
    source_geneids
  ),
  "FINAL SourceGeneID matches authoritative source GeneID"
)


# ============================================================
# 4. Extract FINAL keep set
# ============================================================

keep_idx <- which(
  decision$IncludedLogical
)

keep <- decision[
  keep_idx
]

check(
  length(keep_idx) == 24656,
  "24,656 source rows selected"
)

check(
  all(
    keep$SourceRow ==
      keep_idx
  ),
  "Selected source rows exactly match FINAL SourceRow"
)

check(
  !anyNA(keep$FinalGeneID),
  "No retained FinalGeneID is missing"
)

check(
  !anyNA(keep$FinalSymbol),
  "No retained FinalSymbol is missing"
)

check(
  all(nzchar(keep$FinalGeneID)),
  "All retained FinalGeneID values are non-empty"
)

check(
  all(nzchar(keep$FinalSymbol)),
  "All retained FinalSymbol values are non-empty"
)

check(
  uniqueN(keep$FinalGeneID) == 24656,
  "FinalGeneID values are globally unique"
)

check(
  uniqueN(keep$FinalSymbol) == 24656,
  "FinalSymbol values are globally unique"
)

check(
  !any(grepl("^REVIEW_", keep$FinalDecision)),
  "No retained row requires further review"
)


# ============================================================
# 5. Revalidate frozen decision composition
# ============================================================

decision_summary <- decision[
  ,
  .N,
  by = FinalDecision
][order(-N)]

cat("\n============================================================\n")
cat("FINAL GENE DECISION SUMMARY\n")
cat("============================================================\n")

print(decision_summary)


check(
  decision[
    FinalDecision == "KEEP_DIRECT_CURRENT",
    .N
  ] == 24555,
  "KEEP_DIRECT_CURRENT = 24,555"
)

check(
  decision[
    FinalDecision == "KEEP_HISTORY_RESCUE_UNIQUE",
    .N
  ] == 101,
  "KEEP_HISTORY_RESCUE_UNIQUE = 101"
)

check(
  decision[
    FinalDecision == "EXCLUDE_RETIRED_NO_REPLACEMENT",
    .N
  ] == 1528,
  "EXCLUDE_RETIRED_NO_REPLACEMENT = 1,528"
)

check(
  decision[
    FinalDecision == "EXCLUDE_HISTORY_TARGET_ALREADY_PRESENT",
    .N
  ] == 278,
  "EXCLUDE_HISTORY_TARGET_ALREADY_PRESENT = 278"
)

check(
  decision[
    FinalDecision == "EXCLUDE_HISTORY_TARGET_COLLISION",
    .N
  ] == 6,
  "EXCLUDE_HISTORY_TARGET_COLLISION = 6"
)


# ============================================================
# 6. History-rescue audit
# ============================================================

history_keep <- keep[
  FinalDecision ==
    "KEEP_HISTORY_RESCUE_UNIQUE"
]

cat("\n============================================================\n")
cat("HISTORY RESCUE AUDIT\n")
cat("============================================================\n")

check(
  nrow(history_keep) == 101,
  "Unique history-rescued rows = 101"
)

check(
  all(
    history_keep$ResolveStatus ==
      "HISTORY_TO_CURRENT"
  ),
  "All retained history rescues resolve to current NCBI genes"
)

check(
  all(
    history_keep$SourceGeneID !=
      history_keep$FinalGeneID
  ),
  "All history rescues represent genuine GeneID replacements"
)


# ============================================================
# 7. Extract retained ORIGINAL source rows
# ============================================================

source_subset <- source[
  keep_idx,
  ,
  drop = FALSE
]

# Normalize dimnames to simple character vectors.
rownames(source_subset) <- source_names[keep_idx]
colnames(source_subset) <- sample_ids

check(
  nrow(source_subset) == 24656,
  "Retained source expression has 24,656 rows"
)

check(
  ncol(source_subset) == 73,
  "Retained source expression has 73 columns"
)

check(
  same_text(
    rownames(source_subset),
    keep$SourceName
  ),
  "Retained source rows exactly match FINAL keep set"
)

check(
  same_text(
    colnames(source_subset),
    sample_ids
  ),
  "Retained source sample order exactly matches authoritative matrix"
)


# ============================================================
# 8. Freeze numeric payload before renaming
# ============================================================

numeric_before <- as.vector(
  source_subset
)

check(
  length(numeric_before) ==
    24656 * 73,
  "Expected number of retained RPKM values present"
)


# ============================================================
# 9. Build canonical expression
#
# Numeric matrix is copied unchanged.
# Only row identity is replaced.
# ============================================================

canonical <- source_subset

rownames(canonical) <- as.character(
  keep$FinalSymbol
)

colnames(canonical) <- sample_ids


check(
  same_text(
    rownames(canonical),
    keep$FinalSymbol
  ),
  "Canonical row names equal FINAL NCBI symbols"
)

check(
  same_text(
    colnames(canonical),
    sample_ids
  ),
  "Canonical sample order preserved"
)


# ============================================================
# 10. Exact numeric preservation
# ============================================================

numeric_after <- as.vector(
  canonical
)

check(
  length(numeric_before) ==
    length(numeric_after),
  "Numeric payload length unchanged"
)

check(
  identical(
    numeric_before,
    numeric_after
  ),
  "Every retained RPKM value is exactly unchanged"
)


# ============================================================
# 11. Matrix integrity
# ============================================================

check(
  nrow(canonical) == 24656,
  "Canonical matrix contains 24,656 genes"
)

check(
  ncol(canonical) == 73,
  "Canonical matrix contains 73 samples"
)

check(
  uniqueN(rownames(canonical)) == 24656,
  "Canonical gene symbols are unique"
)

check(
  uniqueN(colnames(canonical)) == 73,
  "Canonical sample identifiers are unique"
)

check(
  sum(is.na(canonical)) == 0,
  "Canonical matrix contains no NA"
)

check(
  sum(!is.finite(canonical)) == 0,
  "Canonical matrix contains no non-finite values"
)

check(
  sum(canonical < 0) == 0,
  "Canonical matrix contains no negative RPKM values"
)


# ============================================================
# 12. Explicit history-rescue numeric forensic
# ============================================================

history_keep_positions <- which(
  keep$FinalDecision ==
    "KEEP_HISTORY_RESCUE_UNIQUE"
)

history_source_rows <- keep$SourceRow[
  history_keep_positions
]

history_source_matrix <- source[
  history_source_rows,
  ,
  drop = FALSE
]

# Normalize sample dimnames only.
colnames(history_source_matrix) <- sample_ids

history_canonical_matrix <- canonical[
  history_keep_positions,
  ,
  drop = FALSE
]

check(
  identical(
    as.vector(history_source_matrix),
    as.vector(history_canonical_matrix)
  ),
  "All 101 history-rescued genes retain original source-row RPKM exactly"
)


# ============================================================
# 13. Feature manifest
# ============================================================

feature_manifest <- data.table(
  CanonicalRow =
    seq_len(nrow(canonical)),

  SourceRow =
    keep$SourceRow,

  SourceName =
    keep$SourceName,

  SourceGeneID =
    keep$SourceGeneID,

  ResolveStatus =
    keep$ResolveStatus,

  GeneIDChain =
    keep$GeneIDChain,

  FinalGeneID =
    keep$FinalGeneID,

  FinalSymbol =
    keep$FinalSymbol,

  FinalDecision =
    keep$FinalDecision
)


check(
  same_text(
    feature_manifest$FinalSymbol,
    rownames(canonical)
  ),
  "Feature manifest aligns exactly with canonical matrix"
)


# ============================================================
# 14. Excluded-feature manifest
# ============================================================

excluded_manifest <- decision[
  IncludedLogical == FALSE,
  .(
    SourceRow,
    SourceName,
    SourceGeneID,
    ResolveStatus,
    GeneIDChain,
    FinalGeneID,
    FinalSymbol,
    FinalDecision
  )
]

check(
  nrow(excluded_manifest) == 1812,
  "Excluded-feature manifest contains 1,812 rows"
)


# ============================================================
# 15. Matrix-wide audit
# ============================================================

matrix_audit <- data.table(
  Metric = c(
    "Genes",
    "Samples",
    "Min",
    "Max",
    "Zero_N",
    "Zero_fraction",
    "Positive_N",
    "Positive_fraction",
    "Negative_N",
    "NA_N",
    "Nonfinite_N",
    "Unique_FinalGeneID",
    "Unique_FinalSymbol",
    "Unique_samples",
    "Direct_current_rows",
    "History_rescue_rows",
    "Values_exactly_unchanged"
  ),

  Value = c(
    as.character(nrow(canonical)),
    as.character(ncol(canonical)),
    format(min(canonical), digits = 16),
    format(max(canonical), digits = 16),
    as.character(sum(canonical == 0)),
    format(mean(canonical == 0), digits = 16),
    as.character(sum(canonical > 0)),
    format(mean(canonical > 0), digits = 16),
    as.character(sum(canonical < 0)),
    as.character(sum(is.na(canonical))),
    as.character(sum(!is.finite(canonical))),
    as.character(uniqueN(keep$FinalGeneID)),
    as.character(uniqueN(keep$FinalSymbol)),
    as.character(uniqueN(sample_ids)),
    as.character(
      sum(
        keep$FinalDecision ==
          "KEEP_DIRECT_CURRENT"
      )
    ),
    as.character(
      sum(
        keep$FinalDecision ==
          "KEEP_HISTORY_RESCUE_UNIQUE"
      )
    ),
    as.character(
      identical(
        numeric_before,
        numeric_after
      )
    )
  )
)


cat("\n============================================================\n")
cat("CANONICAL MATRIX AUDIT\n")
cat("============================================================\n")

print(matrix_audit)


# ============================================================
# 16. Per-sample audit
# ============================================================

sample_audit <- data.table(
  GSM = sample_ids,

  Min = apply(
    canonical,
    2,
    min
  ),

  Median = apply(
    canonical,
    2,
    median
  ),

  Mean = colMeans(
    canonical
  ),

  Max = apply(
    canonical,
    2,
    max
  ),

  SD = apply(
    canonical,
    2,
    sd
  ),

  DetectedGenes = colSums(
    canonical > 0
  ),

  ZeroGenes = colSums(
    canonical == 0
  ),

  ZeroFraction = colMeans(
    canonical == 0
  )
)


sample_summary <- sample_audit[
  ,
  .(
    Samples = .N,

    Median_detected_genes =
      median(DetectedGenes),

    Min_detected_genes =
      min(DetectedGenes),

    Max_detected_genes =
      max(DetectedGenes),

    Median_zero_fraction =
      median(ZeroFraction),

    Min_zero_fraction =
      min(ZeroFraction),

    Max_zero_fraction =
      max(ZeroFraction)
  )
]


cat("\n============================================================\n")
cat("PER-SAMPLE EXPRESSION SUMMARY\n")
cat("============================================================\n")

print(sample_summary)


# ============================================================
# 17. Save outputs
# ============================================================

canonical_rds <- file.path(
  out_dir,
  "GSE72819_canonical_RPKM_24656x73.rds"
)

saveRDS(
  canonical,
  canonical_rds,
  compress = FALSE
)


canonical_tsv <- file.path(
  out_dir,
  "GSE72819_canonical_RPKM_24656x73.tsv.gz"
)

canonical_dt <- data.table(
  GeneSymbol = rownames(canonical)
)

canonical_dt <- cbind(
  canonical_dt,
  as.data.table(
    canonical,
    keep.rownames = FALSE
  )
)

setnames(
  canonical_dt,
  c(
    "GeneSymbol",
    sample_ids
  )
)

fwrite(
  canonical_dt,
  canonical_tsv,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  feature_manifest,
  file.path(
    out_dir,
    "GSE72819_canonical_feature_manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  excluded_manifest,
  file.path(
    out_dir,
    "GSE72819_excluded_feature_manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  matrix_audit,
  file.path(
    out_dir,
    "GSE72819_canonical_matrix_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  sample_audit,
  file.path(
    out_dir,
    "GSE72819_canonical_per_sample_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  sample_summary,
  file.path(
    out_dir,
    "GSE72819_canonical_sample_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  decision_summary,
  file.path(
    out_dir,
    "GSE72819_FINAL_decision_summary_revalidated.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 18. Exact RDS round-trip
# ============================================================

reload <- readRDS(
  canonical_rds
)

# Normalize dimnames for robust comparison.
rownames(reload) <- as.character(
  rownames(reload)
)

colnames(reload) <- as.character(
  colnames(reload)
)

check(
  identical(
    as.vector(reload),
    as.vector(canonical)
  ),
  "Canonical RDS numeric round-trip is exact"
)

check(
  same_text(
    rownames(reload),
    rownames(canonical)
  ),
  "Canonical RDS gene order round-trip is exact"
)

check(
  same_text(
    colnames(reload),
    colnames(canonical)
  ),
  "Canonical RDS sample order round-trip is exact"
)


# ============================================================
# 19. Hash artifacts
# ============================================================

md5 <- tools::md5sum(
  c(
    source_file,
    decision_file,
    canonical_rds
  )
)

artifact_md5 <- data.table(
  File = names(md5),
  MD5 = unname(md5)
)

fwrite(
  artifact_md5,
  file.path(
    out_dir,
    "GSE72819_stage3_artifact_md5.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 20. Final gate
# ============================================================

critical_pass <- all(c(
  nrow(canonical) == 24656,
  ncol(canonical) == 73,
  uniqueN(keep$FinalGeneID) == 24656,
  uniqueN(rownames(canonical)) == 24656,
  uniqueN(colnames(canonical)) == 73,
  sum(is.na(canonical)) == 0,
  sum(!is.finite(canonical)) == 0,
  sum(canonical < 0) == 0,
  identical(
    numeric_before,
    numeric_after
  ),
  identical(
    as.vector(source_subset),
    as.vector(canonical)
  ),
  same_text(
    colnames(canonical),
    sample_ids
  ),
  same_text(
    rownames(canonical),
    keep$FinalSymbol
  )
))


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] 24,656 canonical genes retained\n")
  cat("[PASS] 73 samples retained\n")
  cat("[PASS] 24,656 unique FinalGeneID values\n")
  cat("[PASS] 24,656 unique FinalSymbol values\n")
  cat("[PASS] 101 historical GeneID rescues retained\n")
  cat("[PASS] Historical rescues preserve original source RPKM\n")
  cat("[PASS] Sample identity/order preserved\n")
  cat("[PASS] Expression values preserved exactly\n")
  cat("[PASS] No normalization performed\n")
  cat("[PASS] No log transformation performed\n")
  cat("[PASS] No batch correction performed\n")
  cat("[PASS] No aggregation performed\n")
  cat("[PASS] No imputation performed\n")
  cat("[PASS] Canonical RDS round-trip exact\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE3_CANONICAL_EXPRESSION_BUILT\n")

} else {

  cat("[REVIEW] Stage 3 final validation did not pass\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE3_CANONICAL_EXPRESSION_BUILT\n")
}

cat("============================================================\n")
