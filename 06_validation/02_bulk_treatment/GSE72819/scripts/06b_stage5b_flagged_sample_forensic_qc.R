suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE72819 Stage 5B
# Lightweight forensic QC for 3 diagnostic low-correlation flags
#
# No sample exclusion is performed here.
# No expression values are modified.
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

canonical_file <- file.path(
  root,
  "prepared/03_canonical_expression",
  "GSE72819_canonical_RPKM_24656x73.rds"
)

counts_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE72819_source_counts_26468x73.rds"
)

source_rpkm_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE72819_authoritative_RPKM_26468x73.rds"
)

meta_file <- file.path(
  root,
  "prepared/04_sample_metadata",
  "GSE72819_canonical_sample_metadata_73.rds"
)

out_dir <- file.path(
  root,
  "prepared/05_sample_expression_qc"
)

flagged <- c(
  "GSM1872905",
  "GSM1872911",
  "GSM1872914"
)

check <- function(x, msg) {
  if (!isTRUE(x)) {
    stop(paste0("[FAIL] ", msg), call. = FALSE)
  }
  cat("[PASS] ", msg, "\n", sep = "")
}

cat("============================================================\n")
cat("GSE72819 Stage 5B: flagged-sample forensic QC\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load inputs
# ============================================================

expr <- readRDS(canonical_file)
counts <- readRDS(counts_file)
source_rpkm <- readRDS(source_rpkm_file)
meta <- as.data.table(readRDS(meta_file))

check(is.matrix(expr), "Canonical RPKM loaded")
check(is.matrix(counts), "Source counts loaded")
check(is.matrix(source_rpkm), "Source RPKM loaded")

check(
  all(flagged %in% colnames(expr)),
  "All 3 flagged samples occur in canonical expression"
)

check(
  all(flagged %in% colnames(counts)),
  "All 3 flagged samples occur in source counts"
)

check(
  all(flagged %in% colnames(source_rpkm)),
  "All 3 flagged samples occur in source RPKM"
)

check(
  all(flagged %in% meta$GSM),
  "All 3 flagged samples occur in canonical metadata"
)


# Reorder source matrices to canonical GSM order.
gsm <- as.character(colnames(expr))

counts <- counts[
  ,
  match(gsm, colnames(counts)),
  drop = FALSE
]

source_rpkm <- source_rpkm[
  ,
  match(gsm, colnames(source_rpkm)),
  drop = FALSE
]

colnames(counts) <- gsm
colnames(source_rpkm) <- gsm

check(
  identical(colnames(counts), gsm),
  "Source counts aligned to canonical sample order"
)

check(
  identical(colnames(source_rpkm), gsm),
  "Source RPKM aligned to canonical sample order"
)


# ============================================================
# 2. QC-only correlations
# ============================================================

log_expr <- log2(expr + 1)

median_profile <- apply(
  log_expr,
  1,
  median
)

pearson_to_median <- as.numeric(
  cor(
    log_expr,
    median_profile,
    method = "pearson"
  )
)

spearman_to_median <- as.numeric(
  cor(
    log_expr,
    median_profile,
    method = "spearman"
  )
)

pearson_pair <- cor(
  log_expr,
  method = "pearson"
)

spearman_pair <- cor(
  log_expr,
  method = "spearman"
)

diag(pearson_pair) <- -Inf
diag(spearman_pair) <- -Inf

pearson_nn <- apply(
  pearson_pair,
  1,
  max
)

spearman_nn <- apply(
  spearman_pair,
  1,
  max
)


# ============================================================
# 3. Source-level QC metrics
# ============================================================

count_sum <- colSums(counts)
count_detected <- colSums(counts > 0)

source_rpkm_sum <- colSums(source_rpkm)
source_rpkm_detected <- colSums(source_rpkm > 0)
source_rpkm_zero_fraction <- colMeans(source_rpkm == 0)


# ============================================================
# 4. Percentile helper
# ============================================================

pct <- function(x) {

  r <- rank(
    x,
    ties.method = "average"
  )

  if (length(x) == 1L) {
    return(1)
  }

  (r - 1) / (length(x) - 1)
}


count_sum_pct <- pct(count_sum)
count_detected_pct <- pct(count_detected)
pearson_pct <- pct(pearson_to_median)
spearman_pct <- pct(spearman_to_median)


# ============================================================
# 5. Build targeted forensic table
# ============================================================

idx <- match(flagged, gsm)

meta_flagged <- meta[
  match(flagged, GSM)
]

result <- data.table(
  GSM = flagged,

  PatientID =
    as.character(meta_flagged$PatientID),

  Run =
    as.character(meta_flagged$Run),

  Lane =
    as.character(meta_flagged$Lane),

  ResponseWeek10 =
    as.character(meta_flagged$ResponseWeek10),

  PriorAntiTNF =
    as.character(meta_flagged$PriorAntiTNF),

  CountSum =
    count_sum[idx],

  CountSumPercentile =
    count_sum_pct[idx],

  CountDetectedGenes =
    count_detected[idx],

  CountDetectedPercentile =
    count_detected_pct[idx],

  SourceRPKMSum =
    source_rpkm_sum[idx],

  SourceRPKMDetectedGenes =
    source_rpkm_detected[idx],

  SourceRPKMZeroFraction =
    source_rpkm_zero_fraction[idx],

  PearsonToMedian =
    pearson_to_median[idx],

  PearsonPercentile =
    pearson_pct[idx],

  SpearmanToMedian =
    spearman_to_median[idx],

  SpearmanPercentile =
    spearman_pct[idx],

  PearsonNearestNeighbor =
    pearson_nn[idx],

  SpearmanNearestNeighbor =
    spearman_nn[idx]
)


# ============================================================
# 6. Simple technical-warning screen
#
# Deliberately conservative:
# only obvious source-level extremes are warnings.
# Correlation alone does NOT trigger technical exclusion.
# ============================================================

result[
  ,
  SourceCountWarning :=
    CountSumPercentile <= 0.01 |
    CountSumPercentile >= 0.99 |
    CountDetectedPercentile <= 0.01 |
    CountDetectedPercentile >= 0.99
]

result[
  ,
  StrongCorrelationWarning :=
    (
      PearsonToMedian < 0.85 &
      SpearmanToMedian < 0.85
    )
]

result[
  ,
  TechnicalConcern :=
    SourceCountWarning |
    StrongCorrelationWarning
]


# ============================================================
# 7. Cohort reference summary
# ============================================================

cohort_summary <- data.table(
  Metric = c(
    "CountSum_median",
    "CountSum_min",
    "CountSum_max",
    "CountDetected_median",
    "CountDetected_min",
    "CountDetected_max",
    "PearsonToMedian_median",
    "PearsonToMedian_min",
    "SpearmanToMedian_median",
    "SpearmanToMedian_min",
    "PearsonNN_median",
    "PearsonNN_min",
    "SpearmanNN_median",
    "SpearmanNN_min"
  ),

  Value = c(
    median(count_sum),
    min(count_sum),
    max(count_sum),

    median(count_detected),
    min(count_detected),
    max(count_detected),

    median(pearson_to_median),
    min(pearson_to_median),

    median(spearman_to_median),
    min(spearman_to_median),

    median(pearson_nn),
    min(pearson_nn),

    median(spearman_nn),
    min(spearman_nn)
  )
)


cat("============================================================\n")
cat("FLAGGED SAMPLE FORENSIC TABLE\n")
cat("============================================================\n")

print(result)

cat("\n============================================================\n")
cat("COHORT REFERENCE\n")
cat("============================================================\n")

print(cohort_summary)


# ============================================================
# 8. Save
# ============================================================

fwrite(
  result,
  file.path(
    out_dir,
    "GSE72819_stage5b_flagged_sample_forensic_qc.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  cohort_summary,
  file.path(
    out_dir,
    "GSE72819_stage5b_cohort_reference.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 9. Final status
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] 3 diagnostic samples inspected\n")
cat("[PASS] Source counts inspected\n")
cat("[PASS] Source RPKM inspected\n")
cat("[PASS] Pearson and Spearman correlations compared\n")
cat("[PASS] No expression data modified\n")

cat(
  "[INFO] Samples with targeted technical concern = ",
  sum(result$TechnicalConcern),
  "\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE72819_STAGE5B_FLAGGED_SAMPLE_FORENSIC_QC\n")

cat("============================================================\n")
