suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 5
#
# Sample-level expression QC
#
# Purpose:
#   - detect gross technical anomalies
#   - inspect sample distributions
#   - identify duplicate / near-duplicate profiles
#   - quantify sample similarity to cohort reference
#
# IMPORTANT:
#   - diagnostic only
#   - NO automatic sample exclusion
#   - NO normalization
#   - NO transformation
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

expr_file <- file.path(
  base,
  "prepared/03_canonical_expression",
  "GSE93624_canonical_expression_13323x245.rds"
)

meta_file <- file.path(
  base,
  "prepared/04_sample_metadata",
  "GSE93624_canonical_sample_metadata.tsv"
)

out_dir <- file.path(
  base,
  "prepared/05_sample_expression_qc"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE93624 sample-level expression QC\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load canonical expression + metadata
# ============================================================

obj <- readRDS(expr_file)
expr <- obj$expression

meta <- fread(
  meta_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  is.matrix(expr),
  nrow(expr) == 13323,
  ncol(expr) == 245,
  nrow(meta) == 245,
  identical(
    colnames(expr),
    meta$ExpressionSampleID
  ),
  !anyNA(expr),
  all(is.finite(expr))
)

cat(
  "[PASS] Expression:",
  nrow(expr),
  "x",
  ncol(expr),
  "\n"
)

cat(
  "[PASS] Metadata aligned exactly to expression columns\n"
)


# ============================================================
# 2. Per-sample expression metrics
# ============================================================

sample_qc <- data.table(
  ExpressionSampleID =
    colnames(expr),

  Min =
    apply(expr, 2, min),

  Q01 =
    apply(
      expr,
      2,
      quantile,
      probs = 0.01,
      names = FALSE,
      type = 7
    ),

  Q05 =
    apply(
      expr,
      2,
      quantile,
      probs = 0.05,
      names = FALSE,
      type = 7
    ),

  Median =
    apply(expr, 2, median),

  Mean =
    colMeans(expr),

  Q95 =
    apply(
      expr,
      2,
      quantile,
      probs = 0.95,
      names = FALSE,
      type = 7
    ),

  Q99 =
    apply(
      expr,
      2,
      quantile,
      probs = 0.99,
      names = FALSE,
      type = 7
    ),

  Max =
    apply(expr, 2, max),

  SD =
    apply(expr, 2, sd),

  IQR =
    apply(expr, 2, IQR),

  MAD =
    apply(
      expr,
      2,
      mad
    ),

  Negative_N =
    colSums(expr < 0),

  Zero_N =
    colSums(expr == 0),

  Positive_N =
    colSums(expr > 0),

  NegativeFraction =
    colMeans(expr < 0),

  ZeroFraction =
    colMeans(expr == 0),

  PositiveFraction =
    colMeans(expr > 0)
)


# ============================================================
# 3. Cohort-reference expression profile
#
# Reference = per-gene median across 245 samples.
# No biological interpretation is attached to this profile.
# ============================================================

reference_profile <- apply(
  expr,
  1,
  median
)

stopifnot(
  length(reference_profile) == 13323
)


sample_qc[
  ,
  PearsonToMedianProfile :=
    as.numeric(
      cor(
        expr,
        reference_profile,
        method = "pearson"
      )
    )
]

sample_qc[
  ,
  SpearmanToMedianProfile :=
    as.numeric(
      cor(
        expr,
        reference_profile,
        method = "spearman"
      )
    )
]

sample_qc[
  ,
  MeanAbsDeviationFromMedianProfile :=
    colMeans(
      abs(
        expr -
        reference_profile
      )
    )
]


# ============================================================
# 4. Pairwise sample correlations
# ============================================================

cat(
  "\nComputing pairwise sample correlations...\n"
)

pearson_cor <- cor(
  expr,
  method = "pearson"
)

spearman_cor <- cor(
  expr,
  method = "spearman"
)

stopifnot(
  all(
    dim(pearson_cor) ==
      c(245, 245)
  ),
  all(
    dim(spearman_cor) ==
      c(245, 245)
  )
)


# ============================================================
# 5. Nearest-neighbour similarity
# ============================================================

pearson_nn <- pearson_cor
spearman_nn <- spearman_cor

diag(pearson_nn) <- NA_real_
diag(spearman_nn) <- NA_real_


get_top_two <- function(x) {

  x <- as.numeric(x)

  x[
    is.na(x)
  ] <- -Inf

  ord <- order(
    x,
    decreasing = TRUE
  )

  c(
    first = ord[1],
    second = ord[2]
  )
}


top_idx <- t(
  apply(
    pearson_nn,
    1,
    get_top_two
  )
)

sample_qc[
  ,
  NearestNeighbor :=
    colnames(expr)[
      top_idx[, "first"]
    ]
]

sample_qc[
  ,
  PearsonNearest :=
    vapply(
      seq_len(nrow(sample_qc)),
      function(i) {
        pearson_nn[
          i,
          top_idx[i, "first"]
        ]
      },
      numeric(1)
    )
]

sample_qc[
  ,
  PearsonSecondNearest :=
    vapply(
      seq_len(nrow(sample_qc)),
      function(i) {
        pearson_nn[
          i,
          top_idx[i, "second"]
        ]
      },
      numeric(1)
    )
]

sample_qc[
  ,
  PearsonNNGap :=
    PearsonNearest -
    PearsonSecondNearest
]

sample_qc[
  ,
  SpearmanNearest :=
    vapply(
      seq_len(nrow(sample_qc)),
      function(i) {

        j <- top_idx[
          i,
          "first"
        ]

        spearman_nn[
          i,
          j
        ]

      },
      numeric(1)
    )
]


# ============================================================
# 6. High-similarity pair forensic
#
# Only very highly correlated pairs receive expensive
# element-wise difference checks.
# ============================================================

pair_idx <- which(
  upper.tri(
    pearson_cor
  ) &
    pearson_cor >= 0.9995,
  arr.ind = TRUE
)

if (nrow(pair_idx) > 0) {

  high_pairs <- rbindlist(
    lapply(
      seq_len(nrow(pair_idx)),
      function(k) {

        i <- pair_idx[k, 1]
        j <- pair_idx[k, 2]

        x <- expr[, i]
        y <- expr[, j]

        data.table(
          Sample1 =
            colnames(expr)[i],

          Sample2 =
            colnames(expr)[j],

          Pearson =
            pearson_cor[i, j],

          Spearman =
            spearman_cor[i, j],

          ExactDuplicate =
            identical(
              as.numeric(x),
              as.numeric(y)
            ),

          MeanAbsDiff =
            mean(
              abs(x - y)
            ),

          MaxAbsDiff =
            max(
              abs(x - y)
            )
        )
      }
    )
  )

} else {

  high_pairs <- data.table(
    Sample1 = character(),
    Sample2 = character(),
    Pearson = numeric(),
    Spearman = numeric(),
    ExactDuplicate = logical(),
    MeanAbsDiff = numeric(),
    MaxAbsDiff = numeric()
  )
}


# ============================================================
# 7. Exact duplicate audit
# ============================================================

exact_pairs <- high_pairs[
  ExactDuplicate == TRUE
]

cat("\n============================================================\n")
cat("DUPLICATE / NEAR-DUPLICATE AUDIT\n")
cat("============================================================\n")

cat(
  "Pairs with Pearson >= 0.9995:",
  nrow(high_pairs),
  "\n"
)

cat(
  "Exact duplicate pairs:",
  nrow(exact_pairs),
  "\n"
)

if (nrow(high_pairs) > 0) {

  print(
    high_pairs[
      order(
        -Pearson,
        MaxAbsDiff
      )
    ]
  )
}


# ============================================================
# 8. Robust diagnostic flags
#
# Flags are intentionally conservative.
# They are NOT sample-exclusion rules.
# ============================================================

robust_z <- function(x) {

  med <- median(
    x,
    na.rm = TRUE
  )

  s <- mad(
    x,
    center = med,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (
    !is.finite(s) ||
    s == 0
  ) {

    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  (x - med) / s
}


sample_qc[
  ,
  rz_Mean :=
    robust_z(Mean)
]

sample_qc[
  ,
  rz_Median :=
    robust_z(Median)
]

sample_qc[
  ,
  rz_SD :=
    robust_z(SD)
]

sample_qc[
  ,
  rz_MAD :=
    robust_z(MAD)
]

sample_qc[
  ,
  rz_PearsonToMedian :=
    robust_z(
      PearsonToMedianProfile
    )
]

sample_qc[
  ,
  rz_DeviationFromMedian :=
    robust_z(
      MeanAbsDeviationFromMedianProfile
    )
]


sample_qc[
  ,
  Flag_ExtremeMean :=
    !is.na(rz_Mean) &
    abs(rz_Mean) > 5
]

sample_qc[
  ,
  Flag_ExtremeMedian :=
    !is.na(rz_Median) &
    abs(rz_Median) > 5
]

sample_qc[
  ,
  Flag_ExtremeSD :=
    !is.na(rz_SD) &
    abs(rz_SD) > 5
]

sample_qc[
  ,
  Flag_ExtremeMAD :=
    !is.na(rz_MAD) &
    abs(rz_MAD) > 5
]

sample_qc[
  ,
  Flag_LowCorrelation :=
    !is.na(rz_PearsonToMedian) &
    rz_PearsonToMedian < -5
]

sample_qc[
  ,
  Flag_HighProfileDeviation :=
    !is.na(rz_DeviationFromMedian) &
    rz_DeviationFromMedian > 5
]


flag_columns <- c(
  "Flag_ExtremeMean",
  "Flag_ExtremeMedian",
  "Flag_ExtremeSD",
  "Flag_ExtremeMAD",
  "Flag_LowCorrelation",
  "Flag_HighProfileDeviation"
)

sample_qc[
  ,
  DiagnosticFlagN :=
    rowSums(
      .SD,
      na.rm = TRUE
    ),
  .SDcols = flag_columns
]

sample_qc[
  ,
  AnyDiagnosticFlag :=
    DiagnosticFlagN > 0
]


# ============================================================
# 9. Attach biological metadata
# ============================================================

meta_keep <- meta[
  ,
  .(
    ExpressionSampleID,
    GSM,
    SampleTitle,
    Diagnosis,
    Progression3yr,
    Sex,
    Tissue,
    AgeAtDiagnosis,
    ParisAge,
    Ancestry
  )
]

sample_qc <- merge(
  sample_qc,
  meta_keep,
  by = "ExpressionSampleID",
  all.x = TRUE,
  sort = FALSE
)

sample_qc <- sample_qc[
  match(
    colnames(expr),
    ExpressionSampleID
  )
]

stopifnot(
  identical(
    sample_qc$ExpressionSampleID,
    colnames(expr)
  )
)


# ============================================================
# 10. Summary
# ============================================================

cat("\n============================================================\n")
cat("SAMPLE QC SUMMARY\n")
cat("============================================================\n")

summary_metrics <- data.table(
  Metric = c(
    "Samples",
    "DiagnosticFlaggedSamples",
    "ExactDuplicatePairs",
    "HighSimilarityPairs_Pearson_ge_0.9995",
    "Median_PearsonToMedianProfile",
    "Min_PearsonToMedianProfile",
    "Median_NearestNeighborPearson",
    "Min_NearestNeighborPearson"
  ),

  Value = c(
    245,
    sum(
      sample_qc$AnyDiagnosticFlag
    ),
    nrow(exact_pairs),
    nrow(high_pairs),
    median(
      sample_qc$PearsonToMedianProfile
    ),
    min(
      sample_qc$PearsonToMedianProfile
    ),
    median(
      sample_qc$PearsonNearest
    ),
    min(
      sample_qc$PearsonNearest
    )
  )
)

print(
  summary_metrics
)


cat("\n============================================================\n")
cat("DIAGNOSTIC FLAG SUMMARY\n")
cat("============================================================\n")

flag_summary <- rbindlist(
  lapply(
    flag_columns,
    function(v) {

      data.table(
        Flag =
          v,

        N =
          sum(
            sample_qc[[v]],
            na.rm = TRUE
          )
      )
    }
  )
)

print(
  flag_summary
)


flagged <- sample_qc[
  AnyDiagnosticFlag == TRUE
]

cat("\nFlagged samples:", nrow(flagged), "\n")

if (nrow(flagged) > 0) {

  print(
    flagged[
      ,
      .(
        ExpressionSampleID,
        GSM,
        Diagnosis,
        Progression3yr,
        DiagnosticFlagN,
        Mean,
        Median,
        SD,
        PearsonToMedianProfile,
        MeanAbsDeviationFromMedianProfile,
        NearestNeighbor,
        PearsonNearest
      )
    ]
  )
}


# ============================================================
# 11. Group-level descriptive QC
#
# Descriptive only; biological group differences are not
# interpreted as technical effects.
# ============================================================

group_qc <- sample_qc[
  ,
  .(
    N = .N,

    MeanExpression_median =
      median(Mean),

    MedianExpression_median =
      median(Median),

    SD_median =
      median(SD),

    PearsonToMedian_median =
      median(
        PearsonToMedianProfile
      ),

    NearestNeighborPearson_median =
      median(
        PearsonNearest
      ),

    Flagged_N =
      sum(
        AnyDiagnosticFlag
      )
  ),
  by = .(
    Diagnosis,
    Progression3yr
  )
]


cat("\n============================================================\n")
cat("GROUP-DESCRIPTIVE QC\n")
cat("============================================================\n")

print(
  group_qc
)


# ============================================================
# 12. Save
# ============================================================

fwrite(
  sample_qc,
  file.path(
    out_dir,
    "GSE93624_sample_level_expression_qc.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  summary_metrics,
  file.path(
    out_dir,
    "GSE93624_sample_qc_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  flag_summary,
  file.path(
    out_dir,
    "GSE93624_sample_qc_flag_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  flagged,
  file.path(
    out_dir,
    "GSE93624_sample_qc_flagged_samples.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  high_pairs,
  file.path(
    out_dir,
    "GSE93624_high_similarity_sample_pairs.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  group_qc,
  file.path(
    out_dir,
    "GSE93624_group_descriptive_expression_qc.tsv"
  ),
  sep = "\t",
  na = ""
)


saveRDS(
  list(
    pearson =
      pearson_cor,

    spearman =
      spearman_cor
  ),
  file.path(
    out_dir,
    "GSE93624_pairwise_sample_correlations.rds"
  ),
  compress = FALSE
)


# ============================================================
# 13. Final
# ============================================================

cat("\n============================================================\n")
cat("FINAL QC STATUS\n")
cat("============================================================\n")

cat("[PASS] All 245 samples evaluated\n")
cat("[PASS] No expression transformation performed\n")
cat("[PASS] QC flags are diagnostic only\n")

if (nrow(exact_pairs) == 0) {
  cat("[PASS] No exact duplicate sample profiles detected\n")
} else {
  cat(
    "[REVIEW] Exact duplicate pairs detected:",
    nrow(exact_pairs),
    "\n"
  )
}

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_SAMPLE_LEVEL_EXPRESSION_QC_COMPLETED\n")
cat("============================================================\n")
