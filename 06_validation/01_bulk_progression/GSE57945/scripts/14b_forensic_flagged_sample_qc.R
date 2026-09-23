suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE57945
# Stage 3B: flagged-sample forensic QC
#
# Purpose:
#   Investigate diagnostic QC flags from Stage 13.
#
# This script DOES NOT:
#   - remove samples
#   - normalize expression
#   - batch-correct expression
#   - modify the canonical expression matrix
#
# Main audits:
#   1. expression-composition metrics
#   2. exact duplicate samples
#   3. very-high-correlation sample pairs
#   4. nearest-neighbour expression correlation
#   5. within-disease PCA outlier reassessment
#   6. descriptive 254/68 provenance-segment comparison
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

expr_file <- file.path(
  base,
  "prepared/03_canonical_expression",
  "GSE57945_canonical_expression_34368x322.rds"
)

meta_file <- file.path(
  base,
  "prepared/04_sample_metadata",
  "GSE57945_canonical_sample_metadata_322.rds"
)

qc_file <- file.path(
  base,
  "prepared/05_expression_qc",
  "GSE57945_sample_expression_qc.tsv"
)

pca_file <- file.path(
  base,
  "prepared/05_expression_qc",
  "GSE57945_expression_PCA_scores.tsv"
)

out_dir <- file.path(
  base,
  "prepared/05_expression_qc"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE57945 flagged-sample forensic QC\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load inputs
# ============================================================

expr <- readRDS(expr_file)
meta <- as.data.table(readRDS(meta_file))

qc <- fread(
  qc_file,
  colClasses = NULL
)

pca <- fread(
  pca_file,
  colClasses = NULL
)

stopifnot(
  nrow(expr) == 34368,
  ncol(expr) == 322,
  nrow(meta) == 322,
  nrow(qc) == 322,
  nrow(pca) == 322,
  !anyDuplicated(colnames(expr)),
  !anyDuplicated(meta$SampleID),
  !anyDuplicated(qc$SampleID),
  !anyDuplicated(pca$SampleID),
  all(is.finite(expr)),
  all(expr >= 0)
)

sample_ids <- colnames(expr)

stopifnot(
  setequal(sample_ids, meta$SampleID),
  setequal(sample_ids, qc$SampleID),
  setequal(sample_ids, pca$SampleID)
)


# ------------------------------------------------------------
# Enforce exactly the same sample order everywhere
# ------------------------------------------------------------

meta <- meta[
  match(sample_ids, SampleID)
]

qc <- qc[
  match(sample_ids, SampleID)
]

pca <- pca[
  match(sample_ids, SampleID)
]

stopifnot(
  identical(meta$SampleID, sample_ids),
  identical(qc$SampleID, sample_ids),
  identical(pca$SampleID, sample_ids)
)

cat("[PASS] All four inputs aligned to identical 322-sample order\n")


# ============================================================
# 2. Required-column checks
# ============================================================

qc_required <- c(
  "SampleID",
  "GSM",
  "Disease",
  "DiseaseSubtype",
  "DetectedGenes",
  "ZeroFraction",
  "Log2TotalExpression",
  "UnivariateQCFlag",
  "PCAQCFlag",
  "AnyQCFlag"
)

pca_required <- c(
  "SampleID",
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5",
  "RobustPCDistance",
  "PCAQCFlag"
)

meta_required <- c(
  "SampleID",
  "Histopathology",
  "DeepUlcer",
  "ParisStage"
)

stopifnot(
  all(qc_required %in% names(qc)),
  all(pca_required %in% names(pca)),
  all(meta_required %in% names(meta))
)


# ============================================================
# 3. Confirm Stage 13 PCAQCFlag consistency
# ============================================================

stopifnot(
  identical(
    as.logical(qc$PCAQCFlag),
    as.logical(pca$PCAQCFlag)
  )
)

cat(
  "[PASS] Stage 13 PCAQCFlag identical between QC and PCA tables\n\n"
)


# ============================================================
# 4. Expression-composition metrics
# ============================================================

cat("Computing expression-composition metrics...\n")

total_expression <- colSums(expr)

max_expression <- apply(
  expr,
  2,
  max
)

top1_fraction <- max_expression / total_expression


# ------------------------------------------------------------
# Top-10 contribution
#
# Use full decreasing sort.
# This avoids the unsupported:
#   partial + decreasing=TRUE
# combination in sort.int().
# ------------------------------------------------------------

top10_fraction <- vapply(
  seq_len(ncol(expr)),
  function(j) {

    x <- expr[, j]

    total_x <- sum(x)

    if (
      !is.finite(total_x) ||
      total_x <= 0
    ) {
      return(NA_real_)
    }

    n_top <- min(
      10L,
      length(x)
    )

    top_values <- sort(
      x,
      decreasing = TRUE
    )[seq_len(n_top)]

    sum(top_values) / total_x
  },
  numeric(1)
)


composition_qc <- data.table(
  SampleID = sample_ids,
  Top1Fraction = top1_fraction,
  Top10Fraction = top10_fraction
)

stopifnot(
  !anyNA(composition_qc$Top1Fraction),
  !anyNA(composition_qc$Top10Fraction)
)

cat("[PASS] Expression-composition metrics calculated\n\n")


# ============================================================
# 5. Exact duplicate sample audit
# ============================================================

cat("Auditing exact duplicate sample profiles...\n")

# First use deterministic summaries to identify candidate groups.
# Only candidate groups are then compared value-by-value.

signature <- data.table(
  SampleID = sample_ids,

  DetectedGenes =
    colSums(expr > 0),

  SumExpression =
    colSums(expr),

  SumSquares =
    colSums(expr * expr),

  MaxExpression =
    max_expression
)

signature[
  ,
  SignatureKey := paste(
    DetectedGenes,
    sprintf("%.17g", SumExpression),
    sprintf("%.17g", SumSquares),
    sprintf("%.17g", MaxExpression),
    sep = "|"
  )
]

candidate_keys <- signature[
  ,
  .N,
  by = SignatureKey
][
  N > 1,
  SignatureKey
]

duplicate_list <- list()
dup_counter <- 0L

if (length(candidate_keys) > 0) {

  for (key in candidate_keys) {

    ids <- signature[
      SignatureKey == key,
      SampleID
    ]

    if (length(ids) < 2)
      next

    pairs <- combn(
      ids,
      2,
      simplify = FALSE
    )

    for (pair in pairs) {

      identical_profile <- identical(
        as.numeric(expr[, pair[1]]),
        as.numeric(expr[, pair[2]])
      )

      if (identical_profile) {

        dup_counter <- dup_counter + 1L

        duplicate_list[[dup_counter]] <- data.table(
          Sample1 = pair[1],
          Sample2 = pair[2],
          ExactDuplicate = TRUE
        )
      }
    }
  }
}

if (length(duplicate_list) > 0) {

  exact_dup <- rbindlist(
    duplicate_list
  )

} else {

  exact_dup <- data.table(
    Sample1 = character(),
    Sample2 = character(),
    ExactDuplicate = logical()
  )
}

cat(
  "[PASS] Exact duplicate audit completed\n\n"
)


# ============================================================
# 6. Nearest-neighbour correlation audit
# ============================================================

cat("Computing nearest-neighbour correlations...\n")

log_expr <- log2(
  expr + 1
)

gene_var <- apply(
  log_expr,
  1,
  var
)

valid_gene_idx <- which(
  is.finite(gene_var) &
  gene_var > 0
)

top_n <- min(
  5000L,
  length(valid_gene_idx)
)

top_idx <- valid_gene_idx[
  order(
    gene_var[valid_gene_idx],
    decreasing = TRUE
  )[seq_len(top_n)]
]

cat(
  "Variable genes used for correlation:",
  top_n,
  "\n"
)

cor_mat <- cor(
  log_expr[
    top_idx,
    ,
    drop = FALSE
  ],
  method = "pearson",
  use = "pairwise.complete.obs"
)

rownames(cor_mat) <- sample_ids
colnames(cor_mat) <- sample_ids

diag(cor_mat) <- NA_real_


# ------------------------------------------------------------
# Safe nearest-neighbour extraction
# ------------------------------------------------------------

nearest_index <- vapply(
  seq_len(ncol(cor_mat)),
  function(i) {

    x <- cor_mat[, i]

    valid <- which(
      is.finite(x)
    )

    if (length(valid) == 0) {
      return(NA_integer_)
    }

    valid[
      which.max(
        x[valid]
      )
    ]
  },
  integer(1)
)


nearest_cor <- vapply(
  seq_len(ncol(cor_mat)),
  function(i) {

    idx <- nearest_index[i]

    if (is.na(idx))
      return(NA_real_)

    cor_mat[
      idx,
      i
    ]
  },
  numeric(1)
)


second_nearest_cor <- vapply(
  seq_len(ncol(cor_mat)),
  function(i) {

    x <- cor_mat[, i]

    x <- sort(
      x[
        is.finite(x)
      ],
      decreasing = TRUE
    )

    if (length(x) >= 2) {
      x[2]
    } else {
      NA_real_
    }
  },
  numeric(1)
)


nearest_sample <- rep(
  NA_character_,
  length(nearest_index)
)

valid_nn <- !is.na(
  nearest_index
)

nearest_sample[
  valid_nn
] <- sample_ids[
  nearest_index[
    valid_nn
  ]
]


neighbor_qc <- data.table(
  SampleID = sample_ids,
  NearestSample = nearest_sample,
  NearestCorrelation = nearest_cor,
  SecondNearestCorrelation = second_nearest_cor
)

cat(
  "[PASS] Nearest-neighbour correlation audit completed\n\n"
)


# ============================================================
# 7. Very-high-correlation pair audit
# ============================================================

high_mask <- upper.tri(
  cor_mat
) &
  is.finite(
    cor_mat
  ) &
  cor_mat >= 0.9999

high_idx <- which(
  high_mask,
  arr.ind = TRUE
)

if (nrow(high_idx) > 0) {

  high_corr <- data.table(
    Sample1 =
      rownames(cor_mat)[
        high_idx[, 1]
      ],

    Sample2 =
      colnames(cor_mat)[
        high_idx[, 2]
      ],

    Pearson =
      cor_mat[
        high_idx
      ]
  )

  setorder(
    high_corr,
    -Pearson
  )

} else {

  high_corr <- data.table(
    Sample1 = character(),
    Sample2 = character(),
    Pearson = numeric()
  )
}


# ============================================================
# 8. Within-disease PCA reassessment
# ============================================================

cat("Computing within-disease PCA distances...\n")

pca2 <- copy(
  pca
)

pc_cols <- c(
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5"
)

pca2[
  ,
  WithinDiseasePCDistance :=
    NA_real_
]

pca2[
  ,
  WithinDiseaseThreshold :=
    NA_real_
]

pca2[
  ,
  WithinDiseasePCAFlag :=
    FALSE
]


for (grp in unique(pca2$Disease)) {

  idx <- which(
    pca2$Disease == grp
  )

  if (length(idx) < 3) {
    next
  }

  x <- as.matrix(
    pca2[
      idx,
      ..pc_cols
    ]
  )

  centers <- apply(
    x,
    2,
    median,
    na.rm = TRUE
  )

  scales <- apply(
    x,
    2,
    mad,
    constant = 1.4826,
    na.rm = TRUE
  )

  scales[
    !is.finite(scales) |
    scales == 0
  ] <- 1

  z <- sweep(
    x,
    2,
    centers,
    "-"
  )

  z <- sweep(
    z,
    2,
    scales,
    "/"
  )

  distances <- sqrt(
    rowSums(
      z^2
    )
  )

  dist_median <- median(
    distances,
    na.rm = TRUE
  )

  dist_mad <- mad(
    distances,
    center = dist_median,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (
    !is.finite(dist_mad) ||
    dist_mad == 0
  ) {

    threshold <- Inf

  } else {

    threshold <-
      dist_median +
      5 * dist_mad
  }

  pca2$WithinDiseasePCDistance[
    idx
  ] <- distances

  pca2$WithinDiseaseThreshold[
    idx
  ] <- threshold

  pca2$WithinDiseasePCAFlag[
    idx
  ] <- distances > threshold
}

cat(
  "[PASS] Within-disease PCA reassessment completed\n\n"
)


# ============================================================
# 9. Build forensic master table
#
# Avoid merge() wherever possible so existing column names
# such as PCAQCFlag cannot become PCAQCFlag.x/.y.
# ============================================================

forensic <- copy(
  qc
)

stopifnot(
  identical(
    forensic$SampleID,
    sample_ids
  )
)


# ------------------------------------------------------------
# Add composition metrics
# ------------------------------------------------------------

idx_comp <- match(
  forensic$SampleID,
  composition_qc$SampleID
)

stopifnot(
  !anyNA(idx_comp)
)

forensic[
  ,
  `:=`(
    Top1Fraction =
      composition_qc$Top1Fraction[
        idx_comp
      ],

    Top10Fraction =
      composition_qc$Top10Fraction[
        idx_comp
      ]
  )
]


# ------------------------------------------------------------
# Add neighbour metrics
# ------------------------------------------------------------

idx_nn <- match(
  forensic$SampleID,
  neighbor_qc$SampleID
)

stopifnot(
  !anyNA(idx_nn)
)

forensic[
  ,
  `:=`(
    NearestSample =
      neighbor_qc$NearestSample[
        idx_nn
      ],

    NearestCorrelation =
      neighbor_qc$NearestCorrelation[
        idx_nn
      ],

    SecondNearestCorrelation =
      neighbor_qc$SecondNearestCorrelation[
        idx_nn
      ]
  )
]


# ------------------------------------------------------------
# Add PCA forensic metrics
#
# Keep the original PCAQCFlag from Stage 13.
# ------------------------------------------------------------

idx_pca <- match(
  forensic$SampleID,
  pca2$SampleID
)

stopifnot(
  !anyNA(idx_pca)
)

stopifnot(
  identical(
    as.logical(forensic$PCAQCFlag),
    as.logical(
      pca2$PCAQCFlag[
        idx_pca
      ]
    )
  )
)

forensic[
  ,
  `:=`(
    RobustPCDistance =
      pca2$RobustPCDistance[
        idx_pca
      ],

    WithinDiseasePCDistance =
      pca2$WithinDiseasePCDistance[
        idx_pca
      ],

    WithinDiseaseThreshold =
      pca2$WithinDiseaseThreshold[
        idx_pca
      ],

    WithinDiseasePCAFlag =
      pca2$WithinDiseasePCAFlag[
        idx_pca
      ]
  )
]


# ------------------------------------------------------------
# Add phenotype context
# ------------------------------------------------------------

idx_meta <- match(
  forensic$SampleID,
  meta$SampleID
)

stopifnot(
  !anyNA(idx_meta)
)

forensic[
  ,
  `:=`(
    Histopathology =
      meta$Histopathology[
        idx_meta
      ],

    DeepUlcer =
      meta$DeepUlcer[
        idx_meta
      ],

    ParisStage =
      meta$ParisStage[
        idx_meta
      ]
  )
]


# ============================================================
# 10. Additional nearest-correlation lower-outlier diagnostic
# ============================================================

nn_values <- forensic$NearestCorrelation

nn_med <- median(
  nn_values,
  na.rm = TRUE
)

nn_mad <- mad(
  nn_values,
  center = nn_med,
  constant = 1.4826,
  na.rm = TRUE
)

if (
  is.finite(nn_mad) &&
  nn_mad > 0
) {

  nn_lower_threshold <-
    nn_med -
    5 * nn_mad

} else {

  nn_lower_threshold <- -Inf
}

forensic[
  ,
  LowNearestCorrelationFlag :=
    is.finite(NearestCorrelation) &
    NearestCorrelation <
      nn_lower_threshold
]


# ============================================================
# 11. 254/68 descriptive provenance-segment audit
#
# This is NOT treated as a batch variable.
# ============================================================

risk_number <- suppressWarnings(
  as.integer(
    sub(
      "^CCFA_Risk_",
      "",
      forensic$SampleID
    )
  )
)

stopifnot(
  !anyNA(risk_number),
  all(risk_number >= 1),
  all(risk_number <= 322)
)

forensic[
  ,
  RiskNumber :=
    risk_number
]

forensic[
  ,
  ProvenanceSegment :=
    fifelse(
      RiskNumber <= 254,
      "Risk001_254",
      "Risk255_322"
    )
]


segment_summary <- forensic[
  ,
  .(
    N = .N,

    MedianDetectedGenes =
      median(
        DetectedGenes
      ),

    MedianZeroFraction =
      median(
        ZeroFraction
      ),

    MedianLog2TotalExpression =
      median(
        Log2TotalExpression
      ),

    MedianNearestCorrelation =
      median(
        NearestCorrelation,
        na.rm = TRUE
      ),

    OriginalQCFlagN =
      sum(
        AnyQCFlag
      ),

    OriginalQCFlagFraction =
      mean(
        AnyQCFlag
      ),

    GlobalPCAFlagN =
      sum(
        PCAQCFlag
      ),

    WithinDiseasePCAFlagN =
      sum(
        WithinDiseasePCAFlag
      ),

    LowNearestCorrelationFlagN =
      sum(
        LowNearestCorrelationFlag
      )
  ),
  by = ProvenanceSegment
]


# ============================================================
# 12. Disease-level QC summary
# ============================================================

disease_summary <- forensic[
  ,
  .(
    N = .N,

    OriginalQCFlagN =
      sum(
        AnyQCFlag
      ),

    OriginalQCFlagFraction =
      mean(
        AnyQCFlag
      ),

    UnivariateFlagN =
      sum(
        UnivariateQCFlag
      ),

    GlobalPCAFlagN =
      sum(
        PCAQCFlag
      ),

    WithinDiseasePCAFlagN =
      sum(
        WithinDiseasePCAFlag
      ),

    LowNearestCorrelationFlagN =
      sum(
        LowNearestCorrelationFlag
      ),

    MedianDetectedGenes =
      median(
        DetectedGenes
      ),

    MedianZeroFraction =
      median(
        ZeroFraction
      ),

    MedianNearestCorrelation =
      median(
        NearestCorrelation,
        na.rm = TRUE
      )
  ),
  by = Disease
]


# ============================================================
# 13. Console report
# ============================================================

cat("============================================================\n")
cat("EXACT DUPLICATES\n")
cat("============================================================\n")

cat(
  "Exact duplicate pairs:",
  nrow(exact_dup),
  "\n"
)

if (nrow(exact_dup) > 0) {
  print(exact_dup)
}


cat("\n============================================================\n")
cat("VERY HIGH CORRELATION PAIRS (>=0.9999)\n")
cat("============================================================\n")

cat(
  "Pairs:",
  nrow(high_corr),
  "\n"
)

if (nrow(high_corr) > 0) {

  print(
    head(
      high_corr,
      50
    )
  )
}


cat("\n============================================================\n")
cat("NEAREST-NEIGHBOUR CORRELATION\n")
cat("============================================================\n")

cat(
  "Median nearest correlation:",
  sprintf(
    "%.6f",
    median(
      forensic$NearestCorrelation,
      na.rm = TRUE
    )
  ),
  "\n"
)

cat(
  "Minimum nearest correlation:",
  sprintf(
    "%.6f",
    min(
      forensic$NearestCorrelation,
      na.rm = TRUE
    )
  ),
  "\n"
)

cat(
  "Lower outlier threshold:",
  sprintf(
    "%.6f",
    nn_lower_threshold
  ),
  "\n"
)

cat(
  "Low-nearest-correlation flags:",
  sum(
    forensic$LowNearestCorrelationFlag
  ),
  "\n"
)


cat("\n============================================================\n")
cat("FLAGGED SAMPLE FORENSICS\n")
cat("============================================================\n")

flagged_forensic <- forensic[
  AnyQCFlag == TRUE,
  .(
    SampleID,
    GSM,
    Disease,
    DiseaseSubtype,
    Histopathology,
    DeepUlcer,

    DetectedGenes,
    ZeroFraction,
    Log2TotalExpression,

    Top1Fraction,
    Top10Fraction,

    NearestSample,
    NearestCorrelation,
    SecondNearestCorrelation,

    UnivariateQCFlag,
    PCAQCFlag,
    WithinDiseasePCAFlag,
    LowNearestCorrelationFlag,

    RobustPCDistance,
    WithinDiseasePCDistance,

    RiskNumber,
    ProvenanceSegment
  )
]

print(
  flagged_forensic
)


cat("\n============================================================\n")
cat("FLAG COUNTS\n")
cat("============================================================\n")

cat(
  "Stage 13 AnyQCFlag:",
  sum(
    forensic$AnyQCFlag
  ),
  "\n"
)

cat(
  "Univariate flags:",
  sum(
    forensic$UnivariateQCFlag
  ),
  "\n"
)

cat(
  "Global PCA flags:",
  sum(
    forensic$PCAQCFlag
  ),
  "\n"
)

cat(
  "Within-disease PCA flags:",
  sum(
    forensic$WithinDiseasePCAFlag
  ),
  "\n"
)

cat(
  "Low nearest-correlation flags:",
  sum(
    forensic$LowNearestCorrelationFlag
  ),
  "\n"
)


cat("\n============================================================\n")
cat("FLAG DISTRIBUTION BY DISEASE\n")
cat("============================================================\n")

print(
  disease_summary
)


cat("\n============================================================\n")
cat("254/68 SEGMENT DESCRIPTIVE CHECK\n")
cat("============================================================\n")

print(
  segment_summary
)


# ============================================================
# 14. Save outputs
# ============================================================

fwrite(
  forensic,
  file.path(
    out_dir,
    "GSE57945_flagged_sample_forensic_qc.tsv"
  ),
  sep = "\t"
)

fwrite(
  flagged_forensic,
  file.path(
    out_dir,
    "GSE57945_flagged_samples_detail.tsv"
  ),
  sep = "\t"
)

fwrite(
  exact_dup,
  file.path(
    out_dir,
    "GSE57945_exact_duplicate_sample_pairs.tsv"
  ),
  sep = "\t"
)

fwrite(
  high_corr,
  file.path(
    out_dir,
    "GSE57945_very_high_correlation_sample_pairs.tsv"
  ),
  sep = "\t"
)

fwrite(
  disease_summary,
  file.path(
    out_dir,
    "GSE57945_forensic_qc_by_disease.tsv"
  ),
  sep = "\t"
)

fwrite(
  segment_summary,
  file.path(
    out_dir,
    "GSE57945_254_68_segment_qc_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 15. Final integrity checks
# ============================================================

stopifnot(
  nrow(forensic) == 322,
  identical(
    forensic$SampleID,
    sample_ids
  ),
  !anyDuplicated(
    forensic$SampleID
  ),
  sum(
    forensic$AnyQCFlag
  ) == 13
)

cat("\n============================================================\n")
cat("INTEGRITY CHECKS\n")
cat("============================================================\n")

cat(
  "[PASS] Forensic table contains exactly 322 ordered samples\n"
)

cat(
  "[PASS] Stage 13 original AnyQCFlag count retained at 13\n"
)

cat(
  "[PASS] No expression values were modified\n"
)


cat("\nIMPORTANT:\n")
cat(
  "- QC flags are diagnostic only.\n",
  "- No sample was removed.\n",
  "- No normalization or batch correction was applied.\n",
  "- Correlations use the top 5000 variable genes after log2(x + 1).\n",
  "- Within-disease PCA uses the existing Stage 13 PC1-PC5 coordinates.\n",
  "- The 254/68 split is descriptive only and is NOT treated as a batch.\n",
  "- Global biological extremes must not be interpreted automatically as technical failures.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_FLAGGED_SAMPLE_FORENSIC_QC_COMPLETED\n")
cat("============================================================\n")
