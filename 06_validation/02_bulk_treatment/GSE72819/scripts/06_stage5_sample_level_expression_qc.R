suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE72819 Stage 5
# Sample-level expression QC
#
# Frozen inputs:
#   Stage 3 canonical RPKM:
#       24,656 genes x 73 samples
#
#   Stage 4 canonical sample metadata:
#       73 samples / 73 patients
#
# IMPORTANT:
#   - canonical RPKM is NEVER modified
#   - log2(RPKM + 1) is used ONLY as an in-memory QC copy
#   - diagnostic flag != technical exclusion
#   - no sample is automatically removed
#   - no batch correction is performed
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

expr_file <- file.path(
  root,
  "prepared/03_canonical_expression",
  "GSE72819_canonical_RPKM_24656x73.rds"
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

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Helpers
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


same_chr <- function(x, y) {

  x <- as.character(x)
  y <- as.character(y)

  if (length(x) != length(y)) {
    return(FALSE)
  }

  all(
    !is.na(x) &
      !is.na(y) &
      x == y
  )
}


robust_center_scale <- function(x) {

  x <- as.numeric(x)

  med <- median(
    x,
    na.rm = TRUE
  )

  sc <- mad(
    x,
    center = med,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (!is.finite(sc) || sc <= 0) {

    sc <- sd(
      x,
      na.rm = TRUE
    )
  }

  if (!is.finite(sc) || sc <= 0) {

    sc <- 1
  }

  list(
    center = med,
    scale = sc
  )
}


robust_z <- function(x) {

  pars <- robust_center_scale(x)

  (
    as.numeric(x) -
      pars$center
  ) / pars$scale
}


cat("============================================================\n")
cat("GSE72819 Stage 5: sample-level expression QC\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load frozen canonical expression
# ============================================================

expr <- readRDS(
  expr_file
)

check(
  is.matrix(expr),
  "Canonical expression object is a matrix"
)

check(
  nrow(expr) == 24656,
  "Canonical expression contains 24,656 genes"
)

check(
  ncol(expr) == 73,
  "Canonical expression contains 73 samples"
)

check(
  !is.null(rownames(expr)),
  "Canonical expression has gene symbols"
)

check(
  !is.null(colnames(expr)),
  "Canonical expression has sample identifiers"
)

check(
  uniqueN(rownames(expr)) == 24656,
  "Canonical gene symbols are unique"
)

check(
  uniqueN(colnames(expr)) == 73,
  "Canonical sample identifiers are unique"
)

check(
  sum(is.na(expr)) == 0,
  "Canonical expression contains no NA"
)

check(
  sum(!is.finite(expr)) == 0,
  "Canonical expression contains no non-finite values"
)

check(
  sum(expr < 0) == 0,
  "Canonical expression contains no negative values"
)

sample_ids <- as.character(
  colnames(expr)
)


# ============================================================
# 2. Load frozen canonical metadata
# ============================================================

meta <- readRDS(
  meta_file
)

setDT(meta)

required_meta <- c(
  "GSM",
  "PatientID",
  "Disease",
  "Tissue",
  "ResponseWeek10",
  "PriorAntiTNF",
  "Run",
  "Lane",
  "RecoveryEligible",
  "ResponseAnalysisEligible"
)

check(
  all(required_meta %in% names(meta)),
  "Canonical metadata contains all required fields"
)

check(
  nrow(meta) == 73,
  "Canonical metadata contains 73 samples"
)

check(
  uniqueN(meta$GSM) == 73,
  "Canonical metadata contains 73 unique GSM IDs"
)

check(
  same_chr(
    meta$GSM,
    sample_ids
  ),
  "Metadata order exactly matches expression columns"
)

check(
  all(meta$RecoveryEligible),
  "All 73 samples enter Stage-5 QC"
)


# ============================================================
# 3. Freeze source numeric payload
# ============================================================

raw_payload_before <- as.vector(
  expr
)

check(
  length(raw_payload_before) ==
    24656 * 73,
  "Expected canonical RPKM payload size confirmed"
)


# ============================================================
# 4. Basic raw-RPKM per-sample metrics
# ============================================================

rpkm_sum <- colSums(
  expr
)

detected_genes <- colSums(
  expr > 0
)

zero_genes <- colSums(
  expr == 0
)

zero_fraction <- colMeans(
  expr == 0
)

median_rpkm <- apply(
  expr,
  2,
  median
)

mean_rpkm <- colMeans(
  expr
)

max_rpkm <- apply(
  expr,
  2,
  max
)


check(
  all(is.finite(rpkm_sum)),
  "Per-sample RPKM sums are finite"
)

check(
  all(detected_genes > 0),
  "Every sample contains detected genes"
)


# ============================================================
# 5. QC-only log2(RPKM + 1) copy
#
# This object is NOT saved as the canonical expression matrix.
# ============================================================

qc_log <- log2(
  expr + 1
)

check(
  is.matrix(qc_log),
  "QC-only log2(RPKM+1) matrix created"
)

check(
  all(is.finite(qc_log)),
  "QC-only log matrix contains only finite values"
)

check(
  nrow(qc_log) == 24656 &&
    ncol(qc_log) == 73,
  "QC-only log matrix dimensions unchanged"
)


median_log_rpkm <- apply(
  qc_log,
  2,
  median
)

mean_log_rpkm <- colMeans(
  qc_log
)

sd_log_rpkm <- apply(
  qc_log,
  2,
  sd
)

check(
  all(
    is.finite(sd_log_rpkm) &
      sd_log_rpkm > 0
  ),
  "Every sample has non-zero log-expression variance"
)


# ============================================================
# 6. Correlation to cohort median profile
# ============================================================

median_profile <- apply(
  qc_log,
  1,
  median
)

check(
  length(median_profile) == 24656,
  "Median expression profile contains 24,656 genes"
)

cor_to_median <- as.numeric(
  cor(
    qc_log,
    median_profile,
    method = "pearson"
  )
)

check(
  length(cor_to_median) == 73,
  "Correlation-to-median calculated for all 73 samples"
)

check(
  all(is.finite(cor_to_median)),
  "All correlation-to-median values are finite"
)


# ============================================================
# 7. Full pairwise sample correlation
# ============================================================

pair_cor <- cor(
  qc_log,
  method = "pearson"
)

check(
  is.matrix(pair_cor),
  "Pairwise sample correlation matrix created"
)

check(
  nrow(pair_cor) == 73 &&
    ncol(pair_cor) == 73,
  "Pairwise correlation matrix is 73 x 73"
)

check(
  all(is.finite(pair_cor)),
  "All pairwise sample correlations are finite"
)


# ============================================================
# 8. Nearest-neighbour correlations
# ============================================================

nn_matrix <- pair_cor

diag(nn_matrix) <- -Inf

nearest_index <- max.col(
  nn_matrix,
  ties.method = "first"
)

nearest_corr <- nn_matrix[
  cbind(
    seq_len(nrow(nn_matrix)),
    nearest_index
  )
]

nearest_gsm <- sample_ids[
  nearest_index
]

check(
  length(nearest_corr) == 73,
  "Nearest-neighbour correlation calculated for all samples"
)

check(
  all(is.finite(nearest_corr)),
  "All nearest-neighbour correlations are finite"
)


# ============================================================
# 9. Exact duplicate audit
#
# Exact equality is tested on the frozen raw canonical RPKM.
# ============================================================

pair_index <- t(
  combn(
    seq_len(73),
    2
  )
)

exact_duplicate <- vapply(
  seq_len(nrow(pair_index)),
  function(k) {

    i <- pair_index[k, 1]
    j <- pair_index[k, 2]

    isTRUE(
      all(
        expr[, i] ==
          expr[, j]
      )
    )
  },
  logical(1)
)


if (any(exact_duplicate)) {

  exact_idx <- pair_index[
    exact_duplicate,
    ,
    drop = FALSE
  ]

  exact_pairs <- data.table(
    Sample1 =
      sample_ids[
        exact_idx[, 1]
      ],

    Sample2 =
      sample_ids[
        exact_idx[, 2]
      ],

    Pearson_log2_RPKM =
      vapply(
        seq_len(nrow(exact_idx)),
        function(k) {

          pair_cor[
            exact_idx[k, 1],
            exact_idx[k, 2]
          ]
        },
        numeric(1)
      ),

    ExactDuplicate =
      TRUE
  )

} else {

  exact_pairs <- data.table(
    Sample1 = character(),
    Sample2 = character(),
    Pearson_log2_RPKM = numeric(),
    ExactDuplicate = logical()
  )
}


cat("\n============================================================\n")
cat("EXACT DUPLICATE AUDIT\n")
cat("============================================================\n")

cat(
  "Exact duplicate pairs:",
  nrow(exact_pairs),
  "\n"
)


# ============================================================
# 10. Near-duplicate audit
#
# Threshold follows the same diagnostic criterion used in the
# previous validation audits:
# Pearson(log2(RPKM+1)) >= 0.9995
# ============================================================

near_duplicate_threshold <- 0.9995

near_idx <- which(
  upper.tri(pair_cor) &
    pair_cor >= near_duplicate_threshold,
  arr.ind = TRUE
)


if (nrow(near_idx) > 0) {

  near_all <- data.table(
    Sample1 =
      sample_ids[
        near_idx[, 1]
      ],

    Sample2 =
      sample_ids[
        near_idx[, 2]
      ],

    Pearson_log2_RPKM =
      pair_cor[
        near_idx
      ]
  )

  exact_keys <- if (nrow(exact_pairs) > 0) {

    paste(
      exact_pairs$Sample1,
      exact_pairs$Sample2,
      sep = "||"
    )

  } else {

    character()
  }

  near_all[
    ,
    PairKey :=
      paste(
        Sample1,
        Sample2,
        sep = "||"
      )
  ]

  near_all[
    ,
    ExactDuplicate :=
      PairKey %in%
        exact_keys
  ]

  near_pairs <- near_all[
    ExactDuplicate == FALSE,
    .(
      Sample1,
      Sample2,
      Pearson_log2_RPKM
    )
  ]

} else {

  near_pairs <- data.table(
    Sample1 = character(),
    Sample2 = character(),
    Pearson_log2_RPKM = numeric()
  )
}


cat("\n============================================================\n")
cat("NEAR-DUPLICATE AUDIT\n")
cat("============================================================\n")

cat(
  "Near-duplicate threshold:",
  near_duplicate_threshold,
  "\n"
)

cat(
  "Near-duplicate non-exact pairs:",
  nrow(near_pairs),
  "\n"
)


# ============================================================
# 11. PCA diagnostic
#
# Top 3,000 variable genes from QC-only log matrix.
# PCA is diagnostic only and cannot trigger automatic exclusion.
# ============================================================

gene_variance <- apply(
  qc_log,
  1,
  var
)

variable_gene_idx <- which(
  is.finite(gene_variance) &
    gene_variance > 0
)

check(
  length(variable_gene_idx) >= 100,
  "Sufficient variable genes available for PCA"
)

n_pca_genes <- min(
  3000L,
  length(variable_gene_idx)
)

variable_gene_idx <- variable_gene_idx[
  order(
    gene_variance[
      variable_gene_idx
    ],
    decreasing = TRUE
  )
]

top_gene_idx <- variable_gene_idx[
  seq_len(n_pca_genes)
]

pca_input <- t(
  qc_log[
    top_gene_idx,
    ,
    drop = FALSE
  ]
)

check(
  nrow(pca_input) == 73,
  "PCA input contains 73 samples"
)

check(
  ncol(pca_input) == n_pca_genes,
  "PCA input contains expected number of variable genes"
)


pca_rank <- min(
  10L,
  nrow(pca_input) - 1L,
  ncol(pca_input)
)

pca <- prcomp(
  pca_input,
  center = TRUE,
  scale. = TRUE,
  rank. = pca_rank
)

pc_scores <- pca$x[
  ,
  seq_len(pca_rank),
  drop = FALSE
]

check(
  nrow(pc_scores) == 73,
  "PCA scores generated for all 73 samples"
)

check(
  ncol(pc_scores) == pca_rank,
  "Expected PCA components generated"
)

check(
  all(is.finite(pc_scores)),
  "All PCA scores are finite"
)


pca_sdev <- pca$sdev[
  seq_len(pca_rank)
]

check(
  all(
    is.finite(pca_sdev) &
      pca_sdev > 0
  ),
  "Retained PCA components have positive variance"
)


pc_standardized <- sweep(
  pc_scores,
  2,
  pca_sdev,
  FUN = "/"
)

pca_distance <- sqrt(
  rowSums(
    pc_standardized ^ 2
  )
)

check(
  all(is.finite(pca_distance)),
  "PCA diagnostic distances are finite"
)


# ============================================================
# 12. Robust diagnostic scores
#
# cutoff = |robust Z| > 5
# Correlation metrics are flagged only on LOW side.
# ============================================================

robust_z_cutoff <- 5

z_detected <- robust_z(
  detected_genes
)

z_zero_fraction <- robust_z(
  zero_fraction
)

z_log10_sum <- robust_z(
  log10(
    rpkm_sum + 1
  )
)

z_median_log <- robust_z(
  median_log_rpkm
)

z_cor_median <- robust_z(
  cor_to_median
)

z_nearest <- robust_z(
  nearest_corr
)

z_pca_distance <- robust_z(
  pca_distance
)


distribution_flag <- (
  abs(z_detected) > robust_z_cutoff |
    abs(z_zero_fraction) > robust_z_cutoff |
    abs(z_log10_sum) > robust_z_cutoff |
    abs(z_median_log) > robust_z_cutoff
)

correlation_flag <- (
  z_cor_median < -robust_z_cutoff |
    z_nearest < -robust_z_cutoff
)

pca_flag <- (
  z_pca_distance > robust_z_cutoff
)


# ============================================================
# 13. Duplicate flags per sample
# ============================================================

exact_duplicate_samples <- unique(
  c(
    exact_pairs$Sample1,
    exact_pairs$Sample2
  )
)

near_duplicate_samples <- unique(
  c(
    near_pairs$Sample1,
    near_pairs$Sample2
  )
)

exact_duplicate_flag <- (
  sample_ids %in%
    exact_duplicate_samples
)

near_duplicate_flag <- (
  sample_ids %in%
    near_duplicate_samples
)


# ============================================================
# 14. Build sample-level QC table
# ============================================================

qc <- data.table(
  GSM = sample_ids,

  PatientID =
    as.character(
      meta$PatientID
    ),

  Run =
    as.character(
      meta$Run
    ),

  Lane =
    as.character(
      meta$Lane
    ),

  ResponseWeek10 =
    as.character(
      meta$ResponseWeek10
    ),

  PriorAntiTNF =
    as.character(
      meta$PriorAntiTNF
    ),

  RPKMSum =
    rpkm_sum,

  Log10RPKMSum =
    log10(
      rpkm_sum + 1
    ),

  DetectedGenes =
    detected_genes,

  ZeroGenes =
    zero_genes,

  ZeroFraction =
    zero_fraction,

  MedianRPKM =
    median_rpkm,

  MeanRPKM =
    mean_rpkm,

  MaxRPKM =
    max_rpkm,

  MedianLog2RPKM =
    median_log_rpkm,

  MeanLog2RPKM =
    mean_log_rpkm,

  SDLog2RPKM =
    sd_log_rpkm,

  CorrelationToMedian =
    cor_to_median,

  NearestNeighborGSM =
    nearest_gsm,

  NearestNeighborCorrelation =
    nearest_corr,

  PCADistance =
    pca_distance,

  Z_DetectedGenes =
    z_detected,

  Z_ZeroFraction =
    z_zero_fraction,

  Z_Log10RPKMSum =
    z_log10_sum,

  Z_MedianLog2RPKM =
    z_median_log,

  Z_CorrelationToMedian =
    z_cor_median,

  Z_NearestNeighborCorrelation =
    z_nearest,

  Z_PCADistance =
    z_pca_distance,

  DistributionFlag =
    distribution_flag,

  CorrelationFlag =
    correlation_flag,

  PCAFlag =
    pca_flag,

  ExactDuplicateFlag =
    exact_duplicate_flag,

  NearDuplicateFlag =
    near_duplicate_flag
)


qc[
  ,
  DiagnosticFlag :=
    DistributionFlag |
    CorrelationFlag |
    PCAFlag |
    ExactDuplicateFlag |
    NearDuplicateFlag
]


# ============================================================
# 15. Human-readable flag reasons
# ============================================================

qc[
  ,
  FlagReasons :=
    vapply(
      seq_len(.N),
      function(i) {

        reasons <- character()

        if (DistributionFlag[i]) {
          reasons <- c(
            reasons,
            "expression_distribution"
          )
        }

        if (CorrelationFlag[i]) {
          reasons <- c(
            reasons,
            "low_correlation"
          )
        }

        if (PCAFlag[i]) {
          reasons <- c(
            reasons,
            "PCA_outlier"
          )
        }

        if (ExactDuplicateFlag[i]) {
          reasons <- c(
            reasons,
            "exact_duplicate"
          )
        }

        if (NearDuplicateFlag[i]) {
          reasons <- c(
            reasons,
            "near_duplicate"
          )
        }

        if (length(reasons) == 0) {
          return("")
        }

        paste(
          reasons,
          collapse = ";"
        )
      },
      character(1)
    )
]


# ============================================================
# 16. NO automatic exclusions
# ============================================================

qc[
  ,
  TechnicalExclusion :=
    FALSE
]

qc[
  ,
  Retained :=
    TRUE
]


check(
  sum(qc$TechnicalExclusion) == 0,
  "No sample automatically assigned technical exclusion"
)

check(
  sum(qc$Retained) == 73,
  "All 73 samples remain retained pending review"
)


# ============================================================
# 17. PCA score table
# ============================================================

pca_scores_dt <- data.table(
  GSM = sample_ids,
  PatientID = as.character(meta$PatientID),
  Run = as.character(meta$Run),
  ResponseWeek10 = as.character(meta$ResponseWeek10),
  PriorAntiTNF = as.character(meta$PriorAntiTNF)
)

pca_scores_dt <- cbind(
  pca_scores_dt,
  as.data.table(
    pc_scores
  )
)

setnames(
  pca_scores_dt,
  c(
    "GSM",
    "PatientID",
    "Run",
    "ResponseWeek10",
    "PriorAntiTNF",
    paste0(
      "PC",
      seq_len(pca_rank)
    )
  )
)


# ============================================================
# 18. PCA association with Run
#
# Diagnostic only.
# No adjustment or batch correction is performed.
# ============================================================

run_levels <- sort(
  unique(
    as.character(
      meta$Run
    )
  )
)

check(
  length(run_levels) == 2,
  "Exactly two sequencing runs present"
)


pca_run_assoc <- rbindlist(
  lapply(
    seq_len(pca_rank),
    function(k) {

      score <- pc_scores[, k]

      run <- as.character(
        meta$Run
      )

      p_value <- tryCatch(
        wilcox.test(
          score ~ run,
          exact = FALSE
        )$p.value,
        error = function(e) {
          NA_real_
        }
      )

      data.table(
        PC =
          paste0(
            "PC",
            k
          ),

        Median_R1510 =
          median(
            score[
              run == "R1510"
            ]
          ),

        Median_R1518 =
          median(
            score[
              run == "R1518"
            ]
          ),

        MedianDifference_R1518_minus_R1510 =
          median(
            score[
              run == "R1518"
            ]
          ) -
          median(
            score[
              run == "R1510"
            ]
          ),

        WilcoxonP =
          p_value
      )
    }
  )
)


pca_run_assoc[
  ,
  FDR_BH :=
    p.adjust(
      WilcoxonP,
      method = "BH"
    )
]


# ============================================================
# 19. Descriptive QC by Run
# ============================================================

qc_by_run <- qc[
  ,
  .(
    N = .N,

    Median_RPKMSum =
      median(RPKMSum),

    Median_DetectedGenes =
      median(DetectedGenes),

    Median_ZeroFraction =
      median(ZeroFraction),

    Median_CorrelationToMedian =
      median(CorrelationToMedian),

    Median_NearestNeighborCorrelation =
      median(NearestNeighborCorrelation),

    Median_PCADistance =
      median(PCADistance),

    DiagnosticFlags =
      sum(DiagnosticFlag),

    ExactDuplicateFlags =
      sum(ExactDuplicateFlag),

    NearDuplicateFlags =
      sum(NearDuplicateFlag)
  ),
  by = Run
][order(Run)]


# ============================================================
# 20. Flag summaries by clinical/technical strata
# ============================================================

flag_by_run <- qc[
  ,
  .(
    N = .N,
    Flagged = sum(DiagnosticFlag)
  ),
  by = Run
][order(Run)]


flag_by_response <- qc[
  ,
  .(
    N = .N,
    Flagged = sum(DiagnosticFlag)
  ),
  by = ResponseWeek10
][order(ResponseWeek10)]


flag_by_prior_tnf <- qc[
  ,
  .(
    N = .N,
    Flagged = sum(DiagnosticFlag)
  ),
  by = PriorAntiTNF
][order(PriorAntiTNF)]


# ============================================================
# 21. Threshold audit
# ============================================================

cor_pars <- robust_center_scale(
  cor_to_median
)

nn_pars <- robust_center_scale(
  nearest_corr
)

pca_pars <- robust_center_scale(
  pca_distance
)


thresholds <- data.table(
  Metric = c(
    "Robust_Z_cutoff",
    "Near_duplicate_Pearson_cutoff",
    "PCA_variable_genes",
    "PCA_components",
    "CorrelationToMedian_center",
    "CorrelationToMedian_scale",
    "CorrelationToMedian_lower_flag_threshold",
    "NearestNeighbor_center",
    "NearestNeighbor_scale",
    "NearestNeighbor_lower_flag_threshold",
    "PCADistance_center",
    "PCADistance_scale",
    "PCADistance_upper_flag_threshold"
  ),

  Value = c(
    robust_z_cutoff,
    near_duplicate_threshold,
    n_pca_genes,
    pca_rank,

    cor_pars$center,
    cor_pars$scale,
    cor_pars$center -
      robust_z_cutoff *
      cor_pars$scale,

    nn_pars$center,
    nn_pars$scale,
    nn_pars$center -
      robust_z_cutoff *
      nn_pars$scale,

    pca_pars$center,
    pca_pars$scale,
    pca_pars$center +
      robust_z_cutoff *
      pca_pars$scale
  )
)


# ============================================================
# 22. Console summaries
# ============================================================

cat("\n============================================================\n")
cat("SAMPLE-LEVEL QC SUMMARY\n")
cat("============================================================\n")

qc_summary <- data.table(
  Metric = c(
    "Samples_evaluated",
    "Samples_retained",
    "Technical_exclusions",
    "Diagnostic_flags",
    "Distribution_flags",
    "Correlation_flags",
    "PCA_flags",
    "Exact_duplicate_flagged_samples",
    "Near_duplicate_flagged_samples",
    "Exact_duplicate_pairs",
    "Near_duplicate_pairs"
  ),

  N = c(
    nrow(qc),
    sum(qc$Retained),
    sum(qc$TechnicalExclusion),
    sum(qc$DiagnosticFlag),
    sum(qc$DistributionFlag),
    sum(qc$CorrelationFlag),
    sum(qc$PCAFlag),
    sum(qc$ExactDuplicateFlag),
    sum(qc$NearDuplicateFlag),
    nrow(exact_pairs),
    nrow(near_pairs)
  )
)

print(qc_summary)


cat("\n============================================================\n")
cat("CORE QC DISTRIBUTIONS\n")
cat("============================================================\n")

core_summary <- data.table(
  Metric = c(
    "DetectedGenes_median",
    "DetectedGenes_min",
    "DetectedGenes_max",
    "ZeroFraction_median",
    "ZeroFraction_min",
    "ZeroFraction_max",
    "CorrelationToMedian_median",
    "CorrelationToMedian_min",
    "CorrelationToMedian_max",
    "NearestNeighborCorrelation_median",
    "NearestNeighborCorrelation_min",
    "NearestNeighborCorrelation_max",
    "PCADistance_median",
    "PCADistance_min",
    "PCADistance_max"
  ),

  Value = c(
    median(qc$DetectedGenes),
    min(qc$DetectedGenes),
    max(qc$DetectedGenes),

    median(qc$ZeroFraction),
    min(qc$ZeroFraction),
    max(qc$ZeroFraction),

    median(qc$CorrelationToMedian),
    min(qc$CorrelationToMedian),
    max(qc$CorrelationToMedian),

    median(qc$NearestNeighborCorrelation),
    min(qc$NearestNeighborCorrelation),
    max(qc$NearestNeighborCorrelation),

    median(qc$PCADistance),
    min(qc$PCADistance),
    max(qc$PCADistance)
  )
)

print(core_summary)


cat("\n============================================================\n")
cat("QC BY RUN\n")
cat("============================================================\n")

print(qc_by_run)


cat("\n============================================================\n")
cat("PCA ASSOCIATION WITH RUN\n")
cat("============================================================\n")

print(pca_run_assoc)


cat("\n============================================================\n")
cat("FLAGGED SAMPLES\n")
cat("============================================================\n")

if (any(qc$DiagnosticFlag)) {

  print(
    qc[
      DiagnosticFlag == TRUE,
      .(
        GSM,
        PatientID,
        Run,
        ResponseWeek10,
        PriorAntiTNF,
        DetectedGenes,
        ZeroFraction,
        CorrelationToMedian,
        NearestNeighborGSM,
        NearestNeighborCorrelation,
        PCADistance,
        FlagReasons
      )
    ]
  )

} else {

  cat("None\n")
}


# ============================================================
# 23. Verify canonical matrix was never modified
# ============================================================

raw_payload_after <- as.vector(
  expr
)

check(
  identical(
    raw_payload_before,
    raw_payload_after
  ),
  "Frozen canonical RPKM remained numerically unchanged throughout QC"
)


# ============================================================
# 24. Save outputs
# ============================================================

fwrite(
  qc,
  file.path(
    out_dir,
    "GSE72819_sample_level_qc.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  qc[
    DiagnosticFlag == TRUE
  ],
  file.path(
    out_dir,
    "GSE72819_sample_qc_flags.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  exact_pairs,
  file.path(
    out_dir,
    "GSE72819_exact_duplicate_pairs.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  near_pairs,
  file.path(
    out_dir,
    "GSE72819_near_duplicate_pairs.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  pca_scores_dt,
  file.path(
    out_dir,
    "GSE72819_qc_pca_scores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  pca_run_assoc,
  file.path(
    out_dir,
    "GSE72819_qc_pca_run_association.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  qc_by_run,
  file.path(
    out_dir,
    "GSE72819_qc_by_run.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  flag_by_run,
  file.path(
    out_dir,
    "GSE72819_flags_by_run.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  flag_by_response,
  file.path(
    out_dir,
    "GSE72819_flags_by_response.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  flag_by_prior_tnf,
  file.path(
    out_dir,
    "GSE72819_flags_by_prior_antiTNF.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  thresholds,
  file.path(
    out_dir,
    "GSE72819_sample_qc_thresholds.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  qc_summary,
  file.path(
    out_dir,
    "GSE72819_stage5_qc_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  core_summary,
  file.path(
    out_dir,
    "GSE72819_stage5_core_qc_distributions.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 25. Artifact hashes
# ============================================================

artifact_md5 <- tools::md5sum(
  c(
    expr_file,
    meta_file
  )
)

fwrite(
  data.table(
    File = names(artifact_md5),
    MD5 = unname(artifact_md5)
  ),
  file.path(
    out_dir,
    "GSE72819_stage5_input_md5.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 26. Final structural gate
#
# Flags do NOT cause failure.
# ============================================================

critical_pass <- all(c(
  nrow(qc) == 73,
  uniqueN(qc$GSM) == 73,
  sum(qc$TechnicalExclusion) == 0,
  sum(qc$Retained) == 73,
  all(is.finite(qc$CorrelationToMedian)),
  all(is.finite(qc$NearestNeighborCorrelation)),
  all(is.finite(qc$PCADistance)),
  identical(
    raw_payload_before,
    raw_payload_after
  )
))


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] 73 samples evaluated\n")
  cat("[PASS] 73 samples retained pending interpretation\n")
  cat("[PASS] Technical exclusions = 0\n")
  cat(
    "[INFO] Diagnostic flagged samples = ",
    sum(qc$DiagnosticFlag),
    "\n",
    sep = ""
  )
  cat(
    "[INFO] Exact duplicate pairs = ",
    nrow(exact_pairs),
    "\n",
    sep = ""
  )
  cat(
    "[INFO] Near-duplicate pairs = ",
    nrow(near_pairs),
    "\n",
    sep = ""
  )
  cat("[PASS] Canonical RPKM was never modified\n")
  cat("[PASS] No normalization performed on canonical matrix\n")
  cat("[PASS] No batch correction performed\n")
  cat("[PASS] No automatic sample exclusion performed\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE5_SAMPLE_LEVEL_QC_COMPLETED\n")

} else {

  cat("[REVIEW] Structural Stage-5 QC checks failed\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE5_SAMPLE_LEVEL_QC\n")
}

cat("============================================================\n")
