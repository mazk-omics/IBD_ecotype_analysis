suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 2E v3
# Multi-transcript-cluster expression concordance audit
#
# Frozen inputs:
#   Stage 1 authoritative author-provided log2 RMA matrix
#   Stage 2D FINAL probe -> current-gene decision table
#
# Goals:
#   - audit 618 genes represented by >1 retained transcript cluster
#   - quantify pairwise Pearson / Spearman concordance
#   - quantify baseline offsets and sample-wise expression spread
#   - record diagnostic representative-cluster candidates
#
# IMPORTANT:
#   NO expression normalization
#   NO transcript-cluster selection
#   NO transcript-cluster collapse
#   NO gene-level matrix generated
# ============================================================


# ============================================================
# 0. Paths
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

expr_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE73661_authoritative_log2_RMA_33252x178.rds"
)

decision_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_canonical_probe_gene_decision_table_FINAL.tsv"
)

out_dir <- file.path(
  root,
  "prepared/02_gene_harmonization"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. Helpers
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


safe_cor <- function(x, y, method) {

  sx <- sd(x)
  sy <- sd(y)

  if (
    !is.finite(sx) ||
    !is.finite(sy) ||
    sx == 0 ||
    sy == 0
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x,
      y,
      method = method
    )
  )
}


safe_quantile <- function(x, p) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    quantile(
      x,
      probs = p,
      names = FALSE,
      type = 7
    )
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 2E v3: multi-cluster expression audit\n")
cat("============================================================\n\n")

t0 <- proc.time()[["elapsed"]]


# ============================================================
# 2. Load frozen inputs
# ============================================================

check(
  file.exists(expr_file),
  "Frozen Stage-1 authoritative expression exists"
)

check(
  file.exists(decision_file),
  "Frozen Stage-2D FINAL decision table exists"
)


expr <- readRDS(
  expr_file
)

decision <- fread(
  decision_file,
  showProgress = FALSE
)


check(
  identical(
    dim(expr),
    c(33252L, 178L)
  ),
  "Authoritative expression = 33,252 x 178"
)

check(
  nrow(decision) == 33252L,
  "Stage-2D FINAL decision table contains 33,252 rows"
)


# ============================================================
# 3. Freeze identifier types
# ============================================================

decision[
  ,
  SourceProbeID := as.character(SourceProbeID)
]

decision[
  ,
  FinalGeneID := as.character(FinalGeneID)
]

decision[
  ,
  FinalSymbol := as.character(FinalSymbol)
]


check(
  is.character(decision$SourceProbeID),
  "SourceProbeID frozen as character"
)

check(
  is.character(decision$FinalGeneID),
  "FinalGeneID frozen as character"
)

check(
  is.character(decision$FinalSymbol),
  "FinalSymbol frozen as character"
)


check(
  identical(
    decision$SourceProbeID,
    as.character(rownames(expr))
  ),
  "Stage-2D source order matches expression exactly"
)


# ============================================================
# 4. IncludeFinal
# ============================================================

if (!is.logical(decision$IncludeFinal)) {

  decision[
    ,
    IncludeFinal := as.logical(IncludeFinal)
  ]
}


check(
  is.logical(decision$IncludeFinal),
  "IncludeFinal is logical"
)

check(
  !anyNA(decision$IncludeFinal),
  "IncludeFinal contains no NA"
)


# ============================================================
# 5. Retained mappings
# ============================================================

retained <- decision[
  IncludeFinal == TRUE
]


check(
  nrow(retained) == 19229L,
  "19,229 retained transcript clusters recovered"
)

check(
  uniqueN(retained$FinalGeneID) == 18422L,
  "18,422 unique final genes recovered"
)

check(
  !anyNA(retained$FinalGeneID),
  "Retained FinalGeneID contains no NA"
)

check(
  !anyNA(retained$FinalSymbol),
  "Retained FinalSymbol contains no NA"
)


# ============================================================
# 6. Align expression
# ============================================================

expr_probe_ids <- as.character(
  rownames(expr)
)

expr_idx <- match(
  retained$SourceProbeID,
  expr_probe_ids
)


check(
  !anyNA(expr_idx),
  "Every retained transcript cluster found in expression matrix"
)


retained_expr <- expr[
  expr_idx,
  ,
  drop = FALSE
]

rownames(retained_expr) <- retained$SourceProbeID


check(
  identical(
    dim(retained_expr),
    c(19229L, 178L)
  ),
  "Retained expression matrix = 19,229 x 178"
)

check(
  all(is.finite(retained_expr)),
  "Retained expression contains only finite values"
)


# ============================================================
# 7. Gene multiplicity
# ============================================================

gene_counts <- retained[
  ,
  .(
    N_SourceTranscriptClusters = .N
  ),
  by = .(
    FinalGeneID,
    FinalSymbol
  )
]


single_genes <- gene_counts[
  N_SourceTranscriptClusters == 1L
]

multi_genes <- gene_counts[
  N_SourceTranscriptClusters > 1L
]


check(
  nrow(single_genes) == 17804L,
  "17,804 single-cluster genes confirmed"
)

check(
  nrow(multi_genes) == 618L,
  "618 multi-cluster genes confirmed"
)


multi_gene_ids <- multi_genes$FinalGeneID

multi_retained <- retained[
  FinalGeneID %in% multi_gene_ids
]


check(
  nrow(multi_retained) == 1425L,
  "1,425 transcript clusters belong to multi-cluster genes"
)


# ============================================================
# 8. Per-transcript-cluster statistics
# ============================================================

cat(
  "Calculating per-transcript-cluster expression statistics...\n"
)


cluster_stats <- data.table(
  SourceProbeID =
    retained$SourceProbeID,

  FinalGeneID =
    retained$FinalGeneID,

  FinalSymbol =
    retained$FinalSymbol,

  MeanExpression =
    rowMeans(retained_expr),

  MedianExpression =
    apply(
      retained_expr,
      1L,
      median
    ),

  SDExpression =
    apply(
      retained_expr,
      1L,
      sd
    ),

  IQRExpression =
    apply(
      retained_expr,
      1L,
      IQR
    ),

  MADExpression =
    apply(
      retained_expr,
      1L,
      mad
    ),

  MinExpression =
    apply(
      retained_expr,
      1L,
      min
    ),

  MaxExpression =
    apply(
      retained_expr,
      1L,
      max
    )
)


# Avoid data.table join entirely.
multiplicity_lookup <- setNames(
  gene_counts$N_SourceTranscriptClusters,
  gene_counts$FinalGeneID
)


cluster_stats[
  ,
  N_SourceTranscriptClusters :=
    as.integer(
      multiplicity_lookup[
        FinalGeneID
      ]
    )
]


check(
  nrow(cluster_stats) == 19229L,
  "Per-cluster statistics generated for all retained features"
)

check(
  !anyNA(cluster_stats$N_SourceTranscriptClusters),
  "Every retained cluster assigned gene multiplicity"
)


# ============================================================
# 9. Construct multi-cluster groups
# ============================================================

multi_row_idx <- which(
  retained$FinalGeneID %in%
    multi_gene_ids
)


multi_indices <- split(
  multi_row_idx,
  retained$FinalGeneID[
    multi_row_idx
  ]
)


check(
  length(multi_indices) == 618L,
  "618 multi-cluster expression groups constructed"
)


# ============================================================
# 10. Prepare output containers
# ============================================================

pairwise_list <- vector(
  "list",
  length(multi_indices)
)

gene_summary_list <- vector(
  "list",
  length(multi_indices)
)

candidate_list <- vector(
  "list",
  length(multi_indices)
)


cat(
  "Auditing pairwise concordance across 618 genes...\n"
)


# ============================================================
# 11. Audit each multi-cluster gene
# ============================================================

counter <- 0L


for (gene_id in names(multi_indices)) {

  counter <- counter + 1L


  # IMPORTANT:
  # list extraction must use [[gene_id]]
  idx <- multi_indices[[gene_id]]


  if (
    is.null(idx) ||
    length(idx) < 2L
  ) {
    stop(
      paste0(
        "[FAIL] Invalid multi-cluster index for GeneID ",
        gene_id
      ),
      call. = FALSE
    )
  }


  x <- retained_expr[
    idx,
    ,
    drop = FALSE
  ]


  probe_ids <- retained$SourceProbeID[
    idx
  ]


  gene_ids_check <- unique(
    retained$FinalGeneID[
      idx
    ]
  )


  if (
    length(gene_ids_check) != 1L ||
    gene_ids_check != gene_id
  ) {
    stop(
      paste0(
        "[FAIL] Gene grouping mismatch for GeneID ",
        gene_id
      ),
      call. = FALSE
    )
  }


  symbol <- retained$FinalSymbol[
    idx[1L]
  ]


  k <- nrow(x)


  if (k < 2L) {
    stop(
      paste0(
        "[FAIL] Expected >=2 transcript clusters for GeneID ",
        gene_id
      ),
      call. = FALSE
    )
  }


  # ----------------------------------------------------------
  # Cluster statistics in exact probe order
  # ----------------------------------------------------------

  cs_idx <- match(
    probe_ids,
    cluster_stats$SourceProbeID
  )


  if (anyNA(cs_idx)) {
    stop(
      paste0(
        "[FAIL] Missing cluster statistics for GeneID ",
        gene_id
      ),
      call. = FALSE
    )
  }


  gene_cluster_stats <- cluster_stats[
    cs_idx
  ]


  if (
    !identical(
      gene_cluster_stats$SourceProbeID,
      probe_ids
    )
  ) {
    stop(
      paste0(
        "[FAIL] Cluster-stat order mismatch for GeneID ",
        gene_id
      ),
      call. = FALSE
    )
  }


  # ----------------------------------------------------------
  # Pair combinations
  # ----------------------------------------------------------

  comb <- combn(
    seq_len(k),
    2L
  )

  n_pairs <- ncol(comb)


  pearson <- rep(
    NA_real_,
    n_pairs
  )

  spearman <- rep(
    NA_real_,
    n_pairs
  )

  mean_abs_diff <- rep(
    NA_real_,
    n_pairs
  )

  median_abs_diff <- rep(
    NA_real_,
    n_pairs
  )

  max_abs_diff <- rep(
    NA_real_,
    n_pairs
  )

  exact_equal <- rep(
    FALSE,
    n_pairs
  )


  for (j in seq_len(n_pairs)) {

    r1 <- comb[1L, j]
    r2 <- comb[2L, j]


    a <- as.numeric(
      x[
        r1,
        ,
        drop = TRUE
      ]
    )

    b <- as.numeric(
      x[
        r2,
        ,
        drop = TRUE
      ]
    )


    d <- abs(
      a - b
    )


    pearson[j] <- safe_cor(
      a,
      b,
      method = "pearson"
    )


    spearman[j] <- safe_cor(
      a,
      b,
      method = "spearman"
    )


    mean_abs_diff[j] <- mean(d)

    median_abs_diff[j] <- median(d)

    max_abs_diff[j] <- max(d)

    exact_equal[j] <- identical(
      a,
      b
    )
  }


  pairwise_list[[counter]] <- data.table(
    FinalGeneID =
      rep(
        gene_id,
        n_pairs
      ),

    FinalSymbol =
      rep(
        symbol,
        n_pairs
      ),

    N_SourceTranscriptClusters =
      rep(
        k,
        n_pairs
      ),

    Probe1 =
      probe_ids[
        comb[1L, ]
      ],

    Probe2 =
      probe_ids[
        comb[2L, ]
      ],

    Pearson =
      pearson,

    Spearman =
      spearman,

    MeanAbsoluteDifference =
      mean_abs_diff,

    MedianAbsoluteDifference =
      median_abs_diff,

    MaxAbsoluteDifference =
      max_abs_diff,

    ExactNumericEquality =
      exact_equal
  )


  # ----------------------------------------------------------
  # Sample-wise within-gene spread
  # ----------------------------------------------------------

  sample_spread <- apply(
    x,
    2L,
    function(v) {
      max(v) - min(v)
    }
  )


  # ----------------------------------------------------------
  # Candidate 1: highest mean
  # Diagnostic only
  # ----------------------------------------------------------

  order_max_mean <- order(
    -gene_cluster_stats$MeanExpression,
    -gene_cluster_stats$IQRExpression,
    gene_cluster_stats$SourceProbeID
  )


  candidate_max_mean <- gene_cluster_stats$SourceProbeID[
    order_max_mean[1L]
  ]


  # ----------------------------------------------------------
  # Candidate 2: highest IQR
  # Diagnostic only
  # ----------------------------------------------------------

  order_max_iqr <- order(
    -gene_cluster_stats$IQRExpression,
    -gene_cluster_stats$MeanExpression,
    gene_cluster_stats$SourceProbeID
  )


  candidate_max_iqr <- gene_cluster_stats$SourceProbeID[
    order_max_iqr[1L]
  ]


  # ----------------------------------------------------------
  # Candidate 3: Spearman medoid
  # Diagnostic only
  # ----------------------------------------------------------

  spearman_matrix <- matrix(
    NA_real_,
    nrow = k,
    ncol = k
  )

  diag(spearman_matrix) <- 1


  for (j in seq_len(n_pairs)) {

    r1 <- comb[1L, j]
    r2 <- comb[2L, j]

    spearman_matrix[
      r1,
      r2
    ] <- spearman[j]

    spearman_matrix[
      r2,
      r1
    ] <- spearman[j]
  }


  median_spearman_to_others <- vapply(
    seq_len(k),
    function(r) {

      vals <- spearman_matrix[
        r,
        -r
      ]

      vals <- vals[
        is.finite(vals)
      ]

      if (length(vals) == 0L) {
        return(NA_real_)
      }

      median(vals)
    },
    numeric(1)
  )


  medoid_df <- data.table(
    SourceProbeID =
      probe_ids,

    MedianSpearmanToOthers =
      median_spearman_to_others,

    MeanExpression =
      gene_cluster_stats$MeanExpression,

    IQRExpression =
      gene_cluster_stats$IQRExpression
  )


  medoid_order <- order(
    -medoid_df$MedianSpearmanToOthers,
    -medoid_df$MeanExpression,
    -medoid_df$IQRExpression,
    medoid_df$SourceProbeID,
    na.last = TRUE
  )


  candidate_medoid <- medoid_df$SourceProbeID[
    medoid_order[1L]
  ]


  # ----------------------------------------------------------
  # Gene-level summary
  # ----------------------------------------------------------

  finite_pearson <- pearson[
    is.finite(pearson)
  ]

  finite_spearman <- spearman[
    is.finite(spearman)
  ]


  gene_summary_list[[counter]] <- data.table(
    FinalGeneID =
      gene_id,

    FinalSymbol =
      symbol,

    N_SourceTranscriptClusters =
      k,

    N_Pairs =
      n_pairs,

    N_FinitePearsonPairs =
      length(finite_pearson),

    N_FiniteSpearmanPairs =
      length(finite_spearman),

    MedianPearson =
      safe_quantile(
        pearson,
        0.50
      ),

    MinPearson =
      if (
        length(finite_pearson) > 0L
      ) {
        min(finite_pearson)
      } else {
        NA_real_
      },

    MaxPearson =
      if (
        length(finite_pearson) > 0L
      ) {
        max(finite_pearson)
      } else {
        NA_real_
      },

    MedianSpearman =
      safe_quantile(
        spearman,
        0.50
      ),

    MinSpearman =
      if (
        length(finite_spearman) > 0L
      ) {
        min(finite_spearman)
      } else {
        NA_real_
      },

    MaxSpearman =
      if (
        length(finite_spearman) > 0L
      ) {
        max(finite_spearman)
      } else {
        NA_real_
      },

    MedianPairMeanAbsDiff =
      median(mean_abs_diff),

    MaxPairMeanAbsDiff =
      max(mean_abs_diff),

    MedianPairMaxAbsDiff =
      median(max_abs_diff),

    MaxPairMaxAbsDiff =
      max(max_abs_diff),

    ExactDuplicatePairs =
      sum(exact_equal),

    MinClusterMean =
      min(
        gene_cluster_stats$MeanExpression
      ),

    MaxClusterMean =
      max(
        gene_cluster_stats$MeanExpression
      ),

    ClusterMeanRange =
      max(
        gene_cluster_stats$MeanExpression
      ) -
      min(
        gene_cluster_stats$MeanExpression
      ),

    MinClusterIQR =
      min(
        gene_cluster_stats$IQRExpression
      ),

    MaxClusterIQR =
      max(
        gene_cluster_stats$IQRExpression
      ),

    MedianSampleSpread =
      median(sample_spread),

    Q95SampleSpread =
      safe_quantile(
        sample_spread,
        0.95
      ),

    MaxSampleSpread =
      max(sample_spread)
  )


  candidate_list[[counter]] <- data.table(
    FinalGeneID =
      gene_id,

    FinalSymbol =
      symbol,

    N_SourceTranscriptClusters =
      k,

    Candidate_MaxMean =
      candidate_max_mean,

    Candidate_MaxIQR =
      candidate_max_iqr,

    Candidate_SpearmanMedoid =
      candidate_medoid
  )


  if (
    counter %% 100L == 0L ||
    counter == length(multi_indices)
  ) {

    cat(
      "[INFO] audited ",
      counter,
      " / ",
      length(multi_indices),
      " multi-cluster genes\n",
      sep = ""
    )
  }
}


# ============================================================
# 12. Combine outputs
# ============================================================

pairwise <- rbindlist(
  pairwise_list,
  use.names = TRUE,
  fill = TRUE
)

gene_summary <- rbindlist(
  gene_summary_list,
  use.names = TRUE,
  fill = TRUE
)

candidate_summary <- rbindlist(
  candidate_list,
  use.names = TRUE,
  fill = TRUE
)


check(
  nrow(gene_summary) == 618L,
  "Gene-level concordance summary contains 618 genes"
)

check(
  nrow(candidate_summary) == 618L,
  "Candidate representative summary contains 618 genes"
)


expected_pairs <- sum(
  choose(
    multi_genes$N_SourceTranscriptClusters,
    2L
  )
)


check(
  nrow(pairwise) == expected_pairs,
  "All expected transcript-cluster pairs audited"
)


cat(
  "[INFO] Total pairwise comparisons: ",
  nrow(pairwise),
  "\n",
  sep = ""
)


# ============================================================
# 13. Global concordance summary
# ============================================================

spearman_gene <- gene_summary$MedianSpearman
pearson_gene <- gene_summary$MedianPearson


global_summary <- data.table(
  Metric = c(
    "Multi_cluster_genes",
    "Multi_cluster_source_features",
    "Pairwise_comparisons",

    "Median_gene_median_Spearman",
    "Q05_gene_median_Spearman",
    "Q25_gene_median_Spearman",
    "Q75_gene_median_Spearman",
    "Q95_gene_median_Spearman",

    "Median_gene_median_Pearson",
    "Q05_gene_median_Pearson",
    "Q25_gene_median_Pearson",
    "Q75_gene_median_Pearson",
    "Q95_gene_median_Pearson",

    "Genes_median_Spearman_ge_0.90",
    "Genes_median_Spearman_ge_0.75",
    "Genes_median_Spearman_ge_0.50",
    "Genes_median_Spearman_lt_0.50",

    "Genes_with_any_negative_Spearman_pair",
    "Exact_duplicate_pairs"
  ),

  Value = c(
    618,
    nrow(multi_retained),
    nrow(pairwise),

    safe_quantile(
      spearman_gene,
      0.50
    ),

    safe_quantile(
      spearman_gene,
      0.05
    ),

    safe_quantile(
      spearman_gene,
      0.25
    ),

    safe_quantile(
      spearman_gene,
      0.75
    ),

    safe_quantile(
      spearman_gene,
      0.95
    ),

    safe_quantile(
      pearson_gene,
      0.50
    ),

    safe_quantile(
      pearson_gene,
      0.05
    ),

    safe_quantile(
      pearson_gene,
      0.25
    ),

    safe_quantile(
      pearson_gene,
      0.75
    ),

    safe_quantile(
      pearson_gene,
      0.95
    ),

    sum(
      spearman_gene >= 0.90,
      na.rm = TRUE
    ),

    sum(
      spearman_gene >= 0.75,
      na.rm = TRUE
    ),

    sum(
      spearman_gene >= 0.50,
      na.rm = TRUE
    ),

    sum(
      spearman_gene < 0.50,
      na.rm = TRUE
    ),

    sum(
      gene_summary$MinSpearman < 0,
      na.rm = TRUE
    ),

    sum(
      pairwise$ExactNumericEquality
    )
  )
)


cat("\n============================================================\n")
cat("GLOBAL MULTI-CLUSTER CONCORDANCE SUMMARY\n")
cat("============================================================\n")

print(global_summary)


# ============================================================
# 14. Pairwise Spearman bins
# ============================================================

pairwise[
  ,
  SpearmanDiagnosticBin :=
    fcase(

      is.na(Spearman),
      "NA",

      Spearman >= 0.90,
      ">=0.90",

      Spearman >= 0.75,
      "0.75-0.90",

      Spearman >= 0.50,
      "0.50-0.75",

      Spearman >= 0,
      "0-0.50",

      default =
        "<0"
    )
]


bin_levels <- c(
  ">=0.90",
  "0.75-0.90",
  "0.50-0.75",
  "0-0.50",
  "<0",
  "NA"
)


pairwise_bins <- pairwise[
  ,
  .N,
  by = SpearmanDiagnosticBin
]


pairwise_bins[
  ,
  Order :=
    match(
      SpearmanDiagnosticBin,
      bin_levels
    )
]


setorder(
  pairwise_bins,
  Order
)


pairwise_bins[
  ,
  Order := NULL
]


cat("\n============================================================\n")
cat("PAIRWISE SPEARMAN DIAGNOSTIC DISTRIBUTION\n")
cat("============================================================\n")

print(pairwise_bins)


# ============================================================
# 15. Lowest-concordance genes
# ============================================================

lowest_concordance <- gene_summary[
  order(
    MedianSpearman,
    MinSpearman,
    na.last = TRUE
  )
][
  1:min(.N, 50L)
]


cat("\n============================================================\n")
cat("LOWEST-CONCORDANCE MULTI-CLUSTER GENES — TOP 20\n")
cat("============================================================\n")

print(
  lowest_concordance[
    1:min(.N, 20L),
    .(
      FinalGeneID,
      FinalSymbol,
      N_SourceTranscriptClusters,
      MedianSpearman,
      MinSpearman,
      MedianPearson,
      ClusterMeanRange,
      MedianSampleSpread,
      MaxSampleSpread
    )
  ]
)


# ============================================================
# 16. Candidate-strategy agreement
# ============================================================

candidate_summary[
  ,
  MaxMean_equals_MaxIQR :=
    Candidate_MaxMean ==
      Candidate_MaxIQR
]

candidate_summary[
  ,
  MaxMean_equals_Medoid :=
    Candidate_MaxMean ==
      Candidate_SpearmanMedoid
]

candidate_summary[
  ,
  MaxIQR_equals_Medoid :=
    Candidate_MaxIQR ==
      Candidate_SpearmanMedoid
]

candidate_summary[
  ,
  AllThreeAgree :=
    MaxMean_equals_MaxIQR &
    MaxMean_equals_Medoid
]


candidate_agreement <- data.table(
  Metric = c(
    "Multi_cluster_genes",
    "MaxMean_equals_MaxIQR",
    "MaxMean_equals_Medoid",
    "MaxIQR_equals_Medoid",
    "All_three_agree"
  ),

  N = c(
    nrow(candidate_summary),

    sum(
      candidate_summary$MaxMean_equals_MaxIQR,
      na.rm = TRUE
    ),

    sum(
      candidate_summary$MaxMean_equals_Medoid,
      na.rm = TRUE
    ),

    sum(
      candidate_summary$MaxIQR_equals_Medoid,
      na.rm = TRUE
    ),

    sum(
      candidate_summary$AllThreeAgree,
      na.rm = TRUE
    )
  )
)


cat("\n============================================================\n")
cat("CANDIDATE REPRESENTATIVE-CLUSTER AGREEMENT\n")
cat("============================================================\n")

print(candidate_agreement)


# ============================================================
# 17. Save outputs
# ============================================================

fwrite(
  cluster_stats,
  file.path(
    out_dir,
    "GSE73661_stage2e_transcript_cluster_expression_stats.tsv"
  ),
  sep = "\t"
)


fwrite(
  pairwise,
  file.path(
    out_dir,
    "GSE73661_stage2e_multicluster_pairwise_concordance.tsv"
  ),
  sep = "\t"
)


fwrite(
  gene_summary,
  file.path(
    out_dir,
    "GSE73661_stage2e_multicluster_gene_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  candidate_summary,
  file.path(
    out_dir,
    "GSE73661_stage2e_candidate_representative_clusters.tsv"
  ),
  sep = "\t"
)


fwrite(
  global_summary,
  file.path(
    out_dir,
    "GSE73661_stage2e_global_concordance_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  pairwise_bins,
  file.path(
    out_dir,
    "GSE73661_stage2e_pairwise_spearman_bins.tsv"
  ),
  sep = "\t"
)


fwrite(
  lowest_concordance,
  file.path(
    out_dir,
    "GSE73661_stage2e_lowest_concordance_genes.tsv"
  ),
  sep = "\t"
)


fwrite(
  candidate_agreement,
  file.path(
    out_dir,
    "GSE73661_stage2e_candidate_strategy_agreement.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 18. Final structural checks
# ============================================================

check(
  nrow(cluster_stats) == 19229L,
  "All retained transcript clusters have expression statistics"
)

check(
  nrow(gene_summary) == 618L,
  "All multi-cluster genes have concordance statistics"
)

check(
  nrow(pairwise) == expected_pairs,
  "All expected pairwise comparisons represented"
)

check(
  all(
    gene_summary$N_SourceTranscriptClusters > 1L
  ),
  "Gene summary contains multi-cluster genes only"
)

check(
  identical(
    dim(expr),
    c(33252L, 178L)
  ),
  "Authoritative expression remains unchanged"
)


# ============================================================
# 19. Stage summary
# ============================================================

stage_summary <- data.table(
  Metric = c(
    "Retained_transcript_clusters",
    "Unique_final_genes",
    "Single_cluster_genes",
    "Multi_cluster_genes",
    "Multi_cluster_source_features",
    "Pairwise_comparisons",
    "Exact_duplicate_pairs",
    "Elapsed_seconds"
  ),

  Value = c(
    nrow(retained),
    uniqueN(retained$FinalGeneID),
    nrow(single_genes),
    nrow(multi_genes),
    nrow(multi_retained),
    nrow(pairwise),
    sum(pairwise$ExactNumericEquality),
    round(
      proc.time()[["elapsed"]] - t0,
      2
    )
  )
)


cat("\n============================================================\n")
cat("STAGE 2E SUMMARY\n")
cat("============================================================\n")

print(stage_summary)


fwrite(
  stage_summary,
  file.path(
    out_dir,
    "GSE73661_stage2e_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 20. FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] 19,229 retained transcript clusters audited\n")
cat("[PASS] 17,804 single-cluster genes require no collapse decision\n")
cat("[PASS] 618 multi-cluster genes audited\n")
cat("[PASS] 1,425 transcript clusters evaluated within multi-cluster genes\n")
cat("[PASS] All expected pairwise concordance comparisons completed\n")
cat("[PASS] Pearson and Spearman concordance quantified\n")
cat("[PASS] Baseline-expression offsets quantified\n")
cat("[PASS] Sample-wise within-gene spread quantified\n")
cat("[PASS] MaxMean / MaxIQR / Spearman-medoid candidates recorded\n")
cat("[PASS] No transcript cluster selected\n")
cat("[PASS] No expression values modified\n")
cat("[PASS] No transcript-cluster collapse performed\n")

cat(
  sprintf(
    "[INFO] Total elapsed time: %.2f seconds\n",
    proc.time()[["elapsed"]] - t0
  )
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE73661_STAGE2E_MULTICLUSTER_EXPRESSION_AUDIT\n")
cat("============================================================\n")
