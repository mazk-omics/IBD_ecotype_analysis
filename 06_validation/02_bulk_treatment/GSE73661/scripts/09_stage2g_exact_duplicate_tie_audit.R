suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 2G
# Exact-duplicate audit for the 120 remaining Stage-2F ties
#
# Frozen inputs:
#   Stage 1 authoritative log2-RMA expression
#   Stage 2F annotation/design candidate table
#   Stage 2F multi-cluster annotation/design feature table
#
# Purpose:
#   Determine whether tied annotation/design candidates produce
#   EXACTLY identical expression vectors across all 178 samples.
#
# Rule:
#   EXACT identical tied candidates:
#       deterministic canonical candidate may be defined as the
#       numerically smallest SourceProbeID because the resulting
#       expression vector is identical regardless of choice.
#
#   NON-identical tied candidates:
#       remain unresolved.
#
# IMPORTANT:
#   - correlation does NOT resolve ties
#   - mean/IQR does NOT resolve ties
#   - no expression values are changed
#   - no final representative mapping is frozen here
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

candidate_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2f_annotation_design_candidates.tsv"
)

feature_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2f_multicluster_annotation_design_features.tsv"
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

  cat(
    "[PASS] ",
    message,
    "\n",
    sep = ""
  )
}


safe_cor <- function(x, y, method = "spearman") {

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


split_candidates <- function(x) {

  if (
    length(x) != 1L ||
    is.na(x) ||
    !nzchar(x)
  ) {
    return(character())
  }

  y <- strsplit(
    x,
    ";",
    fixed = TRUE
  )[[1]]

  y <- trimws(y)

  y[
    nzchar(y)
  ]
}


deterministic_probe <- function(ids) {

  ids <- as.character(ids)

  numeric_ids <- suppressWarnings(
    as.numeric(ids)
  )

  if (
    length(ids) == length(numeric_ids) &&
    all(is.finite(numeric_ids))
  ) {

    return(
      ids[
        order(
          numeric_ids,
          ids
        )[1L]
      ]
    )
  }

  sort(ids)[1L]
}


cat("============================================================\n")
cat("GSE73661 Stage 2G: exact-duplicate tie audit\n")
cat("============================================================\n\n")

t0 <- proc.time()[["elapsed"]]


# ============================================================
# 2. Load frozen inputs
# ============================================================

for (f in c(
  expr_file,
  candidate_file,
  feature_file
)) {

  check(
    file.exists(f),
    paste0(
      "Input exists: ",
      basename(f)
    )
  )
}


expr <- readRDS(
  expr_file
)

candidates <- fread(
  candidate_file,
  showProgress = FALSE
)

features <- fread(
  feature_file,
  showProgress = FALSE
)


check(
  identical(
    dim(expr),
    c(
      33252L,
      178L
    )
  ),
  "Authoritative expression = 33,252 x 178"
)


# ============================================================
# 3. Freeze identifier types
# ============================================================

candidates[
  ,
  FinalGeneID :=
    as.character(FinalGeneID)
]

candidates[
  ,
  FinalSymbol :=
    as.character(FinalSymbol)
]

candidates[
  ,
  TiedCandidates :=
    as.character(TiedCandidates)
]


features[
  ,
  FinalGeneID :=
    as.character(FinalGeneID)
]

features[
  ,
  FinalSymbol :=
    as.character(FinalSymbol)
]

features[
  ,
  SourceProbeID :=
    as.character(SourceProbeID)
]


expr_probe_ids <- as.character(
  rownames(expr)
)


# ============================================================
# 4. Recover the 120 Stage-2F ties
# ============================================================

ties <- candidates[
  N_AnnotationDesignCandidates > 1L
]


check(
  nrow(candidates) == 618L,
  "Stage-2F candidate table contains 618 genes"
)


check(
  nrow(ties) == 120L,
  "Exactly 120 annotation/design ties recovered"
)


# ============================================================
# 5. Containers
# ============================================================

tie_summary_list <- vector(
  "list",
  nrow(ties)
)

candidate_long_list <- vector(
  "list",
  nrow(ties)
)


# ============================================================
# 6. Audit each tied gene
# ============================================================

cat(
  "Auditing exact expression identity across 120 tied genes...\n"
)


for (i in seq_len(nrow(ties))) {

  gene_id <- ties$FinalGeneID[i]

  symbol <- ties$FinalSymbol[i]

  tied_ids <- split_candidates(
    ties$TiedCandidates[i]
  )


  expected_n <- as.integer(
    ties$N_AnnotationDesignCandidates[i]
  )


  if (
    length(tied_ids) != expected_n
  ) {
    stop(
      paste0(
        "[FAIL] Candidate-count mismatch for ",
        gene_id,
        " / ",
        symbol
      ),
      call. = FALSE
    )
  }


  if (
    length(tied_ids) < 2L
  ) {
    stop(
      paste0(
        "[FAIL] Expected >=2 tied candidates for ",
        gene_id
      ),
      call. = FALSE
    )
  }


  idx <- match(
    tied_ids,
    expr_probe_ids
  )


  if (anyNA(idx)) {
    stop(
      paste0(
        "[FAIL] Missing tied probe in authoritative expression for ",
        gene_id
      ),
      call. = FALSE
    )
  }


  x <- expr[
    idx,
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Exact vector equality
  # ----------------------------------------------------------

  reference_vector <- as.numeric(
    x[
      1L,
      ,
      drop = TRUE
    ]
  )


  identical_to_reference <- vapply(
    seq_len(nrow(x)),
    function(j) {

      identical(
        as.numeric(
          x[
            j,
            ,
            drop = TRUE
          ]
        ),
        reference_vector
      )
    },
    logical(1)
  )


  all_exact <- all(
    identical_to_reference
  )


  # ----------------------------------------------------------
  # Pairwise diagnostics
  # ----------------------------------------------------------

  comb <- combn(
    seq_len(nrow(x)),
    2L
  )


  n_pairs <- ncol(comb)


  pair_spearman <- rep(
    NA_real_,
    n_pairs
  )

  pair_pearson <- rep(
    NA_real_,
    n_pairs
  )

  pair_mean_abs_diff <- rep(
    NA_real_,
    n_pairs
  )

  pair_max_abs_diff <- rep(
    NA_real_,
    n_pairs
  )

  pair_exact <- rep(
    FALSE,
    n_pairs
  )


  for (j in seq_len(n_pairs)) {

    a <- as.numeric(
      x[
        comb[1L, j],
        ,
        drop = TRUE
      ]
    )

    b <- as.numeric(
      x[
        comb[2L, j],
        ,
        drop = TRUE
      ]
    )


    d <- abs(
      a - b
    )


    pair_spearman[j] <- safe_cor(
      a,
      b,
      method = "spearman"
    )


    pair_pearson[j] <- safe_cor(
      a,
      b,
      method = "pearson"
    )


    pair_mean_abs_diff[j] <- mean(d)

    pair_max_abs_diff[j] <- max(d)

    pair_exact[j] <- identical(
      a,
      b
    )
  }


  finite_spearman <- pair_spearman[
    is.finite(
      pair_spearman
    )
  ]


  finite_pearson <- pair_pearson[
    is.finite(
      pair_pearson
    )
  ]


  min_spearman <- if (
    length(finite_spearman) > 0L
  ) {
    min(finite_spearman)
  } else {
    NA_real_
  }


  median_spearman <- if (
    length(finite_spearman) > 0L
  ) {
    median(finite_spearman)
  } else {
    NA_real_
  }


  min_pearson <- if (
    length(finite_pearson) > 0L
  ) {
    min(finite_pearson)
  } else {
    NA_real_
  }


  # ----------------------------------------------------------
  # Deterministic candidate ONLY if values are exactly equal
  # ----------------------------------------------------------

  deterministic_candidate <- if (
    all_exact
  ) {

    deterministic_probe(
      tied_ids
    )

  } else {

    NA_character_
  }


  # ----------------------------------------------------------
  # Resolution class
  # ----------------------------------------------------------

  resolution_class <- if (
    all_exact
  ) {

    "EXACT_NUMERIC_DUPLICATE_TIE"

  } else if (
    is.na(min_spearman)
  ) {

    "NONIDENTICAL_CORRELATION_UNDEFINED"

  } else if (
    min_spearman >= 0.90
  ) {

    "NONIDENTICAL_ALL_PAIRS_SPEARMAN_GE_0.90"

  } else if (
    min_spearman >= 0.75
  ) {

    "NONIDENTICAL_MIN_SPEARMAN_0.75_0.90"

  } else if (
    min_spearman >= 0.50
  ) {

    "NONIDENTICAL_MIN_SPEARMAN_0.50_0.75"

  } else if (
    min_spearman >= 0
  ) {

    "NONIDENTICAL_MIN_SPEARMAN_0_0.50"

  } else {

    "NONIDENTICAL_NEGATIVE_SPEARMAN_PAIR"
  }


  # ----------------------------------------------------------
  # Gene-level output
  # ----------------------------------------------------------

  tie_summary_list[[i]] <- data.table(
    FinalGeneID =
      gene_id,

    FinalSymbol =
      symbol,

    N_TiedCandidates =
      length(tied_ids),

    N_Pairs =
      n_pairs,

    TiedCandidates =
      paste(
        tied_ids,
        collapse = ";"
      ),

    AllTiedCandidatesExact =
      all_exact,

    ExactPairs =
      sum(pair_exact),

    MinSpearmanTied =
      min_spearman,

    MedianSpearmanTied =
      median_spearman,

    MinPearsonTied =
      min_pearson,

    MedianPairMeanAbsDiff =
      median(
        pair_mean_abs_diff
      ),

    MaxPairMeanAbsDiff =
      max(
        pair_mean_abs_diff
      ),

    MaxAbsoluteDifference =
      max(
        pair_max_abs_diff
      ),

    TieResolutionClass =
      resolution_class,

    DeterministicCandidateIfExact =
      deterministic_candidate
  )


  # ----------------------------------------------------------
  # Long candidate table
  # ----------------------------------------------------------

  candidate_long_list[[i]] <- data.table(
    FinalGeneID =
      rep(
        gene_id,
        length(tied_ids)
      ),

    FinalSymbol =
      rep(
        symbol,
        length(tied_ids)
      ),

    SourceProbeID =
      tied_ids,

    TiedCandidate =
      TRUE,

    GeneAllTiedCandidatesExact =
      all_exact,

    DeterministicCandidateIfExact =
      deterministic_candidate
  )


  if (
    i %% 25L == 0L ||
    i == nrow(ties)
  ) {

    cat(
      "[INFO] audited ",
      i,
      " / ",
      nrow(ties),
      " tied genes\n",
      sep = ""
    )
  }
}


# ============================================================
# 7. Combine outputs
# ============================================================

tie_summary <- rbindlist(
  tie_summary_list,
  use.names = TRUE
)

candidate_long <- rbindlist(
  candidate_long_list,
  use.names = TRUE
)


check(
  nrow(tie_summary) == 120L,
  "Tie summary contains all 120 genes"
)


check(
  uniqueN(
    candidate_long$FinalGeneID
  ) == 120L,
  "Candidate-long table covers all 120 tied genes"
)


# ============================================================
# 8. Exact-equality guards
# ============================================================

exact_rows <- tie_summary[
  AllTiedCandidatesExact == TRUE
]


nonexact_rows <- tie_summary[
  AllTiedCandidatesExact == FALSE
]


check(
  all(
    exact_rows$MaxAbsoluteDifference == 0
  ),
  "Every exact tie has maximum absolute difference = 0"
)


check(
  all(
    exact_rows$ExactPairs ==
      exact_rows$N_Pairs
  ),
  "Every exact tie has all pairs numerically identical"
)


check(
  !anyNA(
    exact_rows$DeterministicCandidateIfExact
  ),
  "Every exact tie has deterministic canonical candidate"
)


check(
  all(
    is.na(
      nonexact_rows$DeterministicCandidateIfExact
    )
  ),
  "No non-exact tie receives deterministic candidate"
)


# ============================================================
# 9. Resolution summary
# ============================================================

resolution_summary <- tie_summary[
  ,
  .N,
  by = TieResolutionClass
][
  order(
    -N
  )
]


n_exact <- nrow(
  exact_rows
)

n_nonexact <- nrow(
  nonexact_rows
)


summary <- data.table(
  Metric = c(
    "Stage2F_remaining_ties",
    "Exact_numeric_duplicate_ties",
    "Nonidentical_remaining_ties",
    "Exact_ties_fraction",
    "Nonexact_ties_fraction"
  ),

  Value = c(
    120,
    n_exact,
    n_nonexact,
    n_exact / 120,
    n_nonexact / 120
  )
)


cat("\n============================================================\n")
cat("EXACT-DUPLICATE TIE AUDIT SUMMARY\n")
cat("============================================================\n")

print(
  summary
)


cat("\n============================================================\n")
cat("NON-EXACT TIE CONCORDANCE CLASSES\n")
cat("============================================================\n")

print(
  resolution_summary
)


# ============================================================
# 10. Detailed unresolved feature table
# ============================================================

unresolved_candidate_long <- candidate_long[
  GeneAllTiedCandidatesExact == FALSE
]


feature_idx <- match(
  paste(
    unresolved_candidate_long$FinalGeneID,
    unresolved_candidate_long$SourceProbeID,
    sep = "||"
  ),
  paste(
    features$FinalGeneID,
    features$SourceProbeID,
    sep = "||"
  )
)


check(
  !anyNA(feature_idx),
  "Every unresolved tied candidate found in Stage-2F feature table"
)


unresolved_details <- cbind(
  unresolved_candidate_long,
  features[
    feature_idx,
    setdiff(
      names(features),
      c(
        "FinalGeneID",
        "FinalSymbol",
        "SourceProbeID"
      )
    ),
    with = FALSE
  ]
)


# ============================================================
# 11. Print lowest-concordance unresolved ties
# ============================================================

lowest_unresolved <- nonexact_rows[
  order(
    MinSpearmanTied,
    MedianSpearmanTied,
    na.last = TRUE
  )
]


cat("\n============================================================\n")
cat("LOWEST-CONCORDANCE NON-EXACT TIES — TOP 30\n")
cat("============================================================\n")

print(
  head(
    lowest_unresolved[
      ,
      .(
        FinalGeneID,
        FinalSymbol,
        N_TiedCandidates,
        MinSpearmanTied,
        MedianSpearmanTied,
        MinPearsonTied,
        MedianPairMeanAbsDiff,
        MaxPairMeanAbsDiff,
        MaxAbsoluteDifference,
        TiedCandidates
      )
    ],
    30L
  )
)


# ============================================================
# 12. Save outputs
# ============================================================

fwrite(
  tie_summary,
  file.path(
    out_dir,
    "GSE73661_stage2g_exact_duplicate_tie_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  candidate_long,
  file.path(
    out_dir,
    "GSE73661_stage2g_tied_candidates_long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  unresolved_details,
  file.path(
    out_dir,
    "GSE73661_stage2g_nonexact_tied_candidate_details.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  resolution_summary,
  file.path(
    out_dir,
    "GSE73661_stage2g_resolution_class_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  summary,
  file.path(
    out_dir,
    "GSE73661_stage2g_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  lowest_unresolved,
  file.path(
    out_dir,
    "GSE73661_stage2g_lowest_concordance_nonexact_ties.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


# ============================================================
# 13. Final structural checks
# ============================================================

check(
  n_exact + n_nonexact == 120L,
  "Exact + non-exact ties account for all 120 genes"
)


check(
  nrow(tie_summary) == 120L,
  "All Stage-2F ties audited"
)


check(
  identical(
    dim(expr),
    c(
      33252L,
      178L
    )
  ),
  "Authoritative expression remains unchanged"
)


# ============================================================
# 14. FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] All 120 Stage-2F annotation/design ties audited\n")
cat("[PASS] Exact numeric equality assessed across all 178 samples\n")
cat("[PASS] Exact ties distinguished from merely correlated ties\n")
cat("[PASS] Deterministic SourceProbeID assigned only to exact ties\n")
cat("[PASS] Non-exact ties remain unresolved\n")
cat("[PASS] No MaxMean/MaxIQR/correlation-based selection performed\n")
cat("[PASS] No expression values modified\n")
cat("[PASS] No final representative mapping frozen\n")

cat(
  sprintf(
    "[INFO] Total elapsed time: %.2f seconds\n",
    proc.time()[["elapsed"]] - t0
  )
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE73661_STAGE2G_EXACT_DUPLICATE_TIE_AUDIT\n")
cat("============================================================\n")
