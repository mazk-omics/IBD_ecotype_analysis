suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE57945 <-> GSE93624
# Expression fingerprint crosswalk
#
# Purpose:
#   empirically determine sample/patient overlap between the
#   two RISK-derived datasets.
#
# IMPORTANT:
#   - NO expression normalization is written back to disk
#   - standardized matrices exist only for fingerprinting
#   - NO sample identity is assumed in advance
#   - NO match is accepted solely because it is top-1
#
# Primary fingerprint metrics:
#
#   A. Gene-wise Z fingerprint
#      Each common gene standardized across samples separately
#      within each cohort.
#
#   B. Gene-wise rank fingerprint
#      Each common gene rank-standardized across samples
#      separately within each cohort.
#
# This removes fixed gene-abundance differences and reduces
# sensitivity to the different deposited expression scales.
# ============================================================


root <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression"

base57945 <- file.path(
  root,
  "GSE57945"
)

base93624 <- file.path(
  root,
  "GSE93624"
)

out_dir <- file.path(
  root,
  "99_cross_cohort",
  "GSE57945_vs_GSE93624"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


file57945 <- file.path(
  base57945,
  "prepared/03_canonical_expression",
  "GSE57945_canonical_expression_34368x322.rds"
)

file93624 <- file.path(
  base93624,
  "prepared/03_canonical_expression",
  "GSE93624_canonical_expression_13323x245.rds"
)

meta93624_file <- file.path(
  base93624,
  "prepared/04_sample_metadata",
  "GSE93624_canonical_sample_metadata.tsv"
)

qc93624_file <- file.path(
  base93624,
  "prepared/05_sample_expression_qc",
  "GSE93624_sample_level_expression_qc.tsv"
)


cat("============================================================\n")
cat("GSE57945 <-> GSE93624 expression fingerprint crosswalk\n")
cat("============================================================\n\n")


# ============================================================
# 1. Flexible expression-object loader
# ============================================================

get_expression <- function(x) {

  if (is.matrix(x)) {
    return(x)
  }

  if (
    is.list(x) &&
    !is.null(x$expression)
  ) {
    return(
      as.matrix(
        x$expression
      )
    )
  }

  stop(
    "RDS is neither a matrix nor a list containing $expression"
  )
}


obj57945 <- readRDS(
  file57945
)

obj93624 <- readRDS(
  file93624
)

x57945 <- get_expression(
  obj57945
)

x93624 <- get_expression(
  obj93624
)


stopifnot(
  nrow(x57945) == 34368,
  ncol(x57945) == 322,
  nrow(x93624) == 13323,
  ncol(x93624) == 245,

  !anyDuplicated(
    rownames(x57945)
  ),

  !anyDuplicated(
    rownames(x93624)
  ),

  !anyDuplicated(
    colnames(x57945)
  ),

  !anyDuplicated(
    colnames(x93624)
  ),

  !anyNA(x57945),
  !anyNA(x93624),

  all(is.finite(x57945)),
  all(is.finite(x93624))
)


cat(
  "[PASS] GSE57945:",
  nrow(x57945),
  "x",
  ncol(x57945),
  "\n"
)

cat(
  "[PASS] GSE93624:",
  nrow(x93624),
  "x",
  ncol(x93624),
  "\n"
)


# ============================================================
# 2. Common canonical genes
# ============================================================

common_genes <- intersect(
  rownames(x93624),
  rownames(x57945)
)

stopifnot(
  length(common_genes) > 0
)


# Preserve GSE93624 canonical order.
common_genes <- rownames(x93624)[
  rownames(x93624) %chin%
    common_genes
]


a <- x57945[
  common_genes,
  ,
  drop = FALSE
]

b <- x93624[
  common_genes,
  ,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(a),
    rownames(b)
  )
)


cat("\n============================================================\n")
cat("COMMON GENE UNIVERSE\n")
cat("============================================================\n")

cat(
  "Common canonical symbols:",
  length(common_genes),
  "\n"
)


# ============================================================
# 3. Gene variability audit
#
# For gene-wise standardization, constant genes cannot
# contribute to fingerprinting and are removed.
# ============================================================

sd57945 <- apply(
  a,
  1,
  sd
)

sd93624 <- apply(
  b,
  1,
  sd
)


usable <- (
  is.finite(sd57945) &
  is.finite(sd93624) &
  sd57945 > 0 &
  sd93624 > 0
)


fingerprint_genes <- common_genes[
  usable
]

a_fp <- a[
  usable,
  ,
  drop = FALSE
]

b_fp <- b[
  usable,
  ,
  drop = FALSE
]


cat(
  "Non-variable genes removed:",
  sum(!usable),
  "\n"
)

cat(
  "Fingerprint genes:",
  length(fingerprint_genes),
  "\n"
)


stopifnot(
  nrow(a_fp) ==
    nrow(b_fp),

  nrow(a_fp) > 1000
)


# ============================================================
# 4. Metric A:
#    gene-wise z-standardization within each cohort
# ============================================================

gene_zscore <- function(m) {

  mu <- rowMeans(m)

  s <- apply(
    m,
    1,
    sd
  )

  stopifnot(
    all(
      is.finite(s)
    ),
    all(
      s > 0
    )
  )

  z <- sweep(
    m,
    1,
    mu,
    FUN = "-"
  )

  z <- sweep(
    z,
    1,
    s,
    FUN = "/"
  )

  z
}


cat("\nBuilding gene-wise z fingerprints...\n")

za <- gene_zscore(
  a_fp
)

zb <- gene_zscore(
  b_fp
)


stopifnot(
  all(is.finite(za)),
  all(is.finite(zb))
)


cat(
  "Computing 322 x 245 z-fingerprint correlation matrix...\n"
)

z_cor <- cor(
  za,
  zb,
  method = "pearson"
)


stopifnot(
  all(
    dim(z_cor) ==
      c(322, 245)
  )
)


# ============================================================
# 5. Metric B:
#    gene-wise rank-standardization
#
# For every gene:
#   rank samples within that cohort,
#   then standardize those ranks to mean 0 / sd 1.
#
# This is invariant to monotonic transformations within gene
# and therefore complements the z-based fingerprint.
# ============================================================

gene_rank_standardize <- function(m) {

  out <- t(
    apply(
      m,
      1,
      function(x) {

        r <- rank(
          x,
          ties.method = "average"
        )

        s <- sd(r)

        if (
          !is.finite(s) ||
          s == 0
        ) {
          return(
            rep(
              0,
              length(r)
            )
          )
        }

        (
          r -
          mean(r)
        ) / s
      }
    )
  )

  rownames(out) <- rownames(m)
  colnames(out) <- colnames(m)

  out
}


cat(
  "\nBuilding gene-wise rank fingerprints...\n"
)

ra <- gene_rank_standardize(
  a_fp
)

rb <- gene_rank_standardize(
  b_fp
)


cat(
  "Computing 322 x 245 rank-fingerprint correlation matrix...\n"
)

rank_cor <- cor(
  ra,
  rb,
  method = "pearson"
)


stopifnot(
  all(
    dim(rank_cor) ==
      c(322, 245)
  )
)


# ============================================================
# 6. Extract top-1 and top-2 match for every GSE93624 sample
# ============================================================

top_two_by_column <- function(m) {

  out <- lapply(
    seq_len(ncol(m)),
    function(j) {

      v <- m[, j]

      ord <- order(
        v,
        decreasing = TRUE,
        na.last = NA
      )

      stopifnot(
        length(ord) >= 2
      )

      data.table(
        ColIndex =
          j,

        Top1Index =
          ord[1],

        Top1 =
          v[
            ord[1]
          ],

        Top2Index =
          ord[2],

        Top2 =
          v[
            ord[2]
          ],

        Gap =
          v[
            ord[1]
          ] -
          v[
            ord[2]
          ]
      )
    }
  )

  rbindlist(out)
}


z_top <- top_two_by_column(
  z_cor
)

r_top <- top_two_by_column(
  rank_cor
)


# ============================================================
# 7. Reciprocal nearest-neighbour audit
# ============================================================

# For each GSE57945 sample, identify its best GSE93624 match.
z_reverse_best <- apply(
  z_cor,
  1,
  which.max
)

rank_reverse_best <- apply(
  rank_cor,
  1,
  which.max
)


crosswalk <- data.table(
  GSE93624_Sample =
    colnames(b_fp),

  GSE57945_Top1_Z =
    rownames(z_cor)[
      z_top$Top1Index
    ],

  Z_Top1 =
    z_top$Top1,

  GSE57945_Top2_Z =
    rownames(z_cor)[
      z_top$Top2Index
    ],

  Z_Top2 =
    z_top$Top2,

  Z_Gap =
    z_top$Gap,

  GSE57945_Top1_Rank =
    rownames(rank_cor)[
      r_top$Top1Index
    ],

  Rank_Top1 =
    r_top$Top1,

  GSE57945_Top2_Rank =
    rownames(rank_cor)[
      r_top$Top2Index
    ],

  Rank_Top2 =
    r_top$Top2,

  Rank_Gap =
    r_top$Gap
)


crosswalk[
  ,
  SameTop1AcrossMetrics :=
    GSE57945_Top1_Z ==
      GSE57945_Top1_Rank
]


crosswalk[
  ,
  Reciprocal_Z :=
    vapply(
      seq_len(.N),
      function(j) {

        i <- z_top$Top1Index[j]

        z_reverse_best[i] == j

      },
      logical(1)
    )
]


crosswalk[
  ,
  Reciprocal_Rank :=
    vapply(
      seq_len(.N),
      function(j) {

        i <- r_top$Top1Index[j]

        rank_reverse_best[i] == j

      },
      logical(1)
    )
]


# ------------------------------------------------------------
# Diagnostic consensus candidate.
#
# This is NOT yet a declared patient identity match.
# ------------------------------------------------------------

crosswalk[
  ,
  ExpressionConcordantCandidate :=
    SameTop1AcrossMetrics &
    Reciprocal_Z &
    Reciprocal_Rank
]


# ============================================================
# 8. Top-1 target multiplicity
#
# If GSE93624 is mostly a subset of GSE57945, a strong
# one-to-one fingerprint should usually produce unique targets.
# ============================================================

crosswalk[
  ,
  Z_Top1_TargetMultiplicity :=
    .N,
  by = GSE57945_Top1_Z
]

crosswalk[
  ,
  Rank_Top1_TargetMultiplicity :=
    .N,
  by = GSE57945_Top1_Rank
]


# ============================================================
# 9. Attach GSE93624 metadata + QC flags
# ============================================================

meta <- fread(
  meta93624_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

qc <- fread(
  qc93624_file,
  colClasses = "character",
  na.strings = c("", "NA")
)


meta_idx <- match(
  crosswalk$GSE93624_Sample,
  meta$ExpressionSampleID
)

qc_idx <- match(
  crosswalk$GSE93624_Sample,
  qc$ExpressionSampleID
)


stopifnot(
  !anyNA(meta_idx),
  !anyNA(qc_idx)
)


crosswalk[
  ,
  `:=`(
    GSE93624_GSM =
      meta$GSM[
        meta_idx
      ],

    Diagnosis =
      meta$Diagnosis[
        meta_idx
      ],

    Progression3yr =
      meta$Progression3yr[
        meta_idx
      ],

    Sex =
      meta$Sex[
        meta_idx
      ],

    AgeAtDiagnosis =
      meta$AgeAtDiagnosis[
        meta_idx
      ],

    ParisAge =
      meta$ParisAge[
        meta_idx
      ],

    Ancestry =
      meta$Ancestry[
        meta_idx
      ],

    GSE93624_QCFlag =
      qc$AnyDiagnosticFlag[
        qc_idx
      ]
  )
]


# ============================================================
# 10. Summary
# ============================================================

cat("\n============================================================\n")
cat("FINGERPRINT CROSSWALK SUMMARY\n")
cat("============================================================\n")


summary <- data.table(
  Metric = c(
    "GSE93624_samples",
    "GSE57945_samples",
    "Common_canonical_genes",
    "Fingerprint_genes",
    "Unique_Z_top1_targets",
    "Unique_rank_top1_targets",
    "Same_top1_across_metrics",
    "Reciprocal_Z",
    "Reciprocal_rank",
    "Expression_concordant_candidates",
    "GSE93624_QC_flagged_candidates"
  ),

  Value = c(
    ncol(b_fp),
    ncol(a_fp),
    length(common_genes),
    length(fingerprint_genes),

    uniqueN(
      crosswalk$GSE57945_Top1_Z
    ),

    uniqueN(
      crosswalk$GSE57945_Top1_Rank
    ),

    sum(
      crosswalk$SameTop1AcrossMetrics
    ),

    sum(
      crosswalk$Reciprocal_Z
    ),

    sum(
      crosswalk$Reciprocal_Rank
    ),

    sum(
      crosswalk$ExpressionConcordantCandidate
    ),

    sum(
      crosswalk$ExpressionConcordantCandidate &
      crosswalk$GSE93624_QCFlag ==
        "TRUE",
      na.rm = TRUE
    )
  )
)

print(
  summary
)


# ============================================================
# 11. Score distributions
# ============================================================

score_summary <- data.table(
  Metric = c(
    "Z_Top1",
    "Z_Gap",
    "Rank_Top1",
    "Rank_Gap"
  ),

  Min = c(
    min(crosswalk$Z_Top1),
    min(crosswalk$Z_Gap),
    min(crosswalk$Rank_Top1),
    min(crosswalk$Rank_Gap)
  ),

  Q05 = c(
    quantile(
      crosswalk$Z_Top1,
      0.05,
      names = FALSE
    ),

    quantile(
      crosswalk$Z_Gap,
      0.05,
      names = FALSE
    ),

    quantile(
      crosswalk$Rank_Top1,
      0.05,
      names = FALSE
    ),

    quantile(
      crosswalk$Rank_Gap,
      0.05,
      names = FALSE
    )
  ),

  Median = c(
    median(crosswalk$Z_Top1),
    median(crosswalk$Z_Gap),
    median(crosswalk$Rank_Top1),
    median(crosswalk$Rank_Gap)
  ),

  Q95 = c(
    quantile(
      crosswalk$Z_Top1,
      0.95,
      names = FALSE
    ),

    quantile(
      crosswalk$Z_Gap,
      0.95,
      names = FALSE
    ),

    quantile(
      crosswalk$Rank_Top1,
      0.95,
      names = FALSE
    ),

    quantile(
      crosswalk$Rank_Gap,
      0.95,
      names = FALSE
    )
  ),

  Max = c(
    max(crosswalk$Z_Top1),
    max(crosswalk$Z_Gap),
    max(crosswalk$Rank_Top1),
    max(crosswalk$Rank_Gap)
  )
)


cat("\n============================================================\n")
cat("FINGERPRINT SCORE DISTRIBUTIONS\n")
cat("============================================================\n")

print(
  score_summary
)


# ============================================================
# 12. Non-consensus / non-reciprocal samples
# ============================================================

non_consensus <- crosswalk[
  ExpressionConcordantCandidate == FALSE
]


cat("\n============================================================\n")
cat("NON-CONSENSUS / NON-RECIPROCAL SAMPLES\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(non_consensus),
  "\n"
)

if (nrow(non_consensus) > 0) {

  print(
    non_consensus[
      order(
        Z_Gap,
        Rank_Gap
      ),
      .(
        GSE93624_Sample,
        GSE93624_GSM,
        Diagnosis,
        Progression3yr,

        GSE57945_Top1_Z,
        Z_Top1,
        Z_Gap,
        Reciprocal_Z,

        GSE57945_Top1_Rank,
        Rank_Top1,
        Rank_Gap,
        Reciprocal_Rank,

        SameTop1AcrossMetrics,
        Z_Top1_TargetMultiplicity,
        Rank_Top1_TargetMultiplicity,

        GSE93624_QCFlag
      )
    ]
  )
}


# ============================================================
# 13. Weakest consensus candidates
#
# Useful for seeing whether there is a clean separation or
# whether some putative matches are only marginal.
# ============================================================

consensus <- crosswalk[
  ExpressionConcordantCandidate == TRUE
]


cat("\n============================================================\n")
cat("20 WEAKEST EXPRESSION-CONCORDANT CANDIDATES\n")
cat("============================================================\n")

if (nrow(consensus) > 0) {

  print(
    head(
      consensus[
        order(
          Z_Gap,
          Rank_Gap
        ),
        .(
          GSE93624_Sample,
          GSE93624_GSM,
          Diagnosis,
          Progression3yr,

          GSE57945_Top1_Z,

          Z_Top1,
          Z_Gap,

          Rank_Top1,
          Rank_Gap,

          GSE93624_QCFlag
        )
      ],
      20
    )
  )
}


# ============================================================
# 14. Save outputs
# ============================================================

fwrite(
  crosswalk,
  file.path(
    out_dir,
    "GSE57945_GSE93624_expression_fingerprint_crosswalk.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  summary,
  file.path(
    out_dir,
    "GSE57945_GSE93624_fingerprint_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  score_summary,
  file.path(
    out_dir,
    "GSE57945_GSE93624_fingerprint_score_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  non_consensus,
  file.path(
    out_dir,
    "GSE57945_GSE93624_non_consensus_samples.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  consensus,
  file.path(
    out_dir,
    "GSE57945_GSE93624_expression_concordant_candidates.tsv"
  ),
  sep = "\t",
  na = ""
)


saveRDS(
  list(
    common_genes =
      common_genes,

    fingerprint_genes =
      fingerprint_genes,

    z_correlation =
      z_cor,

    rank_correlation =
      rank_cor
  ),
  file.path(
    out_dir,
    "GSE57945_GSE93624_fingerprint_correlation_matrices.rds"
  ),
  compress = FALSE
)


# ============================================================
# 15. Final
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] Cross-cohort expression fingerprints computed\n")
cat("[PASS] No sample identity assumed from accession membership\n")
cat("[PASS] No expression matrix modified\n")
cat("[PASS] Candidate matches require two-metric concordance + reciprocity\n")

cat("\nFINAL STATUS:\n")
cat("PASS_GSE57945_GSE93624_EXPRESSION_FINGERPRINT_CROSSWALK\n")
cat("============================================================\n")
