suppressPackageStartupMessages({
  library(data.table)
})

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
cat("GSE57945 sample-level expression QC\n")
cat("============================================================\n\n")

expr <- readRDS(expr_file)
meta <- readRDS(meta_file)

stopifnot(
  nrow(expr) == 34368,
  ncol(expr) == 322,
  identical(colnames(expr), meta$SampleID),
  all(is.finite(expr)),
  all(expr >= 0)
)

# ------------------------------------------------------------
# 1. Per-sample descriptive QC
# ------------------------------------------------------------

sample_qc <- data.table(
  SampleID = colnames(expr),

  DetectedGenes =
    colSums(expr > 0),

  ZeroFraction =
    colMeans(expr == 0),

  TotalExpression =
    colSums(expr),

  MedianExpression =
    apply(expr, 2, median),

  MeanExpression =
    colMeans(expr),

  MaxExpression =
    apply(expr, 2, max)
)

sample_qc[
  ,
  Log2TotalExpression :=
    log2(TotalExpression + 1)
]

sample_qc <- merge(
  sample_qc,
  meta[
    ,
    .(
      SampleID,
      GSM,
      Disease,
      DiseaseSubtype,
      Sex,
      AgeAtDiagnosis
    )
  ],
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

sample_qc <- sample_qc[
  match(colnames(expr), SampleID)
]

# ------------------------------------------------------------
# 2. Robust univariate outlier flags
# ------------------------------------------------------------

robust_z <- function(x) {

  med <- median(x, na.rm = TRUE)
  madv <- mad(
    x,
    center = med,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (is.na(madv) || madv == 0) {
    return(rep(0, length(x)))
  }

  (x - med) / madv
}

sample_qc[
  ,
  Z_DetectedGenes :=
    robust_z(DetectedGenes)
]

sample_qc[
  ,
  Z_ZeroFraction :=
    robust_z(ZeroFraction)
]

sample_qc[
  ,
  Z_Log2Total :=
    robust_z(Log2TotalExpression)
]

sample_qc[
  ,
  UnivariateQCFlag :=
    abs(Z_DetectedGenes) > 5 |
    abs(Z_ZeroFraction) > 5 |
    abs(Z_Log2Total) > 5
]


# ------------------------------------------------------------
# 3. PCA QC on log2(expression + 1)
#
# Diagnostic only:
# no normalization, batch correction, or sample filtering.
# ------------------------------------------------------------

log_expr <- log2(expr + 1)

gene_var <- apply(
  log_expr,
  1,
  var
)

top_n <- min(
  2000L,
  sum(
    is.finite(gene_var) &
    gene_var > 0
  )
)

top_genes <- names(
  sort(
    gene_var,
    decreasing = TRUE
  )
)[seq_len(top_n)]

pca <- prcomp(
  t(log_expr[top_genes, , drop = FALSE]),
  center = TRUE,
  scale. = FALSE
)

pcs_to_keep <- min(
  10L,
  ncol(pca$x)
)

pca_dt <- data.table(
  SampleID = rownames(pca$x),
  pca$x[, seq_len(pcs_to_keep), drop = FALSE]
)

pca_dt <- merge(
  pca_dt,
  meta[
    ,
    .(
      SampleID,
      GSM,
      Disease,
      DiseaseSubtype,
      Sex,
      AgeAtDiagnosis
    )
  ],
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

pca_dt <- pca_dt[
  match(colnames(expr), SampleID)
]


# ------------------------------------------------------------
# 4. PCA robust-distance diagnostic
# ------------------------------------------------------------

pc_cols <- paste0(
  "PC",
  seq_len(
    min(5L, pcs_to_keep)
  )
)

pc_mat <- as.matrix(
  pca_dt[
    ,
    ..pc_cols
  ]
)

pc_center <- apply(
  pc_mat,
  2,
  median
)

pc_scale <- apply(
  pc_mat,
  2,
  mad,
  constant = 1.4826
)

pc_scale[
  is.na(pc_scale) |
  pc_scale == 0
] <- 1

pc_z <- sweep(
  pc_mat,
  2,
  pc_center,
  "-"
)

pc_z <- sweep(
  pc_z,
  2,
  pc_scale,
  "/"
)

pca_dt[
  ,
  RobustPCDistance :=
    sqrt(
      rowSums(pc_z^2)
    )
]

pca_threshold <- median(
  pca_dt$RobustPCDistance
) + 5 * mad(
  pca_dt$RobustPCDistance,
  constant = 1.4826
)

pca_dt[
  ,
  PCAQCFlag :=
    RobustPCDistance > pca_threshold
]


# ------------------------------------------------------------
# 5. Combine flags
# ------------------------------------------------------------

sample_qc[
  ,
  PCAQCFlag :=
    pca_dt$PCAQCFlag
]

sample_qc[
  ,
  AnyQCFlag :=
    UnivariateQCFlag |
    PCAQCFlag
]


# ------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("SAMPLE QC SUMMARY\n")
cat("============================================================\n")

print(
  sample_qc[
    ,
    .(
      Metric = c(
        "DetectedGenes",
        "ZeroFraction",
        "Log2TotalExpression"
      ),

      Min = c(
        min(DetectedGenes),
        min(ZeroFraction),
        min(Log2TotalExpression)
      ),

      Q1 = c(
        quantile(DetectedGenes, 0.25),
        quantile(ZeroFraction, 0.25),
        quantile(Log2TotalExpression, 0.25)
      ),

      Median = c(
        median(DetectedGenes),
        median(ZeroFraction),
        median(Log2TotalExpression)
      ),

      Q3 = c(
        quantile(DetectedGenes, 0.75),
        quantile(ZeroFraction, 0.75),
        quantile(Log2TotalExpression, 0.75)
      ),

      Max = c(
        max(DetectedGenes),
        max(ZeroFraction),
        max(Log2TotalExpression)
      )
    )
  ]
)

cat("\n============================================================\n")
cat("QC FLAGS\n")
cat("============================================================\n")

cat(
  "Univariate flagged:",
  sum(sample_qc$UnivariateQCFlag),
  "\n"
)

cat(
  "PCA flagged:",
  sum(sample_qc$PCAQCFlag),
  "\n"
)

cat(
  "Any QC flag:",
  sum(sample_qc$AnyQCFlag),
  "\n"
)

if (any(sample_qc$AnyQCFlag)) {

  cat("\nFlagged samples:\n")

  print(
    sample_qc[
      AnyQCFlag == TRUE,
      .(
        SampleID,
        GSM,
        Disease,
        DiseaseSubtype,
        DetectedGenes,
        ZeroFraction,
        Log2TotalExpression,
        Z_DetectedGenes,
        Z_ZeroFraction,
        Z_Log2Total,
        UnivariateQCFlag,
        PCAQCFlag
      )
    ]
  )
}


# ------------------------------------------------------------
# 7. PCA variance explained
# ------------------------------------------------------------

variance_explained <-
  pca$sdev^2 /
  sum(pca$sdev^2)

pca_variance <- data.table(
  PC = paste0(
    "PC",
    seq_along(variance_explained)
  ),
  VarianceExplained =
    variance_explained,
  CumulativeVariance =
    cumsum(variance_explained)
)

cat("\n============================================================\n")
cat("PCA VARIANCE EXPLAINED\n")
cat("============================================================\n")

print(
  head(
    pca_variance,
    10
  )
)


# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------

fwrite(
  sample_qc,
  file.path(
    out_dir,
    "GSE57945_sample_expression_qc.tsv"
  ),
  sep = "\t"
)

fwrite(
  pca_dt,
  file.path(
    out_dir,
    "GSE57945_expression_PCA_scores.tsv"
  ),
  sep = "\t"
)

fwrite(
  pca_variance,
  file.path(
    out_dir,
    "GSE57945_expression_PCA_variance.tsv"
  ),
  sep = "\t"
)

saveRDS(
  pca,
  file.path(
    out_dir,
    "GSE57945_expression_PCA.rds"
  )
)


cat("\nIMPORTANT:\n")
cat(
  "- QC flags are diagnostic only.\n",
  "- No sample was removed.\n",
  "- PCA used log2(expression + 1) only for QC visualization/diagnostics.\n",
  "- No normalization or batch correction was applied to the canonical matrix.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_SAMPLE_LEVEL_EXPRESSION_QC_COMPLETED\n")
cat("============================================================\n")
