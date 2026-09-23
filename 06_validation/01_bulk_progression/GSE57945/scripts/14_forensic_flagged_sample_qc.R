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

cat("============================================================\n")
cat("GSE57945 flagged-sample forensic QC\n")
cat("============================================================\n\n")

expr <- readRDS(expr_file)
meta <- readRDS(meta_file)

qc <- fread(qc_file)
pca <- fread(pca_file)

stopifnot(
  identical(colnames(expr), meta$SampleID),
  setequal(qc$SampleID, colnames(expr)),
  setequal(pca$SampleID, colnames(expr))
)

setorder(meta, MatrixOrder)

qc <- qc[
  match(colnames(expr), SampleID)
]

pca <- pca[
  match(colnames(expr), SampleID)
]


# ============================================================
# 1. Detailed phenotype context for flagged samples
# ============================================================

flagged_ids <- qc[
  AnyQCFlag == TRUE,
  SampleID
]

flag_meta <- merge(
  qc[
    SampleID %chin% flagged_ids
  ],
  meta[
    ,
    .(
      SampleID,
      Histopathology,
      DeepUlcer,
      ParisStage,
      Tissue,
      DiseaseRaw
    )
  ],
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

flag_meta <- flag_meta[
  match(flagged_ids, SampleID)
]


# ============================================================
# 2. Additional expression-composition metrics
# ============================================================

total <- colSums(expr)

max_expr <- apply(
  expr,
  2,
  max
)

top1_fraction <- max_expr / total

top10_fraction <- apply(
  expr,
  2,
  function(x) {

    total_x <- sum(x)

    if (!is.finite(total_x) || total_x <= 0) {
      return(NA_real_)
    }

    n <- length(x)

    top10 <- sort(
      x,
      partial = (n - 9):n
    )[
      (n - 9):n
    ]

    sum(top10) / total_x
  }
)

composition_qc <- data.table(
  SampleID = colnames(expr),
  Top1Fraction = top1_fraction,
  Top10Fraction = top10_fraction
)

flag_meta <- merge(
  flag_meta,
  composition_qc,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

flag_meta <- flag_meta[
  match(flagged_ids, SampleID)
]


# ============================================================
# 3. Exact duplicate-column audit
# ============================================================

# Exact duplicates must also share deterministic summaries.
signature <- data.table(
  SampleID = colnames(expr),

  Detected =
    colSums(expr > 0),

  Sum =
    colSums(expr),

  SumSq =
    colSums(expr * expr),

  Max =
    apply(expr, 2, max)
)

signature[
  ,
  SignatureKey := paste(
    Detected,
    sprintf("%.17g", Sum),
    sprintf("%.17g", SumSq),
    sprintf("%.17g", Max),
    sep = "|"
  )
]

candidate_groups <- signature[
  ,
  .N,
  by = SignatureKey
][
  N > 1
]

duplicate_pairs <- list()
k <- 0L

if (nrow(candidate_groups) > 0) {

  for (key in candidate_groups$SignatureKey) {

    ids <- signature[
      SignatureKey == key,
      SampleID
    ]

    if (length(ids) < 2)
      next

    cmb <- combn(
      ids,
      2,
      simplify = FALSE
    )

    for (pair in cmb) {

      same <- identical(
        as.numeric(expr[, pair[1]]),
        as.numeric(expr[, pair[2]])
      )

      if (same) {

        k <- k + 1L

        duplicate_pairs[[k]] <- data.table(
          Sample1 = pair[1],
          Sample2 = pair[2],
          ExactDuplicate = TRUE
        )
      }
    }
  }
}

exact_dup <- if (length(duplicate_pairs)) {
  rbindlist(duplicate_pairs)
} else {
  data.table(
    Sample1 = character(),
    Sample2 = character(),
    ExactDuplicate = logical()
  )
}


# ============================================================
# 4. Nearest-neighbour correlation audit
#
# Use top 5000 variable genes on log2(expr + 1)
# ============================================================

log_expr <- log2(expr + 1)

gene_var <- apply(
  log_expr,
  1,
  var
)

valid <- which(
  is.finite(gene_var) &
  gene_var > 0
)

top_n <- min(
  5000L,
  length(valid)
)

top_idx <- valid[
  order(
    gene_var[valid],
    decreasing = TRUE
  )[seq_len(top_n)]
]

cor_mat <- cor(
  log_expr[
    top_idx,
    ,
    drop = FALSE
  ],
  method = "pearson"
)

diag(cor_mat) <- NA_real_

nearest_index <- apply(
  cor_mat,
  2,
  which.max
)

nearest_cor <- vapply(
  seq_len(ncol(cor_mat)),
  function(i)
    cor_mat[
      nearest_index[i],
      i
    ],
  numeric(1)
)

second_cor <- vapply(
  seq_len(ncol(cor_mat)),
  function(i) {

    x <- cor_mat[, i]

    x <- sort(
      x[
        is.finite(x)
      ],
      decreasing = TRUE
    )

    if (length(x) >= 2)
      x[2]
    else
      NA_real_
  },
  numeric(1)
)

neighbor_qc <- data.table(
  SampleID = colnames(expr),
  NearestSample =
    colnames(expr)[nearest_index],
  NearestCorrelation =
    nearest_cor,
  SecondNearestCorrelation =
    second_cor
)


# ============================================================
# 5. Very-high-correlation pairs
# ============================================================

high_pairs <- which(
  upper.tri(cor_mat) &
  cor_mat >= 0.9999,
  arr.ind = TRUE
)

if (nrow(high_pairs) > 0) {

  high_corr <- data.table(
    Sample1 =
      rownames(cor_mat)[
        high_pairs[, 1]
      ],

    Sample2 =
      colnames(cor_mat)[
        high_pairs[, 2]
      ],

    Pearson =
      cor_mat[high_pairs]
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
# 6. Disease-stratified robust PCA distance
#
# Global PCA can flag biologically distinct disease groups.
# Reassess each sample relative to its own disease category.
# ============================================================

pc_cols <- paste0(
  "PC",
  1:5
)

pca2 <- copy(pca)

pca2[
  ,
  WithinDiseasePCDistance := NA_real_
]

for (grp in unique(pca2$Disease)) {

  idx <- which(
    pca2$Disease == grp
  )

  x <- as.matrix(
    pca2[
      idx,
      ..pc_cols
    ]
  )

  centers <- apply(
    x,
    2,
    median
  )

  scales <- apply(
    x,
    2,
    mad,
    constant = 1.4826
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

  pca2$WithinDiseasePCDistance[idx] <-
    sqrt(
      rowSums(z^2)
    )
}


pca2[
  ,
  WithinDiseaseThreshold :=
    median(
      WithinDiseasePCDistance
    ) +
    5 *
    mad(
      WithinDiseasePCDistance,
      constant = 1.4826
    ),
  by = Disease
]

pca2[
  ,
  WithinDiseasePCAFlag :=
    WithinDiseasePCDistance >
    WithinDiseaseThreshold
]


# ============================================================
# 7. Combine forensic information
# ============================================================

forensic <- merge(
  qc,
  neighbor_qc,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

forensic <- merge(
  forensic,
  composition_qc,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

forensic <- merge(
  forensic,
  pca2[
    ,
    .(
      SampleID,
      RobustPCDistance,
      PCAQCFlag,
      WithinDiseasePCDistance,
      WithinDiseaseThreshold,
      WithinDiseasePCAFlag
    )
  ],
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

forensic <- merge(
  forensic,
  meta[
    ,
    .(
      SampleID,
      Histopathology,
      DeepUlcer,
      ParisStage
    )
  ],
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

forensic <- forensic[
  match(
    colnames(expr),
    SampleID
  )
]


# ============================================================
# 8. 254/68 boundary descriptive audit
#
# Descriptive only — NOT treated as a batch variable.
# ============================================================

forensic[
  ,
  ProvenanceSegment :=
    fifelse(
      as.integer(
        sub(
          "CCFA_Risk_",
          "",
          SampleID
        )
      ) <= 254,
      "Risk001_254",
      "Risk255_322"
    )
]

segment_summary <- forensic[
  ,
  .(
    N = .N,

    MedianDetected =
      median(DetectedGenes),

    MedianZeroFraction =
      median(ZeroFraction),

    MedianLog2Total =
      median(Log2TotalExpression),

    AnyFlagN =
      sum(AnyQCFlag),

    WithinDiseasePCAFlagN =
      sum(WithinDiseasePCAFlag)
  ),
  by = ProvenanceSegment
]


# ============================================================
# 9. Console report
# ============================================================

cat("============================================================\n")
cat("EXACT DUPLICATES\n")
cat("============================================================\n")

cat(
  "Exact duplicate pairs:",
  nrow(exact_dup),
  "\n"
)

if (nrow(exact_dup) > 0)
  print(exact_dup)


cat("\n============================================================\n")
cat("VERY HIGH CORRELATION PAIRS (>=0.9999)\n")
cat("============================================================\n")

cat(
  "Pairs:",
  nrow(high_corr),
  "\n"
)

if (nrow(high_corr) > 0)
  print(
    head(
      high_corr,
      50
    )
  )


cat("\n============================================================\n")
cat("FLAGGED SAMPLE FORENSICS\n")
cat("============================================================\n")

print(
  forensic[
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

      UnivariateQCFlag,
      PCAQCFlag,
      WithinDiseasePCAFlag,

      RobustPCDistance,
      WithinDiseasePCDistance
    )
  ]
)


cat("\n============================================================\n")
cat("FLAG COUNTS\n")
cat("============================================================\n")

cat(
  "Original AnyQCFlag:",
  sum(forensic$AnyQCFlag),
  "\n"
)

cat(
  "Global PCA flags:",
  sum(forensic$PCAQCFlag),
  "\n"
)

cat(
  "Within-disease PCA flags:",
  sum(forensic$WithinDiseasePCAFlag),
  "\n"
)


cat("\n============================================================\n")
cat("FLAG DISTRIBUTION BY DISEASE\n")
cat("============================================================\n")

print(
  forensic[
    ,
    .(
      N = .N,
      OriginalFlagged =
        sum(AnyQCFlag),
      GlobalPCAFlagged =
        sum(PCAQCFlag),
      WithinDiseasePCAFlagged =
        sum(WithinDiseasePCAFlag)
    ),
    by = Disease
  ]
)


cat("\n============================================================\n")
cat("254/68 SEGMENT DESCRIPTIVE CHECK\n")
cat("============================================================\n")

print(segment_summary)


# ============================================================
# 10. Save
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
  segment_summary,
  file.path(
    out_dir,
    "GSE57945_254_68_segment_qc_summary.tsv"
  ),
  sep = "\t"
)


cat("\nIMPORTANT:\n")
cat(
  "- No sample was removed.\n",
  "- Nearest-neighbour correlations use top 5000 variable genes after log2(x+1).\n",
  "- Within-disease PCA flags are diagnostic and account partly for disease-level biology.\n",
  "- The 254/68 split is audited descriptively only and is NOT modeled as a batch.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_FLAGGED_SAMPLE_FORENSIC_QC_COMPLETED\n")
cat("============================================================\n")
