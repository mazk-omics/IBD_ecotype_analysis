suppressPackageStartupMessages({
  library(data.table)
})

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

qc_file <- file.path(
  base,
  "prepared/05_expression_qc",
  "GSE57945_flagged_sample_forensic_qc.tsv"
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

cat("============================================================\n")
cat("GSE57945 freeze sample-QC inclusion decision\n")
cat("============================================================\n\n")

qc <- fread(qc_file)
meta <- as.data.table(readRDS(meta_file))

stopifnot(
  nrow(qc) == 322,
  nrow(meta) == 322,
  !anyDuplicated(qc$SampleID),
  !anyDuplicated(meta$SampleID),
  setequal(qc$SampleID, meta$SampleID)
)

meta <- meta[
  match(qc$SampleID, SampleID)
]

stopifnot(
  identical(
    qc$SampleID,
    meta$SampleID
  )
)


# ============================================================
# Primary decision
#
# No sample currently meets sufficient evidence for a
# technical exclusion.
# ============================================================

manifest <- data.table(
  MatrixOrder = meta$MatrixOrder,
  SampleID = qc$SampleID,
  GSM = qc$GSM,
  Disease = qc$Disease,
  DiseaseSubtype = qc$DiseaseSubtype,

  IncludePrimary = TRUE,
  QCDecision = "KEEP",

  Stage13AnyQCFlag =
    qc$AnyQCFlag,

  UnivariateQCFlag =
    qc$UnivariateQCFlag,

  GlobalPCAQCFlag =
    qc$PCAQCFlag,

  WithinDiseasePCAFlag =
    qc$WithinDiseasePCAFlag,

  LowNearestCorrelationFlag =
    qc$LowNearestCorrelationFlag,

  NearestCorrelation =
    qc$NearestCorrelation,

  DetectedGenes =
    qc$DetectedGenes,

  ZeroFraction =
    qc$ZeroFraction,

  ProvenanceSegment =
    qc$ProvenanceSegment
)


# ------------------------------------------------------------
# Diagnostic severity label
#
# IMPORTANT:
# This is NOT an exclusion rule.
# It is only an audit annotation.
# ------------------------------------------------------------

manifest[
  ,
  QCAnnotation := fcase(

    Stage13AnyQCFlag &
      WithinDiseasePCAFlag &
      LowNearestCorrelationFlag,
    "MULTI_METRIC_FLAG",

    Stage13AnyQCFlag,
    "DIAGNOSTIC_FLAG",

    LowNearestCorrelationFlag |
      WithinDiseasePCAFlag,
    "SECONDARY_DIAGNOSTIC_FLAG",

    default =
      "NO_MAJOR_FLAG"
  )
]


# ============================================================
# Integrity checks
# ============================================================

stopifnot(
  nrow(manifest) == 322,
  all(manifest$IncludePrimary),
  all(manifest$QCDecision == "KEEP"),
  sum(manifest$Stage13AnyQCFlag) == 13
)

cat("Primary included samples:", sum(manifest$IncludePrimary), "\n")
cat("Primary excluded samples:", sum(!manifest$IncludePrimary), "\n")

cat("\nQC annotation counts:\n")

print(
  manifest[
    ,
    .N,
    by = QCAnnotation
  ][order(-N)]
)


cat("\nMulti-metric flagged samples:\n")

print(
  manifest[
    QCAnnotation == "MULTI_METRIC_FLAG",
    .(
      SampleID,
      GSM,
      Disease,
      NearestCorrelation,
      DetectedGenes,
      ZeroFraction
    )
  ]
)


# ============================================================
# Save
# ============================================================

fwrite(
  manifest,
  file.path(
    out_dir,
    "GSE57945_FINAL_sample_inclusion_manifest.tsv"
  ),
  sep = "\t"
)

saveRDS(
  manifest,
  file.path(
    out_dir,
    "GSE57945_FINAL_sample_inclusion_manifest.rds"
  )
)


decision_summary <- data.table(
  Metric = c(
    "Total_samples",
    "Primary_included",
    "Primary_excluded",
    "Stage13_diagnostic_flags",
    "Exact_duplicate_pairs",
    "Very_high_correlation_pairs"
  ),
  Value = c(
    322,
    322,
    0,
    13,
    0,
    0
  )
)

fwrite(
  decision_summary,
  file.path(
    out_dir,
    "GSE57945_FINAL_sample_qc_decision_summary.tsv"
  ),
  sep = "\t"
)


cat("\nDECISION:\n")
cat(
  "All 322 GSE57945 samples are retained for primary downstream analysis.\n"
)

cat(
  "QC flags are retained as diagnostic annotations only.\n"
)

cat(
  "No sample is excluded on the basis of current expression QC.\n"
)

cat(
  "The Risk001-254 / Risk255-322 provenance split is not used as a batch covariate.\n"
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE57945_SAMPLE_QC_DECISION_FROZEN\n")
cat("============================================================\n")
