suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

expr_file <- file.path(
  base,
  "prepared/03_canonical_expression",
  "GSE57945_canonical_expression_34368x322.rds"
)

meta_file <- file.path(
  base,
  "metadata",
  "RISK_metadata_from_GEO.xlsx"
)

out_dir <- file.path(
  base,
  "prepared/04_sample_metadata"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE57945 canonical sample metadata build\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Load canonical expression matrix
# ------------------------------------------------------------

expr <- readRDS(expr_file)

stopifnot(
  nrow(expr) == 34368,
  ncol(expr) == 322,
  !anyDuplicated(colnames(expr))
)

sample_ids <- colnames(expr)

cat("Expression samples:", length(sample_ids), "\n")
cat("First sample:", sample_ids[1], "\n")
cat("Last sample :", sample_ids[length(sample_ids)], "\n\n")


# ------------------------------------------------------------
# 2. Load source metadata
# ------------------------------------------------------------

meta_raw <- as.data.table(
  read_excel(meta_file)
)

cat("Metadata rows:", nrow(meta_raw), "\n")
cat("Metadata columns:", ncol(meta_raw), "\n\n")

required_cols <- c(
  "GSM编号",
  "患者编号",
  "取样部位",
  "性别",
  "年龄（诊断入组时）",
  "Paris分期",
  "IBD疾病类型",
  "IBD疾病细分类型",
  "组织病理学特征",
  "深部溃疡"
)

missing_cols <- setdiff(
  required_cols,
  names(meta_raw)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

stopifnot(nrow(meta_raw) == 322)


# ------------------------------------------------------------
# 3. Construct canonical metadata
# ------------------------------------------------------------

meta <- meta_raw[
  ,
  .(
    GSM = trimws(as.character(`GSM编号`)),
    SampleID = trimws(as.character(`患者编号`)),
    Tissue = trimws(as.character(`取样部位`)),
    Sex = trimws(as.character(`性别`)),
    AgeAtDiagnosis = as.numeric(`年龄（诊断入组时）`),
    ParisStage = trimws(as.character(`Paris分期`)),
    Disease = trimws(as.character(`IBD疾病类型`)),
    DiseaseSubtype = trimws(as.character(`IBD疾病细分类型`)),
    Histopathology = trimws(as.character(`组织病理学特征`)),
    DeepUlcer = trimws(as.character(`深部溃疡`))
  )
]

# clean literal NA-like strings
char_cols <- names(meta)[
  vapply(meta, is.character, logical(1))
]

for (cc in char_cols) {
  set(
    meta,
    which(meta[[cc]] %chin% c("", "NA", "N/A", "na")),
    cc,
    NA_character_
  )
}


# ------------------------------------------------------------
# 4. Identity QC
# ------------------------------------------------------------

stopifnot(
  uniqueN(meta$GSM) == 322,
  uniqueN(meta$SampleID) == 322
)

if (!setequal(sample_ids, meta$SampleID)) {

  cat("Samples in expression but not metadata:\n")
  print(setdiff(sample_ids, meta$SampleID))

  cat("\nSamples in metadata but not expression:\n")
  print(setdiff(meta$SampleID, sample_ids))

  stop("Expression/metadata SampleID mismatch.")
}

cat("[PASS] Expression and metadata contain identical 322 SampleIDs\n")


# ------------------------------------------------------------
# 5. Reorder metadata to expression matrix
# ------------------------------------------------------------

idx <- match(
  sample_ids,
  meta$SampleID
)

stopifnot(!anyNA(idx))

meta <- meta[idx]

stopifnot(
  identical(
    meta$SampleID,
    colnames(expr)
  )
)

meta[
  ,
  MatrixOrder := seq_len(.N)
]

setcolorder(
  meta,
  c(
    "MatrixOrder",
    "SampleID",
    "GSM",
    "Disease",
    "DiseaseSubtype",
    "Tissue",
    "Histopathology",
    "DeepUlcer",
    "Sex",
    "AgeAtDiagnosis",
    "ParisStage"
  )
)

cat("[PASS] Metadata reordered exactly to expression columns\n\n")


# ------------------------------------------------------------
# 6. Basic phenotype QC
# ------------------------------------------------------------

cat("============================================================\n")
cat("DISEASE\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = Disease
  ][order(-N)]
)

cat("\n============================================================\n")
cat("DISEASE SUBTYPE\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = .(
      Disease,
      DiseaseSubtype
    )
  ][order(Disease, -N)]
)

cat("\n============================================================\n")
cat("TISSUE\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = Tissue
  ]
)

cat("\n============================================================\n")
cat("HISTOPATHOLOGY\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = .(
      Disease,
      Histopathology
    )
  ][order(Disease, -N)]
)

cat("\n============================================================\n")
cat("DEEP ULCER\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = .(
      Disease,
      DeepUlcer
    )
  ][order(Disease, -N)]
)

cat("\n============================================================\n")
cat("SEX\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = Sex
  ][order(-N)]
)

cat("\n============================================================\n")
cat("AGE\n")
cat("============================================================\n")

print(
  meta[
    ,
    .(
      N = .N,
      Missing = sum(is.na(AgeAtDiagnosis)),
      Mean = mean(AgeAtDiagnosis, na.rm = TRUE),
      SD = sd(AgeAtDiagnosis, na.rm = TRUE),
      Median = median(AgeAtDiagnosis, na.rm = TRUE),
      Min = min(AgeAtDiagnosis, na.rm = TRUE),
      Max = max(AgeAtDiagnosis, na.rm = TRUE)
    )
  ]
)


# ------------------------------------------------------------
# 7. Known cohort consistency checks
# ------------------------------------------------------------

disease_counts <- meta[, .N, by = Disease]

get_n <- function(group) {
  x <- disease_counts[
    Disease == group,
    N
  ]
  if (length(x) == 0) 0L else x
}

stopifnot(
  get_n("CD") == 218,
  get_n("UC") == 62,
  get_n("Non-IBD") == 42
)

cat("\n[PASS] Disease composition matches expected 322-sample RISK cohort\n")


# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------

fwrite(
  meta,
  file.path(
    out_dir,
    "GSE57945_canonical_sample_metadata_322.tsv"
  ),
  sep = "\t",
  na = "NA"
)

saveRDS(
  meta,
  file.path(
    out_dir,
    "GSE57945_canonical_sample_metadata_322.rds"
  )
)

# Expression ↔ metadata identity manifest
identity <- meta[
  ,
  .(
    MatrixOrder,
    SampleID,
    GSM
  )
]

fwrite(
  identity,
  file.path(
    out_dir,
    "GSE57945_expression_metadata_identity_manifest.tsv"
  ),
  sep = "\t"
)


cat("\n============================================================\n")
cat("FINAL SAMPLE METADATA\n")
cat("============================================================\n")

cat("Rows:", nrow(meta), "\n")
cat("Unique SampleIDs:", uniqueN(meta$SampleID), "\n")
cat("Unique GSMs:", uniqueN(meta$GSM), "\n")

cat(
  "Exact matrix-column alignment:",
  identical(meta$SampleID, colnames(expr)),
  "\n"
)

cat("\nIMPORTANT:\n")
cat(
  "- No phenotype values were inferred or imputed.\n",
  "- Metadata order exactly matches canonical expression columns.\n",
  "- Original GEO/RISK phenotype labels were preserved.\n",
  "- No progression endpoint was created for GSE57945.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_CANONICAL_SAMPLE_METADATA_BUILT\n")
cat("============================================================\n")
