suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# ============================================================
# GSE57945
# Canonical sample metadata build
#
# Important:
#   - Raw metadata labels are preserved.
#   - Analysis labels are normalized separately.
#   - No phenotype values are inferred or imputed.
#   - Metadata rows are reordered exactly to expression columns.
# ============================================================

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

meta_sheet <- "RISK进展队列GEO回末测序数据整理"

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


# ============================================================
# 1. Load canonical expression matrix
# ============================================================

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


# ============================================================
# 2. Load the correct metadata sheet
# ============================================================

meta_raw <- as.data.table(
  read_excel(
    meta_file,
    sheet = meta_sheet
  )
)

setnames(
  meta_raw,
  trimws(names(meta_raw))
)

cat("Metadata sheet:", meta_sheet, "\n")
cat("Raw metadata rows:", nrow(meta_raw), "\n")
cat("Raw metadata columns:", ncol(meta_raw), "\n\n")


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

cat("[PASS] All required metadata columns detected\n")


# ============================================================
# 3. Keep valid sample rows
# ============================================================

meta_raw[
  ,
  `:=`(
    `GSM编号` = trimws(as.character(`GSM编号`)),
    `患者编号` = trimws(as.character(`患者编号`))
  )
]

meta_raw <- meta_raw[
  !is.na(`GSM编号`) &
    `GSM编号` != "" &
    !is.na(`患者编号`) &
    `患者编号` != ""
]

cat("Sample rows after filtering:", nrow(meta_raw), "\n")

stopifnot(
  nrow(meta_raw) == 322
)


# ============================================================
# 4. Construct provenance-preserving metadata
# ============================================================

meta <- meta_raw[
  ,
  .(
    GSM =
      trimws(as.character(`GSM编号`)),

    SampleID =
      trimws(as.character(`患者编号`)),

    TissueRaw =
      trimws(as.character(`取样部位`)),

    SexRaw =
      trimws(as.character(`性别`)),

    AgeAtDiagnosis =
      suppressWarnings(
        as.numeric(`年龄（诊断入组时）`)
      ),

    ParisStageRaw =
      trimws(as.character(`Paris分期`)),

    DiseaseRaw =
      trimws(as.character(`IBD疾病类型`)),

    DiseaseSubtypeRaw =
      trimws(as.character(`IBD疾病细分类型`)),

    HistopathologyRaw =
      trimws(as.character(`组织病理学特征`)),

    DeepUlcerRaw =
      trimws(as.character(`深部溃疡`))
  )
]


# ============================================================
# 5. Standardize literal missing values
# ============================================================

char_cols <- names(meta)[
  vapply(
    meta,
    is.character,
    logical(1)
  )
]

for (cc in char_cols) {

  idx <- which(
    is.na(meta[[cc]]) |
      meta[[cc]] %chin% c(
        "",
        "NA",
        "N/A",
        "na",
        "NaN",
        "NULL"
      )
  )

  if (length(idx) > 0) {
    set(
      meta,
      i = idx,
      j = cc,
      value = NA_character_
    )
  }
}


# ============================================================
# 6. Audit raw disease labels BEFORE normalization
# ============================================================

cat("\n============================================================\n")
cat("RAW DISEASE LABELS\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = DiseaseRaw
  ][order(-N)]
)

cat("\nRaw disease subtype labels:\n")

print(
  meta[
    ,
    .N,
    by = DiseaseSubtypeRaw
  ][order(-N)]
)


# Normalize only for matching.
# Raw columns remain unchanged.
disease_key <- tolower(
  trimws(meta$DiseaseRaw)
)

subtype_key <- tolower(
  trimws(meta$DiseaseSubtypeRaw)
)


# Hard-check that no unexpected disease label exists.
allowed_disease_keys <- c(
  "cd",
  "uc",
  "not ibd"
)

unexpected_disease <- unique(
  disease_key[
    !is.na(disease_key) &
      !disease_key %chin% allowed_disease_keys
  ]
)

if (length(unexpected_disease) > 0) {
  stop(
    "Unexpected raw disease labels: ",
    paste(unexpected_disease, collapse = ", ")
  )
}


# ============================================================
# 7. Create normalized analysis variables
# ============================================================

# ------------------------------------------------------------
# Disease
# Case-insensitive handling:
# "Not IBD" and "not IBD" both -> "Non-IBD"
# ------------------------------------------------------------

meta[
  ,
  Disease := fcase(
    disease_key == "cd", "CD",
    disease_key == "uc", "UC",
    disease_key == "not ibd", "Non-IBD",
    default = NA_character_
  )
]


# ------------------------------------------------------------
# Disease subtype
#
# Disease status takes precedence for UC / Non-IBD.
# CD subtype is normalized case-insensitively.
# ------------------------------------------------------------

meta[
  ,
  DiseaseSubtype := fcase(
    Disease == "Non-IBD", "Non-IBD",
    Disease == "UC", "UC",
    Disease == "CD" & subtype_key == "icd", "iCD",
    Disease == "CD" & subtype_key == "ccd", "cCD",
    default = DiseaseSubtypeRaw
  )
]


# ------------------------------------------------------------
# Sex
# ------------------------------------------------------------

sex_key <- tolower(
  trimws(meta$SexRaw)
)

meta[
  ,
  Sex := fcase(
    SexRaw == "男", "Male",
    SexRaw == "女", "Female",
    sex_key %chin% c("male", "m"), "Male",
    sex_key %chin% c("female", "f"), "Female",
    default = NA_character_
  )
]


# ------------------------------------------------------------
# Tissue
# ------------------------------------------------------------

meta[
  ,
  Tissue := fcase(
    grepl(
      "ileum|回肠",
      TissueRaw,
      ignore.case = TRUE
    ),
    "Ileum",
    default = TissueRaw
  )
]


# ------------------------------------------------------------
# Histopathology
# Preserve source wording
# ------------------------------------------------------------

meta[
  ,
  Histopathology := HistopathologyRaw
]


# ------------------------------------------------------------
# Deep ulcer
# ------------------------------------------------------------

ulcer_key <- tolower(
  trimws(meta$DeepUlcerRaw)
)

meta[
  ,
  DeepUlcer := fcase(
    ulcer_key %chin% c("1", "yes"), "Yes",
    ulcer_key %chin% c("0", "no"), "No",
    default = NA_character_
  )
]


# ------------------------------------------------------------
# Paris stage
# ------------------------------------------------------------

meta[
  ,
  ParisStage := ParisStageRaw
]


# ============================================================
# 8. Identity QC
# ============================================================

stopifnot(
  nrow(meta) == 322,
  uniqueN(meta$GSM) == 322,
  uniqueN(meta$SampleID) == 322
)

if (!setequal(sample_ids, meta$SampleID)) {

  cat("\nSamples in expression but not metadata:\n")
  print(
    setdiff(
      sample_ids,
      meta$SampleID
    )
  )

  cat("\nSamples in metadata but not expression:\n")
  print(
    setdiff(
      meta$SampleID,
      sample_ids
    )
  )

  stop(
    "Expression/metadata SampleID mismatch."
  )
}

cat(
  "[PASS] Expression and metadata contain identical 322 SampleIDs\n"
)


# ============================================================
# 9. Reorder metadata exactly to expression columns
# ============================================================

idx <- match(
  sample_ids,
  meta$SampleID
)

stopifnot(
  !anyNA(idx)
)

meta <- meta[idx]

meta[
  ,
  MatrixOrder := seq_len(.N)
]

stopifnot(
  identical(
    meta$SampleID,
    colnames(expr)
  )
)

cat(
  "[PASS] Metadata order exactly matches expression columns\n\n"
)


# ============================================================
# 10. Reorder output columns
# ============================================================

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
    "ParisStage",

    "DiseaseRaw",
    "DiseaseSubtypeRaw",
    "TissueRaw",
    "HistopathologyRaw",
    "DeepUlcerRaw",
    "SexRaw",
    "ParisStageRaw"
  )
)


# ============================================================
# 11. Cohort phenotype QC
# ============================================================

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
  ][
    order(
      Disease,
      -N
    )
  ]
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
  ][
    order(
      Disease,
      -N
    )
  ]
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
  ][
    order(
      Disease,
      -N
    )
  ]
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

      Missing =
        sum(
          is.na(AgeAtDiagnosis)
        ),

      Mean =
        mean(
          AgeAtDiagnosis,
          na.rm = TRUE
        ),

      SD =
        sd(
          AgeAtDiagnosis,
          na.rm = TRUE
        ),

      Median =
        median(
          AgeAtDiagnosis,
          na.rm = TRUE
        ),

      Min =
        min(
          AgeAtDiagnosis,
          na.rm = TRUE
        ),

      Max =
        max(
          AgeAtDiagnosis,
          na.rm = TRUE
        )
    )
  ]
)


# ============================================================
# 12. Hard cohort checks
# ============================================================

count_value <- function(x, value) {
  sum(
    !is.na(x) &
      x == value
  )
}

# ------------------------------------------------------------
# Disease composition
# ------------------------------------------------------------

stopifnot(
  count_value(meta$Disease, "CD") == 218,
  count_value(meta$Disease, "UC") == 62,
  count_value(meta$Disease, "Non-IBD") == 42,
  sum(is.na(meta$Disease)) == 0
)

cat(
  "\n[PASS] Disease composition = CD 218 / UC 62 / Non-IBD 42\n"
)


# ------------------------------------------------------------
# Disease subtype
# ------------------------------------------------------------

stopifnot(
  count_value(meta$DiseaseSubtype, "iCD") == 163,
  count_value(meta$DiseaseSubtype, "cCD") == 55,
  count_value(meta$DiseaseSubtype, "UC") == 62,
  count_value(meta$DiseaseSubtype, "Non-IBD") == 42
)

cat(
  "[PASS] Disease subtype composition matches expected cohort\n"
)


# ------------------------------------------------------------
# All tissue = ileum
# ------------------------------------------------------------

stopifnot(
  count_value(meta$Tissue, "Ileum") == 322
)

cat(
  "[PASS] All 322 samples normalized to Ileum\n"
)


# ------------------------------------------------------------
# CD-specific histopathology
# ------------------------------------------------------------

cd <- meta[
  Disease == "CD"
]

stopifnot(
  nrow(cd) == 218,

  count_value(
    cd$Histopathology,
    "Macroscopic inflammation"
  ) == 162,

  count_value(
    cd$Histopathology,
    "Microscopic inflammation"
  ) == 29,

  count_value(
    cd$Histopathology,
    "Normal"
  ) == 24,

  count_value(
    cd$Histopathology,
    "Undetermined"
  ) == 3
)

cat(
  "[PASS] CD histopathology = 162 / 29 / 24 / 3\n"
)


# ------------------------------------------------------------
# CD deep ulcer
# ------------------------------------------------------------

stopifnot(
  count_value(cd$DeepUlcer, "Yes") == 76,
  count_value(cd$DeepUlcer, "No") == 142,
  sum(is.na(cd$DeepUlcer)) == 0
)

cat(
  "[PASS] CD deep ulcer = Yes 76 / No 142\n"
)


# ------------------------------------------------------------
# Sex
# ------------------------------------------------------------

stopifnot(
  count_value(meta$Sex, "Male") == 187,
  count_value(meta$Sex, "Female") == 135,
  sum(is.na(meta$Sex)) == 0
)

cat(
  "[PASS] Sex = Male 187 / Female 135\n"
)


# ------------------------------------------------------------
# Age
# ------------------------------------------------------------

stopifnot(
  sum(is.na(meta$AgeAtDiagnosis)) == 0
)

cat(
  "[PASS] Age available for all 322 samples\n"
)


# ============================================================
# 13. Explicit audit of raw Not IBD capitalization variant
# ============================================================

cat("\n============================================================\n")
cat("NON-IBD RAW-LABEL NORMALIZATION AUDIT\n")
cat("============================================================\n")

non_ibd_audit <- meta[
  Disease == "Non-IBD",
  .(
    GSM,
    SampleID,
    DiseaseRaw,
    DiseaseSubtypeRaw,
    Disease,
    DiseaseSubtype
  )
]

print(
  non_ibd_audit[
    ,
    .N,
    by = .(
      DiseaseRaw,
      DiseaseSubtypeRaw,
      Disease,
      DiseaseSubtype
    )
  ]
)

# Confirm the known lower-case record is preserved in Raw fields
stopifnot(
  non_ibd_audit[
    DiseaseRaw == "not IBD",
    .N
  ] == 1
)

stopifnot(
  non_ibd_audit[
    DiseaseRaw == "Not IBD",
    .N
  ] == 41
)

cat(
  "[PASS] Raw Not IBD capitalization preserved: 41 'Not IBD' + 1 'not IBD'\n"
)


# ============================================================
# 14. Missingness audit
# ============================================================

missingness <- data.table(
  Variable = names(meta),

  MissingN = vapply(
    meta,
    function(x)
      sum(is.na(x)),
    integer(1)
  )
)

missingness[
  ,
  MissingFraction :=
    MissingN / nrow(meta)
]

setorder(
  missingness,
  -MissingN,
  Variable
)

cat("\n============================================================\n")
cat("MISSINGNESS\n")
cat("============================================================\n")

print(missingness)


# ============================================================
# 15. Save outputs
# ============================================================

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


fwrite(
  missingness,
  file.path(
    out_dir,
    "GSE57945_sample_metadata_missingness.tsv"
  ),
  sep = "\t"
)


fwrite(
  non_ibd_audit,
  file.path(
    out_dir,
    "GSE57945_nonIBD_label_normalization_audit.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 16. Final validation
# ============================================================

cat("\n============================================================\n")
cat("FINAL SAMPLE METADATA\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(meta),
  "\n"
)

cat(
  "Unique SampleIDs:",
  uniqueN(meta$SampleID),
  "\n"
)

cat(
  "Unique GSMs:",
  uniqueN(meta$GSM),
  "\n"
)

cat(
  "Exact matrix-column alignment:",
  identical(
    meta$SampleID,
    colnames(expr)
  ),
  "\n"
)


cat("\nIMPORTANT:\n")
cat(
  "- Raw GEO/RISK phenotype labels are preserved in *Raw columns.\n",
  "- 'Not IBD' and 'not IBD' are normalized only in analysis variables.\n",
  "- No phenotype value was inferred or imputed.\n",
  "- Metadata order exactly matches expression columns.\n",
  "- No progression endpoint was inferred for GSE57945.\n",
  sep = ""
)


cat("\nFINAL STATUS:\n")
cat("PASS_CANONICAL_SAMPLE_METADATA_BUILT\n")
cat("============================================================\n")
