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
# 2. Load the correct metadata worksheet
# ============================================================

meta_raw <- as.data.table(
  read_excel(
    meta_file,
    sheet = meta_sheet
  )
)

# trim column names
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
# 3. Remove completely blank / non-sample rows if any
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
# 4. Build metadata table
#
# Preserve original labels AND create normalized variables.
# ============================================================

meta <- meta_raw[
  ,
  .(
    GSM = trimws(as.character(`GSM编号`)),
    SampleID = trimws(as.character(`患者编号`)),

    TissueRaw = trimws(as.character(`取样部位`)),
    SexRaw = trimws(as.character(`性别`)),

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


# ------------------------------------------------------------
# Normalize literal missing values
# ------------------------------------------------------------

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
      "NaN"
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
# 5. Create normalized analysis variables
# ============================================================

meta[
  ,
  Disease := fcase(
    DiseaseRaw == "CD", "CD",
    DiseaseRaw == "UC", "UC",
    DiseaseRaw == "Not IBD", "Non-IBD",
    default = NA_character_
  )
]

meta[
  ,
  DiseaseSubtype := fcase(
    DiseaseSubtypeRaw == "iCD", "iCD",
    DiseaseSubtypeRaw == "cCD", "cCD",
    DiseaseRaw == "UC", "UC",
    DiseaseRaw == "Not IBD", "Non-IBD",
    default = DiseaseSubtypeRaw
  )
]

meta[
  ,
  Sex := fcase(
    SexRaw == "男", "Male",
    SexRaw == "女", "Female",
    SexRaw %chin% c("Male", "M"), "Male",
    SexRaw %chin% c("Female", "F"), "Female",
    default = NA_character_
  )
]

meta[
  ,
  Tissue := fcase(
    grepl(
      "Ileum|回肠",
      TissueRaw,
      ignore.case = TRUE
    ),
    "Ileum",
    default = TissueRaw
  )
]

meta[
  ,
  Histopathology := HistopathologyRaw
]

meta[
  ,
  DeepUlcer := fcase(
    DeepUlcerRaw == "1", "Yes",
    DeepUlcerRaw == "0", "No",
    DeepUlcerRaw %chin% c("Yes", "YES", "yes"), "Yes",
    DeepUlcerRaw %chin% c("No", "NO", "no"), "No",
    default = NA_character_
  )
]

meta[
  ,
  ParisStage := ParisStageRaw
]


# ============================================================
# 6. Identity QC
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
# 7. Reorder metadata exactly to expression columns
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
# 8. Reorder columns
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
# 9. Cohort composition audit
# ============================================================

cat("============================================================\n")
cat("DISEASE\n")
cat("============================================================\n")

disease_tab <- meta[
  ,
  .N,
  by = Disease
][order(-N)]

print(disease_tab)


cat("\n============================================================\n")
cat("DISEASE SUBTYPE\n")
cat("============================================================\n")

subtype_tab <- meta[
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

print(subtype_tab)


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

hist_tab <- meta[
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

print(hist_tab)


cat("\n============================================================\n")
cat("DEEP ULCER\n")
cat("============================================================\n")

ulcer_tab <- meta[
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

print(ulcer_tab)


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
      Missing = sum(
        is.na(AgeAtDiagnosis)
      ),
      Mean = mean(
        AgeAtDiagnosis,
        na.rm = TRUE
      ),
      SD = sd(
        AgeAtDiagnosis,
        na.rm = TRUE
      ),
      Median = median(
        AgeAtDiagnosis,
        na.rm = TRUE
      ),
      Min = min(
        AgeAtDiagnosis,
        na.rm = TRUE
      ),
      Max = max(
        AgeAtDiagnosis,
        na.rm = TRUE
      )
    )
  ]
)


# ============================================================
# 10. Hard cohort consistency checks
# ============================================================

get_n <- function(dt, col, value) {

  x <- dt[
    get(col) == value,
    .N
  ]

  if (length(x) == 0)
    return(0L)

  x
}

stopifnot(
  get_n(meta, "Disease", "CD") == 218,
  get_n(meta, "Disease", "UC") == 62,
  get_n(meta, "Disease", "Non-IBD") == 42
)

stopifnot(
  get_n(meta, "DiseaseSubtype", "iCD") == 163,
  get_n(meta, "DiseaseSubtype", "cCD") == 55
)

cat(
  "\n[PASS] Disease and subtype composition match expected RISK cohort\n"
)


# ------------------------------------------------------------
# CD histopathology checks
# ------------------------------------------------------------

cd <- meta[
  Disease == "CD"
]

stopifnot(
  cd[
    Histopathology ==
      "Macroscopic inflammation",
    .N
  ] == 162,

  cd[
    Histopathology ==
      "Microscopic inflammation",
    .N
  ] == 29,

  cd[
    Histopathology ==
      "Normal",
    .N
  ] == 24
)

undetermined_n <- cd[
  is.na(Histopathology) |
  Histopathology %chin% c(
    "Undetermined",
    "undetermined"
  ),
  .N
]

cat(
  "CD histopathology unresolved/undetermined:",
  undetermined_n,
  "\n"
)

stopifnot(
  undetermined_n == 3
)


# ------------------------------------------------------------
# CD deep ulcer checks
# ------------------------------------------------------------

stopifnot(
  cd[
    DeepUlcer == "Yes",
    .N
  ] == 76,

  cd[
    DeepUlcer == "No",
    .N
  ] == 142
)

cat(
  "[PASS] CD histopathology and deep-ulcer counts match expected cohort\n"
)


# ============================================================
# 11. Missingness audit
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
# 12. Save
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


# ============================================================
# 13. Final validation
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
  "- Raw GEO/RISK phenotype labels were preserved in *Raw columns.\n",
  "- Normalized analysis variables were created separately.\n",
  "- No phenotype value was imputed.\n",
  "- Metadata order exactly matches expression columns.\n",
  "- No progression endpoint was inferred for GSE57945.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_CANONICAL_SAMPLE_METADATA_BUILT\n")
cat("============================================================\n")
