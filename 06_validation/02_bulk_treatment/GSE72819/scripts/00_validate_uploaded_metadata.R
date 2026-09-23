suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

meta_file <- file.path(
  root,
  "01_source_metadata/original",
  "metadata_geo_GSE72819.xlsx"
)

out_dir <- file.path(
  root,
  "00_provenance"
)

cat("============================================================\n")
cat("GSE72819 uploaded metadata validation\n")
cat("============================================================\n\n")

stopifnot(file.exists(meta_file))

# ------------------------------------------------------------
# 1. Workbook sheets
# ------------------------------------------------------------

sheets <- excel_sheets(meta_file)

cat("Workbook sheets:\n")
print(sheets)

stopifnot(
  "GEO数据整理" %in% sheets
)

# ------------------------------------------------------------
# 2. Read authoritative metadata sheet
# ------------------------------------------------------------

x <- read_excel(
  meta_file,
  sheet = "GEO数据整理"
)

dt <- as.data.table(x)

cat("\nRows:", nrow(dt), "\n")
cat("Columns:", ncol(dt), "\n")

cat("\nColumn names:\n")
print(names(dt))

expected_cols <- c(
  "GSM编号",
  "疾病",
  "样本名",
  "组织诊断",
  "缓解情况w10",
  "先前抗tnf治疗",
  "tnf ir",
  "组织"
)

stopifnot(
  identical(
    names(dt),
    expected_cols
  )
)

# ------------------------------------------------------------
# 3. Core sample identity checks
# ------------------------------------------------------------

stopifnot(
  nrow(dt) == 73,
  uniqueN(dt$GSM编号) == 73,
  uniqueN(dt$样本名) == 73,
  !anyNA(dt$GSM编号),
  !anyNA(dt$样本名)
)

cat("\n[PASS] 73 rows\n")
cat("[PASS] 73 unique GSM\n")
cat("[PASS] 73 unique patient/sample IDs\n")

# ------------------------------------------------------------
# 4. Cohort structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DISEASE\n")
cat("============================================================\n")
print(dt[, .N, by = 疾病][order(-N)])

cat("\n============================================================\n")
cat("TISSUE\n")
cat("============================================================\n")
print(dt[, .N, by = 组织][order(-N)])

cat("\n============================================================\n")
cat("WEEK-10 REMISSION\n")
cat("============================================================\n")
print(dt[, .N, by = 缓解情况w10][order(-N)])

cat("\n============================================================\n")
cat("PRIOR ANTI-TNF\n")
cat("============================================================\n")
print(dt[, .N, by = 先前抗tnf治疗][order(-N)])

cat("\n============================================================\n")
cat("TNF IR\n")
cat("============================================================\n")
print(dt[, .N, by = `tnf ir`][order(-N)])

cat("\n============================================================\n")
cat("TISSUE DIAGNOSIS\n")
cat("============================================================\n")
print(dt[, .N, by = 组织诊断][order(-N)])

# ------------------------------------------------------------
# 5. Response eligibility
# ------------------------------------------------------------

dt[
  ,
  ResponseWeek10 :=
    fifelse(
      缓解情况w10 == "Remitter",
      "Response",
      fifelse(
        缓解情况w10 == "Non-remitter",
        "Non-response",
        NA_character_
      )
    )
]

cat("\n============================================================\n")
cat("RESPONSE ANALYSIS ELIGIBILITY\n")
cat("============================================================\n")

print(
  dt[
    ,
    .N,
    by = ResponseWeek10
  ][order(-N)]
)

cat(
  "\nEvaluable:",
  sum(!is.na(dt$ResponseWeek10)),
  "\n"
)

cat(
  "Not evaluable:",
  sum(is.na(dt$ResponseWeek10)),
  "\n"
)

stopifnot(
  sum(dt$ResponseWeek10 == "Response", na.rm = TRUE) == 12,
  sum(dt$ResponseWeek10 == "Non-response", na.rm = TRUE) == 58,
  sum(is.na(dt$ResponseWeek10)) == 3
)

# ------------------------------------------------------------
# 6. Missingness
# ------------------------------------------------------------

missing_summary <- data.table(
  Field = names(dt),
  Missing_N = vapply(
    dt,
    function(z) sum(is.na(z) | trimws(as.character(z)) == ""),
    integer(1)
  )
)

cat("\n============================================================\n")
cat("MISSINGNESS\n")
cat("============================================================\n")

print(missing_summary)

# ------------------------------------------------------------
# 7. Save immutable audit summaries only
# ------------------------------------------------------------

fwrite(
  missing_summary,
  file.path(
    out_dir,
    "GSE72819_uploaded_metadata_missingness.tsv"
  ),
  sep = "\t"
)

cohort_summary <- data.table(
  Metric = c(
    "Total_samples",
    "Unique_GSM",
    "Unique_patient_ids",
    "UC_samples",
    "Colonic_biopsy_samples",
    "Week10_response",
    "Week10_nonresponse",
    "Week10_missing",
    "Week10_evaluable"
  ),
  Value = c(
    nrow(dt),
    uniqueN(dt$GSM编号),
    uniqueN(dt$样本名),
    sum(dt$疾病 == "UC"),
    sum(dt$组织 == "colonic biopsy"),
    sum(dt$ResponseWeek10 == "Response", na.rm = TRUE),
    sum(dt$ResponseWeek10 == "Non-response", na.rm = TRUE),
    sum(is.na(dt$ResponseWeek10)),
    sum(!is.na(dt$ResponseWeek10))
  )
)

fwrite(
  cohort_summary,
  file.path(
    out_dir,
    "GSE72819_uploaded_metadata_cohort_summary.tsv"
  ),
  sep = "\t"
)

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] Workbook readable\n")
cat("[PASS] GEO数据整理 sheet present\n")
cat("[PASS] Metadata schema matches expected structure\n")
cat("[PASS] 73 unique GSM / 73 unique IDs\n")
cat("[PASS] Week-10 response structure: 12 / 58 / 3 missing\n")

cat("\nFINAL STATUS:\n")
cat("PASS_GSE72819_UPLOADED_METADATA_VALIDATED\n")
cat("============================================================\n")
