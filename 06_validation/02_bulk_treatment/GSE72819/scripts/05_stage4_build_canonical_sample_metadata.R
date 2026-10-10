suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# ============================================================
# GSE72819 Stage 4
# Build canonical sample metadata
#
# Frozen inputs:
#   Stage 3 canonical expression:
#       24,656 genes x 73 samples
#
#   Stage 1 source sample manifest:
#       GSM + technical provenance + normalized clinical fields
#
#   Stage 0 four-way GSM crosswalk:
#       Local metadata / SOFT / Series Matrix / RAW
#
#   Original uploaded Excel:
#       independent source-level clinical metadata verification
#
# Important:
#   - All 73 expression samples are retained.
#   - 3 samples with missing Week-10 outcome remain eligible for
#     expression/recovery, but are NOT eligible for response analysis.
#   - No expression matrix is modified.
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

canonical_file <- file.path(
  root,
  "prepared/03_canonical_expression",
  "GSE72819_canonical_RPKM_24656x73.rds"
)

manifest_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE72819_source_sample_manifest.tsv"
)

crosswalk_file <- file.path(
  root,
  "00_provenance",
  "GSE72819_four_way_GSM_crosswalk.tsv"
)

excel_file <- file.path(
  root,
  "01_source_metadata/original",
  "metadata_geo_GSE72819.xlsx"
)

out_dir <- file.path(
  root,
  "prepared/04_sample_metadata"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Helper functions
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


clean_chr <- function(x) {

  x <- trimws(as.character(x))

  x[
    is.na(x) |
    x %in% c(
      "",
      "NA",
      "N/A",
      "na",
      "n/a",
      "-"
    )
  ] <- NA_character_

  x
}


same_chr <- function(x, y) {

  x <- clean_chr(x)
  y <- clean_chr(y)

  if (length(x) != length(y)) {
    return(FALSE)
  }

  all(
    (is.na(x) & is.na(y)) |
      (!is.na(x) & !is.na(y) & x == y)
  )
}


cat("============================================================\n")
cat("GSE72819 Stage 4: canonical sample metadata\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load canonical expression and freeze sample order
# ============================================================

expr <- readRDS(
  canonical_file
)

check(
  is.matrix(expr),
  "Canonical expression is a matrix"
)

check(
  nrow(expr) == 24656,
  "Canonical expression contains 24,656 genes"
)

check(
  ncol(expr) == 73,
  "Canonical expression contains 73 samples"
)

sample_order <- as.character(
  colnames(expr)
)

check(
  length(unique(sample_order)) == 73,
  "All 73 canonical expression sample IDs are unique"
)


# ============================================================
# 2. Load source sample manifest
# ============================================================

manifest <- fread(
  manifest_file
)

required_manifest <- c(
  "GSM",
  "Run",
  "Lane",
  "Library",
  "SAM",
  "NXG",
  "PatientID",
  "Disease",
  "Tissue",
  "TissueDiagnosis",
  "ResponseWeek10",
  "PriorAntiTNF",
  "TNF_IR"
)

check(
  all(required_manifest %in% names(manifest)),
  "Source sample manifest contains all required fields"
)

check(
  nrow(manifest) == 73,
  "Source sample manifest contains 73 rows"
)

manifest[
  ,
  (required_manifest) :=
    lapply(
      .SD,
      clean_chr
    ),
  .SDcols = required_manifest
]

check(
  uniqueN(manifest$GSM) == 73,
  "Source manifest contains 73 unique GSM IDs"
)

check(
  uniqueN(manifest$PatientID) == 73,
  "Source manifest contains 73 unique patient IDs"
)


# ============================================================
# 3. Load four-way provenance crosswalk
# ============================================================

crosswalk <- fread(
  crosswalk_file
)

required_crosswalk <- c(
  "GSM",
  "In_LocalMetadata",
  "In_SOFT",
  "In_SeriesMatrix",
  "In_RAW",
  "AllFour"
)

check(
  all(required_crosswalk %in% names(crosswalk)),
  "Four-way GSM crosswalk has required columns"
)

check(
  nrow(crosswalk) == 73,
  "Four-way GSM crosswalk contains 73 rows"
)

check(
  uniqueN(crosswalk$GSM) == 73,
  "Four-way GSM crosswalk contains 73 unique GSM IDs"
)

check(
  all(crosswalk$In_LocalMetadata),
  "All GSM IDs occur in uploaded metadata"
)

check(
  all(crosswalk$In_SOFT),
  "All GSM IDs occur in GEO SOFT"
)

check(
  all(crosswalk$In_SeriesMatrix),
  "All GSM IDs occur in GEO Series Matrix"
)

check(
  all(crosswalk$In_RAW),
  "All GSM IDs occur in GEO RAW processed files"
)

check(
  all(crosswalk$AllFour),
  "All 73 GSM IDs pass four-way provenance concordance"
)


# ============================================================
# 4. Load original Excel metadata
# ============================================================

excel_raw <- as.data.table(
  read_excel(
    excel_file,
    sheet = "GEO数据整理"
  )
)

required_excel <- c(
  "GSM编号",
  "疾病",
  "样本名",
  "组织诊断",
  "缓解情况w10",
  "先前抗tnf治疗",
  "tnf ir",
  "组织"
)

check(
  all(required_excel %in% names(excel_raw)),
  "Original Excel contains all required clinical fields"
)

check(
  nrow(excel_raw) == 73,
  "Original Excel GEO metadata contains 73 rows"
)


excel <- data.table(
  GSM = clean_chr(
    excel_raw[["GSM编号"]]
  ),

  Disease = clean_chr(
    excel_raw[["疾病"]]
  ),

  PatientID = clean_chr(
    excel_raw[["样本名"]]
  ),

  TissueDiagnosis = clean_chr(
    excel_raw[["组织诊断"]]
  ),

  RawResponseWeek10 = clean_chr(
    excel_raw[["缓解情况w10"]]
  ),

  PriorAntiTNF = clean_chr(
    excel_raw[["先前抗tnf治疗"]]
  ),

  TNF_IR = clean_chr(
    excel_raw[["tnf ir"]]
  ),

  Tissue = clean_chr(
    excel_raw[["组织"]]
  )
)


check(
  uniqueN(excel$GSM) == 73,
  "Original Excel contains 73 unique GSM IDs"
)

check(
  uniqueN(excel$PatientID) == 73,
  "Original Excel contains 73 unique patient IDs"
)


# ============================================================
# 5. Normalize Week-10 response from original Excel
# ============================================================

excel[
  ,
  ResponseWeek10 :=
    fcase(

      RawResponseWeek10 %in%
        c(
          "Remitter",
          "Response"
        ),
      "Response",

      RawResponseWeek10 %in%
        c(
          "Non-remitter",
          "Non-response"
        ),
      "Non-response",

      is.na(RawResponseWeek10),
      NA_character_,

      default =
        paste0(
          "REVIEW:",
          RawResponseWeek10
        )
    )
]


check(
  !any(
    grepl(
      "^REVIEW:",
      excel$ResponseWeek10
    ),
    na.rm = TRUE
  ),
  "All original Week-10 response labels normalize unambiguously"
)


# ============================================================
# 6. Verify identical GSM universe across all resources
# ============================================================

gsm_expr <- sort(sample_order)
gsm_manifest <- sort(manifest$GSM)
gsm_crosswalk <- sort(crosswalk$GSM)
gsm_excel <- sort(excel$GSM)

check(
  identical(gsm_expr, gsm_manifest),
  "Expression and source manifest contain identical GSM universe"
)

check(
  identical(gsm_expr, gsm_crosswalk),
  "Expression and provenance crosswalk contain identical GSM universe"
)

check(
  identical(gsm_expr, gsm_excel),
  "Expression and original Excel contain identical GSM universe"
)


# ============================================================
# 7. Align source resources to canonical expression order
# ============================================================

manifest <- manifest[
  match(
    sample_order,
    GSM
  )
]

excel <- excel[
  match(
    sample_order,
    GSM
  )
]

crosswalk <- crosswalk[
  match(
    sample_order,
    GSM
  )
]


check(
  identical(
    manifest$GSM,
    sample_order
  ),
  "Source manifest aligned exactly to canonical expression order"
)

check(
  identical(
    excel$GSM,
    sample_order
  ),
  "Original Excel aligned exactly to canonical expression order"
)

check(
  identical(
    crosswalk$GSM,
    sample_order
  ),
  "Four-way crosswalk aligned exactly to canonical expression order"
)


# ============================================================
# 8. Independent Excel vs source-manifest clinical audit
# ============================================================

clinical_concordance <- data.table(
  Field = c(
    "PatientID",
    "Disease",
    "Tissue",
    "TissueDiagnosis",
    "ResponseWeek10",
    "PriorAntiTNF",
    "TNF_IR"
  ),

  Concordant = c(
    same_chr(
      manifest$PatientID,
      excel$PatientID
    ),

    same_chr(
      manifest$Disease,
      excel$Disease
    ),

    same_chr(
      manifest$Tissue,
      excel$Tissue
    ),

    same_chr(
      manifest$TissueDiagnosis,
      excel$TissueDiagnosis
    ),

    same_chr(
      manifest$ResponseWeek10,
      excel$ResponseWeek10
    ),

    same_chr(
      manifest$PriorAntiTNF,
      excel$PriorAntiTNF
    ),

    same_chr(
      manifest$TNF_IR,
      excel$TNF_IR
    )
  )
)


cat("\n============================================================\n")
cat("EXCEL vs SOURCE MANIFEST CONCORDANCE\n")
cat("============================================================\n")

print(clinical_concordance)


check(
  all(clinical_concordance$Concordant),
  "Original Excel and source manifest clinical fields are fully concordant"
)


# ============================================================
# 9. Build canonical metadata
# ============================================================

meta <- data.table(
  SampleID = sample_order,

  GSM = manifest$GSM,

  PatientID = manifest$PatientID,

  Disease = manifest$Disease,

  Tissue = manifest$Tissue,

  TissueDiagnosis =
    manifest$TissueDiagnosis,

  Timepoint =
    rep(
      "Baseline",
      73
    ),

  Treatment =
    rep(
      "Etrolizumab",
      73
    ),

  DoseArm =
    rep(
      NA_character_,
      73
    ),

  RawResponseWeek10 =
    excel$RawResponseWeek10,

  ResponseWeek10 =
    manifest$ResponseWeek10,

  ResponseEvaluable =
    !is.na(
      manifest$ResponseWeek10
    ),

  PriorAntiTNF =
    manifest$PriorAntiTNF,

  TNF_IR =
    manifest$TNF_IR,

  Run =
    manifest$Run,

  Lane =
    manifest$Lane,

  Library =
    manifest$Library,

  SAM =
    manifest$SAM,

  NXG =
    manifest$NXG,

  FourWayProvenancePASS =
    crosswalk$AllFour,

  ExpressionRetained =
    rep(
      TRUE,
      73
    ),

  RecoveryEligible =
    rep(
      TRUE,
      73
    ),

  ResponseAnalysisEligible =
    !is.na(
      manifest$ResponseWeek10
    )
)


# ============================================================
# 10. Core identity validation
# ============================================================

check(
  nrow(meta) == 73,
  "Canonical metadata contains exactly 73 rows"
)

check(
  uniqueN(meta$SampleID) == 73,
  "Canonical SampleID values are unique"
)

check(
  uniqueN(meta$GSM) == 73,
  "Canonical GSM values are unique"
)

check(
  uniqueN(meta$PatientID) == 73,
  "Canonical PatientID values are unique"
)

check(
  identical(
    meta$SampleID,
    sample_order
  ),
  "Canonical metadata order exactly matches expression columns"
)

check(
  identical(
    meta$GSM,
    sample_order
  ),
  "GSM order exactly matches expression columns"
)


# ============================================================
# 11. Frozen cohort checks
# ============================================================

check(
  all(meta$Disease == "UC"),
  "All 73 samples are UC"
)

check(
  all(meta$Tissue == "colonic biopsy"),
  "All 73 samples are colonic biopsy"
)

check(
  all(meta$Timepoint == "Baseline"),
  "All 73 samples are baseline"
)

check(
  all(meta$Treatment == "Etrolizumab"),
  "All 73 samples are assigned to etrolizumab validation cohort"
)

check(
  all(meta$FourWayProvenancePASS),
  "All 73 samples retain four-way provenance PASS"
)

check(
  all(meta$ExpressionRetained),
  "All 73 samples are retained for expression analysis"
)

check(
  all(meta$RecoveryEligible),
  "All 73 samples are retained for EcoTyper recovery"
)


# ============================================================
# 12. Outcome checks
# ============================================================

response_counts <- meta[
  ,
  .N,
  by = ResponseWeek10
][order(ResponseWeek10)]

cat("\n============================================================\n")
cat("WEEK-10 RESPONSE SUMMARY\n")
cat("============================================================\n")

print(response_counts)


check(
  meta[
    ResponseWeek10 == "Response",
    .N
  ] == 12,
  "Week-10 Response = 12"
)

check(
  meta[
    ResponseWeek10 == "Non-response",
    .N
  ] == 58,
  "Week-10 Non-response = 58"
)

check(
  sum(is.na(meta$ResponseWeek10)) == 3,
  "Week-10 outcome missing = 3"
)

check(
  sum(meta$ResponseEvaluable) == 70,
  "Week-10 response-evaluable samples = 70"
)

check(
  sum(meta$ResponseAnalysisEligible) == 70,
  "Response-analysis eligible samples = 70"
)

check(
  sum(
    !meta$ResponseAnalysisEligible
  ) == 3,
  "Exactly 3 samples are excluded only from response association"
)


# ============================================================
# 13. Prior anti-TNF checks
# ============================================================

prior_tnf_counts <- meta[
  ,
  .N,
  by = PriorAntiTNF
][order(PriorAntiTNF)]

cat("\n============================================================\n")
cat("PRIOR ANTI-TNF SUMMARY\n")
cat("============================================================\n")

print(prior_tnf_counts)


check(
  meta[
    PriorAntiTNF == "Prior anti-TNF",
    .N
  ] == 46,
  "Prior anti-TNF = 46"
)

check(
  meta[
    PriorAntiTNF == "Naive",
    .N
  ] == 27,
  "Anti-TNF naive = 27"
)


# ============================================================
# 14. Technical provenance checks
# ============================================================

run_counts <- meta[
  ,
  .N,
  by = Run
][order(Run)]

lane_counts <- meta[
  ,
  .N,
  by = Lane
][order(Lane)]


cat("\n============================================================\n")
cat("RUN SUMMARY\n")
cat("============================================================\n")

print(run_counts)

cat("\n============================================================\n")
cat("LANE SUMMARY\n")
cat("============================================================\n")

print(lane_counts)


check(
  meta[
    Run == "R1510",
    .N
  ] == 31,
  "Run R1510 = 31"
)

check(
  meta[
    Run == "R1518",
    .N
  ] == 42,
  "Run R1518 = 42"
)


expected_lane <- data.table(
  Lane = c(
    "L1",
    "L2",
    "L3",
    "L4",
    "L5",
    "L6"
  ),

  ExpectedN = c(
    19,
    19,
    11,
    11,
    6,
    7
  )
)

observed_lane <- meta[
  ,
  .N,
  by = Lane
]

lane_check <- merge(
  expected_lane,
  observed_lane,
  by = "Lane",
  all = TRUE
)

check(
  all(
    lane_check$ExpectedN ==
      lane_check$N
  ),
  "Lane counts reproduce frozen Stage-0 provenance"
)


# ============================================================
# 15. Run × outcome table
# ============================================================

run_response <- meta[
  !is.na(ResponseWeek10),
  .N,
  by = .(
    Run,
    ResponseWeek10
  )
][order(Run, ResponseWeek10)]


cat("\n============================================================\n")
cat("RUN x WEEK-10 RESPONSE\n")
cat("============================================================\n")

print(run_response)


# ============================================================
# 16. Run × prior anti-TNF table
# ============================================================

run_prior_tnf <- meta[
  ,
  .N,
  by = .(
    Run,
    PriorAntiTNF
  )
][order(Run, PriorAntiTNF)]


cat("\n============================================================\n")
cat("RUN x PRIOR ANTI-TNF\n")
cat("============================================================\n")

print(run_prior_tnf)


# ============================================================
# 17. Missingness audit
# ============================================================

missingness <- data.table(
  Field = names(meta),

  Missing_N = vapply(
    meta,
    function(x) {
      sum(
        is.na(x) |
        (
          is.character(x) &
          trimws(x) == ""
        )
      )
    },
    integer(1)
  )
)


cat("\n============================================================\n")
cat("CANONICAL METADATA MISSINGNESS\n")
cat("============================================================\n")

print(
  missingness[
    Missing_N > 0
  ]
)


# DoseArm is intentionally unavailable.
check(
  missingness[
    Field == "DoseArm",
    Missing_N
  ] == 73,
  "DoseArm is explicitly unavailable for all 73 public samples"
)

check(
  missingness[
    Field == "ResponseWeek10",
    Missing_N
  ] == 3,
  "Only 3 Week-10 outcomes are missing"
)


# ============================================================
# 18. Save canonical metadata
# ============================================================

metadata_tsv <- file.path(
  out_dir,
  "GSE72819_canonical_sample_metadata_73.tsv"
)

metadata_rds <- file.path(
  out_dir,
  "GSE72819_canonical_sample_metadata_73.rds"
)

fwrite(
  meta,
  metadata_tsv,
  sep = "\t",
  quote = FALSE,
  na = ""
)

saveRDS(
  meta,
  metadata_rds,
  compress = FALSE
)


# ============================================================
# 19. Save audit tables
# ============================================================

fwrite(
  clinical_concordance,
  file.path(
    out_dir,
    "GSE72819_excel_vs_manifest_concordance.tsv"
  ),
  sep = "\t"
)

fwrite(
  response_counts,
  file.path(
    out_dir,
    "GSE72819_week10_response_summary.tsv"
  ),
  sep = "\t",
  na = "NA"
)

fwrite(
  prior_tnf_counts,
  file.path(
    out_dir,
    "GSE72819_prior_antiTNF_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  run_counts,
  file.path(
    out_dir,
    "GSE72819_run_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  lane_counts,
  file.path(
    out_dir,
    "GSE72819_lane_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  run_response,
  file.path(
    out_dir,
    "GSE72819_run_by_response.tsv"
  ),
  sep = "\t"
)

fwrite(
  run_prior_tnf,
  file.path(
    out_dir,
    "GSE72819_run_by_prior_antiTNF.tsv"
  ),
  sep = "\t"
)

fwrite(
  missingness,
  file.path(
    out_dir,
    "GSE72819_canonical_metadata_missingness.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 20. Round-trip validation
# ============================================================

reload <- readRDS(
  metadata_rds
)

check(
  nrow(reload) == 73,
  "Canonical metadata RDS round-trip retains 73 rows"
)

check(
  identical(
    as.character(reload$GSM),
    sample_order
  ),
  "Canonical metadata RDS sample order round-trip is exact"
)


# ============================================================
# 21. Artifact hashes
# ============================================================

md5 <- tools::md5sum(
  c(
    canonical_file,
    manifest_file,
    crosswalk_file,
    excel_file,
    metadata_tsv,
    metadata_rds
  )
)

artifact_md5 <- data.table(
  File = names(md5),
  MD5 = unname(md5)
)

fwrite(
  artifact_md5,
  file.path(
    out_dir,
    "GSE72819_stage4_artifact_md5.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 22. Final summary
# ============================================================

summary <- data.table(
  Metric = c(
    "Total_samples",
    "Unique_GSM",
    "Unique_patients",
    "UC_samples",
    "Colonic_biopsy_samples",
    "Baseline_samples",
    "Etrolizumab_samples",
    "Week10_response",
    "Week10_nonresponse",
    "Week10_missing",
    "Week10_evaluable",
    "Prior_antiTNF",
    "AntiTNF_naive",
    "Run_R1510",
    "Run_R1518",
    "Expression_retained",
    "Recovery_eligible",
    "Response_analysis_eligible",
    "Response_analysis_ineligible"
  ),

  N = c(
    nrow(meta),
    uniqueN(meta$GSM),
    uniqueN(meta$PatientID),
    meta[Disease == "UC", .N],
    meta[Tissue == "colonic biopsy", .N],
    meta[Timepoint == "Baseline", .N],
    meta[Treatment == "Etrolizumab", .N],
    meta[ResponseWeek10 == "Response", .N],
    meta[ResponseWeek10 == "Non-response", .N],
    sum(is.na(meta$ResponseWeek10)),
    sum(meta$ResponseEvaluable),
    meta[PriorAntiTNF == "Prior anti-TNF", .N],
    meta[PriorAntiTNF == "Naive", .N],
    meta[Run == "R1510", .N],
    meta[Run == "R1518", .N],
    sum(meta$ExpressionRetained),
    sum(meta$RecoveryEligible),
    sum(meta$ResponseAnalysisEligible),
    sum(!meta$ResponseAnalysisEligible)
  )
)


cat("\n============================================================\n")
cat("FINAL SAMPLE METADATA SUMMARY\n")
cat("============================================================\n")

print(summary)


fwrite(
  summary,
  file.path(
    out_dir,
    "GSE72819_stage4_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 23. Final gate
# ============================================================

critical_pass <- all(c(
  nrow(meta) == 73,
  uniqueN(meta$GSM) == 73,
  uniqueN(meta$PatientID) == 73,

  identical(
    meta$GSM,
    sample_order
  ),

  all(meta$FourWayProvenancePASS),

  meta[
    ResponseWeek10 == "Response",
    .N
  ] == 12,

  meta[
    ResponseWeek10 == "Non-response",
    .N
  ] == 58,

  sum(
    is.na(meta$ResponseWeek10)
  ) == 3,

  sum(
    meta$ResponseAnalysisEligible
  ) == 70,

  sum(
    meta$RecoveryEligible
  ) == 73,

  meta[
    PriorAntiTNF == "Prior anti-TNF",
    .N
  ] == 46,

  meta[
    PriorAntiTNF == "Naive",
    .N
  ] == 27,

  meta[
    Run == "R1510",
    .N
  ] == 31,

  meta[
    Run == "R1518",
    .N
  ] == 42,

  all(
    clinical_concordance$Concordant
  )
))


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] 73 canonical metadata rows retained\n")
  cat("[PASS] 73 unique GSM IDs\n")
  cat("[PASS] 73 unique patient IDs\n")
  cat("[PASS] Metadata order exactly matches canonical expression\n")
  cat("[PASS] Original Excel agrees with source manifest\n")
  cat("[PASS] Four-way GEO provenance retained for all samples\n")
  cat("[PASS] Week-10 outcome: 12 response / 58 non-response / 3 missing\n")
  cat("[PASS] 70 samples eligible for response association\n")
  cat("[PASS] All 73 samples remain eligible for EcoTyper recovery\n")
  cat("[PASS] Prior anti-TNF: 46 prior / 27 naive\n")
  cat("[PASS] Run structure retained without batch correction\n")
  cat("[PASS] No sample excluded from expression or recovery\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE4_CANONICAL_SAMPLE_METADATA_BUILT\n")

} else {

  cat("[REVIEW] Stage 4 sample metadata requires review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE4_CANONICAL_SAMPLE_METADATA_BUILT\n")
}

cat("============================================================\n")
