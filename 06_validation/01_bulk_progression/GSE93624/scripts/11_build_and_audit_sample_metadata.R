suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 4
#
# Build canonical sample metadata from GEO SOFT.
#
# Rules:
#   - preserve all raw GEO fields
#   - normalize only after raw cross-tab is verified
#   - progression "-" is NOT globally "No"
#   - Non-IBD progression is structural NA
#   - sample order must exactly match canonical expression
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

soft_file <- file.path(
  base,
  "metadata",
  "GSE93624_family.soft.gz"
)

expr_file <- file.path(
  base,
  "prepared/03_canonical_expression",
  "GSE93624_canonical_expression_13323x245.rds"
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
cat("GSE93624 sample metadata build + audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load canonical expression
# ============================================================

obj <- readRDS(
  expr_file
)

expr <- obj$expression

stopifnot(
  nrow(expr) == 13323,
  ncol(expr) == 245
)

expr_samples <- colnames(expr)

stopifnot(
  length(expr_samples) == 245,
  !anyDuplicated(expr_samples)
)

cat(
  "[PASS] Canonical expression:",
  nrow(expr),
  "x",
  ncol(expr),
  "\n"
)

cat(
  "Expression sample labels, first 10:\n"
)

print(
  head(
    expr_samples,
    10
  )
)


# ============================================================
# 2. Read GEO SOFT
# ============================================================

soft <- readLines(
  soft_file,
  warn = FALSE,
  encoding = "UTF-8"
)

cat(
  "[PASS] SOFT lines loaded:",
  length(soft),
  "\n"
)


# ============================================================
# 3. Identify SAMPLE blocks
# ============================================================

sample_start <- grep(
  "^\\^SAMPLE = ",
  soft
)

stopifnot(
  length(sample_start) == 245
)

sample_end <- c(
  sample_start[-1] - 1L,
  length(soft)
)

cat(
  "[PASS] GEO SAMPLE blocks:",
  length(sample_start),
  "\n"
)


# ============================================================
# 4. Parse one sample block
# ============================================================

extract_one <- function(lines) {

  gsm_line <- grep(
    "^\\^SAMPLE = ",
    lines,
    value = TRUE
  )

  title_line <- grep(
    "^!Sample_title = ",
    lines,
    value = TRUE
  )

  source_line <- grep(
    "^!Sample_source_name_ch1 = ",
    lines,
    value = TRUE
  )

  characteristics <- grep(
    "^!Sample_characteristics_ch1 = ",
    lines,
    value = TRUE
  )

  gsm <- sub(
    "^\\^SAMPLE = ",
    "",
    gsm_line[1]
  )

  title <- sub(
    "^!Sample_title = ",
    "",
    title_line[1]
  )

  source <- if (
    length(source_line) > 0
  ) {
    sub(
      "^!Sample_source_name_ch1 = ",
      "",
      source_line[1]
    )
  } else {
    NA_character_
  }


  # ----------------------------------------------------------
  # Characteristics parser:
  #   "diagnosis: Crohn's disease"
  #   "gender: Male"
  # etc.
  # ----------------------------------------------------------

  ch <- sub(
    "^!Sample_characteristics_ch1 = ",
    "",
    characteristics
  )

  keys <- trimws(
    sub(
      ":.*$",
      "",
      ch
    )
  )

  vals <- trimws(
    sub(
      "^[^:]+:",
      "",
      ch
    )
  )

  kv <- setNames(
    vals,
    keys
  )


  get_chr <- function(key) {

    if (key %in% names(kv)) {
      unname(kv[[key]])
    } else {
      NA_character_
    }
  }


  data.table(
    GSM =
      gsm,

    SampleTitle =
      title,

    SourceName_raw =
      source,

    Tissue_raw =
      get_chr("tissue"),

    Progression_raw =
      get_chr("progression to complication"),

    ParisAge_raw =
      get_chr("paris age"),

    Gender_raw =
      get_chr("gender"),

    Diagnosis_raw =
      get_chr("diagnosis"),

    Ancestry_raw =
      get_chr("ancestry"),

    AgeAtDiagnosis_raw =
      get_chr("age at diagnosis")
  )
}


# ============================================================
# 5. Parse all 245 samples
# ============================================================

meta_list <- vector(
  "list",
  245
)

for (i in seq_len(245)) {

  meta_list[[i]] <- extract_one(
    soft[
      sample_start[i]:
      sample_end[i]
    ]
  )
}

meta <- rbindlist(
  meta_list,
  use.names = TRUE,
  fill = TRUE
)


stopifnot(
  nrow(meta) == 245,
  uniqueN(meta$GSM) == 245,
  uniqueN(meta$SampleTitle) == 245
)

cat(
  "[PASS] Parsed metadata rows:",
  nrow(meta),
  "\n"
)


# ============================================================
# 6. Raw field completeness
# ============================================================

cat("\n============================================================\n")
cat("RAW FIELD COMPLETENESS\n")
cat("============================================================\n")

raw_fields <- c(
  "GSM",
  "SampleTitle",
  "Tissue_raw",
  "Progression_raw",
  "ParisAge_raw",
  "Gender_raw",
  "Diagnosis_raw",
  "Ancestry_raw",
  "AgeAtDiagnosis_raw"
)

completeness <- rbindlist(
  lapply(
    raw_fields,
    function(v) {

      x <- meta[[v]]

      data.table(
        Field =
          v,

        NonMissing =
          sum(
            !is.na(x) &
            x != ""
          ),

        Missing =
          sum(
            is.na(x) |
            x == ""
          ),

        UniqueValues =
          uniqueN(
            x[
              !is.na(x) &
              x != ""
            ]
          )
      )
    }
  )
)

print(
  completeness
)


# ============================================================
# 7. Raw distributions
# ============================================================

cat("\n============================================================\n")
cat("RAW DIAGNOSIS\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = Diagnosis_raw
  ][order(-N)]
)


cat("\n============================================================\n")
cat("RAW PROGRESSION\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = Progression_raw
  ][order(-N)]
)


cat("\n============================================================\n")
cat("RAW DIAGNOSIS x PROGRESSION\n")
cat("============================================================\n")

dx_prog <- meta[
  ,
  .N,
  by = .(
    Diagnosis_raw,
    Progression_raw
  )
][
  order(
    Diagnosis_raw,
    Progression_raw
  )
]

print(
  dx_prog
)


cat("\n============================================================\n")
cat("RAW TISSUE\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = Tissue_raw
  ][order(-N)]
)


cat("\n============================================================\n")
cat("RAW GENDER\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = Gender_raw
  ][order(-N)]
)


cat("\n============================================================\n")
cat("RAW PARIS AGE\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = ParisAge_raw
  ][order(-N)]
)


cat("\n============================================================\n")
cat("RAW ANCESTRY\n")
cat("============================================================\n")

print(
  meta[
    ,
    .N,
    by = Ancestry_raw
  ][order(-N)]
)


# ============================================================
# 8. Hard provenance expectations
# ============================================================

stopifnot(
  meta[
    Diagnosis_raw == "Crohn's disease",
    .N
  ] == 210,

  meta[
    Diagnosis_raw == "Non-IBD control",
    .N
  ] == 35,

  meta[
    Progression_raw == "yes",
    .N
  ] == 27,

  meta[
    Progression_raw == "-",
    .N
  ] == 218
)

cat(
  "\n[PASS] Raw diagnosis counts = CD 210 / Non-IBD 35\n"
)

cat(
  "[PASS] Raw progression counts = yes 27 / '-' 218\n"
)


# ============================================================
# 9. Critical diagnosis × progression validation
# ============================================================

n_cd_yes <- meta[
  Diagnosis_raw == "Crohn's disease" &
    Progression_raw == "yes",
  .N
]

n_cd_dash <- meta[
  Diagnosis_raw == "Crohn's disease" &
    Progression_raw == "-",
  .N
]

n_ctrl_yes <- meta[
  Diagnosis_raw == "Non-IBD control" &
    Progression_raw == "yes",
  .N
]

n_ctrl_dash <- meta[
  Diagnosis_raw == "Non-IBD control" &
    Progression_raw == "-",
  .N
]


cat("\n============================================================\n")
cat("CRITICAL PROGRESSION STRUCTURE AUDIT\n")
cat("============================================================\n")

cat(
  "CD + yes :",
  n_cd_yes,
  "\n"
)

cat(
  "CD + '-' :",
  n_cd_dash,
  "\n"
)

cat(
  "Non-IBD + yes :",
  n_ctrl_yes,
  "\n"
)

cat(
  "Non-IBD + '-' :",
  n_ctrl_dash,
  "\n"
)


stopifnot(
  n_cd_yes == 27,
  n_cd_dash == 183,
  n_ctrl_yes == 0,
  n_ctrl_dash == 35
)

cat(
  "[PASS] Progression field is structurally interpretable within CD\n"
)

cat(
  "[PASS] Non-IBD '-' is structural N/A, not progression='No'\n"
)


# ============================================================
# 10. Normalize metadata
# ============================================================

meta[
  ,
  Diagnosis :=
    fcase(

      Diagnosis_raw ==
        "Crohn's disease",
      "CD",

      Diagnosis_raw ==
        "Non-IBD control",
      "Non-IBD",

      default =
        NA_character_
    )
]


meta[
  ,
  Progression3yr :=
    fcase(

      Diagnosis ==
        "CD" &
      Progression_raw ==
        "yes",
      "Yes",

      Diagnosis ==
        "CD" &
      Progression_raw ==
        "-",
      "No",

      Diagnosis ==
        "Non-IBD",
      NA_character_,

      default =
        NA_character_
    )
]


# ------------------------------------------------------------
# Sex normalization:
# preserve raw value, normalize common forms only.
# ------------------------------------------------------------

meta[
  ,
  Sex :=
    fcase(

      tolower(Gender_raw) %chin%
        c(
          "male",
          "m"
        ),
      "Male",

      tolower(Gender_raw) %chin%
        c(
          "female",
          "f"
        ),
      "Female",

      default =
        NA_character_
    )
]


# ------------------------------------------------------------
# Tissue:
# preserve GEO wording while providing broad normalized tissue.
# ------------------------------------------------------------

meta[
  ,
  Tissue :=
    fifelse(
      grepl(
        "ile",
        Tissue_raw,
        ignore.case = TRUE
      ),
      "Ileum",
      Tissue_raw
    )
]


# ------------------------------------------------------------
# Age numeric parsing.
#
# Do NOT discard raw value if parsing fails.
# ------------------------------------------------------------

meta[
  ,
  AgeAtDiagnosis :=
    suppressWarnings(
      as.numeric(
        AgeAtDiagnosis_raw
      )
    )
]


# ============================================================
# 11. Normalized metadata summary
# ============================================================

cat("\n============================================================\n")
cat("NORMALIZED METADATA SUMMARY\n")
cat("============================================================\n")

cat("\nDiagnosis:\n")

print(
  meta[
    ,
    .N,
    by = Diagnosis
  ]
)


cat("\nProgression3yr:\n")

print(
  meta[
    ,
    .N,
    by = .(
      Diagnosis,
      Progression3yr
    )
  ]
)


cat("\nSex:\n")

print(
  meta[
    ,
    .N,
    by = Sex
  ]
)


cat("\nTissue:\n")

print(
  meta[
    ,
    .N,
    by = Tissue
  ]
)


# ============================================================
# 12. Expression-column alignment
# ============================================================

# ------------------------------------------------------------
# First determine whether expression uses SampleTitle or GSM.
# ------------------------------------------------------------

title_match <- setequal(
  expr_samples,
  meta$SampleTitle
)

gsm_match <- setequal(
  expr_samples,
  meta$GSM
)


cat("\n============================================================\n")
cat("EXPRESSION / METADATA ALIGNMENT\n")
cat("============================================================\n")

cat(
  "Expression columns match SampleTitle set:",
  title_match,
  "\n"
)

cat(
  "Expression columns match GSM set:",
  gsm_match,
  "\n"
)


if (title_match) {

  meta[
    ,
    ExpressionSampleID :=
      SampleTitle
  ]

} else if (gsm_match) {

  meta[
    ,
    ExpressionSampleID :=
      GSM
  ]

} else {

  stop(
    paste0(
      "Expression columns match neither SampleTitle nor GSM. ",
      "Inspect Stage 1 column labels before continuing."
    )
  )
}


# ------------------------------------------------------------
# Reorder metadata to exact expression-column order
# ------------------------------------------------------------

idx <- match(
  expr_samples,
  meta$ExpressionSampleID
)

stopifnot(
  !anyNA(idx)
)

meta_aligned <- meta[
  idx
]

stopifnot(
  identical(
    meta_aligned$ExpressionSampleID,
    expr_samples
  )
)

cat(
  "[PASS] Metadata reordered to exact expression-column order\n"
)


# ============================================================
# 13. Sample-number audit
# ============================================================

meta_aligned[
  ,
  SampleNumber :=
    suppressWarnings(
      as.integer(
        sub(
          "^Sample",
          "",
          SampleTitle
        )
      )
    )
]


stopifnot(
  all(
    !is.na(
      meta_aligned$SampleNumber
    )
  ),
  setequal(
    meta_aligned$SampleNumber,
    1:245
  )
)

cat(
  "[PASS] SampleTitle maps exactly to Sample1-Sample245\n"
)


# ============================================================
# 14. Final canonical metadata ordering
# ============================================================

canonical_meta <- meta_aligned[
  ,
  .(
    ExpressionSampleID,
    GSM,
    SampleTitle,
    SampleNumber,

    Diagnosis,
    Progression3yr,
    Sex,
    Tissue,
    AgeAtDiagnosis,
    ParisAge =
      ParisAge_raw,
    Ancestry =
      Ancestry_raw,

    Diagnosis_raw,
    Progression_raw,
    Gender_raw,
    Tissue_raw,
    ParisAge_raw,
    Ancestry_raw,
    AgeAtDiagnosis_raw,
    SourceName_raw
  )
]


# ============================================================
# 15. Final hard checks
# ============================================================

stopifnot(
  nrow(canonical_meta) == 245,
  uniqueN(
    canonical_meta$ExpressionSampleID
  ) == 245,
  uniqueN(
    canonical_meta$GSM
  ) == 245,
  uniqueN(
    canonical_meta$SampleTitle
  ) == 245,

  canonical_meta[
    Diagnosis == "CD",
    .N
  ] == 210,

  canonical_meta[
    Diagnosis == "Non-IBD",
    .N
  ] == 35,

  canonical_meta[
    Progression3yr == "Yes",
    .N
  ] == 27,

  canonical_meta[
    Progression3yr == "No",
    .N
  ] == 183,

  canonical_meta[
    Diagnosis == "Non-IBD" &
      is.na(Progression3yr),
    .N
  ] == 35,

  identical(
    canonical_meta$ExpressionSampleID,
    colnames(expr)
  )
)


# ============================================================
# 16. Save
# ============================================================

fwrite(
  canonical_meta,
  file.path(
    out_dir,
    "GSE93624_canonical_sample_metadata.tsv"
  ),
  sep = "\t",
  na = ""
)

saveRDS(
  canonical_meta,
  file.path(
    out_dir,
    "GSE93624_canonical_sample_metadata.rds"
  ),
  compress = FALSE
)

fwrite(
  dx_prog,
  file.path(
    out_dir,
    "GSE93624_raw_diagnosis_progression_crosstab.tsv"
  ),
  sep = "\t"
)

fwrite(
  completeness,
  file.path(
    out_dir,
    "GSE93624_metadata_field_completeness.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 17. Freeze audit
# ============================================================

freeze <- data.table(
  Metric = c(
    "Samples",
    "CD",
    "NonIBD",
    "CD_ProgressionYes",
    "CD_ProgressionNo",
    "NonIBD_ProgressionNA",
    "ExpressionMetadataOrderExact"
  ),

  Value = c(
    "245",
    "210",
    "35",
    "27",
    "183",
    "35",
    "TRUE"
  )
)

fwrite(
  freeze,
  file.path(
    out_dir,
    "GSE93624_metadata_freeze_summary.tsv"
  ),
  sep = "\t"
)


cat("\n============================================================\n")
cat("FINAL METADATA INTEGRITY CHECKS\n")
cat("============================================================\n")

cat("[PASS] Samples = 245\n")
cat("[PASS] CD = 210\n")
cat("[PASS] Non-IBD = 35\n")
cat("[PASS] CD progression Yes = 27\n")
cat("[PASS] CD progression No = 183\n")
cat("[PASS] Non-IBD progression = structural NA for all 35 controls\n")
cat("[PASS] Metadata order exactly matches canonical expression columns\n")

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_CANONICAL_SAMPLE_METADATA_BUILT\n")
cat("============================================================\n")
