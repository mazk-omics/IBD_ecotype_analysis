suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 0B v3
# Robust longitudinal cohort reconstruction
#
# Primary metadata:
#   GEO Sample characteristics
#
# Secondary/fallback metadata:
#   GEO Sample title
#
# Controls are explicitly assigned:
#   Timepoint = CONTROL
#
# No expression processing occurs here.
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

soft_file <- file.path(
  root,
  "00_provenance/geo_downloads",
  "GSE73661_family.soft.gz"
)

matrix_file <- file.path(
  root,
  "00_provenance/geo_downloads",
  "GSE73661_series_matrix.txt.gz"
)

out_dir <- file.path(
  root,
  "prepared/00_cohort_reconstruction"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Helpers
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
        "-",
        "null",
        "NULL"
      )
  ] <- NA_character_

  x
}


soft_value <- function(x) {

  x <- sub(
    "^[^=]+=[[:space:]]*",
    "",
    x
  )

  x <- trimws(x)

  x <- sub(
    '^"',
    "",
    x
  )

  x <- sub(
    '"$',
    "",
    x
  )

  x
}


normalize_key <- function(x) {

  x <- tolower(
    trimws(
      as.character(x)
    )
  )

  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )

  x
}


extract_characteristic_regex <- function(char_vec, patterns) {

  if (length(char_vec) == 0L) {
    return(NA_character_)
  }

  for (z in char_vec) {

    z <- clean_chr(z)

    if (is.na(z)) {
      next
    }

    pos <- regexpr(
      ":",
      z,
      fixed = TRUE
    )

    if (pos[1] < 1L) {
      next
    }

    key <- substr(
      z,
      1L,
      pos[1] - 1L
    )

    value <- substr(
      z,
      pos[1] + 1L,
      nchar(z)
    )

    key_norm <- normalize_key(key)

    matched <- any(
      vapply(
        patterns,
        function(p) {
          grepl(
            p,
            key_norm,
            perl = TRUE,
            ignore.case = TRUE
          )
        },
        logical(1)
      )
    )

    if (matched) {
      return(
        clean_chr(value)
      )
    }
  }

  NA_character_
}


normalize_week_value <- function(x) {

  x <- clean_chr(x)

  if (is.na(x)) {
    return(NA_character_)
  }

  y <- toupper(x)

  y <- gsub(
    "[[:space:]]+",
    "",
    y
  )

  y <- gsub(
    "\\u2013|\\u2014",
    "-",
    y
  )

  # Explicitly handle interval-style IFX labels such as:
  # W4_W6 / W4-W6 / W4/W6 / 4_6 / 4-6 / 4/6
  if (grepl(
    "^W?4[_/-]W?6$",
    y,
    perl = TRUE
  )) {
    return("W4_6")
  }

  # Remove common textual week prefixes.
  y <- gsub(
    "^WEEKS?",
    "",
    y,
    perl = TRUE
  )

  y <- gsub(
    "^WKS?",
    "",
    y,
    perl = TRUE
  )

  y <- gsub(
    "^WK",
    "",
    y,
    perl = TRUE
  )

  # Standard single-week values:
  # 0 / W0 / 6 / W6 / 12 / W12 / 52 / W52 ...
  if (grepl(
    "^W?[0-9]+$",
    y,
    perl = TRUE
  )) {

    n <- sub(
      "^W",
      "",
      y
    )

    return(
      paste0(
        "W",
        as.integer(n)
      )
    )
  }

  NA_character_
}

extract_week_from_title <- function(x) {

  x <- clean_chr(x)

  if (is.na(x)) {
    return(NA_character_)
  }

  y <- tolower(x)

  if (
    grepl(
      "w4/6",
      y,
      fixed = TRUE
    ) ||
    grepl(
      "w4-6",
      y,
      fixed = TRUE
    ) ||
    grepl(
      "w4_6",
      y,
      fixed = TRUE
    )
  ) {
    return("W4_6")
  }

  m <- regexpr(
    "w[0-9]+",
    y,
    perl = TRUE
  )

  if (m[1] >= 0) {

    hit <- regmatches(
      y,
      m
    )

    return(
      toupper(hit)
    )
  }

  NA_character_
}


mayo_to_response <- function(x) {

  x <- clean_chr(x)

  if (is.na(x)) {
    return(NA_character_)
  }

  suppressWarnings(
    value <- as.numeric(x)
  )

  if (!is.finite(value)) {
    return(NA_character_)
  }

  if (value %in% c(0, 1)) {
    return("Response")
  }

  if (value %in% c(2, 3)) {
    return("Non-response")
  }

  NA_character_
}


title_response <- function(x) {

  y <- tolower(
    as.character(x)
  )

  if (
    grepl(
      "uc nr",
      y,
      fixed = TRUE
    )
  ) {
    return("Non-response")
  }

  if (
    grepl(
      "uc r",
      y,
      fixed = TRUE
    )
  ) {
    return("Response")
  }

  NA_character_
}


collapse_status <- function(x) {

  x <- unique(
    x[
      !is.na(x) &
        x %in%
        c(
          "Response",
          "Non-response"
        )
    ]
  )

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(
    sort(x),
    collapse = ";"
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 0B v3: longitudinal cohort reconstruction\n")
cat("============================================================\n\n")


# ============================================================
# 1. Inputs
# ============================================================

check(
  file.exists(soft_file),
  "GEO SOFT exists"
)

check(
  file.exists(matrix_file),
  "Series Matrix exists"
)


# ============================================================
# 2. Extract sample metadata only
# ============================================================

pattern <- paste0(
  "^\\^SAMPLE =",
  "|^!Sample_geo_accession =",
  "|^!Sample_title =",
  "|^!Sample_source_name_ch1 =",
  "|^!Sample_organism_ch1 =",
  "|^!Sample_description =",
  "|^!Sample_characteristics_ch1 =",
  "|^!Sample_supplementary_file"
)

cmd <- paste(
  "zgrep -E",
  shQuote(pattern),
  shQuote(soft_file)
)

lines <- system(
  cmd,
  intern = TRUE
)

check(
  length(lines) > 0,
  "SOFT sample metadata extracted"
)


# ============================================================
# 3. Parse sample blocks
# ============================================================

rows <- list()
row_counter <- 0L
current <- NULL


flush_sample <- function(x) {

  if (is.null(x)) {
    return(NULL)
  }

  patient <- extract_characteristic_regex(
    x$Characteristics,
    c(
      "study.*individual.*number",
      "individual.*number",
      "^patient.*number$"
    )
  )

  week <- extract_characteristic_regex(
    x$Characteristics,
    c(
      "^week",
      "^timepoint",
      "^time point"
    )
  )

  treatment <- extract_characteristic_regex(
    x$Characteristics,
    c(
      "induction.*therapy.*maintenance.*therapy",
      "induction.*maintenance",
      "treatment.*schedule"
    )
  )

  mayo <- extract_characteristic_regex(
    x$Characteristics,
    c(
      "mayo.*endoscopic.*subscore",
      "endoscopic.*subscore"
    )
  )

  tissue <- extract_characteristic_regex(
    x$Characteristics,
    c(
      "^tissue$",
      "tissue.*type"
    )
  )

  data.table(
    GSM =
      x$GSM,

    GEOAccession =
      x$GEOAccession,

    Title =
      x$Title,

    SourceName =
      x$SourceName,

    Organism =
      x$Organism,

    PatientNumberRaw =
      patient,

    WeekRaw =
      week,

    TreatmentScheduleRaw =
      treatment,

    MayoEndoscopicSubscoreRaw =
      mayo,

    TissueRaw =
      tissue,

    Characteristics =
      paste(
        x$Characteristics,
        collapse = " || "
      ),

    Description =
      paste(
        x$Description,
        collapse = " || "
      ),

    SupplementaryFile =
      paste(
        x$SupplementaryFile,
        collapse = " || "
      )
  )
}


for (line in lines) {

  if (
    grepl(
      "^\\^SAMPLE =",
      line
    )
  ) {

    if (!is.null(current)) {

      row_counter <-
        row_counter + 1L

      rows[[row_counter]] <-
        flush_sample(current)
    }

    current <- list(
      GSM =
        trimws(
          sub(
            "^\\^SAMPLE =[[:space:]]*",
            "",
            line
          )
        ),

      GEOAccession =
        NA_character_,

      Title =
        NA_character_,

      SourceName =
        NA_character_,

      Organism =
        NA_character_,

      Characteristics =
        character(),

      Description =
        character(),

      SupplementaryFile =
        character()
    )

    next
  }


  if (is.null(current)) {
    next
  }


  if (
    grepl(
      "^!Sample_geo_accession =",
      line
    )
  ) {

    current$GEOAccession <-
      soft_value(line)

  } else if (
    grepl(
      "^!Sample_title =",
      line
    )
  ) {

    current$Title <-
      soft_value(line)

  } else if (
    grepl(
      "^!Sample_source_name_ch1 =",
      line
    )
  ) {

    current$SourceName <-
      soft_value(line)

  } else if (
    grepl(
      "^!Sample_organism_ch1 =",
      line
    )
  ) {

    current$Organism <-
      soft_value(line)

  } else if (
    grepl(
      "^!Sample_characteristics_ch1 =",
      line
    )
  ) {

    current$Characteristics <- c(
      current$Characteristics,
      soft_value(line)
    )

  } else if (
    grepl(
      "^!Sample_description =",
      line
    )
  ) {

    current$Description <- c(
      current$Description,
      soft_value(line)
    )

  } else if (
    grepl(
      "^!Sample_supplementary_file",
      line
    )
  ) {

    current$SupplementaryFile <- c(
      current$SupplementaryFile,
      soft_value(line)
    )
  }
}


if (!is.null(current)) {

  row_counter <-
    row_counter + 1L

  rows[[row_counter]] <-
    flush_sample(current)
}


samples <- rbindlist(
  rows,
  fill = TRUE
)


check(
  nrow(samples) == 178,
  "178 SOFT samples parsed"
)

check(
  uniqueN(samples$GSM) == 178,
  "178 unique GSM IDs"
)

check(
  all(
    samples$GSM ==
      samples$GEOAccession
  ),
  "GSM and Sample_geo_accession agree"
)

check(
  all(
    samples$Organism ==
      "Homo sapiens"
  ),
  "All samples are Homo sapiens"
)


# ============================================================
# 4. Series Matrix concordance
# ============================================================

matrix_line <- system(
  paste(
    "zgrep -m1",
    shQuote(
      "^!Sample_geo_accession"
    ),
    shQuote(
      matrix_file
    )
  ),
  intern = TRUE
)


check(
  length(matrix_line) == 1L,
  "Series Matrix GSM row found"
)


matrix_gsm <- strsplit(
  matrix_line,
  "\t",
  fixed = TRUE
)[[1]][-1]


matrix_gsm <- gsub(
  '^"|"$',
  "",
  matrix_gsm
)


check(
  length(matrix_gsm) == 178,
  "Series Matrix has 178 GSMs"
)

check(
  identical(
    samples$GSM,
    matrix_gsm
  ),
  "SOFT and Series Matrix GSM order identical"
)


# ============================================================
# 5. CEL provenance
# ============================================================

samples[
  ,
  CELFile :=
    sub(
      ".*/",
      "",
      sub(
        "[[:space:]]*\\|\\|.*$",
        "",
        SupplementaryFile
      )
    )
]


check(
  all(
    grepl(
      "\\.CEL\\.gz$",
      samples$CELFile,
      ignore.case = TRUE
    )
  ),
  "All 178 CEL.gz references recovered"
)


check(
  uniqueN(
    samples$CELFile
  ) == 178,
  "All CEL filenames unique"
)


# ============================================================
# 6. Patient identity
# ============================================================

samples[
  ,
  PatientNumber :=
    clean_chr(
      PatientNumberRaw
    )
]


fallback_patient <- ifelse(
  grepl(
    "_[0-9]+$",
    samples$Title
  ),
  sub(
    ".*_([0-9]+)$",
    "\\1",
    samples$Title
  ),
  NA_character_
)


need_patient <- is.na(
  samples$PatientNumber
)


samples[
  need_patient,
  PatientNumber :=
    fallback_patient[
      need_patient
    ]
]


check(
  !anyNA(
    samples$PatientNumber
  ),
  "Patient number recovered for all samples"
)


# ============================================================
# 7. Study branch
# ============================================================

title_lower <- tolower(
  samples$Title
)


samples[
  ,
  Branch :=
    fcase(

      grepl(
        "control",
        title_lower,
        fixed = TRUE
      ),
      "CONTROL",

      grepl(
        "ifx",
        title_lower,
        fixed = TRUE
      ),
      "IFX",

      grepl(
        "vdz",
        title_lower,
        fixed = TRUE
      ) |
        grepl(
          "plac",
          title_lower,
          fixed = TRUE
        ),
      "VDZ_TRIAL",

      default =
        "REVIEW"
    )
]


check(
  sum(
    samples$Branch ==
      "REVIEW"
  ) == 0,
  "All 178 samples assigned to study branch"
)


samples[
  ,
  PatientKey :=
    paste(
      Branch,
      PatientNumber,
      sep = "_"
    )
]


samples[
  ,
  Disease :=
    fifelse(
      Branch ==
        "CONTROL",
      "Non-IBD",
      "UC"
    )
]


# ============================================================
# 8. Timepoint reconstruction
# ============================================================

samples[
  ,
  TimepointFromCharacteristics :=
    vapply(
      WeekRaw,
      normalize_week_value,
      character(1)
    )
]


samples[
  ,
  TimepointFromTitle :=
    vapply(
      Title,
      extract_week_from_title,
      character(1)
    )
]


samples[
  ,
  Timepoint :=
    TimepointFromCharacteristics
]


need_time <- is.na(
  samples$Timepoint
)


samples[
  need_time,
  Timepoint :=
    TimepointFromTitle[
      need_time
    ]
]


# Controls are not longitudinal treatment visits.
samples[
  Branch == "CONTROL" &
    is.na(Timepoint),
  Timepoint :=
    "CONTROL"
]


cat("\n============================================================\n")
cat("RAW / RECONSTRUCTED TIMEPOINT AUDIT\n")
cat("============================================================\n")

print(
  samples[
    ,
    .N,
    by = .(
      Branch,
      WeekRaw,
      Timepoint
    )
  ][
    order(
      Branch,
      Timepoint,
      WeekRaw
    )
  ]
)


unresolved_time <- samples[
  Branch != "CONTROL" &
    is.na(Timepoint),
  .(
    GSM,
    Title,
    Branch,
    WeekRaw,
    Characteristics
  )
]


if (
  nrow(
    unresolved_time
  ) > 0
) {

  cat("\nUNRESOLVED TREATMENT TIMEPOINTS:\n")
  print(
    unresolved_time
  )
}


check(
  nrow(
    unresolved_time
  ) == 0,
  "All treatment samples have interpretable timepoints"
)


check(
  all(
    samples[
      Branch == "CONTROL",
      Timepoint
    ] ==
      "CONTROL"
  ),
  "Controls explicitly assigned Timepoint=CONTROL"
)


# ============================================================
# 9. Treatment schedule
# ============================================================

samples[
  ,
  TreatmentSchedule :=
    tolower(
      clean_chr(
        TreatmentScheduleRaw
      )
    )
]


cat("\n============================================================\n")
cat("TREATMENT SCHEDULE VALUES\n")
cat("============================================================\n")

print(
  samples[
    ,
    .N,
    by = .(
      Branch,
      TreatmentSchedule
    )
  ][
    order(
      Branch,
      TreatmentSchedule
    )
  ]
)


# ============================================================
# 10. Mayo and response
# ============================================================

samples[
  ,
  MayoEndoscopicSubscore :=
    suppressWarnings(
      as.numeric(
        MayoEndoscopicSubscoreRaw
      )
    )
]


samples[
  ,
  ResponseFromMayo :=
    vapply(
      MayoEndoscopicSubscoreRaw,
      mayo_to_response,
      character(1)
    )
]


samples[
  ,
  ResponseFromTitle :=
    vapply(
      Title,
      title_response,
      character(1)
    )
]


samples[
  ,
  ResponseAtVisit :=
    fifelse(
      Branch %in%
        c(
          "VDZ_TRIAL",
          "IFX"
        ) &
        !Timepoint %in%
          c(
            "W0",
            "CONTROL"
          ),
      ResponseFromMayo,
      NA_character_
    )
]


# ============================================================
# 11. Mayo vs title response
# ============================================================

response_compare <- samples[
  !is.na(
    ResponseFromMayo
  ) &
    !is.na(
      ResponseFromTitle
    ),
  .(
    GSM,
    Title,
    Branch,
    Timepoint,
    MayoEndoscopicSubscore,
    ResponseFromMayo,
    ResponseFromTitle,
    Concordant =
      ResponseFromMayo ==
        ResponseFromTitle
  )
]


cat("\n============================================================\n")
cat("MAYO vs TITLE RESPONSE CROSS-CHECK\n")
cat("============================================================\n")

cat(
  "Comparable samples:",
  nrow(
    response_compare
  ),
  "\n"
)

cat(
  "Discordant samples:",
  sum(
    !response_compare$Concordant
  ),
  "\n"
)


check(
  sum(
    !response_compare$Concordant
  ) == 0,
  "Mayo-derived and GEO-title response labels agree"
)


# ============================================================
# 12. Branch summary
# ============================================================

branch_summary <- samples[
  ,
  .(
    Samples = .N,
    UniquePatients =
      uniqueN(
        PatientKey
      )
  ),
  by = Branch
][
  order(
    Branch
  )
]


cat("\n============================================================\n")
cat("STUDY BRANCH SUMMARY\n")
cat("============================================================\n")

print(
  branch_summary
)


check(
  samples[
    Branch == "CONTROL",
    .N
  ] == 12,
  "Control samples = 12"
)

check(
  samples[
    Branch == "IFX",
    .N
  ] == 46,
  "IFX samples = 46"
)

check(
  samples[
    Branch == "VDZ_TRIAL",
    .N
  ] == 120,
  "VDZ-trial samples = 120"
)

check(
  uniqueN(
    samples[
      Branch == "IFX",
      PatientKey
    ]
  ) == 23,
  "IFX patients = 23"
)

check(
  uniqueN(
    samples[
      Branch == "VDZ_TRIAL",
      PatientKey
    ]
  ) == 44,
  "VDZ-trial patients = 44"
)


# ============================================================
# 13. IFX longitudinal structure
# ============================================================

ifx <- samples[
  Branch == "IFX"
]


check(
  ifx[
    Timepoint == "W0",
    .N
  ] == 23,
  "IFX W0 = 23"
)

check(
  ifx[
    Timepoint == "W4_6",
    .N
  ] == 23,
  "IFX W4/6 = 23"
)


ifx_pairs <- ifx[
  ,
  .(
    N = .N,

    W0 =
      sum(
        Timepoint ==
          "W0"
      ),

    W4_6 =
      sum(
        Timepoint ==
          "W4_6"
      )
  ),
  by = PatientKey
]


check(
  all(
    ifx_pairs$N == 2 &
      ifx_pairs$W0 == 1 &
      ifx_pairs$W4_6 == 1
  ),
  "All IFX patients have one W0 and one W4/6 sample"
)


# ============================================================
# 14. VDZ baseline
# ============================================================

vdz_baseline <- samples[
  Branch == "VDZ_TRIAL" &
    Timepoint == "W0"
]


check(
  nrow(
    vdz_baseline
  ) == 44,
  "VDZ-trial W0 samples = 44"
)

check(
  uniqueN(
    vdz_baseline$PatientKey
  ) == 44,
  "Exactly one W0 sample per VDZ-trial patient"
)


# ============================================================
# 15. Baseline VDZ exposure
# ============================================================

vdz_baseline[
  ,
  ExposureClass :=
    fcase(

      !is.na(
        TreatmentSchedule
      ) &
        grepl(
          "plac.*plac",
          TreatmentSchedule,
          perl = TRUE
        ),
      "Placebo_only",

      !is.na(
        TreatmentSchedule
      ) &
        grepl(
          "vdz",
          TreatmentSchedule,
          fixed = TRUE
        ),
      "VDZ_exposed",

      default =
        "REVIEW"
    )
]


cat("\n============================================================\n")
cat("VDZ BASELINE EXPOSURE SUMMARY\n")
cat("============================================================\n")

print(
  vdz_baseline[
    ,
    .N,
    by = .(
      TreatmentSchedule,
      ExposureClass
    )
  ][
    order(
      ExposureClass,
      TreatmentSchedule
    )
  ]
)


check(
  vdz_baseline[
    ExposureClass ==
      "VDZ_exposed",
    .N
  ] == 41,
  "VDZ-exposed baseline patients = 41"
)

check(
  vdz_baseline[
    ExposureClass ==
      "Placebo_only",
    .N
  ] == 3,
  "Placebo-only baseline patients = 3"
)

check(
  vdz_baseline[
    ExposureClass ==
      "REVIEW",
    .N
  ] == 0,
  "No baseline exposure requires review"
)


# ============================================================
# 16. Propagate baseline class
# ============================================================

baseline_class <- vdz_baseline[
  ,
  .(
    PatientKey,

    BaselineGSM =
      GSM,

    BaselineTitle =
      Title,

    BaselineTreatmentSchedule =
      TreatmentSchedule,

    ExposureClass
  )
]


samples[
  baseline_class,
  on = "PatientKey",
  `:=`(
    BaselineGSM =
      i.BaselineGSM,

    BaselineTitle =
      i.BaselineTitle,

    BaselineTreatmentSchedule =
      i.BaselineTreatmentSchedule,

    ExposureClass =
      i.ExposureClass
  )
]


# ============================================================
# 17. Placebo-only list
# ============================================================

placebo_only <- unique(
  samples[
    Branch == "VDZ_TRIAL" &
      ExposureClass ==
        "Placebo_only",
    .(
      PatientKey,
      PatientNumber,
      BaselineGSM,
      BaselineTreatmentSchedule,
      BaselineTitle
    )
  ]
)


cat("\n============================================================\n")
cat("PLACEBO-ONLY PATIENTS\n")
cat("============================================================\n")

print(
  placebo_only
)


# ============================================================
# 18. VDZ response
# ============================================================

vdz_response <- samples[
  Branch == "VDZ_TRIAL" &
    ExposureClass ==
      "VDZ_exposed" &
    !is.na(
      ResponseAtVisit
    )
]


endpoint_sample_counts <- vdz_response[
  ,
  .N,
  by = .(
    Timepoint,
    ResponseAtVisit
  )
][
  order(
    Timepoint,
    ResponseAtVisit
  )
]


endpoint_patient_counts <- unique(
  vdz_response[
    ,
    .(
      PatientKey,
      Timepoint,
      ResponseAtVisit
    )
  ]
)[
  ,
  .N,
  by = .(
    Timepoint,
    ResponseAtVisit
  )
][
  order(
    Timepoint,
    ResponseAtVisit
  )
]


cat("\n============================================================\n")
cat("VDZ RESPONSE — SAMPLE COUNTS\n")
cat("============================================================\n")

print(
  endpoint_sample_counts
)


cat("\n============================================================\n")
cat("VDZ RESPONSE — UNIQUE PATIENT COUNTS\n")
cat("============================================================\n")

print(
  endpoint_patient_counts
)


# ============================================================
# 19. Response conflicts
# ============================================================

response_conflicts <- vdz_response[
  ,
  .(
    N_samples =
      .N,

    N_labels =
      uniqueN(
        ResponseAtVisit
      ),

    Labels =
      paste(
        sort(
          unique(
            ResponseAtVisit
          )
        ),
        collapse = ";"
      ),

    GSMs =
      paste(
        GSM,
        collapse = ";"
      )
  ),
  by = .(
    PatientKey,
    Timepoint
  )
][
  N_labels > 1
]


check(
  nrow(
    response_conflicts
  ) == 0,
  "No conflicting VDZ response labels within patient/timepoint"
)


# ============================================================
# 20. Repeated patient x timepoint
# ============================================================

repeats <- samples[
  Branch %in%
    c(
      "VDZ_TRIAL",
      "IFX"
    ),
  .(
    N_samples = .N,

    GSMs =
      paste(
        GSM,
        collapse = ";"
      ),

    Titles =
      paste(
        Title,
        collapse = " || "
      )
  ),
  by = .(
    Branch,
    PatientKey,
    Timepoint
  )
][
  N_samples > 1
]


cat("\n============================================================\n")
cat("REPEATED PATIENT x TIMEPOINT COMBINATIONS\n")
cat("============================================================\n")

if (
  nrow(
    repeats
  ) == 0
) {

  cat("None\n")

} else {

  print(
    repeats
  )
}


# ============================================================
# 21. Patient-level VDZ table
# ============================================================

vdz_patient <- samples[
  Branch == "VDZ_TRIAL",
  .(
    PatientNumber =
      PatientNumber[1],

    BaselineGSM =
      BaselineGSM[1],

    BaselineTreatmentSchedule =
      BaselineTreatmentSchedule[1],

    ExposureClass =
      ExposureClass[1],

    N_samples =
      .N,

    W6Response =
      collapse_status(
        ResponseAtVisit[
          Timepoint ==
            "W6"
        ]
      ),

    W12Response =
      collapse_status(
        ResponseAtVisit[
          Timepoint ==
            "W12"
        ]
      ),

    W52Response =
      collapse_status(
        ResponseAtVisit[
          Timepoint ==
            "W52"
        ]
      )
  ),
  by = PatientKey
]


check(
  nrow(
    vdz_patient
  ) == 44,
  "VDZ patient table contains 44 patients"
)

check(
  vdz_patient[
    ExposureClass ==
      "VDZ_exposed",
    .N
  ] == 41,
  "VDZ patient table contains 41 exposed patients"
)


candidate <- vdz_patient[
  ExposureClass ==
    "VDZ_exposed"
]


check(
  nrow(
    candidate
  ) == 41,
  "Candidate VDZ baseline cohort contains 41 patients"
)

check(
  uniqueN(
    candidate$BaselineGSM
  ) == 41,
  "Candidate VDZ baseline cohort contains 41 unique GSMs"
)


# ============================================================
# 22. Endpoint availability
# ============================================================

endpoint_availability <- data.table(
  Endpoint =
    c(
      "W6",
      "W12",
      "W52"
    ),

  PatientsWithInterpretableResponse =
    c(
      sum(
        !is.na(
          candidate$W6Response
        )
      ),

      sum(
        !is.na(
          candidate$W12Response
        )
      ),

      sum(
        !is.na(
          candidate$W52Response
        )
      )
    ),

  Response =
    c(
      sum(
        candidate$W6Response ==
          "Response",
        na.rm = TRUE
      ),

      sum(
        candidate$W12Response ==
          "Response",
        na.rm = TRUE
      ),

      sum(
        candidate$W52Response ==
          "Response",
        na.rm = TRUE
      )
    ),

  NonResponse =
    c(
      sum(
        candidate$W6Response ==
          "Non-response",
        na.rm = TRUE
      ),

      sum(
        candidate$W12Response ==
          "Non-response",
        na.rm = TRUE
      ),

      sum(
        candidate$W52Response ==
          "Non-response",
        na.rm = TRUE
      )
    )
)


cat("\n============================================================\n")
cat("ENDPOINT AVAILABILITY AMONG 41 VDZ-EXPOSED PATIENTS\n")
cat("============================================================\n")

print(
  endpoint_availability
)


# ============================================================
# 23. Save outputs
# ============================================================

fwrite(
  samples,
  file.path(
    out_dir,
    "GSE73661_cohort_reconstruction_178_samples.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  branch_summary,
  file.path(
    out_dir,
    "GSE73661_study_branch_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  vdz_baseline,
  file.path(
    out_dir,
    "GSE73661_vdz_baseline_44_samples.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  placebo_only,
  file.path(
    out_dir,
    "GSE73661_placebo_only_patients.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  vdz_patient,
  file.path(
    out_dir,
    "GSE73661_vdz_patient_reconstruction_44.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  candidate,
  file.path(
    out_dir,
    "GSE73661_vdz_exposed_baseline_candidate_41.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  endpoint_sample_counts,
  file.path(
    out_dir,
    "GSE73661_vdz_endpoint_sample_counts.tsv"
  ),
  sep = "\t"
)

fwrite(
  endpoint_patient_counts,
  file.path(
    out_dir,
    "GSE73661_vdz_endpoint_patient_counts.tsv"
  ),
  sep = "\t"
)

fwrite(
  endpoint_availability,
  file.path(
    out_dir,
    "GSE73661_vdz_endpoint_availability.tsv"
  ),
  sep = "\t"
)

fwrite(
  repeats,
  file.path(
    out_dir,
    "GSE73661_repeated_patient_timepoint_samples.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

fwrite(
  response_conflicts,
  file.path(
    out_dir,
    "GSE73661_response_label_conflicts.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

fwrite(
  response_compare,
  file.path(
    out_dir,
    "GSE73661_mayo_vs_title_response_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 24. Summary
# ============================================================

summary <- data.table(
  Metric =
    c(
      "Total_samples",
      "Control_samples",
      "IFX_samples",
      "VDZ_trial_samples",
      "IFX_patients",
      "VDZ_trial_patients",
      "VDZ_baseline_samples",
      "VDZ_exposed_baseline_patients",
      "Placebo_only_baseline_patients",
      "Candidate_VDZ_baseline_samples",
      "Mayo_title_response_discordance",
      "VDZ_response_label_conflicts",
      "Repeated_patient_timepoint_combinations"
    ),

  N =
    c(
      nrow(samples),

      samples[
        Branch == "CONTROL",
        .N
      ],

      samples[
        Branch == "IFX",
        .N
      ],

      samples[
        Branch == "VDZ_TRIAL",
        .N
      ],

      uniqueN(
        samples[
          Branch == "IFX",
          PatientKey
        ]
      ),

      uniqueN(
        samples[
          Branch == "VDZ_TRIAL",
          PatientKey
        ]
      ),

      nrow(
        vdz_baseline
      ),

      vdz_baseline[
        ExposureClass ==
          "VDZ_exposed",
        .N
      ],

      vdz_baseline[
        ExposureClass ==
          "Placebo_only",
        .N
      ],

      nrow(
        candidate
      ),

      sum(
        !response_compare$Concordant
      ),

      nrow(
        response_conflicts
      ),

      nrow(
        repeats
      )
    )
)


cat("\n============================================================\n")
cat("FINAL COHORT RECONSTRUCTION SUMMARY\n")
cat("============================================================\n")

print(
  summary
)


fwrite(
  summary,
  file.path(
    out_dir,
    "GSE73661_stage0b_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 25. Final gate
# ============================================================

critical_pass <- all(
  c(
    nrow(samples) == 178,

    uniqueN(
      samples$GSM
    ) == 178,

    samples[
      Branch == "CONTROL",
      .N
    ] == 12,

    samples[
      Branch == "IFX",
      .N
    ] == 46,

    samples[
      Branch == "VDZ_TRIAL",
      .N
    ] == 120,

    uniqueN(
      samples[
        Branch == "IFX",
        PatientKey
      ]
    ) == 23,

    uniqueN(
      samples[
        Branch == "VDZ_TRIAL",
        PatientKey
      ]
    ) == 44,

    nrow(
      unresolved_time
    ) == 0,

    nrow(
      vdz_baseline
    ) == 44,

    vdz_baseline[
      ExposureClass ==
        "VDZ_exposed",
      .N
    ] == 41,

    vdz_baseline[
      ExposureClass ==
        "Placebo_only",
      .N
    ] == 3,

    nrow(
      candidate
    ) == 41,

    nrow(
      response_conflicts
    ) == 0,

    sum(
      !response_compare$Concordant
    ) == 0
  )
)


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] 178 GEO samples reconstructed\n")
  cat("[PASS] Controls explicitly treated as non-longitudinal CONTROL timepoint\n")
  cat("[PASS] All treatment samples have interpretable timepoints\n")
  cat("[PASS] 12 controls identified\n")
  cat("[PASS] 46 IFX samples / 23 IFX patients identified\n")
  cat("[PASS] 120 VDZ-trial samples / 44 patients identified\n")
  cat("[PASS] 44 VDZ-trial W0 samples identified\n")
  cat("[PASS] 41 VDZ-exposed baseline patients identified\n")
  cat("[PASS] 3 placebo-only baseline patients identified\n")
  cat("[PASS] Candidate VDZ baseline cohort = 41\n")
  cat("[PASS] Mayo-derived response agrees with title response\n")
  cat("[PASS] No VDZ response-label conflicts\n")
  cat("[INFO] W6/W12/W52 endpoint remains unfrozen\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE73661_STAGE0B_LONGITUDINAL_COHORT_RECONSTRUCTED\n")

} else {

  cat("[REVIEW] Stage 0B requires review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE73661_STAGE0B_LONGITUDINAL_COHORT\n")
}

cat("============================================================\n")
