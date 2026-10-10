suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 4
# FINAL canonical metadata construction
#
# Frozen upstream sources:
#
#   Stage 0B:
#     longitudinal cohort reconstruction
#
#   Stage 3:
#     canonical expression
#     18,407 genes x 178 samples
#
# Principles:
#
#   1. GSE73661_cohort_reconstruction_178_samples.tsv
#      is the master sample-level metadata source.
#
#   2. Canonical expression column order defines FINAL
#      metadata sample order.
#
#   3. VDZ patient-level W6 / W12 / W52 responses are taken
#      directly from frozen Stage-0B patient reconstruction.
#
#   4. Clinical endpoint labels are attached ONLY to
#      VDZ-exposed baseline W0 predictor samples.
#
#   5. Endpoint definitions remain frozen:
#
#        W6  = PRIMARY
#        W52 = SECONDARY
#        W12 = EXPLORATORY
#
#   6. No sample is removed in Stage 4.
#      Stage 5 will perform technical sample QC.
#
#   7. No expression values are touched.
# ============================================================


# ============================================================
# 0. Paths
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

cohort_dir <- file.path(
  root,
  "prepared/00_cohort_reconstruction"
)

expr_file <- file.path(
  root,
  "prepared/03_canonical_expression",
  "GSE73661_canonical_log2_RMA_18407x178.rds"
)

master_file <- file.path(
  cohort_dir,
  "GSE73661_cohort_reconstruction_178_samples.tsv"
)

vdz_baseline_file <- file.path(
  cohort_dir,
  "GSE73661_vdz_baseline_44_samples.tsv"
)

vdz_exposed_file <- file.path(
  cohort_dir,
  "GSE73661_vdz_exposed_baseline_candidate_41.tsv"
)

vdz_patient_file <- file.path(
  cohort_dir,
  "GSE73661_vdz_patient_reconstruction_44.tsv"
)

endpoint_availability_file <- file.path(
  cohort_dir,
  "GSE73661_vdz_endpoint_availability.tsv"
)

out_dir <- file.path(
  root,
  "prepared/04_canonical_metadata"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. Helpers
# ============================================================

check <- function(condition, message) {

  if (!isTRUE(condition)) {
    stop(
      paste0("[FAIL] ", message),
      call. = FALSE
    )
  }

  cat(
    "[PASS] ",
    message,
    "\n",
    sep = ""
  )
}


clean_chr <- function(x) {

  x <- trimws(
    as.character(x)
  )

  bad <- is.na(x) |
    x %in% c(
      "",
      "-",
      "NA",
      "N/A",
      "NULL",
      "null"
    )

  x[bad] <- NA_character_

  x
}


clean_character_columns <- function(dt) {

  char_cols <- names(dt)[
    vapply(
      dt,
      is.character,
      logical(1)
    )
  ]

  for (cc in char_cols) {

    set(
      dt,
      j = cc,
      value = clean_chr(
        dt[[cc]]
      )
    )
  }

  invisible(dt)
}


same_with_na <- function(x, y) {

  x <- clean_chr(x)
  y <- clean_chr(y)

  if (length(x) != length(y)) {
    return(FALSE)
  }

  all(
    (is.na(x) & is.na(y)) |
      (
        !is.na(x) &
        !is.na(y) &
        x == y
      )
  )
}


response_binary <- function(x) {

  x <- clean_chr(x)

  out <- rep(
    NA_integer_,
    length(x)
  )

  out[
    x == "Response"
  ] <- 1L

  out[
    x == "Non-response"
  ] <- 0L

  out
}


check_response_counts <- function(
  dt,
  column,
  expected_total,
  expected_response,
  expected_nonresponse,
  label
) {

  x <- clean_chr(
    dt[[column]]
  )

  n_total <- sum(
    !is.na(x)
  )

  n_response <- sum(
    x == "Response",
    na.rm = TRUE
  )

  n_nonresponse <- sum(
    x == "Non-response",
    na.rm = TRUE
  )

  check(
    n_total == expected_total,
    paste0(
      label,
      " interpretable N = ",
      expected_total
    )
  )

  check(
    n_response == expected_response,
    paste0(
      label,
      " Response N = ",
      expected_response
    )
  )

  check(
    n_nonresponse == expected_nonresponse,
    paste0(
      label,
      " Non-response N = ",
      expected_nonresponse
    )
  )

  check(
    n_response +
      n_nonresponse ==
      n_total,
    paste0(
      label,
      " contains only Response / Non-response labels"
    )
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 4: FINAL canonical metadata\n")
cat("============================================================\n\n")


# ============================================================
# 2. Input existence
# ============================================================

for (f in c(
  expr_file,
  master_file,
  vdz_baseline_file,
  vdz_exposed_file,
  vdz_patient_file,
  endpoint_availability_file
)) {

  check(
    file.exists(f),
    paste0(
      "Input exists: ",
      basename(f)
    )
  )
}


# ============================================================
# 3. Load frozen inputs
# ============================================================

expr <- readRDS(
  expr_file
)

master <- fread(
  master_file,
  showProgress = FALSE
)

vdz_baseline <- fread(
  vdz_baseline_file,
  showProgress = FALSE
)

vdz_exposed <- fread(
  vdz_exposed_file,
  showProgress = FALSE
)

vdz_patient <- fread(
  vdz_patient_file,
  showProgress = FALSE
)

endpoint_availability <- fread(
  endpoint_availability_file,
  showProgress = FALSE
)


check(
  identical(
    dim(expr),
    c(
      18407L,
      178L
    )
  ),
  "Canonical expression = 18,407 x 178"
)


check(
  nrow(master) == 178L,
  "Master Stage-0B metadata contains 178 samples"
)


# ============================================================
# 4. Clean character fields
# ============================================================

clean_character_columns(master)
clean_character_columns(vdz_baseline)
clean_character_columns(vdz_exposed)
clean_character_columns(vdz_patient)
clean_character_columns(endpoint_availability)


# ============================================================
# 5. Freeze identifier types
# ============================================================

master[
  ,
  GSM :=
    as.character(GSM)
]

master[
  ,
  GEOAccession :=
    as.character(GEOAccession)
]

master[
  ,
  PatientKey :=
    as.character(PatientKey)
]


vdz_baseline[
  ,
  GSM :=
    as.character(GSM)
]

vdz_baseline[
  ,
  PatientKey :=
    as.character(PatientKey)
]


vdz_exposed[
  ,
  PatientKey :=
    as.character(PatientKey)
]

vdz_exposed[
  ,
  BaselineGSM :=
    as.character(BaselineGSM)
]


vdz_patient[
  ,
  PatientKey :=
    as.character(PatientKey)
]

vdz_patient[
  ,
  BaselineGSM :=
    as.character(BaselineGSM)
]


# ============================================================
# 6. Validate master sample universe
# ============================================================

expr_samples <- as.character(
  colnames(expr)
)


check(
  length(expr_samples) == 178L,
  "Canonical expression contains 178 sample IDs"
)


check(
  uniqueN(expr_samples) == 178L,
  "Canonical expression sample IDs are unique"
)


check(
  uniqueN(master$GSM) == 178L,
  "Master metadata GSM IDs are unique"
)


check(
  setequal(
    expr_samples,
    master$GSM
  ),
  "Expression and master metadata contain identical GSM universe"
)


# ============================================================
# 7. Align metadata EXACTLY to canonical expression
# ============================================================

metadata_idx <- match(
  expr_samples,
  master$GSM
)


check(
  !anyNA(metadata_idx),
  "Every expression sample maps to master metadata"
)


meta <- copy(
  master[
    metadata_idx
  ]
)


check(
  identical(
    meta$GSM,
    expr_samples
  ),
  "Canonical metadata order exactly matches expression columns"
)


meta[
  ,
  SampleOrder :=
    seq_len(.N)
]


# ============================================================
# 8. Core cohort structure
# ============================================================

check(
  meta[
    Branch == "CONTROL",
    .N
  ] == 12L,
  "CONTROL branch contains 12 samples"
)


check(
  meta[
    Branch == "IFX",
    .N
  ] == 46L,
  "IFX branch contains 46 samples"
)


check(
  meta[
    Branch == "VDZ_TRIAL",
    .N
  ] == 120L,
  "VDZ_TRIAL branch contains 120 samples"
)


check(
  uniqueN(
    meta[
      Branch == "CONTROL",
      PatientKey
    ]
  ) == 11L,
  "CONTROL branch contains 11 patient labels"
)


check(
  uniqueN(
    meta[
      Branch == "IFX",
      PatientKey
    ]
  ) == 23L,
  "IFX branch contains 23 patients"
)


check(
  uniqueN(
    meta[
      Branch == "VDZ_TRIAL",
      PatientKey
    ]
  ) == 44L,
  "VDZ_TRIAL branch contains 44 patients"
)


check(
  all(
    meta[
      Branch == "IFX",
      Disease
    ] ==
      "UC"
  ),
  "All IFX samples are UC"
)


check(
  all(
    meta[
      Branch == "VDZ_TRIAL",
      Disease
    ] ==
      "UC"
  ),
  "All VDZ trial samples are UC"
)


check(
  all(
    meta[
      Branch == "CONTROL",
      Disease
    ] ==
      "Non-IBD"
  ),
  "All CONTROL samples are Non-IBD"
)


# ============================================================
# 9. IFX schedule validation
# ============================================================

check(
  meta[
    Branch == "IFX" &
      Timepoint == "W0",
    .N
  ] == 23L,
  "IFX contains 23 W0 samples"
)


check(
  meta[
    Branch == "IFX" &
      Timepoint == "W4_6",
    .N
  ] == 23L,
  "IFX contains 23 W4/6 samples"
)


# ============================================================
# 10. VDZ longitudinal timepoint validation
# ============================================================

vdz_expected_timepoints <- data.table(
  Timepoint = c(
    "W0",
    "W6",
    "W8",
    "W9",
    "W10",
    "W12",
    "W14",
    "W16",
    "W18",
    "W26",
    "W28",
    "W32",
    "W38",
    "W52"
  ),

  ExpectedN = c(
    44L,
    30L,
    4L,
    1L,
    1L,
    13L,
    1L,
    1L,
    1L,
    1L,
    2L,
    1L,
    1L,
    19L
  )
)


vdz_observed_timepoints <- meta[
  Branch == "VDZ_TRIAL",
  .(
    ObservedN = .N
  ),
  by = Timepoint
]


vdz_timepoint_check <- merge(
  vdz_expected_timepoints,
  vdz_observed_timepoints,
  by = "Timepoint",
  all = TRUE
)


vdz_timepoint_check[
  is.na(ExpectedN),
  ExpectedN := 0L
]

vdz_timepoint_check[
  is.na(ObservedN),
  ObservedN := 0L
]


vdz_timepoint_check[
  ,
  Match :=
    ExpectedN ==
      ObservedN
]


check(
  all(
    vdz_timepoint_check$Match
  ),
  "VDZ longitudinal timepoint counts reproduce frozen Stage-0B structure"
)


# ============================================================
# 11. Validate VDZ baseline-44 source
# ============================================================

check(
  nrow(vdz_baseline) == 44L,
  "Frozen VDZ baseline table contains 44 samples"
)


check(
  uniqueN(
    vdz_baseline$GSM
  ) == 44L,
  "VDZ baseline GSM IDs are unique"
)


master_vdz_w0 <- meta[
  Branch == "VDZ_TRIAL" &
    Timepoint == "W0"
]


check(
  nrow(master_vdz_w0) == 44L,
  "Master metadata contains exactly 44 VDZ W0 samples"
)


check(
  setequal(
    master_vdz_w0$GSM,
    vdz_baseline$GSM
  ),
  "VDZ W0 GSM universe matches frozen 44-sample baseline table"
)


# ============================================================
# 12. Validate patient-level VDZ reconstruction
# ============================================================

check(
  nrow(vdz_patient) == 44L,
  "VDZ patient reconstruction contains 44 patients"
)


check(
  uniqueN(
    vdz_patient$PatientKey
  ) == 44L,
  "VDZ patient keys are unique"
)


check(
  uniqueN(
    vdz_patient$BaselineGSM
  ) == 44L,
  "VDZ patient baseline GSM IDs are unique"
)


check(
  setequal(
    vdz_patient$BaselineGSM,
    vdz_baseline$GSM
  ),
  "VDZ patient BaselineGSM matches frozen baseline-44 sample universe"
)


check(
  vdz_patient[
    ExposureClass == "VDZ_exposed",
    .N
  ] == 41L,
  "41 VDZ-exposed patients confirmed"
)


check(
  vdz_patient[
    ExposureClass == "Placebo_only",
    .N
  ] == 3L,
  "3 placebo-only patients confirmed"
)


# ============================================================
# 13. Cross-check exposed-41 frozen table
# ============================================================

check(
  nrow(vdz_exposed) == 41L,
  "Frozen exposed-baseline candidate table contains 41 patients"
)


check(
  uniqueN(
    vdz_exposed$PatientKey
  ) == 41L,
  "Exposed-baseline PatientKeys are unique"
)


check(
  all(
    vdz_exposed$ExposureClass ==
      "VDZ_exposed"
  ),
  "All exposed-baseline candidates are VDZ_exposed"
)


check(
  setequal(
    vdz_exposed$BaselineGSM,
    vdz_patient[
      ExposureClass == "VDZ_exposed",
      BaselineGSM
    ]
  ),
  "Exposed-41 BaselineGSM universe matches patient reconstruction"
)


# Direct patient-level response concordance
p1 <- copy(
  vdz_exposed[
    order(PatientKey)
  ]
)

p2 <- copy(
  vdz_patient[
    ExposureClass == "VDZ_exposed"
  ][
    order(PatientKey)
  ]
)


check(
  identical(
    p1$PatientKey,
    p2$PatientKey
  ),
  "Exposed-41 PatientKey order aligns to patient reconstruction"
)


for (cc in c(
  "BaselineGSM",
  "BaselineTreatmentSchedule",
  "W6Response",
  "W12Response",
  "W52Response"
)) {

  check(
    same_with_na(
      p1[[cc]],
      p2[[cc]]
    ),
    paste0(
      "Exposed-41 ",
      cc,
      " matches patient reconstruction"
    )
  )
}


# ============================================================
# 14. Frozen endpoint counts
# ============================================================

vdz_exposed_patients <- vdz_patient[
  ExposureClass == "VDZ_exposed"
]


check_response_counts(
  dt = vdz_exposed_patients,
  column = "W6Response",
  expected_total = 27L,
  expected_response = 6L,
  expected_nonresponse = 21L,
  label = "W6 primary endpoint"
)


check_response_counts(
  dt = vdz_exposed_patients,
  column = "W12Response",
  expected_total = 13L,
  expected_response = 3L,
  expected_nonresponse = 10L,
  label = "W12 exploratory endpoint"
)


check_response_counts(
  dt = vdz_exposed_patients,
  column = "W52Response",
  expected_total = 19L,
  expected_response = 13L,
  expected_nonresponse = 6L,
  label = "W52 secondary endpoint"
)


# ============================================================
# 15. Cross-check endpoint availability artifact
# ============================================================

endpoint_expected <- data.table(
  Endpoint = c(
    "W6",
    "W12",
    "W52"
  ),

  PatientsWithInterpretableResponse = c(
    27L,
    13L,
    19L
  ),

  Response = c(
    6L,
    3L,
    13L
  ),

  NonResponse = c(
    21L,
    10L,
    6L
  )
)


endpoint_check <- merge(
  endpoint_expected,
  endpoint_availability,
  by = "Endpoint",
  suffixes = c(
    "_Expected",
    "_Observed"
  )
)


check(
  nrow(endpoint_check) == 3L,
  "Endpoint availability artifact contains W6/W12/W52"
)


check(
  all(
    endpoint_check$PatientsWithInterpretableResponse_Expected ==
      endpoint_check$PatientsWithInterpretableResponse_Observed
  ),
  "Endpoint interpretable-patient counts match frozen artifact"
)


check(
  all(
    endpoint_check$Response_Expected ==
      endpoint_check$Response_Observed
  ),
  "Endpoint Response counts match frozen artifact"
)


check(
  all(
    endpoint_check$NonResponse_Expected ==
      endpoint_check$NonResponse_Observed
  ),
  "Endpoint Non-response counts match frozen artifact"
)


# ============================================================
# 16. Validate baseline ExposureClass against patient table
# ============================================================

exposure_lookup <- setNames(
  vdz_patient$ExposureClass,
  vdz_patient$PatientKey
)


baseline_gsm_lookup <- setNames(
  vdz_patient$BaselineGSM,
  vdz_patient$PatientKey
)


baseline_schedule_lookup <- setNames(
  vdz_patient$BaselineTreatmentSchedule,
  vdz_patient$PatientKey
)


vdz_w0_idx <- which(
  meta$Branch == "VDZ_TRIAL" &
    meta$Timepoint == "W0"
)


check(
  length(vdz_w0_idx) == 44L,
  "44 VDZ baseline rows located in canonical metadata"
)


check(
  same_with_na(
    meta$ExposureClass[
      vdz_w0_idx
    ],
    unname(
      exposure_lookup[
        meta$PatientKey[
          vdz_w0_idx
        ]
      ]
    )
  ),
  "VDZ baseline ExposureClass matches patient reconstruction"
)


check(
  identical(
    meta$GSM[
      vdz_w0_idx
    ],
    unname(
      baseline_gsm_lookup[
        meta$PatientKey[
          vdz_w0_idx
        ]
      ]
    )
  ),
  "VDZ baseline GSM matches patient reconstruction patient-by-patient"
)


# ============================================================
# 17. Add PRE-QC analysis-role flags
#
# These flags do NOT represent Stage-5 QC retention.
# ============================================================

meta[
  ,
  PreQCRecoveryCandidate :=
    TRUE
]


meta[
  ,
  IsControlSample :=
    Branch == "CONTROL"
]


meta[
  ,
  IsIFXSample :=
    Branch == "IFX"
]


meta[
  ,
  IsVDZTrialSample :=
    Branch == "VDZ_TRIAL"
]


meta[
  ,
  IsIFXBaseline :=
    Branch == "IFX" &
    Timepoint == "W0"
]


meta[
  ,
  IsVDZBaseline :=
    Branch == "VDZ_TRIAL" &
    Timepoint == "W0"
]


meta[
  ,
  IsVDZExposedBaseline :=
    IsVDZBaseline &
    ExposureClass == "VDZ_exposed"
]


meta[
  ,
  IsVDZPlaceboOnlyBaseline :=
    IsVDZBaseline &
    ExposureClass == "Placebo_only"
]


check(
  sum(
    meta$IsVDZBaseline
  ) == 44L,
  "44 VDZ baseline predictor samples flagged"
)


check(
  sum(
    meta$IsVDZExposedBaseline
  ) == 41L,
  "41 VDZ-exposed baseline predictor samples flagged"
)


check(
  sum(
    meta$IsVDZPlaceboOnlyBaseline
  ) == 3L,
  "3 placebo-only baseline samples flagged"
)


check(
  setequal(
    meta[
      IsVDZExposedBaseline == TRUE,
      GSM
    ],
    vdz_exposed$BaselineGSM
  ),
  "VDZ-exposed baseline flags reproduce frozen 41-sample universe"
)


# ============================================================
# 18. Attach patient-level endpoints ONLY to exposed W0
#     predictor samples
# ============================================================

w6_lookup <- setNames(
  vdz_patient$W6Response,
  vdz_patient$PatientKey
)

w12_lookup <- setNames(
  vdz_patient$W12Response,
  vdz_patient$PatientKey
)

w52_lookup <- setNames(
  vdz_patient$W52Response,
  vdz_patient$PatientKey
)


meta[
  ,
  W6EndoscopicHealing :=
    NA_character_
]

meta[
  ,
  W12EndoscopicHealing :=
    NA_character_
]

meta[
  ,
  W52EndoscopicHealing :=
    NA_character_
]


exposed_baseline_idx <- which(
  meta$IsVDZExposedBaseline
)


meta[
  exposed_baseline_idx,
  W6EndoscopicHealing :=
    unname(
      w6_lookup[
        PatientKey
      ]
    )
]


meta[
  exposed_baseline_idx,
  W12EndoscopicHealing :=
    unname(
      w12_lookup[
        PatientKey
      ]
    )
]


meta[
  exposed_baseline_idx,
  W52EndoscopicHealing :=
    unname(
      w52_lookup[
        PatientKey
      ]
    )
]


# ============================================================
# 19. Binary endpoint encoding
#
# Response     = 1
# Non-response = 0
# Missing      = NA
# ============================================================

meta[
  ,
  W6EndoscopicHealingBinary :=
    response_binary(
      W6EndoscopicHealing
    )
]


meta[
  ,
  W12EndoscopicHealingBinary :=
    response_binary(
      W12EndoscopicHealing
    )
]


meta[
  ,
  W52EndoscopicHealingBinary :=
    response_binary(
      W52EndoscopicHealing
    )
]


# ============================================================
# 20. Endpoint eligibility flags
# ============================================================

meta[
  ,
  PrimaryW6Eligible :=
    IsVDZExposedBaseline &
    !is.na(
      W6EndoscopicHealing
    )
]


meta[
  ,
  ExploratoryW12Eligible :=
    IsVDZExposedBaseline &
    !is.na(
      W12EndoscopicHealing
    )
]


meta[
  ,
  SecondaryW52Eligible :=
    IsVDZExposedBaseline &
    !is.na(
      W52EndoscopicHealing
    )
]


# ============================================================
# 21. Validate canonical endpoint cohorts
# ============================================================

primary_w6 <- meta[
  PrimaryW6Eligible == TRUE
]


exploratory_w12 <- meta[
  ExploratoryW12Eligible == TRUE
]


secondary_w52 <- meta[
  SecondaryW52Eligible == TRUE
]


check_response_counts(
  dt = primary_w6,
  column = "W6EndoscopicHealing",
  expected_total = 27L,
  expected_response = 6L,
  expected_nonresponse = 21L,
  label = "Canonical W6 primary baseline cohort"
)


check_response_counts(
  dt = exploratory_w12,
  column = "W12EndoscopicHealing",
  expected_total = 13L,
  expected_response = 3L,
  expected_nonresponse = 10L,
  label = "Canonical W12 exploratory baseline cohort"
)


check_response_counts(
  dt = secondary_w52,
  column = "W52EndoscopicHealing",
  expected_total = 19L,
  expected_response = 13L,
  expected_nonresponse = 6L,
  label = "Canonical W52 secondary baseline cohort"
)


check(
  all(
    primary_w6$Timepoint ==
      "W0"
  ),
  "All W6 primary predictors are baseline W0 samples"
)


check(
  all(
    secondary_w52$Timepoint ==
      "W0"
  ),
  "All W52 secondary predictors are baseline W0 samples"
)


check(
  all(
    exploratory_w12$Timepoint ==
      "W0"
  ),
  "All W12 exploratory predictors are baseline W0 samples"
)


check(
  all(
    primary_w6$ExposureClass ==
      "VDZ_exposed"
  ),
  "All W6 primary predictors are VDZ-exposed"
)


check(
  all(
    secondary_w52$ExposureClass ==
      "VDZ_exposed"
  ),
  "All W52 secondary predictors are VDZ-exposed"
)


# ============================================================
# 22. Ensure endpoints are NOT attached elsewhere
# ============================================================

check(
  meta[
    IsVDZExposedBaseline == FALSE &
      !is.na(
        W6EndoscopicHealing
      ),
    .N
  ] == 0L,
  "W6 endpoint is attached only to exposed baseline predictors"
)


check(
  meta[
    IsVDZExposedBaseline == FALSE &
      !is.na(
        W12EndoscopicHealing
      ),
    .N
  ] == 0L,
  "W12 endpoint is attached only to exposed baseline predictors"
)


check(
  meta[
    IsVDZExposedBaseline == FALSE &
      !is.na(
        W52EndoscopicHealing
      ),
    .N
  ] == 0L,
  "W52 endpoint is attached only to exposed baseline predictors"
)


# ============================================================
# 23. Analysis-role label
#
# Descriptive role only.
# Multiple endpoint eligibility remains represented by
# dedicated logical flags.
# ============================================================

meta[
  ,
  AnalysisRole :=
    fcase(

      IsVDZExposedBaseline,
      "VDZ_EXPOSED_BASELINE_PREDICTOR",

      IsVDZPlaceboOnlyBaseline,
      "VDZ_PLACEBO_ONLY_BASELINE",

      Branch == "VDZ_TRIAL",
      "VDZ_LONGITUDINAL_FOLLOWUP",

      IsIFXBaseline,
      "IFX_BASELINE",

      Branch == "IFX",
      "IFX_FOLLOWUP",

      Branch == "CONTROL",
      "NON_IBD_CONTROL",

      default =
        "OTHER"
    )
]


# ============================================================
# 24. Core analysis metadata
# ============================================================

analysis_meta <- meta[
  ,
  .(
    SampleOrder,
    GSM,
    GEOAccession,
    Title,

    PatientNumber,
    PatientKey,

    Branch,
    Disease,
    Timepoint,
    TreatmentSchedule,

    MayoEndoscopicSubscore,
    ResponseAtVisit,

    ExposureClass,

    PreQCRecoveryCandidate,

    IsControlSample,
    IsIFXSample,
    IsVDZTrialSample,
    IsIFXBaseline,
    IsVDZBaseline,
    IsVDZExposedBaseline,
    IsVDZPlaceboOnlyBaseline,

    W6EndoscopicHealing,
    W6EndoscopicHealingBinary,

    W12EndoscopicHealing,
    W12EndoscopicHealingBinary,

    W52EndoscopicHealing,
    W52EndoscopicHealingBinary,

    PrimaryW6Eligible,
    ExploratoryW12Eligible,
    SecondaryW52Eligible,

    AnalysisRole
  )
]


check(
  identical(
    analysis_meta$GSM,
    expr_samples
  ),
  "Analysis metadata sample order matches canonical expression"
)


check(
  nrow(analysis_meta) == 178L,
  "Analysis metadata contains 178 samples"
)


# ============================================================
# 25. Endpoint definition table
# ============================================================

endpoint_definitions <- data.table(
  Endpoint = c(
    "W6",
    "W12",
    "W52"
  ),

  AnalysisRole = c(
    "PRIMARY",
    "EXPLORATORY",
    "SECONDARY"
  ),

  PredictorSample = c(
    "VDZ-exposed baseline W0 colonic biopsy",
    "VDZ-exposed baseline W0 colonic biopsy",
    "VDZ-exposed baseline W0 colonic biopsy"
  ),

  Outcome = c(
    "Endoscopic healing at W6: Mayo endoscopic subscore 0-1",
    "Endoscopic healing at W12: Mayo endoscopic subscore 0-1",
    "Endoscopic healing at W52: Mayo endoscopic subscore 0-1"
  ),

  EvaluableN = c(
    27L,
    13L,
    19L
  ),

  ResponseN = c(
    6L,
    3L,
    13L
  ),

  NonResponseN = c(
    21L,
    10L,
    6L
  )
)


# ============================================================
# 26. Cohort summaries
# ============================================================

branch_summary <- meta[
  ,
  .(
    Samples =
      .N,

    Patients =
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


role_summary <- analysis_meta[
  ,
  .N,
  by = AnalysisRole
][
  order(
    -N
  )
]


stage_summary <- data.table(
  Metric = c(
    "Canonical_expression_samples",
    "Canonical_metadata_samples",
    "PreQC_recovery_candidates",
    "CONTROL_samples",
    "IFX_samples",
    "VDZ_TRIAL_samples",
    "VDZ_baseline_samples",
    "VDZ_exposed_baseline_samples",
    "VDZ_placebo_only_baseline_samples",
    "Primary_W6_evaluable",
    "Primary_W6_response",
    "Primary_W6_nonresponse",
    "Exploratory_W12_evaluable",
    "Exploratory_W12_response",
    "Exploratory_W12_nonresponse",
    "Secondary_W52_evaluable",
    "Secondary_W52_response",
    "Secondary_W52_nonresponse"
  ),

  N = c(
    178L,
    nrow(meta),
    sum(meta$PreQCRecoveryCandidate),

    meta[
      Branch == "CONTROL",
      .N
    ],

    meta[
      Branch == "IFX",
      .N
    ],

    meta[
      Branch == "VDZ_TRIAL",
      .N
    ],

    sum(
      meta$IsVDZBaseline
    ),

    sum(
      meta$IsVDZExposedBaseline
    ),

    sum(
      meta$IsVDZPlaceboOnlyBaseline
    ),

    nrow(primary_w6),

    primary_w6[
      W6EndoscopicHealing ==
        "Response",
      .N
    ],

    primary_w6[
      W6EndoscopicHealing ==
        "Non-response",
      .N
    ],

    nrow(exploratory_w12),

    exploratory_w12[
      W12EndoscopicHealing ==
        "Response",
      .N
    ],

    exploratory_w12[
      W12EndoscopicHealing ==
        "Non-response",
      .N
    ],

    nrow(secondary_w52),

    secondary_w52[
      W52EndoscopicHealing ==
        "Response",
      .N
    ],

    secondary_w52[
      W52EndoscopicHealing ==
        "Non-response",
      .N
    ]
  )
)


cat("\n============================================================\n")
cat("STAGE 4 SUMMARY\n")
cat("============================================================\n")

print(
  stage_summary
)


cat("\n============================================================\n")
cat("BRANCH SUMMARY\n")
cat("============================================================\n")

print(
  branch_summary
)


cat("\n============================================================\n")
cat("ANALYSIS ROLE SUMMARY\n")
cat("============================================================\n")

print(
  role_summary
)


cat("\n============================================================\n")
cat("FROZEN ENDPOINT DEFINITIONS\n")
cat("============================================================\n")

print(
  endpoint_definitions
)


# ============================================================
# 27. Save artifacts
# ============================================================

canonical_tsv <- file.path(
  out_dir,
  "GSE73661_canonical_metadata_178_samples.tsv"
)

canonical_rds <- file.path(
  out_dir,
  "GSE73661_canonical_metadata_178_samples.rds"
)


fwrite(
  meta,
  canonical_tsv,
  sep = "\t",
  quote = FALSE,
  na = ""
)


saveRDS(
  meta,
  canonical_rds,
  compress = FALSE
)


fwrite(
  analysis_meta,
  file.path(
    out_dir,
    "GSE73661_analysis_metadata_178_samples.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  meta[
    IsVDZExposedBaseline == TRUE
  ],
  file.path(
    out_dir,
    "GSE73661_vdz_exposed_baseline_41.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  primary_w6,
  file.path(
    out_dir,
    "GSE73661_primary_W6_baseline_predictor_27.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  exploratory_w12,
  file.path(
    out_dir,
    "GSE73661_exploratory_W12_baseline_predictor_13.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  secondary_w52,
  file.path(
    out_dir,
    "GSE73661_secondary_W52_baseline_predictor_19.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  endpoint_definitions,
  file.path(
    out_dir,
    "GSE73661_endpoint_definitions_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  vdz_timepoint_check,
  file.path(
    out_dir,
    "GSE73661_stage4_vdz_timepoint_validation.tsv"
  ),
  sep = "\t"
)


fwrite(
  branch_summary,
  file.path(
    out_dir,
    "GSE73661_stage4_branch_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  role_summary,
  file.path(
    out_dir,
    "GSE73661_stage4_analysis_role_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  stage_summary,
  file.path(
    out_dir,
    "GSE73661_stage4_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 28. RDS round-trip
# ============================================================

roundtrip <- readRDS(
  canonical_rds
)


check(
  identical(
    meta,
    roundtrip
  ),
  "Canonical metadata RDS round-trip is exact"
)


# ============================================================
# 29. Expression remains untouched
# ============================================================

check(
  identical(
    dim(expr),
    c(
      18407L,
      178L
    )
  ),
  "Canonical expression remains 18,407 x 178"
)


check(
  identical(
    colnames(expr),
    meta$GSM
  ),
  "Final metadata remains exactly aligned to expression columns"
)


# ============================================================
# 30. Checksums
# ============================================================

hash_files <- c(
  expr_file,
  master_file,
  vdz_patient_file,
  canonical_tsv,
  canonical_rds
)


hashes <- tools::md5sum(
  hash_files
)


fwrite(
  data.table(
    File =
      names(hashes),

    MD5 =
      as.character(
        unname(hashes)
      )
  ),
  file.path(
    out_dir,
    "GSE73661_stage4_artifact_md5.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 31. FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] Canonical metadata contains all 178 samples\n")
cat("[PASS] Metadata order exactly matches canonical expression columns\n")
cat("[PASS] CONTROL / IFX / VDZ branch structure reproduced\n")
cat("[PASS] VDZ longitudinal timepoint structure reproduced\n")
cat("[PASS] Frozen VDZ baseline-44 universe reproduced\n")
cat("[PASS] 41 VDZ-exposed baseline predictors reproduced\n")
cat("[PASS] 3 placebo-only baseline samples reproduced\n")
cat("[PASS] W6 primary endpoint = 27 evaluable: 6 Response / 21 Non-response\n")
cat("[PASS] W12 exploratory endpoint = 13 evaluable: 3 Response / 10 Non-response\n")
cat("[PASS] W52 secondary endpoint = 19 evaluable: 13 Response / 6 Non-response\n")
cat("[PASS] Clinical endpoints attached only to exposed W0 predictor samples\n")
cat("[PASS] All 178 samples remain PRE-QC recovery candidates\n")
cat("[PASS] No sample removed\n")
cat("[PASS] No expression values modified\n")
cat("[PASS] Canonical metadata RDS round-trip exact\n")

cat("\nFINAL STATUS:\n")
cat("PASS_GSE73661_STAGE4_CANONICAL_METADATA_FINALIZATION\n")
cat("============================================================\n")
