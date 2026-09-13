library(rhdf5)
library(data.table)

# ============================================================
# GSE282122 reference census
#
# Purpose:
#   1. Count cells by Patient × final_analysis
#   2. Apply validated reference mapping
#   3. Generate Patient × KEEP identity census
#   4. Generate candidate census (Glial / Tuft)
#   5. Audit REVIEW / EXCLUDE labels
#   6. Summarize donor coverage for later balanced sampling
#
# IMPORTANT:
#   - Does NOT load expression matrix X
#   - Uses categorical integer codes only
# ============================================================


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_root <- "/home/mazekai/IBD_EcoTyper"

dataset_dir <- file.path(
  project_root,
  "03_reference",
  "GSE282122"
)

build_dir <- file.path(
  dataset_dir,
  "processed",
  "reference_build_v1"
)

mapping_dir <- file.path(
  build_dir,
  "00_mapping"
)

output_dir <- file.path(
  build_dir,
  "01_census"
)

log_dir <- file.path(
  build_dir,
  "logs"
)


# ------------------------------------------------------------
# Raw H5AD
# ------------------------------------------------------------

h5ad_file <- paste0(
  "/home/mazekai/enteric_glia/data_scRNA/GSE282122/",
  "TAURUS_raw_counts_annotated_final.h5ad"
)


# ------------------------------------------------------------
# Validated mapping
# ------------------------------------------------------------

mapping_file <- file.path(
  mapping_dir,
  "GSE282122_celltype_mapping.tsv"
)


# ------------------------------------------------------------
# Create directories if absent
# ------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Check required input files
# ------------------------------------------------------------

if (!file.exists(h5ad_file)) {
  stop(
    "ERROR: H5AD file not found:\n",
    h5ad_file
  )
}

if (!file.exists(mapping_file)) {
  stop(
    "ERROR: validated mapping file not found:\n",
    mapping_file
  )
}

cat("Project root : ", project_root, "\n", sep = "")
cat("Dataset dir  : ", dataset_dir, "\n", sep = "")
cat("Build dir    : ", build_dir, "\n", sep = "")
cat("Mapping file : ", mapping_file, "\n", sep = "")
cat("Output dir   : ", output_dir, "\n\n", sep = "")

# ------------------------------------------------------------
# 2. Frozen KEEP ontology
# ------------------------------------------------------------

keep_identities <- c(
  "B",
  "Plasma",
  "CD4_T",
  "CD8_T",
  "NK",
  "ILC",
  "Monocyte_Macrophage",
  "DC",
  "Mast",
  "Fibroblast",
  "Pericyte",
  "Endothelial",
  "Ileal_absorptive",
  "Colonic_absorptive",
  "Goblet",
  "Stem_TA_progenitor",
  "Paneth",
  "Enteroendocrine"
)

candidate_identities <- c(
  "Glial",
  "Tuft"
)


# ------------------------------------------------------------
# 3. Read validated mapping
# ------------------------------------------------------------

mapping <- fread(
  mapping_file,
  na.strings = c("NA", "")
)

cat("========================================\n")
cat("GSE282122 REFERENCE CENSUS\n")
cat("========================================\n\n")

cat("Mapping rows:", nrow(mapping), "\n")


# ------------------------------------------------------------
# 4. Mapping sanity checks
# ------------------------------------------------------------

if (nrow(mapping) != 109L) {
  stop(
    "ERROR: mapping file must contain exactly 109 rows."
  )
}

if (anyDuplicated(mapping$final_analysis)) {
  stop(
    "ERROR: duplicated final_analysis labels in mapping file."
  )
}

actual_keep <- sort(
  unique(
    mapping[
      status == "KEEP",
      reference_identity
    ]
  )
)

if (!setequal(actual_keep, keep_identities)) {

  cat("\nExpected KEEP identities:\n")
  print(sort(keep_identities))

  cat("\nObserved KEEP identities:\n")
  print(actual_keep)

  stop(
    "ERROR: KEEP ontology differs from frozen 18-class ontology."
  )
}

actual_candidate <- sort(
  unique(
    mapping[
      status == "CANDIDATE",
      reference_identity
    ]
  )
)

if (!setequal(actual_candidate, candidate_identities)) {
  stop(
    "ERROR: candidate identities differ from Glial/Tuft."
  )
}

cat("PASS: mapping ontology validated\n")


# ------------------------------------------------------------
# 5. Read H5AD categorical metadata only
# ------------------------------------------------------------

cat("\nReading categorical metadata from H5AD...\n")

patient_categories <- h5read(
  h5ad_file,
  "/obs/Patient/categories"
)

patient_codes <- h5read(
  h5ad_file,
  "/obs/Patient/codes"
)

fa_categories <- h5read(
  h5ad_file,
  "/obs/final_analysis/categories"
)

fa_codes <- h5read(
  h5ad_file,
  "/obs/final_analysis/codes"
)


# ------------------------------------------------------------
# 6. Basic H5AD checks
# ------------------------------------------------------------

n_cells <- length(patient_codes)

cat("Cells:", n_cells, "\n")
cat("Patients:", length(patient_categories), "\n")
cat("final_analysis categories:", length(fa_categories), "\n")

if (length(fa_codes) != n_cells) {
  stop(
    "ERROR: Patient and final_analysis code vectors differ in length."
  )
}

if (n_cells != 987743L) {
  stop(
    "ERROR: unexpected number of cells: ",
    n_cells
  )
}

if (length(fa_categories) != 109L) {
  stop(
    "ERROR: expected 109 final_analysis categories."
  )
}

if (any(patient_codes < 0L)) {
  stop(
    "ERROR: missing Patient values detected."
  )
}

if (any(fa_codes < 0L)) {
  stop(
    "ERROR: missing final_analysis values detected."
  )
}

cat("PASS: H5AD metadata dimensions validated\n")


# ------------------------------------------------------------
# 7. Validate mapping labels against H5AD again
# ------------------------------------------------------------

missing_mapping <- setdiff(
  fa_categories,
  mapping$final_analysis
)

extra_mapping <- setdiff(
  mapping$final_analysis,
  fa_categories
)

if (length(missing_mapping) > 0L) {

  cat("\nMissing mapping labels:\n")
  print(missing_mapping)

  stop("ERROR: incomplete mapping.")
}

if (length(extra_mapping) > 0L) {

  cat("\nMapping labels absent from H5AD:\n")
  print(extra_mapping)

  stop("ERROR: invalid mapping labels.")
}

cat("PASS: H5AD labels match mapping 109/109\n")


# ------------------------------------------------------------
# 8. Validate H5AD categorical codes against mapping code column
# ------------------------------------------------------------

expected_codes <- match(
  mapping$final_analysis,
  fa_categories
) - 1L

if (!all(mapping$code == expected_codes)) {

  bad <- mapping[
    code != expected_codes,
    .(
      final_analysis,
      mapping_code = code,
      expected_code = expected_codes
    )
  ]

  print(bad)

  stop(
    "ERROR: mapping code values do not match H5AD categorical codes."
  )
}

cat("PASS: mapping codes match H5AD categorical codes\n")


# ------------------------------------------------------------
# 9. Count Patient × final_analysis using integer codes
# ------------------------------------------------------------

cat("\nCounting Patient × final_analysis...\n")

cell_codes <- data.table(
  patient_code = as.integer(patient_codes),
  final_code = as.integer(fa_codes)
)

pair_counts <- cell_codes[
  ,
  .(n_cells = .N),
  by = .(
    patient_code,
    final_code
  )
]

rm(cell_codes)
gc()


# ------------------------------------------------------------
# 10. Decode categorical codes
# ------------------------------------------------------------

pair_counts[
  ,
  Patient := patient_categories[
    patient_code + 1L
  ]
]

pair_counts[
  ,
  final_analysis := fa_categories[
    final_code + 1L
  ]
]


# ------------------------------------------------------------
# 11. Check source-level counts against mapping template
# ------------------------------------------------------------

source_counts <- pair_counts[
  ,
  .(
    counted_cells = sum(n_cells)
  ),
  by = final_analysis
]

source_check <- merge(
  mapping[
    ,
    .(
      final_analysis,
      expected_cells = n_cells
    )
  ],
  source_counts,
  by = "final_analysis",
  all = TRUE
)

bad_counts <- source_check[
  expected_cells != counted_cells
]

if (nrow(bad_counts) > 0L) {

  print(bad_counts)

  stop(
    "ERROR: source-label cell counts differ from mapping template."
  )
}

cat(
  "PASS: all final_analysis counts reproduce original template\n"
)


# ------------------------------------------------------------
# 12. Attach mapping
# ------------------------------------------------------------

pair_mapped <- merge(
  pair_counts,
  mapping[
    ,
    .(
      final_analysis,
      reference_identity,
      status,
      action
    )
  ],
  by = "final_analysis",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(pair_mapped$status)) {

  print(
    pair_mapped[
      is.na(status)
    ]
  )

  stop(
    "ERROR: unmapped final_analysis labels remain."
  )
}


# ------------------------------------------------------------
# 13. Global accounting check
# ------------------------------------------------------------

total_counted <- sum(pair_mapped$n_cells)

if (total_counted != n_cells) {
  stop(
    "ERROR: total cell accounting failed."
  )
}

cat(
  "PASS: all ",
  format(n_cells, big.mark = ","),
  " cells accounted for\n",
  sep = ""
)


# ============================================================
# 14. KEEP census
# ============================================================

keep_long <- pair_mapped[
  status == "KEEP",
  .(
    n_cells = sum(n_cells)
  ),
  by = .(
    Patient,
    reference_identity
  )
]


# ------------------------------------------------------------
# 15. Complete KEEP grid
#     Include explicit zeroes for absent patient × identity pairs
# ------------------------------------------------------------

keep_complete <- CJ(
  Patient = patient_categories,
  reference_identity = keep_identities,
  unique = TRUE
)

keep_complete <- merge(
  keep_complete,
  keep_long,
  by = c(
    "Patient",
    "reference_identity"
  ),
  all.x = TRUE,
  sort = FALSE
)

keep_complete[
  is.na(n_cells),
  n_cells := 0L
]

keep_complete[
  ,
  Patient := factor(
    Patient,
    levels = patient_categories
  )
]

keep_complete[
  ,
  reference_identity := factor(
    reference_identity,
    levels = keep_identities
  )
]

setorder(
  keep_complete,
  Patient,
  reference_identity
)

keep_complete[
  ,
  Patient := as.character(Patient)
]

keep_complete[
  ,
  reference_identity := as.character(reference_identity)
]


# ------------------------------------------------------------
# 16. Wide Patient × CellType matrix
# ------------------------------------------------------------

keep_wide <- dcast(
  keep_complete,
  Patient ~ reference_identity,
  value.var = "n_cells"
)


# ============================================================
# 17. Donor coverage statistics
# ============================================================

safe_quantile <- function(x, prob) {

  x <- x[x > 0]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    quantile(
      x,
      probs = prob,
      names = FALSE,
      type = 7
    )
  )
}


keep_donor_summary <- keep_complete[
  ,
  {

    x <- n_cells
    xpos <- x[x > 0]

    list(
      total_cells = sum(x),

      n_patients_total = .N,

      n_patients_with_cells = sum(x > 0),

      pct_patients_with_cells =
        100 * sum(x > 0) / .N,

     min_nonzero =
        if (length(xpos) > 0L)
          as.numeric(min(xpos))
        else
          NA_real_,

      q25_nonzero =
        safe_quantile(x, 0.25),

      median_nonzero =
        safe_quantile(x, 0.50),

      q75_nonzero =
        safe_quantile(x, 0.75),

      max_nonzero =
        if (length(xpos) > 0L)
          as.numeric(max(xpos))
        else
          NA_real_,

      mean_nonzero =
        if (length(xpos) > 0)
          mean(xpos)
        else
          NA_real_
    )
  },
  by = reference_identity
]

setorder(
  keep_donor_summary,
  -total_cells
)


# ============================================================
# 18. Candidate census
# ============================================================

candidate_long <- pair_mapped[
  status == "CANDIDATE",
  .(
    n_cells = sum(n_cells)
  ),
  by = .(
    Patient,
    reference_identity
  )
]


candidate_complete <- CJ(
  Patient = patient_categories,
  reference_identity = candidate_identities,
  unique = TRUE
)

candidate_complete <- merge(
  candidate_complete,
  candidate_long,
  by = c(
    "Patient",
    "reference_identity"
  ),
  all.x = TRUE,
  sort = FALSE
)

candidate_complete[
  is.na(n_cells),
  n_cells := 0L
]


candidate_donor_summary <- candidate_complete[
  ,
  {

    x <- n_cells
    xpos <- x[x > 0]

    list(
      total_cells = as.integer(sum(x)),

      n_patients_total = as.integer(.N),

      n_patients_with_cells =
        as.integer(sum(x > 0)),

      pct_patients_with_cells =
        as.numeric(
          100 * sum(x > 0) / .N
        ),

      min_nonzero =
        if (length(xpos) > 0L)
          as.numeric(min(xpos))
        else
          NA_real_,

      median_nonzero =
        if (length(xpos) > 0L)
          as.numeric(median(xpos))
        else
          NA_real_,

      max_nonzero =
        if (length(xpos) > 0L)
          as.numeric(max(xpos))
        else
          NA_real_
    )
  },
  by = reference_identity
]

# ============================================================
# 19. REVIEW / EXCLUDE audit
# ============================================================

review_exclude <- pair_mapped[
  status %in% c(
    "REVIEW",
    "EXCLUDE"
  ),
  .(
    n_cells = sum(n_cells)
  ),
  by = .(
    Patient,
    final_analysis,
    status,
    action
  )
]

setorder(
  review_exclude,
  status,
  final_analysis,
  Patient
)


# ============================================================
# 20. Dataset-level identity summary
# ============================================================

identity_summary <- pair_mapped[
  ,
  .(
    n_cells = sum(n_cells)
  ),
  by = .(
    status,
    reference_identity
  )
]

setorder(
  identity_summary,
  status,
  -n_cells
)


# ============================================================
# 21. Status summary
# ============================================================

status_summary <- pair_mapped[
  ,
  .(
    n_cells = sum(n_cells)
  ),
  by = status
]

setorder(
  status_summary,
  status
)


# ============================================================
# 22. Save outputs
# ============================================================

fwrite(
  keep_complete,
  file.path(
    output_dir,
    "GSE282122_patient_KEEP_celltype_census_long.tsv"
  ),
  sep = "\t"
)

fwrite(
  keep_wide,
  file.path(
    output_dir,
    "GSE282122_patient_KEEP_celltype_census_wide.tsv"
  ),
  sep = "\t"
)

fwrite(
  keep_donor_summary,
  file.path(
    output_dir,
    "GSE282122_KEEP_donor_coverage_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  candidate_complete,
  file.path(
    output_dir,
    "GSE282122_patient_CANDIDATE_census.tsv"
  ),
  sep = "\t"
)

fwrite(
  candidate_donor_summary,
  file.path(
    output_dir,
    "GSE282122_CANDIDATE_donor_coverage_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  review_exclude,
  file.path(
    output_dir,
    "GSE282122_REVIEW_EXCLUDE_patient_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  identity_summary,
  file.path(
    output_dir,
    "GSE282122_identity_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  status_summary,
  file.path(
    output_dir,
    "GSE282122_status_summary.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 23. Save compact R object for later sampling
# ------------------------------------------------------------

saveRDS(
  list(
    keep_census = keep_complete,
    keep_summary = keep_donor_summary,
    candidate_census = candidate_complete,
    candidate_summary = candidate_donor_summary,
    review_exclude = review_exclude,
    mapping = mapping
  ),
  file.path(
    output_dir,
    "GSE282122_reference_census.rds"
  )
)


# ============================================================
# 24. Console report
# ============================================================

cat("\n========================================\n")
cat("STATUS SUMMARY\n")
cat("========================================\n")

print(status_summary)


cat("\n========================================\n")
cat("KEEP DONOR COVERAGE SUMMARY\n")
cat("========================================\n")

print(keep_donor_summary)


cat("\n========================================\n")
cat("CANDIDATE DONOR COVERAGE SUMMARY\n")
cat("========================================\n")

print(candidate_donor_summary)


cat("\n========================================\n")
cat("OUTPUT DIRECTORY\n")
cat("========================================\n")

cat(output_dir, "\n")


cat("\n========================================\n")
cat("SUCCESS\n")
cat("========================================\n")

cat(
  "GSE282122 Patient × CellType census completed.\n"
)
