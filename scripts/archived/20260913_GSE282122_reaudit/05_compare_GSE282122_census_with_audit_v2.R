library(data.table)

# ============================================================
# GSE282122 census cross-validation
#
# Compare:
#   OLD Python audit_v2 donor-level census
#       vs
#   NEW R reference_build_v1 KEEP census
#
# IMPORTANT:
#   - Read-only comparison
#   - Does NOT modify either source
#   - Does NOT silently remap biological identities
# ============================================================


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_root <- "/home/mazekai/IBD_EcoTyper"

old_file <- file.path(
  project_root,
  "03_reference",
  "GSE282122",
  "audit_v2",
  "output",
  "04_identity_coverage_v1",
  "identity_counts_by_donor.tsv"
)

new_file <- file.path(
  project_root,
  "03_reference",
  "GSE282122",
  "processed",
  "reference_build_v1",
  "01_census",
  "GSE282122_patient_KEEP_celltype_census_long.tsv"
)

output_dir <- file.path(
  project_root,
  "03_reference",
  "GSE282122",
  "processed",
  "reference_build_v1",
  "04_qc",
  "census_crosscheck_v1"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Check inputs
# ------------------------------------------------------------

if (!file.exists(old_file)) {
  stop(
    "ERROR: old Python audit file not found:\n",
    old_file
  )
}

if (!file.exists(new_file)) {
  stop(
    "ERROR: new R census file not found:\n",
    new_file
  )
}


cat("========================================\n")
cat("GSE282122 CENSUS CROSS-VALIDATION\n")
cat("========================================\n\n")

cat("OLD:", old_file, "\n")
cat("NEW:", new_file, "\n\n")


# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------

old <- fread(
  old_file,
  na.strings = c("", "NA")
)

new <- fread(
  new_file,
  na.strings = c("", "NA")
)


cat("Old rows:", nrow(old), "\n")
cat("New rows:", nrow(new), "\n\n")

cat("OLD columns:\n")
print(names(old))

cat("\nNEW columns:\n")
print(names(new))


# ------------------------------------------------------------
# 4. Helper: detect a column conservatively
# ------------------------------------------------------------

detect_column <- function(
    dt,
    candidates,
    role
) {

  hit <- candidates[
    candidates %in% names(dt)
  ]

  if (length(hit) == 0L) {

    cat(
      "\nCould not identify ",
      role,
      " column.\n",
      sep = ""
    )

    cat("Candidate names tried:\n")
    print(candidates)

    cat("\nActual columns:\n")
    print(names(dt))

    stop(
      "ERROR: column detection failed for ",
      role,
      "."
    )
  }

  if (length(hit) > 1L) {

    cat(
      "\nMultiple possible ",
      role,
      " columns detected:\n",
      sep = ""
    )

    print(hit)

    stop(
      "ERROR: ambiguous column detection for ",
      role,
      "."
    )
  }

  hit
}


# ------------------------------------------------------------
# 5. OLD Python audit columns
#    Fixed according to actual file schema
# ------------------------------------------------------------

old_donor_col <- "Patient"
old_identity_col <- "coarse_identity"
old_count_col <- "n_cells"

required_old <- c(
  old_donor_col,
  old_identity_col,
  old_count_col
)

missing_old <- setdiff(
  required_old,
  names(old)
)

if (length(missing_old) > 0L) {
  stop(
    "ERROR: OLD audit missing required columns: ",
    paste(missing_old, collapse = ", ")
  )
}

cat("\nValidated OLD columns:\n")
cat("  donor    :", old_donor_col, "\n")
cat("  identity :", old_identity_col, "\n")
cat("  count    :", old_count_col, "\n")


# ------------------------------------------------------------
# 6. Validate NEW R census columns
# ------------------------------------------------------------

required_new <- c(
  "Patient",
  "reference_identity",
  "n_cells"
)

missing_new <- setdiff(
  required_new,
  names(new)
)

if (length(missing_new) > 0L) {

  stop(
    "ERROR: new census missing columns: ",
    paste(missing_new, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 7. Standardize only structural column names
#
# NO biological identity remapping is performed here.
# ------------------------------------------------------------

old_std <- old[
  ,
  .(
    Patient = as.character(
      get(old_donor_col)
    ),

    reference_identity = as.character(
      get(old_identity_col)
    ),

    n_cells_old = as.numeric(
      get(old_count_col)
    )
  )
]

new_std <- new[
  ,
  .(
    Patient = as.character(Patient),
    reference_identity =
      as.character(reference_identity),
    n_cells_new =
      as.numeric(n_cells)
  )
]


# ------------------------------------------------------------
# 8. Basic cleanup: whitespace only
# ------------------------------------------------------------

old_std[
  ,
  Patient := trimws(Patient)
]

old_std[
  ,
  reference_identity :=
    trimws(reference_identity)
]

new_std[
  ,
  Patient := trimws(Patient)
]

new_std[
  ,
  reference_identity :=
    trimws(reference_identity)
]


# ------------------------------------------------------------
# 9. Reject missing structural fields
# ------------------------------------------------------------

bad_old <- old_std[
  is.na(Patient) |
    Patient == "" |
    is.na(reference_identity) |
    reference_identity == "" |
    is.na(n_cells_old)
]

if (nrow(bad_old) > 0L) {

  cat("\nInvalid rows in OLD audit:\n")
  print(bad_old)

  stop(
    "ERROR: OLD audit contains missing structural values."
  )
}


bad_new <- new_std[
  is.na(Patient) |
    Patient == "" |
    is.na(reference_identity) |
    reference_identity == "" |
    is.na(n_cells_new)
]

if (nrow(bad_new) > 0L) {

  cat("\nInvalid rows in NEW census:\n")
  print(bad_new)

  stop(
    "ERROR: NEW census contains missing structural values."
  )
}


# ------------------------------------------------------------
# 10. Aggregate in case OLD file has multiple rows
#     per donor × identity
# ------------------------------------------------------------

old_agg <- old_std[
  ,
  .(
    n_cells_old = sum(n_cells_old)
  ),
  by = .(
    Patient,
    reference_identity
  )
]

new_agg <- new_std[
  ,
  .(
    n_cells_new = sum(n_cells_new)
  ),
  by = .(
    Patient,
    reference_identity
  )
]


# ------------------------------------------------------------
# 11. Define frozen NEW 18-class KEEP ontology
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


new_identity_set <- sort(
  unique(new_agg$reference_identity)
)

if (!setequal(
  new_identity_set,
  keep_identities
)) {

  cat("\nExpected NEW KEEP ontology:\n")
  print(sort(keep_identities))

  cat("\nObserved NEW identities:\n")
  print(new_identity_set)

  stop(
    "ERROR: NEW census does not contain expected 18-class ontology."
  )
}

cat(
  "\nPASS: NEW census contains frozen 18-class ontology\n"
)


# ------------------------------------------------------------
# 12. Compare identity sets before cell-count comparison
# ------------------------------------------------------------

old_identity_set <- sort(
  unique(old_agg$reference_identity)
)

only_old_identity <- setdiff(
  old_identity_set,
  new_identity_set
)

only_new_identity <- setdiff(
  new_identity_set,
  old_identity_set
)


identity_set_report <- data.table(
  reference_identity = sort(
    unique(
      c(
        old_identity_set,
        new_identity_set
      )
    )
  )
)

identity_set_report[
  ,
  in_old := reference_identity %in%
    old_identity_set
]

identity_set_report[
  ,
  in_new := reference_identity %in%
    new_identity_set
]

fwrite(
  identity_set_report,
  file.path(
    output_dir,
    "identity_set_comparison.tsv"
  ),
  sep = "\t"
)


cat("\n========================================\n")
cat("IDENTITY SET CHECK\n")
cat("========================================\n")

cat(
  "OLD identities:",
  length(old_identity_set),
  "\n"
)

cat(
  "NEW identities:",
  length(new_identity_set),
  "\n"
)


if (length(only_old_identity) > 0L) {

  cat("\nIdentities only in OLD audit:\n")
  print(only_old_identity)
}

if (length(only_new_identity) > 0L) {

  cat("\nIdentities only in NEW census:\n")
  print(only_new_identity)
}


# ------------------------------------------------------------
# 13. Compare donor sets
# ------------------------------------------------------------

old_donor_set <- sort(
  unique(old_agg$Patient)
)

new_donor_set <- sort(
  unique(new_agg$Patient)
)

only_old_donor <- setdiff(
  old_donor_set,
  new_donor_set
)

only_new_donor <- setdiff(
  new_donor_set,
  old_donor_set
)


donor_set_report <- data.table(
  Patient = sort(
    unique(
      c(
        old_donor_set,
        new_donor_set
      )
    )
  )
)

donor_set_report[
  ,
  in_old := Patient %in% old_donor_set
]

donor_set_report[
  ,
  in_new := Patient %in% new_donor_set
]

fwrite(
  donor_set_report,
  file.path(
    output_dir,
    "donor_set_comparison.tsv"
  ),
  sep = "\t"
)


cat("\n========================================\n")
cat("DONOR SET CHECK\n")
cat("========================================\n")

cat(
  "OLD donors:",
  length(old_donor_set),
  "\n"
)

cat(
  "NEW donors:",
  length(new_donor_set),
  "\n"
)

if (length(only_old_donor) > 0L) {

  cat("\nDonors only in OLD audit:\n")
  print(only_old_donor)
}

if (length(only_new_donor) > 0L) {

  cat("\nDonors only in NEW census:\n")
  print(only_new_donor)
}


# ------------------------------------------------------------
# 14. Restrict comparison to identities shared exactly
#
# Important:
# no alias conversion is done automatically.
# ------------------------------------------------------------

shared_identities <- intersect(
  old_identity_set,
  new_identity_set
)

shared_donors <- intersect(
  old_donor_set,
  new_donor_set
)

cat("\nShared identities:", length(shared_identities), "\n")
cat("Shared donors:", length(shared_donors), "\n")


# ------------------------------------------------------------
# 15. Build complete donor × identity grid
# ------------------------------------------------------------

comparison_grid <- CJ(
  Patient = shared_donors,
  reference_identity = shared_identities,
  unique = TRUE
)


old_shared <- old_agg[
  Patient %in% shared_donors &
    reference_identity %in% shared_identities
]

new_shared <- new_agg[
  Patient %in% shared_donors &
    reference_identity %in% shared_identities
]


comparison <- merge(
  comparison_grid,
  old_shared,
  by = c(
    "Patient",
    "reference_identity"
  ),
  all.x = TRUE
)

comparison <- merge(
  comparison,
  new_shared,
  by = c(
    "Patient",
    "reference_identity"
  ),
  all.x = TRUE
)


# Missing donor × identity means zero cells
comparison[
  is.na(n_cells_old),
  n_cells_old := 0
]

comparison[
  is.na(n_cells_new),
  n_cells_new := 0
]


# ------------------------------------------------------------
# 16. Compute differences
# ------------------------------------------------------------

comparison[
  ,
  diff := n_cells_new - n_cells_old
]

comparison[
  ,
  abs_diff := abs(diff)
]

comparison[
  ,
  exact_match := diff == 0
]


# ------------------------------------------------------------
# 17. Identity-level comparison summary
# ------------------------------------------------------------

identity_comparison <- comparison[
  ,
  .(
    old_total = sum(n_cells_old),
    new_total = sum(n_cells_new),

    total_diff =
      sum(n_cells_new) -
      sum(n_cells_old),

    sum_abs_donor_diff =
      sum(abs_diff),

    max_abs_donor_diff =
      max(abs_diff),

    n_donors_compared = .N,

    n_donors_exact =
      sum(exact_match),

    n_donors_different =
      sum(!exact_match)
  ),
  by = reference_identity
]

identity_comparison[
  ,
  exact_identity_match :=
    (
      total_diff == 0 &
      sum_abs_donor_diff == 0
    )
]

setorder(
  identity_comparison,
  -sum_abs_donor_diff,
  reference_identity
)


# ------------------------------------------------------------
# 18. Global summary
# ------------------------------------------------------------

global_summary <- data.table(
  old_donors = length(old_donor_set),
  new_donors = length(new_donor_set),

  old_identities = length(old_identity_set),
  new_identities = length(new_identity_set),

  shared_donors = length(shared_donors),
  shared_identities = length(shared_identities),

  compared_pairs = nrow(comparison),

  exact_pairs =
    sum(comparison$exact_match),

  different_pairs =
    sum(!comparison$exact_match),

  total_abs_difference =
    sum(comparison$abs_diff),

  max_pair_difference =
    if (nrow(comparison) > 0L)
      max(comparison$abs_diff)
    else
      NA_real_
)


# ------------------------------------------------------------
# 19. Determine validation status
# ------------------------------------------------------------

same_identity_set <- setequal(
  old_identity_set,
  new_identity_set
)

same_donor_set <- setequal(
  old_donor_set,
  new_donor_set
)

all_exact <- (
  nrow(comparison) > 0L &&
  all(comparison$exact_match)
)


if (
  same_identity_set &&
  same_donor_set &&
  all_exact
) {

  validation_status <-
    "PASS_EXACT_MATCH"

} else if (
  same_identity_set &&
  same_donor_set &&
  !all_exact
) {

  validation_status <-
    "FAIL_COUNT_DIFFERENCES"

} else {

  validation_status <-
    "NOT_DIRECTLY_COMPARABLE_SET_DIFFERENCES"
}


global_summary[
  ,
  validation_status :=
    validation_status
]


# ------------------------------------------------------------
# 20. Write outputs
# ------------------------------------------------------------

fwrite(
  comparison,
  file.path(
    output_dir,
    "donor_identity_pair_comparison.tsv"
  ),
  sep = "\t"
)

fwrite(
  comparison[
    exact_match == FALSE
  ],
  file.path(
    output_dir,
    "donor_identity_differences_only.tsv"
  ),
  sep = "\t"
)

fwrite(
  identity_comparison,
  file.path(
    output_dir,
    "identity_level_comparison.tsv"
  ),
  sep = "\t"
)

fwrite(
  global_summary,
  file.path(
    output_dir,
    "global_comparison_summary.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 21. Console report
# ------------------------------------------------------------

cat("\n========================================\n")
cat("IDENTITY-LEVEL COMPARISON\n")
cat("========================================\n")

print(identity_comparison)


cat("\n========================================\n")
cat("GLOBAL SUMMARY\n")
cat("========================================\n")

print(global_summary)


cat("\n========================================\n")
cat("VALIDATION STATUS\n")
cat("========================================\n")

cat(validation_status, "\n")


cat("\nOutput directory:\n")
cat(output_dir, "\n")


# ------------------------------------------------------------
# 22. Final interpretation
# ------------------------------------------------------------

if (validation_status == "PASS_EXACT_MATCH") {

  cat(
    "\nSUCCESS: old Python audit and new R census ",
    "are exactly identical at donor × identity level.\n",
    sep = ""
  )

} else if (
  validation_status ==
    "FAIL_COUNT_DIFFERENCES"
) {

  cat(
    "\nWARNING: donor and identity sets are identical, ",
    "but cell counts differ.\n",
    sep = ""
  )

  cat(
    "Inspect donor_identity_differences_only.tsv ",
    "before proceeding to balanced sampling.\n"
  )

} else {

  cat(
    "\nNOTE: old and new pipelines use different ",
    "donor and/or identity sets.\n",
    sep = ""
  )

  cat(
    "Do NOT interpret count differences until ",
    "the set differences are understood.\n"
  )
}