#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

ROOT <- "/home/mazekai/IBD_EcoTyper"

MASS_FILE <- file.path(
  ROOT,
  "03_reference/combined_reference/output/23_ecotyper_interface_audit",
  "Z_allocated_mass_by_sample.tsv"
)

META_FILE <- file.path(
  ROOT,
  "03_reference/combined_reference/output/21_bayesprism_epithelial_qc",
  "GSE193677_metadata_used_for_stage2_qc.tsv"
)

TECH_EXCLUSION_FILE <- file.path(
  ROOT,
  "03_reference/combined_reference/output/26_panibd_panel_classification",
  "technical_exclusion_samples.tsv"
)

OUT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/26b_control_region_coverage"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(MASS_FILE),
  file.exists(META_FILE),
  file.exists(TECH_EXCLUSION_FILE)
)

# ==============================================================
# Load
# ==============================================================

mass <- read.delim(
  MASS_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- read.delim(
  META_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tech <- read.delim(
  TECH_EXCLUSION_FILE,
  stringsAsFactors = FALSE
)

required_meta <- c(
  "sample_id",
  "region_broad",
  "disease",
  "inflammation"
)

stopifnot(
  all(required_meta %in% colnames(meta)),
  "sample_id" %in% colnames(mass),
  "sample_id" %in% colnames(tech)
)

idx <- match(
  mass$sample_id,
  meta$sample_id
)

stopifnot(
  !anyNA(idx)
)

meta <- meta[
  idx,
  ,
  drop = FALSE
]

dat <- cbind(
  meta[, required_meta, drop = FALSE],
  mass[, setdiff(colnames(mass), "sample_id"), drop = FALSE]
)

# Apply same technical exclusion policy
dat <- dat[
  !(dat$sample_id %in% tech$sample_id),
  ,
  drop = FALSE
]

# ==============================================================
# Locked CORE12
# ==============================================================

core12 <- c(
  "B",
  "CD4_T",
  "CD8_T",
  "Endothelial",
  "Fibroblast",
  "Glial",
  "ILC",
  "Mast",
  "Monocyte_Macrophage",
  "Pericyte",
  "Goblet",
  "Enteroendocrine"
)

special5 <- c(
  "DC",
  "NK",
  "Plasma",
  "Colonic_absorptive",
  "Ileal_absorptive"
)

audit_types <- c(
  core12,
  special5
)

stopifnot(
  all(audit_types %in% colnames(dat))
)

# ==============================================================
# Control cohort
# ==============================================================

control <- dat[
  dat$disease == "Control" &
    dat$region_broad %in% c(
      "Colon_rectum",
      "Ileum"
    ),
  ,
  drop = FALSE
]

control$stratum <- paste(
  "Control",
  control$region_broad,
  sep = "_"
)

expected_strata <- c(
  "Control_Colon_rectum",
  "Control_Ileum"
)

cat("============================================\n")
cat("CONTROL × REGION COVERAGE AUDIT\n")
cat("============================================\n\n")

cat("Total Control samples:", nrow(control), "\n\n")

cat("Control stratum sizes:\n")
print(
  table(control$stratum)
)

stopifnot(
  all(expected_strata %in% unique(control$stratum))
)

# ==============================================================
# Summaries
# ==============================================================

rows <- list()
counter <- 1L

for (ct in audit_types) {

  for (st in expected_strata) {

    x <- control[
      control$stratum == st,
      ct,
      drop = TRUE
    ]

    rows[[counter]] <- data.frame(
      cell_type = ct,
      stratum = st,
      n = length(x),

      mass_min = min(x),
      mass_q05 = unname(
        quantile(x, 0.05)
      ),
      mass_q25 = unname(
        quantile(x, 0.25)
      ),
      mass_median = median(x),
      mass_q75 = unname(
        quantile(x, 0.75)
      ),
      mass_q95 = unname(
        quantile(x, 0.95)
      ),

      pct_ge10 =
        mean(x >= 10) * 100,

      pct_ge100 =
        mean(x >= 100) * 100,

      pct_ge1000 =
        mean(x >= 1000) * 100,

      stringsAsFactors = FALSE
    )

    counter <- counter + 1L
  }
}

summary_table <- do.call(
  rbind,
  rows
)

write.table(
  summary_table,
  file.path(
    OUT,
    "control_region_mass_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================
# CORE12 qualification
# ==============================================================

qualification <- lapply(
  core12,
  function(ct) {

    x <- summary_table[
      summary_table$cell_type == ct,
      ,
      drop = FALSE
    ]

    stopifnot(
      nrow(x) == 2
    )

    robust_each <- (
      x$mass_median >= 1000 &
        x$pct_ge100 >= 90 &
        x$pct_ge1000 >= 80
    )

    supported_each <- (
      x$mass_median >= 100 &
        x$pct_ge100 >= 80
    )

    data.frame(
      cell_type = ct,

      control_colon_median =
        x$mass_median[
          x$stratum == "Control_Colon_rectum"
        ],

      control_ileum_median =
        x$mass_median[
          x$stratum == "Control_Ileum"
        ],

      min_control_pct_ge100 =
        min(x$pct_ge100),

      min_control_pct_ge1000 =
        min(x$pct_ge1000),

      control_classification =
        if (
          all(robust_each)
        ) {
          "COMMON_ROBUST"
        } else if (
          all(supported_each)
        ) {
          "COMMON_SUPPORTED"
        } else {
          "REVIEW"
        },

      stringsAsFactors = FALSE
    )
  }
)

qualification <- do.call(
  rbind,
  qualification
)

write.table(
  qualification,
  file.path(
    OUT,
    "core12_control_qualification.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================
# Print
# ==============================================================

cat("\n============================================\n")
cat("CORE12 CONTROL QUALIFICATION\n")
cat("============================================\n\n")

print(
  qualification,
  row.names = FALSE
)

cat("\n============================================\n")
cat("SPECIAL CELL TYPES IN CONTROL\n")
cat("============================================\n\n")

print(
  summary_table[
    summary_table$cell_type %in% special5,
    c(
      "cell_type",
      "stratum",
      "n",
      "mass_median",
      "pct_ge100",
      "pct_ge1000"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

# ==============================================================
# Global decision
# ==============================================================

core12_all_robust <- all(
  qualification$control_classification ==
    "COMMON_ROBUST"
)

core12_all_usable <- all(
  qualification$control_classification %in%
    c(
      "COMMON_ROBUST",
      "COMMON_SUPPORTED"
    )
)

cat("\n============================================\n")
cat("GLOBAL DECISION\n")
cat("============================================\n\n")

cat(
  "CORE12 all COMMON_ROBUST in Control:",
  core12_all_robust,
  "\n"
)

cat(
  "CORE12 all at least COMMON_SUPPORTED in Control:",
  core12_all_usable,
  "\n"
)

if (core12_all_robust) {

  cat(
    "\nFINAL STATUS: PASS_CORE12_SIX_STRATA_ROBUST\n"
  )

} else if (core12_all_usable) {

  cat(
    "\nFINAL STATUS: PASS_CORE12_SIX_STRATA_SUPPORTED\n"
  )

} else {

  cat(
    "\nFINAL STATUS: REVIEW_CORE12_CONTROL_COVERAGE\n"
  )
}

cat(
  "\nOutput:",
  OUT,
  "\n"
)

