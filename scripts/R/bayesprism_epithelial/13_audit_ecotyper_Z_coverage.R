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

OUT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/24_ecotyper_Z_coverage_audit"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(MASS_FILE),
  file.exists(META_FILE)
)

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

required_meta <- c(
  "sample_id",
  "region_broad",
  "disease",
  "inflammation"
)

stopifnot(
  all(required_meta %in% colnames(meta))
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

cell_types <- setdiff(
  colnames(mass),
  "sample_id"
)

cat("Samples:", nrow(dat), "\n")
cat("Cell types:", length(cell_types), "\n\n")


# ==============================================================
# 1. Zero-mass epithelial samples
# ==============================================================

epi4 <- c(
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine"
)

zero_mat <- sapply(
  epi4,
  function(ct) dat[[ct]] <= 0
)

common_zero <- apply(
  zero_mat,
  1,
  all
)

any_zero <- apply(
  zero_mat,
  1,
  any
)

cat("============================================\n")
cat("EPITHELIAL ZERO-MASS AUDIT\n")
cat("============================================\n\n")

cat(
  "Samples zero in ALL four epithelial profiles:",
  sum(common_zero),
  "\n"
)

cat(
  "Samples zero in ANY epithelial profile:",
  sum(any_zero),
  "\n\n"
)

if (any(any_zero)) {

  zero_table <- dat[
    any_zero,
    c(
      "sample_id",
      "region_broad",
      "disease",
      "inflammation",
      epi4
    ),
    drop = FALSE
  ]

  print(
    zero_table,
    row.names = FALSE
  )

  write.table(
    zero_table,
    file.path(
      OUT,
      "epithelial_zero_mass_samples.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


# ==============================================================
# 2. Cell-type mass by anatomical region
# ==============================================================

thresholds <- c(
  10,
  100,
  1000
)

rows <- list()
counter <- 1L

for (region in sort(unique(dat$region_broad))) {

  sub <- dat[
    dat$region_broad == region,
    ,
    drop = FALSE
  ]

  for (ct in cell_types) {

    x <- sub[[ct]]

    rows[[counter]] <- data.frame(
      region_broad = region,
      cell_type = ct,
      n = length(x),

      mass_min = min(x),
      mass_q05 = unname(quantile(x, 0.05)),
      mass_q25 = unname(quantile(x, 0.25)),
      mass_median = median(x),
      mass_q75 = unname(quantile(x, 0.75)),
      mass_q95 = unname(quantile(x, 0.95)),

      n_zero = sum(x <= 0),
      pct_lt10 = mean(x < 10) * 100,
      pct_lt100 = mean(x < 100) * 100,
      pct_lt1000 = mean(x < 1000) * 100,

      stringsAsFactors = FALSE
    )

    counter <- counter + 1L
  }
}

region_summary <- do.call(
  rbind,
  rows
)

write.table(
  region_summary,
  file.path(
    OUT,
    "celltype_mass_by_region.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\n============================================\n")
cat("KEY CELL TYPES BY REGION\n")
cat("============================================\n\n")

key_types <- c(
  "DC",
  "NK",
  "Plasma",
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine"
)

print(
  region_summary[
    region_summary$cell_type %in% key_types,
    c(
      "region_broad",
      "cell_type",
      "n",
      "mass_median",
      "pct_lt10",
      "pct_lt100",
      "pct_lt1000"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


# ==============================================================
# 3. Common-cohort retention diagnostics
#
# These thresholds are DIAGNOSTIC ONLY.
# They are not yet exclusion criteria.
# ==============================================================

core_non_epi <- c(
  "B",
  "CD4_T",
  "CD8_T",
  "Endothelial",
  "Fibroblast",
  "Glial",
  "ILC",
  "Mast",
  "Monocyte_Macrophage",
  "Pericyte"
)

core12 <- c(
  core_non_epi,
  "Goblet",
  "Enteroendocrine"
)

colon13 <- c(
  core_non_epi,
  "Colonic_absorptive",
  "Goblet",
  "Enteroendocrine"
)

ileum13 <- c(
  core_non_epi,
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine"
)

all17 <- cell_types


calc_retention <- function(
  data,
  types,
  cohort_name
) {

  output <- list()

  ts <- c(
    0,
    10,
    100,
    1000
  )

  for (i in seq_along(ts)) {

    threshold <- ts[i]

    m <- as.matrix(
      data[, types, drop = FALSE]
    )

    if (threshold == 0) {

      keep <- apply(
        m,
        1,
        function(x) all(x > 0)
      )

    } else {

      keep <- apply(
        m,
        1,
        function(x) all(x >= threshold)
      )
    }

    output[[i]] <- data.frame(
      cohort = cohort_name,
      n_cell_types = length(types),
      threshold = threshold,
      n_total = nrow(data),
      n_retained = sum(keep),
      pct_retained = mean(keep) * 100,
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    output
  )
}


retention <- list()

retention[[1]] <- calc_retention(
  dat,
  all17,
  "ALL_17_ALL_SAMPLES"
)

retention[[2]] <- calc_retention(
  dat,
  core12,
  "CORE_12_ALL_SAMPLES"
)

colon_dat <- dat[
  dat$region_broad == "Colon_rectum",
  ,
  drop = FALSE
]

ileum_dat <- dat[
  dat$region_broad == "Ileum",
  ,
  drop = FALSE
]

retention[[3]] <- calc_retention(
  colon_dat,
  colon13,
  "COLON_13"
)

retention[[4]] <- calc_retention(
  ileum_dat,
  ileum13,
  "ILEUM_13"
)

retention <- do.call(
  rbind,
  retention
)

write.table(
  retention,
  file.path(
    OUT,
    "common_cohort_retention.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\n============================================\n")
cat("COMMON-COHORT RETENTION\n")
cat("Diagnostic thresholds only\n")
cat("============================================\n\n")

print(
  retention,
  row.names = FALSE
)


# ==============================================================
# 4. Rare-risk types by disease / inflammation
# ==============================================================

rare_types <- c(
  "DC",
  "NK",
  "Plasma"
)

risk_rows <- list()
counter <- 1L

for (ct in rare_types) {

  for (disease in sort(unique(dat$disease))) {

    x <- dat[
      dat$disease == disease,
      ct,
      drop = TRUE
    ]

    risk_rows[[counter]] <- data.frame(
      cell_type = ct,
      stratum = "disease",
      level = disease,
      n = length(x),
      median_mass = median(x),
      pct_lt100 = mean(x < 100) * 100,
      pct_lt1000 = mean(x < 1000) * 100,
      stringsAsFactors = FALSE
    )

    counter <- counter + 1L
  }

  for (inflam in sort(unique(dat$inflammation))) {

    x <- dat[
      dat$inflammation == inflam,
      ct,
      drop = TRUE
    ]

    risk_rows[[counter]] <- data.frame(
      cell_type = ct,
      stratum = "inflammation",
      level = inflam,
      n = length(x),
      median_mass = median(x),
      pct_lt100 = mean(x < 100) * 100,
      pct_lt1000 = mean(x < 1000) * 100,
      stringsAsFactors = FALSE
    )

    counter <- counter + 1L
  }
}

rare_summary <- do.call(
  rbind,
  risk_rows
)

write.table(
  rare_summary,
  file.path(
    OUT,
    "rare_celltype_mass_by_phenotype.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================\n")
cat("DC / NK / PLASMA STRATIFIED SUMMARY\n")
cat("============================================\n\n")

print(
  rare_summary,
  row.names = FALSE
)

cat(
  "\nFINAL STATUS: COVERAGE_AUDIT_COMPLETE\n"
)

cat(
  "Output:",
  OUT,
  "\n"
)

