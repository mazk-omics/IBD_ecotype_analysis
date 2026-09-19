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
  "03_reference/combined_reference/output/26_panibd_panel_classification"
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

epi4 <- c(
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine"
)

# ==============================================================
# 1. Remove the four common Stage2 zero-mass samples globally
# ==============================================================

zero_epi <- apply(
  as.matrix(dat[, epi4, drop = FALSE]),
  1,
  function(x) all(x <= 0)
)

technical_exclusions <- dat$sample_id[
  zero_epi
]

cat("Technical epithelial-zero exclusions:", sum(zero_epi), "\n")

if (length(technical_exclusions) > 0) {
  cat(
    paste(
      technical_exclusions,
      collapse = "\n"
    ),
    "\n"
  )
}

dat_use <- dat[
  !zero_epi,
  ,
  drop = FALSE
]

stopifnot(
  nrow(dat_use) == 2486
)

# ==============================================================
# 2. Main IBD strata
# ==============================================================

ibd <- dat_use[
  dat_use$disease %in% c("CD", "UC") &
    dat_use$region_broad %in% c("Colon_rectum", "Ileum"),
  ,
  drop = FALSE
]

ibd$stratum <- paste(
  ibd$disease,
  ibd$region_broad,
  sep = "_"
)

expected_strata <- c(
  "CD_Colon_rectum",
  "CD_Ileum",
  "UC_Colon_rectum",
  "UC_Ileum"
)

stopifnot(
  all(expected_strata %in% unique(ibd$stratum))
)

cat("\nIBD stratum sizes:\n")
print(
  table(ibd$stratum)
)

# ==============================================================
# 3. Cell type × disease × region summary
# ==============================================================

rows <- list()
counter <- 1L

for (ct in cell_types) {

  for (st in expected_strata) {

    x <- ibd[
      ibd$stratum == st,
      ct,
      drop = TRUE
    ]

    rows[[counter]] <- data.frame(
      cell_type = ct,
      stratum = st,
      n = length(x),

      mass_min = min(x),
      mass_q05 = unname(quantile(x, 0.05)),
      mass_q25 = unname(quantile(x, 0.25)),
      mass_median = median(x),
      mass_q75 = unname(quantile(x, 0.75)),
      mass_q95 = unname(quantile(x, 0.95)),

      pct_ge10 = mean(x >= 10) * 100,
      pct_ge100 = mean(x >= 100) * 100,
      pct_ge1000 = mean(x >= 1000) * 100,

      stringsAsFactors = FALSE
    )

    counter <- counter + 1L
  }
}

stratum_summary <- do.call(
  rbind,
  rows
)

write.table(
  stratum_summary,
  file.path(
    OUT,
    "celltype_by_IBD_region_stratum.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================
# 4. Region-level summaries across CD + UC
# ==============================================================

region_rows <- list()
counter <- 1L

for (ct in cell_types) {

  for (region in c("Colon_rectum", "Ileum")) {

    x <- ibd[
      ibd$region_broad == region,
      ct,
      drop = TRUE
    ]

    region_rows[[counter]] <- data.frame(
      cell_type = ct,
      region = region,
      n = length(x),
      median_mass = median(x),
      pct_ge100 = mean(x >= 100) * 100,
      pct_ge1000 = mean(x >= 1000) * 100,
      stringsAsFactors = FALSE
    )

    counter <- counter + 1L
  }
}

region_summary <- do.call(
  rbind,
  region_rows
)

# ==============================================================
# 5. Transparent technical classification
# ==============================================================

class_rows <- list()

for (i in seq_along(cell_types)) {

  ct <- cell_types[i]

  ss <- stratum_summary[
    stratum_summary$cell_type == ct,
    ,
    drop = FALSE
  ]

  ss <- ss[
    match(
      expected_strata,
      ss$stratum
    ),
    ,
    drop = FALSE
  ]

  rs <- region_summary[
    region_summary$cell_type == ct,
    ,
    drop = FALSE
  ]

  colon_med <- rs$median_mass[
    rs$region == "Colon_rectum"
  ]

  ileum_med <- rs$median_mass[
    rs$region == "Ileum"
  ]

  min_med <- min(
    ss$mass_median
  )

  min_pct100 <- min(
    ss$pct_ge100
  )

  min_pct1000 <- min(
    ss$pct_ge1000
  )

  n_strata_med_lt100 <- sum(
    ss$mass_median < 100
  )

  # Strong directional regional restriction.
  region_ratio <- max(
    colon_med,
    ileum_med
  ) /
    max(
      min(
        colon_med,
        ileum_med
      ),
      1e-12
    )

  one_region_high_other_low <- (
    max(colon_med, ileum_med) >= 1000 &&
      min(colon_med, ileum_med) < 100 &&
      region_ratio >= 50
  )

  classification <- "REVIEW"

  if (n_strata_med_lt100 >= 3) {

    classification <- "LOW_INFORMATION"

  } else if (one_region_high_other_low) {

    classification <- "REGION_RESTRICTED"

  } else if (
    min_med >= 1000 &&
      min_pct100 >= 90 &&
      min_pct1000 >= 80
  ) {

    classification <- "COMMON_ROBUST"

  } else if (
    min_med >= 100 &&
      min_pct100 >= 80
  ) {

    classification <- "COMMON_SUPPORTED"
  }

  class_rows[[i]] <- data.frame(
    cell_type = ct,
    classification = classification,

    min_stratum_median_mass = min_med,
    min_stratum_pct_ge100 = min_pct100,
    min_stratum_pct_ge1000 = min_pct1000,

    n_of_4_strata_median_lt100 =
      n_strata_med_lt100,

    colon_median_mass = colon_med,
    ileum_median_mass = ileum_med,
    colon_to_ileum_or_inverse_ratio =
      region_ratio,

    stringsAsFactors = FALSE
  )
}

classification <- do.call(
  rbind,
  class_rows
)

classification <- classification[
  order(
    factor(
      classification$classification,
      levels = c(
        "COMMON_ROBUST",
        "COMMON_SUPPORTED",
        "REGION_RESTRICTED",
        "LOW_INFORMATION",
        "REVIEW"
      )
    ),
    classification$cell_type
  ),
  ,
  drop = FALSE
]

write.table(
  classification,
  file.path(
    OUT,
    "panibd_ecotyper_panel_classification.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    sample_id = technical_exclusions
  ),
  file.path(
    OUT,
    "technical_exclusion_samples.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================
# 6. Print concise evidence table
# ==============================================================

cat("\n============================================\n")
cat("PAN-IBD ECOTYPER PANEL CLASSIFICATION\n")
cat("============================================\n\n")

print(
  classification,
  row.names = FALSE
)

cat("\n============================================\n")
cat("FOUR IBD STRATA — MEDIAN MASS / % >= 100\n")
cat("============================================\n\n")

for (ct in cell_types) {

  x <- stratum_summary[
    stratum_summary$cell_type == ct,
    c(
      "stratum",
      "mass_median",
      "pct_ge100",
      "pct_ge1000"
    ),
    drop = FALSE
  ]

  cat("\n--- ", ct, " ---\n", sep = "")
  print(
    x,
    row.names = FALSE
  )
}

cat("\nFINAL STATUS: CLASSIFICATION_COMPLETE\n")
cat("Output:", OUT, "\n")

