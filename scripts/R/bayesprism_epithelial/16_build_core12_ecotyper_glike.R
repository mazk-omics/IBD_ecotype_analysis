#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/home/mazekai/IBD_EcoTyper"

S1_DIR <- file.path(
  ROOT,
  "03_reference/combined_reference/output/14_bayesprism_broad_run",
  "celltype_expression"
)

S2_DIR <- file.path(
  ROOT,
  "03_reference/combined_reference/output/20_bayesprism_epithelial_full_run",
  "celltype_expression"
)

META_FILE <- file.path(
  ROOT,
  "03_reference/combined_reference/output/21_bayesprism_epithelial_qc",
  "GSE193677_metadata_used_for_stage2_qc.tsv"
)

TECH_FILE <- file.path(
  ROOT,
  "03_reference/combined_reference/output/26_panibd_panel_classification",
  "technical_exclusion_samples.tsv"
)

OUT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/27_ecotyper_core12_glike"
)

MATRIX_DIR <- file.path(
  OUT,
  "expression_matrices"
)

dir.create(
  MATRIX_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==============================================================
# Locked Common Discovery Panel
# ==============================================================

stage1_types <- c(
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

stage2_types <- c(
  "Goblet",
  "Enteroendocrine"
)

core12 <- c(
  stage1_types,
  stage2_types
)

source_stage <- c(
  setNames(
    rep("Stage1", length(stage1_types)),
    stage1_types
  ),
  setNames(
    rep("Stage2", length(stage2_types)),
    stage2_types
  )
)

get_z_path <- function(ct) {

  base <- if (
    source_stage[[ct]] == "Stage1"
  ) {
    S1_DIR
  } else {
    S2_DIR
  }

  file.path(
    base,
    paste0(
      "bayesprism_expression_",
      ct,
      ".rds"
    )
  )
}

# ==============================================================
# Unified discovery cohort:
# Control + CD + UC
# minus four locked technical exclusions
# ==============================================================

meta <- read.delim(
  META_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tech <- read.delim(
  TECH_FILE,
  stringsAsFactors = FALSE
)

stopifnot(
  all(
    c(
      "sample_id",
      "disease",
      "region_broad",
      "inflammation"
    ) %in% colnames(meta)
  ),
  "sample_id" %in% colnames(tech)
)

keep_meta <- meta[
  meta$disease %in% c(
    "CD",
    "UC",
    "Control"
  ) &
    !(meta$sample_id %in% tech$sample_id),
  ,
  drop = FALSE
]

stopifnot(
  nrow(keep_meta) == 2486,
  length(unique(keep_meta$sample_id)) == 2486
)

sample_ids <- keep_meta$sample_id

cat("============================================\n")
cat("CORE12 G-LIKE CONSTRUCTION\n")
cat("============================================\n\n")

cat("Samples:", length(sample_ids), "\n")
cat("CD:", sum(keep_meta$disease == "CD"), "\n")
cat("UC:", sum(keep_meta$disease == "UC"), "\n")
cat("Control:", sum(keep_meta$disease == "Control"), "\n")
cat(
  "Colon_rectum:",
  sum(keep_meta$region_broad == "Colon_rectum"),
  "\n"
)
cat(
  "Ileum:",
  sum(keep_meta$region_broad == "Ileum"),
  "\n\n"
)

# ==============================================================
# Save immutable manifests
# ==============================================================

write.table(
  keep_meta,
  file.path(
    OUT,
    "discovery_metadata_2486.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  sample_ids,
  file.path(
    OUT,
    "discovery_sample_ids.txt"
  )
)

writeLines(
  core12,
  file.path(
    OUT,
    "core12_cell_types.txt"
  )
)

# ==============================================================
# Construct G-like matrices
#
# For each cell type:
#
# G_like[n,g] =
#     Z[n,g] / sum_g Z[n,g] * 1e6
#
# Then transpose:
#     genes x samples
#
# This is a CPM-like normalized BayesPrism posterior profile.
# It is NOT TPM and is NOT claimed to be numerically identical
# to CIBERSORTx HiRes.
# ==============================================================

qc_rows <- list()
canonical_genes <- NULL

for (i in seq_along(core12)) {

  ct <- core12[i]

  cat(
    "[",
    i,
    "/",
    length(core12),
    "] ",
    ct,
    "\n",
    sep = ""
  )

  path <- get_z_path(ct)

  stopifnot(
    file.exists(path)
  )

  z <- readRDS(path)

  stopifnot(
    is.matrix(z) ||
      inherits(z, "Matrix")
  )

  # Z contains all original 2490 samples.
  # sample_ids contains the 2486 discovery samples after removing
  # the four locked technical exclusions.

  stopifnot(
    all(sample_ids %in% rownames(z))
  )

  extra_samples <- setdiff(
    rownames(z),
    sample_ids
  )

  stopifnot(
    setequal(
      extra_samples,
      tech$sample_id
    )
  )

  stopifnot(
    length(extra_samples) == 4
  )

  z <- z[
    sample_ids,
    ,
    drop = FALSE
  ]

  stopifnot(
    identical(
      rownames(z),
      sample_ids
    )
  )

  if (is.null(canonical_genes)) {

    canonical_genes <- colnames(z)

  } else {

    stopifnot(
      identical(
        colnames(z),
        canonical_genes
      )
    )
  }

  stopifnot(
    all(is.finite(z)),
    min(z) >= -1e-12
  )

  raw_mass <- rowSums(z)

  if (any(raw_mass <= 0)) {
    stop(
      ct,
      ": zero/non-positive posterior mass remains."
    )
  }

  # Row-wise normalization.
  g <- sweep(
    z,
    MARGIN = 1,
    STATS = raw_mass,
    FUN = "/"
  )

  g <- g * 1e6

  norm_mass <- rowSums(g)

  norm_error <- max(
    abs(
      norm_mass - 1e6
    )
  )

  stopifnot(
    all(is.finite(g)),
    min(g) >= -1e-8
  )

  detected_genes <- rowSums(
    z > 0
  )

  # EcoTyper expects genes x samples.
  g_t <- t(g)

  stopifnot(
    identical(
      rownames(g_t),
      canonical_genes
    ),
    identical(
      colnames(g_t),
      sample_ids
    )
  )

  # ============================================================
  # Write EcoTyper-compatible matrix:
  #
  # Gene   sample1 sample2 ...
  # gene1  ...
  # gene2  ...
  # ============================================================

  dt <- as.data.table(
    g_t,
    keep.rownames = "Gene"
  )

  outfile <- file.path(
    MATRIX_DIR,
    paste0(
      ct,
      ".txt"
    )
  )

  fwrite(
    dt,
    outfile,
    sep = "\t",
    quote = FALSE
  )

  qc_rows[[i]] <- data.frame(
    cell_type = ct,
    source_stage = source_stage[[ct]],

    n_samples = nrow(z),
    n_genes = ncol(z),

    raw_mass_min = min(raw_mass),
    raw_mass_q01 = unname(
      quantile(raw_mass, 0.01)
    ),
    raw_mass_q05 = unname(
      quantile(raw_mass, 0.05)
    ),
    raw_mass_median = median(raw_mass),

    pct_mass_lt100 =
      mean(raw_mass < 100) * 100,

    pct_mass_lt1000 =
      mean(raw_mass < 1000) * 100,

    normalized_mass_min =
      min(norm_mass),

    normalized_mass_max =
      max(norm_mass),

    normalized_mass_max_abs_error =
      norm_error,

    min_detected_genes =
      min(detected_genes),

    q05_detected_genes =
      unname(
        quantile(
          detected_genes,
          0.05
        )
      ),

    median_detected_genes =
      median(detected_genes),

    max_detected_genes =
      max(detected_genes),

    matrix_file =
      outfile,

    stringsAsFactors = FALSE
  )

  cat(
    "  raw mass median:",
    signif(
      median(raw_mass),
      6
    ),
    "\n"
  )

  cat(
    "  detected genes median:",
    median(detected_genes),
    "\n"
  )

  cat(
    "  normalization max error:",
    format(
      norm_error,
      scientific = TRUE
    ),
    "\n"
  )

  cat(
    "  written:",
    outfile,
    "\n\n"
  )

  rm(
    z,
    g,
    g_t,
    dt,
    raw_mass,
    norm_mass,
    detected_genes
  )

  gc()
}

# ==============================================================
# Final QC
# ==============================================================

qc <- do.call(
  rbind,
  qc_rows
)

write.table(
  qc,
  file.path(
    OUT,
    "core12_Glike_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  canonical_genes,
  file.path(
    OUT,
    "gene_universe.txt"
  )
)

expected_files <- file.path(
  MATRIX_DIR,
  paste0(
    core12,
    ".txt"
  )
)

files_ok <- all(
  file.exists(
    expected_files
  )
)

dimensions_ok <- (
  all(qc$n_samples == 2486) &&
    all(qc$n_genes == 13668)
)

normalization_ok <- all(
  qc$normalized_mass_max_abs_error <
    1e-5
)

cat("\n============================================\n")
cat("CORE12 G-LIKE QC\n")
cat("============================================\n\n")

print(
  qc[
    ,
    c(
      "cell_type",
      "source_stage",
      "n_samples",
      "n_genes",
      "raw_mass_median",
      "pct_mass_lt100",
      "pct_mass_lt1000",
      "median_detected_genes",
      "normalized_mass_max_abs_error"
    )
  ],
  row.names = FALSE
)

cat("\nFiles present:", files_ok, "\n")
cat("Dimensions correct:", dimensions_ok, "\n")
cat("Normalization correct:", normalization_ok, "\n")

if (
  files_ok &&
    dimensions_ok &&
    normalization_ok
) {

  cat(
    "\nFINAL STATUS: PASS_READY_FOR_ECOTYPER_PRESORTED\n"
  )

} else {

  cat(
    "\nFINAL STATUS: REVIEW\n"
  )

  quit(
    save = "no",
    status = 1
  )
}

cat(
  "\nOutput:",
  OUT,
  "\n"
)

