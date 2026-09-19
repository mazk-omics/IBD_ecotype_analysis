#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

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

OUT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/23_ecotyper_interface_audit"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

stage1_types <- c(
  "B",
  "CD4_T",
  "CD8_T",
  "DC",
  "Endothelial",
  "Fibroblast",
  "Glial",
  "ILC",
  "Mast",
  "Monocyte_Macrophage",
  "NK",
  "Pericyte",
  "Plasma"
)

stage2_types <- c(
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine"
)

all_types <- c(
  stage1_types,
  stage2_types
)

source_stage <- c(
  rep("Stage1", length(stage1_types)),
  rep("Stage2", length(stage2_types))
)

names(source_stage) <- all_types


get_path <- function(cell_type) {

  base <- if (
    source_stage[[cell_type]] == "Stage1"
  ) {
    S1_DIR
  } else {
    S2_DIR
  }

  file.path(
    base,
    paste0(
      "bayesprism_expression_",
      cell_type,
      ".rds"
    )
  )
}


cat("============================================\n")
cat("EcoTyper Z-interface audit\n")
cat("============================================\n\n")

# ------------------------------------------------------------------
# Establish canonical sample/gene universe from first matrix
# ------------------------------------------------------------------

first_path <- get_path(all_types[1])

stopifnot(
  file.exists(first_path)
)

first <- readRDS(first_path)

stopifnot(
  is.matrix(first) ||
    inherits(first, "Matrix")
)

canonical_samples <- rownames(first)
canonical_genes <- colnames(first)

cat(
  "Canonical dimensions:",
  nrow(first),
  "samples x",
  ncol(first),
  "genes\n"
)

cat(
  "Canonical samples:",
  length(canonical_samples),
  "\n"
)

cat(
  "Canonical genes:",
  length(canonical_genes),
  "\n\n"
)

rm(first)
gc()

summary_rows <- list()
mass_table <- data.frame(
  sample_id = canonical_samples,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# Sequential audit
# ------------------------------------------------------------------

for (i in seq_along(all_types)) {

  ct <- all_types[i]
  path <- get_path(ct)

  cat(
    "[",
    i,
    "/",
    length(all_types),
    "] ",
    ct,
    "\n",
    sep = ""
  )

  if (!file.exists(path)) {
    stop(
      "Missing Z file: ",
      path
    )
  }

  z <- readRDS(path)

  if (!(
    is.matrix(z) ||
      inherits(z, "Matrix")
  )) {
    stop(
      ct,
      ": object is not matrix-like."
    )
  }

  sample_set_match <- setequal(
    rownames(z),
    canonical_samples
  )

  gene_set_match <- setequal(
    colnames(z),
    canonical_genes
  )

  sample_order_match <- identical(
    rownames(z),
    canonical_samples
  )

  gene_order_match <- identical(
    colnames(z),
    canonical_genes
  )

  finite_ok <- all(
    is.finite(z)
  )

  min_value <- min(z)

  nonnegative_ok <- (
    min_value >= -1e-12
  )

  rs <- rowSums(z)

  zero_mass <- sum(
    rs <= 0
  )

  mass_lt_1 <- sum(
    rs < 1
  )

  mass_lt_10 <- sum(
    rs < 10
  )

  mass_lt_100 <- sum(
    rs < 100
  )

  mass_lt_1000 <- sum(
    rs < 1000
  )

  q <- quantile(
    rs,
    probs = c(
      0,
      0.01,
      0.05,
      0.25,
      0.50,
      0.75,
      0.95,
      0.99,
      1
    ),
    names = FALSE
  )

  summary_rows[[i]] <- data.frame(
    cell_type = ct,
    source_stage = source_stage[[ct]],
    n_samples = nrow(z),
    n_genes = ncol(z),
    sample_set_match = sample_set_match,
    sample_order_match = sample_order_match,
    gene_set_match = gene_set_match,
    gene_order_match = gene_order_match,
    finite_ok = finite_ok,
    nonnegative_ok = nonnegative_ok,
    min_value = min_value,
    zero_mass_samples = zero_mass,
    mass_lt_1 = mass_lt_1,
    mass_lt_10 = mass_lt_10,
    mass_lt_100 = mass_lt_100,
    mass_lt_1000 = mass_lt_1000,
    mass_min = q[1],
    mass_q01 = q[2],
    mass_q05 = q[3],
    mass_q25 = q[4],
    mass_median = q[5],
    mass_q75 = q[6],
    mass_q95 = q[7],
    mass_q99 = q[8],
    mass_max = q[9],
    stringsAsFactors = FALSE
  )

  # Align only for the diagnostic mass table.
  # No samples are removed here.
  if (sample_set_match) {

    rs <- rs[
      match(
        canonical_samples,
        rownames(z)
      )
    ]

    mass_table[[ct]] <- rs
  }

  cat(
    "  dims:",
    nrow(z), "x", ncol(z),
    "\n"
  )

  cat(
    "  sample set/order:",
    sample_set_match,
    "/",
    sample_order_match,
    "\n"
  )

  cat(
    "  gene set/order:",
    gene_set_match,
    "/",
    gene_order_match,
    "\n"
  )

  cat(
    "  finite/nonnegative:",
    finite_ok,
    "/",
    nonnegative_ok,
    "\n"
  )

  cat(
    "  zero mass:",
    zero_mass,
    "\n"
  )

  cat(
    "  mass median:",
    signif(median(rs), 6),
    "\n\n"
  )

  rm(z, rs)
  gc()
}

audit <- do.call(
  rbind,
  summary_rows
)

# ------------------------------------------------------------------
# Global interface gate
# ------------------------------------------------------------------

schema_pass <- all(
  audit$n_samples == length(canonical_samples)
) &&
  all(
    audit$n_genes == length(canonical_genes)
  ) &&
  all(
    audit$sample_set_match
  ) &&
  all(
    audit$gene_set_match
  ) &&
  all(
    audit$finite_ok
  ) &&
  all(
    audit$nonnegative_ok
  )

normalization_directly_possible <- (
  audit$zero_mass_samples == 0
)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------

write.table(
  audit,
  file.path(
    OUT,
    "Z_interface_audit_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  mass_table,
  file.path(
    OUT,
    "Z_allocated_mass_by_sample.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  canonical_samples,
  file.path(
    OUT,
    "canonical_sample_ids.txt"
  )
)

writeLines(
  canonical_genes,
  file.path(
    OUT,
    "canonical_gene_ids.txt"
  )
)

cat("\n============================================\n")
cat("GLOBAL SUMMARY\n")
cat("============================================\n\n")

print(
  audit[
    ,
    c(
      "cell_type",
      "source_stage",
      "n_samples",
      "n_genes",
      "zero_mass_samples",
      "mass_lt_100",
      "mass_lt_1000",
      "mass_median"
    )
  ],
  row.names = FALSE
)

cat("\nSchema/interface PASS:", schema_pass, "\n")

cat(
  "All 17 types directly normalizable without zero-mass handling:",
  all(normalization_directly_possible),
  "\n"
)

if (!schema_pass) {

  cat("\nFINAL STATUS: FAIL_SCHEMA\n")
  quit(
    save = "no",
    status = 1
  )
}

if (
  any(
    audit$zero_mass_samples > 0
  )
) {

  cat("\nFINAL STATUS: PASS_SCHEMA_ZERO_MASS_REQUIRES_POLICY\n")

} else {

  cat("\nFINAL STATUS: PASS_READY_FOR_NORMALIZATION\n")
}

cat(
  "\nOutput:",
  OUT,
  "\n"
)

