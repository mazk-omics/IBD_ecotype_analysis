#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# ==============================================================================
# Stage 2 epithelial BayesPrism
# QC-03: Tuft discriminative-gene support audit
#
# Goal:
# Determine whether the large Tuft posterior allocation is supported by
# genes that actually distinguish Tuft from the other 6 epithelial references.
#
# NO model fitting.
# NO modification of existing results.
# ==============================================================================


# ==============================================================================
# 1. Paths
# ==============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

RUN_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "20_bayesprism_epithelial_full_run"
)

QC_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "21_bayesprism_epithelial_qc"
)

RESULT_RDS <- file.path(
  RUN_DIR,
  "bayesprism_epithelial7_result.rds"
)

Z_DIR <- file.path(
  RUN_DIR,
  "celltype_expression"
)

dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


OUT_GENE <- file.path(
  QC_DIR,
  "stage2_tuft_discriminative_support_by_gene.tsv"
)

OUT_CLASS <- file.path(
  QC_DIR,
  "stage2_tuft_discriminative_support_by_class.tsv"
)

OUT_TOP <- file.path(
  QC_DIR,
  "stage2_tuft_strong_specific_genes.tsv"
)

OUT_CANONICAL <- file.path(
  QC_DIR,
  "stage2_tuft_canonical_marker_posterior_support.tsv"
)

OUT_LOG <- file.path(
  QC_DIR,
  "stage2_tuft_discriminative_support_audit.log"
)

OUT_SESSION <- file.path(
  QC_DIR,
  "tuft_discriminative_support_sessionInfo.txt"
)


# ==============================================================================
# 2. Logging
# ==============================================================================

log_con <- file(
  OUT_LOG,
  open = "wt"
)

sink(
  log_con,
  type = "output",
  split = TRUE
)

sink(
  log_con,
  type = "message"
)

on.exit({

  while (sink.number(type = "message") > 0) {
    sink(type = "message")
  }

  while (sink.number(type = "output") > 0) {
    sink(type = "output")
  }

  close(log_con)

}, add = TRUE)


write_tsv <- function(x, path) {

  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


assert_true <- function(x, msg) {

  if (!isTRUE(x)) {
    stop(msg, call. = FALSE)
  }
}


cat("============================================================\n")
cat("Stage 2 Tuft discriminative-support audit\n")
cat("============================================================\n\n")

cat(
  "START:",
  format(Sys.time()),
  "\n\n"
)


# ==============================================================================
# 3. Load BayesPrism result
# ==============================================================================

suppressPackageStartupMessages(
  library(BayesPrism)
)

assert_true(
  file.exists(RESULT_RDS),
  "Stage 2 BayesPrism result missing."
)

bp <- readRDS(
  RESULT_RDS
)

assert_true(
  is(bp, "BayesPrism"),
  "Result object is not BayesPrism."
)


# IMPORTANT:
# posterior.initial.cellType@Z was generated using the ORIGINAL reference.
# Therefore reference specificity for this audit must use phi0.

phi <- bp@prism@phi_cellType@phi


assert_true(
  "Tuft" %in% rownames(phi),
  "Tuft absent from original reference."
)


components <- rownames(phi)
genes <- colnames(phi)

other_components <- setdiff(
  components,
  "Tuft"
)


cat(
  "[PASS] Original reference:",
  nrow(phi),
  "components x",
  ncol(phi),
  "genes\n\n"
)


# ==============================================================================
# 4. Reference-based specificity
# ==============================================================================

phi_tuft <- phi[
  "Tuft",
  ,
  drop = TRUE
]

phi_other <- phi[
  other_components,
  ,
  drop = FALSE
]


max_other <- apply(
  phi_other,
  2,
  max
)


second_component <- other_components[
  max.col(
    t(phi_other),
    ties.method = "first"
  )
]


reference_sum <- colSums(phi)

EPS <- 1e-12


tuft_reference_share <- phi_tuft /
  reference_sum


tuft_vs_max_other <- (
  phi_tuft + EPS
) / (
  max_other + EPS
)


all_top_component <- components[
  max.col(
    t(phi),
    ties.method = "first"
  )
]


# ==============================================================================
# 5. Define data-driven specificity classes
# ==============================================================================

specificity_class <- rep(
  "non_Tuft",
  length(genes)
)


# Tuft is the highest-reference component but not strongly dominant
specificity_class[
  all_top_component == "Tuft"
] <- "Tuft_top"


# At least 2-fold over every other epithelial reference
specificity_class[
  all_top_component == "Tuft" &
    tuft_vs_max_other >= 2
] <- "Tuft_enriched_2x"


# At least 5-fold over every other epithelial reference
specificity_class[
  all_top_component == "Tuft" &
    tuft_vs_max_other >= 5
] <- "Tuft_specific_5x"


# Very strong specificity
specificity_class[
  all_top_component == "Tuft" &
    tuft_vs_max_other >= 10
] <- "Tuft_specific_10x"


# Nearly unique to Tuft in normalized reference
specificity_class[
  all_top_component == "Tuft" &
    tuft_reference_share >= 0.80
] <- "Tuft_reference_share_ge80"
  

# ==============================================================================
# 6. Aggregate posterior Z gene masses
# ==============================================================================

cat("Aggregating posterior Z across 7 components...\n\n")


total_Z_gene <- rep(
  0,
  length(genes)
)

names(total_Z_gene) <- genes


component_Z_gene <- matrix(
  0,
  nrow = length(components),
  ncol = length(genes),
  dimnames = list(
    components,
    genes
  )
)


for (component in components) {

  z_file <- file.path(
    Z_DIR,
    paste0(
      "bayesprism_expression_",
      component,
      ".rds"
    )
  )

  assert_true(
    file.exists(z_file),
    paste(
      "Missing Z matrix:",
      component
    )
  )


  cat(
    "Loading:",
    component,
    "\n"
  )


  z <- readRDS(
    z_file
  )


  assert_true(
    identical(
      colnames(z),
      genes
    ),
    paste(
      component,
      "gene order mismatch."
    )
  )


  gene_mass <- colSums(
    z
  )


  component_Z_gene[
    component,
  ] <- gene_mass


  total_Z_gene <- total_Z_gene +
    gene_mass


  rm(
    z,
    gene_mass
  )

  gc(verbose = FALSE)
}


cat(
  "\n[PASS] Posterior Z aggregation complete\n\n"
)


# ==============================================================================
# 7. Posterior Tuft allocation per gene
# ==============================================================================

tuft_Z_gene <- component_Z_gene[
  "Tuft",
  ,
  drop = TRUE
]


posterior_tuft_share <- rep(
  NA_real_,
  length(genes)
)


positive_total <- total_Z_gene > 0


posterior_tuft_share[
  positive_total
] <- tuft_Z_gene[
  positive_total
] / total_Z_gene[
  positive_total
]


# ==============================================================================
# 8. Gene-level audit table
# ==============================================================================

gene_table <- data.frame(

  gene = genes,

  reference_Tuft =
    phi_tuft,

  reference_max_other =
    max_other,

  max_other_component =
    second_component,

  reference_top_component =
    all_top_component,

  tuft_reference_share =
    tuft_reference_share,

  tuft_vs_max_other_fold =
    tuft_vs_max_other,

  specificity_class =
    specificity_class,

  posterior_Tuft_mass =
    tuft_Z_gene,

  posterior_all_epithelial_mass =
    total_Z_gene,

  posterior_Tuft_share =
    posterior_tuft_share,

  stringsAsFactors = FALSE
)


# Ratio of observed posterior allocation to reference share.
# Values >>1 mean the gene is allocated to Tuft more strongly than
# suggested by the original normalized reference proportions.

gene_table$posterior_vs_reference_share <- (
  gene_table$posterior_Tuft_share + EPS
) / (
  gene_table$tuft_reference_share + EPS
)


gene_table <- gene_table[
  order(
    -gene_table$posterior_Tuft_mass
  ),
  ,
  drop = FALSE
]


write_tsv(
  gene_table,
  OUT_GENE
)


# ==============================================================================
# 9. Class-level support
# ==============================================================================

classes_to_report <- c(
  "Tuft_reference_share_ge80",
  "Tuft_specific_10x",
  "Tuft_specific_5x",
  "Tuft_enriched_2x",
  "Tuft_top",
  "non_Tuft"
)


class_rows <- lapply(
  classes_to_report,
  function(cls) {

    idx <- specificity_class == cls

    if (!any(idx)) {

      return(
        data.frame(
          class = cls,
          n_genes = 0,
          reference_Tuft_mass_fraction = NA_real_,
          Tuft_Z_mass_fraction = NA_real_,
          median_reference_share = NA_real_,
          median_posterior_Tuft_share = NA_real_,
          weighted_posterior_Tuft_share = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }


    reference_fraction <- sum(
      phi_tuft[idx]
    ) / sum(
      phi_tuft
    )


    tuft_Z_fraction <- sum(
      tuft_Z_gene[idx]
    ) / sum(
      tuft_Z_gene
    )


    weighted_posterior_share <- sum(
      tuft_Z_gene[idx]
    ) / sum(
      total_Z_gene[idx]
    )


    data.frame(

      class =
        cls,

      n_genes =
        sum(idx),

      reference_Tuft_mass_fraction =
        reference_fraction,

      Tuft_Z_mass_fraction =
        tuft_Z_fraction,

      median_reference_share =
        median(
          tuft_reference_share[idx],
          na.rm = TRUE
        ),

      median_posterior_Tuft_share =
        median(
          posterior_tuft_share[idx],
          na.rm = TRUE
        ),

      weighted_posterior_Tuft_share =
        weighted_posterior_share,

      stringsAsFactors = FALSE
    )
  }
)


class_table <- do.call(
  rbind,
  class_rows
)


write_tsv(
  class_table,
  OUT_CLASS
)


# ==============================================================================
# 10. Strong Tuft-specific genes
# ==============================================================================

strong_idx <- (
  all_top_component == "Tuft" &
    tuft_vs_max_other >= 5
)


strong_table <- gene_table[
  gene_table$gene %in% genes[
    strong_idx
  ],
  ,
  drop = FALSE
]


strong_table <- strong_table[
  order(
    -strong_table$posterior_Tuft_mass
  ),
  ,
  drop = FALSE
]


write_tsv(
  strong_table,
  OUT_TOP
)


# ==============================================================================
# 11. Canonical marker audit
# ==============================================================================

canonical <- c(
  "POU2F3",
  "TRPM5",
  "PLCB2",
  "GNAT3",
  "GFI1B",
  "AVIL",
  "IL25",
  "PTGS1",
  "DCLK1",
  "SH2D6"
)


canonical_table <- gene_table[
  match(
    canonical,
    gene_table$gene
  ),
  ,
  drop = FALSE
]


canonical_table$canonical_marker <- canonical


canonical_table <- canonical_table[
  ,
  c(
    "canonical_marker",
    "gene",
    "reference_Tuft",
    "reference_max_other",
    "tuft_reference_share",
    "tuft_vs_max_other_fold",
    "specificity_class",
    "posterior_Tuft_mass",
    "posterior_all_epithelial_mass",
    "posterior_Tuft_share",
    "posterior_vs_reference_share"
  )
]


write_tsv(
  canonical_table,
  OUT_CANONICAL
)


# ==============================================================================
# 12. Additional biologically useful support sets
# ==============================================================================

idx_top <- all_top_component == "Tuft"

idx_2x <- idx_top &
  tuft_vs_max_other >= 2

idx_5x <- idx_top &
  tuft_vs_max_other >= 5

idx_10x <- idx_top &
  tuft_vs_max_other >= 10

idx_share80 <- idx_top &
  tuft_reference_share >= 0.80


support_summary <- function(idx) {

  if (!any(idx)) {

    return(
      c(
        n = 0,
        ref_mass_fraction = NA,
        Z_mass_fraction = NA,
        epithelial_reads_allocated_to_Tuft = NA
      )
    )
  }


  c(

    n =
      sum(idx),

    ref_mass_fraction =
      sum(phi_tuft[idx]) /
        sum(phi_tuft),

    Z_mass_fraction =
      sum(tuft_Z_gene[idx]) /
        sum(tuft_Z_gene),

    epithelial_reads_allocated_to_Tuft =
      sum(tuft_Z_gene[idx]) /
        sum(total_Z_gene[idx])
  )
}


support_sets <- rbind(

  Tuft_top =
    support_summary(
      idx_top
    ),

  Tuft_2x =
    support_summary(
      idx_2x
    ),

  Tuft_5x =
    support_summary(
      idx_5x
    ),

  Tuft_10x =
    support_summary(
      idx_10x
    ),

  Tuft_reference_share_ge80 =
    support_summary(
      idx_share80
    )
)


# ==============================================================================
# 13. Console report
# ==============================================================================

cat("============================================================\n")
cat("TUFT DISCRIMINATIVE SUPPORT SETS\n")
cat("============================================================\n\n")

print(
  support_sets
)


cat("\n")
cat("============================================================\n")
cat("SPECIFICITY CLASS SUMMARY\n")
cat("============================================================\n\n")

print(
  class_table,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("TOP 30 STRONG TUFT-SPECIFIC GENES\n")
cat("Reference criterion: Tuft >= 5x every other epithelial type\n")
cat("============================================================\n\n")

print(
  head(
    strong_table[
      ,
      c(
        "gene",
        "tuft_reference_share",
        "tuft_vs_max_other_fold",
        "posterior_Tuft_mass",
        "posterior_all_epithelial_mass",
        "posterior_Tuft_share"
      )
    ],
    30
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("CANONICAL TUFT MARKER POSTERIOR SUPPORT\n")
cat("============================================================\n\n")

print(
  canonical_table,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("KEY INTERPRETATION\n")
cat("============================================================\n\n")

cat(
  "For truly Tuft-discriminative genes, posterior_Tuft_share should be\n",
  "substantially enriched toward Tuft.\n\n",
  sep = ""
)

cat(
  "If canonical / >=5x Tuft-specific genes are overwhelmingly assigned\n",
  "to Tuft, the reference identity itself is valid.\n\n",
  sep = ""
)

cat(
  "If these genes support Tuft strongly but account for only a small\n",
  "fraction of total Tuft Z mass, the high global Tuft theta is likely\n",
  "driven by broad transcriptomic similarity / identifiability rather\n",
  "than lack of canonical Tuft markers.\n",
  sep = ""
)


capture.output(
  sessionInfo(),
  file = OUT_SESSION
)


cat("\nOutput directory:\n")
cat(QC_DIR, "\n")

cat(
  "\nEND:",
  format(Sys.time()),
  "\n"
)