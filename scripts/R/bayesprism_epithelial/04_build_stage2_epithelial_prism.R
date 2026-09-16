#!/usr/bin/env Rscript

# ==============================================================================
# IBD EcoTyper
# Stage 2 epithelial-only BayesPrism
#
# S2-B2b: Construct and audit formal Stage 2 prism object
#
# DESIGN:
#   - 7 epithelial reference components
#   - cell.type = cell.state (1:1)
#   - Stage 1 Epithelial Z is the Stage 2 pseudo-mixture
#   - retain the inherited 13,668-gene Stage 1 gene universe
#   - DISABLE second-stage new.prism outlier filtering
#
# IMPORTANT:
#   This script DOES NOT run run.prism().
# ==============================================================================


# ==============================================================================
# 0. Settings
# ==============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(Matrix)
  library(BayesPrism)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "17_bayesprism_epithelial_input"
)

PREPROCESS_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "18_bayesprism_epithelial_preprocessing"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "19_bayesprism_epithelial_prism"
)

SCRIPT_DIR <- file.path(
  PROJECT_ROOT,
  "scripts",
  "R",
  "bayesprism_epithelial"
)

dir.create(
  SCRIPT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Inputs
# ==============================================================================

CLEAN_REFERENCE_RDS <- file.path(
  PREPROCESS_DIR,
  "bayesprism_epithelial7_reference_cleaned.rds"
)

MIXTURE_RDS <- file.path(
  INPUT_DIR,
  "bayesprism_stage2_epithelial_pseudomixture.rds"
)

OUTLIER_AUDIT_RDS <- file.path(
  PREPROCESS_DIR,
  "stage2_newprism_default_outliers.tsv"
)


# ==============================================================================
# 3. Outputs
# ==============================================================================

OUT_PRISM <- file.path(
  OUTPUT_DIR,
  "bayesprism_epithelial7_prism_object.rds"
)

OUT_AUDIT <- file.path(
  OUTPUT_DIR,
  "stage2_prism_construction_audit.tsv"
)

OUT_GENE_UNIVERSE <- file.path(
  OUTPUT_DIR,
  "stage2_prism_gene_universe.txt"
)

OUT_MAP <- file.path(
  OUTPUT_DIR,
  "stage2_prism_type_state_map.tsv"
)

OUT_MANIFEST <- file.path(
  OUTPUT_DIR,
  "stage2_prism_manifest.tsv"
)

OUT_LOG <- file.path(
  OUTPUT_DIR,
  "stage2_prism_construction.log"
)

OUT_SESSIONINFO <- file.path(
  OUTPUT_DIR,
  "sessionInfo.txt"
)


# ==============================================================================
# 4. Logging
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


cat("============================================================\n")
cat("Stage 2 epithelial BayesPrism\n")
cat("S2-B2b: formal prism construction\n")
cat("============================================================\n\n")

cat("START:", format(Sys.time()), "\n")

cat(
  "BayesPrism version:",
  as.character(packageVersion("BayesPrism")),
  "\n\n"
)


# ==============================================================================
# 5. Helpers
# ==============================================================================

assert_true <- function(condition, message) {

  if (!isTRUE(condition)) {
    stop(
      message,
      call. = FALSE
    )
  }
}


audit_row <- function(
  item,
  observed,
  expected = NA_character_,
  status
) {

  data.frame(
    item = as.character(item),
    observed = as.character(observed),
    expected = as.character(expected),
    status = as.character(status),
    stringsAsFactors = FALSE
  )
}


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


# ==============================================================================
# 6. Validate required inputs
# ==============================================================================

assert_true(
  file.exists(CLEAN_REFERENCE_RDS),
  paste(
    "Missing cleaned Stage 2 reference:",
    CLEAN_REFERENCE_RDS
  )
)

assert_true(
  file.exists(MIXTURE_RDS),
  paste(
    "Missing Stage 2 pseudo-mixture:",
    MIXTURE_RDS
  )
)

assert_true(
  file.exists(OUTLIER_AUDIT_RDS),
  paste(
    "Missing Stage 2 outlier audit:",
    OUTLIER_AUDIT_RDS
  )
)

cat("[PASS] Required S2 inputs exist\n\n")


# ==============================================================================
# 7. Load cleaned epithelial reference
# ==============================================================================

cat("Loading cleaned epithelial reference...\n")

reference_obj <- readRDS(
  CLEAN_REFERENCE_RDS
)

assert_true(
  is.list(reference_obj),
  "Cleaned reference RDS is not a list."
)

reference <- reference_obj$reference

cell_type_labels <- as.character(
  reference_obj$cell.type.labels
)

cell_state_labels <- as.character(
  reference_obj$cell.state.labels
)


assert_true(
  nrow(reference) == 93386L,
  paste0(
    "Unexpected reference cell count: ",
    nrow(reference)
  )
)

assert_true(
  ncol(reference) == 13744L,
  paste0(
    "Unexpected cleaned reference gene count: ",
    ncol(reference)
  )
)

assert_true(
  length(cell_type_labels) == nrow(reference),
  "cell.type.labels length mismatch."
)

assert_true(
  length(cell_state_labels) == nrow(reference),
  "cell.state.labels length mismatch."
)

assert_true(
  identical(
    cell_type_labels,
    cell_state_labels
  ),
  "Stage 2 requires cell.type.labels == cell.state.labels."
)

assert_true(
  length(unique(cell_type_labels)) == 7L,
  "Reference does not contain exactly seven epithelial identities."
)

cat(
  "[PASS] Cleaned reference:",
  nrow(reference),
  "cells x",
  ncol(reference),
  "genes\n\n"
)


# ==============================================================================
# 8. Load Stage 1 Epithelial Z pseudo-mixture
# ==============================================================================

cat("Loading Stage 1 Epithelial-Z pseudo-mixture...\n")

mixture <- readRDS(
  MIXTURE_RDS
)

assert_true(
  is.matrix(mixture) ||
    inherits(mixture, "Matrix"),
  "Pseudo-mixture is not matrix-like."
)

assert_true(
  nrow(mixture) == 2490L,
  paste0(
    "Unexpected mixture sample count: ",
    nrow(mixture)
  )
)

assert_true(
  ncol(mixture) == 13668L,
  paste0(
    "Unexpected mixture gene count: ",
    ncol(mixture)
  )
)

assert_true(
  all(colnames(mixture) %in% colnames(reference)),
  paste0(
    "One or more Stage 1 Epithelial-Z genes ",
    "are absent from cleaned epithelial reference."
  )
)

assert_true(
  anyDuplicated(rownames(mixture)) == 0L,
  "Duplicated mixture sample IDs."
)

assert_true(
  anyDuplicated(colnames(mixture)) == 0L,
  "Duplicated mixture genes."
)

cat(
  "[PASS] Pseudo-mixture:",
  nrow(mixture),
  "samples x",
  ncol(mixture),
  "genes\n\n"
)


# ==============================================================================
# 9. Confirm S2-B2a decision
# ==============================================================================

outlier_table <- read.delim(
  OUTLIER_AUDIT_RDS,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

assert_true(
  nrow(outlier_table) == 7L,
  paste0(
    "Expected seven default Stage 2 outlier genes, observed ",
    nrow(outlier_table),
    "."
  )
)

expected_default_outliers <- c(
  "FCGBP",
  "CLCA1",
  "ALDOB",
  "SLC26A3",
  "SI",
  "OLFM4",
  "SLC26A2"
)

assert_true(
  setequal(
    outlier_table$gene,
    expected_default_outliers
  ),
  paste0(
    "Default outlier gene set differs from ",
    "the audited S2-B2a result."
  )
)

assert_true(
  all(
    outlier_table$gene %in%
      colnames(reference)
  ),
  "One or more audited outlier genes are absent from the reference."
)

cat(
  "[PASS] Reproduced S2-B2a audited seven-gene ",
  "default outlier set\n\n"
)


# ==============================================================================
# 10. Lock Stage 2 gene universe
# ==============================================================================

# Stage 2 intentionally inherits the Stage 1 Epithelial-Z gene universe.
#
# We do NOT manually subset or alter expression values here.
# new.prism() will perform its standard reference-mixture intersection.
#
# However, its SECOND outlier filtering step is intentionally disabled because
# the pseudo-mixture is conditional epithelial posterior mass, not independent
# whole-biopsy bulk RNA-seq.

EXPECTED_GENE_UNIVERSE <- colnames(
  mixture
)

assert_true(
  length(EXPECTED_GENE_UNIVERSE) == 13668L,
  "Expected Stage 2 gene universe is not 13,668 genes."
)


# ==============================================================================
# 11. new.prism parameters
# ==============================================================================

INPUT_TYPE <- "count.matrix"

KEY <- NULL

PSEUDO_MIN <- 1e-8


# ------------------------------------------------------------------------------
# IMPORTANT PROJECT-SPECIFIC ADAPTATION
#
# Default BayesPrism:
#   outlier.cut      = 0.01
#   outlier.fraction = 0.10
#
# Stage 2:
#   outlier.cut      = 1
#   outlier.fraction = 1
#
# This prevents re-filtering of lineage-informative epithelial genes that
# become relatively abundant after conditioning on Stage 1 Epithelial Z.
# ------------------------------------------------------------------------------

OUTLIER_CUT <- 1

OUTLIER_FRACTION <- 1


cat("new.prism configuration:\n")
cat("  input.type       :", INPUT_TYPE, "\n")
cat("  key              : NULL\n")
cat("  pseudo.min       :", PSEUDO_MIN, "\n")
cat("  outlier.cut      :", OUTLIER_CUT, "\n")
cat("  outlier.fraction :", OUTLIER_FRACTION, "\n\n")

cat(
  "Second-stage mixture outlier filtering is ",
  "INTENTIONALLY DISABLED.\n\n"
)


# ==============================================================================
# 12. Construct formal Stage 2 prism
# ==============================================================================

cat("Constructing Stage 2 prism object...\n\n")

stage2_prism <- new.prism(

  reference = reference,

  input.type = INPUT_TYPE,

  cell.type.labels = cell_type_labels,

  cell.state.labels = cell_state_labels,

  key = KEY,

  mixture = mixture,

  outlier.cut = OUTLIER_CUT,

  outlier.fraction = OUTLIER_FRACTION,

  pseudo.min = PSEUDO_MIN
)

cat("\nnew.prism() completed.\n\n")


# ==============================================================================
# 13. Extract prism components for audit
# ==============================================================================

assert_true(
  inherits(stage2_prism, "prism"),
  "new.prism() did not return a prism object."
)

prism_mixture <- stage2_prism@mixture

phi_state <- stage2_prism@phi_cellState@phi

phi_type <- stage2_prism@phi_cellType@phi

prism_map <- stage2_prism@map


# ==============================================================================
# 14. Expected gene order
# ==============================================================================

# new.prism uses:
#
# gene.shared <- intersect(
#     colnames(reference),
#     colnames(mixture)
# )
#
# Therefore calculate the same expected ordering.

expected_gene_order <- intersect(
  colnames(reference),
  colnames(mixture)
)

assert_true(
  length(expected_gene_order) == 13668L,
  paste0(
    "Expected shared gene count is not 13,668: ",
    length(expected_gene_order)
  )
)


# ==============================================================================
# 15. Formal prism structural audit
# ==============================================================================

assert_true(
  nrow(prism_mixture) == 2490L,
  "Prism mixture sample count changed."
)

assert_true(
  ncol(prism_mixture) == 13668L,
  paste0(
    "Prism gene count is ",
    ncol(prism_mixture),
    ", expected 13668."
  )
)

assert_true(
  nrow(phi_type) == 7L,
  paste0(
    "phi_cellType has ",
    nrow(phi_type),
    " rows; expected 7."
  )
)

assert_true(
  nrow(phi_state) == 7L,
  paste0(
    "phi_cellState has ",
    nrow(phi_state),
    " rows; expected 7."
  )
)

assert_true(
  ncol(phi_type) == 13668L,
  "phi_cellType gene count is not 13,668."
)

assert_true(
  ncol(phi_state) == 13668L,
  "phi_cellState gene count is not 13,668."
)


# ==============================================================================
# 16. Gene-order audit
# ==============================================================================

assert_true(
  identical(
    colnames(prism_mixture),
    expected_gene_order
  ),
  "Prism mixture gene order differs from expected new.prism intersection."
)

assert_true(
  identical(
    colnames(phi_type),
    expected_gene_order
  ),
  "phi_cellType gene order differs from expected gene order."
)

assert_true(
  identical(
    colnames(phi_state),
    expected_gene_order
  ),
  "phi_cellState gene order differs from expected gene order."
)

cat(
  "[PASS] Prism reference and mixture gene ordering is identical\n"
)


# ==============================================================================
# 17. Sample-order audit
# ==============================================================================

assert_true(
  identical(
    rownames(prism_mixture),
    rownames(mixture)
  ),
  "Prism mixture sample order changed."
)

cat(
  "[PASS] Prism sample IDs/order preserved\n"
)


# ==============================================================================
# 18. Verify mixture values were not altered
# ==============================================================================

mixture_expected <- mixture[
  ,
  expected_gene_order,
  drop = FALSE
]

mixture_concordance <- all.equal(
  unname(prism_mixture),
  unname(mixture_expected),
  tolerance = 0
)

assert_true(
  isTRUE(mixture_concordance),
  paste(
    "Prism mixture differs numerically from input pseudo-mixture:",
    mixture_concordance
  )
)

cat(
  "[PASS] Stage 2 pseudo-mixture values preserved exactly\n"
)


# ==============================================================================
# 19. Cell type / state audit
# ==============================================================================

expected_identities <- c(
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Goblet",
  "Enteroendocrine",
  "Tuft",
  "Paneth",
  "Stem_TA_progenitor"
)

assert_true(
  setequal(
    rownames(phi_type),
    expected_identities
  ),
  "phi_cellType identities do not match expected epithelial identities."
)

assert_true(
  setequal(
    rownames(phi_state),
    expected_identities
  ),
  "phi_cellState identities do not match expected epithelial identities."
)


# ==============================================================================
# 20. Validate 1:1 type/state map
# ==============================================================================

assert_true(
  length(prism_map) == 7L,
  paste0(
    "Prism map contains ",
    length(prism_map),
    " entries; expected 7."
  )
)

assert_true(
  setequal(
    names(prism_map),
    expected_identities
  ),
  "Prism map names do not match expected identities."
)

map_1to1 <- all(
  vapply(
    names(prism_map),
    function(x) {

      length(prism_map[[x]]) == 1L &&
        identical(
          as.character(prism_map[[x]]),
          x
        )
    },
    logical(1)
  )
)

assert_true(
  map_1to1,
  "Stage 2 prism type/state mapping is not strictly 1:1."
)

cat(
  "[PASS] Stage 2 cell.type = cell.state mapping is strictly 1:1\n"
)


# ==============================================================================
# 21. phi normalization audit
# ==============================================================================

phi_type_rowsums <- rowSums(
  phi_type
)

phi_state_rowsums <- rowSums(
  phi_state
)

assert_true(
  max(
    abs(phi_type_rowsums - 1)
  ) < 1e-8,
  "phi_cellType rows do not sum to ~1."
)

assert_true(
  max(
    abs(phi_state_rowsums - 1)
  ) < 1e-8,
  "phi_cellState rows do not sum to ~1."
)

cat(
  "[PASS] phi_cellType and phi_cellState rows normalized to 1\n\n"
)


# ==============================================================================
# 22. Save type/state map
# ==============================================================================

map_table <- data.frame(

  cell_type = names(prism_map),

  cell_state = vapply(
    prism_map,
    function(x) {
      paste(
        as.character(x),
        collapse = ";"
      )
    },
    character(1)
  ),

  stringsAsFactors = FALSE
)

write_tsv(
  map_table,
  OUT_MAP
)


# ==============================================================================
# 23. Formal audit table
# ==============================================================================

audit <- do.call(
  rbind,
  list(

    audit_row(
      "reference_cells",
      nrow(reference),
      93386,
      "PASS"
    ),

    audit_row(
      "cleaned_reference_genes",
      ncol(reference),
      13744,
      "PASS"
    ),

    audit_row(
      "input_pseudomixture_samples",
      nrow(mixture),
      2490,
      "PASS"
    ),

    audit_row(
      "input_pseudomixture_genes",
      ncol(mixture),
      13668,
      "PASS"
    ),

    audit_row(
      "default_stage2_outlier_genes_audited",
      nrow(outlier_table),
      7,
      "PASS"
    ),

    audit_row(
      "stage2_outlier_filtering",
      "disabled",
      "disabled",
      "PASS"
    ),

    audit_row(
      "newprism_outlier_cut",
      OUTLIER_CUT,
      1,
      "PASS"
    ),

    audit_row(
      "newprism_outlier_fraction",
      OUTLIER_FRACTION,
      1,
      "PASS"
    ),

    audit_row(
      "prism_samples",
      nrow(prism_mixture),
      2490,
      ifelse(
        nrow(prism_mixture) == 2490L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "prism_genes",
      ncol(prism_mixture),
      13668,
      ifelse(
        ncol(prism_mixture) == 13668L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "phi_celltypes",
      nrow(phi_type),
      7,
      ifelse(
        nrow(phi_type) == 7L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "phi_cellstates",
      nrow(phi_state),
      7,
      ifelse(
        nrow(phi_state) == 7L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "type_state_mapping_1to1",
      map_1to1,
      TRUE,
      ifelse(
        map_1to1,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "sample_order_preserved",
      identical(
        rownames(prism_mixture),
        rownames(mixture)
      ),
      TRUE,
      ifelse(
        identical(
          rownames(prism_mixture),
          rownames(mixture)
        ),
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "gene_order_consistent",
      identical(
        colnames(prism_mixture),
        colnames(phi_type)
      ) &&
        identical(
          colnames(prism_mixture),
          colnames(phi_state)
        ),
      TRUE,
      "PASS"
    ),

    audit_row(
      "mixture_values_preserved",
      isTRUE(mixture_concordance),
      TRUE,
      ifelse(
        isTRUE(mixture_concordance),
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "max_phi_type_rowsum_error",
      max(
        abs(phi_type_rowsums - 1)
      ),
      "<1e-8",
      ifelse(
        max(
          abs(phi_type_rowsums - 1)
        ) < 1e-8,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "max_phi_state_rowsum_error",
      max(
        abs(phi_state_rowsums - 1)
      ),
      "<1e-8",
      ifelse(
        max(
          abs(phi_state_rowsums - 1)
        ) < 1e-8,
        "PASS",
        "FAIL"
      )
    )
  )
)

write_tsv(
  audit,
  OUT_AUDIT
)

cat("Formal prism audit:\n")
print(audit)

assert_true(
  !any(audit$status == "FAIL"),
  "Stage 2 prism construction audit contains FAIL."
)


# ==============================================================================
# 24. Save gene universe
# ==============================================================================

writeLines(
  colnames(prism_mixture),
  OUT_GENE_UNIVERSE
)


# ==============================================================================
# 25. Save formal Stage 2 prism checkpoint
# ==============================================================================

saveRDS(
  stage2_prism,
  OUT_PRISM,
  compress = TRUE
)

cat(
  "\n[PASS] Saved formal Stage 2 prism object:\n",
  OUT_PRISM,
  "\n\n"
)


# ==============================================================================
# 26. Manifest
# ==============================================================================

manifest <- data.frame(

  role = c(
    "SOURCE_CLEANED_REFERENCE",
    "SOURCE_EPITHELIAL_PSEUDOMIXTURE",
    "SOURCE_OUTLIER_AUDIT",
    "STAGE2_PRISM",
    "PRISM_AUDIT",
    "PRISM_GENE_UNIVERSE",
    "TYPE_STATE_MAP",
    "LOG",
    "SESSION_INFO"
  ),

  path = c(
    CLEAN_REFERENCE_RDS,
    MIXTURE_RDS,
    OUTLIER_AUDIT_RDS,
    OUT_PRISM,
    OUT_AUDIT,
    OUT_GENE_UNIVERSE,
    OUT_MAP,
    OUT_LOG,
    OUT_SESSIONINFO
  ),

  stringsAsFactors = FALSE
)

write_tsv(
  manifest,
  OUT_MANIFEST
)


# ==============================================================================
# 27. sessionInfo
# ==============================================================================

capture.output(
  sessionInfo(),
  file = OUT_SESSIONINFO
)


# ==============================================================================
# 28. Final summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("S2-B2b SUMMARY\n")
cat("============================================================\n")

cat(
  "Reference cells                  :",
  nrow(reference),
  "\n"
)

cat(
  "Cleaned reference genes          :",
  ncol(reference),
  "\n"
)

cat(
  "Input pseudo-mixture             :",
  nrow(mixture),
  "samples x",
  ncol(mixture),
  "genes\n"
)

cat(
  "Default Stage2 outliers audited  :",
  nrow(outlier_table),
  "\n"
)

cat(
  "Second-stage outlier filtering   : DISABLED\n"
)

cat(
  "Prism mixture                    :",
  nrow(prism_mixture),
  "samples x",
  ncol(prism_mixture),
  "genes\n"
)

cat(
  "Prism cell types                 :",
  nrow(phi_type),
  "\n"
)

cat(
  "Prism cell states                :",
  nrow(phi_state),
  "\n"
)

cat(
  "Type/state map                   : 1:1\n"
)

cat(
  "Mixture values preserved         :",
  isTRUE(mixture_concordance),
  "\n"
)

cat(
  "Gene universe preserved          :",
  ncol(prism_mixture) == 13668L,
  "\n"
)

cat("\n")

cat(
  "IMPORTANT:\n",
  "Formal Stage 2 prism has been constructed.\n",
  "run.prism() has NOT been executed.\n"
)

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\n")
cat("============================================================\n")
cat("S2-B2b PRISM CONSTRUCTION COMPLETED\n")
cat("END:", format(Sys.time()), "\n")
cat("============================================================\n")