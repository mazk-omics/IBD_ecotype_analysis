#!/usr/bin/env Rscript

# ==============================================================================
# IBD EcoTyper
# Stage 2 epithelial-only BayesPrism
#
# S2-B1: epithelial reference cleanup + gene-space audit
#
# INPUTS:
#   S2-A epithelial raw reference
#   S2-A Stage1 Epithelial-Z pseudo-mixture
#
# THIS SCRIPT:
#   1. validates frozen S2-A inputs
#   2. applies the SAME cleanup.genes() policy as Broad14
#   3. audits the cleaned epithelial reference gene universe
#   4. calculates overlap with Stage1 Epithelial Z
#   5. saves the cleaned reference and gene-space audit
#
# THIS SCRIPT DOES NOT:
#   - run new.prism()
#   - run run.prism()
#   - normalize mixture
#   - log-transform mixture
#   - round mixture
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
# 1. Project paths
# ==============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "17_bayesprism_epithelial_input"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "18_bayesprism_epithelial_preprocessing"
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
# 2. Input files
# ==============================================================================

REFERENCE_RDS <- file.path(
  INPUT_DIR,
  "bayesprism_epithelial7_reference_raw.rds"
)

MIXTURE_RDS <- file.path(
  INPUT_DIR,
  "bayesprism_stage2_epithelial_pseudomixture.rds"
)


# ==============================================================================
# 3. Output files
# ==============================================================================

OUT_CLEAN_MATRIX <- file.path(
  OUTPUT_DIR,
  "bayesprism_epithelial7_reference_cleaned_matrix.rds"
)

OUT_CLEAN_OBJECT <- file.path(
  OUTPUT_DIR,
  "bayesprism_epithelial7_reference_cleaned.rds"
)

OUT_REMOVED_GENES <- file.path(
  OUTPUT_DIR,
  "stage2_cleanup_removed_genes.tsv"
)

OUT_GENE_OVERLAP <- file.path(
  OUTPUT_DIR,
  "stage2_postcleanup_gene_overlap.tsv"
)

OUT_GENE_SUMMARY <- file.path(
  OUTPUT_DIR,
  "stage2_postcleanup_gene_summary.tsv"
)

OUT_AUDIT <- file.path(
  OUTPUT_DIR,
  "stage2_cleanup_audit.tsv"
)

OUT_COMMON_GENES <- file.path(
  OUTPUT_DIR,
  "stage2_reference_mixture_common_genes.txt"
)

OUT_REFERENCE_ONLY <- file.path(
  OUTPUT_DIR,
  "stage2_postcleanup_reference_only_genes.txt"
)

OUT_MIXTURE_ONLY <- file.path(
  OUTPUT_DIR,
  "stage2_postcleanup_mixture_only_genes.txt"
)

OUT_MANIFEST <- file.path(
  OUTPUT_DIR,
  "stage2_preprocessing_manifest.tsv"
)

OUT_SESSIONINFO <- file.path(
  OUTPUT_DIR,
  "sessionInfo.txt"
)

OUT_LOG <- file.path(
  OUTPUT_DIR,
  "stage2_cleanup.log"
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
cat("Stage 2 epithelial BayesPrism preprocessing\n")
cat("S2-B1: cleanup + gene-space audit\n")
cat("============================================================\n\n")

cat("START:", format(Sys.time()), "\n\n")

cat(
  "BayesPrism version:",
  as.character(packageVersion("BayesPrism")),
  "\n"
)

cat(
  "Matrix version:",
  as.character(packageVersion("Matrix")),
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


# ==============================================================================
# 6. Check S2-A inputs
# ==============================================================================

assert_true(
  file.exists(REFERENCE_RDS),
  paste(
    "Missing S2-A reference:",
    REFERENCE_RDS
  )
)

assert_true(
  file.exists(MIXTURE_RDS),
  paste(
    "Missing S2-A pseudo-mixture:",
    MIXTURE_RDS
  )
)

cat("[PASS] S2-A inputs exist\n\n")


# ==============================================================================
# 7. Load S2-A frozen reference
# ==============================================================================

cat("Loading S2-A epithelial reference...\n")

reference_obj <- readRDS(
  REFERENCE_RDS
)

assert_true(
  is.list(reference_obj),
  "Stage 2 reference RDS is not a list."
)

required_fields <- c(
  "reference",
  "input.type",
  "cell.type.labels",
  "cell.state.labels",
  "metadata"
)

missing_fields <- setdiff(
  required_fields,
  names(reference_obj)
)

assert_true(
  length(missing_fields) == 0L,
  paste(
    "Missing reference fields:",
    paste(
      missing_fields,
      collapse = ", "
    )
  )
)

reference_raw <- reference_obj$reference

cell_type_labels <- as.character(
  reference_obj$cell.type.labels
)

cell_state_labels <- as.character(
  reference_obj$cell.state.labels
)


# ------------------------------------------------------------------------------
# Locked S2-A expectations
# ------------------------------------------------------------------------------

assert_true(
  nrow(reference_raw) == 93386L,
  paste0(
    "Unexpected epithelial cell count: ",
    nrow(reference_raw)
  )
)

assert_true(
  ncol(reference_raw) == 14431L,
  paste0(
    "Unexpected raw epithelial gene count: ",
    ncol(reference_raw)
  )
)

assert_true(
  length(cell_type_labels) == 93386L,
  "cell.type.labels length mismatch."
)

assert_true(
  length(cell_state_labels) == 93386L,
  "cell.state.labels length mismatch."
)

assert_true(
  identical(
    cell_type_labels,
    cell_state_labels
  ),
  paste0(
    "Stage 2 design requires ",
    "cell.type.labels == cell.state.labels."
  )
)

assert_true(
  length(unique(cell_type_labels)) == 7L,
  "Stage 2 reference does not contain exactly 7 identities."
)

assert_true(
  !is.null(colnames(reference_raw)),
  "Reference has no gene names."
)

assert_true(
  anyDuplicated(colnames(reference_raw)) == 0L,
  "Duplicated raw reference genes found."
)

cat(
  "[PASS] S2-A epithelial reference:",
  nrow(reference_raw),
  "cells x",
  ncol(reference_raw),
  "genes\n\n"
)


# ==============================================================================
# 8. Load frozen Stage 1 Epithelial Z
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
  !is.null(colnames(mixture)),
  "Pseudo-mixture has no gene names."
)

assert_true(
  anyDuplicated(colnames(mixture)) == 0L,
  "Duplicated mixture genes found."
)

cat(
  "[PASS] Stage 1 Epithelial Z:",
  nrow(mixture),
  "samples x",
  ncol(mixture),
  "genes\n\n"
)


# ==============================================================================
# 9. Define EXACT Broad14 cleanup policy
# ==============================================================================

gene_groups_remove <- c(
  "Rb",
  "Mrp",
  "other_Rb",
  "chrM",
  "MALAT1",
  "chrX",
  "chrY"
)

cat("Cleanup policy:\n")

cat(
  paste(
    "  gene groups:",
    paste(
      gene_groups_remove,
      collapse = ", "
    )
  ),
  "\n"
)

cat("  species: hs\n")
cat("  input.type: count.matrix\n")
cat("  exp.cells: 1\n\n")


# ==============================================================================
# 10. Run epithelial-reference cleanup
# ==============================================================================

cat("Running cleanup.genes() on epithelial raw-count reference...\n")
cat("This is NOT a BayesPrism deconvolution run.\n\n")

reference_clean <- cleanup.genes(

  input = reference_raw,

  input.type = "count.matrix",

  species = "hs",

  gene.group = gene_groups_remove,

  exp.cells = 1
)


cat("\ncleanup.genes() completed.\n\n")

cat(
  "Reference genes before cleanup:",
  ncol(reference_raw),
  "\n"
)

cat(
  "Reference genes after cleanup:",
  ncol(reference_clean),
  "\n"
)

cat(
  "Genes removed:",
  ncol(reference_raw) - ncol(reference_clean),
  "\n\n"
)


# ==============================================================================
# 11. Basic cleaned-reference validation
# ==============================================================================

assert_true(
  nrow(reference_clean) == 93386L,
  paste0(
    "Cell count changed after cleanup: ",
    nrow(reference_clean)
  )
)

assert_true(
  ncol(reference_clean) <= ncol(reference_raw),
  "Gene count increased after cleanup, which is unexpected."
)

assert_true(
  !is.null(colnames(reference_clean)),
  "Cleaned reference has no gene names."
)

assert_true(
  anyDuplicated(colnames(reference_clean)) == 0L,
  "Cleaned reference contains duplicated genes."
)

assert_true(
  all(
    colnames(reference_clean) %in%
      colnames(reference_raw)
  ),
  paste0(
    "Cleaned reference contains genes absent ",
    "from the raw reference."
  )
)

assert_true(
  length(cell_type_labels) ==
    nrow(reference_clean),
  "Labels no longer align with cleaned reference rows."
)

cat("[PASS] Cleaned epithelial reference structure\n\n")


# ==============================================================================
# 12. Audit removed genes
# ==============================================================================

raw_genes <- colnames(
  reference_raw
)

clean_genes <- colnames(
  reference_clean
)

mixture_genes <- colnames(
  mixture
)

removed_genes <- setdiff(
  raw_genes,
  clean_genes
)

removed_table <- data.frame(
  gene = removed_genes,
  present_in_stage1_epithelial_Z =
    removed_genes %in% mixture_genes,
  stringsAsFactors = FALSE
)

write_tsv(
  removed_table,
  OUT_REMOVED_GENES
)


# ==============================================================================
# 13. Post-cleanup reference/mixture gene-space audit
# ==============================================================================

common_genes <- intersect(
  clean_genes,
  mixture_genes
)

reference_only_genes <- setdiff(
  clean_genes,
  mixture_genes
)

mixture_only_genes <- setdiff(
  mixture_genes,
  clean_genes
)


cat("Post-cleanup gene-space audit:\n")

cat(
  "  cleaned reference genes :",
  length(clean_genes),
  "\n"
)

cat(
  "  pseudo-mixture genes    :",
  length(mixture_genes),
  "\n"
)

cat(
  "  common genes            :",
  length(common_genes),
  "\n"
)

cat(
  "  reference-only genes    :",
  length(reference_only_genes),
  "\n"
)

cat(
  "  mixture-only genes      :",
  length(mixture_only_genes),
  "\n\n"
)


# ==============================================================================
# 14. Save gene lists
# ==============================================================================

writeLines(
  common_genes,
  OUT_COMMON_GENES
)

writeLines(
  reference_only_genes,
  OUT_REFERENCE_ONLY
)

writeLines(
  mixture_only_genes,
  OUT_MIXTURE_ONLY
)


# ==============================================================================
# 15. Detailed gene overlap table
# ==============================================================================

gene_union <- union(
  clean_genes,
  mixture_genes
)

gene_overlap <- data.frame(

  gene = gene_union,

  in_cleaned_epithelial_reference =
    gene_union %in% clean_genes,

  in_stage1_epithelial_Z =
    gene_union %in% mixture_genes,

  stringsAsFactors = FALSE
)

gene_overlap$category <- ifelse(

  gene_overlap$in_cleaned_epithelial_reference &
    gene_overlap$in_stage1_epithelial_Z,

  "overlap",

  ifelse(
    gene_overlap$in_cleaned_epithelial_reference,
    "reference_only",
    "stage1_Z_only"
  )
)

write_tsv(
  gene_overlap,
  OUT_GENE_OVERLAP
)


# ==============================================================================
# 16. Summary table
# ==============================================================================

gene_summary <- data.frame(

  metric = c(
    "raw_epithelial_reference_genes",
    "cleaned_epithelial_reference_genes",
    "genes_removed_by_cleanup",
    "stage1_epithelial_Z_genes",
    "postcleanup_common_genes",
    "postcleanup_reference_only_genes",
    "postcleanup_stage1_Z_only_genes",
    "removed_genes_present_in_stage1_Z"
  ),

  value = c(
    length(raw_genes),
    length(clean_genes),
    length(removed_genes),
    length(mixture_genes),
    length(common_genes),
    length(reference_only_genes),
    length(mixture_only_genes),
    sum(
      removed_genes %in%
        mixture_genes
    )
  ),

  stringsAsFactors = FALSE
)

write_tsv(
  gene_summary,
  OUT_GENE_SUMMARY
)

print(gene_summary)


# ==============================================================================
# 17. Formal audit
# ==============================================================================

audit <- do.call(
  rbind,
  list(

    audit_row(
      "raw_reference_cells",
      nrow(reference_raw),
      93386,
      ifelse(
        nrow(reference_raw) == 93386L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "raw_reference_genes",
      ncol(reference_raw),
      14431,
      ifelse(
        ncol(reference_raw) == 14431L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "cleaned_reference_cells",
      nrow(reference_clean),
      93386,
      ifelse(
        nrow(reference_clean) == 93386L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "cleaned_reference_genes",
      ncol(reference_clean),
      "<=14431",
      ifelse(
        ncol(reference_clean) <= 14431L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "mixture_samples",
      nrow(mixture),
      2490,
      ifelse(
        nrow(mixture) == 2490L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "mixture_genes",
      ncol(mixture),
      13668,
      ifelse(
        ncol(mixture) == 13668L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "duplicated_cleaned_reference_genes",
      anyDuplicated(clean_genes),
      0,
      ifelse(
        anyDuplicated(clean_genes) == 0L,
        "PASS",
        "FAIL"
      )
    ),

    audit_row(
      "postcleanup_common_genes",
      length(common_genes),
      "informational",
      "INFO"
    ),

    audit_row(
      "postcleanup_reference_only_genes",
      length(reference_only_genes),
      "informational",
      "INFO"
    ),

    audit_row(
      "postcleanup_stage1_Z_only_genes",
      length(mixture_only_genes),
      "expected_possible",
      "INFO"
    ),

    audit_row(
      "cell_type_state_labels_identical",
      identical(
        cell_type_labels,
        cell_state_labels
      ),
      TRUE,
      ifelse(
        identical(
          cell_type_labels,
          cell_state_labels
        ),
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

cat("\nFormal audit:\n")
print(audit)

assert_true(
  !any(audit$status == "FAIL"),
  "S2-B1 audit contains FAIL."
)


# ==============================================================================
# 18. Save cleaned matrix
# ==============================================================================

saveRDS(
  reference_clean,
  OUT_CLEAN_MATRIX,
  compress = TRUE
)


# ==============================================================================
# 19. Save formal cleaned Stage 2 reference object
# ==============================================================================

clean_reference_obj <- list(

  reference = reference_clean,

  input.type = "count.matrix",

  cell.type.labels = cell_type_labels,

  cell.state.labels = cell_state_labels,

  metadata = reference_obj$metadata,

  provenance = list(

    project = "IBD_EcoTyper",

    stage =
      "Stage2_epithelial_BayesPrism_S2B1",

    source =
      REFERENCE_RDS,

    source_cells =
      nrow(reference_raw),

    source_genes =
      ncol(reference_raw),

    cleaned_cells =
      nrow(reference_clean),

    cleaned_genes =
      ncol(reference_clean),

    cleanup_species =
      "hs",

    cleanup_gene_groups =
      gene_groups_remove,

    cleanup_exp_cells =
      1,

    stage1_Z_genes =
      length(mixture_genes),

    common_genes =
      length(common_genes),

    reference_only_genes =
      length(reference_only_genes),

    stage1_Z_only_genes =
      length(mixture_only_genes),

    BayesPrism_version =
      as.character(
        packageVersion("BayesPrism")
      ),

    created_at =
      as.character(Sys.time())
  )
)

saveRDS(
  clean_reference_obj,
  OUT_CLEAN_OBJECT,
  compress = TRUE
)

cat(
  "\n[PASS] Saved cleaned epithelial reference:\n",
  OUT_CLEAN_OBJECT,
  "\n"
)


# ==============================================================================
# 20. Manifest
# ==============================================================================

manifest <- data.frame(

  role = c(
    "SOURCE_STAGE2_RAW_REFERENCE",
    "SOURCE_STAGE2_PSEUDOMIXTURE",
    "CLEANED_REFERENCE_MATRIX",
    "CLEANED_REFERENCE_OBJECT",
    "REMOVED_GENES",
    "COMMON_GENES",
    "REFERENCE_ONLY_GENES",
    "STAGE1_Z_ONLY_GENES",
    "GENE_OVERLAP",
    "GENE_SUMMARY",
    "AUDIT",
    "LOG",
    "SESSION_INFO"
  ),

  path = c(
    REFERENCE_RDS,
    MIXTURE_RDS,
    OUT_CLEAN_MATRIX,
    OUT_CLEAN_OBJECT,
    OUT_REMOVED_GENES,
    OUT_COMMON_GENES,
    OUT_REFERENCE_ONLY,
    OUT_MIXTURE_ONLY,
    OUT_GENE_OVERLAP,
    OUT_GENE_SUMMARY,
    OUT_AUDIT,
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
# 21. sessionInfo
# ==============================================================================

capture.output(
  sessionInfo(),
  file = OUT_SESSIONINFO
)


# ==============================================================================
# 22. Final summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("S2-B1 SUMMARY\n")
cat("============================================================\n")

cat(
  "Raw epithelial reference genes     :",
  length(raw_genes),
  "\n"
)

cat(
  "Cleaned epithelial reference genes :",
  length(clean_genes),
  "\n"
)

cat(
  "Genes removed by cleanup           :",
  length(removed_genes),
  "\n"
)

cat(
  "Stage1 Epithelial Z genes          :",
  length(mixture_genes),
  "\n"
)

cat(
  "Post-cleanup common genes          :",
  length(common_genes),
  "\n"
)

cat(
  "Post-cleanup reference-only genes  :",
  length(reference_only_genes),
  "\n"
)

cat(
  "Post-cleanup Stage1-Z-only genes   :",
  length(mixture_only_genes),
  "\n"
)

cat(
  "Removed genes present in Stage1 Z  :",
  sum(
    removed_genes %in%
      mixture_genes
  ),
  "\n"
)

cat("\n")

cat(
  "IMPORTANT:\n",
  "cleanup.genes() has completed.\n",
  "new.prism() has NOT been run.\n",
  "run.prism() has NOT been run.\n"
)

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\n")
cat("============================================================\n")
cat("S2-B1 PREPROCESSING COMPLETED\n")
cat("END:", format(Sys.time()), "\n")
cat("============================================================\n")