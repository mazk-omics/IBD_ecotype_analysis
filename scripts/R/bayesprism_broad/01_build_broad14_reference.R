suppressPackageStartupMessages({
  library(Matrix)
  library(arrow)
})

options(stringsAsFactors = FALSE)

# ==============================================================================
# IBD ECOTYPER / BAYESPRISM
# Build 14-type broad reference from the validated 20-type cell-level reference
#
# IMPORTANT:
# - No cell is removed.
# - No raw count is modified.
# - No downsampling is performed.
# - Only analysis labels are remapped.
# ==============================================================================

root <- "/home/mazekai/IBD_EcoTyper"

source_reference <- file.path(
  root,
  "03_reference/combined_reference/output/11_bayesprism_full_input",
  "bayesprism_full_reference_count_matrix.rds"
)

mixture_source <- file.path(
  root,
  "03_reference/combined_reference/output/11_bayesprism_full_input",
  "bayesprism_full_msccr_mixture.rds"
)

manifest_path <- file.path(
  root,
  "03_reference/combined_reference/output/05_final_cell_selection",
  "selected_cell_manifest.parquet"
)

out_dir <- file.path(
  root,
  "03_reference/combined_reference/output/13_bayesprism_broad_input"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_reference <- file.path(
  out_dir,
  "bayesprism_broad14_reference_count_matrix.rds"
)

cat("================================================================================\n")
cat("BUILD BAYESPRISM 14-TYPE BROAD REFERENCE\n")
cat("================================================================================\n\n")

cat("[1] Loading validated reference\n")

obj <- readRDS(source_reference)
manifest <- as.data.frame(
  arrow::read_parquet(manifest_path)
)

stopifnot(
  is.list(obj),
  "reference" %in% names(obj),
  "cell.type.labels" %in% names(obj),
  "cell.state.labels" %in% names(obj)
)

ref <- obj$reference

stopifnot(
  inherits(ref, "dgCMatrix"),
  nrow(ref) == 323144L,
  ncol(ref) == 14431L,
  nrow(manifest) == nrow(ref)
)

cat(
  "Reference:",
  nrow(ref),
  "cells x",
  ncol(ref),
  "genes\n"
)

# ==============================================================================
# Structural concordance
# ==============================================================================

cat("\n[2] Structural concordance\n")

cell_order_ok <- identical(
  as.character(rownames(ref)),
  as.character(manifest$global_cell_id)
)

type_label_ok <- identical(
  as.character(obj$cell.type.labels),
  as.character(manifest$identity)
)

state_label_ok <- identical(
  as.character(obj$cell.state.labels),
  as.character(manifest$identity)
)

cat("Cell order identical:", cell_order_ok, "\n")
cat("Type labels identical:", type_label_ok, "\n")
cat("State labels identical:", state_label_ok, "\n")

if (!cell_order_ok || !type_label_ok || !state_label_ok) {
  stop("Reference / manifest structural concordance FAILED.")
}

# ==============================================================================
# Raw-count integrity
# ==============================================================================

cat("\n[3] Raw-count integrity\n")

if (anyNA(ref@x)) {
  stop("Reference contains NA counts.")
}

if (any(ref@x < 0)) {
  stop("Reference contains negative counts.")
}

if (any(ref@x != floor(ref@x))) {
  stop("Reference contains non-integer counts.")
}

cat("NA counts: 0\n")
cat("Negative counts: 0\n")
cat("Non-integer counts: 0\n")

# ==============================================================================
# Define broad taxonomy
# ==============================================================================

cat("\n[4] Constructing broad taxonomy\n")

fine_labels <- as.character(
  obj$cell.type.labels
)

epithelial_fine <- c(
  "Colonic_absorptive",
  "Ileal_absorptive",
  "Enteroendocrine",
  "Goblet",
  "Paneth",
  "Stem_TA_progenitor",
  "Tuft"
)

broad_labels <- fine_labels

broad_labels[
  broad_labels %in% epithelial_fine
] <- "Epithelial"

expected_broad <- c(
  "B",
  "CD4_T",
  "CD8_T",
  "DC",
  "Endothelial",
  "Epithelial",
  "Fibroblast",
  "Glial",
  "ILC",
  "Mast",
  "Monocyte_Macrophage",
  "NK",
  "Pericyte",
  "Plasma"
)

if (!setequal(
  unique(broad_labels),
  expected_broad
)) {
  stop("Unexpected broad taxonomy.")
}

if (length(unique(broad_labels)) != 14L) {
  stop("Broad reference does not contain exactly 14 cell types.")
}

epi_n <- sum(
  broad_labels == "Epithelial"
)

if (epi_n != 93386L) {
  stop(
    paste0(
      "Expected 93,386 epithelial cells; found ",
      epi_n
    )
  )
}

broad_counts <- sort(
  table(broad_labels),
  decreasing = TRUE
)

print(broad_counts)

cat(
  "\nNumber of broad identities:",
  length(broad_counts),
  "\n"
)

cat(
  "Epithelial cells:",
  epi_n,
  "\n"
)

# ==============================================================================
# Construct new BayesPrism input object
# ==============================================================================

cat("\n[5] Building formal broad reference object\n")

broad_obj <- obj

# First-stage deconvolution is deliberately broad.
broad_obj$cell.type.labels <- broad_labels
broad_obj$cell.state.labels <- broad_labels

if (
  is.null(broad_obj$metadata) ||
  !is.list(broad_obj$metadata)
) {
  broad_obj$metadata <- list()
}

broad_obj$metadata$broad14_reconstruction <- list(
  created = format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S %Z"
  ),
  parent_reference = source_reference,
  manifest = manifest_path,
  strategy =
    "14-type broad BayesPrism reference; seven epithelial fine identities collapsed to Epithelial",
  n_cells = nrow(ref),
  n_genes = ncol(ref),
  n_broad_types = 14L,
  epithelial_cells = epi_n,
  epithelial_fine_identities =
    epithelial_fine,
  downstream_note =
    paste(
      "Paneth and Stem_TA_progenitor remain represented",
      "inside broad Epithelial at stage 1.",
      "They are not intended as downstream EcoTyper target populations."
    )
)

# ==============================================================================
# Atomic write
# ==============================================================================

cat("\n[6] Writing formal broad reference\n")

tmp_output <- paste0(
  output_reference,
  ".tmp"
)

if (file.exists(tmp_output)) {
  unlink(tmp_output)
}

saveRDS(
  broad_obj,
  tmp_output
)

if (!file.rename(
  tmp_output,
  output_reference
)) {
  stop("Atomic rename of reference RDS failed.")
}

cat(
  "Saved:",
  output_reference,
  "\n"
)

# ==============================================================================
# Write taxonomy and audit files
# ==============================================================================

mapping <- unique(
  data.frame(
    fine_identity = fine_labels,
    broad_identity = broad_labels
  )
)

mapping <- mapping[
  order(
    mapping$broad_identity,
    mapping$fine_identity
  ),
]

write.table(
  mapping,
  file.path(
    out_dir,
    "broad14_taxonomy_mapping.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

count_df <- data.frame(
  broad_identity = names(broad_counts),
  n_cells = as.integer(broad_counts)
)

write.table(
  count_df,
  file.path(
    out_dir,
    "broad14_cell_counts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  c(
    paste0("source_reference\t", source_reference),
    paste0("mixture_source\t", mixture_source),
    paste0("manifest\t", manifest_path),
    paste0("n_cells\t", nrow(ref)),
    paste0("n_genes\t", ncol(ref)),
    paste0("n_broad_types\t", length(unique(broad_labels))),
    paste0("epithelial_cells\t", epi_n),
    paste0("cell_order_identical\t", cell_order_ok),
    paste0("type_label_source_identical\t", type_label_ok),
    paste0("state_label_source_identical\t", state_label_ok),
    "raw_count_NA\t0",
    "raw_count_negative\t0",
    "raw_count_noninteger\t0",
    "STATUS\tPASS"
  ),
  file.path(
    out_dir,
    "broad14_build_audit.tsv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    out_dir,
    "sessionInfo_build.txt"
  )
)

cat("\n================================================================================\n")
cat("BROAD14 REFERENCE BUILD: PASS\n")
cat("================================================================================\n")
