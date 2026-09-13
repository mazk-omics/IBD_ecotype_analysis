library(rhdf5)
library(data.table)

h5ad_file <- paste0(
  "/home/mazekai/enteric_glia/data_scRNA/GSE282122/",
  "TAURUS_raw_counts_annotated_final.h5ad"
)

out_file <- "GSE282122_celltype_mapping_template.tsv"

# -----------------------------
# Read final_analysis categories
# -----------------------------

categories <- h5read(
  h5ad_file,
  "/obs/final_analysis/categories"
)

codes <- h5read(
  h5ad_file,
  "/obs/final_analysis/codes"
)

stopifnot(length(codes) == 987743)

# AnnData categorical codes are zero-based
valid <- codes >= 0

tab <- table(codes[valid])

mapping <- data.table(
  code = as.integer(names(tab)),
  final_analysis = categories[
    as.integer(names(tab)) + 1L
  ],
  n_cells = as.integer(tab)
)

setorder(mapping, -n_cells)

mapping[, reference_identity := ""]
mapping[, status := ""]
mapping[, action := ""]
mapping[, note := ""]

fwrite(
  mapping,
  out_file,
  sep = "\t",
  quote = FALSE
)

cat("Mapping template generated:\n")
cat(out_file, "\n\n")

cat("Number of categories:", nrow(mapping), "\n")
cat("Mapped cells:", sum(mapping$n_cells), "\n")
cat("Expected cells:", length(codes), "\n")

stopifnot(nrow(mapping) == length(categories))
stopifnot(sum(mapping$n_cells) == length(codes))

cat("\nPASS: source labels extracted directly from H5AD.\n")
