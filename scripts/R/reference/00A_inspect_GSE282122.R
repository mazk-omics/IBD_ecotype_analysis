library(rhdf5)

h5ad_file <- paste0(
    "/home/mazekai/enteric_glia/data_scRNA/GSE282122/",
    "TAURUS_raw_counts_annotated_final.h5ad"
)

cat("===== H5AD ROOT STRUCTURE =====\n\n")

h5_structure <- h5ls(
    h5ad_file,
    recursive = FALSE
)

print(h5_structure)

cat("\n===== FILE INFO =====\n")
print(file.info(h5ad_file)[, c("size", "mtime")])
