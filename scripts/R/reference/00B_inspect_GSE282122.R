library(rhdf5)

h5ad_file <- paste0(
    "/home/mazekai/enteric_glia/data_scRNA/GSE282122/",
    "TAURUS_raw_counts_annotated_final.h5ad"
)

cat("========================================\n")
cat("GSE282122 H5AD STRUCTURE AUDIT\n")
cat("========================================\n\n")

## --------------------------------------------------
## 1. Basic file information
## --------------------------------------------------

fi <- file.info(h5ad_file)

cat("===== FILE =====\n")
cat("Path:", h5ad_file, "\n")
cat("Size:", round(fi$size / 1024^3, 2), "GiB\n\n")


## --------------------------------------------------
## 2. Read HDF5 object index only
## --------------------------------------------------

h5 <- h5ls(
    h5ad_file,
    recursive = TRUE
)

cat("===== X ATTRIBUTES =====\n")
print(h5readAttributes(h5ad_file, "/X"))

cat("\n===== OBS ATTRIBUTES =====\n")
print(h5readAttributes(h5ad_file, "/obs"))

cat("\n===== VAR ATTRIBUTES =====\n")
print(h5readAttributes(h5ad_file, "/var"))


## --------------------------------------------------
## 3. X internal structure
## --------------------------------------------------

cat("\n===== X INTERNAL STRUCTURE =====\n")

x_structure <- h5[
    h5$group == "/X",
    c("group", "name", "otype", "dclass", "dim")
]

print(x_structure)


## --------------------------------------------------
## 4. obs top-level fields
## --------------------------------------------------

cat("\n===== OBS TOP-LEVEL FIELDS =====\n")

obs_structure <- h5[
    h5$group == "/obs",
    c("name", "otype", "dclass", "dim")
]

print(obs_structure, row.names = FALSE)


## --------------------------------------------------
## 5. var top-level fields
## --------------------------------------------------

cat("\n===== VAR TOP-LEVEL FIELDS =====\n")

var_structure <- h5[
    h5$group == "/var",
    c("name", "otype", "dclass", "dim")
]

print(var_structure, row.names = FALSE)


## --------------------------------------------------
## 6. layers
## --------------------------------------------------

cat("\n===== LAYERS =====\n")

layers_structure <- h5[
    h5$group == "/layers",
    c("name", "otype", "dclass", "dim")
]

if (nrow(layers_structure) == 0) {
    cat("No layers found.\n")
} else {
    print(layers_structure, row.names = FALSE)
}


cat("\n========================================\n")
cat("Audit completed without loading X data.\n")
cat("========================================\n")