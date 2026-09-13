library(rhdf5)
library(data.table)

h5ad_file <- paste0(
    "/home/mazekai/enteric_glia/data_scRNA/GSE282122/",
    "TAURUS_raw_counts_annotated_final.h5ad"
)

cat("========================================\n")
cat("GSE282122 METADATA AUDIT\n")
cat("========================================\n")

## --------------------------------------------------
## Helper: inspect AnnData categorical column
## --------------------------------------------------

inspect_categorical <- function(file, column) {

    path <- paste0("/obs/", column)

    cat("\n----------------------------------------\n")
    cat("COLUMN:", column, "\n")
    cat("----------------------------------------\n")

    attrs <- h5readAttributes(file, path)

    cat("Attributes:\n")
    print(attrs)

    structure <- h5ls(file, recursive = TRUE)

    sub <- structure[
        structure$group == path,
        c("name", "otype", "dclass", "dim")
    ]

    cat("\nInternal structure:\n")
    print(sub, row.names = FALSE)

    ## Standard AnnData categorical encoding
    categories_path <- paste0(path, "/categories")
    codes_path      <- paste0(path, "/codes")

    categories <- h5read(file, categories_path)
    codes      <- h5read(file, codes_path)

    ## Python categorical codes are zero-based.
    ## -1 indicates missing.
    valid <- codes >= 0

    counts <- table(codes[valid])

    result <- data.frame(
        code     = as.integer(names(counts)),
        category = categories[
            as.integer(names(counts)) + 1L
        ],
        n_cells  = as.integer(counts)
    )

    result <- result[
        order(result$n_cells, decreasing = TRUE),
    ]

    cat("\nCategory counts:\n")
    print(result, row.names = FALSE)

    cat("\nMissing:", sum(!valid), "\n")
    cat("Total:", length(codes), "\n")

    invisible(result)
}


## --------------------------------------------------
## Annotation candidates
## --------------------------------------------------

annotation_fields <- c(
    "final_analysis",
    "major",
    "minor",
    "bucket",
    "sub_bucket"
)

annotation_results <- lapply(
    annotation_fields,
    function(x) inspect_categorical(h5ad_file, x)
)

names(annotation_results) <- annotation_fields


## --------------------------------------------------
## Biological/sample metadata
## --------------------------------------------------

metadata_fields <- c(
    "Patient",
    "sample_id",
    "Disease",
    "Site",
    "Inflammation",
    "Ileum_vs_Colon"
)

metadata_results <- lapply(
    metadata_fields,
    function(x) inspect_categorical(h5ad_file, x)
)

names(metadata_results) <- metadata_fields


cat("\n========================================\n")
cat("METADATA AUDIT COMPLETE\n")
cat("========================================\n")
