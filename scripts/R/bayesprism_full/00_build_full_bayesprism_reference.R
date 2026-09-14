
suppressPackageStartupMessages({
    library(rhdf5)
    library(arrow)
    library(Matrix)
})

options(stringsAsFactors = FALSE)

project_root <- "/home/mazekai/IBD_EcoTyper"

ref_h5 <- file.path(
    project_root,
    "03_reference/combined_reference/output/07_unified_reference",
    "unified_reference_raw_counts.h5"
)

manifest_file <- file.path(
    project_root,
    "03_reference/combined_reference/output/07_unified_reference",
    "unified_reference_cell_manifest.parquet"
)

outdir <- file.path(
    project_root,
    "03_reference/combined_reference/output/11_bayesprism_full_input"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

out_rds <- file.path(
    outdir,
    "bayesprism_full_reference_count_matrix.rds"
)

audit_file <- file.path(
    outdir,
    "bayesprism_full_reference_audit.txt"
)

session_file <- file.path(
    outdir,
    "bayesprism_full_reference_sessionInfo.txt"
)


cat("\n============================================================\n")
cat("FULL BAYESPRISM CELL-LEVEL REFERENCE\n")
cat("============================================================\n")


# ============================================================
# 1. File check
# ============================================================

cat("\n[1] INPUT FILE CHECK\n")

stopifnot(
    file.exists(ref_h5),
    file.exists(manifest_file)
)

cat("Reference HDF5 : PASS\n")
cat("Manifest       : PASS\n")


# ============================================================
# 2. Read HDF5 metadata
# ============================================================

cat("\n[2] HDF5 METADATA\n")

genes <- as.character(
    h5read(
        ref_h5,
        "/raw_counts/genes"
    )
)

barcodes <- as.character(
    h5read(
        ref_h5,
        "/raw_counts/barcodes"
    )
)

shape <- as.integer(
    h5read(
        ref_h5,
        "/raw_counts/shape"
    )
)

indptr <- as.integer(
    h5read(
        ref_h5,
        "/raw_counts/indptr"
    )
)

n_genes <- length(genes)
n_cells <- length(barcodes)
expected_nnz <- as.numeric(tail(indptr, 1))


cat("Genes :", n_genes, "\n")
cat("Cells :", n_cells, "\n")
cat("Shape :", paste(shape, collapse = " x "), "\n")
cat(
    "NNZ   :",
    format(
        expected_nnz,
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)


stopifnot(
    n_genes == 14431L,
    n_cells == 323144L,
    identical(
        shape,
        c(14431L, 323144L)
    ),
    length(indptr) == 323145L,
    expected_nnz == 393754562
)

if (anyDuplicated(genes)) {
    stop("Duplicated genes detected.")
}

if (anyDuplicated(barcodes)) {
    stop("Duplicated barcodes detected.")
}

cat("HDF5 metadata: PASS\n")


# ============================================================
# 3. Manifest
# ============================================================

cat("\n[3] MANIFEST ALIGNMENT\n")

manifest <- as.data.frame(
    read_parquet(
        manifest_file,
        as_data_frame = TRUE
    )
)

required_cols <- c(
    "global_cell_id",
    "reference_identity"
)

if (!all(required_cols %in% names(manifest))) {

    stop(
        "Missing manifest columns: ",
        paste(
            setdiff(
                required_cols,
                names(manifest)
            ),
            collapse = ", "
        )
    )
}

if (nrow(manifest) != n_cells) {

    stop(
        "Manifest cell count does not match HDF5."
    )
}

if (!identical(
    as.character(manifest$global_cell_id),
    barcodes
)) {

    stop(
        "Manifest ordering does not match HDF5 barcode ordering."
    )
}

cell_type_labels <- as.character(
    manifest$reference_identity
)

cell_state_labels <- cell_type_labels


if (
    anyNA(cell_type_labels) ||
    any(cell_type_labels == "")
) {

    stop(
        "Missing reference_identity labels."
    )
}

identity_table <- sort(
    table(cell_type_labels),
    decreasing = TRUE
)

print(identity_table)

if (length(identity_table) != 20L) {

    stop(
        "Expected 20 identities, observed ",
        length(identity_table)
    )
}

cat("Manifest ordering : PASS\n")
cat("Identity labels    : PASS\n")


# ============================================================
# 4. Read complete CSC payload
#
# HDF5 currently stores:
# genes x cells
#
# scipy CSC:
# data
# indices = zero-based gene row indices
# indptr  = column pointers
# ============================================================

cat("\n[4] LOAD FULL CSC PAYLOAD\n")

cat("Reading indices...\n")

indices <- as.integer(
    h5read(
        ref_h5,
        "/raw_counts/indices"
    )
)

cat(
    "Indices loaded:",
    format(
        length(indices),
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)


cat("Reading data...\n")

data_values <- as.numeric(
    h5read(
        ref_h5,
        "/raw_counts/data"
    )
)

cat(
    "Data loaded:",
    format(
        length(data_values),
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)


if (
    length(indices) != expected_nnz ||
    length(data_values) != expected_nnz
) {

    stop(
        "CSC payload length does not match indptr."
    )
}

cat("CSC payload length: PASS\n")


# ============================================================
# 5. Numeric validation
# ============================================================

cat("\n[5] RAW COUNT VALIDATION\n")

nonfinite_n <- sum(
    !is.finite(data_values)
)

negative_n <- sum(
    data_values < 0,
    na.rm = TRUE
)

fractional_n <- sum(
    abs(
        data_values -
            round(data_values)
    ) > 1e-8,
    na.rm = TRUE
)

cat(
    "Non-finite entries :",
    nonfinite_n,
    "\n"
)

cat(
    "Negative entries   :",
    negative_n,
    "\n"
)

cat(
    "Fractional entries :",
    fractional_n,
    "\n"
)


if (
    nonfinite_n != 0L ||
    negative_n != 0L ||
    fractional_n != 0L
) {

    stop(
        "Unexpected numeric properties in scRNA raw counts."
    )
}

cat("Raw-count validation: PASS\n")


# ============================================================
# 6. Reconstruct genes x cells dgCMatrix
#
# IMPORTANT:
# scipy CSC and R dgCMatrix use the same conceptual storage:
#
# i = zero-based row indices
# p = zero-based column pointers
# x = values
# ============================================================

cat("\n[6] RECONSTRUCT genes x cells dgCMatrix\n")

ref_gxc <- new(
    "dgCMatrix",
    i = indices,
    p = indptr,
    x = data_values,
    Dim = as.integer(
        c(
            n_genes,
            n_cells
        )
    ),
    Dimnames = list(
        genes,
        barcodes
    )
)

cat(
    "Matrix:",
    nrow(ref_gxc),
    "genes x",
    ncol(ref_gxc),
    "cells\n"
)

cat(
    "NNZ:",
    format(
        length(ref_gxc@x),
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)


stopifnot(
    inherits(ref_gxc, "dgCMatrix"),
    identical(
        dim(ref_gxc),
        c(14431L, 323144L)
    ),
    length(ref_gxc@x) == expected_nnz
)

cat("genes x cells reconstruction: PASS\n")


# ============================================================
# 7. Transpose to BayesPrism orientation
#
# BayesPrism count.matrix:
# cells x genes
# ============================================================

cat("\n[7] TRANSPOSE TO BAYESPRISM ORIENTATION\n")

reference <- as(
    Matrix::t(ref_gxc),
    "dgCMatrix"
)

cat(
    "BayesPrism reference:",
    nrow(reference),
    "cells x",
    ncol(reference),
    "genes\n"
)

cat(
    "Reference class:",
    class(reference),
    "\n"
)


stopifnot(
    inherits(reference, "dgCMatrix"),
    identical(
        dim(reference),
        c(323144L, 14431L)
    )
)

if (!identical(
    rownames(reference),
    barcodes
)) {

    stop(
        "Reference cell ordering changed during transpose."
    )
}

if (!identical(
    colnames(reference),
    genes
)) {

    stop(
        "Reference gene ordering changed during transpose."
    )
}

if (!identical(
    rownames(reference),
    as.character(
        manifest$global_cell_id
    )
)) {

    stop(
        "Final reference rows do not match manifest."
    )
}

cat("Orientation/order validation: PASS\n")


# ============================================================
# 8. Cell/gene sanity checks
# ============================================================

cat("\n[8] MATRIX SANITY CHECKS\n")

zero_cells <- sum(
    Matrix::rowSums(reference) == 0
)

zero_genes <- sum(
    Matrix::colSums(reference) == 0
)

cat(
    "Zero-count cells:",
    zero_cells,
    "\n"
)

cat(
    "Zero-count genes:",
    zero_genes,
    "\n"
)


if (zero_cells != 0L) {
    stop(
        "Zero-count cells detected in final reference."
    )
}

if (zero_genes != 0L) {
    stop(
        "Zero-count genes detected in final reference."
    )
}


cat(
    "Reference total raw counts:",
    format(
        sum(reference@x),
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

cat("Matrix sanity checks: PASS\n")


# ============================================================
# 9. Release redundant genes x cells object
# ============================================================

cat("\n[9] MEMORY CLEANUP\n")

rm(
    ref_gxc,
    indices,
    indptr,
    data_values
)

gc(verbose = FALSE)

cat("Temporary objects released.\n")


# ============================================================
# 10. Package formal BayesPrism reference
# ============================================================

cat("\n[10] BUILD FORMAL REFERENCE OBJECT\n")

formal_reference <- list(

    reference = reference,

    input.type = "count.matrix",

    cell.type.labels = cell_type_labels,

    cell.state.labels = cell_state_labels,

    metadata = list(

        source = basename(ref_h5),

        reference_gene_count = ncol(reference),

        reference_cell_count = nrow(reference),

        reference_identity_count =
            length(
                unique(
                    cell_type_labels
                )
            ),

        storage_class =
            class(reference)[1],

        matrix_orientation =
            "cells_x_genes",

        raw_counts = TRUE,

        downsampled_for_bayesprism = FALSE
    )
)


stopifnot(
    nrow(formal_reference$reference) ==
        length(
            formal_reference$cell.type.labels
        ),

    nrow(formal_reference$reference) ==
        length(
            formal_reference$cell.state.labels
        )
)

cat("Formal object validation: PASS\n")


# ============================================================
# 11. Save
# ============================================================

cat("\n[11] SAVE FORMAL REFERENCE\n")

saveRDS(
    formal_reference,
    out_rds,
    compress = FALSE
)

cat(
    "Saved:",
    out_rds,
    "\n"
)


# ============================================================
# 12. RDS round-trip
#
# Do not compare 393M values element-by-element twice.
# Verify metadata, dimensions, ordering and sparse payload.
# ============================================================

cat("\n[12] RDS ROUND-TRIP CHECK\n")

check <- readRDS(
    out_rds
)

stopifnot(
    inherits(
        check$reference,
        "dgCMatrix"
    ),

    identical(
        dim(check$reference),
        c(323144L, 14431L)
    ),

    length(check$reference@x) ==
        expected_nnz,

    identical(
        rownames(check$reference),
        barcodes
    ),

    identical(
        colnames(check$reference),
        genes
    ),

    identical(
        check$cell.type.labels,
        cell_type_labels
    )
)

cat("RDS round-trip: PASS\n")


# ============================================================
# 13. Audit
# ============================================================

audit_lines <- c(

    "Full BayesPrism cell-level reference",
    "====================================",

    paste0(
        "Reference: ",
        nrow(reference),
        " cells x ",
        ncol(reference),
        " genes"
    ),

    paste0(
        "Sparse class: ",
        class(reference)[1]
    ),

    paste0(
        "NNZ: ",
        expected_nnz
    ),

    paste0(
        "Reference identities: ",
        length(
            unique(
                cell_type_labels
            )
        )
    ),

    paste0(
        "Fractional entries: ",
        fractional_n
    ),

    paste0(
        "Negative entries: ",
        negative_n
    ),

    paste0(
        "Zero-count cells: ",
        zero_cells
    ),

    paste0(
        "Zero-count genes: ",
        zero_genes
    ),

    "input.type: count.matrix",

    "Downsampling for BayesPrism: NO",

    "Manifest ordering: PASS",

    "RDS round-trip: PASS",

    "FINAL STATUS: PASS"
)

writeLines(
    audit_lines,
    audit_file
)

capture.output(
    sessionInfo(),
    file = session_file
)


cat("\n============================================================\n")
cat("FULL BAYESPRISM REFERENCE: PASS\n")
cat("============================================================\n")

cat(
    paste(
        audit_lines,
        collapse = "\n"
    ),
    "\n"
)

cat("\nOutput files:\n")
cat("  ", out_rds, "\n")
cat("  ", audit_file, "\n")
cat("  ", session_file, "\n")

cat("\nDONE\n")
