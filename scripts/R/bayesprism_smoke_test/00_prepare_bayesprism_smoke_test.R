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

bulk_file <- file.path(
    project_root,
    "03_reference/combined_reference/output/09_bayesprism_input_harmonization",
    "GSE193677_MSCCR_Biopsy_counts_GRCh37_Ensembl75_gene_symbol_unique.txt.gz"
)

outdir <- file.path(
    project_root,
    "03_reference/combined_reference/output/10_bayesprism_smoke_test"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(
    outdir,
    "bayesprism_smoke_inputs.rds"
)

out_manifest <- file.path(
    outdir,
    "bayesprism_smoke_reference_manifest.tsv"
)

out_summary <- file.path(
    outdir,
    "bayesprism_smoke_input_summary.txt"
)


cat("\n============================================================\n")
cat("BAYESPRISM SMOKE TEST INPUT PREPARATION\n")
cat("============================================================\n")


# ============================================================
# 1. Input existence
# ============================================================

cat("\n[1] INPUT FILE CHECK\n")

stopifnot(file.exists(ref_h5))
stopifnot(file.exists(manifest_file))
stopifnot(file.exists(bulk_file))

cat("Reference HDF5 : PASS\n")
cat("Manifest       : PASS\n")
cat("Bulk mixture   : PASS\n")


# ============================================================
# 2. Read HDF5 metadata only
# ============================================================

cat("\n[2] READ HDF5 METADATA\n")

genes <- as.character(
    h5read(ref_h5, "/raw_counts/genes")
)

barcodes <- as.character(
    h5read(ref_h5, "/raw_counts/barcodes")
)

shape <- as.integer(
    h5read(ref_h5, "/raw_counts/shape")
)

indptr <- as.integer(
    h5read(ref_h5, "/raw_counts/indptr")
)

cat("Genes    :", length(genes), "\n")
cat("Cells    :", length(barcodes), "\n")
cat("Shape    :", paste(shape, collapse = " x "), "\n")
cat("indptr   :", length(indptr), "\n")

stopifnot(length(genes) == 14431L)
stopifnot(length(barcodes) == 323144L)
stopifnot(length(indptr) == length(barcodes) + 1L)

if (!all(shape == c(length(genes), length(barcodes)))) {
    stop(
        "Unexpected HDF5 shape: ",
        paste(shape, collapse = " x ")
    )
}

if (anyDuplicated(genes)) {
    stop("Duplicated genes detected in final reference.")
}

if (anyDuplicated(barcodes)) {
    stop("Duplicated barcodes detected in final reference.")
}

cat("HDF5 metadata validation: PASS\n")


# ============================================================
# 3. Manifest alignment guard
# ============================================================

cat("\n[3] MANIFEST ALIGNMENT\n")

manifest <- as.data.frame(
    read_parquet(
        manifest_file,
        as_data_frame = TRUE
    )
)

stopifnot(nrow(manifest) == length(barcodes))

required_cols <- c(
    "global_cell_id",
    "reference_identity"
)

if (!all(required_cols %in% names(manifest))) {
    stop(
        "Missing required manifest columns: ",
        paste(
            setdiff(required_cols, names(manifest)),
            collapse = ", "
        )
    )
}

if (!identical(
    as.character(manifest$global_cell_id),
    barcodes
)) {
    stop(
        "HDF5 barcode ordering does not match manifest global_cell_id ordering."
    )
}

cat("Manifest rows       :", nrow(manifest), "\n")
cat(
    "Reference identities:",
    length(unique(manifest$reference_identity)),
    "\n"
)
cat("Barcode/order match : PASS\n")


# ============================================================
# 4. Select smoke-test cells
#
# Engineering smoke test only.
# First 50 cells / identity are used intentionally.
# This is NOT a biological re-sampling procedure.
# ============================================================

cat("\n[4] SELECT SMOKE REFERENCE CELLS\n")

cells_per_identity <- 50L

identity_levels <- sort(
    unique(as.character(manifest$reference_identity))
)

selected_idx <- unlist(
    lapply(
        identity_levels,
        function(id) {

            idx <- which(
                manifest$reference_identity == id
            )

            if (length(idx) < cells_per_identity) {
                stop(
                    "Identity ",
                    id,
                    " contains fewer than ",
                    cells_per_identity,
                    " cells."
                )
            }

            head(idx, cells_per_identity)
        }
    ),
    use.names = FALSE
)

selected_idx <- sort(as.integer(selected_idx))

smoke_manifest <- manifest[selected_idx, , drop = FALSE]

count_table <- table(smoke_manifest$reference_identity)

print(count_table)

stopifnot(length(count_table) == 20L)
stopifnot(all(count_table == cells_per_identity))
stopifnot(length(selected_idx) == 1000L)

cat(
    "Selected reference cells:",
    length(selected_idx),
    "\n"
)


# ============================================================
# 5. Determine contiguous HDF5 column runs
# ============================================================

cat("\n[5] BUILD HDF5 COLUMN RUNS\n")

new_run <- c(
    TRUE,
    diff(selected_idx) != 1L
)

run_id <- cumsum(new_run)

runs <- split(
    selected_idx,
    run_id
)

cat(
    "Selected cells:",
    length(selected_idx),
    "\n"
)

cat(
    "Contiguous HDF5 runs:",
    length(runs),
    "\n"
)


# ============================================================
# 6. Extract ONLY selected columns from CSC HDF5
#
# HDF5 layout:
# genes x cells
#
# scipy CSC:
# data
# indices = zero-based gene row indices
# indptr  = column pointers
# ============================================================

cat("\n[6] EXTRACT SELECTED HDF5 COLUMNS\n")

i_list <- vector("list", length(runs))
j_list <- vector("list", length(runs))
x_list <- vector("list", length(runs))

for (k in seq_along(runs)) {

    cols_global <- runs[[k]]

    c_start <- min(cols_global)
    c_end   <- max(cols_global)

    # scipy CSC indptr values are zero-based offsets.
    start0 <- indptr[c_start]
    end0   <- indptr[c_end + 1L]

    nnz_per_col <- diff(
        indptr[c_start:(c_end + 1L)]
    )

    local_cols <- match(
        cols_global,
        selected_idx
    )

    if (end0 > start0) {

        # rhdf5 indexing is one-based.
        h5_positions <- seq.int(
            start0 + 1L,
            end0
        )

        idx0 <- as.integer(
            h5read(
                ref_h5,
                "/raw_counts/indices",
                index = list(h5_positions)
            )
        )

        vals <- as.numeric(
            h5read(
                ref_h5,
                "/raw_counts/data",
                index = list(h5_positions)
            )
        )

        expected_nnz <- sum(nnz_per_col)

        if (length(vals) != expected_nnz) {
            stop(
                "Unexpected nnz count while reading HDF5 run ",
                k
            )
        }

        i_list[[k]] <- idx0 + 1L

        j_list[[k]] <- rep(
            local_cols,
            times = nnz_per_col
        )

        x_list[[k]] <- vals

    } else {

        i_list[[k]] <- integer(0)
        j_list[[k]] <- integer(0)
        x_list[[k]] <- numeric(0)
    }

    if (
        k %% 10L == 0L ||
        k == length(runs)
    ) {
        cat(
            "  processed run",
            k,
            "/",
            length(runs),
            "\n"
        )
    }
}

i_all <- as.integer(
    unlist(i_list, use.names = FALSE)
)

j_all <- as.integer(
    unlist(j_list, use.names = FALSE)
)

x_all <- as.numeric(
    unlist(x_list, use.names = FALSE)
)

cat(
    "Extracted non-zero entries:",
    length(x_all),
    "\n"
)

if (any(!is.finite(x_all))) {
    stop("Reference contains NA/NaN/Inf.")
}

if (any(x_all < 0)) {
    stop("Reference contains negative counts.")
}


# ============================================================
# 7. Construct genes x cells sparse matrix
# ============================================================

cat("\n[7] CONSTRUCT SPARSE REFERENCE\n")

selected_barcodes <- barcodes[selected_idx]

ref_gxc <- sparseMatrix(
    i = i_all,
    j = j_all,
    x = x_all,
    dims = c(
        length(genes),
        length(selected_idx)
    ),
    dimnames = list(
        genes,
        selected_barcodes
    ),
    giveCsparse = TRUE
)

cat(
    "genes x cells:",
    nrow(ref_gxc),
    "x",
    ncol(ref_gxc),
    "\n"
)


# ============================================================
# 8. Transpose to BayesPrism orientation
#
# BayesPrism:
# cells x genes
# ============================================================

ref_smoke <- as(
    t(ref_gxc),
    "dgCMatrix"
)

stopifnot(
    nrow(ref_smoke) == 1000L,
    ncol(ref_smoke) == 14431L
)

if (!identical(
    rownames(ref_smoke),
    as.character(smoke_manifest$global_cell_id)
)) {
    stop("Reference row ordering does not match smoke manifest.")
}

cell_type_labels <- as.character(
    smoke_manifest$reference_identity
)

cell_state_labels <- cell_type_labels

names(cell_type_labels) <- rownames(ref_smoke)
names(cell_state_labels) <- rownames(ref_smoke)

cat(
    "BayesPrism reference:",
    nrow(ref_smoke),
    "cells x",
    ncol(ref_smoke),
    "genes\n"
)

cat(
    "Reference nnz:",
    length(ref_smoke@x),
    "\n"
)

cat(
    "Reference non-zero genes:",
    sum(Matrix::colSums(ref_smoke) > 0),
    "\n"
)

reference_fractional_n <- sum(
    abs(ref_smoke@x - round(ref_smoke@x)) > 1e-8
)

cat(
    "Reference fractional non-zero entries:",
    reference_fractional_n,
    "\n"
)


# ============================================================
# 9. Extract 10 MSCCR smoke-test samples
#
# MSCCR final matrix may have:
#
# header:
#   sample1 sample2 ... sample2490
#
# data:
#   GENE1 value1 value2 ... value2490
#
# i.e. the gene-symbol column has no header field.
# ============================================================

cat("\n[9] EXTRACT SMOKE BULK SAMPLES\n")


# ------------------------------------------------------------
# Read header + first data row only
# ------------------------------------------------------------

con <- gzfile(
    bulk_file,
    open = "rt"
)

first_two_lines <- readLines(
    con,
    n = 2L
)

close(con)

if (length(first_two_lines) != 2L) {
    stop(
        "Could not read header and first data row from bulk file."
    )
}

header_line <- first_two_lines[1]
first_data_line <- first_two_lines[2]


# ------------------------------------------------------------
# Helper: split whitespace-delimited line
# ------------------------------------------------------------

split_ws <- function(x) {

    strsplit(
        trimws(x),
        "[[:space:]]+"
    )[[1]]
}


strip_quotes <- function(x) {

    sub(
        '^"(.*)"$',
        '\\1',
        x
    )
}


header_fields_raw <- split_ws(
    header_line
)

first_data_fields <- split_ws(
    first_data_line
)

n_header_fields <- length(
    header_fields_raw
)

n_data_fields <- length(
    first_data_fields
)


cat(
    "Header fields:",
    n_header_fields,
    "\n"
)

cat(
    "First data-row fields:",
    n_data_fields,
    "\n"
)

cat(
    "Header preview:\n"
)

print(
    head(
        header_fields_raw,
        5L
    )
)


# ------------------------------------------------------------
# Detect matrix layout
# ------------------------------------------------------------

if (
    n_header_fields == 2490L &&
    n_data_fields == 2491L
) {

    bulk_layout <- "unnamed_gene_column"

    sample_names <- strip_quotes(
        header_fields_raw
    )

} else if (
    n_header_fields == 2491L &&
    n_data_fields == 2491L
) {

    bulk_layout <- "explicit_gene_column"

    sample_names <- strip_quotes(
        header_fields_raw[-1L]
    )

} else {

    stop(
        paste0(
            "Unexpected bulk layout. ",
            "Header fields = ",
            n_header_fields,
            "; first data-row fields = ",
            n_data_fields
        )
    )
}


cat(
    "Detected bulk layout:",
    bulk_layout,
    "\n"
)

cat(
    "Detected sample count:",
    length(sample_names),
    "\n"
)

if (length(sample_names) != 2490L) {

    stop(
        "Expected 2490 MSCCR samples, observed ",
        length(sample_names)
    )
}

if (anyDuplicated(sample_names)) {

    stop(
        "Duplicated MSCCR sample names detected."
    )
}


# ------------------------------------------------------------
# Select 10 samples spanning the full cohort
#
# sample ordinal:
#     1 ... 2490
#
# data-row field:
#     gene = field 1
#     sample k = field k + 1
# ------------------------------------------------------------

sample_ord_idx <- unique(
    as.integer(
        round(
            seq(
                from = 1L,
                to = length(sample_names),
                length.out = 10L
            )
        )
    )
)

if (length(sample_ord_idx) != 10L) {

    stop(
        "Failed to select exactly 10 smoke-test samples."
    )
}

selected_sample_names <- sample_names[
    sample_ord_idx
]

data_field_idx <- c(
    1L,
    sample_ord_idx + 1L
)


cat(
    "Selected sample ordinals:",
    paste(
        sample_ord_idx,
        collapse = ", "
    ),
    "\n"
)

cat(
    "Selected data fields:",
    paste(
        data_field_idx,
        collapse = ", "
    ),
    "\n"
)

cat(
    "Selected MSCCR samples:\n"
)

print(
    selected_sample_names
)


# ------------------------------------------------------------
# Stream data rows only.
#
# Header is deliberately discarded because the gene-symbol
# column has no header in this matrix.
#
# awk:
#   field 1 = gene symbol
#   fields 2:2491 = 2490 samples
# ------------------------------------------------------------

awk_fields <- paste(
    paste0(
        "$",
        data_field_idx
    ),
    collapse = ","
)

awk_program <- sprintf(
    'BEGIN{OFS="\\t"} NR>1 {print %s}',
    awk_fields
)

cmd <- sprintf(
    "gzip -dc %s | awk %s",
    shQuote(bulk_file),
    shQuote(awk_program)
)


bulk_df <- read.table(
    pipe(cmd),
    header = FALSE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
)


names(bulk_df) <- c(
    "gene_symbol",
    selected_sample_names
)


cat(
    "Extracted bulk dimensions:",
    nrow(bulk_df),
    "x",
    ncol(bulk_df),
    "\n"
)


# ------------------------------------------------------------
# Structural validation
# ------------------------------------------------------------

if (nrow(bulk_df) != 53326L) {

    stop(
        "Expected 53,326 bulk genes, observed ",
        nrow(bulk_df)
    )
}

if (ncol(bulk_df) != 11L) {

    stop(
        "Expected gene column + 10 samples, observed ",
        ncol(bulk_df)
    )
}


bulk_genes <- as.character(
    bulk_df[[1]]
)

if (anyNA(bulk_genes)) {

    stop(
        "NA gene symbols detected in smoke bulk extraction."
    )
}

if (anyDuplicated(bulk_genes)) {

    stop(
        "Duplicated gene symbols found in final bulk matrix."
    )
}


# ------------------------------------------------------------
# Convert selected samples to numeric matrix
#
# Input:
# genes x samples
#
# BayesPrism:
# samples x genes
# ------------------------------------------------------------

bulk_numeric <- vapply(
    bulk_df[-1L],
    function(x) {

        out <- as.numeric(x)

        if (anyNA(out)) {
            stop(
                "Non-numeric values detected in bulk matrix."
            )
        }

        out
    },
    numeric(nrow(bulk_df))
)


bulk_smoke <- t(
    bulk_numeric
)

rownames(bulk_smoke) <- selected_sample_names

colnames(bulk_smoke) <- bulk_genes

storage.mode(
    bulk_smoke
) <- "double"


cat(
    "BayesPrism mixture:",
    nrow(bulk_smoke),
    "samples x",
    ncol(bulk_smoke),
    "genes\n"
)


# ------------------------------------------------------------
# Numeric validation
# ------------------------------------------------------------

if (any(!is.finite(bulk_smoke))) {

    stop(
        "Bulk smoke matrix contains NA/NaN/Inf."
    )
}

if (any(bulk_smoke < 0)) {

    stop(
        "Bulk smoke matrix contains negative values."
    )
}


fractional_idx <- abs(
    bulk_smoke -
        round(bulk_smoke)
) > 1e-8

fractional_n <- sum(
    fractional_idx
)

fractional_rate <- fractional_n /
    length(bulk_smoke)


cat(
    "Fractional bulk entries:",
    fractional_n,
    "\n"
)

cat(
    "Fractional bulk rate:",
    sprintf(
        "%.6f%%",
        fractional_rate * 100
    ),
    "\n"
)


if (fractional_n == 0L) {

    warning(
        paste0(
            "Selected 10 samples contain no fractional entries. ",
            "General smoke test remains valid, but fractional ",
            "count compatibility would not be directly exercised."
        )
    )

} else {

    cat(
        "Fractional-count compatibility will be directly tested: YES\n"
    )
}


cat(
    "Bulk streaming extraction: PASS\n"
)


# ============================================================
# 10. Gene-space diagnostics
# ============================================================

cat("\n[10] GENE SPACE\n")

shared_raw <- intersect(
    colnames(ref_smoke),
    colnames(bulk_smoke)
)

cat(
    "Reference genes:",
    ncol(ref_smoke),
    "\n"
)

cat(
    "Bulk genes:",
    ncol(bulk_smoke),
    "\n"
)

cat(
    "Raw name intersection:",
    length(shared_raw),
    "\n"
)

if (length(shared_raw) != 14353L) {

    warning(
        "Expected raw-name intersection = 14,353; observed ",
        length(shared_raw)
    )

} else {

    cat(
        "Expected 14,353-gene raw intersection: PASS\n"
    )
}


# ============================================================
# 11. Save smoke-test inputs
# ============================================================

cat("\n[11] SAVE INPUT OBJECT\n")

smoke_object <- list(

    reference = ref_smoke,

    cell.type.labels = cell_type_labels,

    cell.state.labels = cell_state_labels,

    mixture = bulk_smoke,

    reference_manifest = smoke_manifest,

    metadata = list(
        purpose = "BayesPrism compatibility smoke test only",
        reference_cells_per_identity = cells_per_identity,
        n_reference_identities = length(identity_levels),
        n_reference_cells = nrow(ref_smoke),
        n_reference_genes = ncol(ref_smoke),
        n_bulk_samples = nrow(bulk_smoke),
        n_bulk_genes = ncol(bulk_smoke),
        raw_gene_intersection = length(shared_raw),
        bulk_fractional_entries = fractional_n
    )
)

saveRDS(
    smoke_object,
    out_rds,
    compress = FALSE
)

write.table(
    smoke_manifest,
    out_manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 12. Re-open saved RDS
# ============================================================

cat("\n[12] RDS ROUND-TRIP CHECK\n")

check <- readRDS(out_rds)

stopifnot(
    inherits(check$reference, "dgCMatrix"),
    identical(dim(check$reference), dim(ref_smoke)),
    identical(dim(check$mixture), dim(bulk_smoke)),
    identical(
        check$cell.type.labels,
        cell_type_labels
    )
)

cat("RDS round-trip: PASS\n")


# ============================================================
# 13. Write summary
# ============================================================

summary_lines <- c(
    "BayesPrism smoke-test input preparation",
    "=======================================",
    paste0(
        "Reference: ",
        nrow(ref_smoke),
        " cells x ",
        ncol(ref_smoke),
        " genes"
    ),
    paste0(
        "Reference identities: ",
        length(unique(cell_type_labels))
    ),
    paste0(
        "Cells per identity: ",
        cells_per_identity
    ),
    paste0(
        "Reference nnz: ",
        length(ref_smoke@x)
    ),
    paste0(
        "Reference non-zero genes: ",
        sum(Matrix::colSums(ref_smoke) > 0)
    ),
    paste0(
        "Mixture: ",
        nrow(bulk_smoke),
        " samples x ",
        ncol(bulk_smoke),
        " genes"
    ),
    paste0(
        "Fractional bulk entries: ",
        fractional_n
    ),
    paste0(
        "Raw reference/bulk gene intersection: ",
        length(shared_raw)
    ),
    paste0(
        "Input RDS: ",
        out_rds
    )
)

writeLines(
    summary_lines,
    out_summary
)


cat("\n============================================================\n")
cat("SMOKE INPUT PREPARATION COMPLETE\n")
cat("============================================================\n")

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n")

cat("\nOutput files:\n")
cat("  ", out_rds, "\n")
cat("  ", out_manifest, "\n")
cat("  ", out_summary, "\n")

cat("\nDONE\n")
