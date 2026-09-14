
options(stringsAsFactors = FALSE)

project_root <- "/home/mazekai/IBD_EcoTyper"

bulk_file <- file.path(
    project_root,
    "03_reference/combined_reference/output/09_bayesprism_input_harmonization",
    "GSE193677_MSCCR_Biopsy_counts_GRCh37_Ensembl75_gene_symbol_unique.txt.gz"
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
    "bayesprism_full_msccr_mixture.rds"
)

audit_file <- file.path(
    outdir,
    "bayesprism_full_msccr_mixture_audit.txt"
)

session_file <- file.path(
    outdir,
    "bayesprism_full_msccr_mixture_sessionInfo.txt"
)


cat("\n============================================================\n")
cat("FULL MSCCR BAYESPRISM MIXTURE PREPARATION\n")
cat("============================================================\n")


# ============================================================
# 1. File check
# ============================================================

cat("\n[1] INPUT FILE CHECK\n")

stopifnot(
    file.exists(bulk_file)
)

cat("Bulk matrix exists: PASS\n")
cat("Input:", bulk_file, "\n")


# ============================================================
# 2. Inspect raw file layout
#
# Final MSCCR file layout:
#
# header:
#   sample1 sample2 ... sample2490
#
# data:
#   GENE1 value1 value2 ... value2490
#
# The gene-symbol column has no explicit header.
# ============================================================

cat("\n[2] RAW FILE LAYOUT\n")

con <- gzfile(
    bulk_file,
    open = "rt"
)

first_two <- readLines(
    con,
    n = 2L
)

close(con)

if (length(first_two) != 2L) {
    stop("Could not read bulk header and first data row.")
}

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

header_fields <- split_ws(first_two[1])
data_fields   <- split_ws(first_two[2])

cat(
    "Header fields:",
    length(header_fields),
    "\n"
)

cat(
    "First data-row fields:",
    length(data_fields),
    "\n"
)

if (
    length(header_fields) != 2490L ||
    length(data_fields) != 2491L
) {
    stop(
        "Unexpected final MSCCR file layout."
    )
}

sample_names_expected <- strip_quotes(
    header_fields
)

if (anyDuplicated(sample_names_expected)) {
    stop("Duplicated sample IDs detected in raw header.")
}

cat(
    "Expected samples:",
    length(sample_names_expected),
    "\n"
)

cat("Raw file layout: PASS\n")


# ============================================================
# 3. Read complete matrix
#
# read.table handles the first field of each data row as
# row.names, while the 2490 header values name the 2490
# expression columns.
# ============================================================

cat("\n[3] READ COMPLETE MSCCR MATRIX\n")

cat(
    "Reading 53,326 genes x 2,490 samples...\n"
)

con <- gzfile(
    bulk_file,
    open = "rt"
)

bulk_df <- read.table(
    con,
    header = TRUE,
    sep = "",
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = 1
)

close(con)

cat(
    "Loaded dimensions:",
    nrow(bulk_df),
    "genes x",
    ncol(bulk_df),
    "samples\n"
)


if (
    nrow(bulk_df) != 53326L ||
    ncol(bulk_df) != 2490L
) {
    stop(
        "Unexpected full bulk dimensions."
    )
}


# ============================================================
# 4. Identifier validation
# ============================================================

cat("\n[4] IDENTIFIER VALIDATION\n")

bulk_genes <- rownames(
    bulk_df
)

bulk_samples <- colnames(
    bulk_df
)


if (is.null(bulk_genes)) {
    stop("Bulk gene symbols are missing.")
}

if (is.null(bulk_samples)) {
    stop("Bulk sample IDs are missing.")
}

if (anyNA(bulk_genes) || any(bulk_genes == "")) {
    stop("Missing bulk gene symbols.")
}

if (anyNA(bulk_samples) || any(bulk_samples == "")) {
    stop("Missing bulk sample IDs.")
}

if (anyDuplicated(bulk_genes)) {
    stop("Duplicated gene symbols detected.")
}

if (anyDuplicated(bulk_samples)) {
    stop("Duplicated sample IDs detected.")
}


if (!identical(
    bulk_samples,
    sample_names_expected
)) {
    stop(
        "Parsed sample ordering does not match raw file header."
    )
}


cat(
    "Unique genes   :",
    length(bulk_genes),
    "\n"
)

cat(
    "Unique samples :",
    length(bulk_samples),
    "\n"
)

cat("Identifier/order validation: PASS\n")


# ============================================================
# 5. Convert to numeric matrix
#
# genes x samples
# ============================================================

cat("\n[5] CONVERT TO NUMERIC MATRIX\n")

bulk_gxs <- as.matrix(
    bulk_df
)

storage.mode(
    bulk_gxs
) <- "double"


if (!is.numeric(bulk_gxs)) {
    stop("Bulk payload is not numeric.")
}

cat(
    "Numeric matrix:",
    nrow(bulk_gxs),
    "genes x",
    ncol(bulk_gxs),
    "samples\n"
)

cat("Numeric conversion: PASS\n")


# bulk_df no longer needed
rm(bulk_df)
gc(verbose = FALSE)


# ============================================================
# 6. Numeric integrity
# ============================================================

cat("\n[6] NUMERIC INTEGRITY\n")

nonfinite_n <- sum(
    !is.finite(bulk_gxs)
)

negative_n <- sum(
    bulk_gxs < 0,
    na.rm = TRUE
)

fractional_n <- sum(
    abs(
        bulk_gxs -
            round(bulk_gxs)
    ) > 1e-8,
    na.rm = TRUE
)

fractional_rate <- fractional_n /
    length(bulk_gxs)


cat(
    "Total numeric entries:",
    format(
        length(bulk_gxs),
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

cat(
    "Non-finite entries:",
    nonfinite_n,
    "\n"
)

cat(
    "Negative entries:",
    negative_n,
    "\n"
)

cat(
    "Fractional entries:",
    format(
        fractional_n,
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

cat(
    "Fractional rate:",
    sprintf(
        "%.6f%%",
        fractional_rate * 100
    ),
    "\n"
)


if (nonfinite_n != 0L) {
    stop("Non-finite bulk values detected.")
}

if (negative_n != 0L) {
    stop("Negative bulk values detected.")
}

if (fractional_n == 0L) {
    stop(
        "No fractional values detected; unexpected for MSCCR raw count-space matrix."
    )
}

cat("Numeric integrity: PASS\n")


# ============================================================
# 7. Check sample library sizes
# ============================================================

cat("\n[7] SAMPLE LIBRARY SIZE CHECK\n")

sample_totals <- colSums(
    bulk_gxs
)

if (any(!is.finite(sample_totals))) {
    stop("Invalid sample library size.")
}

if (any(sample_totals <= 0)) {
    stop("One or more MSCCR samples have zero total counts.")
}

cat(
    "Minimum sample total:",
    format(
        min(sample_totals),
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

cat(
    "Median sample total:",
    format(
        median(sample_totals),
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

cat(
    "Maximum sample total:",
    format(
        max(sample_totals),
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

cat("Sample library-size check: PASS\n")


# ============================================================
# 8. Transpose to BayesPrism orientation
#
# BayesPrism mixture:
# samples x genes
# ============================================================

cat("\n[8] TRANSPOSE TO BAYESPRISM ORIENTATION\n")

mixture <- t(
    bulk_gxs
)

storage.mode(
    mixture
) <- "double"


cat(
    "Formal BayesPrism mixture:",
    nrow(mixture),
    "samples x",
    ncol(mixture),
    "genes\n"
)


stopifnot(
    is.matrix(mixture),
    identical(
        dim(mixture),
        c(2490L, 53326L)
    )
)

if (!identical(
    rownames(mixture),
    bulk_samples
)) {
    stop(
        "Sample ordering changed during transpose."
    )
}

if (!identical(
    colnames(mixture),
    bulk_genes
)) {
    stop(
        "Gene ordering changed during transpose."
    )
}

cat("Orientation/order validation: PASS\n")


# original orientation no longer needed
rm(bulk_gxs)
gc(verbose = FALSE)


# ============================================================
# 9. Final BayesPrism input guards
# ============================================================

cat("\n[9] FINAL MIXTURE GUARDS\n")

if (any(!is.finite(mixture))) {
    stop(
        "Formal BayesPrism mixture contains NA/NaN/Inf."
    )
}

if (any(mixture < 0)) {
    stop(
        "Formal BayesPrism mixture contains negative values."
    )
}

if (any(rowSums(mixture) <= 0)) {
    stop(
        "Formal BayesPrism mixture contains zero-depth samples."
    )
}


fractional_final <- sum(
    abs(
        mixture -
            round(mixture)
    ) > 1e-8
)

if (fractional_final != fractional_n) {
    stop(
        "Fractional-value accounting changed during transpose."
    )
}

cat(
    "Fractional values preserved:",
    format(
        fractional_final,
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

cat("Final mixture guards: PASS\n")


# ============================================================
# 10. Package formal mixture
# ============================================================

cat("\n[10] BUILD FORMAL MIXTURE OBJECT\n")

formal_mixture <- list(

    mixture = mixture,

    metadata = list(

        source = basename(
            bulk_file
        ),

        dataset = "GSE193677_MSCCR",

        gene_annotation =
            "GRCh37_Ensembl75_unique_gene_symbol",

        mixture_sample_count =
            nrow(mixture),

        mixture_gene_count =
            ncol(mixture),

        matrix_orientation =
            "samples_x_genes",

        count_space = TRUE,

        rounded = FALSE,

        normalized = FALSE,

        log_transformed = FALSE,

        fractional_entries =
            fractional_final,

        fractional_rate =
            fractional_rate
    )
)

cat("Formal mixture object: PASS\n")


# ============================================================
# 11. Save
# ============================================================

cat("\n[11] SAVE FORMAL MIXTURE\n")

saveRDS(
    formal_mixture,
    out_rds,
    compress = FALSE
)

cat(
    "Saved:",
    out_rds,
    "\n"
)


# ============================================================
# 12. Round-trip validation
# ============================================================

cat("\n[12] RDS ROUND-TRIP CHECK\n")

check <- readRDS(
    out_rds
)

stopifnot(
    is.matrix(
        check$mixture
    ),

    identical(
        dim(check$mixture),
        c(2490L, 53326L)
    ),

    identical(
        rownames(check$mixture),
        bulk_samples
    ),

    identical(
        colnames(check$mixture),
        bulk_genes
    ),

    identical(
        check$mixture,
        mixture
    )
)

cat("RDS exact round-trip: PASS\n")


# ============================================================
# 13. Audit
# ============================================================

cat("\n[13] WRITE AUDIT\n")

audit_lines <- c(

    "Full MSCCR BayesPrism mixture",
    "============================",

    paste0(
        "Mixture: ",
        nrow(mixture),
        " samples x ",
        ncol(mixture),
        " genes"
    ),

    paste0(
        "Unique samples: ",
        length(
            unique(
                rownames(mixture)
            )
        )
    ),

    paste0(
        "Unique genes: ",
        length(
            unique(
                colnames(mixture)
            )
        )
    ),

    paste0(
        "Total entries: ",
        length(mixture)
    ),

    paste0(
        "Fractional entries: ",
        fractional_final
    ),

    paste0(
        "Fractional rate: ",
        sprintf(
            "%.6f%%",
            fractional_rate * 100
        )
    ),

    paste0(
        "Negative entries: ",
        negative_n
    ),

    paste0(
        "Non-finite entries: ",
        nonfinite_n
    ),

    "Rounded: NO",
    "Normalized: NO",
    "Log transformed: NO",
    "Orientation: samples x genes",
    "Sample ordering: PASS",
    "Gene uniqueness: PASS",
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
cat("FULL MSCCR MIXTURE: PASS\n")
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
