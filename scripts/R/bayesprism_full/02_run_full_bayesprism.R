
suppressPackageStartupMessages({
    library(BayesPrism)
    library(Matrix)
    library(arrow)
})

options(stringsAsFactors = FALSE)

project_root <- "/home/mazekai/IBD_EcoTyper"

input_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/11_bayesprism_full_input"
)

output_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/12_bayesprism_full_run"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

expression_dir <- file.path(
    output_dir,
    "celltype_expression"
)

dir.create(
    expression_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


reference_rds <- file.path(
    input_dir,
    "bayesprism_full_reference_count_matrix.rds"
)

mixture_rds <- file.path(
    input_dir,
    "bayesprism_full_msccr_mixture.rds"
)

prism_rds <- file.path(
    output_dir,
    "bayesprism_full_prism_object.rds"
)

bp_rds <- file.path(
    output_dir,
    "bayesprism_full_result.rds"
)

fraction_first_file <- file.path(
    output_dir,
    "bayesprism_fraction_initial.tsv"
)

fraction_final_file <- file.path(
    output_dir,
    "bayesprism_fraction_final.tsv"
)

cv_first_file <- file.path(
    output_dir,
    "bayesprism_fraction_initial_cv.tsv"
)

cv_final_file <- file.path(
    output_dir,
    "bayesprism_fraction_final_cv.tsv"
)

audit_file <- file.path(
    output_dir,
    "bayesprism_full_run_audit.txt"
)

session_file <- file.path(
    output_dir,
    "bayesprism_full_run_sessionInfo.txt"
)


cat("\n============================================================\n")
cat("FULL BAYESPRISM DECONVOLUTION\n")
cat("============================================================\n")

cat(
    "BayesPrism version:",
    as.character(packageVersion("BayesPrism")),
    "\n"
)

cat(
    "R version:",
    R.version.string,
    "\n"
)


# ============================================================
# 1. LOAD FORMAL INPUTS
# ============================================================

cat("\n[1] LOAD FORMAL INPUTS\n")

stopifnot(
    file.exists(reference_rds),
    file.exists(mixture_rds)
)

ref_obj <- readRDS(
    reference_rds
)

mix_obj <- readRDS(
    mixture_rds
)

reference_raw <- ref_obj$reference

cell_type_labels <- as.character(
    ref_obj$cell.type.labels
)

cell_state_labels <- as.character(
    ref_obj$cell.state.labels
)

mixture <- mix_obj$mixture


cat(
    "Reference:",
    nrow(reference_raw),
    "cells x",
    ncol(reference_raw),
    "genes\n"
)

cat(
    "Mixture:",
    nrow(mixture),
    "samples x",
    ncol(mixture),
    "genes\n"
)

cat(
    "Reference identities:",
    length(unique(cell_type_labels)),
    "\n"
)


stopifnot(
    inherits(reference_raw, "dgCMatrix"),

    identical(
        dim(reference_raw),
        c(323144L, 14431L)
    ),

    is.matrix(mixture),

    identical(
        dim(mixture),
        c(2490L, 53326L)
    ),

    length(cell_type_labels) == 323144L,

    length(cell_state_labels) == 323144L,

    identical(
        cell_type_labels,
        cell_state_labels
    ),

    length(unique(cell_type_labels)) == 20L
)


fractional_before <- sum(
    abs(
        mixture - round(mixture)
    ) > 1e-8
)

cat(
    "Fractional bulk entries:",
    format(
        fractional_before,
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

stopifnot(
    fractional_before == 3735371L
)

cat("Formal input validation: PASS\n")


# ============================================================
# 2. RAW GENE SPACE
# ============================================================

cat("\n[2] RAW GENE SPACE\n")

raw_shared <- intersect(
    colnames(reference_raw),
    colnames(mixture)
)

cat(
    "Reference genes:",
    ncol(reference_raw),
    "\n"
)

cat(
    "Bulk genes:",
    ncol(mixture),
    "\n"
)

cat(
    "Raw shared genes:",
    length(raw_shared),
    "\n"
)

stopifnot(
    length(raw_shared) == 14353L
)

cat("Raw 14,353-gene overlap: PASS\n")


# ============================================================
# 3. OFFICIAL BAYESPRISM REFERENCE CLEANUP
# ============================================================

cat("\n[3] cleanup.genes()\n")

gene_groups_remove <- c(
    "Rb",
    "Mrp",
    "other_Rb",
    "chrM",
    "MALAT1",
    "chrX",
    "chrY"
)

cat(
    "Gene groups:",
    paste(
        gene_groups_remove,
        collapse = ", "
    ),
    "\n"
)

cleanup_start <- Sys.time()

reference_clean <- cleanup.genes(
    input = reference_raw,
    input.type = "count.matrix",
    species = "hs",
    gene.group = gene_groups_remove,
    exp.cells = 1
)

cleanup_end <- Sys.time()


if (!inherits(reference_clean, "dgCMatrix")) {

    reference_clean <- as(
        reference_clean,
        "dgCMatrix"
    )
}


genes_before_cleanup <- ncol(
    reference_raw
)

genes_after_cleanup <- ncol(
    reference_clean
)

genes_removed_cleanup <-
    genes_before_cleanup -
    genes_after_cleanup


cat(
    "Reference genes before cleanup:",
    genes_before_cleanup,
    "\n"
)

cat(
    "Reference genes after cleanup:",
    genes_after_cleanup,
    "\n"
)

cat(
    "Genes removed by cleanup:",
    genes_removed_cleanup,
    "\n"
)


clean_shared <- intersect(
    colnames(reference_clean),
    colnames(mixture)
)

cat(
    "Reference/bulk overlap after cleanup:",
    length(clean_shared),
    "\n"
)


if (
    any(!is.finite(reference_clean@x)) ||
    any(reference_clean@x < 0)
) {
    stop(
        "Invalid values in cleaned reference."
    )
}

if (
    any(Matrix::rowSums(reference_clean) <= 0)
) {
    stop(
        "One or more reference cells became zero-depth after cleanup."
    )
}

cat("cleanup.genes(): PASS\n")


# ============================================================
# 4. RELEASE ORIGINAL RAW REFERENCE
# ============================================================

cat("\n[4] RELEASE ORIGINAL REFERENCE\n")

rm(
    reference_raw,
    ref_obj,
    raw_shared
)

gc(verbose = FALSE)

cat("Original unfiltered reference released.\n")


# ============================================================
# 5. BUILD FULL PRISM OBJECT
# ============================================================

cat("\n[5] new.prism()\n")

prism_start <- Sys.time()

prism_obj <- new.prism(
    reference = reference_clean,

    input.type = "count.matrix",

    cell.type.labels =
        cell_type_labels,

    cell.state.labels =
        cell_state_labels,

    key = NULL,

    mixture = mixture,

    outlier.cut = 0.01,

    outlier.fraction = 0.1,

    pseudo.min = 1e-8
)

prism_end <- Sys.time()


if (!is(prism_obj, "prism")) {
    stop(
        "new.prism() failed to produce a prism object."
    )
}


final_gene_n <- ncol(
    prism_obj@mixture
)


cat(
    "Prism mixture:",
    nrow(prism_obj@mixture),
    "samples x",
    final_gene_n,
    "genes\n"
)

cat(
    "Cell-type reference:",
    nrow(prism_obj@phi_cellType@phi),
    "types x",
    ncol(prism_obj@phi_cellType@phi),
    "genes\n"
)

cat(
    "Cell-state reference:",
    nrow(prism_obj@phi_cellState@phi),
    "states x",
    ncol(prism_obj@phi_cellState@phi),
    "genes\n"
)


stopifnot(
    nrow(prism_obj@mixture) == 2490L,

    nrow(
        prism_obj@phi_cellType@phi
    ) == 20L,

    nrow(
        prism_obj@phi_cellState@phi
    ) == 20L,

    ncol(prism_obj@mixture) ==
        ncol(
            prism_obj@phi_cellType@phi
        ),

    ncol(prism_obj@mixture) ==
        ncol(
            prism_obj@phi_cellState@phi
        )
)


if (!identical(
    colnames(prism_obj@mixture),
    colnames(
        prism_obj@phi_cellType@phi
    )
)) {
    stop(
        "Mixture/reference gene alignment failed."
    )
}


fractional_after_prism <- sum(
    abs(
        prism_obj@mixture -
            round(prism_obj@mixture)
    ) > 1e-8
)

cat(
    "Fractional entries retained in prism:",
    format(
        fractional_after_prism,
        big.mark = ",",
        scientific = FALSE
    ),
    "\n"
)

if (fractional_after_prism <= 0L) {
    stop(
        "Fractional mixture values disappeared."
    )
}

cat("new.prism(): PASS\n")


# ============================================================
# 6. SAVE PRISM CHECKPOINT
# ============================================================

cat("\n[6] SAVE PRISM CHECKPOINT\n")

saveRDS(
    prism_obj,
    prism_rds,
    compress = FALSE
)

cat(
    "Saved prism checkpoint:\n",
    prism_rds,
    "\n"
)

cat("Prism checkpoint: PASS\n")


# ============================================================
# 7. RELEASE CELL-LEVEL WORKING OBJECTS
#
# prism_obj now contains collapsed BayesPrism reference
# profiles plus aligned mixture.
# ============================================================

cat("\n[7] MEMORY CLEANUP BEFORE GIBBS\n")

rm(
    reference_clean,
    mixture,
    mix_obj,
    cell_type_labels,
    cell_state_labels
)

gc(verbose = FALSE)

cat("Cell-level working objects released.\n")


# ============================================================
# 8. FULL BAYESPRISM RUN
#
# Official default Gibbs parameters are retained.
#
# Core count can be changed with:
#   export BAYESPRISM_NCORES=32
# ============================================================

cat("\n[8] run.prism()\n")

n_cores <- suppressWarnings(
    as.integer(
        Sys.getenv(
            "BAYESPRISM_NCORES",
            unset = "20"
        )
    )
)

if (
    length(n_cores) != 1L ||
    is.na(n_cores) ||
    n_cores < 1L
) {
    n_cores <- 20L
}

cat(
    "Using cores:",
    n_cores,
    "\n"
)


run_start <- Sys.time()

bp <- run.prism(
    prism = prism_obj,

    n.cores = n_cores,

    update.gibbs = TRUE
)

run_end <- Sys.time()


elapsed_seconds <- as.numeric(
    difftime(
        run_end,
        run_start,
        units = "secs"
    )
)


cat(
    "run.prism elapsed seconds:",
    round(elapsed_seconds, 2),
    "\n"
)

cat(
    "Result class:",
    class(bp),
    "\n"
)


if (!is(bp, "BayesPrism")) {
    stop(
        "run.prism() failed."
    )
}

cat("run.prism(): PASS\n")


# ============================================================
# 9. SAVE COMPLETE BAYESPRISM OBJECT
# ============================================================

cat("\n[9] SAVE BAYESPRISM RESULT\n")

saveRDS(
    bp,
    bp_rds,
    compress = FALSE
)

cat(
    "Saved BayesPrism result:\n",
    bp_rds,
    "\n"
)

cat("BayesPrism result save: PASS\n")


# ============================================================
# 10. EXTRACT INITIAL + FINAL FRACTIONS
# ============================================================

cat("\n[10] EXTRACT CELL-TYPE FRACTIONS\n")

theta_first <- get.fraction(
    bp = bp,
    which.theta = "first",
    state.or.type = "type"
)

theta_final <- get.fraction(
    bp = bp,
    which.theta = "final",
    state.or.type = "type"
)


cat(
    "Initial theta:",
    nrow(theta_first),
    "samples x",
    ncol(theta_first),
    "cell types\n"
)

cat(
    "Final theta:",
    nrow(theta_final),
    "samples x",
    ncol(theta_final),
    "cell types\n"
)


stopifnot(
    identical(
        dim(theta_first),
        c(2490L, 20L)
    ),

    identical(
        dim(theta_final),
        c(2490L, 20L)
    )
)


validate_theta <- function(x, name) {

    if (any(!is.finite(x))) {
        stop(
            name,
            " contains non-finite values."
        )
    }

    if (any(x < -1e-12)) {
        stop(
            name,
            " contains negative values."
        )
    }

    rs <- rowSums(x)

    cat(
        name,
        " row-sum range:",
        sprintf("%.10f", min(rs)),
        "-",
        sprintf("%.10f", max(rs)),
        "\n"
    )

    if (
        max(
            abs(rs - 1)
        ) > 1e-6
    ) {
        stop(
            name,
            " rows do not sum to 1."
        )
    }
}


validate_theta(
    theta_first,
    "Initial theta"
)

validate_theta(
    theta_final,
    "Final theta"
)


write.table(
    data.frame(
        sample_id = rownames(theta_first),
        theta_first,
        check.names = FALSE
    ),
    fraction_first_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    data.frame(
        sample_id = rownames(theta_final),
        theta_final,
        check.names = FALSE
    ),
    fraction_final_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Fraction extraction: PASS\n")


# ============================================================
# 11. SAVE FRACTION CV MATRICES
# ============================================================

cat("\n[11] SAVE FRACTION CV MATRICES\n")

theta_first_cv <-
    bp@posterior.initial.cellType@theta.cv

theta_final_cv <-
    bp@posterior.theta_f@theta.cv


stopifnot(
    identical(
        dim(theta_first_cv),
        c(2490L, 20L)
    ),

    identical(
        dim(theta_final_cv),
        c(2490L, 20L)
    )
)


write.table(
    data.frame(
        sample_id = rownames(theta_first_cv),
        theta_first_cv,
        check.names = FALSE
    ),
    cv_first_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    data.frame(
        sample_id = rownames(theta_final_cv),
        theta_final_cv,
        check.names = FALSE
    ),
    cv_final_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Fraction CV export: PASS\n")


# ============================================================
# 12. EXTRACT SAMPLE-SPECIFIC CELL-TYPE EXPRESSION
#
# IMPORTANT:
# BayesPrism get.exp() extracts sample-specific Z from
# posterior.initial.cellType.
#
# Save each cell type separately to avoid creating a second
# giant 20-element in-memory list.
# ============================================================

cat("\n[12] EXTRACT SAMPLE-SPECIFIC CELL-TYPE EXPRESSION\n")

cell_types <- colnames(
    theta_first
)

stopifnot(
    length(cell_types) == 20L
)


expression_gene_n <- NA_integer_

for (ct in cell_types) {

    cat(
        "Extracting:",
        ct,
        "\n"
    )

    x <- get.exp(
        bp = bp,
        state.or.type = "type",
        cell.name = ct
    )

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }


    if (nrow(x) != 2490L) {
        stop(
            "Unexpected sample number for ",
            ct
        )
    }

    if (
        any(!is.finite(x)) ||
        any(x < -1e-12)
    ) {
        stop(
            "Invalid expression values for ",
            ct
        )
    }


    if (is.na(expression_gene_n)) {

        expression_gene_n <- ncol(x)

    } else if (
        ncol(x) != expression_gene_n
    ) {

        stop(
            "Inconsistent expression gene count for ",
            ct
        )
    }


    safe_ct <- gsub(
        "[^A-Za-z0-9_.-]",
        "_",
        ct
    )

    expression_file <- file.path(
        expression_dir,
        paste0(
            "bayesprism_expression_",
            safe_ct,
            ".rds"
        )
    )

    saveRDS(
        x,
        expression_file,
        compress = FALSE
    )

    cat(
        "  dimensions:",
        nrow(x),
        "x",
        ncol(x),
        "\n"
    )

    rm(x)

    gc(verbose = FALSE)
}


cat(
    "Expression matrices:",
    length(cell_types),
    "\n"
)

cat(
    "Genes per expression matrix:",
    expression_gene_n,
    "\n"
)

cat("Sample-specific expression extraction: PASS\n")


# ============================================================
# 13. FINAL VALIDATION
# ============================================================

cat("\n[13] FINAL VALIDATION\n")

checks <- c(

    BayesPrism_object =
        is(bp, "BayesPrism"),

    bulk_samples =
        nrow(theta_final) == 2490L,

    cell_types =
        ncol(theta_final) == 20L,

    final_fraction_rows_sum_to_1 =
        max(
            abs(
                rowSums(theta_final) - 1
            )
        ) <= 1e-6,

    fractional_bulk_preserved =
        fractional_after_prism > 0,

    expression_gene_space =
        expression_gene_n ==
        final_gene_n
)

print(checks)


if (!all(checks)) {
    stop(
        "One or more final BayesPrism validations failed."
    )
}


# ============================================================
# 14. AUDIT
# ============================================================

cat("\n[14] WRITE AUDIT\n")

audit_lines <- c(

    "Full BayesPrism deconvolution",
    "============================",

    paste0(
        "BayesPrism version: ",
        packageVersion("BayesPrism")
    ),

    "Reference input.type: count.matrix",

    "Reference cells: 323144",

    "Reference genes before cleanup: 14431",

    paste0(
        "Reference genes after cleanup: ",
        genes_after_cleanup
    ),

    paste0(
        "Genes removed by cleanup: ",
        genes_removed_cleanup
    ),

    paste0(
        "Reference/bulk overlap after cleanup: ",
        length(clean_shared)
    ),

    "Bulk samples: 2490",

    "Bulk genes before alignment: 53326",

    paste0(
        "Final prism gene space: ",
        final_gene_n
    ),

    paste0(
        "Fractional bulk entries before new.prism: ",
        fractional_before
    ),

    paste0(
        "Fractional entries retained in prism: ",
        fractional_after_prism
    ),

    paste0(
        "run.prism cores: ",
        n_cores
    ),

    paste0(
        "run.prism elapsed seconds: ",
        round(elapsed_seconds, 2)
    ),

    "Initial Gibbs: PASS",

    "Reference update: PASS",

    "Final Gibbs: PASS",

    "Final fraction matrix: 2490 x 20",

    paste0(
        "Sample-specific expression matrices: ",
        length(cell_types)
    ),

    paste0(
        "Expression genes per cell type: ",
        expression_gene_n
    ),

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
cat("FULL BAYESPRISM RUN: PASS\n")
cat("============================================================\n")

cat(
    paste(
        audit_lines,
        collapse = "\n"
    ),
    "\n"
)

cat("\nMain outputs:\n")
cat("  Prism checkpoint :", prism_rds, "\n")
cat("  BayesPrism result:", bp_rds, "\n")
cat("  Initial fractions:", fraction_first_file, "\n")
cat("  Final fractions  :", fraction_final_file, "\n")
cat("  Expression dir   :", expression_dir, "\n")
cat("  Audit            :", audit_file, "\n")

cat("\nDONE\n")
