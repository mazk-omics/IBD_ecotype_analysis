
suppressPackageStartupMessages({
    library(BayesPrism)
    library(Matrix)
})

options(stringsAsFactors = FALSE)

project_root <- "/home/mazekai/IBD_EcoTyper"

indir <- file.path(
    project_root,
    "03_reference/combined_reference/output/10_bayesprism_smoke_test"
)

input_rds <- file.path(
    indir,
    "bayesprism_smoke_inputs.rds"
)

prism_rds <- file.path(
    indir,
    "bayesprism_smoke_prism_object.rds"
)

bp_rds <- file.path(
    indir,
    "bayesprism_smoke_result.rds"
)

fraction_first_file <- file.path(
    indir,
    "bayesprism_smoke_fraction_first.tsv"
)

fraction_final_file <- file.path(
    indir,
    "bayesprism_smoke_fraction_final.tsv"
)

expression_rds <- file.path(
    indir,
    "bayesprism_smoke_celltype_expression.rds"
)

summary_file <- file.path(
    indir,
    "bayesprism_smoke_run_summary.txt"
)


cat("\n============================================================\n")
cat("BAYESPRISM SMOKE TEST\n")
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
# 1. Load prepared smoke-test input
# ============================================================

cat("\n[1] LOAD SMOKE INPUT\n")

stopifnot(file.exists(input_rds))

obj <- readRDS(input_rds)

required_names <- c(
    "reference",
    "cell.type.labels",
    "cell.state.labels",
    "mixture"
)

if (!all(required_names %in% names(obj))) {
    stop(
        "Smoke-test RDS is missing required objects: ",
        paste(
            setdiff(required_names, names(obj)),
            collapse = ", "
        )
    )
}

reference_raw <- obj$reference
cell_type_labels <- obj$cell.type.labels
cell_state_labels <- obj$cell.state.labels
mixture <- obj$mixture


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

cat(
    "Raw gene intersection:",
    length(
        intersect(
            colnames(reference_raw),
            colnames(mixture)
        )
    ),
    "\n"
)


# ============================================================
# 2. Input guards
# ============================================================

cat("\n[2] INPUT GUARDS\n")

stopifnot(
    inherits(reference_raw, "dgCMatrix"),
    nrow(reference_raw) == length(cell_type_labels),
    nrow(reference_raw) == length(cell_state_labels),
    nrow(reference_raw) == 1000L,
    nrow(mixture) == 10L
)

if (any(!is.finite(reference_raw@x))) {
    stop("Reference contains NA/NaN/Inf.")
}

if (any(reference_raw@x < 0)) {
    stop("Reference contains negative values.")
}

if (any(!is.finite(mixture))) {
    stop("Mixture contains NA/NaN/Inf.")
}

if (any(mixture < 0)) {
    stop("Mixture contains negative values.")
}

if (any(Matrix::rowSums(reference_raw) <= 0)) {
    stop("One or more smoke reference cells have zero total counts.")
}

if (any(rowSums(mixture) <= 0)) {
    stop("One or more mixture samples have zero total counts.")
}

fractional_bulk_n <- sum(
    abs(
        mixture - round(mixture)
    ) > 1e-8
)

cat(
    "Fractional bulk entries:",
    fractional_bulk_n,
    "\n"
)

if (fractional_bulk_n <= 0L) {
    stop(
        "Smoke mixture does not contain fractional values; ",
        "fractional compatibility would not be tested."
    )
}

cat("Input guards: PASS\n")


# ============================================================
# 3. Reference-specific cleanup
#
# Conservative BayesPrism cleanup:
# - ribosomal genes
# - mitochondrial ribosomal genes
# - curated ribosome-related genes
# - mitochondrial chromosome
# - MALAT1
# - sex chromosomes
#
# We intentionally DO NOT remove actin or hemoglobin groups here.
#
# Mixture is not independently transformed. Gene alignment will
# occur inside new.prism().
# ============================================================

cat("\n[3] CLEANUP REFERENCE GENES\n")

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
    "Gene groups to remove:",
    paste(gene_groups_remove, collapse = ", "),
    "\n"
)

reference_clean <- cleanup.genes(
    input = reference_raw,
    input.type = "count.matrix",
    species = "hs",
    gene.group = gene_groups_remove,
    exp.cells = 1
)

reference_clean <- as(
    reference_clean,
    "dgCMatrix"
)

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
    "\n"
)

stopifnot(
    nrow(reference_clean) == nrow(reference_raw),
    ncol(reference_clean) > 1000L
)

if (any(Matrix::rowSums(reference_clean) <= 0)) {
    stop(
        "One or more reference cells became zero-depth after cleanup."
    )
}

shared_after_cleanup <- intersect(
    colnames(reference_clean),
    colnames(mixture)
)

cat(
    "Reference/bulk overlap after cleanup:",
    length(shared_after_cleanup),
    "\n"
)

if (length(shared_after_cleanup) < 1000L) {
    stop(
        "Too few shared genes remain after cleanup."
    )
}

cat("Reference cleanup: PASS\n")


# ============================================================
# 4. Construct prism object
# ============================================================

cat("\n[4] new.prism()\n")

prism_obj <- new.prism(
    reference = reference_clean,
    input.type = "count.matrix",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_state_labels,
    key = NULL,
    mixture = mixture
)

cat(
    "Constructed object class:",
    class(prism_obj),
    "\n"
)

cat(
    "Prism mixture:",
    nrow(prism_obj@mixture),
    "samples x",
    ncol(prism_obj@mixture),
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
    is(prism_obj, "prism"),
    nrow(prism_obj@mixture) == 10L,
    nrow(prism_obj@phi_cellType@phi) == 20L,
    nrow(prism_obj@phi_cellState@phi) == 20L,
    ncol(prism_obj@mixture) ==
        ncol(prism_obj@phi_cellType@phi),
    ncol(prism_obj@mixture) ==
        ncol(prism_obj@phi_cellState@phi)
)

if (!identical(
    colnames(prism_obj@mixture),
    colnames(prism_obj@phi_cellType@phi)
)) {
    stop(
        "Gene ordering differs between prism mixture ",
        "and cell-type reference."
    )
}

if (any(!is.finite(prism_obj@mixture))) {
    stop("Prism mixture contains non-finite values.")
}

if (any(prism_obj@mixture < 0)) {
    stop("Prism mixture contains negative values.")
}

fractional_after_new_prism <- sum(
    abs(
        prism_obj@mixture -
            round(prism_obj@mixture)
    ) > 1e-8
)

cat(
    "Fractional entries retained in prism mixture:",
    fractional_after_new_prism,
    "\n"
)

if (fractional_after_new_prism <= 0L) {
    stop(
        "Fractional entries disappeared before run.prism()."
    )
}

saveRDS(
    prism_obj,
    prism_rds,
    compress = FALSE
)

cat(
    "Saved prism object:",
    prism_rds,
    "\n"
)

cat("new.prism(): PASS\n")


# ============================================================
# 5. Run full BayesPrism smoke test
#
# Keep official Gibbs defaults:
# chain.length = 1000
# burn.in      = 500
# thinning     = 2
#
# Only restrict parallel cores.
# ============================================================

cat("\n[5] run.prism()\n")

n_cores <- suppressWarnings(
    as.integer(
        Sys.getenv(
            "BAYESPRISM_NCORES",
            unset = "2"
        )
    )
)

if (
    length(n_cores) != 1L ||
    is.na(n_cores) ||
    n_cores < 1L
) {
    n_cores <- 2L
}

cat(
    "Using cores:",
    n_cores,
    "\n"
)

start_time <- Sys.time()

bp <- run.prism(
    prism = prism_obj,
    n.cores = n_cores,
    update.gibbs = TRUE
)

end_time <- Sys.time()

elapsed_seconds <- as.numeric(
    difftime(
        end_time,
        start_time,
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
        "run.prism() did not return a BayesPrism object."
    )
}

saveRDS(
    bp,
    bp_rds,
    compress = FALSE
)

cat(
    "Saved BayesPrism object:",
    bp_rds,
    "\n"
)

cat("run.prism(): PASS\n")


# ============================================================
# 6. Extract cell-type fractions
# ============================================================

cat("\n[6] EXTRACT FRACTIONS\n")

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
    "Initial fraction matrix:",
    nrow(theta_first),
    "samples x",
    ncol(theta_first),
    "cell types\n"
)

cat(
    "Final fraction matrix:",
    nrow(theta_final),
    "samples x",
    ncol(theta_final),
    "cell types\n"
)


stopifnot(
    nrow(theta_first) == 10L,
    nrow(theta_final) == 10L,
    ncol(theta_first) == 20L,
    ncol(theta_final) == 20L
)


validate_fraction <- function(x, label) {

    if (any(!is.finite(x))) {
        stop(label, " contains NA/NaN/Inf.")
    }

    if (any(x < -1e-12)) {
        stop(label, " contains negative fractions.")
    }

    rs <- rowSums(x)

    cat(
        label,
        " row-sum range:",
        sprintf("%.10f", min(rs)),
        "-",
        sprintf("%.10f", max(rs)),
        "\n"
    )

    if (max(abs(rs - 1)) > 1e-6) {
        stop(
            label,
            " rows do not sum to 1 within tolerance."
        )
    }
}


validate_fraction(
    theta_first,
    "Initial fractions"
)

validate_fraction(
    theta_final,
    "Final fractions"
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
# 7. Extract sample-specific cell-type expression
#
# get.exp() returns:
# samples x genes
# for the requested cell type.
#
# We extract all 20 types in the smoke test so the downstream
# EcoTyper presorted path is also exercised.
# ============================================================

cat("\n[7] EXTRACT SAMPLE-SPECIFIC CELL-TYPE EXPRESSION\n")

cell_types <- colnames(theta_first)

if (length(cell_types) != 20L) {
    stop(
        "Expected 20 cell types in BayesPrism output."
    )
}

expression_by_celltype <- setNames(
    lapply(
        cell_types,
        function(ct) {

            cat(
                "  extracting:",
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

            if (nrow(x) != 10L) {
                stop(
                    "Unexpected sample count for cell type: ",
                    ct
                )
            }

            if (any(!is.finite(x))) {
                stop(
                    "Non-finite expression values for cell type: ",
                    ct
                )
            }

            if (any(x < -1e-12)) {
                stop(
                    "Negative expression values for cell type: ",
                    ct
                )
            }

            x
        }
    ),
    cell_types
)

expression_dims <- do.call(
    rbind,
    lapply(
        expression_by_celltype,
        dim
    )
)

colnames(expression_dims) <- c(
    "samples",
    "genes"
)

print(expression_dims)

if (
    length(
        unique(
            expression_dims[, "genes"]
        )
    ) != 1L
) {
    stop(
        "Cell-type expression matrices have inconsistent gene counts."
    )
}

saveRDS(
    expression_by_celltype,
    expression_rds,
    compress = FALSE
)

cat(
    "Saved cell-type expression:",
    expression_rds,
    "\n"
)

cat(
    "Sample-specific expression extraction: PASS\n"
)


# ============================================================
# 8. Final smoke-test verdict
# ============================================================

cat("\n[8] FINAL VALIDATION\n")

checks <- c(
    BayesPrism_version_2_2_3 =
        packageVersion("BayesPrism") == "2.2.3",

    fractional_input_present =
        fractional_bulk_n > 0,

    fractional_input_reached_prism =
        fractional_after_new_prism > 0,

    prism_object_valid =
        is(prism_obj, "prism"),

    BayesPrism_run_completed =
        is(bp, "BayesPrism"),

    initial_fraction_valid =
        all(abs(rowSums(theta_first) - 1) <= 1e-6),

    final_fraction_valid =
        all(abs(rowSums(theta_final) - 1) <= 1e-6),

    celltype_expression_valid =
        length(expression_by_celltype) == 20L
)

print(checks)

if (!all(checks)) {
    stop(
        "One or more BayesPrism smoke-test checks failed."
    )
}


summary_lines <- c(
    "BayesPrism smoke-test run",
    "=========================",
    paste0(
        "BayesPrism version: ",
        packageVersion("BayesPrism")
    ),
    paste0(
        "Input reference: ",
        nrow(reference_raw),
        " cells x ",
        ncol(reference_raw),
        " genes"
    ),
    paste0(
        "Reference after cleanup: ",
        nrow(reference_clean),
        " cells x ",
        ncol(reference_clean),
        " genes"
    ),
    paste0(
        "Raw reference/bulk intersection: ",
        length(
            intersect(
                colnames(reference_raw),
                colnames(mixture)
            )
        )
    ),
    paste0(
        "Intersection after reference cleanup: ",
        length(shared_after_cleanup)
    ),
    paste0(
        "Prism gene space after outlier filtering/alignment: ",
        ncol(prism_obj@mixture)
    ),
    paste0(
        "Fractional entries in original smoke mixture: ",
        fractional_bulk_n
    ),
    paste0(
        "Fractional entries retained in prism mixture: ",
        fractional_after_new_prism
    ),
    paste0(
        "Final fraction matrix: ",
        nrow(theta_final),
        " samples x ",
        ncol(theta_final),
        " cell types"
    ),
    paste0(
        "Cell-type expression matrices: ",
        length(expression_by_celltype)
    ),
    paste0(
        "run.prism elapsed seconds: ",
        round(elapsed_seconds, 2)
    ),
    "FINAL STATUS: PASS"
)

writeLines(
    summary_lines,
    summary_file
)


cat("\n============================================================\n")
cat("BAYESPRISM SMOKE TEST: PASS\n")
cat("============================================================\n")

cat(
    paste(
        summary_lines,
        collapse = "\n"
    ),
    "\n"
)

cat("\nOutput files:\n")
cat("  ", prism_rds, "\n")
cat("  ", bp_rds, "\n")
cat("  ", fraction_first_file, "\n")
cat("  ", fraction_final_file, "\n")
cat("  ", expression_rds, "\n")
cat("  ", summary_file, "\n")

cat("\nDONE\n")
