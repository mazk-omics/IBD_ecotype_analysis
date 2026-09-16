suppressPackageStartupMessages({
    library(BayesPrism)
    library(Matrix)
    library(arrow)
})

options(stringsAsFactors = FALSE)

project_root <- "/home/mazekai/IBD_EcoTyper"

# ============================================================
# INPUT / OUTPUT
# ============================================================

reference_rds <- file.path(
    project_root,
    "03_reference/combined_reference/output/13_bayesprism_broad_input",
    "bayesprism_broad14_reference_count_matrix.rds"
)

mixture_rds <- file.path(
    project_root,
    "03_reference/combined_reference/output/11_bayesprism_full_input",
    "bayesprism_full_msccr_mixture.rds"
)

previous_prism_rds <- file.path(
    project_root,
    "03_reference/combined_reference/output/12_bayesprism_full_run",
    "bayesprism_full_prism_object.rds"
)

output_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/14_bayesprism_broad_run"
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

prism_rds <- file.path(
    output_dir,
    "bayesprism_broad14_prism_object.rds"
)

bp_rds <- file.path(
    output_dir,
    "bayesprism_broad14_result.rds"
)

fraction_first_file <- file.path(
    output_dir,
    "bayesprism_broad14_fraction_initial.tsv"
)

fraction_final_file <- file.path(
    output_dir,
    "bayesprism_broad14_fraction_final.tsv"
)

cv_first_file <- file.path(
    output_dir,
    "bayesprism_broad14_fraction_initial_cv.tsv"
)

cv_final_file <- file.path(
    output_dir,
    "bayesprism_broad14_fraction_final_cv.tsv"
)

epi_z_audit_file <- file.path(
    output_dir,
    "Epithelial_Z_initial_fraction_audit.tsv"
)

audit_file <- file.path(
    output_dir,
    "bayesprism_broad14_run_audit.txt"
)

session_file <- file.path(
    output_dir,
    "bayesprism_broad14_run_sessionInfo.txt"
)

cat("\n============================================================\n")
cat("BROAD14 BAYESPRISM DECONVOLUTION\n")
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

expected_types <- c(
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
    "Broad reference identities:",
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

    length(unique(cell_type_labels)) == 14L,

    setequal(
        unique(cell_type_labels),
        expected_types
    ),

    sum(
        cell_type_labels == "Epithelial"
    ) == 93386L
)

cat("\nBroad cell counts:\n")

print(
    sort(
        table(cell_type_labels),
        decreasing = TRUE
    )
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

cat("Formal broad14 input validation: PASS\n")

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
#
# Exact same cleanup settings as the previous successful
# 20-type formal run.
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
# 5. BUILD BROAD14 PRISM OBJECT
#
# STANDARD BayesPrism.
#
# PROJECT-SPECIFIC CHANGE:
# Only taxonomy changed:
# 20 fine types -> 14 broad types.
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
    ) == 14L,

    nrow(
        prism_obj@phi_cellState@phi
    ) == 14L,

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

# ------------------------------------------------------------
# Because only labels changed, the deconvolution gene universe
# should remain identical to the previous successful formal run.
# ------------------------------------------------------------

cat(
    "Final broad14 prism gene space:",
    final_gene_n,
    "\n"
)

stopifnot(
    final_gene_n == 13668L
)

if (file.exists(previous_prism_rds)) {

    cat(
        "Comparing gene universe with previous 20-type formal prism...\n"
    )

    previous_prism <- readRDS(
        previous_prism_rds
    )

    same_gene_space <- identical(
        colnames(prism_obj@mixture),
        colnames(previous_prism@mixture)
    )

    cat(
        "Gene universe identical to previous formal run:",
        same_gene_space,
        "\n"
    )

    if (!same_gene_space) {
        stop(
            "Broad14 gene universe differs from the previous formal run."
        )
    }

    rm(previous_prism)

    gc(verbose = FALSE)
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
# Same settings as previous successful formal run.
#
# Override core count with:
# export BAYESPRISM_NCORES=20
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
        c(2490L, 14L)
    ),

    identical(
        dim(theta_final),
        c(2490L, 14L)
    ),

    setequal(
        colnames(theta_first),
        expected_types
    ),

    setequal(
        colnames(theta_final),
        expected_types
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
# 11. VALIDATE 1:1 CELL STATE -> CELL TYPE MAP
#
# Broad run intentionally has one same-named state per type.
# ============================================================

cat("\n[11] VALIDATE STATE/TYPE MAP\n")

bp_map <- bp@prism@map

one_to_one <- (
    length(bp_map) == 14L &&
    all(lengths(bp_map) == 1L)
)

same_name <- (
    one_to_one &&
    all(
        vapply(
            bp_map,
            function(x) x[[1]],
            character(1)
        ) == names(bp_map)
    )
)

stopifnot(
    one_to_one,
    same_name
)

theta_state_first <-
    bp@posterior.initial.cellState@theta

stopifnot(
    identical(
        dim(theta_state_first),
        c(2490L, 14L)
    )
)

stopifnot(
    isTRUE(
        all.equal(
            theta_state_first,
            theta_first,
            tolerance = 1e-12,
            check.attributes = TRUE
        )
    )
)

cat("One state per type: PASS\n")
cat("Same-name mapping: PASS\n")
cat("Initial state theta == initial type theta: PASS\n")

# ============================================================
# 12. RECOVER FRACTION CV MATRICES
#
# mergeK() does not propagate initial cell-type theta.cv.
# Under strict 1:1 same-name state/type mapping,
# state-level initial CV is the correct corresponding matrix.
# ============================================================

cat("\n[12] RECOVER FRACTION CV MATRICES\n")

theta_first_cv <-
    bp@posterior.initial.cellState@theta.cv

theta_final_cv <-
    bp@posterior.theta_f@theta.cv

stopifnot(
    identical(
        dim(theta_first_cv),
        c(2490L, 14L)
    ),

    identical(
        dim(theta_final_cv),
        c(2490L, 14L)
    ),

    identical(
        rownames(theta_first_cv),
        rownames(theta_first)
    ),

    identical(
        colnames(theta_first_cv),
        colnames(theta_first)
    ),

    identical(
        rownames(theta_final_cv),
        rownames(theta_final)
    ),

    identical(
        colnames(theta_final_cv),
        colnames(theta_final)
    )
)

if (
    any(!is.finite(theta_first_cv)) ||
    any(!is.finite(theta_final_cv))
) {
    stop(
        "Non-finite values detected in CV matrices."
    )
}

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

cat("Fraction CV recovery: PASS\n")

# ============================================================
# 13. EXTRACT INITIAL POSTERIOR CELL-TYPE Z
#
# get.exp() returns posterior.initial.cellType Z.
# This is initial-Gibbs allocated read count,
# NOT final-theta expression.
#
# Save all 14 types for audit/reuse.
# ============================================================

cat("\n[13] EXTRACT SAMPLE-SPECIFIC CELL-TYPE Z\n")

cell_types <- colnames(
    theta_first
)

stopifnot(
    length(cell_types) == 14L
)

final_gene_n <- ncol(
    bp@prism@mixture
)

expression_gene_n <- NA_integer_

epi_z_fraction <- NULL

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

    stopifnot(
        nrow(x) == 2490L
    )

    if (
        any(!is.finite(x)) ||
        any(x < -1e-12)
    ) {
        stop(
            "Invalid Z values for ",
            ct
        )
    }

    if (is.na(expression_gene_n)) {

        expression_gene_n <- ncol(x)

    } else {

        stopifnot(
            ncol(x) == expression_gene_n
        )
    }

    stopifnot(
        ncol(x) == final_gene_n
    )

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

    # --------------------------------------------------------
    # Special audit for Epithelial Z.
    #
    # Z belongs to the initial Gibbs posterior, therefore its
    # read-mass fraction should match theta_first rather than
    # theta_final.
    # --------------------------------------------------------

    if (ct == "Epithelial") {

        stopifnot(
            identical(
                rownames(x),
                rownames(bp@prism@mixture)
            )
        )

        stopifnot(
            identical(
                colnames(x),
                colnames(bp@prism@mixture)
            )
        )

        epi_z_fraction <-
            rowSums(x) /
            rowSums(bp@prism@mixture)

        epi_audit <- data.frame(
            sample_id = rownames(x),
            z_fraction = epi_z_fraction,
            theta_first =
                theta_first[
                    rownames(x),
                    "Epithelial"
                ],
            theta_final =
                theta_final[
                    rownames(x),
                    "Epithelial"
                ],
            check.names = FALSE
        )

        write.table(
            epi_audit,
            epi_z_audit_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        cat(
            "  Epithelial Z vs initial theta Pearson:",
            cor(
                epi_audit$z_fraction,
                epi_audit$theta_first,
                method = "pearson"
            ),
            "\n"
        )

        cat(
            "  Epithelial Z vs initial theta Spearman:",
            cor(
                epi_audit$z_fraction,
                epi_audit$theta_first,
                method = "spearman"
            ),
            "\n"
        )

        cat(
            "  Epithelial Z vs initial theta MAE:",
            mean(
                abs(
                    epi_audit$z_fraction -
                    epi_audit$theta_first
                )
            ),
            "\n"
        )

        cat(
            "  Epithelial Z vs final theta Pearson:",
            cor(
                epi_audit$z_fraction,
                epi_audit$theta_final,
                method = "pearson"
            ),
            "\n"
        )
    }

    rm(x)

    gc(verbose = FALSE)
}

cat(
    "Z matrices:",
    length(cell_types),
    "\n"
)

cat(
    "Genes per Z matrix:",
    expression_gene_n,
    "\n"
)

cat("Sample-specific Z extraction: PASS\n")

# ============================================================
# 14. FRACTION SUMMARY
# ============================================================

cat("\n[14] FRACTION SUMMARY\n")

fraction_summary <- data.frame(
    cell_type = colnames(theta_final),
    initial_mean =
        colMeans(theta_first),
    initial_median =
        apply(
            theta_first,
            2,
            median
        ),
    final_mean =
        colMeans(theta_final),
    final_median =
        apply(
            theta_final,
            2,
            median
        ),
    check.names = FALSE
)

fraction_summary <- fraction_summary[
    order(
        fraction_summary$final_mean,
        decreasing = TRUE
    ),
]

print(
    fraction_summary,
    row.names = FALSE
)

write.table(
    fraction_summary,
    file.path(
        output_dir,
        "bayesprism_broad14_fraction_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 15. FINAL VALIDATION
# ============================================================

cat("\n[15] FINAL VALIDATION\n")

expression_files <- list.files(
    expression_dir,
    pattern = "^bayesprism_expression_.*\\.rds$",
    full.names = TRUE
)

checks <- c(

    BayesPrism_object =
        is(bp, "BayesPrism"),

    bulk_samples =
        nrow(theta_final) == 2490L,

    cell_types =
        ncol(theta_final) == 14L,

    epithelial_present =
        "Epithelial" %in%
        colnames(theta_final),

    initial_fraction_rows_sum_to_1 =
        max(
            abs(
                rowSums(theta_first) - 1
            )
        ) <= 1e-6,

    final_fraction_rows_sum_to_1 =
        max(
            abs(
                rowSums(theta_final) - 1
            )
        ) <= 1e-6,

    state_type_map_one_to_one =
        one_to_one,

    state_type_map_same_name =
        same_name,

    initial_cv =
        identical(
            dim(theta_first_cv),
            c(2490L, 14L)
        ),

    final_cv =
        identical(
            dim(theta_final_cv),
            c(2490L, 14L)
        ),

    expression_matrices =
        length(expression_files) == 14L,

    expression_gene_space =
        expression_gene_n ==
        final_gene_n,

    final_gene_space =
        final_gene_n == 13668L,

    epithelial_Z_audit =
        !is.null(epi_z_fraction)
)

print(checks)

stopifnot(
    all(checks)
)

# ============================================================
# 16. AUDIT
# ============================================================

cat("\n[16] WRITE AUDIT\n")

epi_audit_final <- read.delim(
    epi_z_audit_file,
    check.names = FALSE
)

audit_lines <- c(

    "Broad14 BayesPrism deconvolution",
    "================================",

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
        "Broad cell types: ",
        ncol(theta_final)
    ),

    "Broad epithelial reference cells: 93386",

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

    "Initial fraction matrix: 2490 x 14",

    "Final fraction matrix: 2490 x 14",

    "Initial CV recovery basis: strict 1:1 same-name state-to-type mapping",

    paste0(
        "Sample-specific Z matrices: ",
        length(expression_files)
    ),

    paste0(
        "Z genes per cell type: ",
        expression_gene_n
    ),

    paste0(
        "Epithelial Z vs initial theta Pearson: ",
        cor(
            epi_audit_final$z_fraction,
            epi_audit_final$theta_first,
            method = "pearson"
        )
    ),

    paste0(
        "Epithelial Z vs initial theta Spearman: ",
        cor(
            epi_audit_final$z_fraction,
            epi_audit_final$theta_first,
            method = "spearman"
        )
    ),

    paste0(
        "Epithelial Z vs initial theta MAE: ",
        mean(
            abs(
                epi_audit_final$z_fraction -
                epi_audit_final$theta_first
            )
        )
    ),

    "PROJECT ADAPTATION: 7 epithelial fine identities collapsed to broad Epithelial for first-stage whole-biopsy deconvolution",

    "Epithelial Z is posterior.initial.cellType Z, not final-theta expression",

    "Second-stage epithelial BayesPrism is project-specific and not an official BayesPrism standard workflow",

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
cat("BROAD14 BAYESPRISM RUN: PASS\n")
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
cat("  Z matrices       :", expression_dir, "\n")
cat("  Epithelial audit :", epi_z_audit_file, "\n")
cat("  Audit            :", audit_file, "\n")

cat("\nDONE\n")
