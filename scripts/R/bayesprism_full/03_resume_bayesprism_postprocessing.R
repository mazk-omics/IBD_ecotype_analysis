
suppressPackageStartupMessages({
    library(BayesPrism)
})

options(stringsAsFactors = FALSE)

project_root <- "/home/mazekai/IBD_EcoTyper"

output_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/12_bayesprism_full_run"
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

expression_dir <- file.path(
    output_dir,
    "celltype_expression"
)

audit_file <- file.path(
    output_dir,
    "bayesprism_full_run_audit.txt"
)

session_file <- file.path(
    output_dir,
    "bayesprism_full_run_sessionInfo.txt"
)

dir.create(
    expression_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


cat("\n============================================================\n")
cat("RESUME BAYESPRISM POST-PROCESSING\n")
cat("============================================================\n")


# ============================================================
# 1. LOAD EXISTING COMPLETE BAYESPRISM RESULT
# ============================================================

cat("\n[1] LOAD EXISTING BAYESPRISM RESULT\n")

stopifnot(file.exists(bp_rds))

bp <- readRDS(bp_rds)

stopifnot(
    is(bp, "BayesPrism")
)

cat("BayesPrism object: PASS\n")


# ============================================================
# 2. RECOVER FRACTIONS
# ============================================================

cat("\n[2] RECOVER FRACTIONS\n")

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

stopifnot(
    max(abs(rowSums(theta_first) - 1)) <= 1e-6,
    max(abs(rowSums(theta_final) - 1)) <= 1e-6
)

cat(
    "Initial theta:",
    nrow(theta_first), "x", ncol(theta_first), "\n"
)

cat(
    "Final theta:",
    nrow(theta_final), "x", ncol(theta_final), "\n"
)

cat("Fraction recovery: PASS\n")


# ============================================================
# 3. VALIDATE 1:1 CELL STATE -> CELL TYPE MAP
# ============================================================

cat("\n[3] VALIDATE STATE/TYPE MAP\n")

bp_map <- bp@prism@map

one_to_one <- (
    length(bp_map) == 20L &&
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
        c(2490L, 20L)
    )
)

# Because every type contains exactly one same-named state,
# initial cell-state theta must be identical to initial cell-type theta.

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
# 4. RECOVER CV MATRICES
#
# Important BayesPrism implementation detail:
# posterior.initial.cellType is created by mergeK().
# mergeK() does not propagate theta.cv.
#
# Here state/type mapping is strictly 1:1 and same-name.
# Therefore initial state-level theta.cv is exactly the
# corresponding initial type-level theta.cv.
# ============================================================

cat("\n[4] RECOVER FRACTION CV MATRICES\n")

theta_first_cv <-
    bp@posterior.initial.cellState@theta.cv

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

stopifnot(
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
    stop("Non-finite values detected in CV matrices.")
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

cat(
    "Initial CV:",
    nrow(theta_first_cv),
    "x",
    ncol(theta_first_cv),
    "\n"
)

cat(
    "Final CV:",
    nrow(theta_final_cv),
    "x",
    ncol(theta_final_cv),
    "\n"
)

cat("Fraction CV recovery: PASS\n")


# ============================================================
# 5. CONFIRM FRACTION FILES ALREADY EXIST
# ============================================================

cat("\n[5] CHECK EXISTING FRACTION OUTPUTS\n")

stopifnot(
    file.exists(fraction_first_file),
    file.exists(fraction_final_file)
)

cat("Existing fraction files: PASS\n")


# ============================================================
# 6. EXTRACT SAMPLE-SPECIFIC CELL-TYPE EXPRESSION
# ============================================================

cat("\n[6] EXTRACT SAMPLE-SPECIFIC CELL-TYPE EXPRESSION\n")

cell_types <- colnames(theta_first)

stopifnot(
    length(cell_types) == 20L
)

final_gene_n <- ncol(
    bp@prism@mixture
)

expression_gene_n <- NA_integer_

for (ct in cell_types) {

    cat("Extracting:", ct, "\n")

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
            "Invalid expression values for ",
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

    rm(x)

    gc(verbose = FALSE)
}

cat(
    "Expression matrices:",
    length(cell_types),
    "\n"
)

cat(
    "Genes per matrix:",
    expression_gene_n,
    "\n"
)

cat("Sample-specific expression extraction: PASS\n")


# ============================================================
# 7. FINAL VALIDATION
# ============================================================

cat("\n[7] FINAL VALIDATION\n")

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
        ncol(theta_final) == 20L,

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
            c(2490L, 20L)
        ),

    final_cv =
        identical(
            dim(theta_final_cv),
            c(2490L, 20L)
        ),

    expression_matrices =
        length(expression_files) == 20L,

    expression_gene_space =
        expression_gene_n == final_gene_n
)

print(checks)

stopifnot(
    all(checks)
)


# ============================================================
# 8. AUDIT
# ============================================================

cat("\n[8] WRITE AUDIT\n")

audit_lines <- c(

    "Full BayesPrism deconvolution",
    "============================",

    paste0(
        "BayesPrism version: ",
        packageVersion("BayesPrism")
    ),

    "Bulk samples: 2490",
    "Cell types: 20",

    paste0(
        "Final prism gene space: ",
        final_gene_n
    ),

    "Initial Gibbs: PASS",
    "Reference update: PASS",
    "Final Gibbs: PASS",

    "Initial fraction matrix: 2490 x 20",
    "Final fraction matrix: 2490 x 20",

    paste0(
        "Initial cell-type theta.cv slot after mergeK: ",
        paste(
            dim(
                bp@posterior.initial.cellType@theta.cv
            ),
            collapse = " x "
        )
    ),

    paste0(
        "Initial state CV matrix recovered: ",
        paste(
            dim(theta_first_cv),
            collapse = " x "
        )
    ),

    "Initial CV recovery basis: strict 1:1 same-name state-to-type mapping",

    paste0(
        "Final CV matrix: ",
        paste(
            dim(theta_final_cv),
            collapse = " x "
        )
    ),

    paste0(
        "Sample-specific expression matrices: ",
        length(expression_files)
    ),

    paste0(
        "Expression genes per cell type: ",
        expression_gene_n
    ),

    "POST-PROCESSING RECOVERY: PASS",
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
cat("BAYESPRISM POST-PROCESSING RECOVERY: PASS\n")
cat("============================================================\n")

cat(
    paste(
        audit_lines,
        collapse = "\n"
    ),
    "\n"
)

cat("\nDONE\n")

