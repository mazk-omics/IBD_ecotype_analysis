#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# ==============================================================================
# Stage 2 epithelial BayesPrism
# QC-02: Tuft reference identifiability + posterior assignment audit
#
# NO model fitting.
# NO modification of existing BayesPrism results.
# ==============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

RUN_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference/combined_reference/output/20_bayesprism_epithelial_full_run"
)

QC_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference/combined_reference/output/21_bayesprism_epithelial_qc"
)

RESULT_RDS <- file.path(
  RUN_DIR,
  "bayesprism_epithelial7_result.rds"
)

TUFT_Z_RDS <- file.path(
  RUN_DIR,
  "celltype_expression",
  "bayesprism_expression_Tuft.rds"
)

MASS_QC_TSV <- file.path(
  QC_DIR,
  "stage2_fractional_mass_audit_by_sample.tsv"
)

dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_THETA <- file.path(QC_DIR, "stage2_component_theta_distribution.tsv")
OUT_COR_ORIG <- file.path(QC_DIR, "stage2_reference_spearman_original.tsv")
OUT_COR_UPDATED <- file.path(QC_DIR, "stage2_reference_spearman_updated.tsv")
OUT_SELF_UPDATE <- file.path(QC_DIR, "stage2_reference_original_vs_updated.tsv")
OUT_TUFT_DRIVERS <- file.path(QC_DIR, "stage2_tuft_Z_driver_genes.tsv")
OUT_TUFT_SPEC <- file.path(QC_DIR, "stage2_tuft_reference_specificity.tsv")
OUT_MARKERS <- file.path(QC_DIR, "stage2_tuft_canonical_marker_audit.tsv")
OUT_SAMPLE <- file.path(QC_DIR, "stage2_tuft_sample_distribution.tsv")
OUT_SUMMARY <- file.path(QC_DIR, "stage2_tuft_diagnostic_summary.tsv")
OUT_LOG <- file.path(QC_DIR, "stage2_tuft_reference_assignment_audit.log")
OUT_SESSION <- file.path(QC_DIR, "tuft_audit_sessionInfo.txt")

log_con <- file(OUT_LOG, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

cosine_similarity <- function(x) {
  nr <- sqrt(rowSums(x^2))
  tcrossprod(x) / outer(nr, nr)
}

cat("============================================================\n")
cat("Stage 2 Tuft reference + posterior assignment audit\n")
cat("============================================================\n\n")
cat("START:", format(Sys.time()), "\n\n")


# ==============================================================================
# 1. Load formal BayesPrism result
# ==============================================================================

assert_true(file.exists(RESULT_RDS), "Stage 2 result RDS missing.")

suppressPackageStartupMessages(
  library(BayesPrism)
)

bp <- readRDS(RESULT_RDS)

assert_true(is(bp, "BayesPrism"), "Object is not BayesPrism.")

phi0 <- bp@prism@phi_cellType@phi

assert_true("Tuft" %in% rownames(phi0), "Tuft missing from reference.")

theta_initial <- bp@posterior.initial.cellType@theta
theta_final <- bp@posterior.theta_f@theta

assert_true(
  identical(colnames(theta_initial), colnames(theta_final)),
  "Initial/final theta component names differ."
)

components <- colnames(theta_final)

cat("[PASS] Loaded BayesPrism result\n")
cat("Components:", paste(components, collapse = ", "), "\n")
cat("Reference:", nrow(phi0), "x", ncol(phi0), "\n\n")


# ==============================================================================
# 2. Updated reference
# ==============================================================================

has_updated_reference <- !is.null(bp@reference.update)

if (has_updated_reference) {
  phi1 <- bp@reference.update@phi

  assert_true(
    identical(dimnames(phi0), dimnames(phi1)),
    "Original and updated reference dimensions/names differ."
  )

  cat("[PASS] Updated reference available\n\n")
} else {
  phi1 <- NULL
  cat("[WARN] Updated reference not available\n\n")
}


# ==============================================================================
# 3. Sample reliability from fractional-mass audit
# ==============================================================================

reliable <- rep(TRUE, nrow(theta_final))
names(reliable) <- rownames(theta_final)

mass_qc <- NULL

if (file.exists(MASS_QC_TSV)) {

  mass_qc <- read.delim(
    MASS_QC_TSV,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  idx <- match(rownames(theta_final), mass_qc$sample_id)

  assert_true(!anyNA(idx), "Some theta samples absent from mass QC.")

  rel_loss <- mass_qc$relative_mass_loss[idx]
  allocated_mass <- mass_qc$posterior_allocated_mass[idx]

  # Diagnostic reliability flag only.
  # Not yet a formal exclusion rule.
  reliable <- allocated_mass > 0 & rel_loss <= 0.01

  names(reliable) <- rownames(theta_final)

  cat("Fractional-mass reliability diagnostics:\n")
  cat("  samples total             :", length(reliable), "\n")
  cat("  allocated mass = 0        :", sum(allocated_mass == 0), "\n")
  cat("  relative loss > 1%        :", sum(rel_loss > 0.01), "\n")
  cat("  relative loss > 5%        :", sum(rel_loss > 0.05), "\n")
  cat("  relative loss > 10%       :", sum(rel_loss > 0.10), "\n")
  cat("  diagnostic reliable <=1%  :", sum(reliable), "\n\n")
}


# ==============================================================================
# 4. Component theta distributions
# ==============================================================================

theta_summary_one <- function(mat, label) {

  do.call(
    rbind,
    lapply(colnames(mat), function(component) {

      x <- mat[, component]

      data.frame(
        stage = label,
        component = component,
        mean = mean(x),
        sd = sd(x),
        min = min(x),
        q01 = unname(quantile(x, 0.01)),
        q05 = unname(quantile(x, 0.05)),
        q25 = unname(quantile(x, 0.25)),
        median = median(x),
        q75 = unname(quantile(x, 0.75)),
        q95 = unname(quantile(x, 0.95)),
        q99 = unname(quantile(x, 0.99)),
        max = max(x),
        stringsAsFactors = FALSE
      )
    })
  )
}

theta_summary <- rbind(
  theta_summary_one(theta_initial, "initial"),
  theta_summary_one(theta_final, "final")
)

write_tsv(theta_summary, OUT_THETA)


# ==============================================================================
# 5. Tuft sample-wise distribution
# ==============================================================================

tuft_sample <- data.frame(
  sample_id = rownames(theta_final),
  tuft_initial = theta_initial[, "Tuft"],
  tuft_final = theta_final[, "Tuft"],
  tuft_delta = theta_final[, "Tuft"] - theta_initial[, "Tuft"],
  reliable_mass = reliable,
  stringsAsFactors = FALSE
)

if (!is.null(mass_qc)) {

  idx <- match(tuft_sample$sample_id, mass_qc$sample_id)

  tuft_sample$stage1_epithelial_mass <-
    mass_qc$original_mass[idx]

  tuft_sample$stage2_integerized_mass <-
    mass_qc$posterior_allocated_mass[idx]

  tuft_sample$fractional_mass_loss <-
    mass_qc$relative_mass_loss[idx]
}

write_tsv(tuft_sample, OUT_SAMPLE)


# ==============================================================================
# 6. Reference pairwise similarity
# ==============================================================================

spearman0 <- cor(
  t(phi0),
  method = "spearman"
)

write_tsv(
  data.frame(component = rownames(spearman0), spearman0, check.names = FALSE),
  OUT_COR_ORIG
)

cat("Original-reference Spearman similarity to Tuft:\n")
print(
  sort(
    spearman0["Tuft", setdiff(colnames(spearman0), "Tuft")],
    decreasing = TRUE
  )
)

cat("\n")


if (!is.null(phi1)) {

  spearman1 <- cor(
    t(phi1),
    method = "spearman"
  )

  write_tsv(
    data.frame(component = rownames(spearman1), spearman1, check.names = FALSE),
    OUT_COR_UPDATED
  )

  cat("Updated-reference Spearman similarity to Tuft:\n")
  print(
    sort(
      spearman1["Tuft", setdiff(colnames(spearman1), "Tuft")],
      decreasing = TRUE
    )
  )

  cat("\n")
}


# ==============================================================================
# 7. Original vs updated reference stability
# ==============================================================================

if (!is.null(phi1)) {

  update_stability <- do.call(
    rbind,
    lapply(rownames(phi0), function(component) {

      x <- phi0[component, ]
      y <- phi1[component, ]

      data.frame(
        component = component,
        pearson = cor(x, y, method = "pearson"),
        spearman = cor(x, y, method = "spearman"),
        cosine = sum(x * y) /
          sqrt(sum(x^2) * sum(y^2)),
        stringsAsFactors = FALSE
      )
    })
  )

  write_tsv(update_stability, OUT_SELF_UPDATE)

  cat("Original -> updated reference stability:\n")
  print(update_stability, row.names = FALSE)
  cat("\n")
}


# ==============================================================================
# 8. Load Tuft posterior Z
# ==============================================================================

assert_true(file.exists(TUFT_Z_RDS), "Tuft Z matrix missing.")

cat("Loading Tuft posterior Z...\n")

tuft_z <- readRDS(TUFT_Z_RDS)

assert_true(
  identical(colnames(tuft_z), colnames(phi0)),
  "Tuft Z gene order differs from reference."
)

cat("[PASS] Tuft Z:", nrow(tuft_z), "samples x", ncol(tuft_z), "genes\n\n")


# ==============================================================================
# 9. Tuft posterior Z driver genes
# ==============================================================================

tuft_mass_gene <- colSums(tuft_z)

total_tuft_mass <- sum(tuft_mass_gene)

tuft_mass_fraction <- tuft_mass_gene / total_tuft_mass

tuft_rank <- rank(
  -tuft_mass_gene,
  ties.method = "min"
)

other_components <- setdiff(rownames(phi0), "Tuft")

phi_tuft <- phi0["Tuft", ]
phi_other <- phi0[other_components, , drop = FALSE]

max_other <- apply(phi_other, 2, max)

max_other_component <- other_components[
  max.col(t(phi_other), ties.method = "first")
]

phi_sum <- colSums(phi0)

eps <- 1e-12

tuft_reference_share <- phi_tuft / phi_sum

tuft_vs_max_other <- (phi_tuft + eps) /
  (max_other + eps)

driver_table <- data.frame(
  gene = colnames(phi0),
  tuft_Z_mass = tuft_mass_gene,
  tuft_Z_fraction = tuft_mass_fraction,
  tuft_Z_rank = tuft_rank,
  reference_Tuft = phi_tuft,
  reference_max_other = max_other,
  max_other_component = max_other_component,
  tuft_reference_share = tuft_reference_share,
  tuft_vs_max_other_fold = tuft_vs_max_other,
  stringsAsFactors = FALSE
)

driver_table <- driver_table[
  order(-driver_table$tuft_Z_mass),
  ,
  drop = FALSE
]

driver_table$cumulative_tuft_Z_fraction <-
  cumsum(driver_table$tuft_Z_fraction)

write_tsv(driver_table, OUT_TUFT_DRIVERS)


# ==============================================================================
# 10. Data-driven Tuft reference specificity
# ==============================================================================

specificity_table <- data.frame(
  gene = colnames(phi0),
  reference_Tuft = phi_tuft,
  reference_max_other = max_other,
  max_other_component = max_other_component,
  tuft_reference_share = tuft_reference_share,
  tuft_vs_max_other_fold = tuft_vs_max_other,
  tuft_Z_mass = tuft_mass_gene,
  tuft_Z_fraction = tuft_mass_fraction,
  stringsAsFactors = FALSE
)

specificity_table <- specificity_table[
  order(
    -specificity_table$tuft_reference_share,
    -specificity_table$reference_Tuft
  ),
  ,
  drop = FALSE
]

write_tsv(specificity_table, OUT_TUFT_SPEC)


# ==============================================================================
# 11. Canonical Tuft marker audit
# ==============================================================================

tuft_markers <- c(
  "POU2F3",
  "TRPM5",
  "PLCB2",
  "GNAT3",
  "GFI1B",
  "AVIL",
  "IL25",
  "PTGS1",
  "DCLK1",
  "SH2D6"
)

marker_rows <- lapply(tuft_markers, function(gene) {

  if (!gene %in% colnames(phi0)) {

    return(
      data.frame(
        gene = gene,
        present = FALSE,
        reference_Tuft = NA_real_,
        reference_max_other = NA_real_,
        max_other_component = NA_character_,
        tuft_reference_share = NA_real_,
        tuft_vs_max_other_fold = NA_real_,
        tuft_Z_mass = NA_real_,
        tuft_Z_fraction = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  j <- match(gene, colnames(phi0))

  data.frame(
    gene = gene,
    present = TRUE,
    reference_Tuft = phi_tuft[j],
    reference_max_other = max_other[j],
    max_other_component = max_other_component[j],
    tuft_reference_share = tuft_reference_share[j],
    tuft_vs_max_other_fold = tuft_vs_max_other[j],
    tuft_Z_mass = tuft_mass_gene[j],
    tuft_Z_fraction = tuft_mass_fraction[j],
    stringsAsFactors = FALSE
  )
})

marker_table <- do.call(rbind, marker_rows)

write_tsv(marker_table, OUT_MARKERS)


# ==============================================================================
# 12. Diagnostic summaries
# ==============================================================================

top10_fraction <- sum(
  head(driver_table$tuft_Z_fraction, 10)
)

top50_fraction <- sum(
  head(driver_table$tuft_Z_fraction, 50)
)

top100_fraction <- sum(
  head(driver_table$tuft_Z_fraction, 100)
)

top500_fraction <- sum(
  head(driver_table$tuft_Z_fraction, 500)
)

tuft_final <- theta_final[, "Tuft"]

tuft_initial <- theta_initial[, "Tuft"]

summary_table <- data.frame(
  metric = c(
    "tuft_initial_mean",
    "tuft_initial_median",
    "tuft_final_mean",
    "tuft_final_median",
    "tuft_final_q05",
    "tuft_final_q95",
    "tuft_initial_vs_final_spearman",
    "tuft_fraction_gt_0.10",
    "tuft_fraction_gt_0.20",
    "tuft_fraction_gt_0.30",
    "tuft_Z_top10_gene_fraction",
    "tuft_Z_top50_gene_fraction",
    "tuft_Z_top100_gene_fraction",
    "tuft_Z_top500_gene_fraction",
    "tuft_driver_top100_median_reference_share",
    "tuft_driver_top100_median_vs_max_other_fold"
  ),
  value = c(
    mean(tuft_initial),
    median(tuft_initial),
    mean(tuft_final),
    median(tuft_final),
    unname(quantile(tuft_final, 0.05)),
    unname(quantile(tuft_final, 0.95)),
    cor(tuft_initial, tuft_final, method = "spearman"),
    mean(tuft_final > 0.10),
    mean(tuft_final > 0.20),
    mean(tuft_final > 0.30),
    top10_fraction,
    top50_fraction,
    top100_fraction,
    top500_fraction,
    median(
      head(driver_table$tuft_reference_share, 100),
      na.rm = TRUE
    ),
    median(
      head(driver_table$tuft_vs_max_other_fold, 100),
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

write_tsv(summary_table, OUT_SUMMARY)


# ==============================================================================
# 13. Console report
# ==============================================================================

cat("============================================================\n")
cat("TUFT FRACTION DISTRIBUTION\n")
cat("============================================================\n\n")

cat("Initial mean   :", mean(tuft_initial), "\n")
cat("Initial median :", median(tuft_initial), "\n")
cat("Final mean     :", mean(tuft_final), "\n")
cat("Final median   :", median(tuft_final), "\n")
cat("Final Q05      :", quantile(tuft_final, 0.05), "\n")
cat("Final Q95      :", quantile(tuft_final, 0.95), "\n")
cat("Fraction >10%  :", mean(tuft_final > 0.10), "\n")
cat("Fraction >20%  :", mean(tuft_final > 0.20), "\n")
cat("Fraction >30%  :", mean(tuft_final > 0.30), "\n\n")


cat("============================================================\n")
cat("TOP 30 TUFT POSTERIOR-Z DRIVER GENES\n")
cat("============================================================\n\n")

print(
  head(
    driver_table[
      ,
      c(
        "gene",
        "tuft_Z_fraction",
        "cumulative_tuft_Z_fraction",
        "tuft_reference_share",
        "tuft_vs_max_other_fold",
        "max_other_component"
      )
    ],
    30
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("CANONICAL TUFT MARKERS\n")
cat("============================================================\n\n")

print(
  marker_table,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("TOP-GENE CONCENTRATION\n")
cat("============================================================\n\n")

cat("Top 10 genes :", top10_fraction, "\n")
cat("Top 50 genes :", top50_fraction, "\n")
cat("Top 100 genes:", top100_fraction, "\n")
cat("Top 500 genes:", top500_fraction, "\n\n")


cat("============================================================\n")
cat("INTERPRETATION TARGETS\n")
cat("============================================================\n\n")

cat(
  "1. Is Tuft highly correlated with another epithelial reference?\n",
  "2. Are canonical Tuft markers actually Tuft-specific in phi?\n",
  "3. Are posterior Tuft reads driven by Tuft-specific genes,\n",
  "   or by abundant shared epithelial genes?\n",
  "4. Did reference updating materially alter the Tuft profile?\n",
  sep = ""
)

capture.output(sessionInfo(), file = OUT_SESSION)

cat("\nOutput directory:\n")
cat(QC_DIR, "\n")

cat("\nEND:", format(Sys.time()), "\n")