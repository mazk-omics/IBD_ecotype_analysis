#!/usr/bin/env Rscript
# 本脚本从已完成的 100 次 Cell 2021 shuffle 中截取前 30 次，派生官方仓库口径的敏感性结果。
# 不重跑 NMF，不覆盖 100 次主结果；在服务器 RUN_ROOT 下运行。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取已经验收的 100 次结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04f_derive_GSE214695_repo30_sensitivity.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

input_dir <- file.path(RUN_ROOT, "output/04e_GSE214695_recovery_significance_100")
output_dir <- file.path(RUN_ROOT, "output/04f_GSE214695_recovery_significance_repo30")
jobs <- read.delim(file.path(input_dir, "GSE214695_recovery_significance_100_jobs.tsv"), check.names = FALSE)
summary_100 <- readLines(file.path(input_dir, "GSE214695_recovery_significance_100_summary.txt"), warn = FALSE)
if (!any(summary_100 == "STATUS=PASS_GSE214695_RECOVERY_SIGNIFICANCE_100")) stop("The 100-shuffle run did not pass", call. = FALSE)
if (nrow(jobs) != 18L || any(jobs$SignificanceExitStatus != 0L)) stop("The 100-shuffle job table is incomplete", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step2. 逐细胞类型派生前 30 次的 Z 值
comparison_all <- list()
audit_all <- list()

for (job_index in seq_len(nrow(jobs))) {
  cell_type <- jobs$CellType[[job_index]]
  result_dir <- jobs$ResultDir[[job_index]]
  cell_output_dir <- file.path(output_dir, cell_type)
  dir.create(cell_output_dir, recursive = TRUE, showWarnings = FALSE)

  bootstrap_file <- file.path(result_dir, "bootstrap_marker_fc_100.txt")
  z100_initial_file <- file.path(result_dir, "recovery_z_scores_initial_states.txt")
  z100_file <- file.path(result_dir, "recovery_z_scores.txt")
  mapping_file <- file.path(state_discovery_root, cell_type, jobs$FrozenRank[[job_index]], "mapping_to_initial_states.txt")
  required_files <- c(bootstrap_file, z100_initial_file, z100_file, mapping_file)
  if (!all(file.exists(required_files))) stop("Missing 100-shuffle input for ", cell_type, call. = FALSE)

  bootstrap_100 <- read.delim(bootstrap_file, row.names = 1L, check.names = FALSE)
  z100_initial <- read.delim(z100_initial_file, check.names = FALSE)
  z100 <- read.delim(z100_file, check.names = FALSE)
  mapping <- read.delim(mapping_file, check.names = FALSE)
  if (nrow(bootstrap_100) != 100L) stop("Bootstrap row count changed for ", cell_type, call. = FALSE)
  if (!identical(rownames(bootstrap_100), sprintf("B%03d", seq_len(100L)))) stop("Bootstrap order changed for ", cell_type, call. = FALSE)
  if (!setequal(colnames(bootstrap_100), z100_initial$State)) stop("Initial states do not match for ", cell_type, call. = FALSE)

  bootstrap_30 <- bootstrap_100[seq_len(30L), , drop = FALSE]
  mean_30 <- colMeans(bootstrap_30, na.rm = TRUE)
  sd_30 <- apply(bootstrap_30, 2L, sd, na.rm = TRUE)
  mean_100 <- colMeans(bootstrap_100, na.rm = TRUE)
  sd_100 <- apply(bootstrap_100, 2L, sd, na.rm = TRUE)
  mean_30 <- mean_30[match(z100_initial$State, names(mean_30))]
  sd_30 <- sd_30[match(z100_initial$State, names(sd_30))]
  mean_100 <- mean_100[match(z100_initial$State, names(mean_100))]
  sd_100 <- sd_100[match(z100_initial$State, names(sd_100))]
  z100_rebuilt <- (z100_initial$RefMean - mean_100) / sd_100
  if (!isTRUE(all.equal(z100_initial$Z, z100_rebuilt, tolerance = 1e-12, check.attributes = FALSE))) stop("Cannot reproduce saved 100-shuffle Z for ", cell_type, call. = FALSE)

  z30_initial <- data.frame(
    State = z100_initial$State,
    Z = (z100_initial$RefMean - mean_30) / sd_30,
    RefMean = z100_initial$RefMean,
    BootstrapMean = as.numeric(mean_30),
    BootstrapSD = as.numeric(sd_30)
  )
  z30 <- z30_initial[match(mapping[, 2L], z30_initial$State), , drop = FALSE]
  z30$State <- as.character(mapping[, 1L])
  if (!identical(z30$State, z100$State)) stop("Mapped states changed for ", cell_type, call. = FALSE)

  comparison <- data.frame(
    Dataset = "GSE214695",
    CellType = cell_type,
    State = z100$State,
    Z_repo30 = z30$Z,
    Z_paper100 = z100$Z,
    DeltaZ_repo30_minus_paper100 = z30$Z - z100$Z,
    Positive_repo30 = is.finite(z30$Z) & z30$Z > 1.96,
    Positive_paper100 = is.finite(z100$Z) & z100$Z > 1.96,
    Same_threshold_call = (is.finite(z30$Z) & z30$Z > 1.96) == (is.finite(z100$Z) & z100$Z > 1.96)
  )
  comparison_all[[job_index]] <- comparison
  audit_all[[job_index]] <- data.frame(
    Dataset = "GSE214695",
    CellType = cell_type,
    States = nrow(comparison),
    FiniteZ_repo30 = sum(is.finite(comparison$Z_repo30)),
    FiniteZ_paper100 = sum(is.finite(comparison$Z_paper100)),
    Positive_repo30 = sum(comparison$Positive_repo30, na.rm = TRUE),
    Positive_paper100 = sum(comparison$Positive_paper100, na.rm = TRUE),
    Same_threshold_call = sum(comparison$Same_threshold_call, na.rm = TRUE),
    MaxAbsoluteDeltaZ = max(abs(comparison$DeltaZ_repo30_minus_paper100), na.rm = TRUE)
  )

  write.table(bootstrap_30, file.path(cell_output_dir, "bootstrap_marker_fc_first30.txt"), sep = "\t")
  write.table(z30_initial, file.path(cell_output_dir, "recovery_z_scores_initial_states.txt"), sep = "\t", row.names = FALSE)
  write.table(z30, file.path(cell_output_dir, "recovery_z_scores.txt"), sep = "\t", row.names = FALSE)
  write.table(comparison, file.path(cell_output_dir, "repo30_vs_paper100.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  cat("DERIVED_REPO30", cell_type, nrow(comparison), "states\n")
}

# Step3. 汇总 30 次与 100 次口径差异
comparison_all <- do.call(rbind, comparison_all)
audit_all <- do.call(rbind, audit_all)
write.table(comparison_all, file.path(output_dir, "GSE214695_repo30_vs_paper100_by_state.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(audit_all, file.path(output_dir, "GSE214695_repo30_vs_paper100_by_cell_type.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(
  c(
    "STATUS=PASS_GSE214695_REPO30_SENSITIVITY",
    "PRIMARY_METHOD=PAPER_100_SHUFFLES",
    "SENSITIVITY_METHOD=OFFICIAL_REPOSITORY_30_SHUFFLES",
    "DERIVATION=FIRST_30_OF_FIXED_SEED_100_SHUFFLES",
    "NMF_RERUN=FALSE",
    "PAPER100_RESULTS_OVERWRITTEN=FALSE",
    paste0("PARENT_TYPES=", nrow(audit_all)),
    paste0("STATES=", nrow(comparison_all)),
    paste0("THRESHOLD_CALLS_IDENTICAL=", sum(comparison_all$Same_threshold_call, na.rm = TRUE)),
    paste0("THRESHOLD_CALLS_TOTAL=", nrow(comparison_all))
  ),
  file.path(output_dir, "GSE214695_repo30_sensitivity_summary.txt")
)
cat("GSE214695_REPO30_SENSITIVITY_PASS", nrow(comparison_all), "states\n")
