#!/usr/bin/env Rscript
# 本脚本按 Cell 2021 方法执行 100 次 recovery shuffle，并省略无统计用途的大型中间矩阵。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取通过审计的官方 recovery 任务
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04e_run_GSE214695_recovery_significance_100.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

suppressPackageStartupMessages({
  library(NMF)
  library(data.table)
})
source(file.path(ecotyper_root, "pipeline/lib/misc.R"))
source(file.path(ecotyper_root, "pipeline/lib/ecotyper.R"))

recovery_audit <- readLines(file.path(RUN_ROOT, "output/04d_GSE214695_official_recovery_audit/GSE214695_official_recovery_audit.json"), warn = FALSE)
if (!any(grepl("PASS_GSE214695_OFFICIAL_RECOVERY_AUDIT", recovery_audit, fixed = TRUE))) stop("GSE214695 official recovery audit did not pass", call. = FALSE)
jobs <- fread(file.path(RUN_ROOT, "output/04c_GSE214695_official_recovery/GSE214695_official_recovery_jobs_completed.tsv"), colClasses = "character")
if (nrow(jobs) != 18L) stop("GSE214695 recovery job count changed", call. = FALSE)
if (any(file.exists(file.path(jobs$ResultDir, "coef_null.txt")))) stop("One or more 100-shuffle output files already exist", call. = FALSE)

ecotyper_lib_sha256 <- strsplit(system2("sha256sum", file.path(ecotyper_root, "pipeline/lib/ecotyper.R"), stdout = TRUE), "[[:space:]]+")[[1L]][1L]
if (ecotyper_lib_sha256 != "be8dc381634d025ccec3e9f93a1f1a2d91c1b24c3bb4f35a613c28277e8bf4b9") stop("Official EcoTyper significance code hash changed", call. = FALSE)

log_dir <- file.path(RUN_ROOT, "output/04e_GSE214695_recovery_significance_100")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Step2. 定义官方数据读取和等价的流式统计
read_recovery_bundle <- function(row_index) {
  cell_type <- jobs$CellType[[row_index]]
  discovery_dir <- file.path(state_discovery_root, cell_type, jobs$FrozenRank[[row_index]])
  W_local <- read.delim(file.path(discovery_dir, "W.txt"), check.names = FALSE)
  mapping_local <- read.delim(file.path(discovery_dir, "mapping_to_initial_states.txt"), check.names = FALSE)
  gene_info_local <- read.delim(file.path(discovery_dir, "initial_gene_info.txt"), check.names = FALSE)
  assignment_local <- read.delim(file.path(jobs$ResultDir[[row_index]], "initial_state_assignment.txt"), check.names = FALSE)
  data_local <- read.delim(file.path(jobs$ResultDir[[row_index]], "expression_matrix_scaled.txt"), row.names = 1L, check.names = FALSE)
  data_local <- as.matrix(data_local)
  storage.mode(data_local) <- "double"
  data_local <- data_local[match(gene_info_local$Gene, rownames(data_local)), match(assignment_local$ID, colnames(data_local)), drop = FALSE]
  rownames(data_local) <- gene_info_local$Gene
  data_local[is.na(data_local)] <- 0
  list(data = data_local, W = W_local, mapping = mapping_local, gene_info = gene_info_local, assignment = assignment_local)
}

run_streamed_significance <- function(bundle, output_dir, n_bootstraps) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data_local <- bundle$data
  W_local <- bundle$W
  mapping_local <- bundle$mapping
  gene_info_local <- bundle$gene_info
  assignment_local <- bundle$assignment
  n_states_local <- ncol(W_local)
  state_levels <- levels(as.factor(gene_info_local$State))
  ref_average <- averge_FC_marker_genes(data_local, gene_info_local, assignment_local)
  coef_null <- matrix(NA_real_, nrow = n_states_local * n_bootstraps, ncol = ncol(data_local))
  colnames(coef_null) <- colnames(data_local)
  boot_averages <- matrix(NA_real_, nrow = n_bootstraps, ncol = length(state_levels), dimnames = list(sprintf("B%03d", seq_len(n_bootstraps)), state_levels))

  for (bootstrap_index in seq_len(n_bootstraps)) {
    set.seed(100 + bootstrap_index)
    rand_data <- t(apply(data_local, 1L, sample))
    colnames(rand_data) <- colnames(data_local)
    nmf_null <- NMFpredict(W_local, rand_data)
    coef <- apply(nmf_null@fit@H, 2L, function(x) x / sum(x))
    rownames(coef) <- sprintf("IS%02d.B%02d", seq_len(n_states_local), bootstrap_index)
    coef_rows <- ((bootstrap_index - 1L) * n_states_local + 1L):(bootstrap_index * n_states_local)
    coef_null[coef_rows, ] <- coef
    null_assignment <- data.frame(ID = colnames(rand_data), State = apply(coef, 2L, function(x) colnames(W_local)[which.max(x)]))
    boot_average <- averge_FC_marker_genes(rand_data, gene_info_local, null_assignment)
    boot_averages[bootstrap_index, names(boot_average)] <- boot_average
    if (bootstrap_index %% 10L == 0L || bootstrap_index == n_bootstraps) cat("BOOTSTRAPS_COMPLETED", bootstrap_index, "of", n_bootstraps, "\n")
  }

  rownames(coef_null) <- unlist(lapply(seq_len(n_bootstraps), function(x) sprintf("IS%02d.B%02d", seq_len(n_states_local), x)))
  write.table(coef_null, file.path(output_dir, "coef_null.txt"), sep = "\t")
  write.table(boot_averages, file.path(output_dir, "bootstrap_marker_fc_100.txt"), sep = "\t")
  means <- apply(boot_averages, 2L, function(x) mean(x, na.rm = TRUE))
  sds <- apply(boot_averages, 2L, function(x) sd(x, na.rm = TRUE))
  Z <- (ref_average - means) / sds
  zs_corr <- data.frame(State = names(ref_average), Z = Z, RefMean = as.numeric(ref_average), BootstrapMean = means, BootstrapSD = sds)
  write.table(zs_corr, file.path(output_dir, "recovery_z_scores_initial_states.txt"), sep = "\t", row.names = FALSE)
  filtered_z <- zs_corr[match(mapping_local[, 2L], zs_corr$State), ]
  filtered_z$State <- as.character(mapping_local[, 1L])
  write.table(filtered_z, file.path(output_dir, "recovery_z_scores.txt"), sep = "\t", row.names = FALSE)
  invisible(TRUE)
}

# Step3. 用 EEC 三次 shuffle 验证与官方函数逐值一致
test_index <- which(jobs$CellType == "Enteroendocrine")
test_bundle <- read_recovery_bundle(test_index)
test_official_dir <- file.path(log_dir, "equivalence_test_official_3")
test_streamed_dir <- file.path(log_dir, "equivalence_test_streamed_3")
dir.create(test_official_dir, recursive = TRUE, showWarnings = FALSE)
W <- test_bundle$W
mapping <- test_bundle$mapping
n_states <- ncol(W)
calculate_recovery_significance(test_bundle$data, test_bundle$gene_info, test_bundle$assignment, test_official_dir, n_bootstraps = 3L)
run_streamed_significance(test_bundle, test_streamed_dir, n_bootstraps = 3L)

official_coef <- read.delim(file.path(test_official_dir, "coef_null.txt"), row.names = 1L, check.names = FALSE)
streamed_coef <- read.delim(file.path(test_streamed_dir, "coef_null.txt"), row.names = 1L, check.names = FALSE)
official_z <- read.delim(file.path(test_official_dir, "recovery_z_scores_initial_states.txt"), check.names = FALSE)
streamed_z <- read.delim(file.path(test_streamed_dir, "recovery_z_scores_initial_states.txt"), check.names = FALSE)
coef_equal <- isTRUE(all.equal(official_coef, streamed_coef, tolerance = 0, check.attributes = FALSE))
z_equal <- isTRUE(all.equal(official_z, streamed_z, tolerance = 0, check.attributes = FALSE))
if (!coef_equal || !z_equal) stop("Streamed 100-shuffle implementation differs from official function", call. = FALSE)
writeLines(c("STATUS=PASS_STREAMED_SIGNIFICANCE_EQUIVALENCE", "CELL_TYPE=Enteroendocrine", "BOOTSTRAPS=3", paste0("COEF_EQUAL=", coef_equal), paste0("Z_EQUAL=", z_equal)), file.path(log_dir, "streamed_significance_equivalence.txt"))
rm(test_bundle, W, mapping, official_coef, streamed_coef, official_z, streamed_z)
invisible(gc(FALSE))

# Step4. 四个 worker 分别完成 100 次正式 shuffle
workers <- 4L
job_status <- parallel::mclapply(
  seq_len(nrow(jobs)),
  function(row_index) {
    cell_type <- jobs$CellType[[row_index]]
    job_log <- file.path(log_dir, paste0(cell_type, ".console.txt"))
    job_exit <- file.path(log_dir, paste0(cell_type, ".exit.txt"))
    log_connection <- file(job_log, open = "wt")
    sink(log_connection)
    sink(log_connection, type = "message")
    exit_status <- tryCatch({
      cat("SIGNIFICANCE_START", cell_type, "100 bootstraps\n")
      bundle <- read_recovery_bundle(row_index)
      run_streamed_significance(bundle, jobs$ResultDir[[row_index]], n_bootstraps = 100L)
      cat("SIGNIFICANCE_END", cell_type, "0\n")
      0L
    }, error = function(e) {
      cat("SIGNIFICANCE_ERROR", cell_type, conditionMessage(e), "\n")
      1L
    })
    sink(type = "message")
    sink()
    close(log_connection)
    writeLines(paste0("EXIT_STATUS=", exit_status), job_exit)
    exit_status
  },
  mc.cores = workers,
  mc.preschedule = FALSE
)

# Step5. 保存 100 次 shuffle 执行结论
jobs[, SignificanceExitStatus := as.integer(unlist(job_status))]
fwrite(jobs, file.path(log_dir, "GSE214695_recovery_significance_100_jobs.tsv"), sep = "\t")
if (any(jobs$SignificanceExitStatus != 0L)) stop("One or more 100-shuffle significance jobs failed", call. = FALSE)
writeLines(
  c(
    "STATUS=PASS_GSE214695_RECOVERY_SIGNIFICANCE_100",
    "BOOTSTRAPS=100",
    paste0("PARENT_TYPES=", nrow(jobs)),
    paste0("WORKERS=", workers),
    paste0("OFFICIAL_ECOTYPER_LIB_SHA256=", ecotyper_lib_sha256),
    "FULL_RANDOM_MATRICES_SAVED=FALSE",
    "FULL_NULL_NMF_OBJECTS_SAVED=FALSE",
    "NULL_COEFFICIENTS_SAVED=TRUE",
    "NULL_MARKER_FC_SAVED=TRUE"
  ),
  file.path(log_dir, "GSE214695_recovery_significance_100_summary.txt")
)
cat("GSE214695_RECOVERY_SIGNIFICANCE_100_PASS", nrow(jobs), "parent types\n")
