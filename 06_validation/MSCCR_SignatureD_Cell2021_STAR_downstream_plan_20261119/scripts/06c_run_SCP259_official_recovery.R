#!/usr/bin/env Rscript
# 本脚本并行调用未修改的 EcoTyper 官方脚本，恢复 SCP259 单细胞状态。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取已审计输入并冻结任务表
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 06c_run_SCP259_official_recovery.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing R package: data.table", call. = FALSE)
library(data.table)

prepare_dir <- file.path(RUN_ROOT, "prepared_scRNA/SCP259")
audit_status <- readLines(file.path(prepare_dir, "SCP259_recovery_input_audit.json"), warn = FALSE)
if (!any(grepl("PASS_SCP259_RECOVERY_INPUT_AUDIT", audit_status, fixed = TRUE))) stop("SCP259 recovery input audit did not pass", call. = FALSE)

manifest <- fread(file.path(prepare_dir, "SCP259_recovery_input_manifest.tsv"), colClasses = "character")
manifest[, DatasetID := paste("SCP259_Cell2021_20261121", CellType, sep = "_")]
manifest[, StageDir := file.path(ecotyper_root, "datasets/scRNA", DatasetID)]
manifest[, ResultDir := file.path(model_root, "Cell_States/recovery/scRNA", DatasetID, CellType, FrozenRank)]
manifest[, LogFile := file.path(RUN_ROOT, "output/06c_SCP259_official_recovery", paste0(CellType, ".console.txt"))]
manifest[, ExitFile := file.path(RUN_ROOT, "output/06c_SCP259_official_recovery", paste0(CellType, ".exit.txt"))]
if (any(dir.exists(manifest$StageDir))) stop("One or more SCP259 staging directories already exist", call. = FALSE)
if (any(dir.exists(manifest$ResultDir))) stop("One or more SCP259 recovery result directories already exist", call. = FALSE)

official_script <- file.path(ecotyper_root, "pipeline/state_recovery_scRNA.R")
rscript_bin <- "/home/mazekai/miniconda3/envs/ecotyper/bin/Rscript"
official_sha256 <- strsplit(system2("sha256sum", official_script, stdout = TRUE), "[[:space:]]+")[[1L]][1L]
if (official_sha256 != "ccd7380a6962ac6905bff1c4890440991b130d574efd5bb69c52d7ebb4b2332d") stop("Official state recovery script hash changed", call. = FALSE)

log_dir <- file.path(RUN_ROOT, "output/06c_SCP259_official_recovery")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(manifest, file.path(log_dir, "SCP259_official_recovery_jobs.tsv"), sep = "\t")

# Step2. 建立每个父类型的官方输入目录
for (row_index in seq_len(nrow(manifest))) {
  dir.create(manifest$StageDir[[row_index]], recursive = TRUE, showWarnings = FALSE)
  expression_link <- file.path(manifest$StageDir[[row_index]], "data.txt")
  annotation_link <- file.path(manifest$StageDir[[row_index]], "annotation.txt")
  link_ok <- c(
    file.symlink(manifest$ExpressionFile[[row_index]], expression_link),
    file.symlink(manifest$AnnotationFile[[row_index]], annotation_link)
  )
  if (!all(link_ok)) stop("Failed to create official input links: ", manifest$CellType[[row_index]], call. = FALSE)
}

# Step3. 用四个 worker 调用官方恢复脚本
workers <- 4L
pipeline_dir <- file.path(ecotyper_root, "pipeline")
job_status <- parallel::mclapply(
  seq_len(nrow(manifest)),
  function(row_index) {
    cell_type <- manifest$CellType[[row_index]]
    command_args <- c(
      official_script,
      dataset_name,
      fractions_name,
      cell_type,
      manifest$FrozenRank[[row_index]],
      manifest$DatasetID[[row_index]],
      "FALSE",
      "-1",
      "Disease",
      "Inflammation",
      "Region3"
    )
    setwd(pipeline_dir)
    cat("SCP259_OFFICIAL_RECOVERY_START", cell_type, manifest$DatasetID[[row_index]], "\n")
    exit_status <- system2(
      rscript_bin,
      args = command_args,
      stdout = manifest$LogFile[[row_index]],
      stderr = manifest$LogFile[[row_index]]
    )
    writeLines(paste0("EXIT_STATUS=", exit_status), manifest$ExitFile[[row_index]])
    cat("SCP259_OFFICIAL_RECOVERY_END", cell_type, exit_status, "\n")
    exit_status
  },
  mc.cores = workers,
  mc.preschedule = FALSE
)

# Step4. 保存官方恢复退出状态
manifest[, ExitStatus := as.integer(unlist(job_status))]
fwrite(manifest, file.path(log_dir, "SCP259_official_recovery_jobs_completed.tsv"), sep = "\t")
if (any(manifest$ExitStatus != 0L)) stop("One or more official SCP259 recovery jobs failed", call. = FALSE)
writeLines(
  c(
    "STATUS=PASS_SCP259_OFFICIAL_RECOVERY_EXECUTION",
    paste0("PARENT_TYPES=", nrow(manifest)),
    paste0("WORKERS=", workers),
    paste0("OFFICIAL_SCRIPT_SHA256=", official_sha256),
    "BUILTIN_SIGNIFICANCE=FALSE",
    "RECOVERY_SIGNIFICANCE_BOOTSTRAPS_PENDING=100"
  ),
  file.path(log_dir, "SCP259_official_recovery_execution_summary.txt")
)
cat("SCP259_OFFICIAL_RECOVERY_EXECUTION_PASS", nrow(manifest), "parent types\n")
