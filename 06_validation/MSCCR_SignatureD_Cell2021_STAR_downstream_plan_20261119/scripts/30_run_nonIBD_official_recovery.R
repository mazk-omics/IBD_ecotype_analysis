#!/usr/bin/env Rscript
# 本脚本按 IBD-only 相同参数调用未修改的 EcoTyper 官方脚本，恢复一套 nonIBD scRNA 队列的 cell states。
# 每次只运行一个 dataset；四套队列可在独立 tmux 窗口并行。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取已验收输入并冻结任务表
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 30_run_nonIBD_official_recovery.R RUN_ROOT DATASET", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
dataset <- args[[2]]
source(file.path(RUN_ROOT, "config/run_paths.R"))

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing R package: data.table", call. = FALSE)
library(data.table)

datasets <- c("SCP1884", "GSE214695", "SCP259", "IBDverse_E_MTAB_16999")
manifest_files <- c(
  SCP1884 = "SCP1884_recovery_input_manifest.tsv",
  GSE214695 = "GSE214695_recovery_input_manifest.tsv",
  SCP259 = "SCP259_recovery_input_manifest.tsv",
  IBDverse_E_MTAB_16999 = "IBDverse_recovery_input_manifest.tsv"
)
expected_parent_types <- c(SCP1884 = 20L, GSE214695 = 18L, SCP259 = 18L, IBDverse_E_MTAB_16999 = 16L)
if (!dataset %in% datasets) stop("Unknown nonIBD dataset: ", dataset, call. = FALSE)

input_audit <- readLines(file.path(output_root, "18_nonIBD_scrna_input_audit/nonIBD_recovery_input_audit.json"), warn = FALSE)
if (!any(grepl("PASS_NONIBD_RECOVERY_INPUT_AUDIT", input_audit, fixed = TRUE))) stop("nonIBD recovery input audit did not pass", call. = FALSE)

prepare_dir <- file.path(RUN_ROOT, "prepared_scRNA_nonIBD", dataset)
manifest <- fread(file.path(prepare_dir, manifest_files[[dataset]]), colClasses = "character")
if (nrow(manifest) != expected_parent_types[[dataset]]) stop("nonIBD parent type count changed: ", dataset, call. = FALSE)
manifest[, DatasetID := paste(dataset, "nonIBD_Cell2021_20261123", CellType, sep = "_")]
manifest[, StageDir := file.path(ecotyper_root, "datasets/scRNA", DatasetID)]
manifest[, ResultDir := file.path(model_root, "Cell_States/recovery/scRNA", DatasetID, CellType, FrozenRank)]

log_dir <- file.path(output_root, "19_nonIBD_official_recovery", dataset)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
manifest[, LogFile := file.path(log_dir, paste0(CellType, ".console.txt"))]
manifest[, ExitFile := file.path(log_dir, paste0(CellType, ".exit.txt"))]
if (any(dir.exists(manifest$StageDir))) stop("One or more nonIBD staging directories already exist: ", dataset, call. = FALSE)
if (any(dir.exists(manifest$ResultDir))) stop("One or more nonIBD recovery result directories already exist: ", dataset, call. = FALSE)

official_script <- file.path(ecotyper_root, "pipeline/state_recovery_scRNA.R")
rscript_bin <- "/home/mazekai/miniconda3/envs/ecotyper/bin/Rscript"
official_sha256 <- strsplit(system2("sha256sum", official_script, stdout = TRUE), "[[:space:]]+")[[1L]][1L]
if (official_sha256 != "ccd7380a6962ac6905bff1c4890440991b130d574efd5bb69c52d7ebb4b2332d") stop("Official state recovery script hash changed", call. = FALSE)
fwrite(manifest, file.path(log_dir, "official_recovery_jobs.tsv"), sep = "\t")

# Step2. 建立官方输入链接
for (row_index in seq_len(nrow(manifest))) {
  dir.create(manifest$StageDir[[row_index]], recursive = TRUE, showWarnings = FALSE)
  expression_link <- file.path(manifest$StageDir[[row_index]], "data.txt")
  annotation_link <- file.path(manifest$StageDir[[row_index]], "annotation.txt")
  link_ok <- c(file.symlink(manifest$ExpressionFile[[row_index]], expression_link), file.symlink(manifest$AnnotationFile[[row_index]], annotation_link))
  if (!all(link_ok)) stop("Failed to create nonIBD official input links: ", dataset, " ", manifest$CellType[[row_index]], call. = FALSE)
}

# Step3. 用四个 worker 调用官方 recovery
workers <- 4L
pipeline_dir <- file.path(ecotyper_root, "pipeline")
job_status <- parallel::mclapply(
  seq_len(nrow(manifest)),
  function(row_index) {
    cell_type <- manifest$CellType[[row_index]]
    command_args <- c(
      official_script, dataset_name, fractions_name, cell_type, manifest$FrozenRank[[row_index]], manifest$DatasetID[[row_index]],
      "FALSE", "-1", "Disease", "Inflammation", "Region3"
    )
    setwd(pipeline_dir)
    cat("NONIBD_OFFICIAL_RECOVERY_START", dataset, cell_type, manifest$DatasetID[[row_index]], "\n")
    exit_status <- system2(rscript_bin, args = command_args, stdout = manifest$LogFile[[row_index]], stderr = manifest$LogFile[[row_index]])
    writeLines(paste0("EXIT_STATUS=", exit_status), manifest$ExitFile[[row_index]])
    cat("NONIBD_OFFICIAL_RECOVERY_END", dataset, cell_type, exit_status, "\n")
    exit_status
  },
  mc.cores = workers,
  mc.preschedule = FALSE
)

# Step4. 保存官方 recovery 退出状态
manifest[, ExitStatus := as.integer(unlist(job_status))]
fwrite(manifest, file.path(log_dir, "official_recovery_jobs_completed.tsv"), sep = "\t")
if (any(manifest$ExitStatus != 0L)) stop("One or more official nonIBD recovery jobs failed: ", dataset, call. = FALSE)
writeLines(
  c(
    "STATUS=PASS_NONIBD_OFFICIAL_RECOVERY_EXECUTION", paste0("DATASET=", dataset), paste0("PARENT_TYPES=", nrow(manifest)),
    paste0("WORKERS=", workers), paste0("OFFICIAL_SCRIPT_SHA256=", official_sha256), "BUILTIN_SIGNIFICANCE=FALSE", "RECOVERY_SIGNIFICANCE_BOOTSTRAPS_PENDING=100"
  ),
  file.path(log_dir, "official_recovery_execution_summary.txt")
)
cat("NONIBD_OFFICIAL_RECOVERY_EXECUTION_PASS", dataset, nrow(manifest), "parent types\n")
