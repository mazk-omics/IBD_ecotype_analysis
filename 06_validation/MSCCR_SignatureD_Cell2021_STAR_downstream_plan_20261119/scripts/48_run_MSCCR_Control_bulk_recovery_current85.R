#!/usr/bin/env Rscript
# 本脚本按部位顺序调用 EcoTyper 官方 bulk recovery，仅恢复 MSCCR Control 样本。
# 输入为从原始 counts 按完整文库分母重建的线性 CPM；不重跑 discovery CD/UC 样本。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定任务路径和参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 48_run_MSCCR_Control_bulk_recovery_current85.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

project_root <- normalizePath(file.path(RUN_ROOT, "../.."), mustWork = TRUE)
ecotyper_root <- file.path(project_root, "software/ecotyper")
discovery_name <- "MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101"
model_root <- file.path(ecotyper_root, "EcoTyper", discovery_name, "Custom")
input_root <- file.path(
  RUN_ROOT,
  "output/35_MSCCR_Control_CPM_current85_preflight/reconstructed_current_model_gene_union/input"
)
preflight_root <- file.path(
  RUN_ROOT,
  "output/35_MSCCR_Control_CPM_current85_preflight/reconstructed_current_model_gene_union"
)
job_root <- file.path(RUN_ROOT, "output/36_MSCCR_Control_bulk_recovery_current85_retry1")
dir.create(job_root, recursive = TRUE, showWarnings = FALSE)

regions <- c("Ileum", "Colon", "Rectum")
threads <- 20L
conda_bin <- "/home/mazekai/miniconda3/bin/conda"
official_script <- file.path(ecotyper_root, "EcoTyper_recovery_bulk.R")
pipeline_files <- c(
  official_script,
  file.path(ecotyper_root, "pipeline/state_recovery_bulk.R"),
  file.path(ecotyper_root, "pipeline/ecotypes_recovery.R"),
  file.path(model_root, "Analysis/rank_selection/rank_data.txt")
)
if (any(!file.exists(c(conda_bin, pipeline_files)))) stop("EcoTyper executable, scripts or current model are missing", call. = FALSE)

session_info_text <- capture.output(sessionInfo())
r_subprocess_environment <- c("R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_ENVIRON", "R_PROFILE")
r_subprocess_environment_before <- data.table(
  Variable = r_subprocess_environment,
  ValueBefore = Sys.getenv(r_subprocess_environment, unset = "<unset>")
)
Sys.unsetenv(r_subprocess_environment)
r_subprocess_environment_before[, ValuePassedToConda := Sys.getenv(Variable, unset = "<unset>")]
fwrite(r_subprocess_environment_before, file.path(job_root, "R_subprocess_environment_audit.tsv"), sep = "\t", quote = FALSE)

# Step2. 冻结三个部位的输入和预检查结论
jobs <- data.table(
  Region3 = regions,
  ExpressionFile = file.path(input_root, paste0("MSCCR_Control_CPM_current85_", regions, ".txt")),
  AnnotationFile = file.path(input_root, paste0("MSCCR_Control_CPM_current85_", regions, "_annotation.tsv")),
  PreflightFile = file.path(preflight_root, paste0("preflight_", regions), "preflight_status.txt")
)
if (any(!file.exists(unlist(jobs[, .(ExpressionFile, AnnotationFile, PreflightFile)])))) stop("One or more MSCCR Control inputs or preflight files are missing", call. = FALSE)

jobs[, Preflight := vapply(PreflightFile, function(path) paste(readLines(path, warn = FALSE), collapse = "; "), character(1L))]
if (any(!grepl("status=ACCEPTED", jobs$Preflight, fixed = TRUE))) stop("One or more MSCCR Control preflight checks did not pass", call. = FALSE)

jobs[, RecoveryName := tools::file_path_sans_ext(basename(ExpressionFile))]
jobs[, RegionJobDir := file.path(job_root, Region3)]
jobs[, RecoveryOutputDir := file.path(RegionJobDir, "RecoveryOutput")]
jobs[, ResultRoot := file.path(RecoveryOutputDir, RecoveryName)]
if (any(dir.exists(jobs$ResultRoot))) stop("One or more MSCCR Control recovery result directories already exist", call. = FALSE)

jobs[, ExpressionMD5 := unname(tools::md5sum(ExpressionFile))]
jobs[, AnnotationMD5 := unname(tools::md5sum(AnnotationFile))]
jobs[, Threads := threads]
jobs[, Discovery := discovery_name]
fwrite(jobs, file.path(job_root, "MSCCR_Control_bulk_recovery_jobs.tsv"), sep = "\t", quote = FALSE)

pipeline_identity <- data.table(File = pipeline_files, MD5 = unname(tools::md5sum(pipeline_files)))
fwrite(pipeline_identity, file.path(job_root, "EcoTyper_pipeline_identity.tsv"), sep = "\t", quote = FALSE)

# Step3. 按 Ileum、Colon、Rectum 顺序运行官方 bulk recovery
for (job_index in seq_len(nrow(jobs))) {
  region <- jobs$Region3[[job_index]]
  region_job_dir <- jobs$RegionJobDir[[job_index]]
  dir.create(region_job_dir, recursive = TRUE, showWarnings = FALSE)
  status_file <- file.path(region_job_dir, "status.txt")
  console_file <- file.path(region_job_dir, "recovery.console.txt")
  exit_file <- file.path(region_job_dir, "recovery.exit.txt")

  writeLines(
    c("status=RUNNING", paste0("start=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), paste0("region=", region)),
    status_file
  )
  cat("MSCCR_CONTROL_BULK_RECOVERY_START", region, "\n")

  command_args <- c(
    "run", "--no-capture-output", "-n", "ecotyper", "Rscript", official_script,
    "-d", discovery_name,
    "-m", jobs$ExpressionFile[[job_index]],
    "-a", jobs$AnnotationFile[[job_index]],
    "-c", "Disease,DiseaseSubtype,Tissue,Histopathology",
    "-t", as.character(threads),
    "-o", jobs$RecoveryOutputDir[[job_index]]
  )
  old_dir <- getwd()
  setwd(ecotyper_root)
  exit_status <- system2(conda_bin, args = command_args, stdout = console_file, stderr = console_file)
  setwd(old_dir)
  writeLines(as.character(exit_status), exit_file)

  result_root <- jobs$ResultRoot[[job_index]]
  state_abundance_files <- list.files(result_root, pattern = "^state_abundances[.]txt$", recursive = TRUE, full.names = TRUE)
  state_assignment_files <- list.files(result_root, pattern = "^state_assignment[.]txt$", recursive = TRUE, full.names = TRUE)
  ecotype_abundance_file <- file.path(result_root, "Ecotypes", "ecotype_abundance.txt")
  ecotype_assignment_file <- file.path(result_root, "Ecotypes", "ecotype_assignment.txt")

  audit <- data.table(
    Region3 = region,
    ExitStatus = exit_status,
    StateAbundanceFiles = length(state_abundance_files),
    StateAssignmentFiles = length(state_assignment_files),
    EcotypeAbundancePresent = file.exists(ecotype_abundance_file) && file.info(ecotype_abundance_file)$size > 0,
    EcotypeAssignmentPresent = file.exists(ecotype_assignment_file) && file.info(ecotype_assignment_file)$size > 0
  )
  fwrite(audit, file.path(region_job_dir, "output_audit.tsv"), sep = "\t", quote = FALSE)

  accepted <- exit_status == 0L && audit$StateAbundanceFiles == 20L && audit$StateAssignmentFiles == 20L &&
    audit$EcotypeAbundancePresent && audit$EcotypeAssignmentPresent
  if (!accepted) {
    writeLines(
      c("status=FAILED", paste0("end=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), "detail=Official recovery exit or output audit failed; no automatic rerun"),
      status_file
    )
    stop("MSCCR Control bulk recovery failed for ", region, call. = FALSE)
  }

  writeLines(
    c("status=ACCEPTED", paste0("end=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), "detail=Official bulk recovery completed with 20 cell types and ecotype outputs"),
    status_file
  )
  cat("MSCCR_CONTROL_BULK_RECOVERY_ACCEPTED", region, "\n")
}

# Step4. 保存总体完成状态
region_audits <- rbindlist(lapply(jobs$RegionJobDir, function(path) fread(file.path(path, "output_audit.tsv"))))
fwrite(region_audits, file.path(job_root, "MSCCR_Control_bulk_recovery_output_audit.tsv"), sep = "\t", quote = FALSE)
completion <- list(
  status = "ACCEPTED_MSCCR_CONTROL_BULK_RECOVERY_CURRENT85",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  regions = regions,
  samples = c(Ileum = 121L, Colon = 115L, Rectum = 225L),
  discovery = discovery_name,
  threads_per_region = threads,
  execution = "sequential by region",
  scale = "linear CPM; full raw-library denominator; no TPM rebuild",
  role = "internal nonIBD coverage only; discovery CD/UC samples were not rerun"
)
jsonlite::write_json(completion, file.path(job_root, "MSCCR_Control_bulk_recovery_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(session_info_text, file.path(job_root, "R_sessionInfo.txt"))
print(completion)
cat("ACCEPTED_MSCCR_CONTROL_BULK_RECOVERY_CURRENT85\n")
