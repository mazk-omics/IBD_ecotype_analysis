#!/usr/bin/env Rscript
# 本脚本只读验收已完成的 MSCCR Control bulk recovery，不重复运行 EcoTyper。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定结果目录和验收标准
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 49_finalize_MSCCR_Control_bulk_recovery_current85.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

job_root <- file.path(RUN_ROOT, "output/36_MSCCR_Control_bulk_recovery_current85_retry1")
regions <- c("Ileum", "Colon", "Rectum")
expected_samples <- c(Ileum = 121L, Colon = 115L, Rectum = 225L)
completion_file <- file.path(job_root, "MSCCR_Control_bulk_recovery_complete.json")
if (!file.exists(completion_file)) stop("MSCCR Control bulk recovery completion JSON is missing", call. = FALSE)

# Step2. 逐部位复核退出码、20 个父类型和 ecotype 输出
region_audits <- rbindlist(lapply(regions, function(region) {
  status_file <- file.path(job_root, region, "status.txt")
  audit_file <- file.path(job_root, region, "output_audit.tsv")
  if (!file.exists(status_file) || !file.exists(audit_file)) stop("Missing region audit: ", region, call. = FALSE)
  status_text <- paste(readLines(status_file, warn = FALSE), collapse = "; ")
  audit <- fread(audit_file)
  audit[, StatusText := status_text]
  audit
}))

accepted <- region_audits$ExitStatus == 0L & region_audits$StateAbundanceFiles == 20L &
  region_audits$StateAssignmentFiles == 20L & region_audits$EcotypeAbundancePresent &
  region_audits$EcotypeAssignmentPresent & grepl("status=ACCEPTED", region_audits$StatusText, fixed = TRUE)
if (!all(accepted)) stop("One or more MSCCR Control region audits did not pass", call. = FALSE)

region_audits[, RegionOrder := match(Region3, regions)]
setorder(region_audits, RegionOrder)
region_audits[, RegionOrder := NULL]
fwrite(region_audits, file.path(job_root, "MSCCR_Control_bulk_recovery_final_audit.tsv"), sep = "\t", quote = FALSE)

# Step3. 冻结最终验收，不改动已有 recovery 结果
final_acceptance <- list(
  status = "ACCEPTED_MSCCR_CONTROL_BULK_RECOVERY_CURRENT85",
  finalized_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  regions = regions,
  samples = expected_samples,
  state_abundance_files = sum(region_audits$StateAbundanceFiles),
  state_assignment_files = sum(region_audits$StateAssignmentFiles),
  ecotype_abundance_files = sum(region_audits$EcotypeAbundancePresent),
  ecotype_assignment_files = sum(region_audits$EcotypeAssignmentPresent),
  rerun_performed = FALSE,
  scale = "linear CPM; full raw-library denominator; no TPM rebuild",
  role = "internal nonIBD coverage only; discovery CD/UC samples were not rerun"
)
jsonlite::write_json(final_acceptance, file.path(job_root, "MSCCR_Control_bulk_recovery_final_acceptance.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(job_root, "R_sessionInfo_finalize.txt"))
print(final_acceptance)
cat("ACCEPTED_MSCCR_CONTROL_BULK_RECOVERY_CURRENT85\n")
