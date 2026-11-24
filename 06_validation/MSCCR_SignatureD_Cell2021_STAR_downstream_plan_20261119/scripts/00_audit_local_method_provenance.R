#!/usr/bin/env Rscript
# 本脚本只审计仓库内的正式运行脚本、Cell 2021 原文和 EcoTyper 官方源码。
# 在本地运行；不连接服务器，不读取正式数值结果，输出方法证据与当前阻断项。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定目录和依赖
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 00_audit_local_method_provenance.R RUN_ROOT PROJECT_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
PROJECT_ROOT <- normalizePath(args[[2]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)

audit_dir <- file.path(RUN_ROOT, "output/00_freeze/local_method_provenance")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

formal_script_root <- file.path(PROJECT_ROOT, "06_validation/MSCCR_SignatureD_IBD2025_TPM_20261101/02_ecotyper_bulk_20261101/scripts")
official_root <- normalizePath(file.path(PROJECT_ROOT, "../Code/Ecotyper_Reference_code/ecotyper-master"), mustWork = TRUE)
paper_xml <- file.path(PROJECT_ROOT, "06_validation/paper_aligned_recovery_characterization_20261012/00_audit/PMC8526411_bioc.xml")

source_files <- c(
  formal_prepare = file.path(formal_script_root, "08_prepare_tpm_ecotyper_formal.py"),
  formal_runner = file.path(formal_script_root, "11_run_tpm_ecotyper_formal.sh"),
  formal_resume = file.path(formal_script_root, "14_resume_ecotyper_after_zero_sd_fix.sh"),
  official_wrapper = file.path(official_root, "EcoTyper_discovery_bulk.R"),
  official_ecotypes = file.path(official_root, "pipeline/ecotypes.R"),
  official_assignment = file.path(official_root, "pipeline/ecotypes_assign_samples.R"),
  official_dropout = file.path(official_root, "pipeline/state_discovery_calculate_dropout_score.R"),
  official_second_filter = file.path(official_root, "pipeline/state_discovery_second_filter.R"),
  cell2021_xml = paper_xml
)
if (any(!file.exists(source_files))) stop("Missing method source: ", paste(source_files[!file.exists(source_files)], collapse = "; "), call. = FALSE)

source_text <- lapply(source_files, readLines, warn = FALSE, encoding = "UTF-8")
names(source_text) <- names(source_files)

# Step2. 定位正式配置与官方实现证据
evidence_patterns <- data.frame(
  EvidenceID = c(
    "E01_FORMAL_SKIP7", "E02_FORMAL_REGION3", "E03_FORMAL_NMF50",
    "E04_WRAPPER_STEP7", "E05_DROPOUT_OUTPUT", "E06_SECOND_FILTER_BACKUP",
    "E07_CODE_DISTANCE", "E08_STAR_DISTANCE", "E09_STAR_THREE_QC",
    "E10_ASSIGNMENT_Q025", "E11_RESUME_SKIPPED_ADVANCED_QC"
  ),
  Source = c(
    "formal_prepare", "formal_prepare", "formal_prepare",
    "official_wrapper", "official_dropout", "official_second_filter",
    "official_ecotypes", "cell2021_xml", "cell2021_xml",
    "official_assignment", "formal_resume"
  ),
  Pattern = c(
    "Pipeline steps to skip : [7]", "Annotation file column to scale by : \"Region3\"", "Number of NMF restarts : 50",
    "Step 7 (advanced cell state QC filter)", "dropout_scores_with_outliers.txt", "state_assignment_first_filter.txt",
    "hclust(as.dist(1-jaccard), method = \"average\")", "average linkage with Euclidean distance", "We applied three quality control filters",
    "clinical$AssignmentQ < 0.25", "跳过 Fractions、HiRes 和高级 QC"
  ),
  stringsAsFactors = FALSE
)

evidence_rows <- vector("list", nrow(evidence_patterns))
for (i in seq_len(nrow(evidence_patterns))) {
  source_id <- evidence_patterns$Source[[i]]
  pattern <- evidence_patterns$Pattern[[i]]
  hit <- grep(pattern, source_text[[source_id]], fixed = TRUE, useBytes = TRUE)
  evidence_rows[[i]] <- data.frame(
    EvidenceID = evidence_patterns$EvidenceID[[i]],
    Source = source_id,
    SourcePath = source_files[[source_id]],
    Pattern = pattern,
    Found = length(hit) > 0L,
    Line = if (length(hit) > 0L) hit[[1]] else NA_integer_,
    Text = if (length(hit) > 0L) trimws(source_text[[source_id]][hit[[1]]]) else NA_character_
  )
}
evidence <- data.table::rbindlist(evidence_rows)
print(evidence[, .(EvidenceID, Found, Line)])
if (any(!evidence$Found)) stop("Method evidence pattern missing: ", paste(evidence[Found == FALSE, EvidenceID], collapse = ", "), call. = FALSE)

# Step3. 写出方法判断，不把本地脚本证据冒充服务器数值验收
findings <- data.frame(
  FindingID = c(
    "F01_FORMAL_METHOD_IDENTITY", "F02_ADVANCED_QC_STATUS", "F03_STAR_QC_PARITY",
    "F04_ECOTYPE_DISTANCE_DISCREPANCY", "F05_OFFICIAL_ASSIGNMENT_RULE", "F06_NUMERIC_MODEL_FREEZE"
  ),
  Status = c("CONFIRMED", "CONFIRMED_SKIPPED", "NOT_YET_COMPLETE", "CONFIRMED", "CONFIRMED", "WAITING_SERVER_INPUTS"),
  Conclusion = c(
    "Formal preparation fixes Region3 scaling, 50 NMF restarts and skips pipeline step 7.",
    "Official wrapper defines step 7 as advanced dropout QC; the formal resume also states that advanced QC was skipped.",
    "Cell 2021 specifies three state-QC filters, so the accepted formal result is not yet a literal three-filter STAR Methods result.",
    "STAR Methods states Euclidean distance on filtered Jaccard rows, while official ecotypes.R uses as.dist(1-Jaccard). Both must be audited without overwriting the accepted model.",
    "Official bulk assignment requires winner q<0.25 and at least one valid ecotype-state assignment.",
    "The current local workspace contains plots and run scripts but not the complete accepted formal numeric model; hashes and matrix identities must be frozen through one authorized persistent server session."
  ),
  ChangesModel = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  AutomaticAction = c("NONE", "NONE", "NONE", "NONE", "NONE", "NONE"),
  stringsAsFactors = FALSE
)
print(findings)

# Step4. 保存源码身份和本地审计状态
manifest <- data.frame(
  Source = names(source_files),
  Path = normalizePath(source_files, mustWork = TRUE),
  Bytes = file.info(source_files)$size,
  SHA256 = vapply(source_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)

data.table::fwrite(evidence, file.path(audit_dir, "method_evidence.tsv"), sep = "\t")
data.table::fwrite(findings, file.path(audit_dir, "method_findings.tsv"), sep = "\t")
data.table::fwrite(manifest, file.path(audit_dir, "method_source_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo.txt"))

summary <- list(
  status = "PASS_LOCAL_METHOD_ONLY",
  server_connected = FALSE,
  numeric_model_frozen = FALSE,
  formal_step7_advanced_qc = "SKIPPED",
  literal_star_three_filter_qc = "NOT_YET_COMPLETE",
  model_changed = FALSE,
  automatic_rerun = FALSE,
  evidence_rows = nrow(evidence),
  findings = nrow(findings)
)
jsonlite::write_json(summary, file.path(audit_dir, "local_method_audit.json"), pretty = TRUE, auto_unbox = TRUE)
cat("LOCAL_METHOD_AUDIT_PASS", nrow(evidence), "evidence rows; formal step 7 skipped; no model changes\n")
