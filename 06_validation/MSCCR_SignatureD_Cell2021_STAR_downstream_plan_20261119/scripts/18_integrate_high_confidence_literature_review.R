#!/usr/bin/env Rscript
# 本脚本将首批高置信 state 的逐条文献审阅写入正式审阅表，不改冻结注释表。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 18_integrate_high_confidence_literature_review.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

review_file <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review/state_annotation_literature_database_review.tsv")
evidence_file <- file.path(RUN_ROOT, "inputs/state_annotation_literature_evidence_high_confidence.tsv")
source_file <- file.path(RUN_ROOT, "inputs/state_annotation_source_register.tsv")
output_dir <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review")
output_file <- file.path(output_dir, "state_annotation_literature_database_review_v1.tsv")

input_files <- c(review_file, evidence_file, source_file)
if (any(!file.exists(input_files))) stop("Missing literature review input", call. = FALSE)

review <- fread(review_file, colClasses = "character", check.names = FALSE)
evidence <- fread(evidence_file, colClasses = "character", check.names = FALSE)
sources <- fread(source_file, colClasses = "character", check.names = FALSE)
if (nrow(review) != 85L || anyDuplicated(review$CellState)) stop("Expected 85 unique review states", call. = FALSE)
if (nrow(evidence) != 26L || anyDuplicated(evidence$CellState)) stop("Expected 26 unique high-confidence states", call. = FALSE)
if (!all(evidence$CellState %in% review$CellState)) stop("Literature evidence contains unknown states", call. = FALSE)

# Step2. 检查文献来源和冻结置信度
source_tokens <- unique(unlist(strsplit(evidence$LiteratureSourceIDs, ";", fixed = TRUE)))
source_tokens <- source_tokens[nzchar(source_tokens)]
if (!all(source_tokens %in% sources$SourceID)) stop("Literature evidence contains unknown source IDs", call. = FALSE)
if (!all(review[CellState %in% evidence$CellState, FrozenAnnotationConfidence] == "HIGH")) stop("Evidence set is not the frozen high-confidence set", call. = FALSE)

# Step3. 逐列写入已审阅证据
evidence_columns <- setdiff(names(evidence), "CellState")
for (column_name in evidence_columns) {
  if (!column_name %in% names(review)) stop("Review table lacks evidence column: ", column_name, call. = FALSE)
  state_values <- setNames(evidence[[column_name]], evidence$CellState)
  selected <- review$CellState %in% evidence$CellState
  review[selected, (column_name) := state_values[CellState]]
}

setorder(review, CellType, State)
if (sum(review$ReviewerDecision == "PENDING") != 59L) stop("Unexpected pending-state count", call. = FALSE)
fwrite(review, output_file, sep = "\t", quote = FALSE, na = "")

# Step4. 输出审阅审计
audit <- review[, .N, by = .(FrozenAnnotationConfidence, ReviewerDecision, LiteratureEvidenceGrade, LabelClass)]
setorder(audit, FrozenAnnotationConfidence, ReviewerDecision, LiteratureEvidenceGrade, LabelClass)
fwrite(audit, file.path(output_dir, "state_annotation_literature_review_v1_audit.tsv"), sep = "\t", quote = FALSE)

manifest <- data.table(
  Input = c("ReviewScaffold", "HighConfidenceLiteratureEvidence", "SourceRegister"),
  Path = input_files,
  Bytes = file.info(input_files)$size,
  SHA256 = vapply(input_files, digest::digest, character(1L), file = TRUE, algo = "sha256")
)
fwrite(manifest, file.path(output_dir, "state_annotation_literature_review_v1_input_manifest.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_HIGH_CONFIDENCE_STATE_LITERATURE_REVIEW_V1",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), states = nrow(review), reviewed_states = nrow(evidence),
  keep = sum(review$ReviewerDecision == "KEEP"), revise = sum(review$ReviewerDecision == "REVISE"), pending = sum(review$ReviewerDecision == "PENDING"),
  frozen_annotation_changed = FALSE
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_literature_review_v1_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(audit)
print(completion)
cat("PASS_HIGH_CONFIDENCE_STATE_LITERATURE_REVIEW_V1\n")
