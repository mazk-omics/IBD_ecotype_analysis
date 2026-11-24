#!/usr/bin/env Rscript
# 本脚本将剩余 59 个 state 的文献审阅合并为完整 v2 表，不改冻结注释和模型。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 19_integrate_remaining_literature_review.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

review_file <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review/state_annotation_literature_database_review_v1.tsv")
evidence_file <- file.path(RUN_ROOT, "inputs/state_annotation_literature_evidence_remaining.tsv")
source_file <- file.path(RUN_ROOT, "inputs/state_annotation_source_register.tsv")
output_dir <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review")
output_file <- file.path(output_dir, "state_annotation_literature_database_review_v2.tsv")

input_files <- c(review_file, evidence_file, source_file)
if (any(!file.exists(input_files))) stop("Missing complete literature review input", call. = FALSE)

review <- fread(review_file, colClasses = "character", check.names = FALSE)
evidence <- fread(evidence_file, colClasses = "character", check.names = FALSE)
sources <- fread(source_file, colClasses = "character", check.names = FALSE)
if (nrow(review) != 85L || anyDuplicated(review$CellState)) stop("Expected 85 unique review states", call. = FALSE)
if (nrow(evidence) != 59L || anyDuplicated(evidence$CellState)) stop("Expected 59 unique remaining states", call. = FALSE)
if (!setequal(evidence$CellState, review[ReviewerDecision == "PENDING", CellState])) stop("Remaining evidence does not match pending states", call. = FALSE)

# Step2. 检查来源和证据字段
source_tokens <- unique(unlist(strsplit(evidence$LiteratureSourceIDs, ";", fixed = TRUE)))
source_tokens <- source_tokens[nzchar(source_tokens)]
if (!all(source_tokens %in% sources$SourceID)) stop("Literature evidence contains unknown source IDs", call. = FALSE)
required_values <- c("ProposedLabel", "ProposedLabel_CN", "LabelClass", "LiteratureSourceIDs", "EvidenceDirection", "LiteratureEvidenceGrade", "ReviewerDecision", "ReviewerNote")
if (any(!nzchar(as.matrix(evidence[, ..required_values])))) stop("Remaining evidence contains an empty required field", call. = FALSE)

# Step3. 写入剩余逐 state 证据
evidence_columns <- setdiff(names(evidence), "CellState")
for (column_name in evidence_columns) {
  if (!column_name %in% names(review)) stop("Review table lacks evidence column: ", column_name, call. = FALSE)
  state_values <- setNames(evidence[[column_name]], evidence$CellState)
  selected <- review$CellState %in% evidence$CellState
  review[selected, (column_name) := state_values[CellState]]
}

setorder(review, CellType, State)
if (any(review$ReviewerDecision == "PENDING")) stop("Literature review still contains pending states", call. = FALSE)
if (any(!review$LiteratureEvidenceGrade %in% c("A", "B", "C", "D"))) stop("Unexpected literature evidence grade", call. = FALSE)
fwrite(review, output_file, sep = "\t", quote = FALSE, na = "")

# Step4. 输出总体审阅审计
manuscript_table <- review[, .(
  CellState, CellType, State, Ecotype, EcotypeMember, FrozenLabel,
  ManuscriptLabel = ProposedLabel, ManuscriptLabel_CN = ProposedLabel_CN,
  LabelClass, LiteratureEvidenceGrade, ReviewerDecision,
  KeyMarkers, ExternalValidatedMarkers, DatabaseMatchedMarkers,
  LiteratureMatchedMarkers, DatabaseSourceIDs, LiteratureSourceIDs,
  EvidenceDirection, ConflictSummary, ReviewerNote
)]
manuscript_table[, ManuscriptUseTier := fifelse(
  LiteratureEvidenceGrade %in% c("A", "B"), "PRIMARY_INTERPRETABLE",
  fifelse(LiteratureEvidenceGrade == "C", "DESCRIPTIVE_PROGRAM_ONLY", "SUPPLEMENTARY_UNRESOLVED_OR_CONFLICT")
)]
setcolorder(manuscript_table, c(
  "CellState", "CellType", "State", "Ecotype", "EcotypeMember", "ManuscriptUseTier",
  "ManuscriptLabel", "ManuscriptLabel_CN", "LabelClass", "LiteratureEvidenceGrade", "ReviewerDecision",
  "FrozenLabel", "KeyMarkers", "ExternalValidatedMarkers", "DatabaseMatchedMarkers",
  "LiteratureMatchedMarkers", "DatabaseSourceIDs", "LiteratureSourceIDs", "EvidenceDirection", "ConflictSummary", "ReviewerNote"
))
fwrite(manuscript_table, file.path(output_dir, "state_annotation_manuscript_label_table_v2.tsv"), sep = "\t", quote = FALSE, na = "")

audit <- review[, .N, by = .(ReviewerDecision, LiteratureEvidenceGrade, LabelClass)]
setorder(audit, ReviewerDecision, LiteratureEvidenceGrade, LabelClass)
fwrite(audit, file.path(output_dir, "state_annotation_literature_review_v2_audit.tsv"), sep = "\t", quote = FALSE)

lineage_audit <- review[, .(
  States = .N, Keep = sum(ReviewerDecision == "KEEP"), Revise = sum(ReviewerDecision == "REVISE"),
  GradeA = sum(LiteratureEvidenceGrade == "A"), GradeB = sum(LiteratureEvidenceGrade == "B"),
  GradeC = sum(LiteratureEvidenceGrade == "C"), GradeD = sum(LiteratureEvidenceGrade == "D")
), by = CellType]
setorder(lineage_audit, CellType)
fwrite(lineage_audit, file.path(output_dir, "state_annotation_literature_review_v2_lineage_audit.tsv"), sep = "\t", quote = FALSE)

manifest <- data.table(
  Input = c("LiteratureReviewV1", "RemainingLiteratureEvidence", "SourceRegister"),
  Path = input_files,
  Bytes = file.info(input_files)$size,
  SHA256 = vapply(input_files, digest::digest, character(1L), file = TRUE, algo = "sha256")
)
fwrite(manifest, file.path(output_dir, "state_annotation_literature_review_v2_input_manifest.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_COMPLETE_STATE_LITERATURE_REVIEW_V2",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), states = nrow(review), reviewed_states = nrow(review),
  keep = sum(review$ReviewerDecision == "KEEP"), revise = sum(review$ReviewerDecision == "REVISE"), pending = sum(review$ReviewerDecision == "PENDING"),
  grade_a = sum(review$LiteratureEvidenceGrade == "A"), grade_b = sum(review$LiteratureEvidenceGrade == "B"),
  grade_c = sum(review$LiteratureEvidenceGrade == "C"), grade_d = sum(review$LiteratureEvidenceGrade == "D"), frozen_annotation_changed = FALSE
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_literature_review_v2_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(audit)
print(lineage_audit)
print(completion)
cat("PASS_COMPLETE_STATE_LITERATURE_REVIEW_V2\n")
