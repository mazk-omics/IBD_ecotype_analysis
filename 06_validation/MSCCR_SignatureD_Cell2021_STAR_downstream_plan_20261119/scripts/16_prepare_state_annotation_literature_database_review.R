#!/usr/bin/env Rscript
# 本脚本为 85 个正式 state 建立可追溯的文献与数据库注释证据表。
# 当前冻结标签只作为待审对象；本脚本不改 state、ecotype、marker 或原注释结果。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入和输出
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 16_prepare_state_annotation_literature_database_review.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

annotation_file <- file.path(RUN_ROOT, "output/10_state_annotation/state_annotation_final_table.tsv")
source_file <- file.path(RUN_ROOT, "inputs/state_annotation_source_register.tsv")
database_summary_file <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review/database_overlap/state_annotation_database_overlap_summary.tsv")
output_dir <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(annotation_file)) stop("Missing final state annotation table", call. = FALSE)
if (!file.exists(source_file)) stop("Missing annotation source register", call. = FALSE)
if (!file.exists(database_summary_file)) stop("Run database cross-reference before preparing the review table", call. = FALSE)

annotation <- fread(annotation_file, colClasses = "character", check.names = FALSE)
sources <- fread(source_file, colClasses = "character", check.names = FALSE)
database_summary <- fread(database_summary_file, colClasses = "character", check.names = FALSE)

required_annotation <- c("CellState", "CellType", "State", "Ecotype", "EcotypeMember", "FinalStateLabel", "FinalStateLabel_CN", "AnnotationConfidence", "StateEvidenceTier", "KeyMarkers", "ExternalValidatedMarkers", "KeyPathways", "ValidationSummary", "AnnotationBasis", "AnnotationCaveat")
required_sources <- c("SourceID", "SourceType", "SourceName", "VersionOrYear", "Citation", "DOIorAccession", "URL", "CurrentUse", "Scope", "Limit", "AccessedDate")
if (!all(required_annotation %in% names(annotation))) stop("Final annotation columns changed", call. = FALSE)
if (!all(required_sources %in% names(sources))) stop("Source register columns changed", call. = FALSE)
if (nrow(annotation) != 85L || anyDuplicated(annotation$CellState)) stop("Expected 85 unique states", call. = FALSE)
if (anyDuplicated(sources$SourceID) || any(!nzchar(sources$SourceID))) stop("Source IDs must be unique and non-empty", call. = FALSE)
if (nrow(database_summary) != 255L || uniqueN(database_summary$CellState) != 85L) stop("Expected 85 states by three databases", call. = FALSE)

print(dim(annotation))
print(table(annotation$CellType))
print(table(annotation$AnnotationConfidence))
print(table(sources$CurrentUse))

# Step2. 汇总逐 state 的数据库证据
database_source_ids <- c(CellMarker2 = "DB2", PanglaoDB = "DB3", HPA26 = "DB4")
database_review_parts <- vector("list", nrow(annotation))
for (state_index in seq_len(nrow(annotation))) {
  current_state <- annotation$CellState[[state_index]]
  current_summary <- database_summary[CellState == current_state]
  preferred_text <- current_summary$PreferredContextMarkers[nzchar(current_summary$PreferredContextMarkers)]
  preferred_markers <- unique(unlist(strsplit(preferred_text, ";", fixed = TRUE)))
  preferred_markers <- sort(preferred_markers[nzchar(preferred_markers)])
  used_databases <- current_summary$Database[as.integer(current_summary$DatabaseMatchedMarkerN) > 0L]
  evidence_summary <- paste0(
    current_summary$Database, ": all=", current_summary$DatabaseMatchedMarkerN,
    ", preferred=", current_summary$PreferredContextMarkerN
  )
  database_review_parts[[state_index]] <- data.table(
    CellState = current_state,
    DatabaseMatchedMarkers = paste(preferred_markers, collapse = ";"),
    DatabaseSourceIDs = paste(unname(database_source_ids[used_databases]), collapse = ";"),
    DatabaseEvidenceSummary = paste(evidence_summary, collapse = " | ")
  )
}
database_review <- rbindlist(database_review_parts)

# Step3. 建立逐 state 的人工审阅表；旧标签只作为待审对象
review <- annotation[, .(
  CellState, CellType, State, Ecotype, EcotypeMember,
  FrozenLabel = FinalStateLabel, FrozenLabel_CN = FinalStateLabel_CN,
  FrozenAnnotationConfidence = AnnotationConfidence, StateEvidenceTier,
  KeyMarkers, ExternalValidatedMarkers, KeyPathways, ValidationSummary,
  FrozenAnnotationBasis = AnnotationBasis, FrozenAnnotationCaveat = AnnotationCaveat
)]
review <- merge(review, database_review, by = "CellState", all.x = TRUE, sort = FALSE)

review[, `:=`(
  ProposedLabel = "",
  ProposedLabel_CN = "",
  LabelClass = "",
  CanonicalPositiveMarkers = "",
  CanonicalNegativeOrMissingMarkers = "",
  LiteratureMatchedMarkers = "",
  LiteratureSourceIDs = "",
  ReferenceCellTypeOrState = "",
  TissueAndRegionContext = "",
  DiseaseAndInflammationContext = "",
  EvidenceDirection = "",
  LiteratureEvidenceGrade = "",
  ConflictSummary = "",
  ReviewerDecision = "PENDING",
  ReviewerNote = "",
  DatabaseEvidenceDetailFile = "database_overlap/state_annotation_database_overlap_summary.tsv"
)]

setorder(review, CellType, State)
fwrite(review, file.path(output_dir, "state_annotation_literature_database_review.tsv"), sep = "\t", quote = FALSE, na = "")

# Step4. 保存来源快照和输入校验信息
source_snapshot <- copy(sources)
source_snapshot[, SourceRegisterSHA256 := digest::digest(source_file, file = TRUE, algo = "sha256")]
fwrite(source_snapshot, file.path(output_dir, "state_annotation_source_register_snapshot.tsv"), sep = "\t", quote = FALSE, na = "")

manifest <- data.table(
  Input = c("FrozenStateAnnotation", "AnnotationSourceRegister", "DatabaseOverlapSummary"),
  Path = c(annotation_file, source_file, database_summary_file),
  Bytes = file.info(c(annotation_file, source_file, database_summary_file))$size,
  SHA256 = vapply(c(annotation_file, source_file, database_summary_file), digest::digest, character(1L), file = TRUE, algo = "sha256")
)
fwrite(manifest, file.path(output_dir, "state_annotation_literature_review_input_manifest.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_STATE_ANNOTATION_LITERATURE_REVIEW_SCAFFOLD",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  states = nrow(review),
  sources = nrow(sources),
  database_summary_rows = nrow(database_summary),
  pending_states = sum(review$ReviewerDecision == "PENDING"),
  frozen_annotation_changed = FALSE
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_literature_review_scaffold_complete.json"), pretty = TRUE, auto_unbox = TRUE)

print(completion)
cat("PASS_STATE_ANNOTATION_LITERATURE_REVIEW_SCAFFOLD\n")
