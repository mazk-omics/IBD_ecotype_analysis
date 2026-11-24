#!/usr/bin/env Rscript
# 本脚本按 Cell 2021 的解释口径区分文献锚定状态、marker 驱动功能程序和未解析状态。
# 它不改变 A-D 证据等级、建议标签、冻结 state 或 ecotype。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 20_align_annotation_interpretability_with_cell2021.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

review_file <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review/state_annotation_literature_database_review_v2.tsv")
source_file <- file.path(RUN_ROOT, "inputs/state_annotation_source_register.tsv")
output_dir <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review")
if (!file.exists(review_file) || !file.exists(source_file)) stop("Missing complete review or source register", call. = FALSE)

review <- fread(review_file, colClasses = "character", check.names = FALSE)
if (nrow(review) != 85L || anyDuplicated(review$CellState)) stop("Expected 85 unique reviewed states", call. = FALSE)
if (any(!review$LiteratureEvidenceGrade %in% c("A", "B", "C", "D"))) stop("Unexpected evidence grade", call. = FALSE)

# Step2. 采用与 Cell 2021 一致的解释层次
manuscript <- review[, .(
  CellState, CellType, State, Ecotype, EcotypeMember,
  ManuscriptLabel = ProposedLabel, ManuscriptLabel_CN = ProposedLabel_CN,
  LabelClass, LiteratureEvidenceGrade, ReviewerDecision, FrozenLabel,
  KeyMarkers, ExternalValidatedMarkers, DatabaseMatchedMarkers, LiteratureMatchedMarkers,
  DatabaseSourceIDs, LiteratureSourceIDs, EvidenceDirection, ConflictSummary, ReviewerNote
)]
manuscript[, InterpretationClass := fifelse(
  LiteratureEvidenceGrade %in% c("A", "B"), "LITERATURE_ANCHORED_NAMED_STATE",
  fifelse(LiteratureEvidenceGrade == "C", "MARKER_LED_FUNCTIONAL_PROGRAM", "UNRESOLVED_OR_CROSS_LINEAGE")
)]
manuscript[, InterpretableForManuscript := LiteratureEvidenceGrade %in% c("A", "B", "C")]
manuscript[, RecommendedUse := fifelse(
  LiteratureEvidenceGrade %in% c("A", "B"), "MAIN_TEXT_OR_SUPPLEMENT",
  fifelse(LiteratureEvidenceGrade == "C", "MAIN_TEXT_IF_ECOTYPE_OR_CLINICALLY_RELEVANT", "SUPPLEMENT_ONLY")
)]
setcolorder(manuscript, c(
  "CellState", "CellType", "State", "Ecotype", "EcotypeMember", "InterpretationClass",
  "InterpretableForManuscript", "RecommendedUse", "ManuscriptLabel", "ManuscriptLabel_CN",
  "LabelClass", "LiteratureEvidenceGrade", "ReviewerDecision", "FrozenLabel", "KeyMarkers",
  "ExternalValidatedMarkers", "DatabaseMatchedMarkers", "LiteratureMatchedMarkers",
  "DatabaseSourceIDs", "LiteratureSourceIDs", "EvidenceDirection", "ConflictSummary", "ReviewerNote"
))
fwrite(manuscript, file.path(output_dir, "state_annotation_manuscript_label_table_v3_cell2021_aligned.tsv"), sep = "\t", quote = FALSE, na = "")

# Step3. 保存原文对照，不把不同验证设计混为同一指标
benchmark <- data.table(
  Metric = c(
    "Candidate states before doublet exclusion", "States retained for downstream analysis",
    "States significantly recoverable in scRNA-seq", "States attributable to prior genes or phenotypes",
    "IBD frozen states", "IBD literature-anchored named states", "IBD marker-led functional programs",
    "IBD total interpretable states", "IBD unresolved or cross-lineage states"
  ),
  Cell2021 = c("71", "69", "67/71 (94%)", "Approximately two-thirds of 69", "NA", "NA", "NA", "NA", "NA"),
  IBD_SignatureD = c("NA", "NA", "Not directly comparable here", "NA", "85", "25/85 (29.4%)", "30/85 (35.3%)", "55/85 (64.7%)", "30/85 (35.3%)"),
  Interpretation = c(
    "Cell 2021 evaluated 71 states before excluding two potential-doublet states",
    "All 69 retained states entered downstream clinical and ecotype analyses",
    "Cell 2021 used reference-guided permutation recovery across seven tumor atlases",
    "Cell 2021 allowed known programs, prominent markers and prior scRNA states as annotation evidence",
    "Current frozen IBD model", "Direct or near-direct literature anchoring", "Biologically interpretable program without a definitive subtype identity",
    "Comparable interpretation denominator to the approximately two-thirds statement in Cell 2021", "Not suitable for specific biological naming"
  )
)
fwrite(benchmark, file.path(output_dir, "cell2021_annotation_interpretability_benchmark.tsv"), sep = "\t", quote = FALSE)

summary <- manuscript[, .(States = .N), by = .(InterpretationClass, InterpretableForManuscript, RecommendedUse)]
setorder(summary, -InterpretableForManuscript, InterpretationClass)
fwrite(summary, file.path(output_dir, "state_annotation_manuscript_v3_summary.tsv"), sep = "\t", quote = FALSE)

manifest <- data.table(
  Input = c("CompleteLiteratureReviewV2", "SourceRegister"),
  Path = c(review_file, source_file),
  Bytes = file.info(c(review_file, source_file))$size,
  SHA256 = vapply(c(review_file, source_file), digest::digest, character(1L), file = TRUE, algo = "sha256")
)
fwrite(manifest, file.path(output_dir, "state_annotation_manuscript_v3_input_manifest.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_CELL2021_ALIGNED_STATE_INTERPRETABILITY",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), states = nrow(manuscript),
  literature_anchored = sum(manuscript$InterpretationClass == "LITERATURE_ANCHORED_NAMED_STATE"),
  marker_led_program = sum(manuscript$InterpretationClass == "MARKER_LED_FUNCTIONAL_PROGRAM"),
  interpretable = sum(manuscript$InterpretableForManuscript), unresolved_or_cross_lineage = sum(!manuscript$InterpretableForManuscript),
  evidence_grades_changed = FALSE, frozen_annotation_changed = FALSE
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_manuscript_v3_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(summary)
print(benchmark)
print(completion)
cat("PASS_CELL2021_ALIGNED_STATE_INTERPRETABILITY\n")
