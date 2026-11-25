#!/usr/bin/env Rscript
# 本脚本将首批高优先级人工审阅决定写入 v4 注释表，不修改 state、ecotype 或验证结果。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取优先级主表和人工审阅记录
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 55_integrate_priority_state_annotation_review.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

master_file <- file.path(RUN_ROOT, "output/41_state_annotation_priority_master/state_annotation_evidence_master_v4.tsv")
review_file <- file.path(RUN_ROOT, "inputs/state_annotation_priority_review_v4.tsv")
output_dir <- file.path(RUN_ROOT, "output/42_state_annotation_review_v4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(master_file) || !file.exists(review_file)) stop("Priority master or review file is missing", call. = FALSE)
master <- fread(master_file)
review <- fread(review_file)

if (nrow(master) != 85L || anyDuplicated(master$CellState)) stop("Priority master is not 85 unique states", call. = FALSE)
if (nrow(review) != 10L || anyDuplicated(review$CellState)) stop("Priority review is not 10 unique states", call. = FALSE)
expected <- master[ReviewPriorityTier == "P1_UNRESOLVED_ECOTYPE_MEMBER_HIGH_EVIDENCE", sort(CellState)]
if (!identical(sort(review$CellState), expected)) stop("Priority review does not match the frozen P1 state set", call. = FALSE)
if (any(review$ReviewedGrade != "D")) stop("This review batch must not upgrade unresolved states", call. = FALSE)

# Step2. 合并审阅决定并冻结文章使用边界
reviewed <- merge(master, review, by = "CellState", all.x = TRUE, sort = FALSE)
reviewed[, FinalReviewedLabel := fifelse(!is.na(ReviewedLabel), ReviewedLabel, ManuscriptLabel)]
reviewed[, FinalReviewedLabel_CN := fifelse(!is.na(ReviewedLabel_CN), ReviewedLabel_CN, ManuscriptLabel_CN)]
reviewed[, FinalReviewedGrade := fifelse(!is.na(ReviewedGrade), ReviewedGrade, LiteratureEvidenceGrade)]
reviewed[, FinalInterpretationClass := fifelse(
  ReviewDecision == "REFINE_D_CROSS_LINEAGE",
  "REPRODUCIBLE_CROSS_LINEAGE_SIGNAL",
  fifelse(
    ReviewDecision == "KEEP_D_REFINE_LABEL",
    "UNRESOLVED_PARENT_ASSIGNED_PROGRAM",
    LiteratureInterpretationClass
  )
)]
reviewed[, MainTextNamedStateEligible := FinalReviewedGrade %in% c("A", "B")]
reviewed[, FunctionalProgramEligible := FinalReviewedGrade == "C"]
reviewed[, SupplementOnly := FinalReviewedGrade == "D"]
reviewed[, EcotypeInterpretationUse := fifelse(
  ReviewDecision == "REFINE_D_CROSS_LINEAGE",
  "SENSITIVITY_FLAG_CROSS_LINEAGE",
  fifelse(FinalReviewedGrade == "D", "SUPPLEMENT_ONLY", "STANDARD_ECOTYPE_INTERPRETATION")
)]
reviewed[, ModelDefinitionChanged := FALSE]
setorder(reviewed, CellType, State)

# Step3. 保存 v4 注释表和 ecotype 风险清单
reviewed_batch <- reviewed[!is.na(ReviewDecision)]
ecotype_flags <- reviewed[EcotypeSensitivityFlag %in% TRUE, .(
  CellState, Ecotype, FinalReviewedLabel, FinalInterpretationClass,
  EcotypeInterpretationUse, InflammationDirection, InflamedOR, InflamedFDR,
  NoninflamedOR, NoninflamedFDR, InterpretationBoundary,
  LiteratureSourceIDs = LiteratureSourceIDs.y
)]
summary_table <- reviewed[, .(
  States = .N,
  EcotypeMembers = sum(EcotypeMember %in% TRUE),
  MainTextNamedStateEligible = sum(MainTextNamedStateEligible),
  FunctionalProgramEligible = sum(FunctionalProgramEligible),
  SupplementOnly = sum(SupplementOnly),
  CrossLineageSignals = sum(FinalInterpretationClass == "REPRODUCIBLE_CROSS_LINEAGE_SIGNAL", na.rm = TRUE)
), by = FinalReviewedGrade]

fwrite(reviewed, file.path(output_dir, "state_annotation_evidence_master_v4_reviewed.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(reviewed_batch, file.path(output_dir, "state_annotation_priority_batch1_decisions.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(ecotype_flags, file.path(output_dir, "ecotype_member_annotation_sensitivity_flags.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(summary_table, file.path(output_dir, "state_annotation_v4_review_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

completion <- list(
  status = "PASS_STATE_ANNOTATION_PRIORITY_REVIEW_V4_BATCH1",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  formal_states = nrow(reviewed),
  reviewed_p1_states = nrow(reviewed_batch),
  label_upgrades = 0L,
  cross_lineage_signals_refined = sum(reviewed_batch$FinalInterpretationClass == "REPRODUCIBLE_CROSS_LINEAGE_SIGNAL"),
  model_definition_changed = FALSE
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_review_v4_batch1_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(summary_table)
print(reviewed_batch[, .(
  CellState, Ecotype, FinalReviewedLabel, FinalReviewedGrade,
  FinalInterpretationClass, EcotypeInterpretationUse, InflammationDirection
)])
cat("PASS_STATE_ANNOTATION_PRIORITY_REVIEW_V4_BATCH1\n")
