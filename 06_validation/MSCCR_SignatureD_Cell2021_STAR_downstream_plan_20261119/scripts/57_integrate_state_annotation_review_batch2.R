#!/usr/bin/env Rscript
# 本脚本冻结第二批 D 级未解析状态和 C 级功能程序的使用边界。
# 只写审阅决定，不改标签、证据等级、state 或 ecotype 定义。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取完整证据表和第二批审阅记录
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 57_integrate_state_annotation_review_batch2.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

master_file <- file.path(RUN_ROOT, "output/43_state_annotation_inflammation_sensitivity/state_annotation_evidence_master_v4_with_inflammation_sensitivity.tsv")
review_file <- file.path(RUN_ROOT, "inputs/state_annotation_priority_review_v4_batch2.tsv")
output_dir <- file.path(RUN_ROOT, "output/44_state_annotation_review_v4_batch2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(master_file) || !file.exists(review_file)) stop("Annotation master or batch2 review is missing", call. = FALSE)
master <- fread(master_file)
review <- fread(review_file)

if (nrow(master) != 85L || anyDuplicated(master$CellState)) stop("Annotation master is not 85 unique states", call. = FALSE)
if (nrow(review) != 26L || anyDuplicated(review$CellState)) stop("Batch2 review is not 26 unique states", call. = FALSE)
expected <- master[ReviewPriorityTier %in% c("P2_UNRESOLVED_REPRODUCIBLE", "P3_FUNCTIONAL_PROGRAM_ECOTYPE_REVIEW"), sort(CellState)]
if (!identical(sort(review$CellState), expected)) stop("Batch2 review does not match the frozen P2 and P3 state set", call. = FALSE)

# Step2. 检查等级与决策一致，不允许静默升级
review_check <- merge(master[, .(CellState, FinalReviewedGrade)], review, by = "CellState", all.y = TRUE, sort = FALSE)
if (any(review_check$FinalReviewedGrade == "D" & !grepl("^KEEP_D_", review_check$ReviewDecision))) stop("A D-grade state has an invalid batch2 decision", call. = FALSE)
if (any(review_check$FinalReviewedGrade == "C" & review_check$ReviewDecision != "KEEP_C_FUNCTIONAL")) stop("A C-grade state has an invalid batch2 decision", call. = FALSE)

reviewed <- merge(master, review, by = "CellState", all.x = TRUE, sort = FALSE, suffixes = c("", ".batch2"))
reviewed[, Batch2InterpretationUse := fifelse(!is.na(FinalUse), FinalUse, fifelse(FinalReviewedGrade == "D", "SUPPLEMENT_ONLY", fifelse(FinalReviewedGrade == "C", "FUNCTIONAL_PROGRAM", "STANDARD")))]
reviewed[, Batch2SensitivityFlag := fifelse(!is.na(ModelSensitivityFlag), ModelSensitivityFlag, FALSE)]
reviewed[, AnnotationReviewComplete := TRUE]
reviewed[, ModelDefinitionChanged := FALSE]
setorder(reviewed, CellType, State)

# Step3. 保存完整表、决策表和 ecotype 提醒清单
batch2_decisions <- reviewed[!is.na(ReviewDecision.batch2)]
ecotype_flags <- reviewed[Batch2SensitivityFlag %in% TRUE, .(
  CellState, Ecotype, FinalReviewedGrade, FinalReviewedLabel, ReviewDecision = ReviewDecision.batch2,
  ReviewBoundary, InflammationDirection, PatientAwareSupport
)]
summary_table <- reviewed[, .(
  States = .N,
  EcotypeMembers = sum(EcotypeMember %in% TRUE),
  MainTextNamedStates = sum(MainTextNamedStateEligible %in% TRUE),
  FunctionalPrograms = sum(FunctionalProgramEligible %in% TRUE),
  SupplementOnlyStates = sum(SupplementOnly %in% TRUE),
  Batch2SensitivityFlags = sum(Batch2SensitivityFlag %in% TRUE)
), by = FinalReviewedGrade]

fwrite(reviewed, file.path(output_dir, "state_annotation_evidence_master_v4_reviewed_batch2.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(batch2_decisions, file.path(output_dir, "state_annotation_priority_batch2_decisions.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(ecotype_flags, file.path(output_dir, "ecotype_member_annotation_sensitivity_flags_batch2.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(summary_table, file.path(output_dir, "state_annotation_v4_batch2_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

completion <- list(
  status = "PASS_STATE_ANNOTATION_PRIORITY_REVIEW_V4_BATCH2",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  formal_states = nrow(reviewed),
  reviewed_batch2_states = nrow(batch2_decisions),
  label_or_grade_changes = 0L,
  model_definition_changed = FALSE
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_review_v4_batch2_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(summary_table)
print(batch2_decisions[, .(CellState, Ecotype, FinalReviewedGrade, ReviewDecision = ReviewDecision.batch2, FinalUse, ModelSensitivityFlag)])
cat("PASS_STATE_ANNOTATION_PRIORITY_REVIEW_V4_BATCH2\n")
