#!/usr/bin/env Rscript
# 本脚本审计 IBDverse_E_MTAB_16999 官方 recovery 输出和正式 state 覆盖，不重算状态。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞和已完成任务
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 07d_audit_IBDverse_official_recovery.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

input_audit_dir <- file.path(output_root, "03_scrna_input_audit")
run_log_dir <- file.path(RUN_ROOT, "output/07c_IBDverse_official_recovery")
audit_dir <- file.path(RUN_ROOT, "output/07d_IBDverse_official_recovery_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

selected_cells <- fread(file.path(input_audit_dir, "recovery_selected_cells.tsv.gz"), colClasses = "character")
selected_cells <- selected_cells[Dataset == "IBDverse_E_MTAB_16999"]
ecotype_members <- fread(file.path(output_root, "01_state_qc_audit/ecotype_member_dropout_impact.tsv"), colClasses = "character")
state_evaluability <- fread(file.path(input_audit_dir, "dataset_state_evaluability.tsv"), colClasses = "character")
evaluable_member_ids <- state_evaluability[Dataset == "IBDverse_E_MTAB_16999" & Evaluable == "TRUE", ID]
jobs <- fread(file.path(run_log_dir, "IBDverse_E_MTAB_16999_official_recovery_jobs_completed.tsv"), colClasses = "character")
jobs[, ExitStatus := as.integer(ExitStatus)]
if (nrow(jobs) != 16L || any(jobs$ExitStatus != 0L)) stop("Official IBDverse_E_MTAB_16999 recovery execution did not pass", call. = FALSE)

required_files <- c(
  "H.txt",
  "annotation.txt",
  "expression_matrix.txt",
  "expression_matrix_scaled.txt",
  "heatmap_data.txt",
  "heatmap_top_ann.txt",
  "initial_state_abundances.txt",
  "initial_state_assignment.txt",
  "nmf.RData",
  "state_abundances.txt",
  "state_assignment.txt",
  "state_assignment_heatmap.pdf",
  "state_assignment_heatmap.png"
)

# Step2. 逐父类型核验初始分配和正式分配
audit_parts <- vector("list", nrow(jobs))
state_count_parts <- vector("list", nrow(jobs))
for (row_index in seq_len(nrow(jobs))) {
  cell_type <- jobs$CellType[[row_index]]
  result_dir <- jobs$ResultDir[[row_index]]
  frozen_ids <- selected_cells[CellType == cell_type, CellID]
  result_files <- file.path(result_dir, required_files)
  if (any(!file.exists(result_files)) || any(file.info(result_files)$size <= 0)) stop("Required official recovery files are missing: ", cell_type, call. = FALSE)

  mapping_file <- file.path(state_discovery_root, cell_type, jobs$FrozenRank[[row_index]], "mapping_to_initial_states.txt")
  mapping <- fread(mapping_file, colClasses = "character")
  initial_assignment <- fread(file.path(result_dir, "initial_state_assignment.txt"), colClasses = "character")
  formal_assignment <- fread(file.path(result_dir, "state_assignment.txt"), colClasses = "character")
  initial_abundance <- fread(file.path(result_dir, "initial_state_abundances.txt"), check.names = FALSE)
  formal_abundance <- fread(file.path(result_dir, "state_abundances.txt"), check.names = FALSE)

  if (anyDuplicated(mapping$InitialState) || anyDuplicated(mapping$State)) stop("Formal state mapping is not one-to-one: ", cell_type, call. = FALSE)
  if (anyDuplicated(initial_assignment$ID) || !setequal(initial_assignment$ID, frozen_ids)) stop("Initial assignment does not cover frozen cells exactly: ", cell_type, call. = FALSE)
  if (anyDuplicated(formal_assignment$ID) || !all(formal_assignment$ID %in% frozen_ids)) stop("Formal assignment IDs are invalid: ", cell_type, call. = FALSE)
  expected_state <- mapping$State[match(formal_assignment$InitialState, mapping$InitialState)]
  if (anyNA(expected_state) || !identical(formal_assignment$State, expected_state)) stop("Formal state mapping changed: ", cell_type, call. = FALSE)

  initial_states <- initial_abundance[[1L]]
  initial_values <- as.matrix(initial_abundance[, -1L])
  storage.mode(initial_values) <- "double"
  formal_states <- formal_abundance[[1L]]
  formal_values <- as.matrix(formal_abundance[, -1L])
  storage.mode(formal_values) <- "double"
  if (!setequal(colnames(initial_values), frozen_ids) || !setequal(colnames(formal_values), frozen_ids)) stop("Abundance cell IDs differ from frozen cells: ", cell_type, call. = FALSE)
  if (!identical(formal_states, mapping$State)) stop("Formal abundance rows differ from frozen mapping: ", cell_type, call. = FALSE)
  if (any(!is.finite(initial_values)) || any(initial_values < 0) || max(abs(colSums(initial_values) - 1)) > 1e-6) stop("Initial state abundances are invalid: ", cell_type, call. = FALSE)
  if (any(!is.finite(formal_values)) || any(formal_values < 0) || any(colSums(formal_values) > 1 + 1e-6)) stop("Formal state abundances are invalid: ", cell_type, call. = FALSE)

  formal_counts <- table(factor(formal_assignment$State, levels = mapping$State))
  member_type <- ecotype_members[CellType == cell_type]
  is_ecotype_member <- mapping$State %in% member_type$State
  ecotype_lookup <- setNames(member_type$Ecotype, member_type$State)
  state_count_parts[[row_index]] <- data.frame(
    Dataset = "IBDverse_E_MTAB_16999",
    CellType = cell_type,
    State = mapping$State,
    InitialState = mapping$InitialState,
    EcotypeMember = is_ecotype_member,
    Ecotype = unname(ecotype_lookup[mapping$State]),
    AssignedCells = as.integer(formal_counts),
    stringsAsFactors = FALSE
  )
  audit_parts[[row_index]] <- data.frame(
    Dataset = "IBDverse_E_MTAB_16999",
    CellType = cell_type,
    FrozenRank = jobs$FrozenRank[[row_index]],
    FrozenCells = length(frozen_ids),
    InitialAssignedCells = nrow(initial_assignment),
    FormalAssignedCells = nrow(formal_assignment),
    FormalAssignmentFraction = nrow(formal_assignment) / length(frozen_ids),
    InitialStates = length(initial_states),
    CandidateStates = nrow(mapping),
    RepresentedCandidateStates = sum(formal_counts > 0),
    EvaluableEcotypeMemberStates = sum(is_ecotype_member),
    RepresentedEcotypeMemberStates = sum(formal_counts[is_ecotype_member] > 0),
    ResultFiles = length(required_files),
    StateAssignmentSHA256 = digest::digest(file.path(result_dir, "state_assignment.txt"), algo = "sha256", file = TRUE, serialize = FALSE),
    StateAbundanceSHA256 = digest::digest(file.path(result_dir, "state_abundances.txt"), algo = "sha256", file = TRUE, serialize = FALSE),
    HeatmapPNGSHA256 = digest::digest(file.path(result_dir, "state_assignment_heatmap.png"), algo = "sha256", file = TRUE, serialize = FALSE),
    LowInformationFlag = length(frozen_ids) < 50L || sum(formal_counts > 0) < nrow(mapping),
    stringsAsFactors = FALSE
  )
  cat("IBDverse_E_MTAB_16999_RECOVERY_AUDITED", cell_type, nrow(formal_assignment), "of", length(frozen_ids), "cells", sum(formal_counts > 0), "of", nrow(mapping), "formal states\n")
}

# Step3. 保存结果覆盖和审计结论
recovery_audit <- rbindlist(audit_parts)
state_counts <- rbindlist(state_count_parts)
present_member_ids <- state_counts[EcotypeMember == TRUE, paste(CellType, State, sep = "_")]
if (sum(recovery_audit$FrozenCells) != 35899L) stop("Frozen IBDverse_E_MTAB_16999 cell total changed", call. = FALSE)
if (sum(recovery_audit$CandidateStates) != 72L) stop("IBDverse_E_MTAB_16999 evaluable candidate state total changed", call. = FALSE)
if (nrow(ecotype_members) != 79L || length(evaluable_member_ids) != 68L || !setequal(present_member_ids, evaluable_member_ids)) stop("IBDverse_E_MTAB_16999 evaluable ecotype member set changed", call. = FALSE)
fwrite(recovery_audit, file.path(audit_dir, "IBDverse_E_MTAB_16999_official_recovery_audit.tsv"), sep = "\t")
fwrite(state_counts, file.path(audit_dir, "IBDverse_E_MTAB_16999_official_recovery_state_counts.tsv"), sep = "\t")
summary <- list(
  status = "PASS_IBDverse_E_MTAB_16999_OFFICIAL_RECOVERY_AUDIT",
  dataset = "IBDverse_E_MTAB_16999",
  parent_types = nrow(recovery_audit),
  frozen_cells = sum(recovery_audit$FrozenCells),
  initial_assigned_cells = sum(recovery_audit$InitialAssignedCells),
  formal_assigned_cells = sum(recovery_audit$FormalAssignedCells),
  formal_assignment_fraction = sum(recovery_audit$FormalAssignedCells) / sum(recovery_audit$FrozenCells),
  all_candidate_states = 85,
  evaluable_candidate_states = sum(recovery_audit$CandidateStates),
  represented_candidate_states = sum(recovery_audit$RepresentedCandidateStates),
  all_ecotype_member_states = nrow(ecotype_members),
  evaluable_ecotype_member_states = sum(recovery_audit$EvaluableEcotypeMemberStates),
  represented_ecotype_member_states = sum(recovery_audit$RepresentedEcotypeMemberStates),
  low_information_parent_types = recovery_audit[LowInformationFlag == TRUE, CellType],
  official_script_sha256 = "ccd7380a6962ac6905bff1c4890440991b130d574efd5bb69c52d7ebb4b2332d",
  builtin_significance = FALSE,
  recovery_significance_bootstraps_pending = 100
)
jsonlite::write_json(summary, file.path(audit_dir, "IBDverse_E_MTAB_16999_official_recovery_audit.json"), pretty = TRUE, auto_unbox = TRUE)
cat("IBDverse_E_MTAB_16999_OFFICIAL_RECOVERY_AUDIT_PASS", sum(recovery_audit$FormalAssignedCells), "formal assignments", sum(recovery_audit$RepresentedCandidateStates), "represented candidate states", sum(recovery_audit$RepresentedEcotypeMemberStates), "represented ecotype member states\n")
