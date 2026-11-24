#!/usr/bin/env Rscript
# 本脚本汇总四个单细胞 validation 数据集的 cell-state recovery 和 100 次 shuffle 显著性。
# 在服务器运行，读取各任务的正式输出，结果写入 output/08_scRNA_recovery_summary。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入与判定标准
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 08_summarize_scRNA_recovery.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

z_threshold <- 1.96
output_dir <- file.path(RUN_ROOT, "output/08_scRNA_recovery_summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_specs <- data.table(
  Dataset = c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999"),
  JobsFile = file.path(RUN_ROOT, c(
    "output/04e_GSE214695_recovery_significance_100/GSE214695_recovery_significance_100_jobs.tsv",
    "output/05e_SCP1884_recovery_significance_100/SCP1884_recovery_significance_100_jobs.tsv",
    "output/06e_SCP259_recovery_significance_100/SCP259_recovery_significance_100_jobs.tsv",
    "output/07e_IBDverse_recovery_significance_100/IBDverse_E_MTAB_16999_recovery_significance_100_jobs.tsv"
  )),
  PassFile = file.path(RUN_ROOT, c(
    "output/04e_GSE214695_recovery_significance_100/GSE214695_recovery_significance_100_summary.txt",
    "output/05e_SCP1884_recovery_significance_100/SCP1884_recovery_significance_100_summary.txt",
    "output/06e_SCP259_recovery_significance_100/SCP259_recovery_significance_100_summary.txt",
    "output/07e_IBDverse_recovery_significance_100/IBDverse_E_MTAB_16999_recovery_significance_100_summary.txt"
  )),
  PassPattern = c(
    "STATUS=PASS_GSE214695_RECOVERY_SIGNIFICANCE_100",
    "STATUS=PASS_SCP1884_RECOVERY_SIGNIFICANCE_100",
    "STATUS=PASS_SCP259_RECOVERY_SIGNIFICANCE_100",
    "STATUS=PASS_IBDverse_E_MTAB_16999_RECOVERY_SIGNIFICANCE_100"
  )
)

if (any(!file.exists(dataset_specs$JobsFile))) stop("One or more significance job manifests are missing", call. = FALSE)
if (any(!file.exists(dataset_specs$PassFile))) stop("One or more significance PASS summaries are missing", call. = FALSE)
pass_ok <- vapply(seq_len(nrow(dataset_specs)), function(i) any(readLines(dataset_specs$PassFile[[i]], warn = FALSE) == dataset_specs$PassPattern[[i]]), logical(1L))
if (any(!pass_ok)) stop("One or more significance runs did not pass", call. = FALSE)

# Step2. 逐数据集和 parent type 读取正式 recovery 输出
state_tables <- list()
parent_tables <- list()
table_index <- 0L

for (dataset_index in seq_len(nrow(dataset_specs))) {
  jobs <- fread(dataset_specs$JobsFile[[dataset_index]], colClasses = "character")
  if (any(jobs$ExitStatus != "0" | jobs$SignificanceExitStatus != "0")) stop("A recovery job did not exit successfully", call. = FALSE)
  print(jobs[, .N, by = .(Dataset)])

  for (job_index in seq_len(nrow(jobs))) {
    result_dir <- jobs$ResultDir[[job_index]]
    required_files <- file.path(result_dir, c("initial_state_assignment.txt", "state_assignment.txt", "recovery_z_scores_initial_states.txt", "recovery_z_scores.txt"))
    if (any(!file.exists(required_files))) stop(paste("Missing recovery output in", result_dir), call. = FALSE)

    initial_assignment <- fread(required_files[[1L]])
    formal_assignment <- fread(required_files[[2L]])
    initial_z <- fread(required_files[[3L]])
    formal_z <- fread(required_files[[4L]])
    annotation <- fread(jobs$AnnotationFile[[job_index]])
    if (!"Disease" %in% names(annotation)) stop(paste("Disease is missing from", jobs$AnnotationFile[[job_index]]), call. = FALSE)
    disease_levels <- sort(unique(annotation$Disease))
    healthy_level <- grepl("healthy|non[-_ ]?ibd|control|^hc$", disease_levels, ignore.case = TRUE)
    if (any(healthy_level)) stop(paste("Healthy/nonIBD cells entered the IBD-only recovery:", paste(disease_levels[healthy_level], collapse = ";")), call. = FALSE)

    initial_counts <- initial_assignment[, .(InitialAssignedCells = .N), by = .(InitialState = State)]
    formal_counts <- formal_assignment[, .(FormalAssignedCells = .N), by = .(State, InitialState)]
    if (formal_counts[, anyDuplicated(State)] != 0L) stop(paste("A formal state maps to multiple initial states in", result_dir), call. = FALSE)

    state_result <- merge(formal_z, formal_counts, by = "State", all.x = TRUE, sort = FALSE)
    state_result <- merge(state_result, initial_counts, by = "InitialState", all.x = TRUE, sort = FALSE)
    state_result[is.na(FormalAssignedCells), FormalAssignedCells := 0L]
    state_result[is.na(InitialAssignedCells), InitialAssignedCells := 0L]
    state_result[, Dataset := dataset_specs$Dataset[[dataset_index]]]
    state_result[, CellType := jobs$CellType[[job_index]]]
    state_result[, FrozenRank := as.integer(jobs$FrozenRank[[job_index]])]
    state_result[, InputCells := as.integer(jobs$Cells[[job_index]])]
    state_result[, FormalAssignmentFractionOfInput := FormalAssignedCells / InputCells]
    state_result[, HasFormalAssignment := FormalAssignedCells > 0L]
    state_result[, FiniteZ := is.finite(Z)]
    state_result[, RecoveryEvaluable := HasFormalAssignment & FiniteZ]
    state_result[, Significant_Zgt1.96 := FiniteZ & Z > z_threshold]
    state_result[, P_one_sided := ifelse(FiniteZ, pnorm(Z, lower.tail = FALSE), NA_real_)]
    setcolorder(state_result, c("Dataset", "CellType", "State", "InitialState", "FrozenRank", "InputCells", "InitialAssignedCells", "FormalAssignedCells", "FormalAssignmentFractionOfInput", "HasFormalAssignment", "Z", "RefMean", "BootstrapMean", "BootstrapSD", "FiniteZ", "RecoveryEvaluable", "Significant_Zgt1.96", "P_one_sided"))

    parent_result <- data.table(
      Dataset = dataset_specs$Dataset[[dataset_index]],
      CellType = jobs$CellType[[job_index]],
      DiseaseLevels = paste(disease_levels, collapse = ";"),
      FrozenRank = as.integer(jobs$FrozenRank[[job_index]]),
      InputCells = as.integer(jobs$Cells[[job_index]]),
      InitialAssignedCells = nrow(initial_assignment),
      FormalAssignedCells = nrow(formal_assignment),
      FormalAssignmentFraction = nrow(formal_assignment) / as.integer(jobs$Cells[[job_index]]),
      InitialStatesTested = nrow(initial_z),
      FormalStatesEvaluated = nrow(formal_z),
      StatesWithFormalAssignment = sum(state_result$HasFormalAssignment),
      FiniteZStates = sum(is.finite(formal_z$Z)),
      RecoveryEvaluableStates = sum(state_result$RecoveryEvaluable),
      SignificantStates_Zgt1.96 = sum(is.finite(formal_z$Z) & formal_z$Z > z_threshold),
      MedianZ = median(formal_z$Z[is.finite(formal_z$Z)]),
      MaxZ = max(formal_z$Z[is.finite(formal_z$Z)])
    )
    parent_result[, SignificantFractionFormalStates := SignificantStates_Zgt1.96 / FormalStatesEvaluated]
    parent_result[, SignificantFractionEvaluable := SignificantStates_Zgt1.96 / RecoveryEvaluableStates]

    table_index <- table_index + 1L
    state_tables[[table_index]] <- state_result
    parent_tables[[table_index]] <- parent_result
  }
}

state_by_dataset <- rbindlist(state_tables, use.names = TRUE, fill = TRUE)
parent_by_dataset <- rbindlist(parent_tables, use.names = TRUE, fill = TRUE)
state_by_dataset[, FDR_BH_within_dataset := p.adjust(P_one_sided, method = "BH"), by = Dataset]
state_by_dataset[, Significant_BH_within_dataset := FiniteZ & Z > 0 & is.finite(FDR_BH_within_dataset) & FDR_BH_within_dataset < 0.05]
setorder(state_by_dataset, Dataset, CellType, State)
setorder(parent_by_dataset, Dataset, CellType)

print(dim(state_by_dataset))
print(state_by_dataset[, .(States = .N, Significant = sum(Significant_Zgt1.96)), by = Dataset])

# Step3. 汇总数据集内恢复和显著恢复
dataset_summary <- parent_by_dataset[, .(
  ParentTypes = .N,
  InputCells = sum(InputCells),
  InitialAssignedCells = sum(InitialAssignedCells),
  FormalAssignedCells = sum(FormalAssignedCells),
  FormalAssignmentFraction = sum(FormalAssignedCells) / sum(InputCells),
  InitialStatesTested = sum(InitialStatesTested),
  FormalStatesEvaluated = sum(FormalStatesEvaluated),
  StatesWithFormalAssignment = sum(StatesWithFormalAssignment),
  FiniteZStates = sum(FiniteZStates),
  RecoveryEvaluableStates = sum(RecoveryEvaluableStates),
  SignificantStates_Zgt1.96 = sum(SignificantStates_Zgt1.96)
), by = Dataset]
dataset_summary[, SignificantFractionFormalStates := SignificantStates_Zgt1.96 / FormalStatesEvaluated]
dataset_summary[, SignificantFractionEvaluable := SignificantStates_Zgt1.96 / RecoveryEvaluableStates]

significant_states <- state_by_dataset[Significant_Zgt1.96 == TRUE]
state_across_datasets <- state_by_dataset[, .(
  DatasetsPresent = .N,
  DatasetsWithAssignment = sum(HasFormalAssignment),
  DatasetsEvaluable = sum(RecoveryEvaluable),
  DatasetsFiniteZ = sum(FiniteZ),
  DatasetsSignificant = sum(Significant_Zgt1.96),
  SignificantDatasets = paste(Dataset[Significant_Zgt1.96], collapse = ";"),
  TotalFormalAssignedCells = sum(FormalAssignedCells),
  MedianZ = median(Z[FiniteZ]),
  MaxZ = max(Z[FiniteZ])
), by = .(CellType, State)]
state_across_datasets[, SignificantInAtLeastOneDataset := DatasetsSignificant >= 1L]
state_across_datasets[, SignificantInAtLeastTwoDatasets := DatasetsSignificant >= 2L]
setorder(state_across_datasets, CellType, State)

# Step4. 按 Cell 2021 原法等权合并各数据集 Z
state_presence <- state_by_dataset[, .(
  DatasetsPresent = .N,
  PresentDatasets = paste(Dataset, collapse = ";")
), by = .(CellType, State)]

stouffer_meta <- state_by_dataset[RecoveryEvaluable == TRUE, .(
  DatasetsEvaluable = .N,
  EvaluableDatasets = paste(Dataset, collapse = ";"),
  IndividualZ = paste(sprintf("%s=%.6f", Dataset, Z), collapse = ";"),
  DatasetsIndividuallySignificant = sum(Significant_Zgt1.96),
  SumZ = sum(Z),
  MetaZ = sum(Z) / sqrt(.N)
), by = .(CellType, State)]
stouffer_meta <- merge(state_presence, stouffer_meta, by = c("CellType", "State"), all.x = TRUE, sort = FALSE)
stouffer_meta[, MetaP_one_sided := pnorm(MetaZ, lower.tail = FALSE)]
stouffer_meta[, MetaSignificant_Zgt1.96 := is.finite(MetaZ) & MetaZ > z_threshold]
stouffer_meta[, MetaFDR_BH := p.adjust(MetaP_one_sided, method = "BH")]
stouffer_meta[, MetaSignificant_BH := is.finite(MetaZ) & MetaZ > 0 & is.finite(MetaFDR_BH) & MetaFDR_BH < 0.05]
stouffer_meta[, EvidenceScope := ifelse(DatasetsEvaluable >= 2L, "MULTI_DATASET_STOUFFER", "SINGLE_DATASET_ONLY")]
setcolorder(stouffer_meta, c("CellType", "State", "DatasetsPresent", "PresentDatasets", "DatasetsEvaluable", "EvaluableDatasets", "IndividualZ", "DatasetsIndividuallySignificant", "SumZ", "MetaZ", "MetaP_one_sided", "MetaSignificant_Zgt1.96", "MetaFDR_BH", "MetaSignificant_BH", "EvidenceScope"))
setorder(stouffer_meta, CellType, State)

stouffer_significant <- stouffer_meta[MetaSignificant_Zgt1.96 == TRUE]
stouffer_nonsignificant <- stouffer_meta[MetaSignificant_Zgt1.96 == FALSE]
print(stouffer_meta[, .(States = .N, MetaSignificant = sum(MetaSignificant_Zgt1.96)), by = DatasetsEvaluable])

# Step5. 保存统计表和简要报告
fwrite(state_by_dataset, file.path(output_dir, "scRNA_recovery_state_by_dataset.tsv"), sep = "\t", na = "NA")
fwrite(significant_states, file.path(output_dir, "scRNA_recovery_significant_states.tsv"), sep = "\t", na = "NA")
fwrite(parent_by_dataset, file.path(output_dir, "scRNA_recovery_parent_by_dataset.tsv"), sep = "\t", na = "NA")
fwrite(dataset_summary, file.path(output_dir, "scRNA_recovery_dataset_summary.tsv"), sep = "\t", na = "NA")
fwrite(state_across_datasets, file.path(output_dir, "scRNA_recovery_state_across_datasets_descriptive.tsv"), sep = "\t", na = "NA")
fwrite(stouffer_meta, file.path(output_dir, "scRNA_recovery_stouffer_meta.tsv"), sep = "\t", na = "NA")
fwrite(stouffer_significant, file.path(output_dir, "scRNA_recovery_stouffer_significant_states.tsv"), sep = "\t", na = "NA")
fwrite(stouffer_nonsignificant, file.path(output_dir, "scRNA_recovery_stouffer_nonsignificant_states.tsv"), sep = "\t", na = "NA")

report_lines <- c(
  "# scRNA cell-state recovery 汇总",
  "",
  "本汇总以每个独立数据集内的正式 state recovery 为主。显著恢复按 Cell 2021 原文轨道定义为 100 次 shuffle 后 Z > 1.96；跨数据集重复仅作描述，不作为通过门槛。",
  "",
  "## 数据集汇总",
  "",
  paste0("- ", dataset_summary$Dataset, "：", dataset_summary$RecoveryEvaluableStates, "/", dataset_summary$FormalStatesEvaluated, " states 可评估恢复，其中 ", dataset_summary$SignificantStates_Zgt1.96, "/", dataset_summary$RecoveryEvaluableStates, " 显著（", sprintf("%.1f%%", 100 * dataset_summary$SignificantFractionEvaluable), "）；formal assignment 细胞覆盖 ", sprintf("%.1f%%", 100 * dataset_summary$FormalAssignmentFraction), "。"),
  "",
  "## 跨数据集描述",
  "",
  paste0("- 至少在 1 个数据集中显著恢复：", sum(state_across_datasets$SignificantInAtLeastOneDataset), "/", nrow(state_across_datasets), " 个 discovery states。"),
  paste0("- 至少在 2 个数据集中显著恢复：", sum(state_across_datasets$SignificantInAtLeastTwoDatasets), "/", nrow(state_across_datasets), " 个 discovery states。"),
  "",
  "## Stouffer 总体恢复证据",
  "",
  paste0("- 等权 Stouffer Meta Z > 1.96：", nrow(stouffer_significant), "/", nrow(stouffer_meta), " 个 discovery states。"),
  paste0("- 至少两个数据集可评估的正式 meta-analysis states：", sum(stouffer_meta$DatasetsEvaluable >= 2L), "/", nrow(stouffer_meta), "。"),
  paste0("- 仅单个数据集可评估、保持单队列证据标签：", sum(stouffer_meta$DatasetsEvaluable == 1L), "/", nrow(stouffer_meta), "。"),
  "",
  "BH 校正列作为补充敏感性结果保留；主判定仍为原文 Z > 1.96。"
)
writeLines(report_lines, file.path(output_dir, "SC_RNA_RECOVERY_SUMMARY.md"))

writeLines(c(
  "STATUS=PASS_SC_RNA_RECOVERY_SUMMARY",
  paste0("DATASETS=", nrow(dataset_summary)),
  paste0("PARENT_TYPE_JOBS=", nrow(parent_by_dataset)),
  paste0("STATE_DATASET_TESTS=", nrow(state_by_dataset)),
  paste0("STATE_DATASET_TESTS_WITH_ASSIGNMENT=", sum(state_by_dataset$HasFormalAssignment)),
  paste0("RECOVERY_EVALUABLE_STATE_DATASET_TESTS=", sum(state_by_dataset$RecoveryEvaluable)),
  paste0("SIGNIFICANT_STATE_DATASET_TESTS=", nrow(significant_states)),
  paste0("UNIQUE_DISCOVERY_STATES=", nrow(state_across_datasets)),
  paste0("SIGNIFICANT_IN_AT_LEAST_ONE_DATASET=", sum(state_across_datasets$SignificantInAtLeastOneDataset)),
  paste0("STOUFFER_META_STATES=", nrow(stouffer_meta)),
  paste0("STOUFFER_META_SIGNIFICANT=", nrow(stouffer_significant)),
  paste0("STOUFFER_MULTI_DATASET_STATES=", sum(stouffer_meta$DatasetsEvaluable >= 2L)),
  paste0("STOUFFER_SINGLE_DATASET_STATES=", sum(stouffer_meta$DatasetsEvaluable == 1L)),
  paste0("Z_THRESHOLD=", z_threshold)
), file.path(output_dir, "summary_status.txt"))

cat("SC_RNA_RECOVERY_SUMMARY_PASS", nrow(state_by_dataset), "state-dataset tests\n")
