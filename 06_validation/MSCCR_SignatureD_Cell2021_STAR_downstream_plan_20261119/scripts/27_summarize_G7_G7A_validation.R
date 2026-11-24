#!/usr/bin/env Rscript
# 本脚本汇总 G7/G7A、marker LODO 和 recovery 结果，生成可直接审阅的统计表与中文报告。
# 原文轨道、患者感知扩展和 representation 保持分表，不用单一通过率掩盖分母差异。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读入完成结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 27_summarize_G7_G7A_validation.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

output_dir <- file.path(RUN_ROOT, "output/16_G7_G7A_validation_summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_files <- c(
  recovery = file.path(RUN_ROOT, "output/08_scRNA_recovery_summary/scRNA_recovery_dataset_summary.tsv"),
  recovery_meta = file.path(RUN_ROOT, "output/08_scRNA_recovery_summary/scRNA_recovery_stouffer_meta.tsv"),
  marker = file.path(RUN_ROOT, "output/09_marker_recovery/ranking_and_LODO/marker_LODO_summary.tsv"),
  marker_state = file.path(RUN_ROOT, "output/09_marker_recovery/ranking_and_LODO/marker_LODO_per_state.tsv"),
  exact = file.path(RUN_ROOT, "output/11_known_phenotype_enrichment_exact/state_known_phenotype_enrichment_exact.tsv"),
  patient = file.path(RUN_ROOT, "output/12_known_phenotype_enrichment_patient_aware/state_known_phenotype_enrichment_patient_aware.tsv"),
  representation = file.path(RUN_ROOT, "output/13_cell_state_representation/sample_dominant_state_representation_fisher.tsv"),
  nonibd = file.path(RUN_ROOT, "output/13_cell_state_representation/IBD_nonIBD_representation_gate.tsv")
)
if (any(!file.exists(input_files))) stop("Missing one or more completed validation inputs", call. = FALSE)

recovery <- fread(input_files[["recovery"]])
recovery_meta <- fread(input_files[["recovery_meta"]])
marker <- fread(input_files[["marker"]])
marker_state <- fread(input_files[["marker_state"]])
exact <- fread(input_files[["exact"]])
patient <- fread(input_files[["patient"]])
representation <- fread(input_files[["representation"]])
nonibd <- fread(input_files[["nonibd"]], colClasses = "character")

# Step2. 汇总每个正式轨道
exact_summary <- exact[, .(
  Rows = .N, EvaluableRows = sum(Status == "EVALUABLE"), PositiveRows = sum(PositiveZgt1.96),
  UniqueStatesEvaluable = uniqueN(CellState[Status == "EVALUABLE"]), UniqueStatesPositive = uniqueN(CellState[PositiveZgt1.96 == TRUE])
), by = .(Dataset, Axis)]
exact_summary[, PositiveFractionEvaluableRows := PositiveRows / EvaluableRows]

patient_summary <- patient[, .(
  Rows = .N, EvaluableRows = sum(Status == "EVALUABLE"), PositiveRows = sum(PositiveZgt1.96),
  UniqueStatesEvaluable = uniqueN(CellState[Status == "EVALUABLE"]), UniqueStatesPositive = uniqueN(CellState[PositiveZgt1.96 == TRUE])
), by = .(Dataset, Axis)]
patient_summary[, PositiveFractionEvaluableRows := PositiveRows / EvaluableRows]

comparison <- merge(
  exact[, .(Dataset, CellType, State, CellState, Axis, Level, ExactZ = Z, ExactPositive = PositiveZgt1.96, ExactStatus = Status)],
  patient[, .(Dataset, CellType, State, CellState, Axis, Level, PatientZ = Z, PatientPositive = PositiveZgt1.96, PatientStatus = Status)],
  by = c("Dataset", "CellType", "State", "CellState", "Axis", "Level"), all = TRUE, sort = FALSE
)
comparison[, ExactPositiveFlag := !is.na(ExactPositive) & ExactPositive == TRUE]
comparison[, PatientPositiveFlag := !is.na(PatientPositive) & PatientPositive == TRUE]
comparison_summary <- comparison[, .(
  UnionRows = .N, BothEvaluable = sum(ExactStatus == "EVALUABLE" & PatientStatus == "EVALUABLE", na.rm = TRUE),
  ExactPositive = sum(ExactPositiveFlag), PatientPositive = sum(PatientPositiveFlag),
  PositiveBoth = sum(ExactPositiveFlag & PatientPositiveFlag),
  PositiveExactOnly = sum(ExactPositiveFlag & !PatientPositiveFlag),
  PositivePatientOnly = sum(!ExactPositiveFlag & PatientPositiveFlag)
), by = .(Dataset, Axis)]

representation_summary <- representation[, .(
  Rows = .N, EvaluableRows = sum(Status == "EVALUABLE"), SignificantRows = sum(SignificantEnrichment),
  UniqueStatesSignificant = uniqueN(CellState[SignificantEnrichment == TRUE]), Levels = uniqueN(Level)
), by = .(Scope, DatasetScope, Axis)]
marker_state_summary <- marker_state[, .(
  FormalStates = .N, StatesEvaluableAtLeastOneFold = sum(EvaluableFolds > 0),
  StatesConcordantAtLeastOneFold = sum(ConcordantFolds > 0), MedianStateConcordance = median(Concordance, na.rm = TRUE)
), by = Variant]

setorder(exact_summary, Dataset, Axis)
setorder(patient_summary, Dataset, Axis)
setorder(comparison_summary, Dataset, Axis)
setorder(representation_summary, Scope, DatasetScope, Axis)
fwrite(exact_summary, file.path(output_dir, "known_phenotype_exact_summary.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(patient_summary, file.path(output_dir, "known_phenotype_patient_aware_summary.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(comparison_summary, file.path(output_dir, "exact_vs_patient_aware_summary.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(comparison, file.path(output_dir, "exact_vs_patient_aware_by_state_level.tsv.gz"), sep = "\t", quote = FALSE, na = "")
fwrite(representation_summary, file.path(output_dir, "cell_state_representation_summary.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(marker_state_summary, file.path(output_dir, "marker_LODO_state_summary.tsv"), sep = "\t", quote = FALSE, na = "")

# Step3. 生成中文审阅报告
report <- c(
  "# Signature-D cell-state validation：G7/G7A 汇总",
  "",
  "## 1. State recovery",
  paste0("- 四个 scRNA-seq 队列共检验 ", sum(recovery$FormalStatesEvaluated), " 个 state×dataset 条目；", sum(recovery$RecoveryEvaluableStates), " 个可评价，", sum(recovery$SignificantStates_Zgt1.96), " 个达到 Z>1.96。"),
  paste0("- 等权 Stouffer Meta Z>1.96：", sum(recovery_meta$MetaSignificant_Zgt1.96, na.rm = TRUE), "/", nrow(recovery_meta), " 个 discovery states。"),
  "",
  "## 2. Marker LODO",
  paste0("- Disease 版：", marker[Variant == "Disease", sum(ConcordantStates)], "/", marker[Variant == "Disease", sum(EvaluableStates)], " 个可评价 state×held-out 条目预测一致。"),
  paste0("- Region3 敏感性版：", marker[Variant == "Region3", sum(ConcordantStates)], "/", marker[Variant == "Region3", sum(EvaluableStates)], " 个可评价 state×held-out 条目预测一致。"),
  "- 不可评价条目保留在热图中，不作为错误或成功计入分母。",
  "",
  "## 3. Known phenotype enrichment",
  paste0("- Cell 2021 exact 逐细胞置换：", sum(exact_summary$PositiveRows), "/", sum(exact_summary$EvaluableRows), " 个可评价 state×phenotype 格子 Z>1.96。"),
  paste0("- 患者感知 IBD 扩展：", sum(patient_summary$PositiveRows), "/", sum(patient_summary$EvaluableRows), " 个可评价格子 Z>1.96。"),
  paste0("- 两轨道共同阳性：", sum(comparison_summary$PositiveBoth), " 个格子；exact-only ", sum(comparison_summary$PositiveExactOnly), "；patient-aware-only ", sum(comparison_summary$PositivePatientOnly), "。"),
  "- 两轨道统计单位不同，必须并列解释，不能把患者感知结果替换成原文复现结果。",
  "",
  "## 4. Cell-state representation",
  paste0("- dominant-state Fisher 共得到 ", sum(representation_summary$SignificantRows), " 个 BH FDR<0.05 的 phenotype/dataset representation 格子。"),
  paste0("- IBD/nonIBD 当前状态：", nonibd$Status[[1L]], "；", nonibd$Reason[[1L]], "。"),
  "",
  "## 5. 解释边界",
  "- 当前 G7/G7A 评价的是冻结 state 在四个 IBD 单细胞队列中的恢复、marker 特异性和表型分布，不重新训练 state。",
  "- pooled representation 可能同时包含 cohort、部位和疾病构成差异；逐数据集结果及患者感知扩展是必要的并列证据。",
  "- nonIBD 比较必须等待独立健康/非 IBD 数据进入 G11B，不能从本轮 IBD-only 队列推断。"
)
writeLines(report, file.path(output_dir, "G7_G7A_VALIDATION_SUMMARY.md"))

completion <- list(
  status = "PASS_G7_G7A_VALIDATION_SUMMARY", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  recovery_meta_significant_states = sum(recovery_meta$MetaSignificant_Zgt1.96, na.rm = TRUE),
  exact_positive_rows = sum(exact_summary$PositiveRows), patient_aware_positive_rows = sum(patient_summary$PositiveRows),
  representation_significant_rows = sum(representation_summary$SignificantRows), ibd_nonibd_status = nonibd$Status[[1L]]
)
jsonlite::write_json(completion, file.path(output_dir, "G7_G7A_validation_summary_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_G7_G7A_VALIDATION_SUMMARY\n")
