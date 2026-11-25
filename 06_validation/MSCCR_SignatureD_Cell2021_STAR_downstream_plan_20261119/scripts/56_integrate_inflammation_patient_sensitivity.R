#!/usr/bin/env Rscript
# 本脚本把外部 scRNA inflamed/noninflamed exact 与供者感知结果并入注释证据。
# 该结果不替代 discovery Fig. S2G/S2H；只增加外部敏感性列，不改模型或等级。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取审阅表和表型结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 56_integrate_inflammation_patient_sensitivity.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

annotation_file <- file.path(RUN_ROOT, "output/42_state_annotation_review_v4/state_annotation_evidence_master_v4_reviewed.tsv")
phenotype_file <- file.path(RUN_ROOT, "output/16_G7_G7A_validation_summary/exact_vs_patient_aware_by_state_level.tsv.gz")
output_dir <- file.path(RUN_ROOT, "output/43_state_annotation_inflammation_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(annotation_file) || !file.exists(phenotype_file)) stop("Annotation or phenotype input is missing", call. = FALSE)
annotation <- fread(annotation_file)
phenotype <- fread(phenotype_file)
phenotype <- phenotype[Axis == "Inflammation"]

if (nrow(annotation) != 85L || anyDuplicated(annotation$CellState)) stop("Reviewed annotation is not 85 unique states", call. = FALSE)
if (!all(c("Inflamed", "Non-inflamed") %in% phenotype$Level)) stop("Inflammation levels are incomplete", call. = FALSE)

print(phenotype[, .N, by = .(Dataset, ExactStatus, PatientStatus)][order(Dataset, ExactStatus, PatientStatus)])

# Step2. 按 state 汇总 exact 与供者感知证据
exact_summary <- phenotype[, .(
  ExactDatasetsEvaluable = uniqueN(Dataset[ExactStatus == "EVALUABLE"]),
  ExactDatasetsPositive = uniqueN(Dataset[ExactStatus == "EVALUABLE" & ExactPositive %in% TRUE]),
  ExactInflamedDatasets = uniqueN(Dataset[ExactStatus == "EVALUABLE" & ExactPositive %in% TRUE & Level == "Inflamed"]),
  ExactNoninflamedDatasets = uniqueN(Dataset[ExactStatus == "EVALUABLE" & ExactPositive %in% TRUE & Level == "Non-inflamed"]),
  ExactPositiveDatasetLevels = paste(sort(unique(paste(Dataset[ExactStatus == "EVALUABLE" & ExactPositive %in% TRUE], Level[ExactStatus == "EVALUABLE" & ExactPositive %in% TRUE], sep = ":"))), collapse = ";")
), by = CellState]

patient_summary <- phenotype[, .(
  PatientDatasetsEvaluable = uniqueN(Dataset[PatientStatus == "EVALUABLE"]),
  PatientDatasetsPositive = uniqueN(Dataset[PatientStatus == "EVALUABLE" & PatientPositive %in% TRUE]),
  PatientInflamedDatasets = uniqueN(Dataset[PatientStatus == "EVALUABLE" & PatientPositive %in% TRUE & Level == "Inflamed"]),
  PatientNoninflamedDatasets = uniqueN(Dataset[PatientStatus == "EVALUABLE" & PatientPositive %in% TRUE & Level == "Non-inflamed"]),
  PatientPositiveDatasetLevels = paste(sort(unique(paste(Dataset[PatientStatus == "EVALUABLE" & PatientPositive %in% TRUE], Level[PatientStatus == "EVALUABLE" & PatientPositive %in% TRUE], sep = ":"))), collapse = ";")
), by = CellState]

state_sensitivity <- merge(exact_summary, patient_summary, by = "CellState", all = TRUE, sort = FALSE)
state_sensitivity <- merge(annotation[, .(CellState, CellType, State, Ecotype, EcotypeMember, FinalReviewedGrade, FinalReviewedLabel, InflammationDirection)], state_sensitivity, by = "CellState", all.x = TRUE, sort = FALSE)
state_sensitivity[, ExactDirectionMatchedDatasets := fifelse(
  InflammationDirection == "Inflamed", ExactInflamedDatasets,
  fifelse(InflammationDirection == "Non-inflamed", ExactNoninflamedDatasets, 0L)
)]
state_sensitivity[, PatientDirectionMatchedDatasets := fifelse(
  InflammationDirection == "Inflamed", PatientInflamedDatasets,
  fifelse(InflammationDirection == "Non-inflamed", PatientNoninflamedDatasets, 0L)
)]
state_sensitivity[, PatientAwareSupport := fifelse(
  PatientDirectionMatchedDatasets >= 1L,
  "DIRECTION_MATCHED_IN_AT_LEAST_ONE_DATASET",
  fifelse(PatientDatasetsPositive >= 1L, "PATIENT_POSITIVE_OPPOSITE_OR_DISCOVERY_NULL", "NO_PATIENT_AWARE_POSITIVE")
)]
setorder(state_sensitivity, -PatientDirectionMatchedDatasets, -ExactDirectionMatchedDatasets, CellType, State)

# Step3. 并入注释证据并保存审计表
annotation_with_sensitivity <- merge(annotation, state_sensitivity[, .(
  CellState, ExactDatasetsEvaluable, ExactDatasetsPositive, ExactInflamedDatasets, ExactNoninflamedDatasets,
  ExactPositiveDatasetLevels, PatientDatasetsEvaluable, PatientDatasetsPositive, PatientInflamedDatasets,
  PatientNoninflamedDatasets, PatientPositiveDatasetLevels, ExactDirectionMatchedDatasets,
  PatientDirectionMatchedDatasets, PatientAwareSupport
)], by = "CellState", all.x = TRUE, sort = FALSE)
annotation_with_sensitivity[, ModelDefinitionChanged := FALSE]
setorder(annotation_with_sensitivity, CellType, State)

summary_table <- annotation_with_sensitivity[, .(
  States = .N,
  DiscoveryInflammationAssociated = sum(InflammationDirection %in% c("Inflamed", "Non-inflamed")),
  ExactPositiveAtLeastOneDataset = sum(ExactDatasetsPositive >= 1L),
  ExactDirectionMatchedAtLeastOneDataset = sum(ExactDirectionMatchedDatasets >= 1L),
  PatientPositiveAtLeastOneDataset = sum(PatientDatasetsPositive >= 1L),
  PatientDirectionMatchedAtLeastOneDataset = sum(PatientDirectionMatchedDatasets >= 1L)
), by = FinalReviewedGrade]

fwrite(state_sensitivity, file.path(output_dir, "state_inflammation_exact_patient_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(annotation_with_sensitivity, file.path(output_dir, "state_annotation_evidence_master_v4_with_inflammation_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(summary_table, file.path(output_dir, "state_annotation_inflammation_sensitivity_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

completion <- list(
  status = "PASS_STATE_ANNOTATION_INFLAMMATION_SENSITIVITY",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  formal_states = nrow(annotation_with_sensitivity),
  patient_aware_direction_matched_states = sum(annotation_with_sensitivity$PatientDirectionMatchedDatasets >= 1L),
  labels_changed = 0L,
  model_definition_changed = FALSE
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_inflammation_sensitivity_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(summary_table)
print(state_sensitivity[PatientDatasetsPositive >= 1L, .(
  CellState, Ecotype, InflammationDirection, PatientPositiveDatasetLevels,
  PatientDirectionMatchedDatasets, FinalReviewedGrade, FinalReviewedLabel
)])
cat("PASS_STATE_ANNOTATION_INFLAMMATION_SENSITIVITY\n")
