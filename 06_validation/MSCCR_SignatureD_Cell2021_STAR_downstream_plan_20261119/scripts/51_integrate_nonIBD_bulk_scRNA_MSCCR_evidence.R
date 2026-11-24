#!/usr/bin/env Rscript
# 本脚本整合外部 bulk、四套 scRNA 和 MSCCR Control 的 nonIBD 证据；MSCCR 仅作为覆盖支持，不计入显著性验证。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入和正式 85-state 清单
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 51_integrate_nonIBD_bulk_scRNA_MSCCR_evidence.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

bulk_dir <- file.path(RUN_ROOT, "output/32_GSE57945_bulk_known_phenotype_Cell2021")
concordance_dir <- file.path(RUN_ROOT, "output/33_Cell2021_bulk_scRNA_nonIBD_concordance")
msccr_dir <- file.path(RUN_ROOT, "output/37_MSCCR_Control_bulk_recovery_summary")
output_dir <- file.path(RUN_ROOT, "output/38_integrated_nonIBD_validation_evidence")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

formal_states <- fread(file.path(bulk_dir, "formal_state_universe.tsv"), colClasses = "character")
if (nrow(formal_states) != 85L || anyDuplicated(formal_states$CellState)) stop("Formal 85-state universe changed", call. = FALSE)

# Step2. 提取 GSE57945 nonIBD 表型富集
bulk <- fread(file.path(bulk_dir, "GSE57945_bulk_known_phenotype_enrichment_Cell2021.tsv"))
bulk <- bulk[Level == "nonIBD", .(
  CellType,
  State,
  CellState,
  BulkStatus = Status,
  BulkObserved = ObservedAssignedSamples,
  BulkZ = Z,
  BulkP = P_one_sided,
  BulkEmpiricalP = EmpiricalP_one_sided,
  BulkSignificant = PositiveZgt1.96,
  BulkBHSignificantSensitivity = Significant_BH_sensitivity
)]
if (anyDuplicated(bulk$CellState)) stop("Duplicated GSE57945 nonIBD state rows", call. = FALSE)

# Step3. 汇总四套 scRNA 的 nonIBD 表型富集
scrna <- fread(file.path(concordance_dir, "Cell2021_bulk_scRNA_nonIBD_state_concordance.tsv"))
scrna <- scrna[, .(
  ScRNADataset,
  CellType,
  State,
  CellState,
  ScRNAStatus,
  ScRNAObserved,
  ScRNAZ,
  ScRNAP,
  ScRNASignificant
)]
if (anyDuplicated(scrna[, .(ScRNADataset, CellState)])) stop("Duplicated scRNA state rows", call. = FALSE)

scrna_summary <- scrna[, .(
  ScRNADatasetsEvaluable = sum(ScRNAStatus == "EVALUABLE"),
  ScRNADatasetsSignificant = sum(ScRNAStatus == "EVALUABLE" & ScRNASignificant, na.rm = TRUE),
  ScRNADatasetsPositiveZ = sum(ScRNAStatus == "EVALUABLE" & ScRNAZ > 0, na.rm = TRUE),
  ScRNAMaxZ = if (any(ScRNAStatus == "EVALUABLE")) max(ScRNAZ[ScRNAStatus == "EVALUABLE"], na.rm = TRUE) else NA_real_,
  ScRNAMedianZ = if (any(ScRNAStatus == "EVALUABLE")) median(ScRNAZ[ScRNAStatus == "EVALUABLE"], na.rm = TRUE) else NA_real_,
  ScRNASignificantDatasets = paste(sort(ScRNADataset[ScRNAStatus == "EVALUABLE" & ScRNASignificant]), collapse = ";")
), by = .(CellType, State, CellState)]

scrna_wide <- dcast(
  scrna,
  CellType + State + CellState ~ ScRNADataset,
  value.var = c("ScRNAZ", "ScRNAP", "ScRNASignificant", "ScRNAStatus")
)

# Step4. 汇总 MSCCR Control 三部位覆盖，不生成表型 p 值
msccr <- fread(file.path(msccr_dir, "MSCCR_Control_cell_state_assignment_counts.tsv"))
if (anyDuplicated(msccr[, .(Region3, CellType, State)])) stop("Duplicated MSCCR state coverage rows", call. = FALSE)
msccr[, CellState := paste(CellType, State, sep = "_")]

msccr_summary <- msccr[, .(
  MSCCRRegionsObserved = sum(AssignedSamples > 0L),
  MSCCRRegionsAtLeastFiveSamples = sum(AssignedSamples >= 5L),
  MSCCRMinimumAssignmentFraction = min(AssignmentFraction),
  MSCCRMedianAssignmentFraction = median(AssignmentFraction),
  MSCCRMaximumAssignmentFraction = max(AssignmentFraction),
  MSCCRTotalAssignedSamples = sum(AssignedSamples)
), by = .(CellType, State, CellState)]

msccr_wide <- dcast(
  msccr,
  CellType + State + CellState ~ Region3,
  value.var = c("AssignedSamples", "AssignmentFraction")
)

# Step5. 合并为 state-level 证据总表
evidence <- merge(formal_states, bulk, by = c("CellType", "State", "CellState"), all.x = TRUE)
evidence <- merge(evidence, scrna_summary, by = c("CellType", "State", "CellState"), all.x = TRUE)
evidence <- merge(evidence, scrna_wide, by = c("CellType", "State", "CellState"), all.x = TRUE)
evidence <- merge(evidence, msccr_summary, by = c("CellType", "State", "CellState"), all.x = TRUE)
evidence <- merge(evidence, msccr_wide, by = c("CellType", "State", "CellState"), all.x = TRUE)

evidence[is.na(ScRNADatasetsEvaluable), ScRNADatasetsEvaluable := 0L]
evidence[is.na(ScRNADatasetsSignificant), ScRNADatasetsSignificant := 0L]
evidence[, ExternalBulkValidated := BulkStatus == "EVALUABLE" & BulkSignificant]
evidence[, ExternalScRNAValidated := ScRNADatasetsSignificant >= 1L]
evidence[, ExternalCrossModalityValidated := ExternalBulkValidated & ExternalScRNAValidated]
evidence[, MSCCRAllRegionsObserved := MSCCRRegionsObserved == 3L]
evidence[, EvidenceClass := fifelse(
  ExternalCrossModalityValidated & MSCCRAllRegionsObserved,
  "external bulk + scRNA; MSCCR all-region coverage",
  fifelse(
    ExternalCrossModalityValidated,
    "external bulk + scRNA; partial MSCCR coverage",
    fifelse(
      ExternalBulkValidated & MSCCRAllRegionsObserved,
      "external bulk only; MSCCR all-region coverage",
      fifelse(
        ExternalScRNAValidated & MSCCRAllRegionsObserved,
        "external scRNA only; MSCCR all-region coverage",
        "MSCCR coverage or no external significance"
      )
    )
  )
)]

if (nrow(evidence) != 85L || anyDuplicated(evidence$CellState)) stop("Integrated evidence table is not one row per formal state", call. = FALSE)
setorder(evidence, CellType, State)

evidence_summary <- evidence[, .(
  FormalStates = .N,
  ExternalBulkValidated = sum(ExternalBulkValidated, na.rm = TRUE),
  ExternalScRNAValidated = sum(ExternalScRNAValidated, na.rm = TRUE),
  ExternalCrossModalityValidated = sum(ExternalCrossModalityValidated, na.rm = TRUE),
  MSCCRAllRegionsObserved = sum(MSCCRAllRegionsObserved, na.rm = TRUE),
  MSCCRAtLeastTwoRegionsObserved = sum(MSCCRRegionsObserved >= 2L, na.rm = TRUE)
)]

lineage_summary <- evidence[, .(
  FormalStates = .N,
  ExternalBulkValidated = sum(ExternalBulkValidated, na.rm = TRUE),
  ExternalScRNAValidated = sum(ExternalScRNAValidated, na.rm = TRUE),
  ExternalCrossModalityValidated = sum(ExternalCrossModalityValidated, na.rm = TRUE),
  MSCCRAllRegionsObserved = sum(MSCCRAllRegionsObserved, na.rm = TRUE)
), by = CellType]

# Step6. 保存表格和方法边界
fwrite(evidence, file.path(output_dir, "integrated_nonIBD_state_validation_evidence.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(evidence[ExternalCrossModalityValidated == TRUE], file.path(output_dir, "external_cross_modality_validated_states.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(evidence_summary, file.path(output_dir, "integrated_nonIBD_validation_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(lineage_summary, file.path(output_dir, "integrated_nonIBD_validation_by_cell_type.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_INTEGRATED_NONIBD_VALIDATION_EVIDENCE",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  formal_states = nrow(evidence),
  external_bulk_dataset = "GSE57945",
  external_scrna_datasets = sort(unique(scrna$ScRNADataset)),
  msccr_role = "descriptive internal coverage only; excluded from significance and cross-modality validation calls",
  primary_external_call = "GSE57945 nonIBD Z>1.96 and at least one scRNA nonIBD dataset Z>1.96",
  caveat = "scRNA datasets are summarized separately; no pooled Stouffer call is used as primary evidence"
)
jsonlite::write_json(completion, file.path(output_dir, "integrated_nonIBD_validation_evidence_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(evidence_summary)
print(evidence[ExternalCrossModalityValidated == TRUE, .(CellState, BulkZ, ScRNADatasetsSignificant, ScRNASignificantDatasets, MSCCRRegionsObserved)])
cat("PASS_INTEGRATED_NONIBD_VALIDATION_EVIDENCE\n")
