#!/usr/bin/env Rscript
# 本脚本冻结四套外部 scRNA eligible pools、父类型映射和供者/样本表型。
# 只写样本级清单和审计表；不复制细胞矩阵，不运行 state recovery。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取路径、正式 ecotype 成员和外部资格表
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 03_audit_scrna_inputs.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("arrow", "data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

ecotype_audit_file <- file.path(output_root, "02_ecotype_audit/ecotype_discovery_audit.json")
if (!file.exists(ecotype_audit_file)) stop("Run 02_audit_ecotype_discovery.R first", call. = FALSE)
ecotype_audit <- jsonlite::read_json(ecotype_audit_file, simplifyVector = TRUE)
if (!identical(ecotype_audit$status, "PASS_ECOTYPE_AUDIT")) stop("Ecotype audit did not pass", call. = FALSE)
if (!file.exists(qualification_file)) stop("Missing frozen donor/sample qualification table", call. = FALSE)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("SCP1884", "GSE214695", "SCP259", "IBDverse_E_MTAB_16999")
expected_cells <- c(SCP1884 = 504050L, GSE214695 = 26182L, SCP259 = 254333L, IBDverse_E_MTAB_16999 = 532643L)
expected_donors <- c(SCP1884 = 46L, GSE214695 = 12L, SCP259 = 18L, IBDverse_E_MTAB_16999 = 111L)
expected_samples <- c(SCP1884 = 91L, GSE214695 = 12L, SCP259 = 85L, IBDverse_E_MTAB_16999 = 111L)
expected_types <- c(SCP1884 = 20L, GSE214695 = 18L, SCP259 = 18L, IBDverse_E_MTAB_16999 = 16L)
joint20 <- c(
  "B", "Plasma", "CD4_T", "CD8_T", "NK", "ILC", "Monocyte_Macrophage", "DC", "Mast",
  "Fibroblast", "Pericyte", "Endothelial", "Ileal_absorptive", "Colonic_absorptive", "Goblet",
  "Stem_TA_progenitor", "Paneth", "Enteroendocrine", "Glial", "Tuft"
)

qualification <- data.table::fread(qualification_file, colClasses = "character")
required_qualification <- c("dataset", "donor_id", "disease", "inflammation", "eligibility")
if (!all(required_qualification %in% names(qualification))) stop("Qualification columns changed", call. = FALSE)
ecotypes <- data.table::fread(file.path(ecotype_discovery_root, "ecotypes.txt"), colClasses = "character", check.names = FALSE)
if (nrow(ecotypes) != 79L || length(unique(ecotypes$Ecotype)) != 8L) stop("Frozen ecotype membership changed", call. = FALSE)

# Step2. 逐队列核对细胞、供者、样本和表型映射
pool_manifest <- list()
dataset_summary <- list()
sample_parts <- list()
type_parts <- list()
state_parts <- list()

for (dataset in datasets) {
  current_dataset <- dataset
  pool_file <- file.path(eligible_pool_root, paste0(dataset, "_eligible_ibd_cells.parquet"))
  if (!file.exists(pool_file)) stop("Missing eligible pool: ", dataset, call. = FALSE)
  pool <- as.data.table(arrow::read_parquet(pool_file, as_data_frame = TRUE))
  required_pool <- c("Dataset", "CellID", "CellType", "DonorID", "BiologicalSampleID", "Region", "Inflammation")
  if (!all(required_pool %in% names(pool))) stop("Eligible pool columns changed: ", dataset, call. = FALSE)
  pool[, DonorID := as.character(DonorID)]
  pool[, BiologicalSampleID := as.character(BiologicalSampleID)]
  if (nrow(pool) != expected_cells[[dataset]] || uniqueN(pool$DonorID) != expected_donors[[dataset]] || uniqueN(pool$BiologicalSampleID) != expected_samples[[dataset]]) stop("Eligible pool dimensions changed: ", dataset, call. = FALSE)
  if (anyNA(pool$CellID) || anyDuplicated(pool$CellID) || !all(pool$CellType %in% joint20) || uniqueN(pool$CellType) != expected_types[[dataset]]) stop("Eligible pool identity changed: ", dataset, call. = FALSE)

  donor_disease <- unique(qualification[dataset == current_dataset & eligibility == "INCLUDE_IBD", .(DonorID = donor_id, Disease = disease)])
  if (donor_disease[, .N, by = DonorID][N > 1L, .N] > 0L) stop("Multiple diseases per donor: ", dataset, call. = FALSE)
  pool <- merge(pool, donor_disease, by = "DonorID", all.x = TRUE, sort = FALSE)
  if (anyNA(pool$Disease) || !all(pool$Disease %in% c("CD", "UC"))) stop("Disease mapping incomplete: ", dataset, call. = FALSE)

  pool[, InflammationFrozen := as.character(Inflammation)]
  if (any(is.na(pool$InflammationFrozen) | pool$InflammationFrozen == "")) {
    donor_inflammation <- unique(qualification[dataset == current_dataset & eligibility == "INCLUDE_IBD", .(DonorID = donor_id, QualifiedInflammation = inflammation)])
    if (donor_inflammation[, .N, by = DonorID][N > 1L, .N] > 0L) stop("Missing pool inflammation cannot be filled at donor level: ", dataset, call. = FALSE)
    pool <- merge(pool, donor_inflammation, by = "DonorID", all.x = TRUE, sort = FALSE)
    pool[is.na(InflammationFrozen) | InflammationFrozen == "", InflammationFrozen := QualifiedInflammation]
    pool[, QualifiedInflammation := NULL]
  }
  inflammation_lower <- tolower(pool$InflammationFrozen)
  pool[, InflammationFrozen := ifelse(grepl("non|uninf", inflammation_lower), "Non-inflamed", ifelse(grepl("inflam", inflammation_lower), "Inflamed", NA_character_))]

  region_lower <- tolower(as.character(pool$Region))
  pool[, Region3 := ifelse(region_lower %in% c("co", "colon"), "Colon", ifelse(region_lower == "rectum", "Rectum", ifelse(region_lower %in% c("ti", "sb", "terminal ileum"), "Ileum", NA_character_)))]
  if (anyNA(pool$InflammationFrozen) || anyNA(pool$Region3)) stop("Inflammation/Region3 mapping incomplete: ", dataset, call. = FALSE)

  sample_metadata <- unique(pool[, .(Dataset = dataset, DonorID, BiologicalSampleID, Disease, Inflammation = InflammationFrozen, Region3)])
  sample_conflict <- sample_metadata[, .N, by = .(Dataset, BiologicalSampleID)][N > 1L]
  if (nrow(sample_conflict) > 0L) stop("Biological sample maps to multiple phenotypes: ", dataset, call. = FALSE)
  sample_parts[[dataset]] <- sample_metadata

  type_counts <- pool[, .(EligibleCells = .N, Donors = uniqueN(DonorID), BiologicalSamples = uniqueN(BiologicalSampleID)), by = CellType]
  type_grid <- data.table::data.table(Dataset = dataset, CellType = joint20)
  type_grid <- merge(type_grid, type_counts, by = "CellType", all.x = TRUE, sort = FALSE)
  type_grid[is.na(EligibleCells), `:=`(EligibleCells = 0L, Donors = 0L, BiologicalSamples = 0L)]
  type_grid[, Status := ifelse(EligibleCells > 0L, "SUPPORTED", "SOURCE_TYPE_ABSENT")]
  type_parts[[dataset]] <- type_grid

  state_grid <- merge(ecotypes[, .(ID, CellType, Ecotype)], type_grid[, .(CellType, EligibleCells, Status)], by = "CellType", all.x = TRUE, sort = FALSE)
  state_grid[, Dataset := dataset]
  state_grid[, Evaluable := Status == "SUPPORTED"]
  data.table::setcolorder(state_grid, c("Dataset", "ID", "CellType", "Ecotype", "EligibleCells", "Status", "Evaluable"))
  state_parts[[dataset]] <- state_grid

  pool_manifest[[dataset]] <- data.frame(Dataset = dataset, Path = pool_file, Bytes = file.info(pool_file)$size, SHA256 = digest::digest(pool_file, algo = "sha256", file = TRUE, serialize = FALSE))
  dataset_summary[[dataset]] <- data.frame(
    Dataset = dataset, EligibleCells = nrow(pool), Donors = uniqueN(pool$DonorID), BiologicalSamples = uniqueN(pool$BiologicalSampleID),
    ParentTypes = uniqueN(pool$CellType), EvaluableStates = sum(state_grid$Evaluable), NonEvaluableStates = sum(!state_grid$Evaluable),
    Diseases = paste(sort(unique(pool$Disease)), collapse = ";"), Inflammation = paste(sort(unique(pool$InflammationFrozen)), collapse = ";"), Region3 = paste(sort(unique(pool$Region3)), collapse = ";")
  )
  cat("SCRNA_INPUT_AUDITED", dataset, nrow(pool), "cells", uniqueN(pool$DonorID), "donors", sum(state_grid$Evaluable), "/79 states evaluable\n")
  rm(pool)
  invisible(gc(FALSE))
}

# Step3. 记录表达源路径和仍待闭环的完整原始文件哈希
raw_sources <- data.frame(
  Dataset = datasets,
  Path = c(
    "/home/mazekai/enteric_glia/data_scRNA/SCP1884/expression",
    file.path(project_root, "06_validation/GSE214695_scRNA_recovery_20261009/source"),
    file.path(project_root, "01_raw_data/scRNA/SCP259/raw/SCP259/expression"),
    "/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad"
  ),
  stringsAsFactors = FALSE
)
raw_sources$Exists <- file.exists(raw_sources$Path) | dir.exists(raw_sources$Path)
raw_sources$PriorStructureAudit <- raw_sources$Exists
raw_sources$FullRawSHA256 <- NA_character_
raw_sources$Status <- ifelse(raw_sources$Exists, "PATH_AND_STRUCTURE_VERIFIED_FULL_HASH_PENDING", "MISSING")
if (any(!raw_sources$Exists)) stop("External scRNA expression source missing", call. = FALSE)

# Step4. 保存冻结清单和数据角色声明
pool_manifest_all <- data.table::rbindlist(pool_manifest)
qualification_manifest <- data.frame(Dataset = "ALL", Path = qualification_file, Bytes = file.info(qualification_file)$size, SHA256 = digest::digest(qualification_file, algo = "sha256", file = TRUE, serialize = FALSE))
input_manifest <- data.table::rbindlist(list(pool_manifest_all, qualification_manifest), fill = TRUE)
dataset_summary_all <- data.table::rbindlist(dataset_summary)
sample_metadata_all <- data.table::rbindlist(sample_parts)
type_coverage_all <- data.table::rbindlist(type_parts)
state_evaluability_all <- data.table::rbindlist(state_parts)

data.table::fwrite(input_manifest, file.path(audit_dir, "scrna_metadata_input_manifest.tsv"), sep = "\t")
data.table::fwrite(dataset_summary_all, file.path(audit_dir, "dataset_summary.tsv"), sep = "\t")
data.table::fwrite(sample_metadata_all, file.path(audit_dir, "donor_sample_phenotype_freeze.tsv"), sep = "\t")
data.table::fwrite(type_coverage_all, file.path(audit_dir, "dataset_parent_type_coverage.tsv"), sep = "\t")
data.table::fwrite(state_evaluability_all, file.path(audit_dir, "dataset_state_evaluability.tsv"), sep = "\t")
data.table::fwrite(raw_sources, file.path(audit_dir, "raw_expression_source_status.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo.txt"))

role_statement <- c(
  "# 外部 scRNA 数据角色与泄漏边界",
  "",
  "四套 eligible pools 曾用于 Signature-D 的 expression-purification/signature validation，但未参与本轮 state NMF 或 ecotype discovery。",
  "因此后续 recovery 是未参与 state/ecotype discovery 的正交终点验证，不表述为完全未接触的独立验证。",
  "父类型、供者、样本、Disease、Inflammation 和 Region3 映射在 recovery 前冻结；不得根据 recovery Z、LODO 或 ecotype detection 结果修改。",
  "缺失父类型只标记对应 state 为不可评价，不以 0 补齐，也不把细胞重新分配到相近父类型。"
)
writeLines(role_statement, file.path(audit_dir, "LEAKAGE_AND_ROLE_STATEMENT.md"))

summary <- list(
  status = "PASS_SCRNA_METADATA_FREEZE_RAW_HASH_PENDING",
  datasets = datasets,
  total_eligible_cells = sum(dataset_summary_all$EligibleCells),
  total_donors_by_dataset = as.list(setNames(dataset_summary_all$Donors, dataset_summary_all$Dataset)),
  evaluable_states_by_dataset = as.list(setNames(dataset_summary_all$EvaluableStates, dataset_summary_all$Dataset)),
  ecotype_member_states = nrow(ecotypes),
  metadata_hashes_complete = TRUE,
  raw_expression_paths_verified = TRUE,
  raw_expression_full_sha256_complete = FALSE,
  leakage_statement_complete = TRUE,
  model_changed = FALSE,
  automatic_recovery_started = FALSE,
  release_decision = "G4_PARTIAL_FULL_RAW_SHA256_PENDING"
)
jsonlite::write_json(summary, file.path(audit_dir, "scrna_input_audit.json"), pretty = TRUE, auto_unbox = TRUE)
cat("SCRNA_INPUT_AUDIT_PASS metadata frozen; full raw expression SHA256 pending; recovery not started\n")
