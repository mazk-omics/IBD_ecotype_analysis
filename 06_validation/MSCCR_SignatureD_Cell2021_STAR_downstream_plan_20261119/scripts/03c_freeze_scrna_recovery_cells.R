#!/usr/bin/env Rscript
# 本脚本按官方固定种子冻结外部 scRNA recovery 的实际细胞清单。
# 每个 dataset × parent type 最多 2,500 cells；只保存 ID，不读取表达矩阵。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取路径和前序输入审计
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 03c_freeze_scrna_recovery_cells.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("arrow", "data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
input_audit_file <- file.path(audit_dir, "scrna_input_audit.json")
if (!file.exists(input_audit_file)) stop("Run 03_audit_scrna_inputs.R first", call. = FALSE)
input_audit <- jsonlite::read_json(input_audit_file, simplifyVector = TRUE)
if (!identical(input_audit$status, "PASS_SCRNA_METADATA_FREEZE_RAW_HASH_PENDING")) stop("scRNA metadata audit status changed", call. = FALSE)

# Step2. 按 dataset 和父类型固定抽样
datasets <- c("SCP1884", "GSE214695", "SCP259", "IBDverse_E_MTAB_16999")
max_cells <- 2500L
sampling_seed <- 1234L
selected_parts <- list()
sampling_parts <- list()

for (dataset in datasets) {
  pool_file <- file.path(eligible_pool_root, paste0(dataset, "_eligible_ibd_cells.parquet"))
  pool <- as.data.table(arrow::read_parquet(pool_file, as_data_frame = TRUE))
  required_cols <- c("CellID", "CellType", "DonorID", "BiologicalSampleID")
  if (!all(required_cols %in% names(pool))) stop("Eligible pool columns changed: ", dataset, call. = FALSE)
  if (anyNA(pool$CellID) || anyDuplicated(pool$CellID)) stop("Eligible pool cell identity changed: ", dataset, call. = FALSE)

  parent_types <- sort(unique(pool$CellType))
  for (cell_type in parent_types) {
    type_pool <- pool[CellType == cell_type, .(CellID, CellType, DonorID, BiologicalSampleID)]
    setorder(type_pool, CellID)
    set.seed(sampling_seed)
    selected_index <- if (nrow(type_pool) > max_cells) sample(seq_len(nrow(type_pool)), max_cells, replace = FALSE) else seq_len(nrow(type_pool))
    selected <- type_pool[selected_index]
    selected[, `:=`(Dataset = dataset, SelectionOrder = seq_len(.N), SamplingSeed = sampling_seed)]
    setcolorder(selected, c("Dataset", "CellType", "CellID", "DonorID", "BiologicalSampleID", "SelectionOrder", "SamplingSeed"))

    selected_parts[[paste(dataset, cell_type, sep = "::")]] <- selected
    sampling_parts[[paste(dataset, cell_type, sep = "::")]] <- data.frame(
      Dataset = dataset,
      CellType = cell_type,
      EligibleCells = nrow(type_pool),
      SelectedCells = nrow(selected),
      SamplingApplied = nrow(type_pool) > max_cells,
      MaxCells = max_cells,
      SamplingSeed = sampling_seed,
      SelectedCellIDSHA256 = digest::digest(paste(selected$CellID, collapse = "\n"), algo = "sha256", serialize = FALSE),
      stringsAsFactors = FALSE
    )
    cat("RECOVERY_CELLS_FROZEN", dataset, cell_type, nrow(selected), "/", nrow(type_pool), "\n")
  }
  rm(pool)
  invisible(gc(FALSE))
}

# Step3. 保存冻结细胞和抽样身份
selected_cells <- rbindlist(selected_parts)
sampling_manifest <- rbindlist(sampling_parts)
if (anyDuplicated(selected_cells[, .(Dataset, CellID)])) stop("Selected cell IDs are not unique within dataset", call. = FALSE)
if (any(sampling_manifest$SelectedCells > max_cells)) stop("Selected cell cap exceeded", call. = FALSE)

data.table::fwrite(selected_cells, file.path(audit_dir, "recovery_selected_cells.tsv.gz"), sep = "\t")
data.table::fwrite(sampling_manifest, file.path(audit_dir, "recovery_cell_sampling_manifest.tsv"), sep = "\t")
summary <- list(
  status = "PASS_SCRNA_RECOVERY_CELL_FREEZE",
  datasets = datasets,
  parent_type_groups = nrow(sampling_manifest),
  total_selected_cells = nrow(selected_cells),
  max_cells_per_dataset_parent = max_cells,
  sampling_seed = sampling_seed,
  selected_cells_sha256 = digest::digest(file.path(audit_dir, "recovery_selected_cells.tsv.gz"), algo = "sha256", file = TRUE, serialize = FALSE),
  model_changed = FALSE,
  automatic_recovery_started = FALSE
)
jsonlite::write_json(summary, file.path(audit_dir, "scrna_recovery_cell_freeze.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo_cell_freeze.txt"))
cat("SCRNA_RECOVERY_CELL_FREEZE_PASS", nrow(selected_cells), "cells in", nrow(sampling_manifest), "dataset-parent groups\n")
