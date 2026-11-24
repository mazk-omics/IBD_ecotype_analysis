#!/usr/bin/env Rscript
# 本脚本比较官方 Cell 2021 Carcinoma 资源与本研究的 state QC 后样本纳入比例。
# 只读取正式 initial/state assignment，不重新分配样本或修改模型。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入和输出
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 14d_audit_cell2021_qc_sample_retention.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) stop("Missing audit packages", call. = FALSE)
library(data.table)

result_dir <- file.path(output_root, "10_state_annotation/hstar_umap_cell2021_exact")
carcinoma_root <- file.path(project_root, "software/ecotyper/EcoTyper/Carcinoma/Carcinoma_Fractions")
carcinoma_rank_file <- file.path(carcinoma_root, "Analysis/rank_selection/rank_data.txt")
current_manifest_file <- file.path(result_dir, "hstar_umap_manifest.tsv")
if (!file.exists(carcinoma_rank_file) || !file.exists(current_manifest_file)) stop("QC assignment input is missing", call. = FALSE)

# Step2. 统计官方 Carcinoma 模型各细胞类型的正式纳入数
carcinoma_rank <- fread(carcinoma_rank_file, colClasses = "character", check.names = FALSE)
carcinoma_rows <- vector("list", nrow(carcinoma_rank))
for (cell_index in seq_len(nrow(carcinoma_rank))) {
  cell_type <- carcinoma_rank[[1L]][[cell_index]]
  rank <- trimws(carcinoma_rank[[2L]][[cell_index]])
  type_dir <- file.path(carcinoma_root, "Cell_States/discovery", cell_type, rank)
  initial_file <- file.path(type_dir, "initial_state_assignment.txt")
  final_file <- file.path(type_dir, "state_assignment.txt")
  if (!file.exists(initial_file) || !file.exists(final_file)) stop("Carcinoma assignment missing: ", cell_type, call. = FALSE)

  initial_n <- nrow(fread(initial_file, check.names = FALSE))
  assigned_n <- nrow(fread(final_file, check.names = FALSE))
  carcinoma_rows[[cell_index]] <- data.table(
    Source = "Cell2021_Carcinoma", CellType = cell_type, Rank = rank,
    InitialSamples = initial_n, AssignedSamples = assigned_n,
    UnassignedAfterStateQC = initial_n - assigned_n, AssignedRate = assigned_n / initial_n
  )
}
carcinoma_counts <- rbindlist(carcinoma_rows)
if (uniqueN(carcinoma_counts$InitialSamples) != 1L || unique(carcinoma_counts$InitialSamples) != 6475L) stop("Unexpected Cell 2021 discovery sample count", call. = FALSE)

# Step3. 统计本研究并保存可追溯比较表
current_manifest <- fread(current_manifest_file)
current_counts <- current_manifest[, .(
  Source = "IBD2025_SignatureD", CellType, Rank = as.character(Rank),
  InitialSamples = AllSamples, AssignedSamples = FormalAssigned,
  UnassignedAfterStateQC = AllSamples - FormalAssigned, AssignedRate = FormalAssigned / AllSamples
)]

counts <- rbindlist(list(carcinoma_counts, current_counts), use.names = TRUE)
summary <- counts[, .(
  CellTypes = .N,
  AssignmentOpportunities = sum(InitialSamples),
  AssignedSamples = sum(AssignedSamples),
  UnassignedAfterStateQC = sum(UnassignedAfterStateQC),
  AssignedRate = sum(AssignedSamples) / sum(InitialSamples),
  MinimumCellTypeRate = min(AssignedRate),
  MaximumCellTypeRate = max(AssignedRate)
), by = Source]

fwrite(counts, file.path(result_dir, "cell2021_vs_ibd_state_qc_assignment_counts.tsv"), sep = "\t", quote = FALSE)
fwrite(summary, file.path(result_dir, "cell2021_vs_ibd_state_qc_assignment_summary.tsv"), sep = "\t", quote = FALSE)
jsonlite::write_json(
  list(status = "PASS_CELL2021_VS_IBD_STATE_QC_ASSIGNMENT_AUDIT", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sources = nrow(summary), rows = nrow(counts)),
  file.path(result_dir, "cell2021_vs_ibd_state_qc_assignment_audit_complete.json"), pretty = TRUE, auto_unbox = TRUE
)

print(summary)
cat("PASS_CELL2021_VS_IBD_STATE_QC_ASSIGNMENT_AUDIT\n")
