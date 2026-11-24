#!/usr/bin/env Rscript
# 本脚本比较当前 IBD 模型与 Cell 2021 官方 Carcinoma 模型的 H* state 分离度。
# 只读取冻结的 state abundance 和正式 assignment，不重新运行 UMAP 或修改模型。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入和输出
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 14e_audit_hstar_compactness_vs_cell2021.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

cell2021_root <- file.path(ecotyper_root, "EcoTyper/Carcinoma/Carcinoma_Fractions/Cell_States/discovery")
ibd_root <- state_discovery_root
result_dir <- file.path(output_root, "10_state_annotation/hstar_compactness_audit")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

model_roots <- c(IBD2025 = ibd_root, Cell2021_Carcinoma = cell2021_root)
model_rank_files <- c(
  IBD2025 = rank_data_file,
  Cell2021_Carcinoma = file.path(dirname(dirname(cell2021_root)), "Analysis/rank_selection/rank_data.txt")
)
if (any(!dir.exists(model_roots))) stop("H* compactness model root is missing", call. = FALSE)
if (any(!file.exists(model_rank_files))) stop("H* compactness rank file is missing", call. = FALSE)

# Step2. 逐模型和细胞类型计算 H* 清晰度
cell_results <- list()
sample_results <- list()
result_index <- 1L

for (dataset_name in names(model_roots)) {
  rank_data <- fread(model_rank_files[[dataset_name]], colClasses = "character", check.names = FALSE)
  rank_lookup <- setNames(trimws(rank_data[[2L]]), trimws(rank_data[[1L]]))
  cell_dirs <- list.dirs(model_roots[[dataset_name]], full.names = TRUE, recursive = FALSE)
  cell_dirs <- cell_dirs[order(basename(cell_dirs))]

  for (cell_dir in cell_dirs) {
    selected_rank <- rank_lookup[[basename(cell_dir)]]
    if (is.null(selected_rank) || is.na(selected_rank)) stop("Selected rank is missing: ", cell_dir, call. = FALSE)
    rank_dir <- file.path(cell_dir, selected_rank)
    if (!dir.exists(rank_dir)) stop("Selected rank directory is missing: ", rank_dir, call. = FALSE)

    abundance_file <- file.path(rank_dir, "state_abundances.txt")
    assignment_file <- file.path(rank_dir, "state_assignment.txt")
    if (!file.exists(abundance_file) || !file.exists(assignment_file)) stop("Missing frozen H* input: ", cell_dir, call. = FALSE)

    abundance <- fread(abundance_file, check.names = FALSE)
    state_ids <- abundance[[1L]]
    abundance <- as.matrix(abundance[, -1L, with = FALSE])
    rownames(abundance) <- state_ids
    storage.mode(abundance) <- "double"

    assignment <- fread(assignment_file, header = FALSE, skip = 1L, colClasses = "character")
    if (ncol(assignment) != 4L) stop("Unexpected assignment layout: ", assignment_file, call. = FALSE)
    setnames(assignment, c("RowIndex", "SampleID", "State", "InitialState"))
    if (anyDuplicated(assignment$SampleID) || any(!assignment$SampleID %in% colnames(abundance))) stop("Assignment identity failed: ", cell_dir, call. = FALSE)

    hstar <- t(abundance)[assignment$SampleID, , drop = FALSE]
    if (any(!is.finite(hstar)) || any(hstar < 0) || any(rowSums(hstar) <= 0)) stop("Invalid assigned H*: ", cell_dir, call. = FALSE)
    formal_argmax <- colnames(hstar)[max.col(hstar, ties.method = "first")]
    if (any(formal_argmax != assignment$State)) stop("Formal assignment and H* argmax differ: ", cell_dir, call. = FALSE)

    hstar_fraction <- hstar / rowSums(hstar)
    ordered_fraction <- t(apply(hstar_fraction, 1L, sort, decreasing = TRUE))
    top_share <- ordered_fraction[, 1L]
    second_share <- ordered_fraction[, 2L]
    top_margin <- top_share - second_share
    normalized_entropy <- -rowSums(ifelse(hstar_fraction > 0, hstar_fraction * log(hstar_fraction), 0)) / log(ncol(hstar_fraction))

    centroids <- rowsum(hstar, assignment$State, reorder = FALSE) / as.numeric(table(factor(assignment$State, levels = unique(assignment$State))))
    centroid_distance <- matrix(NA_real_, nrow = nrow(hstar), ncol = nrow(centroids), dimnames = list(rownames(hstar), rownames(centroids)))
    for (state_index in seq_len(nrow(centroids))) {
      centroid_difference <- sweep(hstar, 2L, centroids[state_index, ], "-")
      centroid_distance[, state_index] <- sqrt(rowSums(centroid_difference^2))
    }

    own_index <- match(assignment$State, colnames(centroid_distance))
    own_distance <- centroid_distance[cbind(seq_len(nrow(centroid_distance)), own_index)]
    centroid_distance[cbind(seq_len(nrow(centroid_distance)), own_index)] <- Inf
    nearest_other_distance <- apply(centroid_distance, 1L, min)
    centroid_margin <- (nearest_other_distance - own_distance) / (nearest_other_distance + own_distance + .Machine$double.eps)
    own_centroid_closest <- own_distance < nearest_other_distance

    sample_metrics <- data.table(
      Dataset = dataset_name, CellType = basename(cell_dir), SampleID = assignment$SampleID, State = assignment$State,
      TopStateShare = top_share, TopStateMargin = top_margin, NormalizedEntropy = normalized_entropy,
      OwnCentroidDistance = own_distance, NearestOtherCentroidDistance = nearest_other_distance,
      CentroidMargin = centroid_margin, OwnCentroidClosest = own_centroid_closest
    )
    sample_results[[result_index]] <- sample_metrics

    state_counts_current <- sample_metrics[, .N, by = State]
    cell_results[[result_index]] <- sample_metrics[, .(
      States = uniqueN(State), AssignedSamples = .N,
      MedianTopStateShare = median(TopStateShare), MedianTopStateMargin = median(TopStateMargin),
      MedianNormalizedEntropy = median(NormalizedEntropy), MedianCentroidMargin = median(CentroidMargin),
      OwnCentroidClosestFraction = mean(OwnCentroidClosest)
    ), by = .(Dataset, CellType)]
    cell_results[[result_index]][, MedianSamplesPerState := median(state_counts_current$N)]
    setcolorder(cell_results[[result_index]], c("Dataset", "CellType", "States", "AssignedSamples", "MedianSamplesPerState", "MedianTopStateShare", "MedianTopStateMargin", "MedianNormalizedEntropy", "MedianCentroidMargin", "OwnCentroidClosestFraction"))

    cat("HSTAR_COMPACTNESS", dataset_name, basename(cell_dir), nrow(hstar), "samples", ncol(hstar), "states\n")
    result_index <- result_index + 1L
  }
}

# Step3. 汇总并保存可核对结果
sample_metrics <- rbindlist(sample_results, use.names = TRUE)
cell_metrics <- rbindlist(cell_results, use.names = TRUE)
state_counts <- sample_metrics[, .N, by = .(Dataset, CellType, State)]
dataset_metrics <- sample_metrics[, .(
  CellTypes = uniqueN(CellType), States = uniqueN(paste(CellType, State, sep = "::")), AssignedSamples = .N,
  MedianTopStateShare = median(TopStateShare), MedianTopStateMargin = median(TopStateMargin),
  MedianNormalizedEntropy = median(NormalizedEntropy), MedianCentroidMargin = median(CentroidMargin),
  OwnCentroidClosestFraction = mean(OwnCentroidClosest)
), by = Dataset]
dataset_metrics <- merge(dataset_metrics, state_counts[, .(MedianSamplesPerState = median(N)), by = Dataset], by = "Dataset", sort = FALSE)
setcolorder(dataset_metrics, c("Dataset", "CellTypes", "States", "AssignedSamples", "MedianSamplesPerState", "MedianTopStateShare", "MedianTopStateMargin", "MedianNormalizedEntropy", "MedianCentroidMargin", "OwnCentroidClosestFraction"))

setorder(cell_metrics, Dataset, CellType)
fwrite(sample_metrics, file.path(result_dir, "hstar_compactness_sample_metrics.tsv.gz"), sep = "\t", quote = FALSE)
fwrite(cell_metrics, file.path(result_dir, "hstar_compactness_cell_type_metrics.tsv"), sep = "\t", quote = FALSE)
fwrite(dataset_metrics, file.path(result_dir, "hstar_compactness_dataset_summary.tsv"), sep = "\t", quote = FALSE)
print(dataset_metrics)

completion <- list(
  status = "PASS_HSTAR_COMPACTNESS_VS_CELL2021_AUDIT", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  current_model = ibd_root, comparison_model = cell2021_root, cell_types = nrow(cell_metrics), samples = nrow(sample_metrics)
)
jsonlite::write_json(completion, file.path(result_dir, "hstar_compactness_audit_complete.json"), pretty = TRUE, auto_unbox = TRUE)
