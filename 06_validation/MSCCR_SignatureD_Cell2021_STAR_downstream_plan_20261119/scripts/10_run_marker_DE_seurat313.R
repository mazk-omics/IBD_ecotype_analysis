#!/usr/bin/env Rscript
# 本脚本用 Seurat 3.1.3 对一个 dataset × parent type 做全基因 marker 检验。
# 正式 assignment 固定；每个分析层先排除少于 5 个 assigned cells 的 state 及其细胞。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取一块全转录组 counts
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 10_run_marker_DE_seurat313.R RUN_ROOT DATASET CELL_TYPE", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
dataset <- args[[2]]
cell_type <- args[[3]]
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite", "Matrix", "Seurat")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
if (as.character(packageVersion("Seurat")) != "3.1.3") stop("Marker DE requires Seurat 3.1.3", call. = FALSE)
library(data.table)
library(Matrix)
library(Seurat)

input_dir <- file.path(output_root, "09_marker_recovery/full_counts", dataset, cell_type)
count_file <- file.path(input_dir, "counts.rds")
cell_file <- file.path(input_dir, "cells.tsv.gz")
if (!file.exists(count_file) || !file.exists(cell_file)) stop("Marker count input is incomplete: ", dataset, " ", cell_type, call. = FALSE)

counts <- readRDS(count_file)
cells <- fread(cell_file, colClasses = "character")
cells[, Assigned := Assigned == "TRUE"]
cells[, SelectionOrder := as.integer(SelectionOrder)]
if (!methods::is(counts, "dgCMatrix") || ncol(counts) != nrow(cells) || !identical(colnames(counts), cells$CellID)) stop("Count/cell identity mismatch", call. = FALSE)
if (anyDuplicated(rownames(counts)) || anyDuplicated(cells$CellID) || any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8) || any(Matrix::colSums(counts) <= 0)) stop("Invalid full-transcriptome counts", call. = FALSE)
cat("MARKER_DE_INPUT", dataset, cell_type, nrow(counts), "genes", ncol(counts), "cells", sum(cells$Assigned), "assigned\n")

rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
if (is.na(rank_lookup[[cell_type]])) stop("Frozen rank is missing: ", cell_type, call. = FALSE)
gene_info_file <- file.path(state_discovery_root, cell_type, rank_lookup[[cell_type]], "gene_info.txt")
gene_info <- as.data.table(read.delim(gene_info_file, check.names = FALSE, stringsAsFactors = FALSE))
gene_info[, MaxFC := as.numeric(MaxFC)]
gene_info <- gene_info[State != "-1" & !is.na(State)]
if (!all(c("Gene", "State", "MaxFC") %in% names(gene_info)) || anyDuplicated(gene_info$Gene)) stop("Formal gene_info is invalid: ", cell_type, call. = FALSE)
candidate_genes <- intersect(rownames(counts), gene_info$Gene)
if (!length(candidate_genes)) stop("No formal discovery genes overlap full counts", call. = FALSE)

output_dir <- file.path(output_root, "09_marker_recovery/DE_checkpoints", dataset, cell_type)
completed_file <- file.path(output_dir, "completed.json")
if (file.exists(completed_file)) stop("Marker DE checkpoint already complete; do not overwrite", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step2. 复核 Seurat 3.1.3 数值定义
set.seed(20261122)
toy <- matrix(rpois(120, 4), nrow = 10, dimnames = list(paste0("g", 1:10), paste0("c", 1:12)))
toy[1, 1:6] <- 40
toy_log <- NormalizeData(toy, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
toy_expected <- log1p(sweep(toy, 2, colSums(toy), "/") * 10000)
if (max(abs(toy_log - toy_expected)) >= 1e-10) stop("Seurat NormalizeData numerical smoke failed", call. = FALSE)
toy_de <- suppressWarnings(FindMarkers(toy_log, cells.1 = colnames(toy)[1:6], cells.2 = colnames(toy)[7:12], min.pct = 0.1, logfc.threshold = 0.05, verbose = FALSE))
toy_fc <- log(Matrix::rowMeans(expm1(toy_log[, 1:6])) + 1) - log(Matrix::rowMeans(expm1(toy_log[, 7:12])) + 1)
if (max(abs(toy_de$avg_logFC - toy_fc[rownames(toy_de)])) >= 1e-10) stop("Seurat native lnFC smoke failed", call. = FALSE)
if (max(abs(toy_de$p_val_adj - pmin(toy_de$p_val * nrow(toy), 1))) >= 1e-10) stop("Seurat Bonferroni smoke failed", call. = FALSE)

# Step3. 归一化并保存 n1 所需的平均 CPM
data <- NormalizeData(counts, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
mean_cpm <- Matrix::rowMeans(counts %*% Diagonal(x = 1e6 / Matrix::colSums(counts)))
gene_expression <- data.table(Dataset = dataset, CellType = cell_type, Gene = rownames(counts), MeanCPM = as.numeric(mean_cpm))
fwrite(gene_expression, file.path(output_dir, "gene_expression.tsv.gz"), sep = "\t", quote = FALSE)

scope_table <- rbind(
  data.table(Scope = "Dataset", Stratum = "All"),
  data.table(Scope = "Disease", Stratum = sort(unique(cells$Disease))),
  data.table(Scope = "Region3", Stratum = sort(unique(cells$Region3)))
)
scope_table <- unique(scope_table)
eligibility_parts <- list()

# Step4. 逐层和逐 state 做 one-vs-rest 检验
for (scope_index in seq_len(nrow(scope_table))) {
  scope <- scope_table$Scope[[scope_index]]
  stratum <- scope_table$Stratum[[scope_index]]
  group_cells <- if (scope == "Dataset") cells else cells[get(scope) == stratum]
  state_counts <- table(group_cells$State[group_cells$Assigned])
  allowed <- names(state_counts)[state_counts >= 5L]
  eligible_cells <- group_cells[Assigned == TRUE & State %in% allowed]
  eligibility <- data.table(Dataset = dataset, CellType = cell_type, Scope = scope, Stratum = stratum, State = names(state_counts), AssignedCells = as.integer(state_counts))
  eligibility[, `:=`(EligibleGE5 = State %in% allowed, EvaluableContrast = length(allowed) >= 2L)]
  eligibility_parts[[paste(scope, stratum, sep = "__")]] <- eligibility

  tag <- paste(scope, gsub("[^A-Za-z0-9]+", "_", stratum), sep = "__")
  if (length(allowed) < 2L) {
    writeLines("NOT_EVALUABLE: fewer than two states have at least 5 assigned cells", file.path(output_dir, paste0(tag, "__not_evaluable.txt")))
    next
  }
  cat("MARKER_DE_SCOPE", dataset, cell_type, scope, stratum, nrow(eligible_cells), "cells", length(allowed), "states\n")

  de_parts <- list()
  fc_parts <- list()
  for (state in allowed) {
    own <- eligible_cells[State == state]$CellID
    other <- eligible_cells[State != state]$CellID
    fold <- log(Matrix::rowMeans(expm1(data[candidate_genes, own, drop = FALSE])) + 1) - log(Matrix::rowMeans(expm1(data[candidate_genes, other, drop = FALSE])) + 1)
    fc_parts[[state]] <- data.table(Dataset = dataset, CellType = cell_type, Scope = scope, Stratum = stratum, Gene = candidate_genes, State = state, log2FC = as.numeric(fold / log(2)), OwnCells = length(own), OtherCells = length(other))

    result <- tryCatch(
      suppressWarnings(FindMarkers(data[, eligible_cells$CellID, drop = FALSE], cells.1 = own, cells.2 = other, min.pct = 0.1, logfc.threshold = 0.05, test.use = "wilcox", verbose = FALSE)),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      if (!grepl("No features pass", conditionMessage(result))) stop(result)
      writeLines(conditionMessage(result), file.path(output_dir, paste0(tag, "__", state, "__empty.txt")))
      next
    }
    result <- as.data.table(result, keep.rownames = "Gene")
    result[, `:=`(Dataset = dataset, CellType = cell_type, Scope = scope, Stratum = stratum, State = state, log2FC = avg_logFC / log(2), OwnCells = length(own), OtherCells = length(other))]
    setcolorder(result, c("Dataset", "CellType", "Scope", "Stratum", "State", "Gene", "p_val", "avg_logFC", "log2FC", "pct.1", "pct.2", "p_val_adj", "OwnCells", "OtherCells"))
    de_parts[[state]] <- result
  }
  fwrite(rbindlist(fc_parts), file.path(output_dir, paste0(tag, "__candidate_full_FC.tsv.gz")), sep = "\t", quote = FALSE)
  if (length(de_parts)) fwrite(rbindlist(de_parts), file.path(output_dir, paste0(tag, "__DEG.tsv.gz")), sep = "\t", quote = FALSE)
}

# Step5. 保存 checkpoint 身份和 session
eligibility_all <- rbindlist(eligibility_parts, fill = TRUE)
fwrite(eligibility_all, file.path(output_dir, "state_eligibility.tsv"), sep = "\t", quote = FALSE)
summary <- list(
  status = "PASS_MARKER_DE_SEURAT313",
  dataset = dataset,
  cell_type = cell_type,
  genes = nrow(counts),
  candidate_discovery_genes = length(candidate_genes),
  frozen_cells = nrow(cells),
  assigned_cells = sum(cells$Assigned),
  dataset_states_ge5 = sum(eligibility_all$Scope == "Dataset" & eligibility_all$EligibleGE5),
  seurat_version = as.character(packageVersion("Seurat")),
  normalization = "ln(CP10k+1)",
  native_fold_change = "natural log; log2FC is avg_logFC/log(2)",
  multiple_testing = "Seurat 3.1.3 default Bonferroni over all input genes",
  warning_policy = "expected exact-Wilcoxon tie warnings suppressed; Seurat test and returned values unchanged",
  count_sha256 = digest::digest(count_file, algo = "sha256", file = TRUE, serialize = FALSE),
  cell_sha256 = digest::digest(cell_file, algo = "sha256", file = TRUE, serialize = FALSE)
)
jsonlite::write_json(summary, completed_file, pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
cat("MARKER_DE_SEURAT313_PASS", dataset, cell_type, sum(eligibility_all$Scope == "Dataset" & eligibility_all$EligibleGE5), "dataset states >=5\n")
