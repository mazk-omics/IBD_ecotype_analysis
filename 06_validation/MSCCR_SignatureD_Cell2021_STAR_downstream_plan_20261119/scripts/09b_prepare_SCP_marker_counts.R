#!/usr/bin/env Rscript
# 本脚本从 SCP1884 或 SCP259 原始 UMI bundle 构建 marker validation 的全转录组稀疏输入。
# 每个 bundle 只读一次；不同 bundle 的基因集合在父类型内部取并集并显式补零。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞和正式 assignment
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 09b_prepare_SCP_marker_counts.R RUN_ROOT DATASET", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
dataset <- args[[2]]
if (!dataset %in% c("SCP1884", "SCP259")) stop("DATASET must be SCP1884 or SCP259", call. = FALSE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("arrow", "data.table", "digest", "jsonlite", "Matrix")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(Matrix)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
selected_cells <- fread(file.path(audit_dir, "recovery_selected_cells.tsv.gz"), colClasses = "character")
selected_cells <- selected_cells[Dataset == dataset]
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]

eligible_file <- file.path(eligible_pool_root, paste0(dataset, "_eligible_ibd_cells.parquet"))
eligible <- as.data.table(arrow::read_parquet(eligible_file, as_data_frame = TRUE))
if (!all(c("CellID", "SourceBundle") %in% names(eligible))) stop("Eligible pool lacks SourceBundle: ", dataset, call. = FALSE)
selected_cells <- merge(selected_cells, unique(eligible[, .(CellID, SourceBundle)]), by = "CellID", all.x = TRUE, sort = FALSE)
selected_cells <- selected_cells[match(fread(file.path(audit_dir, "recovery_selected_cells.tsv.gz"), colClasses = "character")[Dataset == dataset]$CellID, CellID)]
if (anyNA(selected_cells$SourceBundle) || anyDuplicated(selected_cells$CellID)) stop("Frozen cells do not map uniquely to source bundles: ", dataset, call. = FALSE)

expected_cells <- if (dataset == "SCP1884") 48047L else 33755L
expected_types <- if (dataset == "SCP1884") 20L else 18L
if (nrow(selected_cells) != expected_cells || uniqueN(selected_cells$CellType) != expected_types) stop("Frozen cell inventory changed: ", dataset, call. = FALSE)

sample_metadata <- fread(file.path(audit_dir, "donor_sample_phenotype_freeze.tsv"), colClasses = "character")
sample_metadata <- sample_metadata[Dataset == dataset]
if (anyDuplicated(sample_metadata[, .(DonorID, BiologicalSampleID)])) stop("Frozen phenotype identity changed: ", dataset, call. = FALSE)

if (dataset == "SCP1884") {
  job_file <- file.path(output_root, "05c_SCP1884_official_recovery/SCP1884_official_recovery_jobs_completed.tsv")
  source_root <- "/home/mazekai/enteric_glia/data_scRNA/SCP1884/expression"
  bundle_config <- data.table(
    SourceBundle = c("CO_EPI", "CO_IMM", "CO_STR", "TI_EPI", "TI_IMM", "TI_STR"),
    Subdir = c("Colon/Epithelial", "Colon/Immune", "Colon/Stromal", "TerminalIleum/Epithelial", "TerminalIleum/Immune", "TerminalIleum/Stromal"),
    Prefix = c("CO_EPI", "CO_IMM", "CO_STR", "TI_EPI", "TI_IMM", "TI_STR")
  )
} else {
  job_file <- file.path(output_root, "06c_SCP259_official_recovery/SCP259_official_recovery_jobs_completed.tsv")
  source_root <- file.path(project_root, "01_raw_data/scRNA/SCP259/raw/SCP259/expression")
  bundle_config <- data.table(
    SourceBundle = c("Epithelial", "Stromal", "Immune"),
    Subdir = c("5cdc540d328cee7a2efc2348", "5cdc540d328cee7a2efc2349", "5cdc540d328cee7a2efc234a"),
    Prefix = c("Epi", "Fib", "Imm")
  )
}

jobs <- fread(job_file, colClasses = "character")
if (nrow(jobs) != expected_types || any(jobs$Dataset != dataset) || any(jobs$ExitStatus != "0") || anyDuplicated(jobs$CellType)) stop("Completed recovery manifest changed: ", dataset, call. = FALSE)

assignment_parts <- list()
for (cell_type in jobs$CellType) {
  assignment_file <- file.path(jobs[CellType == cell_type]$ResultDir, "state_assignment.txt")
  assignment <- fread(assignment_file, colClasses = "character")
  if (!identical(names(assignment), c("ID", "State", "InitialState")) || anyDuplicated(assignment$ID)) stop("Invalid state assignment: ", dataset, " ", cell_type, call. = FALSE)
  if (!all(assignment$ID %in% selected_cells[CellType == cell_type]$CellID)) stop("State assignment contains non-frozen cells: ", dataset, " ", cell_type, call. = FALSE)
  assignment_parts[[cell_type]] <- assignment
}
assignment_all <- rbindlist(assignment_parts)
if (anyDuplicated(assignment_all$ID)) stop("Assignment IDs overlap across parent types: ", dataset, call. = FALSE)

cells <- merge(selected_cells, sample_metadata[, .(DonorID, BiologicalSampleID, Disease, Inflammation, Region3)], by = c("DonorID", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
cells <- merge(cells, assignment_all[, .(CellID = ID, State, InitialState)], by = "CellID", all.x = TRUE, sort = FALSE)
cells <- cells[match(selected_cells$CellID, CellID)]
cells[is.na(State), State := "Unassigned"]
cells[is.na(InitialState), InitialState := "Unassigned"]
cells[, Assigned := State != "Unassigned"]
if (anyNA(cells$Disease) || anyNA(cells$Inflammation) || anyNA(cells$Region3) || !identical(cells$CellID, selected_cells$CellID)) stop("Cell metadata join failed: ", dataset, call. = FALSE)
print(cells[, .(FrozenCells = .N, AssignedCells = sum(Assigned), States = uniqueN(State[Assigned])), by = CellType][order(CellType)])

output_dir <- file.path(output_root, "09_marker_recovery/full_counts", dataset)
manifest_file <- file.path(output_dir, paste0(dataset, "_marker_counts_manifest.tsv"))
if (file.exists(manifest_file)) stop(dataset, " marker counts already exist; do not overwrite", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cell_types <- sort(unique(cells$CellType))
matrix_parts <- setNames(vector("list", length(cell_types)), cell_types)
id_parts <- setNames(vector("list", length(cell_types)), cell_types)
library_parts <- setNames(vector("list", length(cell_types)), cell_types)
bundle_audit_parts <- list()

# Step2. 逐 bundle 抽取冻结细胞的全部基因 counts
for (bundle_index in seq_len(nrow(bundle_config))) {
  source_bundle <- bundle_config$SourceBundle[[bundle_index]]
  prefix <- bundle_config$Prefix[[bundle_index]]
  bundle_dir <- file.path(source_root, bundle_config$Subdir[[bundle_index]])
  selected_bundle <- cells[SourceBundle == source_bundle]
  if (!nrow(selected_bundle)) stop("No frozen cells in source bundle: ", dataset, " ", source_bundle, call. = FALSE)

  if (dataset == "SCP1884") {
    feature_file <- file.path(bundle_dir, paste0(prefix, ".features.tsv"))
    barcode_file <- file.path(bundle_dir, paste0(prefix, ".barcodes.tsv"))
    matrix_file <- file.path(bundle_dir, "raw.mtx")
    features <- fread(feature_file, header = FALSE, select = 1:2, col.names = c("GeneID", "Gene"), na.strings = character())
    barcodes <- fread(barcode_file, header = FALSE, select = 1, col.names = "CellID")
  } else {
    feature_file <- file.path(bundle_dir, paste0(prefix, ".genes.tsv"))
    barcode_file <- file.path(bundle_dir, paste0(prefix, ".barcodes2.tsv"))
    matrix_file <- file.path(bundle_dir, paste0("gene_sorted-", prefix, ".matrix.mtx"))
    features <- fread(feature_file, header = FALSE, select = 1, col.names = "Gene", na.strings = character())
    features[, GeneID := Gene]
    setcolorder(features, c("GeneID", "Gene"))
    barcodes <- fread(barcode_file, header = FALSE, select = 1, col.names = "CellID")
  }
  features[is.na(Gene) | !nzchar(Gene), Gene := paste0("UNMAPPED__", GeneID)]
  features[, FeatureKey := paste(GeneID, Gene, sep = "|")]
  if (anyDuplicated(features$FeatureKey) || anyDuplicated(barcodes$CellID)) stop("Duplicated source identity: ", dataset, " ", source_bundle, call. = FALSE)
  raw_cell_index <- match(selected_bundle$CellID, barcodes$CellID)
  if (anyNA(raw_cell_index) || anyDuplicated(raw_cell_index)) stop("Frozen cells do not map uniquely: ", dataset, " ", source_bundle, call. = FALSE)

  gene_symbols <- unique(features$Gene)
  source_gene_to_symbol <- match(features$Gene, gene_symbols)
  raw <- Matrix::readMM(matrix_file)
  if (!methods::is(raw, "TsparseMatrix")) raw <- methods::as(raw, "dgTMatrix")
  if (nrow(raw) != nrow(features) || ncol(raw) != nrow(barcodes)) stop("Raw matrix dimensions changed: ", dataset, " ", source_bundle, call. = FALSE)
  if (any(!is.finite(raw@x)) || any(raw@x < 0) || any(abs(raw@x - round(raw@x)) > 1e-8)) stop("Raw matrix is not nonnegative integer counts: ", dataset, " ", source_bundle, call. = FALSE)

  if (dataset == "SCP1884") {
    source_cell_to_selected <- integer(nrow(barcodes))
    source_cell_to_selected[raw_cell_index] <- seq_len(nrow(selected_bundle))
    symbol_counts <- sparseMatrix(i = integer(0), j = integer(0), dims = c(length(gene_symbols), nrow(selected_bundle)), dimnames = list(gene_symbols, selected_bundle$CellID))
    nnz_block <- 5000000L
    block_starts <- seq.int(1L, length(raw@x), by = nnz_block)
    for (block_index in seq_along(block_starts)) {
      from <- block_starts[[block_index]]
      to <- min(from + nnz_block - 1L, length(raw@x))
      source_cell <- raw@j[from:to] + 1L
      selected_cell <- source_cell_to_selected[source_cell]
      keep <- selected_cell > 0L
      if (any(keep)) {
        source_gene <- raw@i[from:to][keep] + 1L
        partial <- sparseMatrix(i = source_gene_to_symbol[source_gene], j = selected_cell[keep], x = raw@x[from:to][keep], dims = dim(symbol_counts))
        symbol_counts <- symbol_counts + partial
        rm(partial)
      }
      if (block_index %% 10L == 0L || block_index == length(block_starts)) cat("SCP_MARKER_NNZ_BLOCK", dataset, source_bundle, block_index, "/", length(block_starts), "\n")
    }
  } else {
    selected_raw <- raw[, raw_cell_index, drop = FALSE]
    collapse <- sparseMatrix(i = source_gene_to_symbol, j = seq_len(nrow(features)), x = 1, dims = c(length(gene_symbols), nrow(features)), dimnames = list(gene_symbols, features$FeatureKey))
    symbol_counts <- collapse %*% selected_raw
    colnames(symbol_counts) <- selected_bundle$CellID
    rm(selected_raw, collapse)
  }

  full_library <- Matrix::colSums(symbol_counts)
  if (any(!is.finite(full_library)) || any(full_library <= 0)) stop("Invalid exported library size: ", dataset, " ", source_bundle, call. = FALSE)
  for (cell_type in sort(unique(selected_bundle$CellType))) {
    type_cells <- selected_bundle[CellType == cell_type]
    local_column <- match(type_cells$CellID, selected_bundle$CellID)
    type_counts <- methods::as(symbol_counts[, local_column, drop = FALSE], "dgCMatrix")
    colnames(type_counts) <- type_cells$CellID
    matrix_parts[[cell_type]][[source_bundle]] <- type_counts
    id_parts[[cell_type]][[source_bundle]] <- type_cells$CellID
    library_parts[[cell_type]][[source_bundle]] <- full_library[local_column]
  }
  bundle_audit_parts[[source_bundle]] <- data.table(Dataset = dataset, SourceBundle = source_bundle, SourceGenes = nrow(raw), SourceCells = ncol(raw), SourceNNZ = length(raw@x), ExportGenes = nrow(symbol_counts), SelectedCells = nrow(selected_bundle), MinLibrarySize = min(full_library), MedianLibrarySize = median(full_library), MaxLibrarySize = max(full_library))
  cat("SCP_MARKER_BUNDLE_EXTRACTED", dataset, source_bundle, nrow(selected_bundle), "cells", nrow(symbol_counts), "genes\n")
  rm(raw, symbol_counts, full_library, barcodes, features)
  invisible(gc(FALSE))
}

# Step3. 父类型内对齐 bundle 基因并保存唯一稀疏矩阵
manifest_parts <- list()
state_audit_parts <- list()
for (cell_type in cell_types) {
  frozen_type <- cells[CellType == cell_type][order(SelectionOrder)]
  genes_type <- sort(unique(unlist(lapply(matrix_parts[[cell_type]], rownames), use.names = FALSE)))
  aligned_parts <- list()
  for (source_bundle in names(matrix_parts[[cell_type]])) {
    part <- matrix_parts[[cell_type]][[source_bundle]]
    align <- sparseMatrix(i = match(rownames(part), genes_type), j = seq_len(nrow(part)), x = 1, dims = c(length(genes_type), nrow(part)), dimnames = list(genes_type, rownames(part)))
    aligned_parts[[source_bundle]] <- align %*% part
  }
  counts_type <- methods::as(do.call(cbind, aligned_parts), "dgCMatrix")
  ids_type <- unlist(id_parts[[cell_type]], use.names = FALSE)
  library_type <- unlist(library_parts[[cell_type]], use.names = FALSE)
  frozen_order <- match(frozen_type$CellID, ids_type)
  if (anyNA(frozen_order) || anyDuplicated(ids_type) || ncol(counts_type) != nrow(frozen_type)) stop("Frozen cell assembly failed: ", dataset, " ", cell_type, call. = FALSE)
  counts_type <- counts_type[, frozen_order, drop = FALSE]
  library_type <- library_type[frozen_order]
  colnames(counts_type) <- frozen_type$CellID
  if (!identical(colnames(counts_type), frozen_type$CellID) || !isTRUE(all.equal(as.numeric(Matrix::colSums(counts_type)), as.numeric(library_type), tolerance = 0))) stop("Saved count audit failed: ", dataset, " ", cell_type, call. = FALSE)

  type_dir <- file.path(output_dir, cell_type)
  dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)
  count_file <- file.path(type_dir, "counts.rds")
  cell_file <- file.path(type_dir, "cells.tsv.gz")
  library_file <- file.path(type_dir, "library_size.tsv.gz")
  saveRDS(counts_type, count_file, compress = "gzip")
  fwrite(frozen_type, cell_file, sep = "\t", quote = FALSE)
  fwrite(data.table(CellID = frozen_type$CellID, RawLibrarySize = library_type, ExportedLibrarySize = Matrix::colSums(counts_type)), library_file, sep = "\t", quote = FALSE)

  state_audit <- frozen_type[Assigned == TRUE, .(AssignedCells = .N), by = State]
  state_audit[, `:=`(Dataset = dataset, CellType = cell_type, EligibleGE5 = AssignedCells >= 5L)]
  setcolorder(state_audit, c("Dataset", "CellType", "State", "AssignedCells", "EligibleGE5"))
  state_audit_parts[[cell_type]] <- state_audit
  manifest_parts[[cell_type]] <- data.table(
    Dataset = dataset, CellType = cell_type, FrozenCells = nrow(frozen_type), AssignedCells = sum(frozen_type$Assigned), UnassignedCells = sum(!frozen_type$Assigned),
    AssignedStates = uniqueN(frozen_type[Assigned == TRUE]$State), EligibleStatesGE5 = sum(state_audit$EligibleGE5), Genes = nrow(counts_type), NonzeroGenes = sum(Matrix::rowSums(counts_type) > 0),
    NNZ = length(counts_type@x), MinLibrarySize = min(library_type), MedianLibrarySize = median(library_type), MaxLibrarySize = max(library_type), CountsFile = count_file,
    CountsBytes = file.info(count_file)$size, CountsSHA256 = digest::digest(count_file, algo = "sha256", file = TRUE, serialize = FALSE), CellsFile = cell_file,
    CellsSHA256 = digest::digest(cell_file, algo = "sha256", file = TRUE, serialize = FALSE)
  )
  cat("SCP_MARKER_PARENT_WRITTEN", dataset, cell_type, nrow(frozen_type), "cells", nrow(counts_type), "genes", length(counts_type@x), "nnz\n")
  rm(counts_type, aligned_parts)
  invisible(gc(FALSE))
}

# Step4. 写入审计清单
manifest <- rbindlist(manifest_parts)
state_audit <- rbindlist(state_audit_parts)
bundle_audit <- rbindlist(bundle_audit_parts)
fwrite(manifest, manifest_file, sep = "\t")
fwrite(state_audit, file.path(output_dir, paste0(dataset, "_marker_state_eligibility.tsv")), sep = "\t")
fwrite(bundle_audit, file.path(output_dir, paste0(dataset, "_marker_bundle_audit.tsv")), sep = "\t")

summary <- list(
  status = paste0("PASS_", dataset, "_MARKER_FULL_COUNTS"), dataset = dataset, frozen_cells = nrow(cells), assigned_cells = sum(cells$Assigned), parent_types = nrow(manifest),
  eligible_state_dataset_pairs = sum(state_audit$EligibleGE5), normalization = "raw UMI counts; Seurat NormalizeData applies ln(CP10k+1) in the next step",
  gene_identity = "duplicate gene symbols summed within bundle; missing symbols retain UNMAPPED__GeneID; parent-type bundle union is zero-filled",
  storage = "one gzip-compressed dgCMatrix RDS per dataset x parent type; no duplicate MTX or NPZ", recovery_model_changed = FALSE
)
jsonlite::write_json(summary, file.path(output_dir, paste0(dataset, "_marker_counts_summary.json")), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
cat(dataset, "MARKER_FULL_COUNTS_PASS", nrow(cells), "frozen cells", sum(cells$Assigned), "assigned cells", sum(state_audit$EligibleGE5), "eligible state-dataset pairs\n")
