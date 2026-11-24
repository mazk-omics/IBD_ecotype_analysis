#!/usr/bin/env Rscript
# 本脚本从 GSE214695 原始 10x counts 构建 marker validation 的全转录组稀疏输入。
# 只读取本轮冻结细胞，并拼接正式 EcoTyper state assignment；不修改 recovery 结果。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞和正式 assignment
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 09a_prepare_GSE214695_marker_counts.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite", "Matrix")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(Matrix)

dataset <- "GSE214695"
audit_dir <- file.path(output_root, "03_scrna_input_audit")
selected_cells <- fread(file.path(audit_dir, "recovery_selected_cells.tsv.gz"), colClasses = "character")
selected_cells <- selected_cells[Dataset == dataset]
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]
if (nrow(selected_cells) != 17109L || anyDuplicated(selected_cells$CellID)) stop("GSE214695 frozen cell identity changed", call. = FALSE)

sample_metadata <- fread(file.path(audit_dir, "donor_sample_phenotype_freeze.tsv"), colClasses = "character")
sample_metadata <- sample_metadata[Dataset == dataset]
if (anyDuplicated(sample_metadata[, .(DonorID, BiologicalSampleID)])) stop("GSE214695 phenotype identity changed", call. = FALSE)

job_file <- file.path(output_root, "04c_GSE214695_official_recovery/GSE214695_official_recovery_jobs_completed.tsv")
jobs <- fread(job_file, colClasses = "character")
if (nrow(jobs) != 18L || any(jobs$Dataset != dataset) || any(jobs$ExitStatus != "0") || anyDuplicated(jobs$CellType)) stop("GSE214695 completed recovery manifest changed", call. = FALSE)

assignment_parts <- list()
for (cell_type in jobs$CellType) {
  result_dir <- jobs[CellType == cell_type]$ResultDir
  assignment_file <- file.path(result_dir, "state_assignment.txt")
  assignment <- fread(assignment_file, colClasses = "character")
  if (!identical(names(assignment), c("ID", "State", "InitialState")) || anyDuplicated(assignment$ID)) stop("Invalid state assignment: ", cell_type, call. = FALSE)
  if (!all(assignment$ID %in% selected_cells[CellType == cell_type]$CellID)) stop("State assignment contains non-frozen cells: ", cell_type, call. = FALSE)
  assignment[, CellType := cell_type]
  assignment_parts[[cell_type]] <- assignment
}
assignment_all <- rbindlist(assignment_parts)
if (anyDuplicated(assignment_all$ID)) stop("GSE214695 assignment IDs overlap across parent types", call. = FALSE)

cells <- merge(selected_cells, sample_metadata[, .(DonorID, BiologicalSampleID, Disease, Inflammation, Region3)], by = c("DonorID", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
cells <- merge(cells, assignment_all[, .(CellID = ID, State, InitialState)], by = "CellID", all.x = TRUE, sort = FALSE)
cells <- cells[match(selected_cells$CellID, CellID)]
cells[is.na(State), State := "Unassigned"]
cells[is.na(InitialState), InitialState := "Unassigned"]
cells[, Assigned := State != "Unassigned"]
if (anyNA(cells$Disease) || anyNA(cells$Inflammation) || anyNA(cells$Region3) || !identical(cells$CellID, selected_cells$CellID)) stop("GSE214695 cell metadata join failed", call. = FALSE)
print(cells[, .(FrozenCells = .N, AssignedCells = sum(Assigned), States = uniqueN(State[Assigned])), by = CellType][order(CellType)])

# Step2. 定位原始 10x 文件
source_root <- file.path(project_root, "06_validation/GSE214695_scRNA_recovery_20261009/source")
matrix_files <- sort(list.files(source_root, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE))
sample_labels <- sub("^GSM[0-9]+_(.*)_matrix\\.mtx\\.gz$", "\\1", basename(matrix_files))
sample_ids <- gsub("-", "", sample_labels)
keep_ibd <- sample_ids %in% unique(cells$BiologicalSampleID)
matrix_files <- matrix_files[keep_ibd]
sample_ids <- sample_ids[keep_ibd]
if (!setequal(sample_ids, unique(cells$BiologicalSampleID)) || anyDuplicated(sample_ids)) stop("GSE214695 raw sample identity changed", call. = FALSE)

output_dir <- file.path(output_root, "09_marker_recovery/full_counts", dataset)
manifest_file <- file.path(output_dir, "GSE214695_marker_counts_manifest.tsv")
if (file.exists(manifest_file)) stop("GSE214695 marker counts already exist; do not overwrite", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cell_types <- sort(unique(cells$CellType))
matrix_parts <- setNames(vector("list", length(cell_types)), cell_types)
id_parts <- setNames(vector("list", length(cell_types)), cell_types)
raw_library_parts <- setNames(vector("list", length(cell_types)), cell_types)
symbol_library_parts <- setNames(vector("list", length(cell_types)), cell_types)
feature_reference <- NULL
gene_reference <- NULL

# Step3. 每个样本只读一次并按基因符号合并
for (sample_index in seq_along(matrix_files)) {
  matrix_file <- matrix_files[[sample_index]]
  sample_id <- sample_ids[[sample_index]]
  file_prefix <- sub("_matrix\\.mtx\\.gz$", "", matrix_file)
  feature_file <- paste0(file_prefix, "_features.tsv.gz")
  barcode_file <- paste0(file_prefix, "_barcodes.tsv.gz")

  feature <- fread(feature_file, header = FALSE, select = 1:2, col.names = c("GeneID", "Gene"), na.strings = character())
  barcode <- fread(barcode_file, header = FALSE, select = 1, col.names = "Barcode")
  feature[is.na(Gene) | !nzchar(Gene), Gene := paste0("UNMAPPED__", GeneID)]
  feature[, FeatureKey := paste(GeneID, Gene, sep = "|")]
  if (anyDuplicated(feature$FeatureKey)) stop("Duplicated source feature identity: ", sample_id, call. = FALSE)
  if (is.null(feature_reference)) {
    feature_reference <- feature$FeatureKey
    gene_reference <- unique(feature$Gene)
  }
  if (!identical(feature$FeatureKey, feature_reference)) stop("GSE214695 feature order differs across samples", call. = FALSE)

  sample_cells <- cells[BiologicalSampleID == sample_id]
  sample_cells[, RawBarcode := sub("^.*\\|", "", CellID)]
  column_index <- match(sample_cells$RawBarcode, barcode$Barcode)
  if (anyNA(column_index) || anyDuplicated(column_index)) stop("Frozen barcodes do not map uniquely: ", sample_id, call. = FALSE)

  counts <- Matrix::readMM(gzfile(matrix_file))
  if (nrow(counts) != nrow(feature) || ncol(counts) != nrow(barcode)) stop("Raw matrix dimensions changed: ", sample_id, call. = FALSE)
  selected_counts <- counts[, column_index, drop = FALSE]
  raw_library <- Matrix::colSums(selected_counts)
  collapse <- sparseMatrix(i = match(feature$Gene, gene_reference), j = seq_len(nrow(feature)), x = 1, dims = c(length(gene_reference), nrow(feature)), dimnames = list(gene_reference, feature$FeatureKey))
  symbol_counts <- collapse %*% selected_counts
  symbol_library <- Matrix::colSums(symbol_counts)
  if (any(!is.finite(raw_library)) || any(raw_library <= 0) || !isTRUE(all.equal(raw_library, symbol_library, tolerance = 0))) stop("Gene-symbol collapse changed library size: ", sample_id, call. = FALSE)

  for (cell_type in unique(sample_cells$CellType)) {
    type_cells <- sample_cells[CellType == cell_type]
    local_column <- match(type_cells$CellID, sample_cells$CellID)
    type_counts <- symbol_counts[, local_column, drop = FALSE]
    colnames(type_counts) <- type_cells$CellID
    matrix_parts[[cell_type]][[sample_id]] <- type_counts
    id_parts[[cell_type]][[sample_id]] <- type_cells$CellID
    raw_library_parts[[cell_type]][[sample_id]] <- raw_library[local_column]
    symbol_library_parts[[cell_type]][[sample_id]] <- symbol_library[local_column]
  }
  cat("GSE_MARKER_SAMPLE_EXTRACTED", sample_id, nrow(sample_cells), "cells", nrow(symbol_counts), "genes\n")
  rm(counts, selected_counts, symbol_counts, collapse, barcode)
  invisible(gc(FALSE))
}

# Step4. 按父类型保存唯一一份压缩稀疏矩阵
manifest_parts <- list()
state_audit_parts <- list()
for (cell_type in cell_types) {
  frozen_type <- cells[CellType == cell_type][order(SelectionOrder)]
  counts_type <- do.call(cbind, matrix_parts[[cell_type]])
  ids_type <- unlist(id_parts[[cell_type]], use.names = FALSE)
  raw_library <- unlist(raw_library_parts[[cell_type]], use.names = FALSE)
  symbol_library <- unlist(symbol_library_parts[[cell_type]], use.names = FALSE)
  frozen_order <- match(frozen_type$CellID, ids_type)
  if (anyNA(frozen_order) || anyDuplicated(ids_type) || ncol(counts_type) != nrow(frozen_type)) stop("Frozen cell assembly failed: ", cell_type, call. = FALSE)
  counts_type <- methods::as(counts_type[, frozen_order, drop = FALSE], "dgCMatrix")
  raw_library <- raw_library[frozen_order]
  symbol_library <- symbol_library[frozen_order]
  colnames(counts_type) <- frozen_type$CellID
  if (!identical(colnames(counts_type), frozen_type$CellID) || !isTRUE(all.equal(as.numeric(Matrix::colSums(counts_type)), as.numeric(raw_library), tolerance = 0))) stop("Saved count audit failed: ", cell_type, call. = FALSE)

  type_dir <- file.path(output_dir, cell_type)
  dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)
  count_file <- file.path(type_dir, "counts.rds")
  cell_file <- file.path(type_dir, "cells.tsv.gz")
  library_file <- file.path(type_dir, "library_size.tsv.gz")
  saveRDS(counts_type, count_file, compress = "gzip")
  fwrite(frozen_type, cell_file, sep = "\t", quote = FALSE)
  fwrite(data.table(CellID = frozen_type$CellID, RawLibrarySize = raw_library, ExportedLibrarySize = symbol_library), library_file, sep = "\t", quote = FALSE)

  state_audit <- frozen_type[Assigned == TRUE, .(AssignedCells = .N), by = State]
  state_audit[, `:=`(Dataset = dataset, CellType = cell_type, EligibleGE5 = AssignedCells >= 5L)]
  setcolorder(state_audit, c("Dataset", "CellType", "State", "AssignedCells", "EligibleGE5"))
  state_audit_parts[[cell_type]] <- state_audit

  manifest_parts[[cell_type]] <- data.frame(
    Dataset = dataset,
    CellType = cell_type,
    FrozenCells = nrow(frozen_type),
    AssignedCells = sum(frozen_type$Assigned),
    UnassignedCells = sum(!frozen_type$Assigned),
    AssignedStates = uniqueN(frozen_type[Assigned == TRUE]$State),
    EligibleStatesGE5 = sum(state_audit$EligibleGE5),
    Genes = nrow(counts_type),
    NonzeroGenes = sum(Matrix::rowSums(counts_type) > 0),
    NNZ = length(counts_type@x),
    MinLibrarySize = min(raw_library),
    MedianLibrarySize = median(raw_library),
    MaxLibrarySize = max(raw_library),
    CountsFile = count_file,
    CountsBytes = file.info(count_file)$size,
    CountsSHA256 = digest::digest(count_file, algo = "sha256", file = TRUE, serialize = FALSE),
    CellsFile = cell_file,
    CellsSHA256 = digest::digest(cell_file, algo = "sha256", file = TRUE, serialize = FALSE),
    stringsAsFactors = FALSE
  )
  cat("GSE_MARKER_PARENT_WRITTEN", cell_type, nrow(frozen_type), "cells", nrow(counts_type), "genes", length(counts_type@x), "nnz\n")
  rm(counts_type)
  invisible(gc(FALSE))
}

# Step5. 写入审计清单
manifest <- rbindlist(manifest_parts)
state_audit <- rbindlist(state_audit_parts)
fwrite(manifest, manifest_file, sep = "\t")
fwrite(state_audit, file.path(output_dir, "GSE214695_marker_state_eligibility.tsv"), sep = "\t")

summary <- list(
  status = "PASS_GSE214695_MARKER_FULL_COUNTS",
  dataset = dataset,
  frozen_cells = nrow(cells),
  assigned_cells = sum(cells$Assigned),
  parent_types = nrow(manifest),
  genes = unique(manifest$Genes),
  eligible_state_dataset_pairs = sum(state_audit$EligibleGE5),
  normalization = "raw UMI counts; Seurat NormalizeData applies ln(CP10k+1) in the next step",
  gene_identity = "duplicate gene symbols summed; missing symbols retain UNMAPPED__GeneID and remain in library size",
  storage = "one gzip-compressed dgCMatrix RDS per dataset x parent type; no duplicate MTX or NPZ",
  recovery_model_changed = FALSE
)
jsonlite::write_json(summary, file.path(output_dir, "GSE214695_marker_counts_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
cat("GSE214695_MARKER_FULL_COUNTS_PASS", nrow(cells), "frozen cells", sum(cells$Assigned), "assigned cells", sum(state_audit$EligibleGE5), "eligible state-dataset pairs\n")
