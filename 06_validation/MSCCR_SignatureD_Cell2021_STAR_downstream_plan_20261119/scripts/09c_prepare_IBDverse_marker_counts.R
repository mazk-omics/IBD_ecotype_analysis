#!/usr/bin/env Rscript
# 本脚本从 IBDverse H5AD 原始 counts 层构建 marker validation 的全转录组稀疏输入。
# 逐父类型读取冻结细胞 CSR 行；每个父类型 RDS 同时作为可复用 checkpoint。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞和正式 assignment
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 09c_prepare_IBDverse_marker_counts.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("arrow", "data.table", "digest", "jsonlite", "Matrix", "rhdf5")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(Matrix)
library(rhdf5)

dataset <- "IBDverse_E_MTAB_16999"
audit_dir <- file.path(output_root, "03_scrna_input_audit")
selected_cells_raw <- fread(file.path(audit_dir, "recovery_selected_cells.tsv.gz"), colClasses = "character")
selected_cells <- selected_cells_raw[Dataset == dataset]
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]

eligible <- as.data.table(arrow::read_parquet(file.path(eligible_pool_root, "IBDverse_E_MTAB_16999_eligible_ibd_cells.parquet"), as_data_frame = TRUE))
selected_cells <- merge(selected_cells, unique(eligible[, .(CellID, SourceRow)]), by = "CellID", all.x = TRUE, sort = FALSE)
selected_cells <- selected_cells[match(selected_cells_raw[Dataset == dataset]$CellID, CellID)]
selected_cells[, SourceRow := as.integer(SourceRow)]
if (nrow(selected_cells) != 35899L || uniqueN(selected_cells$CellType) != 16L || anyNA(selected_cells$SourceRow) || any(selected_cells$SourceRow < 0L) || anyDuplicated(selected_cells$CellID) || anyDuplicated(selected_cells$SourceRow)) stop("IBDverse frozen source identity changed", call. = FALSE)

sample_metadata <- fread(file.path(audit_dir, "donor_sample_phenotype_freeze.tsv"), colClasses = "character")
sample_metadata <- sample_metadata[Dataset == dataset]
if (anyDuplicated(sample_metadata[, .(DonorID, BiologicalSampleID)])) stop("IBDverse frozen phenotype identity changed", call. = FALSE)

job_file <- file.path(output_root, "07c_IBDverse_official_recovery/IBDverse_E_MTAB_16999_official_recovery_jobs_completed.tsv")
jobs <- fread(job_file, colClasses = "character")
if (nrow(jobs) != 16L || any(jobs$Dataset != dataset) || any(jobs$ExitStatus != "0") || anyDuplicated(jobs$CellType)) stop("IBDverse completed recovery manifest changed", call. = FALSE)

assignment_parts <- list()
for (cell_type in jobs$CellType) {
  assignment_file <- file.path(jobs[CellType == cell_type]$ResultDir, "state_assignment.txt")
  assignment <- fread(assignment_file, colClasses = "character")
  if (!identical(names(assignment), c("ID", "State", "InitialState")) || anyDuplicated(assignment$ID)) stop("Invalid IBDverse state assignment: ", cell_type, call. = FALSE)
  if (!all(assignment$ID %in% selected_cells[CellType == cell_type]$CellID)) stop("IBDverse assignment contains non-frozen cells: ", cell_type, call. = FALSE)
  assignment_parts[[cell_type]] <- assignment
}
assignment_all <- rbindlist(assignment_parts)
if (anyDuplicated(assignment_all$ID)) stop("IBDverse assignment IDs overlap across parent types", call. = FALSE)

cells <- merge(selected_cells, sample_metadata[, .(DonorID, BiologicalSampleID, Disease, Inflammation, Region3)], by = c("DonorID", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
cells <- merge(cells, assignment_all[, .(CellID = ID, State, InitialState)], by = "CellID", all.x = TRUE, sort = FALSE)
cells <- cells[match(selected_cells$CellID, CellID)]
cells[is.na(State), State := "Unassigned"]
cells[is.na(InitialState), InitialState := "Unassigned"]
cells[, Assigned := State != "Unassigned"]
if (anyNA(cells$Disease) || anyNA(cells$Inflammation) || anyNA(cells$Region3) || !identical(cells$CellID, selected_cells$CellID)) stop("IBDverse cell metadata join failed", call. = FALSE)
print(cells[, .(FrozenCells = .N, AssignedCells = sum(Assigned), States = uniqueN(State[Assigned])), by = CellType][order(CellType)])

# Step2. 固定 H5AD counts 和基因身份
source_file <- "/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad"
source_info_before <- file.info(source_file)[, c("size", "mtime")]
source_gene_id <- as.character(h5read(source_file, "/var/_index"))
source_categories <- as.character(h5read(source_file, "/var/gene_symbols/categories"))
source_codes <- as.integer(h5read(source_file, "/var/gene_symbols/codes"))
source_gene <- rep(NA_character_, length(source_codes))
valid_code <- source_codes >= 0L & source_codes < length(source_categories)
source_gene[valid_code] <- source_categories[source_codes[valid_code] + 1L]
source_gene[is.na(source_gene) | !nzchar(source_gene)] <- paste0("UNMAPPED__", source_gene_id[is.na(source_gene) | !nzchar(source_gene)])
source_key <- paste(source_gene_id, source_gene, sep = "|")
if (anyDuplicated(source_key)) stop("IBDverse source feature identity changed", call. = FALSE)

gene_symbols <- unique(source_gene)
source_to_symbol <- match(source_gene, gene_symbols)
source_cell_id <- as.character(h5read(source_file, "/obs/_index"))
if (max(cells$SourceRow) + 1L > length(source_cell_id) || !identical(source_cell_id[cells$SourceRow + 1L], cells$CellID)) stop("IBDverse SourceRow no longer matches CellID", call. = FALSE)

pointer <- as.numeric(h5read(source_file, "/layers/counts/indptr", bit64conversion = "double"))
count_attributes <- h5readAttributes(source_file, "/layers/counts")
count_shape <- as.numeric(count_attributes$shape)
if (!identical(count_attributes[["encoding-type"]], "csr_matrix") || length(count_shape) != 2L || any(count_shape != c(length(source_cell_id), length(source_gene)))) stop("IBDverse counts layer is no longer cell-by-gene CSR", call. = FALSE)

output_dir <- file.path(output_root, "09_marker_recovery/full_counts", dataset)
manifest_file <- file.path(output_dir, "IBDverse_E_MTAB_16999_marker_counts_manifest.tsv")
if (file.exists(manifest_file)) stop("IBDverse marker counts already exist; do not overwrite", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

file_handle <- H5Fopen(source_file, flags = "H5F_ACC_RDONLY")
data_handle <- H5Dopen(file_handle, "/layers/counts/data")
index_handle <- H5Dopen(file_handle, "/layers/counts/indices")

# Step3. 逐父类型读取 CSR 行并保存 checkpoint
manifest_parts <- list()
state_audit_parts <- list()
for (cell_type in sort(unique(cells$CellType))) {
  frozen_type <- cells[CellType == cell_type][order(SelectionOrder)]
  type_dir <- file.path(output_dir, cell_type)
  count_file <- file.path(type_dir, "counts.rds")
  cell_file <- file.path(type_dir, "cells.tsv.gz")
  library_file <- file.path(type_dir, "library_size.tsv.gz")

  if (file.exists(count_file) && file.exists(cell_file) && file.exists(library_file)) {
    counts_type <- readRDS(count_file)
    checkpoint_cells <- fread(cell_file, colClasses = "character")
    if (!methods::is(counts_type, "dgCMatrix") || !identical(colnames(counts_type), frozen_type$CellID) || !identical(checkpoint_cells$CellID, frozen_type$CellID)) stop("IBDverse checkpoint identity changed: ", cell_type, call. = FALSE)
    library_type <- as.numeric(Matrix::colSums(counts_type))
    cat("IBDVERSE_MARKER_CHECKPOINT_REUSED", cell_type, ncol(counts_type), "cells\n")
  } else {
    triplet_i <- vector("list", nrow(frozen_type))
    triplet_j <- vector("list", nrow(frozen_type))
    triplet_x <- vector("list", nrow(frozen_type))
    library_type <- numeric(nrow(frozen_type))

    for (cell_index in seq_len(nrow(frozen_type))) {
      source_cell_row <- frozen_type$SourceRow[[cell_index]]
      start <- pointer[[source_cell_row + 1L]]
      stop <- pointer[[source_cell_row + 2L]]
      count <- stop - start
      if (!is.finite(count) || count <= 0 || count > 10000000) stop("Invalid IBDverse CSR pointer: ", frozen_type$CellID[[cell_index]], call. = FALSE)

      data_space <- H5Dget_space(data_handle)
      H5Sselect_hyperslab(data_space, start = start + 1, count = 1, block = count)
      memory_space <- H5Screate_simple(count)
      values <- H5Dread(data_handle, h5spaceFile = data_space, h5spaceMem = memory_space, drop = TRUE)
      H5Sclose(memory_space)
      H5Sclose(data_space)

      index_space <- H5Dget_space(index_handle)
      H5Sselect_hyperslab(index_space, start = start + 1, count = 1, block = count)
      memory_space <- H5Screate_simple(count)
      source_index <- as.integer(H5Dread(index_handle, h5spaceFile = index_space, h5spaceMem = memory_space, drop = TRUE)) + 1L
      H5Sclose(memory_space)
      H5Sclose(index_space)

      if (any(!is.finite(values)) || any(values < 0) || any(abs(values - round(values)) > 1e-5) || any(source_index < 1L | source_index > length(source_gene))) stop("Invalid IBDverse raw counts: ", frozen_type$CellID[[cell_index]], call. = FALSE)
      keep <- values > 0
      library_type[[cell_index]] <- sum(values)
      triplet_i[[cell_index]] <- source_to_symbol[source_index[keep]]
      triplet_j[[cell_index]] <- rep.int(cell_index, sum(keep))
      triplet_x[[cell_index]] <- values[keep]
      if (cell_index %% 250L == 0L || cell_index == nrow(frozen_type)) cat("IBDVERSE_MARKER_CSR_ROWS", cell_type, cell_index, "/", nrow(frozen_type), "\n")
    }
    if (any(!is.finite(library_type)) || any(library_type <= 0)) stop("Invalid IBDverse library size: ", cell_type, call. = FALSE)

    counts_type <- sparseMatrix(
      i = unlist(triplet_i, use.names = FALSE), j = unlist(triplet_j, use.names = FALSE), x = unlist(triplet_x, use.names = FALSE),
      dims = c(length(gene_symbols), nrow(frozen_type)), dimnames = list(gene_symbols, frozen_type$CellID)
    )
    counts_type <- methods::as(counts_type, "dgCMatrix")
    if (!isTRUE(all.equal(as.numeric(Matrix::colSums(counts_type)), library_type, tolerance = 0))) stop("IBDverse gene-symbol collapse changed library size: ", cell_type, call. = FALSE)

    dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(counts_type, count_file, compress = "gzip")
    fwrite(frozen_type, cell_file, sep = "\t", quote = FALSE)
    fwrite(data.table(CellID = frozen_type$CellID, RawLibrarySize = library_type, ExportedLibrarySize = Matrix::colSums(counts_type)), library_file, sep = "\t", quote = FALSE)
    cat("IBDVERSE_MARKER_PARENT_WRITTEN", cell_type, nrow(frozen_type), "cells", nrow(counts_type), "genes", length(counts_type@x), "nnz\n")
    rm(triplet_i, triplet_j, triplet_x)
  }

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
  rm(counts_type)
  invisible(gc(FALSE))
}

H5Dclose(index_handle)
H5Dclose(data_handle)
H5Fclose(file_handle)

# Step4. 写入审计清单
source_info_after <- file.info(source_file)[, c("size", "mtime")]
if (!identical(source_info_before, source_info_after)) stop("IBDverse source changed during read", call. = FALSE)
manifest <- rbindlist(manifest_parts)
state_audit <- rbindlist(state_audit_parts)
fwrite(manifest, manifest_file, sep = "\t")
fwrite(state_audit, file.path(output_dir, "IBDverse_E_MTAB_16999_marker_state_eligibility.tsv"), sep = "\t")

summary <- list(
  status = "PASS_IBDVERSE_MARKER_FULL_COUNTS", dataset = dataset, frozen_cells = nrow(cells), assigned_cells = sum(cells$Assigned), parent_types = nrow(manifest),
  genes = unique(manifest$Genes), eligible_state_dataset_pairs = sum(state_audit$EligibleGE5), source_file = source_file, source_bytes = source_info_after$size,
  normalization = "raw UMI counts; Seurat NormalizeData applies ln(CP10k+1) in the next step",
  gene_identity = "duplicate gene symbols summed; missing symbols retain UNMAPPED__GeneID and remain in library size",
  storage = "one gzip-compressed dgCMatrix RDS per dataset x parent type; each final RDS is a restart checkpoint", recovery_model_changed = FALSE
)
jsonlite::write_json(summary, file.path(output_dir, "IBDverse_E_MTAB_16999_marker_counts_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
cat("IBDVERSE_MARKER_FULL_COUNTS_PASS", nrow(cells), "frozen cells", sum(cells$Assigned), "assigned cells", sum(state_audit$EligibleGE5), "eligible state-dataset pairs\n")
