#!/usr/bin/env Rscript
# 本脚本从 IBDverse H5AD 原始 counts 层构建冻结细胞的父类型级 EcoTyper CPM 输入。
# 逐细胞读取 CSR 行并保存父类型 checkpoint；原始 H5AD 只读且完整库大小保留。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞、表型和基因映射
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) stop("Usage: 07a_prepare_IBDverse_recovery_inputs.R RUN_ROOT [IBD|nonIBD]", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))
cohort <- if (length(args) == 2L) args[[2]] else "IBD"
if (!cohort %in% c("IBD", "nonIBD")) stop("Cohort must be IBD or nonIBD", call. = FALSE)

required_pkgs <- c("arrow", "data.table", "digest", "jsonlite", "Matrix", "rhdf5")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(Matrix)
library(rhdf5)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
gene_gate <- jsonlite::read_json(file.path(audit_dir, "scrna_gene_identity_resolution.json"), simplifyVector = TRUE)
if (!grepl("^PASS_SCRNA_GENE_IDENTITY_RESOLUTION", gene_gate$status)) stop("G4 gene identity gate did not pass", call. = FALSE)

cohort_audit_dir <- if (cohort == "IBD") audit_dir else file.path(output_root, "18_nonIBD_scrna_input_audit")
selected_file <- if (cohort == "IBD") "recovery_selected_cells.tsv.gz" else "nonIBD_recovery_selected_cells.tsv.gz"
phenotype_file <- if (cohort == "IBD") "donor_sample_phenotype_freeze.tsv" else "nonIBD_donor_sample_phenotype_freeze.tsv"
pool_file <- if (cohort == "IBD") file.path(eligible_pool_root, "IBDverse_E_MTAB_16999_eligible_ibd_cells.parquet") else file.path(cohort_audit_dir, "eligible_pools/IBDverse_E_MTAB_16999_eligible_nonIBD_cells.parquet")
selected_cells <- fread(file.path(cohort_audit_dir, selected_file), colClasses = "character")
selected_cells <- selected_cells[Dataset == "IBDverse_E_MTAB_16999"]
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]
eligible <- as.data.table(arrow::read_parquet(pool_file, as_data_frame = TRUE))
selected_cells <- merge(selected_cells, unique(eligible[, .(CellID, SourceRow)]), by = "CellID", all.x = TRUE, sort = FALSE)
selected_cells[, SourceRow := as.integer(SourceRow)]
if (anyNA(selected_cells$SourceRow) || any(selected_cells$SourceRow < 0L) || anyDuplicated(selected_cells$CellID) || anyDuplicated(selected_cells$SourceRow)) stop("IBDverse frozen cells do not map uniquely to H5AD rows", call. = FALSE)

sample_metadata <- fread(file.path(cohort_audit_dir, phenotype_file), colClasses = "character")
sample_metadata <- sample_metadata[Dataset == "IBDverse_E_MTAB_16999"]
gene_mapping <- fread(file.path(audit_dir, "scrna_resolved_model_gene_mapping.tsv.gz"), colClasses = "character")
gene_mapping <- gene_mapping[Dataset == "IBDverse_E_MTAB_16999"]
rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
cell_types <- sort(unique(selected_cells$CellType))
if (length(cell_types) != 16L) stop("IBDverse parent type count changed", call. = FALSE)
if (anyNA(rank_lookup[cell_types])) stop("Frozen rank is missing for IBDverse parent types", call. = FALSE)
if (anyDuplicated(sample_metadata[, .(DonorID, BiologicalSampleID)])) stop("IBDverse frozen phenotype contains duplicate samples", call. = FALSE)

model_genes_by_type <- setNames(vector("list", length(cell_types)), cell_types)
for (cell_type in cell_types) {
  W <- read.delim(file.path(state_discovery_root, cell_type, rank_lookup[[cell_type]], "W.txt"), check.names = FALSE)
  model_genes <- sub("__pos$", "", rownames(W)[grepl("__pos$", rownames(W))])
  if (!length(model_genes) || anyDuplicated(model_genes)) stop("Invalid frozen W genes: ", cell_type, call. = FALSE)
  mapping_type <- gene_mapping[CellType == cell_type]
  if (!nrow(mapping_type) || anyDuplicated(mapping_type$SourceKey)) stop("Invalid resolved source genes: ", cell_type, call. = FALSE)
  if (!all(mapping_type$ModelGene %in% model_genes)) stop("Resolved genes are absent from frozen W: ", cell_type, call. = FALSE)
  model_genes_by_type[[cell_type]] <- model_genes
}

source_file <- "/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad"
source_info_before <- file.info(source_file)[, c("size", "mtime")]
source_gene_id <- as.character(h5read(source_file, "/var/_index"))
source_categories <- as.character(h5read(source_file, "/var/gene_symbols/categories"))
source_codes <- as.integer(h5read(source_file, "/var/gene_symbols/codes"))
source_gene <- rep(NA_character_, length(source_codes))
valid_code <- source_codes >= 0L & source_codes < length(source_categories)
source_gene[valid_code] <- source_categories[source_codes[valid_code] + 1L]
source_key <- paste(source_gene_id, source_gene, sep = "|")
if (anyNA(source_gene) || anyDuplicated(source_key)) stop("IBDverse source feature identity changed", call. = FALSE)

source_cell_id <- as.character(h5read(source_file, "/obs/_index"))
if (max(selected_cells$SourceRow) + 1L > length(source_cell_id)) stop("IBDverse frozen source row exceeds H5AD dimensions", call. = FALSE)
if (!identical(source_cell_id[selected_cells$SourceRow + 1L], selected_cells$CellID)) stop("IBDverse frozen SourceRow no longer matches CellID", call. = FALSE)
pointer <- as.numeric(h5read(source_file, "/layers/counts/indptr", bit64conversion = "double"))
count_attributes <- h5readAttributes(source_file, "/layers/counts")
count_shape <- as.numeric(count_attributes$shape)
if (!identical(count_attributes[["encoding-type"]], "csr_matrix") || length(count_shape) != 2L || any(count_shape != c(length(source_cell_id), length(source_gene)))) stop("IBDverse counts layer is no longer cell-by-gene CSR", call. = FALSE)

prepare_root <- if (cohort == "IBD") "prepared_scRNA" else "prepared_scRNA_nonIBD"
prepare_dir <- file.path(RUN_ROOT, prepare_root, "IBDverse_E_MTAB_16999")
checkpoint_dir <- file.path(prepare_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(prepare_dir, "IBDverse_recovery_input_manifest.tsv"))) stop("IBDverse prepared manifest already exists", call. = FALSE)

file_handle <- H5Fopen(source_file, flags = "H5F_ACC_RDONLY")
data_handle <- H5Dopen(file_handle, "/layers/counts/data")
index_handle <- H5Dopen(file_handle, "/layers/counts/indices")

# Step2. 逐父类型读取冻结细胞的 CSR 行
manifest_parts <- list()
for (cell_type in cell_types) {
  frozen_type <- selected_cells[CellType == cell_type][order(SelectionOrder)]
  checkpoint_file <- file.path(checkpoint_dir, paste0(cell_type, ".rds"))

  if (file.exists(checkpoint_file)) {
    type_result <- readRDS(checkpoint_file)
    cat("IBDVERSE_CHECKPOINT_REUSED", cell_type, "\n")
  } else {
    model_genes <- model_genes_by_type[[cell_type]]
    mapping_type <- gene_mapping[CellType == cell_type & SourceKey %in% source_key]
    source_row <- match(mapping_type$SourceKey, source_key)
    source_to_model <- integer(length(source_gene))
    source_to_model[source_row] <- match(mapping_type$ModelGene, model_genes)
    triplet_i <- vector("list", nrow(frozen_type))
    triplet_j <- vector("list", nrow(frozen_type))
    triplet_x <- vector("list", nrow(frozen_type))
    full_library <- numeric(nrow(frozen_type))

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
      full_library[[cell_index]] <- sum(values)
      model_row <- source_to_model[source_index]
      keep_model <- model_row > 0L & values > 0
      triplet_i[[cell_index]] <- model_row[keep_model]
      triplet_j[[cell_index]] <- rep.int(cell_index, sum(keep_model))
      triplet_x[[cell_index]] <- values[keep_model]
      if (cell_index %% 250L == 0L || cell_index == nrow(frozen_type)) cat("IBDVERSE_CSR_ROWS", cell_type, cell_index, "/", nrow(frozen_type), "\n")
    }
    if (any(!is.finite(full_library)) || any(full_library <= 0)) stop("Invalid IBDverse full library size: ", cell_type, call. = FALSE)

    model_counts <- sparseMatrix(
      i = unlist(triplet_i, use.names = FALSE),
      j = unlist(triplet_j, use.names = FALSE),
      x = unlist(triplet_x, use.names = FALSE),
      dims = c(length(model_genes), nrow(frozen_type)),
      dimnames = list(model_genes, frozen_type$CellID)
    )
    type_result <- list(counts = model_counts, ids = frozen_type$CellID, libraries = full_library)
    saveRDS(type_result, checkpoint_file, compress = "gzip")
    cat("IBDVERSE_CHECKPOINT_SAVED", cell_type, nrow(frozen_type), "cells\n")
    rm(model_counts, triplet_i, triplet_j, triplet_x, full_library)
    invisible(gc(FALSE))
  }

  if (!identical(type_result$ids, frozen_type$CellID)) stop("IBDverse checkpoint cell order changed: ", cell_type, call. = FALSE)
  counts_type <- type_result$counts
  library_type <- type_result$libraries
  cpm_type <- counts_type %*% Diagonal(x = 1e6 / library_type)
  type_dir <- file.path(prepare_dir, cell_type)
  dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)
  expression_file <- file.path(type_dir, "data.txt")
  annotation_file_out <- file.path(type_dir, "annotation.txt")
  expression_table <- data.table(Gene = rownames(cpm_type))
  expression_table <- cbind(expression_table, as.data.table(as.matrix(cpm_type)))
  setnames(expression_table, c("Gene", frozen_type$CellID))
  fwrite(expression_table, expression_file, sep = "\t", quote = FALSE)

  annotation <- merge(frozen_type[, .(ID = CellID, CellType, DonorID, BiologicalSampleID)], sample_metadata[, .(DonorID, BiologicalSampleID, Disease, Inflammation, Region3)], by = c("DonorID", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
  annotation <- annotation[match(frozen_type$CellID, ID)]
  if (anyNA(annotation$Disease) || anyNA(annotation$Inflammation) || anyNA(annotation$Region3)) stop("IBDverse frozen phenotype join failed: ", cell_type, call. = FALSE)
  fwrite(annotation, annotation_file_out, sep = "\t", quote = FALSE)

  manifest_parts[[cell_type]] <- data.frame(
    Dataset = "IBDverse_E_MTAB_16999", CellType = cell_type, FrozenRank = rank_lookup[[cell_type]], Cells = nrow(frozen_type),
    ModelGenes = nrow(cpm_type), NonzeroModelGenes = sum(Matrix::rowSums(counts_type) > 0),
    MinLibrarySize = min(library_type), MedianLibrarySize = median(library_type), MaxLibrarySize = max(library_type),
    ExpressionFile = expression_file, ExpressionSHA256 = digest::digest(expression_file, algo = "sha256", file = TRUE, serialize = FALSE),
    AnnotationFile = annotation_file_out, AnnotationSHA256 = digest::digest(annotation_file_out, algo = "sha256", file = TRUE, serialize = FALSE)
  )
  cat("IBDVERSE_PARENT_INPUT_WRITTEN", cell_type, nrow(frozen_type), "cells", nrow(cpm_type), "model genes\n")
  rm(type_result, counts_type, cpm_type, expression_table)
  invisible(gc(FALSE))
}

H5Dclose(index_handle)
H5Dclose(data_handle)
H5Fclose(file_handle)

# Step3. 保存输入清单和源文件审计
source_info_after <- file.info(source_file)[, c("size", "mtime")]
if (!identical(source_info_before, source_info_after)) stop("IBDverse source changed during read", call. = FALSE)
manifest <- rbindlist(manifest_parts)
fwrite(manifest, file.path(prepare_dir, "IBDverse_recovery_input_manifest.tsv"), sep = "\t")
summary <- list(
  status = if (cohort == "IBD") "PASS_IBDVERSE_RECOVERY_INPUT_PREPARATION" else "PASS_IBDVERSE_NONIBD_RECOVERY_INPUT_PREPARATION", dataset = "IBDverse_E_MTAB_16999", cohort = cohort, parent_types = nrow(manifest), cells = sum(manifest$Cells), source_file = source_file,
  source_bytes = source_info_after$size, normalization = "raw counts collapsed by frozen model gene, full-feature raw UMI library CPM; EcoTyper applies log2(CPM+1) and row scaling",
  duplicate_source_rows = "summed before CPM", model_changed = FALSE, recovery_started = FALSE
)
jsonlite::write_json(summary, file.path(prepare_dir, "IBDverse_recovery_input_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(prepare_dir, "R_sessionInfo.txt"))
cat("IBDVERSE_RECOVERY_INPUT_PREPARATION_PASS", sum(manifest$Cells), "cells", nrow(manifest), "parent types\n")
