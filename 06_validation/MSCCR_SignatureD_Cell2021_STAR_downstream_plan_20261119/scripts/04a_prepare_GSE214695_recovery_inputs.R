#!/usr/bin/env Rscript
# 本脚本从 GSE214695 原始 10x counts 构建冻结细胞的父类型级 EcoTyper CPM 输入。
# 每个原始样本只读一次；重复源基因先在 counts 层求和，再按完整库大小计算 CPM。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞、表型和基因映射
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) stop("Usage: 04a_prepare_GSE214695_recovery_inputs.R RUN_ROOT [IBD|nonIBD]", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))
cohort <- if (length(args) == 2L) args[[2]] else "IBD"
if (!cohort %in% c("IBD", "nonIBD")) stop("Cohort must be IBD or nonIBD", call. = FALSE)

required_pkgs <- c("data.table", "digest", "jsonlite", "Matrix")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(Matrix)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
gene_gate <- jsonlite::read_json(file.path(audit_dir, "scrna_gene_identity_resolution.json"), simplifyVector = TRUE)
if (!grepl("^PASS_SCRNA_GENE_IDENTITY_RESOLUTION", gene_gate$status)) stop("G4 gene identity gate did not pass", call. = FALSE)

cohort_audit_dir <- if (cohort == "IBD") audit_dir else file.path(output_root, "18_nonIBD_scrna_input_audit")
selected_file <- if (cohort == "IBD") "recovery_selected_cells.tsv.gz" else "nonIBD_recovery_selected_cells.tsv.gz"
phenotype_file <- if (cohort == "IBD") "donor_sample_phenotype_freeze.tsv" else "nonIBD_donor_sample_phenotype_freeze.tsv"
selected_cells <- fread(file.path(cohort_audit_dir, selected_file), colClasses = "character")
selected_cells <- selected_cells[Dataset == "GSE214695"]
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]
sample_metadata <- fread(file.path(cohort_audit_dir, phenotype_file), colClasses = "character")
sample_metadata <- sample_metadata[Dataset == "GSE214695"]
gene_mapping <- fread(file.path(audit_dir, "scrna_resolved_model_gene_mapping.tsv.gz"), colClasses = "character")
gene_mapping <- gene_mapping[Dataset == "GSE214695"]

rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
cell_types <- sort(unique(selected_cells$CellType))
if (length(cell_types) != 18L) stop("GSE214695 parent type count changed", call. = FALSE)
if (anyNA(rank_lookup[cell_types])) stop("Frozen rank is missing for GSE214695 parent types", call. = FALSE)
if (anyDuplicated(sample_metadata[, .(DonorID, BiologicalSampleID)])) stop("GSE214695 frozen phenotype contains duplicate samples", call. = FALSE)

source_root <- file.path(project_root, "06_validation/GSE214695_scRNA_recovery_20261009/source")
matrix_files <- sort(list.files(source_root, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE))
sample_labels <- sub("^GSM[0-9]+_(.*)_matrix\\.mtx\\.gz$", "\\1", basename(matrix_files))
sample_ids <- gsub("-", "", sample_labels)
keep_ibd <- sample_ids %in% unique(selected_cells$BiologicalSampleID)
matrix_files <- matrix_files[keep_ibd]
sample_ids <- sample_ids[keep_ibd]
if (!setequal(sample_ids, unique(selected_cells$BiologicalSampleID))) stop("GSE214695 raw sample IDs differ from frozen IBD samples", call. = FALSE)
if (anyDuplicated(sample_ids)) stop("GSE214695 raw IBD sample IDs are duplicated", call. = FALSE)

model_genes_by_type <- setNames(vector("list", length(cell_types)), cell_types)
for (cell_type in cell_types) {
  W_file <- file.path(state_discovery_root, cell_type, rank_lookup[[cell_type]], "W.txt")
  W <- read.delim(W_file, check.names = FALSE)
  W_rows <- rownames(W)
  model_genes <- sub("__pos$", "", W_rows[grepl("__pos$", W_rows)])
  if (!length(model_genes) || anyDuplicated(model_genes)) stop("Invalid frozen W genes: ", cell_type, call. = FALSE)
  mapping_type <- gene_mapping[CellType == cell_type]
  if (!nrow(mapping_type) || anyDuplicated(mapping_type$SourceKey)) stop("Invalid resolved source genes: ", cell_type, call. = FALSE)
  if (!all(mapping_type$ModelGene %in% model_genes)) stop("Resolved genes are absent from frozen W: ", cell_type, call. = FALSE)
  model_genes_by_type[[cell_type]] <- model_genes
}

prepare_root <- if (cohort == "IBD") "prepared_scRNA" else "prepared_scRNA_nonIBD"
prepare_dir <- file.path(RUN_ROOT, prepare_root, "GSE214695")
dir.create(prepare_dir, recursive = TRUE, showWarnings = FALSE)
matrix_parts <- setNames(vector("list", length(cell_types)), cell_types)
id_parts <- setNames(vector("list", length(cell_types)), cell_types)
library_parts <- setNames(vector("list", length(cell_types)), cell_types)

# Step2. 逐原始样本抽取冻结细胞并保留完整库大小
for (sample_index in seq_along(matrix_files)) {
  matrix_file <- matrix_files[[sample_index]]
  sample_id <- sample_ids[[sample_index]]
  file_prefix <- sub("_matrix\\.mtx\\.gz$", "", matrix_file)
  feature_file <- paste0(file_prefix, "_features.tsv.gz")
  barcode_file <- paste0(file_prefix, "_barcodes.tsv.gz")

  feature <- fread(feature_file, header = FALSE, select = 1:2, col.names = c("GeneID", "Gene"))
  barcode <- fread(barcode_file, header = FALSE, select = 1, col.names = "Barcode")
  feature[, SourceKey := paste(GeneID, Gene, sep = "|")]
  if (anyDuplicated(feature$SourceKey)) stop("Raw feature keys are duplicated: ", sample_id, call. = FALSE)
  sample_cells <- selected_cells[BiologicalSampleID == sample_id]
  sample_cells[, RawBarcode := sub("^.*\\|", "", CellID)]
  column_index <- match(sample_cells$RawBarcode, barcode$Barcode)
  if (anyNA(column_index) || anyDuplicated(column_index)) stop("Frozen barcodes do not map uniquely: ", sample_id, call. = FALSE)

  counts <- Matrix::readMM(gzfile(matrix_file))
  if (nrow(counts) != nrow(feature) || ncol(counts) != nrow(barcode)) stop("Raw matrix dimensions changed: ", sample_id, call. = FALSE)
  selected_counts <- counts[, column_index, drop = FALSE]
  library_size <- Matrix::colSums(selected_counts)
  if (any(!is.finite(library_size)) || any(library_size <= 0)) stop("Invalid library size: ", sample_id, call. = FALSE)

  for (cell_type in unique(sample_cells$CellType)) {
    type_cells <- sample_cells[CellType == cell_type]
    local_column <- match(type_cells$CellID, sample_cells$CellID)
    mapping_type <- gene_mapping[CellType == cell_type]
    source_row <- match(mapping_type$SourceKey, feature$SourceKey)
    if (anyNA(source_row) || anyDuplicated(mapping_type[, .(ModelGene, SourceKey)])) stop("Resolved gene rows do not map uniquely: ", cell_type, " ", sample_id, call. = FALSE)

    model_genes <- model_genes_by_type[[cell_type]]
    collapse_matrix <- sparseMatrix(
      i = match(mapping_type$ModelGene, model_genes),
      j = seq_len(nrow(mapping_type)),
      x = 1,
      dims = c(length(model_genes), nrow(mapping_type)),
      dimnames = list(model_genes, mapping_type$SourceKey)
    )
    model_counts <- collapse_matrix %*% selected_counts[source_row, local_column, drop = FALSE]
    colnames(model_counts) <- type_cells$CellID

    matrix_parts[[cell_type]][[sample_id]] <- model_counts
    id_parts[[cell_type]][[sample_id]] <- type_cells$CellID
    library_parts[[cell_type]][[sample_id]] <- library_size[local_column]
  }
  cat("GSE_SAMPLE_EXTRACTED", sample_id, nrow(sample_cells), "cells", nrow(feature), "features\n")
  rm(counts, selected_counts, barcode)
  invisible(gc(FALSE))
}

# Step3. 按冻结顺序计算 CPM 并写父类型级输入
manifest_parts <- list()
for (cell_type in cell_types) {
  frozen_type <- selected_cells[CellType == cell_type][order(SelectionOrder)]
  counts_type <- do.call(cbind, matrix_parts[[cell_type]])
  ids_type <- unlist(id_parts[[cell_type]], use.names = FALSE)
  library_type <- unlist(library_parts[[cell_type]], use.names = FALSE)
  frozen_order <- match(frozen_type$CellID, ids_type)
  if (anyNA(frozen_order) || anyDuplicated(ids_type) || ncol(counts_type) != nrow(frozen_type)) stop("Frozen cell assembly failed: ", cell_type, call. = FALSE)
  counts_type <- counts_type[, frozen_order, drop = FALSE]
  library_type <- library_type[frozen_order]
  colnames(counts_type) <- frozen_type$CellID
  cpm_type <- counts_type %*% Diagonal(x = 1e6 / library_type)

  type_dir <- file.path(prepare_dir, cell_type)
  dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)
  expression_file <- file.path(type_dir, "data.txt")
  annotation_file_out <- file.path(type_dir, "annotation.txt")

  expression_table <- data.table(Gene = rownames(cpm_type))
  expression_table <- cbind(expression_table, as.data.table(as.matrix(cpm_type)))
  setnames(expression_table, c("Gene", frozen_type$CellID))
  fwrite(expression_table, expression_file, sep = "\t", quote = FALSE)

  annotation <- merge(
    frozen_type[, .(ID = CellID, CellType, DonorID, BiologicalSampleID)],
    sample_metadata[, .(DonorID, BiologicalSampleID, Disease, Inflammation, Region3)],
    by = c("DonorID", "BiologicalSampleID"),
    all.x = TRUE,
    sort = FALSE
  )
  annotation <- annotation[match(frozen_type$CellID, ID)]
  if (anyNA(annotation$Disease) || anyNA(annotation$Inflammation) || anyNA(annotation$Region3)) stop("Frozen phenotype join failed: ", cell_type, call. = FALSE)
  fwrite(annotation, annotation_file_out, sep = "\t", quote = FALSE)

  manifest_parts[[cell_type]] <- data.frame(
    Dataset = "GSE214695",
    CellType = cell_type,
    FrozenRank = rank_lookup[[cell_type]],
    Cells = nrow(frozen_type),
    ModelGenes = nrow(cpm_type),
    NonzeroModelGenes = sum(Matrix::rowSums(counts_type) > 0),
    MinLibrarySize = min(library_type),
    MedianLibrarySize = median(library_type),
    MaxLibrarySize = max(library_type),
    ExpressionFile = expression_file,
    ExpressionSHA256 = digest::digest(expression_file, algo = "sha256", file = TRUE, serialize = FALSE),
    AnnotationFile = annotation_file_out,
    AnnotationSHA256 = digest::digest(annotation_file_out, algo = "sha256", file = TRUE, serialize = FALSE),
    stringsAsFactors = FALSE
  )
  cat("GSE_PARENT_INPUT_WRITTEN", cell_type, nrow(frozen_type), "cells", nrow(cpm_type), "model genes\n")
}

# Step4. 保存输入清单
manifest <- rbindlist(manifest_parts)
data.table::fwrite(manifest, file.path(prepare_dir, "GSE214695_recovery_input_manifest.tsv"), sep = "\t")
summary <- list(
  status = if (cohort == "IBD") "PASS_GSE214695_RECOVERY_INPUT_PREPARATION" else "PASS_GSE214695_NONIBD_RECOVERY_INPUT_PREPARATION",
  dataset = "GSE214695",
  cohort = cohort,
  parent_types = nrow(manifest),
  cells = sum(manifest$Cells),
  normalization = "raw counts collapsed by frozen model gene, full-library CPM; EcoTyper applies log2(CPM+1) and row scaling",
  duplicate_source_rows = "summed before CPM",
  model_changed = FALSE,
  recovery_started = FALSE
)
jsonlite::write_json(summary, file.path(prepare_dir, "GSE214695_recovery_input_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(prepare_dir, "R_sessionInfo.txt"))
cat("GSE214695_RECOVERY_INPUT_PREPARATION_PASS", sum(manifest$Cells), "cells", nrow(manifest), "parent types\n")
