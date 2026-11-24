#!/usr/bin/env Rscript
# 本脚本从 SCP259 三个原始 UMI bundle 构建冻结细胞的父类型级 EcoTyper CPM 输入。
# 原始矩阵只读；重复源基因在 counts 层求和，CPM 分母使用完整 feature 库。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞、表型和基因映射
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) stop("Usage: 06a_prepare_SCP259_recovery_inputs.R RUN_ROOT [IBD|nonIBD]", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))
cohort <- if (length(args) == 2L) args[[2]] else "IBD"
if (!cohort %in% c("IBD", "nonIBD")) stop("Cohort must be IBD or nonIBD", call. = FALSE)

required_pkgs <- c("arrow", "data.table", "digest", "jsonlite", "Matrix")
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
pool_file <- if (cohort == "IBD") file.path(eligible_pool_root, "SCP259_eligible_ibd_cells.parquet") else file.path(cohort_audit_dir, "eligible_pools/SCP259_eligible_nonIBD_cells.parquet")
selected_cells <- fread(file.path(cohort_audit_dir, selected_file), colClasses = "character")
selected_cells <- selected_cells[Dataset == "SCP259"]
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]
eligible <- as.data.table(arrow::read_parquet(pool_file, as_data_frame = TRUE))
selected_cells <- merge(selected_cells, unique(eligible[, .(CellID, SourceBundle)]), by = "CellID", all.x = TRUE, sort = FALSE)
if (anyNA(selected_cells$SourceBundle) || anyDuplicated(selected_cells$CellID)) stop("SCP259 frozen cells do not map uniquely to bundles", call. = FALSE)

sample_metadata <- fread(file.path(cohort_audit_dir, phenotype_file), colClasses = "character")
sample_metadata <- sample_metadata[Dataset == "SCP259"]
gene_mapping <- fread(file.path(audit_dir, "scrna_resolved_model_gene_mapping.tsv.gz"), colClasses = "character")
gene_mapping <- gene_mapping[Dataset == "SCP259"]
rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
cell_types <- sort(unique(selected_cells$CellType))
if (length(cell_types) != 18L) stop("SCP259 parent type count changed", call. = FALSE)
if (anyNA(rank_lookup[cell_types])) stop("Frozen rank is missing for SCP259 parent types", call. = FALSE)
if (anyDuplicated(sample_metadata[, .(DonorID, BiologicalSampleID)])) stop("SCP259 frozen phenotype contains duplicate samples", call. = FALSE)

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

source_root <- file.path(project_root, "01_raw_data/scRNA/SCP259/raw/SCP259/expression")
bundle_config <- data.table(
  SourceBundle = c("Epithelial", "Stromal", "Immune"),
  Subdir = c("5cdc540d328cee7a2efc2348", "5cdc540d328cee7a2efc2349", "5cdc540d328cee7a2efc234a"),
  Prefix = c("Epi", "Fib", "Imm")
)
prepare_root <- if (cohort == "IBD") "prepared_scRNA" else "prepared_scRNA_nonIBD"
prepare_dir <- file.path(RUN_ROOT, prepare_root, "SCP259")
checkpoint_dir <- file.path(prepare_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(prepare_dir, "SCP259_recovery_input_manifest.tsv"))) stop("SCP259 prepared manifest already exists", call. = FALSE)

counts_parts <- setNames(vector("list", length(cell_types)), cell_types)
id_parts <- setNames(vector("list", length(cell_types)), cell_types)
library_parts <- setNames(vector("list", length(cell_types)), cell_types)
bundle_audit <- list()

# Step2. 逐 bundle 提取冻结细胞和完整库大小
for (bundle_index in seq_len(nrow(bundle_config))) {
  source_bundle <- bundle_config$SourceBundle[[bundle_index]]
  prefix <- bundle_config$Prefix[[bundle_index]]
  bundle_dir <- file.path(source_root, bundle_config$Subdir[[bundle_index]])
  feature_file <- file.path(bundle_dir, paste0(prefix, ".genes.tsv"))
  barcode_file <- file.path(bundle_dir, paste0(prefix, ".barcodes2.tsv"))
  matrix_file <- file.path(bundle_dir, paste0("gene_sorted-", prefix, ".matrix.mtx"))
  checkpoint_file <- file.path(checkpoint_dir, paste0(source_bundle, ".rds"))

  if (file.exists(checkpoint_file)) {
    bundle_result <- readRDS(checkpoint_file)
    cat("SCP259_CHECKPOINT_REUSED", source_bundle, "\n")
  } else {
    selected_bundle <- selected_cells[SourceBundle == source_bundle]
    features <- fread(feature_file, header = FALSE, select = 1L, col.names = "Gene")
    barcodes <- fread(barcode_file, header = FALSE, select = 1L, col.names = "CellID")
    features[, SourceKey := paste(Gene, Gene, sep = "|")]
    if (anyDuplicated(features$SourceKey) || anyDuplicated(barcodes$CellID)) stop("Duplicated SCP259 source identity: ", source_bundle, call. = FALSE)
    raw_cell_index <- match(selected_bundle$CellID, barcodes$CellID)
    if (anyNA(raw_cell_index) || anyDuplicated(raw_cell_index)) stop("SCP259 frozen cells do not map uniquely: ", source_bundle, call. = FALSE)

    raw <- Matrix::readMM(matrix_file)
    if (nrow(raw) != nrow(features) || ncol(raw) != nrow(barcodes)) stop("SCP259 raw dimensions changed: ", source_bundle, call. = FALSE)
    if (any(!is.finite(raw@x)) || any(raw@x < 0) || any(abs(raw@x - round(raw@x)) > 1e-8)) stop("SCP259 raw matrix is not nonnegative integer counts: ", source_bundle, call. = FALSE)
    selected_raw <- raw[, raw_cell_index, drop = FALSE]
    full_library <- Matrix::colSums(selected_raw)
    if (any(!is.finite(full_library)) || any(full_library <= 0)) stop("Invalid SCP259 full library size: ", source_bundle, call. = FALSE)

    bundle_counts <- list()
    bundle_ids <- list()
    bundle_libraries <- list()
    bundle_cell_types <- sort(unique(selected_bundle$CellType))
    for (cell_type in bundle_cell_types) {
      type_cells <- selected_bundle[CellType == cell_type]
      local_column <- match(type_cells$CellID, selected_bundle$CellID)
      mapping_type <- gene_mapping[CellType == cell_type & SourceKey %in% features$SourceKey]
      source_row <- match(mapping_type$SourceKey, features$SourceKey)
      if (!nrow(mapping_type) || anyNA(source_row)) stop("No SCP259 mapped genes are available: ", cell_type, " ", source_bundle, call. = FALSE)
      collapse_matrix <- sparseMatrix(
        i = match(mapping_type$ModelGene, model_genes_by_type[[cell_type]]),
        j = seq_len(nrow(mapping_type)),
        x = 1,
        dims = c(length(model_genes_by_type[[cell_type]]), nrow(mapping_type)),
        dimnames = list(model_genes_by_type[[cell_type]], mapping_type$SourceKey)
      )
      model_counts <- collapse_matrix %*% selected_raw[source_row, local_column, drop = FALSE]
      colnames(model_counts) <- type_cells$CellID
      bundle_counts[[cell_type]] <- model_counts
      bundle_ids[[cell_type]] <- type_cells$CellID
      bundle_libraries[[cell_type]] <- full_library[local_column]
    }

    bundle_result <- list(
      SourceBundle = source_bundle,
      counts = bundle_counts,
      ids = bundle_ids,
      libraries = bundle_libraries,
      audit = data.frame(SourceBundle = source_bundle, SourceGenes = nrow(raw), SourceCells = ncol(raw), SourceNNZ = length(raw@x), SelectedCells = nrow(selected_bundle), MinLibrarySize = min(full_library), MedianLibrarySize = median(full_library), MaxLibrarySize = max(full_library), Status = "PASS_ZERO_FILL_MISSING_BUNDLE_GENES")
    )
    saveRDS(bundle_result, checkpoint_file, compress = "gzip")
    cat("SCP259_CHECKPOINT_SAVED", source_bundle, nrow(selected_bundle), "cells\n")
    rm(raw, selected_raw, full_library, bundle_counts, bundle_ids, bundle_libraries)
    invisible(gc(FALSE))
  }

  for (cell_type in names(bundle_result$counts)) {
    counts_parts[[cell_type]][[source_bundle]] <- bundle_result$counts[[cell_type]]
    id_parts[[cell_type]][[source_bundle]] <- bundle_result$ids[[cell_type]]
    library_parts[[cell_type]][[source_bundle]] <- bundle_result$libraries[[cell_type]]
  }
  bundle_audit[[source_bundle]] <- bundle_result$audit
  rm(bundle_result)
}

# Step3. 按冻结顺序计算 CPM 并写父类型级输入
manifest_parts <- list()
for (cell_type in cell_types) {
  frozen_type <- selected_cells[CellType == cell_type][order(SelectionOrder)]
  counts_type <- do.call(cbind, counts_parts[[cell_type]])
  ids_type <- unlist(id_parts[[cell_type]], use.names = FALSE)
  library_type <- unlist(library_parts[[cell_type]], use.names = FALSE)
  frozen_order <- match(frozen_type$CellID, ids_type)
  if (anyNA(frozen_order) || anyDuplicated(ids_type) || ncol(counts_type) != nrow(frozen_type)) stop("SCP259 frozen cell assembly failed: ", cell_type, call. = FALSE)
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

  annotation <- merge(frozen_type[, .(ID = CellID, CellType, DonorID, BiologicalSampleID)], sample_metadata[, .(DonorID, BiologicalSampleID, Disease, Inflammation, Region3)], by = c("DonorID", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
  annotation <- annotation[match(frozen_type$CellID, ID)]
  if (anyNA(annotation$Disease) || anyNA(annotation$Inflammation) || anyNA(annotation$Region3)) stop("SCP259 frozen phenotype join failed: ", cell_type, call. = FALSE)
  fwrite(annotation, annotation_file_out, sep = "\t", quote = FALSE)

  manifest_parts[[cell_type]] <- data.frame(
    Dataset = "SCP259", CellType = cell_type, FrozenRank = rank_lookup[[cell_type]], Cells = nrow(frozen_type),
    ModelGenes = nrow(cpm_type), NonzeroModelGenes = sum(Matrix::rowSums(counts_type) > 0),
    MinLibrarySize = min(library_type), MedianLibrarySize = median(library_type), MaxLibrarySize = max(library_type),
    ExpressionFile = expression_file, ExpressionSHA256 = digest::digest(expression_file, algo = "sha256", file = TRUE, serialize = FALSE),
    AnnotationFile = annotation_file_out, AnnotationSHA256 = digest::digest(annotation_file_out, algo = "sha256", file = TRUE, serialize = FALSE)
  )
  cat("SCP259_PARENT_INPUT_WRITTEN", cell_type, nrow(frozen_type), "cells", nrow(cpm_type), "model genes\n")
  rm(counts_type, cpm_type, expression_table)
  invisible(gc(FALSE))
}

# Step4. 保存输入清单和 bundle 审计
manifest <- rbindlist(manifest_parts)
fwrite(manifest, file.path(prepare_dir, "SCP259_recovery_input_manifest.tsv"), sep = "\t")
fwrite(rbindlist(bundle_audit), file.path(prepare_dir, "SCP259_bundle_extraction_audit.tsv"), sep = "\t")
summary <- list(
  status = if (cohort == "IBD") "PASS_SCP259_RECOVERY_INPUT_PREPARATION" else "PASS_SCP259_NONIBD_RECOVERY_INPUT_PREPARATION", dataset = "SCP259", cohort = cohort, parent_types = nrow(manifest), cells = sum(manifest$Cells), source_bundles = bundle_config$SourceBundle,
  normalization = "raw counts collapsed by frozen model gene, full-feature raw UMI library CPM; EcoTyper applies log2(CPM+1) and row scaling",
  duplicate_source_rows = "summed before CPM", model_changed = FALSE, recovery_started = FALSE
)
jsonlite::write_json(summary, file.path(prepare_dir, "SCP259_recovery_input_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(prepare_dir, "R_sessionInfo.txt"))
cat("SCP259_RECOVERY_INPUT_PREPARATION_PASS", sum(manifest$Cells), "cells", nrow(manifest), "parent types\n")
