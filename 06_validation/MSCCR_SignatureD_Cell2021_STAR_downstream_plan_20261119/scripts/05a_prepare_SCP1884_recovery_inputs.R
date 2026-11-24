#!/usr/bin/env Rscript
# 本脚本从 SCP1884 六个原始 UMI bundle 构建冻结细胞的父类型级 EcoTyper CPM 输入。
# 逐 bundle 扫描稀疏矩阵并保存 checkpoint；完整 feature 库计算 library size。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞、表型和基因映射
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) stop("Usage: 05a_prepare_SCP1884_recovery_inputs.R RUN_ROOT [IBD|nonIBD]", call. = FALSE)
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
pool_file <- if (cohort == "IBD") file.path(eligible_pool_root, "SCP1884_eligible_ibd_cells.parquet") else file.path(cohort_audit_dir, "eligible_pools/SCP1884_eligible_nonIBD_cells.parquet")
selected_cells <- fread(file.path(cohort_audit_dir, selected_file), colClasses = "character")
selected_cells <- selected_cells[Dataset == "SCP1884"]
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]
eligible <- as.data.table(arrow::read_parquet(pool_file, as_data_frame = TRUE))
if (!all(c("CellID", "SourceBundle") %in% names(eligible))) stop("SCP1884 eligible-pool source columns changed", call. = FALSE)
selected_cells <- merge(selected_cells, unique(eligible[, .(CellID, SourceBundle)]), by = "CellID", all.x = TRUE, sort = FALSE)
if (anyNA(selected_cells$SourceBundle) || anyDuplicated(selected_cells$CellID)) stop("SCP1884 frozen cells do not map uniquely to bundles", call. = FALSE)

sample_metadata <- fread(file.path(cohort_audit_dir, phenotype_file), colClasses = "character")
sample_metadata <- sample_metadata[Dataset == "SCP1884"]
gene_mapping <- fread(file.path(audit_dir, "scrna_resolved_model_gene_mapping.tsv.gz"), colClasses = "character")
gene_mapping <- gene_mapping[Dataset == "SCP1884"]
rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
cell_types <- sort(unique(selected_cells$CellType))
if (length(cell_types) != 20L) stop("SCP1884 parent type count changed", call. = FALSE)
if (anyNA(rank_lookup[cell_types])) stop("Frozen rank is missing for SCP1884 parent types", call. = FALSE)
if (anyDuplicated(sample_metadata[, .(DonorID, BiologicalSampleID)])) stop("SCP1884 frozen phenotype contains duplicate samples", call. = FALSE)

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

source_root <- "/home/mazekai/enteric_glia/data_scRNA/SCP1884/expression"
bundle_config <- data.frame(
  SourceBundle = c("CO_EPI", "CO_IMM", "CO_STR", "TI_EPI", "TI_IMM", "TI_STR"),
  Subdir = c("Colon/Epithelial", "Colon/Immune", "Colon/Stromal", "TerminalIleum/Epithelial", "TerminalIleum/Immune", "TerminalIleum/Stromal")
)
prepare_root <- if (cohort == "IBD") "prepared_scRNA" else "prepared_scRNA_nonIBD"
prepare_dir <- file.path(RUN_ROOT, prepare_root, "SCP1884")
checkpoint_dir <- file.path(prepare_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(prepare_dir, "SCP1884_recovery_input_summary.json"))) stop("SCP1884 recovery input preparation already completed", call. = FALSE)

counts_parts <- setNames(vector("list", length(cell_types)), cell_types)
id_parts <- setNames(vector("list", length(cell_types)), cell_types)
library_parts <- setNames(vector("list", length(cell_types)), cell_types)
bundle_audit <- list()
nnz_block <- 5000000L

# Step2. 逐 bundle 提取冻结细胞和 frozen-W 基因
for (bundle_index in seq_len(nrow(bundle_config))) {
  source_bundle <- bundle_config$SourceBundle[[bundle_index]]
  source_subdir <- bundle_config$Subdir[[bundle_index]]
  checkpoint_file <- file.path(checkpoint_dir, paste0(source_bundle, "_selected_counts.rds"))

  if (file.exists(checkpoint_file)) {
    bundle_result <- readRDS(checkpoint_file)
    if (!identical(bundle_result$SourceBundle, source_bundle)) stop("SCP1884 checkpoint identity changed: ", source_bundle, call. = FALSE)
    cat("SCP1884_CHECKPOINT_REUSED", source_bundle, "\n")
  } else {
    selected_bundle <- selected_cells[SourceBundle == source_bundle]
    if (!nrow(selected_bundle)) stop("No frozen SCP1884 cells in bundle: ", source_bundle, call. = FALSE)
    barcode_file <- file.path(source_root, source_subdir, paste0(source_bundle, ".barcodes.tsv"))
    feature_file <- file.path(source_root, source_subdir, paste0(source_bundle, ".features.tsv"))
    matrix_file <- file.path(source_root, source_subdir, "raw.mtx")
    if (!all(file.exists(c(barcode_file, feature_file, matrix_file)))) stop("Missing SCP1884 raw bundle file: ", source_bundle, call. = FALSE)

    barcodes <- fread(barcode_file, header = FALSE, select = 1)[[1L]]
    features <- fread(feature_file, header = FALSE, select = 1:2, col.names = c("GeneID", "Gene"))
    features[, SourceKey := paste(GeneID, Gene, sep = "|")]
    if (anyDuplicated(barcodes) || anyDuplicated(features$SourceKey)) stop("Duplicated SCP1884 raw identity: ", source_bundle, call. = FALSE)
    raw_cell_index <- match(selected_bundle$CellID, barcodes)
    if (anyNA(raw_cell_index) || anyDuplicated(raw_cell_index)) stop("SCP1884 frozen barcodes do not map uniquely: ", source_bundle, call. = FALSE)

    bundle_cell_types <- sort(unique(selected_bundle$CellType))
    requested_source_keys <- unique(gene_mapping[CellType %in% bundle_cell_types]$SourceKey)
    union_source_keys <- requested_source_keys[requested_source_keys %in% features$SourceKey]
    union_source_row <- match(union_source_keys, features$SourceKey)
    if (!length(union_source_keys) || anyNA(union_source_row)) stop("No resolved SCP1884 source genes are available: ", source_bundle, call. = FALSE)
    source_gene_to_union <- integer(nrow(features))
    source_gene_to_union[union_source_row] <- seq_along(union_source_keys)
    source_cell_to_selected <- integer(length(barcodes))
    source_cell_to_selected[raw_cell_index] <- seq_len(nrow(selected_bundle))

    raw <- Matrix::readMM(matrix_file)
    if (!methods::is(raw, "TsparseMatrix")) raw <- methods::as(raw, "dgTMatrix")
    if (nrow(raw) != nrow(features) || ncol(raw) != length(barcodes)) stop("SCP1884 raw dimensions changed: ", source_bundle, call. = FALSE)
    if (any(!is.finite(raw@x)) || any(raw@x < 0) || any(abs(raw@x - round(raw@x)) > 1e-8)) stop("SCP1884 raw matrix is not nonnegative integer counts: ", source_bundle, call. = FALSE)

    union_counts <- sparseMatrix(i = integer(0), j = integer(0), dims = c(length(union_source_keys), nrow(selected_bundle)), dimnames = list(union_source_keys, selected_bundle$CellID))
    full_library <- numeric(nrow(selected_bundle))
    block_starts <- seq.int(1L, length(raw@x), by = nnz_block)

    for (block_index in seq_along(block_starts)) {
      from <- block_starts[[block_index]]
      to <- min(from + nnz_block - 1L, length(raw@x))
      source_gene <- raw@i[from:to] + 1L
      source_cell <- raw@j[from:to] + 1L
      block_value <- raw@x[from:to]
      selected_cell <- source_cell_to_selected[source_cell]
      keep_cell <- selected_cell > 0L
      if (any(keep_cell)) {
        full_library <- full_library + as.numeric(sparseMatrix(i = selected_cell[keep_cell], j = rep.int(1L, sum(keep_cell)), x = block_value[keep_cell], dims = c(nrow(selected_bundle), 1L)))
      }
      union_gene <- source_gene_to_union[source_gene]
      keep_model <- keep_cell & union_gene > 0L
      if (any(keep_model)) {
        partial <- sparseMatrix(i = union_gene[keep_model], j = selected_cell[keep_model], x = block_value[keep_model], dims = dim(union_counts))
        union_counts <- union_counts + partial
        rm(partial)
      }
      if (block_index %% 10L == 0L || block_index == length(block_starts)) cat("SCP1884_NNZ_BLOCK", source_bundle, block_index, "/", length(block_starts), "\n")
    }
    if (any(!is.finite(full_library)) || any(full_library <= 0)) stop("Invalid SCP1884 full library size: ", source_bundle, call. = FALSE)

    bundle_counts <- list()
    bundle_ids <- list()
    bundle_libraries <- list()
    for (cell_type in bundle_cell_types) {
      type_cells <- selected_bundle[CellType == cell_type]
      local_column <- match(type_cells$CellID, selected_bundle$CellID)
      mapping_type <- gene_mapping[CellType == cell_type & SourceKey %in% union_source_keys]
      if (!nrow(mapping_type)) stop("No SCP1884 mapped genes are available: ", cell_type, " ", source_bundle, call. = FALSE)
      source_column <- match(mapping_type$SourceKey, union_source_keys)
      if (anyNA(source_column)) stop("SCP1884 source collapse failed: ", cell_type, " ", source_bundle, call. = FALSE)
      collapse_matrix <- sparseMatrix(
        i = match(mapping_type$ModelGene, model_genes_by_type[[cell_type]]),
        j = source_column,
        x = 1,
        dims = c(length(model_genes_by_type[[cell_type]]), length(union_source_keys)),
        dimnames = list(model_genes_by_type[[cell_type]], union_source_keys)
      )
      model_counts <- collapse_matrix %*% union_counts[, local_column, drop = FALSE]
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
      audit = data.frame(SourceBundle = source_bundle, SourceGenes = nrow(raw), SourceCells = ncol(raw), SourceNNZ = length(raw@x), SelectedCells = nrow(selected_bundle), RequestedModelSourceGenes = length(requested_source_keys), AvailableModelSourceGenes = length(union_source_keys), MissingModelSourceGenes = length(setdiff(requested_source_keys, union_source_keys)), MinLibrarySize = min(full_library), MedianLibrarySize = median(full_library), MaxLibrarySize = max(full_library), Status = "PASS_ZERO_FILL_MISSING_BUNDLE_GENES")
    )
    saveRDS(bundle_result, checkpoint_file, compress = "gzip")
    cat("SCP1884_CHECKPOINT_SAVED", source_bundle, nrow(selected_bundle), "cells\n")
    rm(raw, union_counts, full_library, bundle_counts, bundle_ids, bundle_libraries)
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
  if (anyNA(frozen_order) || anyDuplicated(ids_type) || ncol(counts_type) != nrow(frozen_type)) stop("SCP1884 frozen cell assembly failed: ", cell_type, call. = FALSE)
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
    by = c("DonorID", "BiologicalSampleID"), all.x = TRUE, sort = FALSE
  )
  annotation <- annotation[match(frozen_type$CellID, ID)]
  if (anyNA(annotation$Disease) || anyNA(annotation$Inflammation) || anyNA(annotation$Region3)) stop("SCP1884 frozen phenotype join failed: ", cell_type, call. = FALSE)
  fwrite(annotation, annotation_file_out, sep = "\t", quote = FALSE)

  manifest_parts[[cell_type]] <- data.frame(
    Dataset = "SCP1884", CellType = cell_type, FrozenRank = rank_lookup[[cell_type]], Cells = nrow(frozen_type),
    ModelGenes = nrow(cpm_type), NonzeroModelGenes = sum(rowSums(counts_type) > 0),
    MinLibrarySize = min(library_type), MedianLibrarySize = median(library_type), MaxLibrarySize = max(library_type),
    ExpressionFile = expression_file, ExpressionSHA256 = digest::digest(expression_file, algo = "sha256", file = TRUE, serialize = FALSE),
    AnnotationFile = annotation_file_out, AnnotationSHA256 = digest::digest(annotation_file_out, algo = "sha256", file = TRUE, serialize = FALSE)
  )
  cat("SCP1884_PARENT_INPUT_WRITTEN", cell_type, nrow(frozen_type), "cells", nrow(cpm_type), "model genes\n")
  rm(counts_type, cpm_type, expression_table)
  invisible(gc(FALSE))
}

# Step4. 保存输入清单和 bundle 审计
manifest <- rbindlist(manifest_parts)
fwrite(manifest, file.path(prepare_dir, "SCP1884_recovery_input_manifest.tsv"), sep = "\t")
fwrite(rbindlist(bundle_audit), file.path(prepare_dir, "SCP1884_bundle_extraction_audit.tsv"), sep = "\t")
summary <- list(
  status = if (cohort == "IBD") "PASS_SCP1884_RECOVERY_INPUT_PREPARATION" else "PASS_SCP1884_NONIBD_RECOVERY_INPUT_PREPARATION",
  dataset = "SCP1884",
  cohort = cohort,
  parent_types = nrow(manifest),
  cells = sum(manifest$Cells),
  source_bundles = bundle_config$SourceBundle,
  normalization = "raw counts collapsed by frozen model gene, full-feature raw UMI library CPM; EcoTyper applies log2(CPM+1) and row scaling",
  duplicate_source_rows = "summed before CPM",
  model_changed = FALSE,
  recovery_started = FALSE
)
jsonlite::write_json(summary, file.path(prepare_dir, "SCP1884_recovery_input_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(prepare_dir, "R_sessionInfo.txt"))
cat("SCP1884_RECOVERY_INPUT_PREPARATION_PASS", sum(manifest$Cells), "cells", nrow(manifest), "parent types\n")
