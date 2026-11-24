#!/usr/bin/env Rscript
# 本脚本统一验收四套 nonIBD scRNA recovery 表达输入，不修改矩阵。
# 核对冻结细胞顺序、注释、数值范围和文件哈希，验收通过后才允许进入 EcoTyper。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结细胞和四套输入清单
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 29_audit_nonIBD_recovery_inputs.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

audit_dir <- file.path(output_root, "18_nonIBD_scrna_input_audit")
prepare_root <- file.path(RUN_ROOT, "prepared_scRNA_nonIBD")
selected_cells <- fread(file.path(audit_dir, "nonIBD_recovery_selected_cells.tsv.gz"), colClasses = "character")
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]
datasets <- c("SCP1884", "GSE214695", "SCP259", "IBDverse_E_MTAB_16999")
manifest_files <- c(
  SCP1884 = "SCP1884_recovery_input_manifest.tsv",
  GSE214695 = "GSE214695_recovery_input_manifest.tsv",
  SCP259 = "SCP259_recovery_input_manifest.tsv",
  IBDverse_E_MTAB_16999 = "IBDverse_recovery_input_manifest.tsv"
)

# Step2. 逐 dataset × parent 核对表达和注释
audit_parts <- list()
for (dataset in datasets) {
  manifest <- fread(file.path(prepare_root, dataset, manifest_files[[dataset]]), colClasses = "character")
  manifest[, `:=`(Cells = as.integer(Cells), ModelGenes = as.integer(ModelGenes))]
  selected_dataset <- selected_cells[Dataset == dataset]
  if (sum(manifest$Cells) != nrow(selected_dataset) || !setequal(manifest$CellType, unique(selected_dataset$CellType))) stop("nonIBD manifest dimensions changed: ", dataset, call. = FALSE)

  for (row_index in seq_len(nrow(manifest))) {
    cell_type <- manifest$CellType[[row_index]]
    expression_file <- manifest$ExpressionFile[[row_index]]
    annotation_file <- manifest$AnnotationFile[[row_index]]
    frozen_type <- selected_dataset[CellType == cell_type][order(SelectionOrder)]
    expression <- fread(expression_file, check.names = FALSE)
    annotation <- fread(annotation_file, colClasses = "character", check.names = FALSE)

    genes <- expression[[1L]]
    values <- as.matrix(expression[, -1L])
    storage.mode(values) <- "double"
    cell_ids <- colnames(expression)[-1L]
    if (nrow(values) != manifest$ModelGenes[[row_index]] || ncol(values) != manifest$Cells[[row_index]]) stop("nonIBD expression dimensions changed: ", dataset, " ", cell_type, call. = FALSE)
    if (anyNA(genes) || anyDuplicated(genes)) stop("nonIBD expression genes are not unique: ", dataset, " ", cell_type, call. = FALSE)
    if (!identical(cell_ids, frozen_type$CellID) || !identical(annotation$ID, frozen_type$CellID)) stop("nonIBD frozen cell order changed: ", dataset, " ", cell_type, call. = FALSE)
    if (any(annotation$CellType != cell_type) || any(annotation$Disease != "nonIBD") || any(annotation$Inflammation != "NonIBD") || anyNA(annotation$Region3)) stop("nonIBD annotation changed: ", dataset, " ", cell_type, call. = FALSE)
    if (any(!is.finite(values)) || any(values < 0)) stop("nonIBD expression contains invalid values: ", dataset, " ", cell_type, call. = FALSE)
    if (digest::digest(expression_file, algo = "sha256", file = TRUE, serialize = FALSE) != manifest$ExpressionSHA256[[row_index]]) stop("nonIBD expression hash changed: ", dataset, " ", cell_type, call. = FALSE)
    if (digest::digest(annotation_file, algo = "sha256", file = TRUE, serialize = FALSE) != manifest$AnnotationSHA256[[row_index]]) stop("nonIBD annotation hash changed: ", dataset, " ", cell_type, call. = FALSE)

    model_cpm_sum <- colSums(values)
    log_values <- log2(values + 1)
    audit_parts[[paste(dataset, cell_type, sep = "::")]] <- data.frame(
      Dataset = dataset, CellType = cell_type, Cells = ncol(values), ModelGenes = nrow(values), DetectedGenes = sum(rowSums(values) > 0),
      ConstantGenesAfterLog2 = sum(apply(log_values, 1L, function(x) max(x) == min(x))),
      ZeroModelSignalCells = sum(model_cpm_sum <= 0),
      MinModelCPMSum = min(model_cpm_sum), MedianModelCPMSum = median(model_cpm_sum), MaxModelCPMSum = max(model_cpm_sum),
      ExactFrozenCellOrder = TRUE, AnnotationFrozen = TRUE, HashesMatch = TRUE, FiniteNonnegative = TRUE
    )
    cat("NONIBD_INPUT_AUDITED", dataset, cell_type, ncol(values), "cells", nrow(values), "genes\n")
    rm(expression, annotation, values, log_values)
    invisible(gc(FALSE))
  }
}

# Step3. 保存统一验收结论
input_audit <- rbindlist(audit_parts)
if (nrow(input_audit) != 72L || sum(input_audit$Cells) != 112001L) stop("nonIBD audited input totals changed", call. = FALSE)
if (any(input_audit$MaxModelCPMSum > 1e6 + 1e-6)) stop("nonIBD model-gene CPM sums exceed the full-library bound", call. = FALSE)
fwrite(input_audit, file.path(audit_dir, "nonIBD_recovery_input_audit.tsv"), sep = "\t")

dataset_summary <- input_audit[, .(
  ParentTypes = .N, Cells = sum(Cells), ZeroModelSignalCells = sum(ZeroModelSignalCells), MinDetectedGenes = min(DetectedGenes), MaxConstantGenesAfterLog2 = max(ConstantGenesAfterLog2),
  MinModelCPMSum = min(MinModelCPMSum), MedianModelCPMSum = median(MedianModelCPMSum), MaxModelCPMSum = max(MaxModelCPMSum)
), by = Dataset]
fwrite(dataset_summary, file.path(audit_dir, "nonIBD_recovery_input_audit_by_dataset.tsv"), sep = "\t")

summary <- list(
  status = "PASS_NONIBD_RECOVERY_INPUT_AUDIT", datasets = datasets, dataset_parent_groups = nrow(input_audit), cells = sum(input_audit$Cells),
  exact_frozen_cell_order = all(input_audit$ExactFrozenCellOrder), annotation_frozen = all(input_audit$AnnotationFrozen),
  hashes_match = all(input_audit$HashesMatch), finite_nonnegative = all(input_audit$FiniteNonnegative),
  zero_model_signal_cells_recorded_not_filtered = sum(input_audit$ZeroModelSignalCells),
  cpm_definition = "full-library CPM restricted to frozen W genes; model-gene column sums may be below 1e6",
  current_model_changed = FALSE, existing_ibd_results_changed = FALSE, official_recovery_started = FALSE
)
jsonlite::write_json(summary, file.path(audit_dir, "nonIBD_recovery_input_audit.json"), pretty = TRUE, auto_unbox = TRUE)
cat("NONIBD_RECOVERY_INPUT_AUDIT_PASS", sum(input_audit$Cells), "cells in", nrow(input_audit), "dataset-parent groups\n")
