#!/usr/bin/env Rscript
# 本脚本审计 GSE214695 父类型级 EcoTyper 输入，不修改表达矩阵。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结清单和输入清单
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04b_audit_GSE214695_recovery_inputs.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
prepare_dir <- file.path(RUN_ROOT, "prepared_scRNA/GSE214695")
selected_cells <- fread(file.path(audit_dir, "recovery_selected_cells.tsv.gz"), colClasses = "character")
selected_cells <- selected_cells[Dataset == "GSE214695"]
selected_cells[, SelectionOrder := as.integer(SelectionOrder)]
manifest <- fread(file.path(prepare_dir, "GSE214695_recovery_input_manifest.tsv"), colClasses = "character")
manifest[, Cells := as.integer(Cells)]
manifest[, ModelGenes := as.integer(ModelGenes)]
if (nrow(manifest) != 18L || sum(manifest$Cells) != nrow(selected_cells)) stop("GSE214695 input manifest dimensions changed", call. = FALSE)

# Step2. 逐父类型核验维度、顺序、数值和哈希
audit_parts <- vector("list", nrow(manifest))
for (row_index in seq_len(nrow(manifest))) {
  cell_type <- manifest$CellType[[row_index]]
  expression_file <- manifest$ExpressionFile[[row_index]]
  annotation_file <- manifest$AnnotationFile[[row_index]]
  frozen_type <- selected_cells[CellType == cell_type][order(SelectionOrder)]
  expression <- fread(expression_file, check.names = FALSE)
  annotation <- fread(annotation_file, colClasses = "character", check.names = FALSE)

  genes <- expression[[1L]]
  values <- as.matrix(expression[, -1L])
  storage.mode(values) <- "double"
  cell_ids <- colnames(expression)[-1L]
  if (nrow(values) != manifest$ModelGenes[[row_index]] || ncol(values) != manifest$Cells[[row_index]]) stop("Expression dimensions changed: ", cell_type, call. = FALSE)
  if (anyNA(genes) || anyDuplicated(genes)) stop("Expression genes are not unique: ", cell_type, call. = FALSE)
  if (!identical(cell_ids, frozen_type$CellID)) stop("Expression cell order differs from frozen order: ", cell_type, call. = FALSE)
  if (!identical(annotation$ID, frozen_type$CellID)) stop("Annotation cell order differs from frozen order: ", cell_type, call. = FALSE)
  if (any(annotation$CellType != cell_type)) stop("Annotation parent type changed: ", cell_type, call. = FALSE)
  if (any(!is.finite(values)) || any(values < 0)) stop("Expression contains invalid values: ", cell_type, call. = FALSE)
  if (digest::digest(expression_file, algo = "sha256", file = TRUE, serialize = FALSE) != manifest$ExpressionSHA256[[row_index]]) stop("Expression hash changed: ", cell_type, call. = FALSE)
  if (digest::digest(annotation_file, algo = "sha256", file = TRUE, serialize = FALSE) != manifest$AnnotationSHA256[[row_index]]) stop("Annotation hash changed: ", cell_type, call. = FALSE)

  model_cpm_sum <- colSums(values)
  log_values <- log2(values + 1)
  constant_genes <- sum(apply(log_values, 1L, function(x) max(x) == min(x)))
  audit_parts[[row_index]] <- data.frame(
    Dataset = "GSE214695",
    CellType = cell_type,
    Cells = ncol(values),
    ModelGenes = nrow(values),
    DetectedGenes = sum(rowSums(values) > 0),
    ConstantGenesAfterLog2 = constant_genes,
    MinModelCPMSum = min(model_cpm_sum),
    MedianModelCPMSum = median(model_cpm_sum),
    MaxModelCPMSum = max(model_cpm_sum),
    ExactFrozenCellOrder = TRUE,
    HashesMatch = TRUE,
    FiniteNonnegative = TRUE,
    stringsAsFactors = FALSE
  )
  cat("GSE_INPUT_AUDITED", cell_type, ncol(values), "cells", nrow(values), "genes\n")
  rm(expression, annotation, values, log_values)
  invisible(gc(FALSE))
}

# Step3. 保存审计结论
input_audit <- rbindlist(audit_parts)
if (any(input_audit$MaxModelCPMSum > 1e6 + 1e-6) || any(input_audit$MinModelCPMSum <= 0)) stop("Model-gene CPM sums are outside expected bounds", call. = FALSE)
fwrite(input_audit, file.path(prepare_dir, "GSE214695_recovery_input_audit.tsv"), sep = "\t")
summary <- list(
  status = "PASS_GSE214695_RECOVERY_INPUT_AUDIT",
  dataset = "GSE214695",
  parent_types = nrow(input_audit),
  cells = sum(input_audit$Cells),
  exact_frozen_cell_order = all(input_audit$ExactFrozenCellOrder),
  hashes_match = all(input_audit$HashesMatch),
  finite_nonnegative = all(input_audit$FiniteNonnegative),
  cpm_definition = "full-library CPM restricted to frozen W genes; model-gene column sums may be below 1e6",
  recovery_started = FALSE
)
jsonlite::write_json(summary, file.path(prepare_dir, "GSE214695_recovery_input_audit.json"), pretty = TRUE, auto_unbox = TRUE)
cat("GSE214695_RECOVERY_INPUT_AUDIT_PASS", sum(input_audit$Cells), "cells", nrow(input_audit), "parent types\n")
