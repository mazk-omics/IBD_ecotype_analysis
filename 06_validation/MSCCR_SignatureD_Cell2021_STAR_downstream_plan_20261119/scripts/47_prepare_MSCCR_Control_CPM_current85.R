#!/usr/bin/env Rscript
# 本脚本从 MSCCR 原始 counts 重建当前 85-state 模型所需基因的 Control CPM 输入，不重建 TPM，也不重复 discovery。
# CPM 分母使用每个样本的全部原始基因行；仅在计算完成后筛出当前模型基因，并按 Ileum、Colon、Rectum 分开保存。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定原始输入、当前模型和 Control 样本
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 47_prepare_MSCCR_Control_CPM_current85.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

project_root <- normalizePath(file.path(RUN_ROOT, "../.."), mustWork = TRUE)
raw_counts_file <- file.path(project_root, "01_raw_data/bulk/MSCCR_GSE193677/GSE193677_MSCCR_Biopsy_counts.txt.gz")
identity_file <- file.path(project_root, "03_reference/web_reference_three_pilots/msccr_gene_identity_20261006/MSCCR_Ensembl75_HGNC_identity.tsv")
annotation_file <- file.path(project_root, "05_EcoTyper/input/MSCCR_CORE12/annotation/MSCCR_CORE12_annotation.tsv")
discovery_name <- "MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101"
model_root <- file.path(project_root, "software/ecotyper/EcoTyper", discovery_name, "Custom")
output_dir <- file.path(RUN_ROOT, "output/35_MSCCR_Control_CPM_current85_preflight/reconstructed_current_model_gene_union")
input_dir <- file.path(output_dir, "input")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(raw_counts_file, identity_file, annotation_file, file.path(model_root, "Analysis/rank_selection/rank_data.txt"))
if (any(!file.exists(required_files))) stop("MSCCR raw counts, identity, annotation or current model is missing", call. = FALSE)

annotation <- fread(annotation_file, colClasses = "character")
required_annotation <- c("ID", "PatientID", "Disease", "Inflammation", "Region3", "RegionDetailed", "RegionBroad")
if (!all(required_annotation %in% names(annotation))) stop("MSCCR annotation columns changed", call. = FALSE)
controls <- annotation[Disease == "Control" & Inflammation == "Non-inflamed" & Region3 %in% c("Ileum", "Colon", "Rectum")]
if (nrow(controls) != 461L || uniqueN(controls$PatientID) != 243L) stop("Frozen MSCCR Control denominator changed", call. = FALSE)
print(controls[, .(Samples = .N, Donors = uniqueN(PatientID)), by = Region3])

# Step2. 冻结当前 20 个父类型的正式模型基因
rank_data <- fread(file.path(model_root, "Analysis/rank_selection/rank_data.txt"))
if (nrow(rank_data) != 20L || !all(c("CellType", "Chosen_Rank") %in% names(rank_data))) stop("Current model rank data changed", call. = FALSE)
model_gene_parts <- vector("list", nrow(rank_data))

for (cell_index in seq_len(nrow(rank_data))) {
  gene_info_file <- file.path(
    model_root, "Cell_States/discovery", rank_data$CellType[[cell_index]],
    as.character(rank_data$Chosen_Rank[[cell_index]]), "initial_gene_info.txt"
  )
  gene_info <- fread(gene_info_file, colClasses = "character")
  model_gene_parts[[cell_index]] <- data.table(CellType = rank_data$CellType[[cell_index]], Gene = unique(gene_info$Gene))
}

model_genes <- rbindlist(model_gene_parts)
model_gene_union <- sort(unique(model_genes$Gene))
cat("CURRENT_MODEL_GENES", length(model_gene_union), "\n")

# Step3. 读取 461 个 Control 的全部原始 counts 并计算全库 CPM 分母
raw_header <- fread(cmd = paste("zcat", shQuote(raw_counts_file)), nrows = 0L, showProgress = FALSE)
raw_sample_names <- names(raw_header)[-1L]
if (!all(controls$ID %in% raw_sample_names)) stop("Control samples are missing from raw counts", call. = FALSE)
selected_columns <- c(1L, match(controls$ID, names(raw_header)))
raw_counts <- fread(cmd = paste("zcat", shQuote(raw_counts_file)), select = selected_columns, showProgress = TRUE)
setnames(raw_counts, 1L, "EnsemblID")
if (!identical(names(raw_counts)[-1L], controls$ID)) stop("Raw count sample order changed", call. = FALSE)
if (anyDuplicated(raw_counts$EnsemblID)) stop("Raw Ensembl IDs are duplicated", call. = FALSE)

count_matrix <- as.matrix(raw_counts[, -1L])
storage.mode(count_matrix) <- "double"
if (any(!is.finite(count_matrix)) || any(count_matrix < 0)) stop("Raw counts contain invalid values", call. = FALSE)
library_size <- colSums(count_matrix)
if (any(!is.finite(library_size)) || any(library_size <= 0)) stop("Invalid full-library denominators", call. = FALSE)
cat("RAW_COUNTS_READY", nrow(count_matrix), "genes", ncol(count_matrix), "samples\n")

# Step4. 映射当前模型基因并在 symbol 层合并重复 Ensembl 行
identity <- fread(identity_file, colClasses = "character")
required_identity <- c("msccr_ensembl_id", "ensembl75_gene_name", "ensembl75_gene_biotype", "mapping_status")
if (!all(required_identity %in% names(identity))) stop("MSCCR gene identity columns changed", call. = FALSE)
identity <- identity[
  mapping_status == "unique_gene_symbol" & ensembl75_gene_biotype == "protein_coding" & ensembl75_gene_name %in% model_gene_union,
  .(EnsemblID = msccr_ensembl_id, Gene = ensembl75_gene_name)
]
identity <- identity[EnsemblID %in% raw_counts$EnsemblID]
if (anyDuplicated(identity$EnsemblID)) stop("Gene identity has duplicate Ensembl IDs", call. = FALSE)

mapped_index <- match(identity$EnsemblID, raw_counts$EnsemblID)
mapped_counts <- count_matrix[mapped_index, , drop = FALSE]
counts_by_gene <- rowsum(mapped_counts, group = identity$Gene, reorder = FALSE)
cpm_by_gene <- sweep(counts_by_gene, 2L, library_size, "/") * 1e6
cpm_by_gene <- cpm_by_gene[order(rownames(cpm_by_gene)), , drop = FALSE]
if (any(!is.finite(cpm_by_gene)) || any(cpm_by_gene < 0) || max(cpm_by_gene) <= 50) stop("Reconstructed CPM failed numerical checks", call. = FALSE)

coverage <- model_genes[, .(
  ModelGenes = uniqueN(Gene), PresentGenes = uniqueN(Gene[Gene %in% rownames(cpm_by_gene)])
), by = CellType]
coverage[, `:=`(MissingGenes = ModelGenes - PresentGenes, PresentFraction = PresentGenes / ModelGenes)]
setorder(coverage, CellType)
print(coverage)
if (min(coverage$PresentFraction) < 0.75) stop("Current model gene coverage is below the prespecified 75% stop line", call. = FALSE)

# Step5. 按部位保存同尺度 CPM 和多列 canonical annotation
regions <- c("Ileum", "Colon", "Rectum")
region_manifest <- vector("list", length(regions))

for (region_index in seq_along(regions)) {
  region <- regions[[region_index]]
  region_annotation <- controls[Region3 == region]
  region_cpm <- cpm_by_gene[, region_annotation$ID, drop = FALSE]
  expression_file <- file.path(input_dir, paste0("MSCCR_Control_CPM_current85_", region, ".txt"))
  canonical_annotation_file <- file.path(input_dir, paste0("MSCCR_Control_CPM_current85_", region, "_annotation.tsv"))

  expression_output <- data.table(GeneSymbol = rownames(region_cpm))
  expression_output <- cbind(expression_output, as.data.table(region_cpm))
  canonical_annotation <- region_annotation[, .(
    ID, PatientID, Disease = "Non-IBD", DiseaseSubtype = "Control", Tissue = Region3,
    Histopathology = Inflammation, Region3, RegionDetailed, RegionBroad
  )]
  if (!identical(names(expression_output)[-1L], canonical_annotation$ID)) stop("Region sample order mismatch", call. = FALSE)
  fwrite(expression_output, expression_file, sep = "\t", quote = FALSE)
  fwrite(canonical_annotation, canonical_annotation_file, sep = "\t", quote = FALSE)

  region_manifest[[region_index]] <- data.table(
    Region3 = region, Samples = nrow(canonical_annotation), Donors = uniqueN(canonical_annotation$PatientID),
    Genes = nrow(expression_output), Scale = "linear CPM; all raw Ensembl rows denominator; no TPM rebuild",
    Role = "internal nonIBD coverage only; not an independent discovery validation",
    ExpressionFile = expression_file, AnnotationFile = canonical_annotation_file,
    ExpressionMD5 = unname(tools::md5sum(expression_file)), AnnotationMD5 = unname(tools::md5sum(canonical_annotation_file))
  )
  cat("MSCCR_CONTROL_INPUT_READY", region, nrow(canonical_annotation), "samples\n")
}

# Step6. 保存输入身份和完成状态
region_manifest <- rbindlist(region_manifest)
library_audit <- data.table(ID = controls$ID, PatientID = controls$PatientID, Region3 = controls$Region3, FullRawLibrary = library_size)
missing_genes <- data.table(Gene = setdiff(model_gene_union, rownames(cpm_by_gene)))
input_identity <- data.table(File = required_files, MD5 = unname(tools::md5sum(required_files)))
fwrite(region_manifest, file.path(output_dir, "MSCCR_Control_CPM_current85_region_manifest.tsv"), sep = "\t", quote = FALSE)
fwrite(coverage, file.path(output_dir, "MSCCR_Control_CPM_current85_gene_coverage.tsv"), sep = "\t", quote = FALSE)
fwrite(missing_genes, file.path(output_dir, "MSCCR_Control_CPM_current85_missing_genes.tsv"), sep = "\t", quote = FALSE)
fwrite(library_audit, file.path(output_dir, "MSCCR_Control_CPM_current85_library_audit.tsv"), sep = "\t", quote = FALSE)
fwrite(input_identity, file.path(output_dir, "MSCCR_Control_CPM_current85_input_identity.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_MSCCR_CONTROL_CPM_CURRENT85_PREPARATION", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  samples = nrow(controls), donors = uniqueN(controls$PatientID), regions = regions, current_model_cell_types = nrow(rank_data),
  model_gene_union = length(model_gene_union), mapped_model_genes = nrow(cpm_by_gene), minimum_parent_coverage = min(coverage$PresentFraction),
  scale = "linear CPM using all raw Ensembl rows as denominator; no TPM rebuild",
  role = "internal nonIBD coverage only; discovery CD/UC samples are not rerun"
)
jsonlite::write_json(completion, file.path(output_dir, "MSCCR_Control_CPM_current85_preparation_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(completion)
cat("PASS_MSCCR_CONTROL_CPM_CURRENT85_PREPARATION\n")
