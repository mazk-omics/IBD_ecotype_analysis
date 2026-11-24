#!/usr/bin/env Rscript
# 本脚本只读计算四套外部 scRNA 原始表达源的逐文件 SHA256。
# 在服务器运行，输出文件清单和每队列聚合哈希，不读取表达矩阵内容到内存。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取运行目录和依赖
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 03b_freeze_scrna_source_hashes.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing R package: data.table", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("Missing R package: digest", call. = FALSE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Missing R package: jsonlite", call. = FALSE)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
input_audit_file <- file.path(audit_dir, "scrna_input_audit.json")
if (!file.exists(input_audit_file)) stop("Run 03_audit_scrna_inputs.R first", call. = FALSE)
input_audit <- jsonlite::read_json(input_audit_file, simplifyVector = TRUE)
if (!identical(input_audit$status, "PASS_SCRNA_METADATA_FREEZE_RAW_HASH_PENDING")) stop("scRNA metadata audit status changed", call. = FALSE)

# Step2. 冻结四套原始表达源
sources <- data.frame(
  Dataset = c("SCP1884", "GSE214695", "SCP259", "IBDverse_E_MTAB_16999"),
  Path = c(
    "/home/mazekai/enteric_glia/data_scRNA/SCP1884/expression",
    file.path(project_root, "06_validation/GSE214695_scRNA_recovery_20261009/source"),
    file.path(project_root, "01_raw_data/scRNA/SCP259/raw/SCP259/expression"),
    "/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad"
  ),
  stringsAsFactors = FALSE
)

file_parts <- list()
source_parts <- list()
for (i in seq_len(nrow(sources))) {
  dataset <- sources$Dataset[[i]]
  source_path <- sources$Path[[i]]
  if (!file.exists(source_path) && !dir.exists(source_path)) stop("Missing raw expression source: ", dataset, call. = FALSE)

  if (dir.exists(source_path)) {
    files <- sort(list.files(source_path, recursive = TRUE, full.names = TRUE, all.files = TRUE, include.dirs = FALSE, no.. = TRUE))
    relative_path <- substring(files, nchar(source_path) + 2L)
  } else {
    files <- source_path
    relative_path <- basename(source_path)
  }
  if (length(files) == 0L) stop("Raw expression source has no files: ", dataset, call. = FALSE)

  info <- file.info(files)
  sha256 <- character(length(files))
  cat("RAW_HASH_START", dataset, length(files), "files", sum(info$size), "bytes\n")
  for (j in seq_along(files)) {
    sha256[[j]] <- digest::digest(files[[j]], algo = "sha256", file = TRUE, serialize = FALSE)
    cat("RAW_HASH_FILE", dataset, j, "/", length(files), relative_path[[j]], info$size[[j]], "bytes\n")
  }

  file_manifest <- data.frame(
    Dataset = dataset,
    SourcePath = source_path,
    RelativePath = relative_path,
    Bytes = info$size,
    MTimeUTC = format(info$mtime, tz = "UTC", usetz = TRUE),
    SHA256 = sha256,
    stringsAsFactors = FALSE
  )
  aggregate_text <- paste(file_manifest$RelativePath, file_manifest$Bytes, file_manifest$SHA256, sep = "\t", collapse = "\n")
  source_manifest <- data.frame(
    Dataset = dataset,
    SourcePath = source_path,
    Files = nrow(file_manifest),
    Bytes = sum(file_manifest$Bytes),
    AggregateSHA256 = digest::digest(aggregate_text, algo = "sha256", serialize = FALSE),
    stringsAsFactors = FALSE
  )
  file_parts[[dataset]] <- file_manifest
  source_parts[[dataset]] <- source_manifest
  cat("RAW_HASH_DONE", dataset, source_manifest$AggregateSHA256, "\n")
}

# Step3. 保存清单并关闭输入 gate
raw_file_manifest <- data.table::rbindlist(file_parts)
raw_source_manifest <- data.table::rbindlist(source_parts)
data.table::fwrite(raw_file_manifest, file.path(audit_dir, "raw_expression_file_manifest.tsv"), sep = "\t")
data.table::fwrite(raw_source_manifest, file.path(audit_dir, "raw_expression_source_manifest.tsv"), sep = "\t")

summary <- list(
  status = "PASS_SCRNA_RAW_SOURCE_HASH_FREEZE",
  datasets = sources$Dataset,
  total_files = nrow(raw_file_manifest),
  total_bytes = sum(raw_file_manifest$Bytes),
  source_aggregate_sha256 = as.list(setNames(raw_source_manifest$AggregateSHA256, raw_source_manifest$Dataset)),
  source_paths = as.list(setNames(raw_source_manifest$SourcePath, raw_source_manifest$Dataset)),
  source_files_read_only = TRUE,
  model_changed = FALSE,
  automatic_recovery_started = FALSE,
  release_decision = "RAW_SOURCE_IDENTITY_COMPLETE_GENE_MAPPING_PENDING"
)
jsonlite::write_json(summary, file.path(audit_dir, "scrna_raw_source_hash_freeze.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo_source_hashes.txt"))
cat("SCRNA_RAW_SOURCE_HASH_FREEZE_PASS", nrow(raw_file_manifest), "files", sum(raw_file_manifest$Bytes), "bytes\n")
