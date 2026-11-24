#!/usr/bin/env Rscript
# 本脚本冻结四套 scRNA 原始基因 ID 到 symbol 的映射，并审计正式 W 的基因覆盖。
# 只读取小型 features/var 注释，不读取大表达矩阵，不启动 recovery。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取路径和前序 gate
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 03d_audit_scrna_gene_mapping.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite", "rhdf5")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
raw_hash_file <- file.path(audit_dir, "scrna_raw_source_hash_freeze.json")
cell_freeze_file <- file.path(audit_dir, "scrna_recovery_cell_freeze.json")
if (!file.exists(raw_hash_file) || !file.exists(cell_freeze_file)) stop("Complete source and cell freeze first", call. = FALSE)
raw_hash <- jsonlite::read_json(raw_hash_file, simplifyVector = TRUE)
cell_freeze <- jsonlite::read_json(cell_freeze_file, simplifyVector = TRUE)
if (!identical(raw_hash$status, "PASS_SCRNA_RAW_SOURCE_HASH_FREEZE")) stop("Raw source hash freeze did not pass", call. = FALSE)
if (!identical(cell_freeze$status, "PASS_SCRNA_RECOVERY_CELL_FREEZE")) stop("Recovery cell freeze did not pass", call. = FALSE)

# Step2. 读取四队列原始 gene ID-symbol 映射
mapping_parts <- list()

scp1884_root <- "/home/mazekai/enteric_glia/data_scRNA/SCP1884/expression"
scp1884_files <- sort(list.files(scp1884_root, pattern = "features\\.tsv$", recursive = TRUE, full.names = TRUE))
for (feature_file in scp1884_files) {
  feature <- fread(feature_file, header = FALSE, select = 1:2, col.names = c("GeneID", "Gene"))
  feature[, `:=`(Dataset = "SCP1884", FeatureFile = substring(feature_file, nchar(scp1884_root) + 2L))]
  mapping_parts[[paste0("SCP1884::", feature_file)]] <- feature
}

gse_root <- file.path(project_root, "06_validation/GSE214695_scRNA_recovery_20261009/source")
gse_files <- sort(list.files(gse_root, pattern = "features\\.tsv\\.gz$", full.names = TRUE))
feature <- fread(gse_files[[1]], header = FALSE, select = 1:2, col.names = c("GeneID", "Gene"))
for (feature_file in gse_files[-1]) {
  feature_check <- fread(feature_file, header = FALSE, select = 1:2, col.names = c("GeneID", "Gene"))
  if (!identical(feature, feature_check)) stop("GSE214695 decompressed feature mappings differ across samples", call. = FALSE)
}
feature[, `:=`(Dataset = "GSE214695", FeatureFile = basename(gse_files[[1]]))]
mapping_parts[["GSE214695"]] <- feature

scp259_root <- file.path(project_root, "01_raw_data/scRNA/SCP259/raw/SCP259/expression")
scp259_files <- sort(list.files(scp259_root, pattern = "\\.genes\\.tsv$", recursive = TRUE, full.names = TRUE))
for (feature_file in scp259_files) {
  feature <- fread(feature_file, header = FALSE, select = 1, col.names = "Gene")
  feature[, GeneID := Gene]
  feature[, `:=`(Dataset = "SCP259", FeatureFile = substring(feature_file, nchar(scp259_root) + 2L))]
  mapping_parts[[paste0("SCP259::", feature_file)]] <- feature
}

ibdverse_file <- "/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad"
ibdverse_gene_id <- as.character(rhdf5::h5read(ibdverse_file, "/var/_index"))
ibdverse_categories <- as.character(rhdf5::h5read(ibdverse_file, "/var/gene_symbols/categories"))
ibdverse_codes <- as.integer(rhdf5::h5read(ibdverse_file, "/var/gene_symbols/codes"))
ibdverse_gene <- rep(NA_character_, length(ibdverse_codes))
valid_code <- ibdverse_codes >= 0L & ibdverse_codes < length(ibdverse_categories)
ibdverse_gene[valid_code] <- ibdverse_categories[ibdverse_codes[valid_code] + 1L]
feature <- data.table(GeneID = ibdverse_gene_id, Gene = ibdverse_gene, Dataset = "IBDverse_E_MTAB_16999", FeatureFile = "Full_cohort.h5ad:/var")
mapping_parts[["IBDverse_E_MTAB_16999"]] <- feature

source_feature_mapping <- rbindlist(mapping_parts, use.names = TRUE)
source_feature_mapping[, GeneID := as.character(GeneID)]
source_feature_mapping[, Gene := as.character(Gene)]
source_feature_mapping <- source_feature_mapping[!is.na(Gene) & Gene != ""]
setcolorder(source_feature_mapping, c("Dataset", "FeatureFile", "GeneID", "Gene"))
source_mapping <- unique(source_feature_mapping[, .(Dataset, GeneID, Gene)])

# Step3. 审计重复 symbol 和正式模型基因覆盖
source_summary <- source_mapping[, .(
  FeatureRows = .N,
  UniqueGeneIDs = uniqueN(GeneID),
  UniqueSymbols = uniqueN(Gene),
  DuplicateSymbolRows = .N - uniqueN(Gene)
), by = Dataset]

rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_cell_type <- rank_data[[1L]]
rank_selected <- trimws(rank_data[[2L]])
if (anyDuplicated(rank_cell_type) || length(rank_cell_type) != 20L) stop("Frozen rank table changed", call. = FALSE)
cell_types <- sort(rank_cell_type)
coverage_parts <- list()
missing_parts <- list()
for (dataset in unique(source_mapping$Dataset)) {
  source_genes <- unique(source_mapping[Dataset == dataset, Gene])
  for (cell_type in cell_types) {
    selected_rank <- rank_selected[match(cell_type, rank_cell_type)]
    rank_dir <- file.path(state_discovery_root, cell_type, selected_rank)
    if (!dir.exists(rank_dir)) stop("Frozen rank directory missing: ", cell_type, call. = FALSE)
    initial_gene_info <- read.delim(file.path(rank_dir, "initial_gene_info.txt"), check.names = FALSE)
    W <- read.delim(file.path(rank_dir, "W.txt"), check.names = FALSE)
    w_genes <- sub("__pos$", "", rownames(W)[grepl("__pos$", rownames(W))])
    if (anyDuplicated(initial_gene_info$Gene) || anyDuplicated(w_genes)) stop("Frozen model gene identity changed: ", cell_type, call. = FALSE)

    initial_detected <- initial_gene_info$Gene %in% source_genes
    w_detected <- w_genes %in% source_genes
    coverage_parts[[paste(dataset, cell_type, sep = "::")]] <- data.frame(
      Dataset = dataset,
      CellType = cell_type,
      FrozenRank = selected_rank,
      InitialGenes = nrow(initial_gene_info),
      InitialGenesDetected = sum(initial_detected),
      InitialGeneCoverage = mean(initial_detected),
      WGenes = length(w_genes),
      WGenesDetected = sum(w_detected),
      WGeneCoverage = mean(w_detected),
      MissingWGenes = sum(!w_detected),
      stringsAsFactors = FALSE
    )
    if (any(!w_detected)) {
      missing_parts[[paste(dataset, cell_type, sep = "::")]] <- data.frame(Dataset = dataset, CellType = cell_type, Gene = w_genes[!w_detected], stringsAsFactors = FALSE)
    }
  }
}

coverage <- rbindlist(coverage_parts)
missing_w_genes <- rbindlist(missing_parts, use.names = TRUE, fill = TRUE)

# Step4. 保存冻结映射和 gate
data.table::fwrite(source_feature_mapping, file.path(audit_dir, "scrna_source_feature_mapping.tsv.gz"), sep = "\t")
data.table::fwrite(source_mapping, file.path(audit_dir, "scrna_source_gene_mapping.tsv.gz"), sep = "\t")
data.table::fwrite(source_summary, file.path(audit_dir, "scrna_source_gene_mapping_summary.tsv"), sep = "\t")
data.table::fwrite(coverage, file.path(audit_dir, "scrna_model_gene_coverage.tsv"), sep = "\t")
data.table::fwrite(missing_w_genes, file.path(audit_dir, "scrna_missing_W_genes.tsv"), sep = "\t")
if (any(coverage$WGeneCoverage < 0.90)) stop("At least one dataset-cell type has <90% frozen W gene coverage; diagnostic tables saved", call. = FALSE)

summary <- list(
  status = "PASS_SCRNA_GENE_MAPPING_FREEZE",
  datasets = sort(unique(source_mapping$Dataset)),
  model_cell_types = length(cell_types),
  minimum_W_gene_coverage = min(coverage$WGeneCoverage),
  median_W_gene_coverage = median(coverage$WGeneCoverage),
  mapping_sha256 = digest::digest(file.path(audit_dir, "scrna_source_gene_mapping.tsv.gz"), algo = "sha256", file = TRUE, serialize = FALSE),
  coverage_sha256 = digest::digest(file.path(audit_dir, "scrna_model_gene_coverage.tsv"), algo = "sha256", file = TRUE, serialize = FALSE),
  duplicate_symbols_will_be_summed_on_raw_counts = TRUE,
  missing_W_genes_will_follow_official_zero_fill = TRUE,
  model_changed = FALSE,
  automatic_recovery_started = FALSE,
  release_decision = "G4_SCRNA_INPUTS_PASS_RECOVERY_NOT_STARTED"
)
jsonlite::write_json(summary, file.path(audit_dir, "scrna_gene_mapping_audit.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo_gene_mapping.txt"))
cat("SCRNA_GENE_MAPPING_FREEZE_PASS", nrow(source_mapping), "mapping rows; minimum W coverage", min(coverage$WGeneCoverage), "\n")
