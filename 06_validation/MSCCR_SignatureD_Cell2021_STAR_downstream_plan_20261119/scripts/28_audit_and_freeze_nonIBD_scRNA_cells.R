#!/usr/bin/env Rscript
# 本脚本从四套外部 scRNA 原始元数据中冻结 nonIBD/健康细胞，用于当前 85-state 模型的独立 recovery。
# 只处理元数据和细胞 ID；不改写既有 IBD 输入、state 模型或 recovery 结果。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取路径和固定参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 28_audit_and_freeze_nonIBD_scRNA_cells.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("arrow", "data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

audit_dir <- file.path(output_root, "18_nonIBD_scrna_input_audit")
pool_dir <- file.path(audit_dir, "eligible_pools")
dir.create(pool_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(audit_dir, "nonIBD_input_freeze.json"))) stop("nonIBD input freeze already completed", call. = FALSE)

joint20 <- c(
  "B", "Plasma", "CD4_T", "CD8_T", "NK", "ILC", "Monocyte_Macrophage", "DC", "Mast",
  "Fibroblast", "Pericyte", "Endothelial", "Ileal_absorptive", "Colonic_absorptive", "Goblet",
  "Stem_TA_progenitor", "Paneth", "Enteroendocrine", "Glial", "Tuft"
)
supported <- list(
  SCP1884 = joint20,
  GSE214695 = setdiff(joint20, c("Ileal_absorptive", "Paneth")),
  SCP259 = setdiff(joint20, c("Ileal_absorptive", "Paneth")),
  IBDverse_E_MTAB_16999 = setdiff(joint20, c("NK", "Colonic_absorptive", "Paneth", "Glial"))
)
max_cells <- 2500L
sampling_seed <- 1234L

scp1884_file <- file.path(project_root, "03_reference/SCP1884/audit_v1/output/02_candidate_taxonomy/scp1884_candidate_cells_provisional_taxonomy.parquet")
gse_file <- file.path(project_root, "06_validation/deconv_fraction_sensitivity_joint20_20261013/08_external_pseudobulk/GSE214695/00_audit/eligible_author_cells.tsv.gz")
scp259_file <- file.path(project_root, "03_reference/SCP259/audit_v1/output/01_candidate_taxonomy_coverage/scp259_candidate_cells_provisional_taxonomy.parquet")
ibdverse_file <- file.path(project_root, "06_validation/deconv_fraction_sensitivity_joint20_20261013/08_external_pseudobulk/IBDverse/00_audit/eligible_author_cells.tsv.gz")
source_files <- c(SCP1884 = scp1884_file, GSE214695 = gse_file, SCP259 = scp259_file, IBDverse_E_MTAB_16999 = ibdverse_file)
if (any(!file.exists(source_files))) stop("Missing nonIBD metadata source: ", paste(names(source_files)[!file.exists(source_files)], collapse = ", "), call. = FALSE)

# Step2. 按上一轮 Joint20 规则建立四套 nonIBD 细胞池
scp1884 <- as.data.table(arrow::read_parquet(scp1884_file, as_data_frame = TRUE))
scp1884 <- scp1884[Type == "Heal"]
scp1884[, CellType := ""]
scp1884[anno2 == "Plasma cells", CellType := "Plasma"]
scp1884[grepl("^B cells", anno2), CellType := "B"]
scp1884[grepl("^T cells CD4", anno2) | anno2 %in% c("Tregs", "T cells Naive CD4+"), CellType := "CD4_T"]
scp1884[grepl("^T cells CD8", anno2), CellType := "CD8_T"]
scp1884[anno2 == "NK cells KLRF1+ CD3G-", CellType := "NK"]
scp1884[anno2 == "ILCs", CellType := "ILC"]
scp1884[grepl("^(Macrophages|Monocytes)", anno2), CellType := "Monocyte_Macrophage"]
scp1884[grepl("^(DC1|DC2)", anno2) | anno2 == "Mature DCs", CellType := "DC"]
scp1884[anno2 == "Mast cells", CellType := "Mast"]
scp1884[grepl("^(Fibroblasts|Myofibroblasts|Activated fibroblasts|Inflammatory fibroblasts)", anno2), CellType := "Fibroblast"]
scp1884[grepl("^Pericytes", anno2), CellType := "Pericyte"]
scp1884[grepl("^Endothelial cells", anno2) | anno2 == "Lymphatics", CellType := "Endothelial"]
scp1884[anno2 == "Glial cells", CellType := "Glial"]
scp1884[grepl("^Enterocytes", anno2) & Site == "CO", CellType := "Colonic_absorptive"]
scp1884[grepl("^Enterocytes", anno2) & Site %in% c("TI", "SB"), CellType := "Ileal_absorptive"]
scp1884[grepl("^Stem cells", anno2), CellType := "Stem_TA_progenitor"]
scp1884[grepl("^Goblet cells", anno2), CellType := "Goblet"]
scp1884[anno2 == "Paneth cells", CellType := "Paneth"]
scp1884[anno2 %in% c("L cells", "Enterochromaffin cells", "Enteroendocrine cells"), CellType := "Enteroendocrine"]
scp1884[anno2 == "Tuft cells", CellType := "Tuft"]
scp1884 <- scp1884[CellType %in% supported$SCP1884]
scp1884_pool <- scp1884[, .(
  Dataset = "SCP1884", CellID = cell_id, CellType, AuthorLabel = anno2, DonorID = as.character(PubID),
  BiologicalSampleID = as.character(biological_sample_id), SourceUnitID = as.character(PubIDSample),
  Region = as.character(Site), Inflammation = "NonIBD", SourceBundle = as.character(expression_bundle), SourceRow = -1L
)]
rm(scp1884)
invisible(gc(FALSE))

gse <- fread(gse_file)
gse <- gse[Disease == "HC" & CellType %in% supported$GSE214695]
gse_pool <- gse[, .(
  Dataset = "GSE214695", CellID, CellType, AuthorLabel, DonorID = as.character(Donor),
  BiologicalSampleID = as.character(Sample), SourceUnitID = as.character(Accession), Region = as.character(Region3),
  Inflammation = "NonIBD", SourceBundle = as.character(Sample), SourceRow = as.integer(SourceRow)
)]
rm(gse)
invisible(gc(FALSE))

scp259 <- as.data.table(arrow::read_parquet(scp259_file, as_data_frame = TRUE))
scp259 <- scp259[disease == "Healthy" & coarse_identity %in% supported$SCP259]
scp259_pool <- scp259[, .(
  Dataset = "SCP259", CellID = cell_id, CellType = coarse_identity, AuthorLabel = author_annotation,
  DonorID = as.character(donor_id), BiologicalSampleID = as.character(sample_id), SourceUnitID = as.character(sample_id),
  Region = as.character(site), Inflammation = "NonIBD", SourceBundle = as.character(expression_bundle), SourceRow = -1L
)]
rm(scp259)
invisible(gc(FALSE))

ibdverse <- fread(ibdverse_file)
ibdverse <- ibdverse[disease_status == "Healthy" & metadata_eligible & counts_QC_eligible & eligible & CellType %in% supported$IBDverse_E_MTAB_16999]
ibdverse_pool <- ibdverse[, .(
  Dataset = "IBDverse_E_MTAB_16999", CellID, CellType, AuthorLabel, DonorID = as.character(Donor),
  BiologicalSampleID = as.character(Sample), SourceUnitID = as.character(sanger_sample_id), Region = as.character(biopsy_type),
  Inflammation = "NonIBD", SourceBundle = "Full_cohort.h5ad", SourceRow = as.integer(SourceRow)
)]
rm(ibdverse)
invisible(gc(FALSE))

# Step3. 统一表型、核对身份并写 eligible pools
pool_objects <- list(SCP1884 = scp1884_pool, GSE214695 = gse_pool, SCP259 = scp259_pool, IBDverse_E_MTAB_16999 = ibdverse_pool)
pool_summary <- list()
type_summary <- list()
sample_parts <- list()

for (dataset in names(pool_objects)) {
  pool <- pool_objects[[dataset]]
  if (!nrow(pool) || anyNA(pool$CellID) || anyDuplicated(pool$CellID)) stop("Invalid nonIBD cell identity: ", dataset, call. = FALSE)
  if (!setequal(unique(pool$CellType), supported[[dataset]])) stop("nonIBD supported parent types changed: ", dataset, call. = FALSE)

  region_lower <- tolower(trimws(pool$Region))
  pool[, Region3 := ifelse(region_lower %in% c("co", "colon"), "Colon", ifelse(region_lower == "rectum", "Rectum", ifelse(region_lower %in% c("ti", "sb", "terminal ileum", "ileum"), "Ileum", NA_character_)))]
  if (anyNA(pool$Region3)) stop("nonIBD Region3 mapping incomplete: ", dataset, call. = FALSE)
  pool[, Disease := "nonIBD"]

  sample_metadata <- unique(pool[, .(Dataset, DonorID, BiologicalSampleID, Disease, Inflammation, Region3)])
  if (sample_metadata[, .N, by = .(Dataset, BiologicalSampleID)][N > 1L, .N] > 0L) stop("nonIBD sample phenotype conflict: ", dataset, call. = FALSE)
  sample_parts[[dataset]] <- sample_metadata

  pool_file <- file.path(pool_dir, paste0(dataset, "_eligible_nonIBD_cells.parquet"))
  arrow::write_parquet(as.data.frame(pool), pool_file, compression = "zstd")
  pool_summary[[dataset]] <- data.frame(
    Dataset = dataset, EligibleCells = nrow(pool), Donors = uniqueN(pool$DonorID), BiologicalSamples = uniqueN(pool$BiologicalSampleID),
    ParentTypes = uniqueN(pool$CellType), Regions = paste(sort(unique(pool$Region3)), collapse = ";"),
    PoolFile = pool_file, PoolSHA256 = digest::digest(pool_file, algo = "sha256", file = TRUE, serialize = FALSE)
  )
  type_summary[[dataset]] <- pool[, .(EligibleCells = .N, Donors = uniqueN(DonorID), BiologicalSamples = uniqueN(BiologicalSampleID)), by = .(Dataset, CellType)]
  pool_objects[[dataset]] <- pool
  cat("NONIBD_POOL_FROZEN", dataset, nrow(pool), "cells", uniqueN(pool$DonorID), "donors", uniqueN(pool$CellType), "parent types\n")
}

# Step4. 固定每个 dataset × parent 的 recovery 细胞
selected_parts <- list()
sampling_parts <- list()
for (dataset in names(pool_objects)) {
  pool <- pool_objects[[dataset]]
  for (cell_type in sort(unique(pool$CellType))) {
    type_pool <- pool[CellType == cell_type, .(CellID, CellType, DonorID, BiologicalSampleID)]
    setorder(type_pool, CellID)
    set.seed(sampling_seed)
    selected_index <- if (nrow(type_pool) > max_cells) sample(seq_len(nrow(type_pool)), max_cells, replace = FALSE) else seq_len(nrow(type_pool))
    selected <- type_pool[selected_index]
    selected[, `:=`(Dataset = dataset, SelectionOrder = seq_len(.N), SamplingSeed = sampling_seed)]
    setcolorder(selected, c("Dataset", "CellType", "CellID", "DonorID", "BiologicalSampleID", "SelectionOrder", "SamplingSeed"))
    key <- paste(dataset, cell_type, sep = "::")
    selected_parts[[key]] <- selected
    sampling_parts[[key]] <- data.frame(
      Dataset = dataset, CellType = cell_type, EligibleCells = nrow(type_pool), SelectedCells = nrow(selected),
      SamplingApplied = nrow(type_pool) > max_cells, MaxCells = max_cells, SamplingSeed = sampling_seed,
      SelectedCellIDSHA256 = digest::digest(paste(selected$CellID, collapse = "\n"), algo = "sha256", serialize = FALSE)
    )
    cat("NONIBD_RECOVERY_CELLS_FROZEN", dataset, cell_type, nrow(selected), "/", nrow(type_pool), "\n")
  }
}

# Step5. 保存审计结果并确认与 IBD recovery 无细胞重叠
selected_cells <- rbindlist(selected_parts)
sampling_manifest <- rbindlist(sampling_parts)
sample_metadata_all <- rbindlist(sample_parts)
pool_summary_all <- rbindlist(pool_summary)
type_summary_all <- rbindlist(type_summary)
ibd_selected <- fread(file.path(output_root, "03_scrna_input_audit/recovery_selected_cells.tsv.gz"), select = c("Dataset", "CellID"), colClasses = "character")
overlap <- merge(selected_cells[, .(Dataset, CellID)], ibd_selected, by = c("Dataset", "CellID"))
if (nrow(overlap) > 0L) stop("nonIBD selected cells overlap the frozen IBD recovery set", call. = FALSE)

fwrite(pool_summary_all, file.path(audit_dir, "nonIBD_dataset_summary.tsv"), sep = "\t")
fwrite(type_summary_all, file.path(audit_dir, "nonIBD_parent_type_coverage.tsv"), sep = "\t")
fwrite(sample_metadata_all, file.path(audit_dir, "nonIBD_donor_sample_phenotype_freeze.tsv"), sep = "\t")
fwrite(selected_cells, file.path(audit_dir, "nonIBD_recovery_selected_cells.tsv.gz"), sep = "\t")
fwrite(sampling_manifest, file.path(audit_dir, "nonIBD_recovery_cell_sampling_manifest.tsv"), sep = "\t")

source_manifest <- data.table(Dataset = names(source_files), Path = unname(source_files))
source_manifest[, `:=`(Bytes = file.info(Path)$size, SHA256 = vapply(Path, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE))]
fwrite(source_manifest, file.path(audit_dir, "nonIBD_metadata_source_manifest.tsv"), sep = "\t")

summary <- list(
  status = "PASS_NONIBD_SCRNA_INPUT_AND_CELL_FREEZE", datasets = names(pool_objects), total_eligible_cells = sum(pool_summary_all$EligibleCells),
  total_selected_cells = nrow(selected_cells), parent_type_groups = nrow(sampling_manifest), max_cells_per_dataset_parent = max_cells,
  sampling_seed = sampling_seed, overlap_with_ibd_selected_cells = nrow(overlap), model_changed = FALSE, existing_ibd_results_changed = FALSE,
  biological_contrast = "nonIBD healthy donors versus frozen IBD recovery assignments after separate current-model recovery"
)
jsonlite::write_json(summary, file.path(audit_dir, "nonIBD_input_freeze.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo.txt"))
cat("NONIBD_INPUT_FREEZE_PASS", nrow(selected_cells), "selected cells in", nrow(sampling_manifest), "dataset-parent groups\n")
