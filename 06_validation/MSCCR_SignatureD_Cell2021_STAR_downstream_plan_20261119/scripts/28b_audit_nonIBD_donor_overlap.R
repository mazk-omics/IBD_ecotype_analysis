#!/usr/bin/env Rscript
# 本脚本冻结四套 nonIBD 队列的供者重叠；用于后续合并时避免把重复健康供者当作独立证据。
# 只读取 28 脚本生成的样本表，不修改 recovery 输入。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取 nonIBD 供者与样本
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 28b_audit_nonIBD_donor_overlap.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing R package: data.table", call. = FALSE)
library(data.table)

audit_dir <- file.path(output_root, "18_nonIBD_scrna_input_audit")
sample_metadata <- fread(file.path(audit_dir, "nonIBD_donor_sample_phenotype_freeze.tsv"), colClasses = "character")
donors <- unique(sample_metadata[, .(Dataset, DonorID)])
datasets <- sort(unique(donors$Dataset))

# Step2. 逐队列对核对相同供者 ID
dataset_pairs <- CJ(DatasetA = datasets, DatasetB = datasets)[DatasetA < DatasetB]
overlap_parts <- vector("list", nrow(dataset_pairs))
summary_parts <- vector("list", nrow(dataset_pairs))

for (pair_index in seq_len(nrow(dataset_pairs))) {
  dataset_a <- dataset_pairs$DatasetA[[pair_index]]
  dataset_b <- dataset_pairs$DatasetB[[pair_index]]
  donor_overlap <- intersect(donors[Dataset == dataset_a, DonorID], donors[Dataset == dataset_b, DonorID])
  overlap_parts[[pair_index]] <- data.table(DatasetA = dataset_a, DatasetB = dataset_b, DonorID = donor_overlap)
  summary_parts[[pair_index]] <- data.table(DatasetA = dataset_a, DatasetB = dataset_b, OverlappingDonors = length(donor_overlap))
}

overlap_detail <- rbindlist(overlap_parts)
overlap_summary <- rbindlist(summary_parts)
unexpected <- overlap_summary[OverlappingDonors > 0L & !(DatasetA == "SCP1884" & DatasetB == "SCP259")]
if (nrow(unexpected) > 0L) stop("Unexpected nonIBD donor overlap across datasets", call. = FALSE)
if (overlap_summary[DatasetA == "SCP1884" & DatasetB == "SCP259", OverlappingDonors] != 12L) stop("SCP1884/SCP259 nonIBD donor overlap changed", call. = FALSE)

# Step3. 保存后续汇总的独立性边界
fwrite(overlap_summary, file.path(audit_dir, "nonIBD_cross_dataset_donor_overlap_summary.tsv"), sep = "\t")
fwrite(overlap_detail, file.path(audit_dir, "nonIBD_cross_dataset_donor_overlap_detail.tsv"), sep = "\t")
writeLines(
  c(
    "STATUS=PASS_NONIBD_DONOR_OVERLAP_AUDIT",
    "SCP1884_SCP259_OVERLAPPING_HEALTHY_DONORS=12",
    "OTHER_CROSS_DATASET_OVERLAPS=0",
    "META_ANALYSIS_RULE=SCP1884 and SCP259 are not independent for the 12 shared healthy donors; report cohort-specific results and donor-deduplicated sensitivity"
  ),
  file.path(audit_dir, "nonIBD_donor_overlap_audit.txt")
)
cat("NONIBD_DONOR_OVERLAP_AUDIT_PASS 12 shared healthy donors in SCP1884/SCP259; all other pairs 0\n")
