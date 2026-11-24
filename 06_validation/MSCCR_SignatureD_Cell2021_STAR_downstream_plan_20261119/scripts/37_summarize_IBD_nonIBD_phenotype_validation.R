#!/usr/bin/env Rscript
# 本脚本汇总 IBD/nonIBD 的四条 cell-state phenotype validation 轨道，并显式保留全部阳性和阴性结果。
# 在服务器 RUN_ROOT 运行；结果写入 output/26_IBD_nonIBD_phenotype_validation_summary。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取四条验证轨道
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 37_summarize_IBD_nonIBD_phenotype_validation.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

permutation_dir <- file.path(RUN_ROOT, "output/23_IBD_nonIBD_known_phenotype_enrichment_Cell2021")
representation_dir <- file.path(RUN_ROOT, "output/24_IBD_nonIBD_state_representation")
abundance_dir <- file.path(RUN_ROOT, "output/25_IBD_nonIBD_sample_donor_abundance")
output_dir <- file.path(RUN_ROOT, "output/26_IBD_nonIBD_phenotype_validation_summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

permutation_files <- list.files(permutation_dir, pattern = "^state_known_phenotype_enrichment_Cell2021\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(permutation_files) != 4L) stop("Expected four Cell2021 permutation result files", call. = FALSE)
permutation <- rbindlist(lapply(permutation_files, fread), use.names = TRUE, fill = TRUE)

representation_file <- file.path(representation_dir, "sample_dominant_state_representation_fisher.tsv")
representation_pass_file <- file.path(representation_dir, "cell_state_representation_complete.json")
abundance_file <- file.path(abundance_dir, "sample_donor_state_abundance_tests.tsv")
abundance_pass_file <- file.path(abundance_dir, "sample_donor_state_abundance_complete.json")
if (!all(file.exists(c(representation_file, representation_pass_file, abundance_file, abundance_pass_file)))) stop("Scripts 35-36 have not completed", call. = FALSE)

representation_pass <- jsonlite::read_json(representation_pass_file, simplifyVector = TRUE)
abundance_pass <- jsonlite::read_json(abundance_pass_file, simplifyVector = TRUE)
if (representation_pass$status != "PASS_IBD_NONIBD_CELL_STATE_REPRESENTATION") stop("Representation audit did not pass", call. = FALSE)
if (abundance_pass$status != "PASS_IBD_NONIBD_SAMPLE_DONOR_STATE_ABUNDANCE") stop("Sample/donor abundance audit did not pass", call. = FALSE)

representation <- fread(representation_file)
abundance <- fread(abundance_file)

cat("Permutation rows:", nrow(permutation), "\n")
cat("Representation rows:", nrow(representation), "\n")
cat("Sample/donor rows:", nrow(abundance), "\n")

# Step2. 汇总 Cell 2021 式逐细胞置换与 dominant-state representation
permutation_summary <- permutation[, .(
  EvaluableStates = uniqueN(CellState[Status == "EVALUABLE"]),
  PositiveZgt1.96States = uniqueN(CellState[Status == "EVALUABLE" & PositiveZgt1.96 == TRUE]),
  SignificantBHStates = uniqueN(CellState[Status == "EVALUABLE" & Significant_BH_sensitivity == TRUE])
), by = .(Dataset, Axis)]
setorder(permutation_summary, Dataset, Axis)

representation_summary <- representation[, .(
  EvaluableStates = uniqueN(CellState[Status == "EVALUABLE"]),
  SignificantFDRStates = uniqueN(CellState[Status == "EVALUABLE" & SignificantEnrichment == TRUE]),
  EligibleSampleParents = max(EligibleSampleParents, na.rm = TRUE),
  Donors = max(Donors, na.rm = TRUE)
), by = .(MinParentCells, Scope, DatasetScope, Axis)]
setorder(representation_summary, MinParentCells, Scope, DatasetScope, Axis)

fwrite(permutation_summary, file.path(output_dir, "cell2021_permutation_summary.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(representation_summary, file.path(output_dir, "dominant_state_representation_summary.tsv"), sep = "\t", quote = FALSE, na = "")

# Step3. 配对 sample 与 donor 结果
abundance_paired <- dcast(
  abundance,
  MinParentCells + Dataset + CellType + State + CellState + Axis + Contrast + Target + Reference ~ Unit,
  value.var = c("TargetUnits", "ReferenceUnits", "DeltaMedian", "RankBiserial", "P_two_sided", "Status", "FDR_BH_within_dataset_contrast", "SignificantFDR", "Direction")
)

required_paired <- c("Status_Sample", "Status_Donor", "DeltaMedian_Sample", "DeltaMedian_Donor", "SignificantFDR_Sample", "SignificantFDR_Donor")
if (!all(required_paired %in% names(abundance_paired))) stop("Sample/donor pairing columns changed", call. = FALSE)

abundance_paired[, EvaluableBoth := Status_Sample == "EVALUABLE" & Status_Donor == "EVALUABLE"]
abundance_paired[, SameEffectDirection := is.finite(DeltaMedian_Sample) & is.finite(DeltaMedian_Donor) & sign(DeltaMedian_Sample) == sign(DeltaMedian_Donor) & sign(DeltaMedian_Sample) != 0]
abundance_paired[, RobustBothFDR := EvaluableBoth & SignificantFDR_Sample == TRUE & SignificantFDR_Donor == TRUE & SameEffectDirection]
abundance_paired[, RobustEnrichedLevel := fifelse(!RobustBothFDR, NA_character_, fifelse(DeltaMedian_Sample > 0, Target, Reference))]
setorder(abundance_paired, MinParentCells, Axis, Contrast, Dataset, CellType, State)
fwrite(abundance_paired, file.path(output_dir, "sample_donor_all_concordance.tsv"), sep = "\t", quote = FALSE, na = "")

main_concordance <- abundance_paired[MinParentCells == 20L & Axis == "BroadPhenotype" & Contrast == "IBD vs nonIBD"]
main_robust <- main_concordance[RobustBothFDR == TRUE]
setorder(main_robust, Dataset, CellType, State)
fwrite(main_concordance, file.path(output_dir, "sample_donor_main_broad_phenotype_concordance.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(main_robust, file.path(output_dir, "sample_donor_robust_state_dataset_pairs.tsv"), sep = "\t", quote = FALSE, na = "")

# Step4. 统计跨数据集复现和方向一致性
cross_dataset <- abundance_paired[, .(
  EvaluableDatasets = sum(EvaluableBoth, na.rm = TRUE),
  SampleSignificantDatasets = sum(SignificantFDR_Sample == TRUE, na.rm = TRUE),
  DonorSignificantDatasets = sum(SignificantFDR_Donor == TRUE, na.rm = TRUE),
  RobustDatasets = sum(RobustBothFDR, na.rm = TRUE),
  TargetHigherRobustDatasets = sum(RobustBothFDR & DeltaMedian_Sample > 0, na.rm = TRUE),
  ReferenceHigherRobustDatasets = sum(RobustBothFDR & DeltaMedian_Sample < 0, na.rm = TRUE),
  RobustDirectionConsistent = sum(RobustBothFDR, na.rm = TRUE) > 0L && !(any(RobustBothFDR & DeltaMedian_Sample > 0) && any(RobustBothFDR & DeltaMedian_Sample < 0)),
  RobustDatasetList = paste(Dataset[RobustBothFDR], collapse = ";")
), by = .(MinParentCells, Axis, Contrast, CellType, State, CellState)]
setorder(cross_dataset, MinParentCells, Axis, Contrast, -RobustDatasets, CellType, State)
fwrite(cross_dataset, file.path(output_dir, "sample_donor_cross_dataset_replication.tsv"), sep = "\t", quote = FALSE, na = "")

# Step5. 对稳健信号核对其他两条方法
permutation_main <- permutation[Axis == "BroadPhenotype" & Level %in% c("IBD", "nonIBD")]
permutation_support <- permutation_main[, .(
  PermutationRaw_IBD = any(Level == "IBD" & PositiveZgt1.96 == TRUE),
  PermutationRaw_nonIBD = any(Level == "nonIBD" & PositiveZgt1.96 == TRUE),
  PermutationBH_IBD = any(Level == "IBD" & Significant_BH_sensitivity == TRUE),
  PermutationBH_nonIBD = any(Level == "nonIBD" & Significant_BH_sensitivity == TRUE),
  Z_IBD = Z[Level == "IBD"][1L],
  Z_nonIBD = Z[Level == "nonIBD"][1L]
), by = .(Dataset, CellType, State, CellState)]

representation_main <- representation[MinParentCells == 20L & Scope == "WITHIN_DATASET" & Axis == "BroadPhenotype" & Level %in% c("IBD", "nonIBD")]
representation_support <- representation_main[, .(
  RepresentationFDR_IBD = any(Level == "IBD" & SignificantEnrichment == TRUE),
  RepresentationFDR_nonIBD = any(Level == "nonIBD" & SignificantEnrichment == TRUE),
  OddsRatio_IBD = OddsRatio[Level == "IBD"][1L],
  OddsRatio_nonIBD = OddsRatio[Level == "nonIBD"][1L]
), by = .(Dataset = DatasetScope, CellType, State, CellState)]

method_concordance <- merge(main_concordance, permutation_support, by = c("Dataset", "CellType", "State", "CellState"), all.x = TRUE, sort = FALSE)
method_concordance <- merge(method_concordance, representation_support, by = c("Dataset", "CellType", "State", "CellState"), all.x = TRUE, sort = FALSE)
logical_columns <- c("PermutationRaw_IBD", "PermutationRaw_nonIBD", "PermutationBH_IBD", "PermutationBH_nonIBD", "RepresentationFDR_IBD", "RepresentationFDR_nonIBD")
for (column_name in logical_columns) set(method_concordance, which(is.na(method_concordance[[column_name]])), column_name, FALSE)
setorder(method_concordance, Dataset, CellType, State)
fwrite(method_concordance, file.path(output_dir, "main_broad_phenotype_method_concordance.tsv"), sep = "\t", quote = FALSE, na = "")

robust_methods <- method_concordance[RobustBothFDR == TRUE]
method_support_long <- rbindlist(list(
  robust_methods[, .(Dataset, CellType, State, CellState, RobustEnrichedLevel, Method = "Cell permutation Z>1.96", SupportSameDirection = fifelse(RobustEnrichedLevel == "IBD", PermutationRaw_IBD, PermutationRaw_nonIBD), SupportOppositeDirection = fifelse(RobustEnrichedLevel == "IBD", PermutationRaw_nonIBD, PermutationRaw_IBD))],
  robust_methods[, .(Dataset, CellType, State, CellState, RobustEnrichedLevel, Method = "Cell permutation BH<0.05", SupportSameDirection = fifelse(RobustEnrichedLevel == "IBD", PermutationBH_IBD, PermutationBH_nonIBD), SupportOppositeDirection = fifelse(RobustEnrichedLevel == "IBD", PermutationBH_nonIBD, PermutationBH_IBD))],
  robust_methods[, .(Dataset, CellType, State, CellState, RobustEnrichedLevel, Method = "Dominant-state Fisher", SupportSameDirection = fifelse(RobustEnrichedLevel == "IBD", RepresentationFDR_IBD, RepresentationFDR_nonIBD), SupportOppositeDirection = fifelse(RobustEnrichedLevel == "IBD", RepresentationFDR_nonIBD, RepresentationFDR_IBD))],
  robust_methods[, .(Dataset, CellType, State, CellState, RobustEnrichedLevel, Method = "Sample Wilcoxon", SupportSameDirection = SignificantFDR_Sample == TRUE, SupportOppositeDirection = FALSE)],
  robust_methods[, .(Dataset, CellType, State, CellState, RobustEnrichedLevel, Method = "Donor Wilcoxon", SupportSameDirection = SignificantFDR_Donor == TRUE, SupportOppositeDirection = FALSE)]
), use.names = TRUE)
method_support_long[, SupportClass := fifelse(SupportSameDirection, "Same direction", fifelse(SupportOppositeDirection, "Opposite direction", "Not significant"))]
fwrite(method_support_long, file.path(output_dir, "robust_state_method_support_long.tsv"), sep = "\t", quote = FALSE, na = "")

# Step6. 给作图准备可比较的显著比例
permutation_counts_raw <- permutation_main[, .(
  EvaluableStates = uniqueN(CellState[Status == "EVALUABLE"]),
  SignificantStates = uniqueN(CellState[Status == "EVALUABLE" & PositiveZgt1.96 == TRUE])
), by = Dataset][, Method := "Cell permutation Z>1.96"]
permutation_counts_bh <- permutation_main[, .(
  EvaluableStates = uniqueN(CellState[Status == "EVALUABLE"]),
  SignificantStates = uniqueN(CellState[Status == "EVALUABLE" & Significant_BH_sensitivity == TRUE])
), by = Dataset][, Method := "Cell permutation BH<0.05"]
representation_counts <- representation_main[, .(
  EvaluableStates = uniqueN(CellState[Status == "EVALUABLE"]),
  SignificantStates = uniqueN(CellState[Status == "EVALUABLE" & SignificantEnrichment == TRUE])
), by = .(Dataset = DatasetScope)][, Method := "Dominant-state Fisher"]
sample_counts <- main_concordance[, .(EvaluableStates = sum(Status_Sample == "EVALUABLE"), SignificantStates = sum(SignificantFDR_Sample == TRUE, na.rm = TRUE)), by = Dataset][, Method := "Sample Wilcoxon"]
donor_counts <- main_concordance[, .(EvaluableStates = sum(Status_Donor == "EVALUABLE"), SignificantStates = sum(SignificantFDR_Donor == TRUE, na.rm = TRUE)), by = Dataset][, Method := "Donor Wilcoxon"]
robust_counts <- main_concordance[, .(EvaluableStates = sum(EvaluableBoth), SignificantStates = sum(RobustBothFDR, na.rm = TRUE)), by = Dataset][, Method := "Sample + donor concordant"]

method_counts <- rbindlist(list(permutation_counts_raw, permutation_counts_bh, representation_counts, sample_counts, donor_counts, robust_counts), use.names = TRUE, fill = TRUE)
method_counts[, SignificantPercent := 100 * SignificantStates / EvaluableStates]
setcolorder(method_counts, c("Dataset", "Method", "SignificantStates", "EvaluableStates", "SignificantPercent"))
setorder(method_counts, Dataset, Method)
fwrite(method_counts, file.path(output_dir, "main_broad_phenotype_method_significant_counts.tsv"), sep = "\t", quote = FALSE, na = "")

# Step7. 完成审计
completion <- list(
  status = "PASS_IBD_NONIBD_PHENOTYPE_VALIDATION_SUMMARY",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = uniqueN(permutation$Dataset), permutation_rows = nrow(permutation), representation_rows = nrow(representation), abundance_rows = nrow(abundance),
  main_broad_phenotype_robust_pairs = nrow(main_robust), robust_cell_states = uniqueN(main_robust$CellState),
  robust_pairs_supported_by_cell_permutation_raw = sum(method_support_long[Method == "Cell permutation Z>1.96", SupportSameDirection]),
  robust_pairs_supported_by_representation_fdr = sum(method_support_long[Method == "Dominant-state Fisher", SupportSameDirection]),
  interpretation_limit = "WITHIN_DATASET_ASSOCIATION_NOT_UNIVERSAL_IBD_EFFECT"
)
jsonlite::write_json(completion, file.path(output_dir, "IBD_nonIBD_phenotype_validation_summary_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(main_robust[, .(Dataset, CellState, RobustEnrichedLevel, DeltaMedian_Sample, DeltaMedian_Donor, FDR_BH_within_dataset_contrast_Sample, FDR_BH_within_dataset_contrast_Donor)])
print(completion)
cat("PASS_IBD_NONIBD_PHENOTYPE_VALIDATION_SUMMARY\n")
