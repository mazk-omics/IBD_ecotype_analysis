#!/usr/bin/env Rscript
# 本脚本合并已通过验收的 IBD 与 nonIBD 单细胞 recovery，冻结后续表型分析的逐细胞、样本和供者输入。
# 在服务器运行，只读取正式 EcoTyper 输出，结果写入 output/22_IBD_nonIBD_phenotype_input。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入和输出
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 33_build_IBD_nonIBD_phenotype_input.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

output_dir <- file.path(RUN_ROOT, "output/22_IBD_nonIBD_phenotype_input")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_specs <- rbindlist(list(
  data.table(
    AnalysisGroup = "IBD",
    Dataset = c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999"),
    JobsFile = file.path(RUN_ROOT, c(
      "output/04e_GSE214695_recovery_significance_100/GSE214695_recovery_significance_100_jobs.tsv",
      "output/05e_SCP1884_recovery_significance_100/SCP1884_recovery_significance_100_jobs.tsv",
      "output/06e_SCP259_recovery_significance_100/SCP259_recovery_significance_100_jobs.tsv",
      "output/07e_IBDverse_recovery_significance_100/IBDverse_E_MTAB_16999_recovery_significance_100_jobs.tsv"
    )),
    PassFile = file.path(RUN_ROOT, c(
      "output/04e_GSE214695_recovery_significance_100/GSE214695_recovery_significance_100_summary.txt",
      "output/05e_SCP1884_recovery_significance_100/SCP1884_recovery_significance_100_summary.txt",
      "output/06e_SCP259_recovery_significance_100/SCP259_recovery_significance_100_summary.txt",
      "output/07e_IBDverse_recovery_significance_100/IBDverse_E_MTAB_16999_recovery_significance_100_summary.txt"
    )),
    PassPattern = c(
      "STATUS=PASS_GSE214695_RECOVERY_SIGNIFICANCE_100",
      "STATUS=PASS_SCP1884_RECOVERY_SIGNIFICANCE_100",
      "STATUS=PASS_SCP259_RECOVERY_SIGNIFICANCE_100",
      "STATUS=PASS_IBDverse_E_MTAB_16999_RECOVERY_SIGNIFICANCE_100"
    )
  ),
  data.table(
    AnalysisGroup = "nonIBD",
    Dataset = c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999"),
    JobsFile = file.path(RUN_ROOT, "output/21_nonIBD_recovery_significance_100", c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999"), "recovery_significance_100_jobs.tsv"),
    PassFile = file.path(RUN_ROOT, "output/21_nonIBD_recovery_significance_100", c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999"), "recovery_significance_100_summary.txt"),
    PassPattern = "STATUS=PASS_NONIBD_RECOVERY_SIGNIFICANCE_100"
  )
))

if (any(!file.exists(dataset_specs$JobsFile)) || any(!file.exists(dataset_specs$PassFile))) stop("Missing recovery job manifest or PASS file", call. = FALSE)
pass_ok <- vapply(seq_len(nrow(dataset_specs)), function(i) any(readLines(dataset_specs$PassFile[[i]], warn = FALSE) == dataset_specs$PassPattern[[i]]), logical(1L))
if (any(!pass_ok)) stop("One or more IBD/nonIBD recovery significance runs did not pass", call. = FALSE)

# Step2. 读取正式分配并保留未分配细胞
cell_parts <- list()
manifest_parts <- list()
part_index <- 0L

for (spec_index in seq_len(nrow(dataset_specs))) {
  jobs <- fread(dataset_specs$JobsFile[[spec_index]], colClasses = "character")
  if (any(jobs$ExitStatus != "0" | jobs$SignificanceExitStatus != "0")) stop("A recovery job did not exit successfully", call. = FALSE)

  for (job_index in seq_len(nrow(jobs))) {
    annotation_file <- jobs$AnnotationFile[[job_index]]
    assignment_file <- file.path(jobs$ResultDir[[job_index]], "state_assignment.txt")
    z_file <- file.path(jobs$ResultDir[[job_index]], "recovery_z_scores.txt")
    input_files <- c(annotation_file, assignment_file, z_file)
    if (any(!file.exists(input_files))) stop("Missing annotation, formal assignment or recovery Z file", call. = FALSE)

    annotation <- fread(annotation_file, colClasses = "character")
    assignment <- fread(assignment_file, colClasses = "character")
    formal_z <- fread(z_file, colClasses = "character")
    required_annotation <- c("DonorID", "BiologicalSampleID", "ID", "CellType", "Disease", "Inflammation", "Region3")
    if (!all(required_annotation %in% names(annotation))) stop("Recovery annotation columns changed", call. = FALSE)
    if (!all(c("ID", "State") %in% names(assignment)) || !"State" %in% names(formal_z)) stop("Recovery output columns changed", call. = FALSE)
    if (anyDuplicated(annotation$ID) || anyDuplicated(assignment$ID)) stop("Duplicate cell IDs within a recovery job", call. = FALSE)
    if (!all(assignment$ID %in% annotation$ID)) stop("Assigned cells are absent from annotation", call. = FALSE)

    cells <- merge(annotation, assignment[, .(ID, State)], by = "ID", all.x = TRUE, sort = FALSE)
    cells <- cells[match(annotation$ID, ID)]
    cells[is.na(State) | !nzchar(State), State := "Unassigned"]
    cells[, `:=`(
      AnalysisGroup = dataset_specs$AnalysisGroup[[spec_index]], Dataset = dataset_specs$Dataset[[spec_index]],
      CellType = jobs$CellType[[job_index]], FormalStateUniverse = paste(sort(unique(formal_z$State)), collapse = ";")
    )]
    cells[, BroadPhenotype := fifelse(AnalysisGroup == "nonIBD", "nonIBD", "IBD")]
    setcolorder(cells, c("AnalysisGroup", "BroadPhenotype", "Dataset", "CellType", "ID", "State", "DonorID", "BiologicalSampleID", "Disease", "Inflammation", "Region3", "FormalStateUniverse"))

    part_index <- part_index + 1L
    cell_parts[[part_index]] <- cells
    manifest_parts[[part_index]] <- data.table(
      AnalysisGroup = dataset_specs$AnalysisGroup[[spec_index]], Dataset = dataset_specs$Dataset[[spec_index]], CellType = jobs$CellType[[job_index]],
      Input = c("RecoveryAnnotation", "FormalAssignment", "FormalRecoveryZ"), Path = input_files,
      Bytes = file.info(input_files)$size, SHA256 = vapply(input_files, digest::digest, character(1L), file = TRUE, algo = "sha256")
    )
  }
}

cells <- rbindlist(cell_parts, use.names = TRUE, fill = TRUE)
input_manifest <- rbindlist(manifest_parts, use.names = TRUE, fill = TRUE)
setorder(cells, AnalysisGroup, Dataset, CellType, ID)
setorder(input_manifest, AnalysisGroup, Dataset, CellType, Input)

if (anyDuplicated(cells[, .(Dataset, CellType, ID)])) stop("IBD and nonIBD inputs contain duplicate cells", call. = FALSE)
if (cells[AnalysisGroup == "IBD", any(BroadPhenotype != "IBD")] || cells[AnalysisGroup == "nonIBD", any(BroadPhenotype != "nonIBD")]) stop("Broad phenotype assignment failed", call. = FALSE)
if (cells[AnalysisGroup == "IBD", any(grepl("non[-_ ]?ibd|healthy|control|^hc$", Disease, ignore.case = TRUE))]) stop("Healthy/nonIBD cells entered the IBD input", call. = FALSE)
if (cells[AnalysisGroup == "nonIBD", any(!grepl("non[-_ ]?ibd|healthy|control|^hc$", Disease, ignore.case = TRUE))]) stop("Disease cells entered the nonIBD input", call. = FALSE)

state_universe_audit <- unique(cells[, .(AnalysisGroup, Dataset, CellType, FormalStateUniverse)])
state_universe_check <- state_universe_audit[, .(StateUniverses = uniqueN(FormalStateUniverse)), by = CellType]
if (any(state_universe_check$StateUniverses != 1L)) stop("Formal state universe differs between IBD and nonIBD inputs", call. = FALSE)

# Step3. 冻结样本和供者层分母
sample_meta <- unique(cells[, .(AnalysisGroup, BroadPhenotype, Dataset, CellType, BiologicalSampleID, DonorID, Disease, Inflammation, Region3)])
if (sample_meta[, anyDuplicated(paste(AnalysisGroup, Dataset, CellType, BiologicalSampleID))] != 0L) stop("One biological sample has conflicting donor or phenotype labels", call. = FALSE)

sample_parent <- cells[, .(
  InputCells = .N, AssignedCells = sum(State != "Unassigned"), UnassignedCells = sum(State == "Unassigned")
), by = .(AnalysisGroup, BroadPhenotype, Dataset, CellType, BiologicalSampleID, DonorID, Disease, Inflammation, Region3)]
sample_parent[, AssignmentRate := AssignedCells / InputCells]
setorder(sample_parent, AnalysisGroup, Dataset, CellType, BiologicalSampleID)

dataset_summary <- cells[, .(
  InputCells = .N, AssignedCells = sum(State != "Unassigned"), BiologicalSamples = uniqueN(BiologicalSampleID), Donors = uniqueN(DonorID), ParentTypes = uniqueN(CellType)
), by = .(AnalysisGroup, BroadPhenotype, Dataset)]
dataset_summary[, AssignmentRate := AssignedCells / InputCells]
setorder(dataset_summary, AnalysisGroup, Dataset)

cell_count_distribution <- sample_parent[, .(
  SampleParentPairs = .N, MinCells = as.numeric(min(InputCells)), Q10Cells = as.numeric(quantile(InputCells, 0.10)),
  Q25Cells = as.numeric(quantile(InputCells, 0.25)), MedianCells = as.numeric(median(InputCells)), Q75Cells = as.numeric(quantile(InputCells, 0.75)),
  Q90Cells = as.numeric(quantile(InputCells, 0.90)), MaxCells = as.numeric(max(InputCells))
), by = .(AnalysisGroup, Dataset, CellType)]
setorder(cell_count_distribution, AnalysisGroup, Dataset, CellType)

thresholds <- c(5L, 10L, 20L, 30L, 50L, 100L)
threshold_retention <- rbindlist(lapply(thresholds, function(threshold) {
  sample_parent[, .(
    SampleParentPairs = .N, RetainedPairs = sum(InputCells >= threshold), RetainedInputCells = sum(InputCells[InputCells >= threshold]),
    TotalInputCells = sum(InputCells), RetainedSamples = uniqueN(BiologicalSampleID[InputCells >= threshold]), TotalSamples = uniqueN(BiologicalSampleID)
  ), by = .(AnalysisGroup, Dataset)][, MinParentCells := threshold]
}))
threshold_retention[, `:=`(RetainedPairFraction = RetainedPairs / SampleParentPairs, RetainedCellFraction = RetainedInputCells / TotalInputCells)]
setcolorder(threshold_retention, c("MinParentCells", "AnalysisGroup", "Dataset", "SampleParentPairs", "RetainedPairs", "RetainedPairFraction", "RetainedInputCells", "TotalInputCells", "RetainedCellFraction", "RetainedSamples", "TotalSamples"))
setorder(threshold_retention, MinParentCells, AnalysisGroup, Dataset)

donors <- unique(cells[, .(AnalysisGroup, Dataset, DonorID, BroadPhenotype)])
dataset_pairs <- CJ(DatasetA = sort(unique(donors$Dataset)), DatasetB = sort(unique(donors$Dataset)))[DatasetA < DatasetB]
donor_overlap <- dataset_pairs[, {
  overlap <- intersect(donors[Dataset == DatasetA, DonorID], donors[Dataset == DatasetB, DonorID])
  .(OverlappingDonors = length(overlap), DonorIDs = paste(sort(overlap), collapse = ";"))
}, by = .(DatasetA, DatasetB)]

cat("Unified cells:", nrow(cells), "\n")
print(dataset_summary)
print(threshold_retention[MinParentCells %in% c(10L, 20L, 30L)])
print(donor_overlap)

# Step4. 保存统一分析输入和验收记录
fwrite(cells, file.path(output_dir, "IBD_nonIBD_formal_cell_assignments.tsv.gz"), sep = "\t", quote = FALSE, na = "")
fwrite(sample_meta, file.path(output_dir, "IBD_nonIBD_sample_metadata.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(sample_parent, file.path(output_dir, "IBD_nonIBD_sample_parent_denominators.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(dataset_summary, file.path(output_dir, "IBD_nonIBD_dataset_summary.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(cell_count_distribution, file.path(output_dir, "IBD_nonIBD_sample_parent_cell_count_distribution.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(threshold_retention, file.path(output_dir, "IBD_nonIBD_minimum_cell_threshold_retention.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(donor_overlap, file.path(output_dir, "IBD_nonIBD_cross_dataset_donor_overlap.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(state_universe_audit, file.path(output_dir, "IBD_nonIBD_state_universe_audit.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(input_manifest, file.path(output_dir, "IBD_nonIBD_phenotype_input_manifest.tsv"), sep = "\t", quote = FALSE, na = "")

completion <- list(
  status = "PASS_IBD_NONIBD_PHENOTYPE_INPUT_AUDIT", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = uniqueN(cells$Dataset), analysis_groups = uniqueN(cells$AnalysisGroup), parent_types = uniqueN(cells$CellType),
  input_cells = nrow(cells), assigned_cells = sum(cells$State != "Unassigned"),
  biological_samples = uniqueN(paste(cells$Dataset, cells$BiologicalSampleID)), donors = uniqueN(cells$DonorID),
  state_universes_consistent = all(state_universe_check$StateUniverses == 1L),
  threshold_status = "AUDITED_NOT_FROZEN"
)
jsonlite::write_json(completion, file.path(output_dir, "IBD_nonIBD_phenotype_input_complete.json"), pretty = TRUE, auto_unbox = TRUE)
cat("PASS_IBD_NONIBD_PHENOTYPE_INPUT_AUDIT\n")
