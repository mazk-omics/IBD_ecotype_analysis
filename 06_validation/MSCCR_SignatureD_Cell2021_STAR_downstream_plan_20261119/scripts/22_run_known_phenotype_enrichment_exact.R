#!/usr/bin/env Rscript
# 本脚本按 Cell 2021 原文轨道，在同数据集、同父类型的正式 assigned cells 内打乱 state labels 1,000 次。
# Disease、Inflammation 和 Region3 分开检验；Unassigned 只进入分母审计，不进入置换。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入和参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 22_run_known_phenotype_enrichment_exact.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

permutations <- 1000L
base_seed <- 20261123L
z_threshold <- 1.96
axes <- c("Disease", "Inflammation", "Region3")
output_dir <- file.path(RUN_ROOT, "output/11_known_phenotype_enrichment_exact")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_specs <- data.table(
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
)

if (any(!file.exists(dataset_specs$JobsFile)) || any(!file.exists(dataset_specs$PassFile))) stop("Missing recovery job manifest or PASS file", call. = FALSE)
pass_ok <- vapply(seq_len(nrow(dataset_specs)), function(i) any(readLines(dataset_specs$PassFile[[i]], warn = FALSE) == dataset_specs$PassPattern[[i]]), logical(1L))
if (any(!pass_ok)) stop("A recovery significance run did not pass", call. = FALSE)

# Step2. 冻结逐细胞正式分配和分母
cell_parts <- list()
denominator_parts <- list()
input_manifest_parts <- list()
part_index <- 0L

for (dataset_index in seq_len(nrow(dataset_specs))) {
  jobs <- fread(dataset_specs$JobsFile[[dataset_index]], colClasses = "character")
  if (any(jobs$ExitStatus != "0" | jobs$SignificanceExitStatus != "0")) stop("A recovery job did not exit successfully", call. = FALSE)

  for (job_index in seq_len(nrow(jobs))) {
    result_dir <- jobs$ResultDir[[job_index]]
    assignment_file <- file.path(result_dir, "state_assignment.txt")
    z_file <- file.path(result_dir, "recovery_z_scores.txt")
    annotation_file <- jobs$AnnotationFile[[job_index]]
    required_files <- c(assignment_file, z_file, annotation_file)
    if (any(!file.exists(required_files))) stop("Missing formal assignment, Z score or annotation file", call. = FALSE)

    assignment <- fread(assignment_file, colClasses = "character")
    formal_z <- fread(z_file, colClasses = "character")
    annotation <- fread(annotation_file, colClasses = "character")
    required_annotation <- c("ID", "CellType", "DonorID", "BiologicalSampleID", axes)
    if (!all(required_annotation %in% names(annotation))) stop("Recovery annotation columns changed", call. = FALSE)
    if (!all(c("ID", "State") %in% names(assignment)) || !"State" %in% names(formal_z)) stop("Recovery output columns changed", call. = FALSE)
    if (anyDuplicated(annotation$ID) || anyDuplicated(assignment$ID)) stop("Duplicate cell IDs in recovery inputs", call. = FALSE)
    if (!all(assignment$ID %in% annotation$ID)) stop("Assigned cells are absent from annotation", call. = FALSE)

    cells <- merge(annotation, assignment[, .(ID, State)], by = "ID", all.x = TRUE, sort = FALSE)
    cells <- cells[match(annotation$ID, ID)]
    cells[is.na(State) | !nzchar(State), State := "Unassigned"]
    cells[, Dataset := dataset_specs$Dataset[[dataset_index]]]
    cells[, CellType := jobs$CellType[[job_index]]]
    cells[, FormalStateUniverse := paste(sort(unique(formal_z$State)), collapse = ";")]
    setcolorder(cells, c("Dataset", "CellType", "ID", "State", "DonorID", "BiologicalSampleID", axes, "FormalStateUniverse"))

    part_index <- part_index + 1L
    cell_parts[[part_index]] <- cells
    input_manifest_parts[[part_index]] <- data.table(
      Dataset = dataset_specs$Dataset[[dataset_index]], CellType = jobs$CellType[[job_index]],
      Input = c("FormalAssignment", "FormalRecoveryZ", "RecoveryAnnotation"), Path = required_files,
      Bytes = file.info(required_files)$size,
      SHA256 = vapply(required_files, digest::digest, character(1L), file = TRUE, algo = "sha256")
    )

    for (axis_name in axes) {
      levels_axis <- sort(unique(cells[[axis_name]][!is.na(cells[[axis_name]]) & nzchar(cells[[axis_name]])]))
      for (level_name in levels_axis) {
        level_cells <- cells[get(axis_name) == level_name]
        denominator_parts[[length(denominator_parts) + 1L]] <- data.table(
          Dataset = dataset_specs$Dataset[[dataset_index]], CellType = jobs$CellType[[job_index]],
          Axis = axis_name, Level = level_name, InputCells = nrow(level_cells),
          AssignedCells = sum(level_cells$State != "Unassigned"), UnassignedCells = sum(level_cells$State == "Unassigned"),
          BiologicalSamples = uniqueN(level_cells$BiologicalSampleID), Donors = uniqueN(level_cells$DonorID)
        )
      }
    }
  }
}

cells_all <- rbindlist(cell_parts, use.names = TRUE, fill = TRUE)
denominators <- rbindlist(denominator_parts, use.names = TRUE, fill = TRUE)
input_manifest <- rbindlist(input_manifest_parts, use.names = TRUE, fill = TRUE)
setorder(cells_all, Dataset, CellType, ID)
setorder(denominators, Dataset, CellType, Axis, Level)
fwrite(cells_all, file.path(output_dir, "formal_cell_assignments_with_phenotypes.tsv.gz"), sep = "\t", quote = FALSE)
fwrite(denominators, file.path(output_dir, "phenotype_denominators.tsv"), sep = "\t", quote = FALSE)
fwrite(input_manifest, file.path(output_dir, "phenotype_enrichment_input_manifest.tsv"), sep = "\t", quote = FALSE)

# Step3. 原文式逐细胞标签置换
result_parts <- list()
result_index <- 0L

job_keys <- unique(cells_all[, .(Dataset, CellType, FormalStateUniverse)])
setorder(job_keys, Dataset, CellType)
for (job_index in seq_len(nrow(job_keys))) {
  dataset_name <- job_keys$Dataset[[job_index]]
  cell_type <- job_keys$CellType[[job_index]]
  state_universe <- strsplit(job_keys$FormalStateUniverse[[job_index]], ";", fixed = TRUE)[[1L]]
  parent_cells <- cells_all[Dataset == dataset_name & CellType == cell_type]

  for (axis_index in seq_along(axes)) {
    axis_name <- axes[[axis_index]]
    observations <- parent_cells[State != "Unassigned" & !is.na(get(axis_name)) & nzchar(get(axis_name))]
    axis_levels <- sort(unique(observations[[axis_name]]))
    seed <- base_seed + job_index * 100L + axis_index

    if (!nrow(observations) || length(axis_levels) < 2L) {
      grid <- CJ(Level = axis_levels, State = state_universe, sorted = FALSE)
      if (!nrow(grid)) grid <- data.table(Level = "", State = state_universe)
      grid[, `:=`(
        Dataset = dataset_name, CellType = cell_type, Axis = axis_name,
        ObservedAssignedCells = 0L, NullMean = NA_real_, NullSD = NA_real_, NullQ025 = NA_real_, NullQ975 = NA_real_,
        Z = NA_real_, PositiveZgt1.96 = FALSE, P_one_sided = NA_real_, EmpiricalP_one_sided = NA_real_,
        Permutations = permutations, Seed = seed, Status = "NOT_EVALUABLE_LT2_LEVELS"
      )]
      result_index <- result_index + 1L
      result_parts[[result_index]] <- grid
      next
    }

    labels <- observations$State
    phenotype <- observations[[axis_name]]
    observed_table <- table(factor(labels, levels = state_universe), factor(phenotype, levels = axis_levels))
    null_counts <- matrix(0L, nrow = permutations, ncol = length(state_universe) * length(axis_levels))
    set.seed(seed)
    for (permutation_index in seq_len(permutations)) {
      shuffled <- sample(labels, length(labels), replace = FALSE)
      permuted_table <- table(factor(shuffled, levels = state_universe), factor(phenotype, levels = axis_levels))
      null_counts[permutation_index, ] <- as.vector(permuted_table)
    }

    grid <- data.table(
      Level = rep(axis_levels, each = length(state_universe)),
      State = rep(state_universe, times = length(axis_levels)),
      ObservedAssignedCells = as.integer(as.vector(observed_table)),
      NullMean = colMeans(null_counts), NullSD = apply(null_counts, 2L, sd),
      NullQ025 = apply(null_counts, 2L, quantile, probs = 0.025, names = FALSE),
      NullQ975 = apply(null_counts, 2L, quantile, probs = 0.975, names = FALSE)
    )
    grid[, Z := fifelse(NullSD > 0, (ObservedAssignedCells - NullMean) / NullSD, NA_real_)]
    grid[, PositiveZgt1.96 := is.finite(Z) & Z > z_threshold]
    grid[, P_one_sided := fifelse(is.finite(Z), pnorm(Z, lower.tail = FALSE), NA_real_)]
    observed_matrix <- matrix(grid$ObservedAssignedCells, nrow = permutations, ncol = nrow(grid), byrow = TRUE)
    grid[, EmpiricalP_one_sided := (1 + colSums(null_counts >= observed_matrix)) / (permutations + 1)]
    grid[, `:=`(
      Dataset = dataset_name, CellType = cell_type, Axis = axis_name,
      Permutations = permutations, Seed = seed,
      Status = fifelse(is.finite(Z), "EVALUABLE", "NOT_EVALUABLE_ZERO_NULL_SD")
    )]
    result_index <- result_index + 1L
    result_parts[[result_index]] <- grid
  }
}

results <- rbindlist(result_parts, use.names = TRUE, fill = TRUE)
results[, CellState := paste(CellType, State, sep = "_")]
results[, FDR_BH_within_dataset_axis := p.adjust(P_one_sided, method = "BH"), by = .(Dataset, Axis)]
results[, Significant_BH_sensitivity := is.finite(Z) & Z > 0 & FDR_BH_within_dataset_axis < 0.05]
setcolorder(results, c(
  "Dataset", "CellType", "State", "CellState", "Axis", "Level", "ObservedAssignedCells",
  "NullMean", "NullSD", "NullQ025", "NullQ975", "Z", "PositiveZgt1.96", "P_one_sided",
  "EmpiricalP_one_sided", "FDR_BH_within_dataset_axis", "Significant_BH_sensitivity",
  "Permutations", "Seed", "Status"
))
setorder(results, Dataset, CellType, Axis, Level, State)
fwrite(results, file.path(output_dir, "state_known_phenotype_enrichment_exact.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(results[PositiveZgt1.96 == TRUE], file.path(output_dir, "state_known_phenotype_enrichment_exact_positive.tsv"), sep = "\t", quote = FALSE, na = "")

# Step4. 输出完成记录
completion <- list(
  status = "PASS_KNOWN_PHENOTYPE_ENRICHMENT_EXACT",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), datasets = uniqueN(cells_all$Dataset),
  parent_jobs = nrow(job_keys), input_cells = nrow(cells_all), assigned_cells = sum(cells_all$State != "Unassigned"),
  permutations = permutations, axes = axes, evaluable_rows = sum(results$Status == "EVALUABLE"),
  positive_z_rows = sum(results$PositiveZgt1.96), method = "CELL2021_EXACT_CELL_LABEL_PERMUTATION"
)
jsonlite::write_json(completion, file.path(output_dir, "known_phenotype_enrichment_exact_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_KNOWN_PHENOTYPE_ENRICHMENT_EXACT\n")
