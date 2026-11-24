#!/usr/bin/env Rscript
# 本脚本按 Cell 2021“known phenotype enrichment”算法，在一个数据集、同一父类型内打乱正式 state labels 1,000 次。
# BroadPhenotype、Disease、Inflammation 和 Region3 是 IBD 场景适配；结果写入 output/23_IBD_nonIBD_known_phenotype_enrichment_Cell2021。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取统一验收输入
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 34_run_IBD_nonIBD_known_phenotype_enrichment_Cell2021.R RUN_ROOT DATASET", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)
dataset <- args[[2L]]

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

datasets <- c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999")
if (!dataset %in% datasets) stop("Unknown dataset: ", dataset, call. = FALSE)

input_dir <- file.path(RUN_ROOT, "output/22_IBD_nonIBD_phenotype_input")
input_file <- file.path(input_dir, "IBD_nonIBD_formal_cell_assignments.tsv.gz")
pass_file <- file.path(input_dir, "IBD_nonIBD_phenotype_input_complete.json")
output_dir <- file.path(RUN_ROOT, "output/23_IBD_nonIBD_known_phenotype_enrichment_Cell2021", dataset)
if (!file.exists(input_file) || !file.exists(pass_file)) stop("Run script 33 before phenotype enrichment", call. = FALSE)
pass <- jsonlite::read_json(pass_file, simplifyVector = TRUE)
if (pass$status != "PASS_IBD_NONIBD_PHENOTYPE_INPUT_AUDIT") stop("Unified phenotype input audit did not pass", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cells <- fread(input_file, colClasses = "character")
cells <- cells[Dataset == dataset]
axes <- c("BroadPhenotype", "Disease", "Inflammation", "Region3")
required_columns <- c("Dataset", "CellType", "ID", "State", "DonorID", "BiologicalSampleID", axes, "FormalStateUniverse")
if (!all(required_columns %in% names(cells))) stop("Unified phenotype input columns changed", call. = FALSE)
if (anyDuplicated(cells[, .(CellType, ID)])) stop("Duplicate cells in the selected dataset", call. = FALSE)

# Step2. 逐父类型和表型执行原文式标签置换
permutations <- 1000L
base_seed <- 20261124L
z_threshold <- 1.96
job_keys <- unique(cells[, .(CellType, FormalStateUniverse)])
setorder(job_keys, CellType)
result_parts <- list()
denominator_parts <- list()
result_index <- 0L

for (job_index in seq_len(nrow(job_keys))) {
  cell_type <- job_keys$CellType[[job_index]]
  state_universe <- strsplit(job_keys$FormalStateUniverse[[job_index]], ";", fixed = TRUE)[[1L]]
  parent_cells <- cells[CellType == cell_type]

  for (axis_index in seq_along(axes)) {
    axis_name <- axes[[axis_index]]
    seed <- base_seed + match(dataset, datasets) * 10000L + job_index * 100L + axis_index
    observations <- parent_cells[State != "Unassigned" & !is.na(get(axis_name)) & nzchar(get(axis_name))]
    axis_levels <- sort(unique(observations[[axis_name]]))

    denominator_parts[[length(denominator_parts) + 1L]] <- parent_cells[!is.na(get(axis_name)) & nzchar(get(axis_name)), .(
      InputCells = .N, AssignedCells = sum(State != "Unassigned"), BiologicalSamples = uniqueN(BiologicalSampleID), Donors = uniqueN(DonorID)
    ), by = .(Level = get(axis_name))][, `:=`(Dataset = dataset, CellType = cell_type, Axis = axis_name)]

    if (!nrow(observations) || length(axis_levels) < 2L) {
      grid <- CJ(Level = axis_levels, State = state_universe, sorted = FALSE)
      if (!nrow(grid)) grid <- data.table(Level = "", State = state_universe)
      grid[, `:=`(
        ObservedAssignedCells = 0L, NullMean = NA_real_, NullSD = NA_real_, NullQ025 = NA_real_, NullQ975 = NA_real_,
        Z = NA_real_, PositiveZgt1.96 = FALSE, P_one_sided = NA_real_, EmpiricalP_one_sided = NA_real_,
        Status = "NOT_EVALUABLE_LT2_LEVELS"
      )]
    } else {
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
        Level = rep(axis_levels, each = length(state_universe)), State = rep(state_universe, times = length(axis_levels)),
        ObservedAssignedCells = as.integer(as.vector(observed_table)), NullMean = colMeans(null_counts), NullSD = apply(null_counts, 2L, sd),
        NullQ025 = apply(null_counts, 2L, quantile, probs = 0.025, names = FALSE), NullQ975 = apply(null_counts, 2L, quantile, probs = 0.975, names = FALSE)
      )
      grid[, Z := fifelse(NullSD > 0, (ObservedAssignedCells - NullMean) / NullSD, NA_real_)]
      grid[, PositiveZgt1.96 := is.finite(Z) & Z > z_threshold]
      grid[, P_one_sided := fifelse(is.finite(Z), pnorm(Z, lower.tail = FALSE), NA_real_)]
      observed_matrix <- matrix(grid$ObservedAssignedCells, nrow = permutations, ncol = nrow(grid), byrow = TRUE)
      grid[, EmpiricalP_one_sided := (1 + colSums(null_counts >= observed_matrix)) / (permutations + 1)]
      grid[, Status := fifelse(is.finite(Z), "EVALUABLE", "NOT_EVALUABLE_ZERO_NULL_SD")]
    }

    grid[, `:=`(
      Dataset = dataset, CellType = cell_type, CellState = paste(cell_type, State, sep = "_"), Axis = axis_name,
      Permutations = permutations, Seed = seed, Method = "CELL2021_STATE_LABEL_PERMUTATION_IBD_NONIBD_EXTENSION"
    )]
    result_index <- result_index + 1L
    result_parts[[result_index]] <- grid
  }
  cat("CELL2021_PHENOTYPE_COMPLETE", dataset, cell_type, "\n")
}

# Step3. 保存每个数据集的正式结果
results <- rbindlist(result_parts, use.names = TRUE, fill = TRUE)
results[, FDR_BH_within_dataset_axis := p.adjust(P_one_sided, method = "BH"), by = Axis]
results[, Significant_BH_sensitivity := is.finite(Z) & Z > 0 & FDR_BH_within_dataset_axis < 0.05]
setcolorder(results, c(
  "Dataset", "CellType", "State", "CellState", "Axis", "Level", "ObservedAssignedCells", "NullMean", "NullSD", "NullQ025", "NullQ975",
  "Z", "PositiveZgt1.96", "P_one_sided", "EmpiricalP_one_sided", "FDR_BH_within_dataset_axis", "Significant_BH_sensitivity",
  "Permutations", "Seed", "Method", "Status"
))
setorder(results, CellType, Axis, Level, State)
denominators <- rbindlist(denominator_parts, use.names = TRUE, fill = TRUE)
setcolorder(denominators, c("Dataset", "CellType", "Axis", "Level", "InputCells", "AssignedCells", "BiologicalSamples", "Donors"))
setorder(denominators, CellType, Axis, Level)

fwrite(results, file.path(output_dir, "state_known_phenotype_enrichment_Cell2021.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(results[PositiveZgt1.96 == TRUE], file.path(output_dir, "state_known_phenotype_enrichment_Cell2021_positive.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(denominators, file.path(output_dir, "state_known_phenotype_enrichment_denominators.tsv"), sep = "\t", quote = FALSE, na = "")

completion <- list(
  status = "PASS_IBD_NONIBD_KNOWN_PHENOTYPE_ENRICHMENT_CELL2021", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  dataset = dataset, parent_types = nrow(job_keys), input_cells = nrow(cells), assigned_cells = sum(cells$State != "Unassigned"),
  permutations = permutations, axes = axes, evaluable_rows = sum(results$Status == "EVALUABLE"),
  positive_z_rows = sum(results$PositiveZgt1.96), method = "CELL2021_STATE_LABEL_PERMUTATION_IBD_NONIBD_EXTENSION"
)
jsonlite::write_json(completion, file.path(output_dir, "known_phenotype_enrichment_Cell2021_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_IBD_NONIBD_KNOWN_PHENOTYPE_ENRICHMENT_CELL2021", dataset, "\n")
