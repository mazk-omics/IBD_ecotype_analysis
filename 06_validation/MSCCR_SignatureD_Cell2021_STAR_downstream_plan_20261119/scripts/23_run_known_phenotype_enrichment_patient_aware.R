#!/usr/bin/env Rscript
# 本脚本在 G7 原文逐细胞轨道之外，补充患者感知的 IBD 表型检验。
# Disease 在供者层打乱标签；Inflammation 和 Region3 仅在同一供者内打乱样本标签。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读入正式逐细胞分配
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 23_run_known_phenotype_enrichment_patient_aware.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

permutations <- 1000L
base_seed <- 20261123L
z_threshold <- 1.96
axes <- c("Disease", "Inflammation", "Region3")
input_file <- file.path(RUN_ROOT, "output/11_known_phenotype_enrichment_exact/formal_cell_assignments_with_phenotypes.tsv.gz")
output_dir <- file.path(RUN_ROOT, "output/12_known_phenotype_enrichment_patient_aware")
if (!file.exists(input_file)) stop("Run script 22 before the patient-aware extension", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cells <- fread(input_file, colClasses = "character")
required_columns <- c("Dataset", "CellType", "ID", "State", "DonorID", "BiologicalSampleID", axes, "FormalStateUniverse")
if (!all(required_columns %in% names(cells))) stop("Formal cell assignment columns changed", call. = FALSE)
if (anyDuplicated(cells[, .(Dataset, CellType, ID)])) stop("Duplicate cells in formal assignment input", call. = FALSE)

# Step2. 逐数据集、父类型和表型建立样本比例
result_parts <- list()
result_index <- 0L
denominator_parts <- list()
job_keys <- unique(cells[, .(Dataset, CellType, FormalStateUniverse)])
setorder(job_keys, Dataset, CellType)

for (job_index in seq_len(nrow(job_keys))) {
  dataset_name <- job_keys$Dataset[[job_index]]
  cell_type <- job_keys$CellType[[job_index]]
  state_universe <- strsplit(job_keys$FormalStateUniverse[[job_index]], ";", fixed = TRUE)[[1L]]
  parent <- cells[Dataset == dataset_name & CellType == cell_type & State != "Unassigned"]

  for (axis_index in seq_along(axes)) {
    axis_name <- axes[[axis_index]]
    seed <- base_seed + job_index * 100L + axis_index
    observations <- parent[!is.na(get(axis_name)) & nzchar(get(axis_name)) & !is.na(DonorID) & nzchar(DonorID) & !is.na(BiologicalSampleID) & nzchar(BiologicalSampleID)]
    sample_meta <- unique(observations[, .(BiologicalSampleID, DonorID, Level = get(axis_name))])
    if (sample_meta[, anyDuplicated(BiologicalSampleID)] != 0L) stop("One sample has multiple donor or phenotype labels: ", dataset_name, " ", cell_type, " ", axis_name, call. = FALSE)

    sample_counts <- observations[, .N, by = .(BiologicalSampleID, State)]
    sample_totals <- sample_counts[, .(AssignedCells = sum(N)), by = BiologicalSampleID]
    sample_grid <- CJ(BiologicalSampleID = sample_meta$BiologicalSampleID, State = state_universe, unique = TRUE)
    sample_grid <- merge(sample_grid, sample_counts, by = c("BiologicalSampleID", "State"), all.x = TRUE, sort = FALSE)
    sample_grid[is.na(N), N := 0L]
    sample_grid <- merge(sample_grid, sample_totals, by = "BiologicalSampleID", all.x = TRUE, sort = FALSE)
    sample_grid <- merge(sample_grid, sample_meta, by = "BiologicalSampleID", all.x = TRUE, sort = FALSE)
    sample_grid[, StateFraction := N / AssignedCells]

    if (axis_name == "Disease") {
      donor_labels <- unique(sample_meta[, .(DonorID, Level)])
      ambiguous_donors <- donor_labels[, .(Levels = uniqueN(Level)), by = DonorID][Levels != 1L]
      if (nrow(ambiguous_donors)) stop("Disease varies within donor: ", dataset_name, " ", cell_type, call. = FALSE)
      donor_state <- sample_grid[, .(StateFraction = mean(StateFraction), BiologicalSamples = uniqueN(BiologicalSampleID)), by = .(DonorID, State)]
      donor_state <- merge(donor_state, donor_labels, by = "DonorID", all.x = TRUE, sort = FALSE)
      levels_axis <- sort(unique(donor_labels$Level))
      denominator_parts[[length(denominator_parts) + 1L]] <- data.table(
        Dataset = dataset_name, CellType = cell_type, Axis = axis_name, EligibleDonors = uniqueN(donor_labels$DonorID),
        EligibleSamples = uniqueN(sample_meta$BiologicalSampleID), Levels = paste(levels_axis, collapse = ";"),
        Restriction = "Disease labels permuted across donors"
      )

      if (length(levels_axis) < 2L || uniqueN(donor_labels$DonorID) < 3L) {
        grid <- CJ(Level = levels_axis, State = state_universe, sorted = FALSE)
        if (!nrow(grid)) grid <- data.table(Level = "", State = state_universe)
        grid[, `:=`(ObservedMeanStateFraction = NA_real_, NullMean = NA_real_, NullSD = NA_real_, NullQ025 = NA_real_, NullQ975 = NA_real_, Z = NA_real_, PositiveZgt1.96 = FALSE, EmpiricalP_one_sided = NA_real_, Status = "NOT_EVALUABLE_DONOR_LEVELS")]
      } else {
        observed <- donor_state[, .(ObservedMeanStateFraction = mean(StateFraction)), by = .(Level, State)]
        grid <- CJ(Level = levels_axis, State = state_universe, sorted = FALSE)
        grid <- merge(grid, observed, by = c("Level", "State"), all.x = TRUE, sort = FALSE)
        grid[is.na(ObservedMeanStateFraction), ObservedMeanStateFraction := 0]
        null_values <- matrix(NA_real_, nrow = permutations, ncol = nrow(grid))
        set.seed(seed)
        for (permutation_index in seq_len(permutations)) {
          shuffled_labels <- copy(donor_labels)
          shuffled_labels[, Level := sample(Level, .N, replace = FALSE)]
          permuted <- merge(donor_state[, .(DonorID, State, StateFraction)], shuffled_labels, by = "DonorID", all.x = TRUE, sort = FALSE)
          permuted_mean <- permuted[, .(Value = mean(StateFraction)), by = .(Level, State)]
          null_values[permutation_index, ] <- permuted_mean$Value[match(paste(grid$Level, grid$State), paste(permuted_mean$Level, permuted_mean$State))]
        }
        grid[, `:=`(NullMean = colMeans(null_values), NullSD = apply(null_values, 2L, sd), NullQ025 = apply(null_values, 2L, quantile, probs = 0.025, names = FALSE), NullQ975 = apply(null_values, 2L, quantile, probs = 0.975, names = FALSE))]
        grid[, Z := fifelse(NullSD > 0, (ObservedMeanStateFraction - NullMean) / NullSD, NA_real_)]
        observed_matrix <- matrix(grid$ObservedMeanStateFraction, nrow = permutations, ncol = nrow(grid), byrow = TRUE)
        grid[, EmpiricalP_one_sided := (1 + colSums(null_values >= observed_matrix)) / (permutations + 1)]
        grid[, `:=`(PositiveZgt1.96 = is.finite(Z) & Z > z_threshold, Status = fifelse(is.finite(Z), "EVALUABLE", "NOT_EVALUABLE_ZERO_NULL_SD"))]
      }
    } else {
      eligible_donors <- sample_meta[, .(Levels = uniqueN(Level)), by = DonorID][Levels >= 2L]$DonorID
      eligible_meta <- sample_meta[DonorID %in% eligible_donors]
      eligible_grid <- sample_grid[DonorID %in% eligible_donors]
      levels_axis <- sort(unique(eligible_meta$Level))
      denominator_parts[[length(denominator_parts) + 1L]] <- data.table(
        Dataset = dataset_name, CellType = cell_type, Axis = axis_name, EligibleDonors = length(eligible_donors),
        EligibleSamples = uniqueN(eligible_meta$BiologicalSampleID), Levels = paste(levels_axis, collapse = ";"),
        Restriction = "Sample phenotype labels permuted within donor"
      )

      if (length(levels_axis) < 2L || length(eligible_donors) < 1L) {
        grid <- CJ(Level = levels_axis, State = state_universe, sorted = FALSE)
        if (!nrow(grid)) grid <- data.table(Level = "", State = state_universe)
        grid[, `:=`(ObservedMeanStateFraction = NA_real_, NullMean = NA_real_, NullSD = NA_real_, NullQ025 = NA_real_, NullQ975 = NA_real_, Z = NA_real_, PositiveZgt1.96 = FALSE, EmpiricalP_one_sided = NA_real_, Status = "NOT_EVALUABLE_WITHIN_DONOR_LEVELS")]
      } else {
        observed_donor <- eligible_grid[, .(DonorMeanStateFraction = mean(StateFraction)), by = .(DonorID, Level, State)]
        observed <- observed_donor[, .(ObservedMeanStateFraction = mean(DonorMeanStateFraction)), by = .(Level, State)]
        grid <- CJ(Level = levels_axis, State = state_universe, sorted = FALSE)
        grid <- merge(grid, observed, by = c("Level", "State"), all.x = TRUE, sort = FALSE)
        null_values <- matrix(NA_real_, nrow = permutations, ncol = nrow(grid))
        set.seed(seed)
        for (permutation_index in seq_len(permutations)) {
          shuffled_meta <- copy(eligible_meta)
          shuffled_meta[, Level := sample(Level, .N, replace = FALSE), by = DonorID]
          permuted <- merge(eligible_grid[, .(BiologicalSampleID, DonorID, State, StateFraction)], shuffled_meta[, .(BiologicalSampleID, Level)], by = "BiologicalSampleID", all.x = TRUE, sort = FALSE)
          permuted_donor <- permuted[, .(DonorMeanStateFraction = mean(StateFraction)), by = .(DonorID, Level, State)]
          permuted_mean <- permuted_donor[, .(Value = mean(DonorMeanStateFraction)), by = .(Level, State)]
          null_values[permutation_index, ] <- permuted_mean$Value[match(paste(grid$Level, grid$State), paste(permuted_mean$Level, permuted_mean$State))]
        }
        grid[, `:=`(NullMean = colMeans(null_values, na.rm = TRUE), NullSD = apply(null_values, 2L, sd, na.rm = TRUE), NullQ025 = apply(null_values, 2L, quantile, probs = 0.025, names = FALSE, na.rm = TRUE), NullQ975 = apply(null_values, 2L, quantile, probs = 0.975, names = FALSE, na.rm = TRUE))]
        grid[, Z := fifelse(NullSD > 0, (ObservedMeanStateFraction - NullMean) / NullSD, NA_real_)]
        observed_matrix <- matrix(grid$ObservedMeanStateFraction, nrow = permutations, ncol = nrow(grid), byrow = TRUE)
        grid[, EmpiricalP_one_sided := (1 + colSums(null_values >= observed_matrix, na.rm = TRUE)) / (permutations + 1)]
        grid[, `:=`(PositiveZgt1.96 = is.finite(Z) & Z > z_threshold, Status = fifelse(is.finite(Z), "EVALUABLE", "NOT_EVALUABLE_ZERO_NULL_SD"))]
      }
    }

    grid[, `:=`(
      Dataset = dataset_name, CellType = cell_type, Axis = axis_name, Permutations = permutations, Seed = seed,
      Method = ifelse(axis_name == "Disease", "DONOR_LABEL_PERMUTATION", "WITHIN_DONOR_SAMPLE_LABEL_PERMUTATION")
    )]
    result_index <- result_index + 1L
    result_parts[[result_index]] <- grid
  }
  cat("PATIENT_AWARE_COMPLETE", dataset_name, cell_type, "\n")
}

# Step3. 保存患者感知扩展结果
results <- rbindlist(result_parts, use.names = TRUE, fill = TRUE)
results[, CellState := paste(CellType, State, sep = "_")]
results[, P_one_sided := fifelse(is.finite(Z), pnorm(Z, lower.tail = FALSE), NA_real_)]
results[, FDR_BH_within_dataset_axis := p.adjust(P_one_sided, method = "BH"), by = .(Dataset, Axis)]
results[, Significant_BH_sensitivity := is.finite(Z) & Z > 0 & FDR_BH_within_dataset_axis < 0.05]
setcolorder(results, c("Dataset", "CellType", "State", "CellState", "Axis", "Level", "ObservedMeanStateFraction", "NullMean", "NullSD", "NullQ025", "NullQ975", "Z", "PositiveZgt1.96", "P_one_sided", "EmpiricalP_one_sided", "FDR_BH_within_dataset_axis", "Significant_BH_sensitivity", "Permutations", "Seed", "Method", "Status"))
setorder(results, Dataset, CellType, Axis, Level, State)
denominators <- rbindlist(denominator_parts, use.names = TRUE, fill = TRUE)
setorder(denominators, Dataset, CellType, Axis)

fwrite(results, file.path(output_dir, "state_known_phenotype_enrichment_patient_aware.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(results[PositiveZgt1.96 == TRUE], file.path(output_dir, "state_known_phenotype_enrichment_patient_aware_positive.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(denominators, file.path(output_dir, "patient_aware_eligible_units.tsv"), sep = "\t", quote = FALSE, na = "")

completion <- list(
  status = "PASS_KNOWN_PHENOTYPE_ENRICHMENT_PATIENT_AWARE", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = uniqueN(results$Dataset), parent_jobs = nrow(job_keys), permutations = permutations,
  evaluable_rows = sum(results$Status == "EVALUABLE"), positive_z_rows = sum(results$PositiveZgt1.96),
  method = "SEPARATE_IBD_PATIENT_AWARE_EXTENSION"
)
jsonlite::write_json(completion, file.path(output_dir, "known_phenotype_enrichment_patient_aware_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_KNOWN_PHENOTYPE_ENRICHMENT_PATIENT_AWARE\n")
