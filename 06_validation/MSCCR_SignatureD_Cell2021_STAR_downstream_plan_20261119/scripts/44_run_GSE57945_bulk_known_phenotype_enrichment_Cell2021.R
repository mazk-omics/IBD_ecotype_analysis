#!/usr/bin/env Rscript
# 本脚本对 GSE57945 bulk recovery 的正式 state assignment 执行 Cell 2021 known phenotype enrichment。
# 在每个父细胞类型内打乱 state labels 1,000 次；Non-IBD 与 IBD 的结果写入 output/32。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取已验收的 bulk recovery
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 44_run_GSE57945_bulk_known_phenotype_enrichment_Cell2021.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

project_root <- normalizePath(file.path(RUN_ROOT, "../.."), mustWork = TRUE)
ecotyper_root <- file.path(project_root, "software/ecotyper")
discovery_name <- "MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101"
model_root <- file.path(ecotyper_root, "EcoTyper", discovery_name, "Custom")
job_root <- file.path(RUN_ROOT, "output/31_bulk_recovery/GSE57945_current85_20261124")
result_root <- file.path(job_root, "RecoveryOutput/GSE57945_EcoTyper_recovery_expression")
annotation_file <- file.path(job_root, "input/GSE57945_current85_annotation_v2.tsv")
status_file <- file.path(job_root, "status.txt")
output_dir <- file.path(RUN_ROOT, "output/32_GSE57945_bulk_known_phenotype_Cell2021")

required_files <- c(status_file, annotation_file, file.path(model_root, "Analysis/rank_selection/rank_data.txt"))
if (any(!file.exists(required_files)) || !dir.exists(result_root)) stop("Accepted GSE57945 recovery inputs are incomplete", call. = FALSE)
if (!any(readLines(status_file, warn = FALSE) == "status=ACCEPTED")) stop("GSE57945 recovery is not ACCEPTED", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation <- fread(annotation_file, colClasses = "character")
rank_data <- fread(file.path(model_root, "Analysis/rank_selection/rank_data.txt"))
if (!all(c("ID", "Disease") %in% names(annotation))) stop("GSE57945 annotation columns changed", call. = FALSE)
if (!all(c("CellType", "Chosen_Rank") %in% names(rank_data))) stop("EcoTyper rank data columns changed", call. = FALSE)
if (anyDuplicated(annotation$ID) || nrow(annotation) != 322L) stop("GSE57945 annotation identity changed", call. = FALSE)
annotation[, BroadPhenotype := fifelse(Disease == "Non-IBD", "nonIBD", "IBD")]
print(annotation[, .N, by = .(Disease, BroadPhenotype)])

# Step2. 冻结 20 个父类型的 assignment、state universe 和分母
cell_parts <- vector("list", nrow(rank_data))
state_universe_parts <- vector("list", nrow(rank_data))
denominator_parts <- vector("list", nrow(rank_data))

for (cell_index in seq_len(nrow(rank_data))) {
  cell_type <- rank_data$CellType[[cell_index]]
  chosen_rank <- rank_data$Chosen_Rank[[cell_index]]
  assignment_file <- file.path(result_root, cell_type, "state_assignment.txt")
  mapping_file <- file.path(model_root, "Cell_States/discovery", cell_type, as.character(chosen_rank), "mapping_to_initial_states.txt")
  if (!file.exists(assignment_file) || !file.exists(mapping_file)) stop("Missing assignment or state mapping for ", cell_type, call. = FALSE)

  assignment <- fread(assignment_file, colClasses = "character")
  mapping <- fread(mapping_file, colClasses = "character")
  if (!all(c("ID", "State") %in% names(assignment)) || !"State" %in% names(mapping)) stop("Recovery output columns changed", call. = FALSE)
  if (anyDuplicated(assignment$ID) || !all(assignment$ID %in% annotation$ID)) stop("Invalid assignment IDs for ", cell_type, call. = FALSE)

  cells <- merge(annotation[, .(ID, Disease, BroadPhenotype)], assignment[, .(ID, State)], by = "ID", all.x = TRUE, sort = FALSE)
  cells <- cells[match(annotation$ID, ID)]
  cells[is.na(State) | !nzchar(State), State := "Unassigned"]
  cells[, CellType := cell_type]
  cell_parts[[cell_index]] <- cells

  state_universe <- sort(unique(mapping$State))
  state_universe_parts[[cell_index]] <- data.table(CellType = cell_type, State = state_universe, CellState = paste(cell_type, state_universe, sep = "_"))
  denominator_parts[[cell_index]] <- cells[, .(
    InputSamples = .N, AssignedSamples = sum(State != "Unassigned"), UnassignedSamples = sum(State == "Unassigned")
  ), by = .(CellType, BroadPhenotype)]
  cat("BULK_INPUT_READY", cell_type, "states", length(state_universe), "assigned", sum(cells$State != "Unassigned"), "\n")
}

cells_all <- rbindlist(cell_parts, use.names = TRUE)
state_universe_all <- rbindlist(state_universe_parts, use.names = TRUE)
denominators <- rbindlist(denominator_parts, use.names = TRUE)
if (uniqueN(cells_all$CellType) != 20L || uniqueN(state_universe_all$CellState) != 85L) stop("Expected 20 cell types and 85 formal states", call. = FALSE)
setorder(denominators, CellType, BroadPhenotype)
fwrite(state_universe_all, file.path(output_dir, "formal_state_universe.tsv"), sep = "\t", quote = FALSE)
fwrite(denominators, file.path(output_dir, "bulk_known_phenotype_denominators.tsv"), sep = "\t", quote = FALSE)

# Step3. 在父类型内执行 1,000 次 state-label shuffle
permutations <- 1000L
base_seed <- 20261124L
z_threshold <- 1.96
phenotype_levels <- c("IBD", "nonIBD")
result_parts <- vector("list", nrow(rank_data))

for (cell_index in seq_len(nrow(rank_data))) {
  cell_type <- rank_data$CellType[[cell_index]]
  state_universe <- state_universe_all[CellType == cell_type, State]
  observations <- cells_all[CellType == cell_type & State != "Unassigned"]
  if (!all(phenotype_levels %in% unique(observations$BroadPhenotype))) stop("Both phenotype levels are required for ", cell_type, call. = FALSE)

  labels <- observations$State
  phenotype <- observations$BroadPhenotype
  observed_table <- table(factor(labels, levels = state_universe), factor(phenotype, levels = phenotype_levels))
  null_counts <- matrix(0L, nrow = permutations, ncol = length(state_universe) * length(phenotype_levels))
  seed <- base_seed + cell_index
  set.seed(seed)

  for (permutation_index in seq_len(permutations)) {
    shuffled <- sample(labels, length(labels), replace = FALSE)
    permuted_table <- table(factor(shuffled, levels = state_universe), factor(phenotype, levels = phenotype_levels))
    null_counts[permutation_index, ] <- as.vector(permuted_table)
  }

  results <- data.table(
    Level = rep(phenotype_levels, each = length(state_universe)), State = rep(state_universe, times = length(phenotype_levels)),
    ObservedAssignedSamples = as.integer(as.vector(observed_table)), NullMean = colMeans(null_counts), NullSD = apply(null_counts, 2L, sd),
    NullQ025 = apply(null_counts, 2L, quantile, probs = 0.025, names = FALSE), NullQ975 = apply(null_counts, 2L, quantile, probs = 0.975, names = FALSE)
  )
  results[, Z := fifelse(NullSD > 0, (ObservedAssignedSamples - NullMean) / NullSD, NA_real_)]
  results[, PositiveZgt1.96 := is.finite(Z) & Z > z_threshold]
  results[, P_one_sided := fifelse(is.finite(Z), pnorm(Z, lower.tail = FALSE), NA_real_)]
  observed_matrix <- matrix(results$ObservedAssignedSamples, nrow = permutations, ncol = nrow(results), byrow = TRUE)
  results[, EmpiricalP_one_sided := (1 + colSums(null_counts >= observed_matrix)) / (permutations + 1)]
  results[, `:=`(
    Dataset = "GSE57945", CellType = cell_type, CellState = paste(cell_type, State, sep = "_"), Axis = "BroadPhenotype",
    Permutations = permutations, Seed = seed, Method = "CELL2021_STATE_LABEL_PERMUTATION_BULK",
    Status = fifelse(is.finite(Z), "EVALUABLE", "NOT_EVALUABLE_ZERO_NULL_SD")
  )]
  result_parts[[cell_index]] <- results
  cat("BULK_PHENOTYPE_COMPLETE", cell_type, "\n")
}

# Step4. 保存正式 bulk phenotype enrichment
results <- rbindlist(result_parts, use.names = TRUE)
results[, FDR_BH_within_axis := p.adjust(P_one_sided, method = "BH")]
results[, Significant_BH_sensitivity := is.finite(Z) & Z > 0 & FDR_BH_within_axis < 0.05]
setcolorder(results, c(
  "Dataset", "CellType", "State", "CellState", "Axis", "Level", "ObservedAssignedSamples", "NullMean", "NullSD", "NullQ025", "NullQ975",
  "Z", "PositiveZgt1.96", "P_one_sided", "EmpiricalP_one_sided", "FDR_BH_within_axis", "Significant_BH_sensitivity",
  "Permutations", "Seed", "Method", "Status"
))
setorder(results, CellType, Level, State)
fwrite(results, file.path(output_dir, "GSE57945_bulk_known_phenotype_enrichment_Cell2021.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(results[PositiveZgt1.96 == TRUE], file.path(output_dir, "GSE57945_bulk_known_phenotype_enrichment_Cell2021_positive.tsv"), sep = "\t", quote = FALSE, na = "")

completion <- list(
  status = "PASS_GSE57945_BULK_KNOWN_PHENOTYPE_ENRICHMENT_CELL2021", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  dataset = "GSE57945", samples = nrow(annotation), nonIBD_samples = sum(annotation$BroadPhenotype == "nonIBD"),
  IBD_samples = sum(annotation$BroadPhenotype == "IBD"), parent_types = uniqueN(results$CellType), formal_states = uniqueN(results$CellState),
  permutations = permutations, evaluable_nonIBD_states = results[Level == "nonIBD" & Status == "EVALUABLE", uniqueN(CellState)],
  significant_nonIBD_states = results[Level == "nonIBD" & PositiveZgt1.96 == TRUE, uniqueN(CellState)]
)
jsonlite::write_json(completion, file.path(output_dir, "GSE57945_bulk_known_phenotype_Cell2021_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_GSE57945_BULK_KNOWN_PHENOTYPE_ENRICHMENT_CELL2021\n")
