#!/usr/bin/env Rscript
# 本脚本参照 Cell 2021 Figure S2G-S2H：给每个样本确定父类型内的 dominant state，再做 Fisher representation 检验。
# 20 个父类型细胞是主分析门槛，10 个是敏感性分析；结果写入 output/24_IBD_nonIBD_state_representation。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取统一逐细胞分配
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 35_run_IBD_nonIBD_state_representation_fisher.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

input_dir <- file.path(RUN_ROOT, "output/22_IBD_nonIBD_phenotype_input")
input_file <- file.path(input_dir, "IBD_nonIBD_formal_cell_assignments.tsv.gz")
pass_file <- file.path(input_dir, "IBD_nonIBD_phenotype_input_complete.json")
output_dir <- file.path(RUN_ROOT, "output/24_IBD_nonIBD_state_representation")
if (!file.exists(input_file) || !file.exists(pass_file)) stop("Run script 33 before state representation", call. = FALSE)
pass <- jsonlite::read_json(pass_file, simplifyVector = TRUE)
if (pass$status != "PASS_IBD_NONIBD_PHENOTYPE_INPUT_AUDIT") stop("Unified phenotype input audit did not pass", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cells <- fread(input_file, colClasses = "character")
axes <- c("BroadPhenotype", "Disease", "Inflammation", "Region3")
required_columns <- c("AnalysisGroup", "Dataset", "CellType", "ID", "State", "DonorID", "BiologicalSampleID", axes, "FormalStateUniverse")
if (!all(required_columns %in% names(cells))) stop("Unified phenotype input columns changed", call. = FALSE)
if (anyDuplicated(cells[, .(Dataset, CellType, ID)])) stop("Duplicate cells in unified phenotype input", call. = FALSE)

# Step2. 冻结每个样本的 dominant state
sample_meta <- unique(cells[, c("AnalysisGroup", "Dataset", "CellType", "BiologicalSampleID", "DonorID", axes), with = FALSE])
if (sample_meta[, anyDuplicated(paste(AnalysisGroup, Dataset, CellType, BiologicalSampleID))] != 0L) stop("A biological sample has conflicting phenotype labels", call. = FALSE)

sample_totals <- cells[, .(
  InputCells = .N, AssignedCells = sum(State != "Unassigned"), UnassignedCells = sum(State == "Unassigned")
), by = .(AnalysisGroup, Dataset, CellType, BiologicalSampleID)]
assigned_counts <- cells[State != "Unassigned", .N, by = .(AnalysisGroup, Dataset, CellType, BiologicalSampleID, State)]
max_counts <- assigned_counts[, .(MaxN = max(N)), by = .(AnalysisGroup, Dataset, CellType, BiologicalSampleID)]
dominant_candidates <- merge(assigned_counts, max_counts, by = c("AnalysisGroup", "Dataset", "CellType", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)[N == MaxN]
dominant <- dominant_candidates[, .(
  DominantState = if (.N == 1L) State else "TiedDominant", DominantAssignedCells = MaxN[[1L]], DominantTies = .N
), by = .(AnalysisGroup, Dataset, CellType, BiologicalSampleID)]

sample_assignments <- merge(sample_meta, sample_totals, by = c("AnalysisGroup", "Dataset", "CellType", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
sample_assignments <- merge(sample_assignments, dominant, by = c("AnalysisGroup", "Dataset", "CellType", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
sample_assignments[is.na(DominantState), `:=`(DominantState = "Unassigned", DominantAssignedCells = 0L, DominantTies = 0L)]
sample_assignments[, AssignmentRate := AssignedCells / InputCells]
setorder(sample_assignments, AnalysisGroup, Dataset, CellType, BiologicalSampleID)
fwrite(sample_assignments, file.path(output_dir, "sample_dominant_state_assignments_all.tsv"), sep = "\t", quote = FALSE, na = "")

# Step3. 按原文结构做 sample-level Fisher 检验
thresholds <- c(20L, 10L)
fisher_parts <- list()
fisher_index <- 0L
datasets <- sort(unique(sample_assignments$Dataset))
analysis_specs <- rbindlist(list(
  data.table(Scope = "POOLED_FOUR_DATASETS", DatasetScope = "ALL", Axis = c(axes, "Dataset")),
  CJ(Scope = "WITHIN_DATASET", DatasetScope = datasets, Axis = axes, unique = TRUE)
))

for (threshold in thresholds) {
  eligible_samples <- sample_assignments[InputCells >= threshold]

  for (spec_index in seq_len(nrow(analysis_specs))) {
    scope_name <- analysis_specs$Scope[[spec_index]]
    dataset_scope <- analysis_specs$DatasetScope[[spec_index]]
    axis_name <- analysis_specs$Axis[[spec_index]]
    scope_samples <- if (scope_name == "WITHIN_DATASET") eligible_samples[Dataset == dataset_scope] else copy(eligible_samples)
    scope_samples[, Level := get(axis_name)]
    scope_samples <- scope_samples[!is.na(Level) & nzchar(Level)]

    for (cell_type in sort(unique(scope_samples$CellType))) {
      parent <- scope_samples[CellType == cell_type]
      formal_universes <- unique(cells[CellType == cell_type, FormalStateUniverse])
      if (length(formal_universes) != 1L) stop("Formal state universe differs across datasets: ", cell_type, call. = FALSE)
      state_universe <- strsplit(formal_universes[[1L]], ";", fixed = TRUE)[[1L]]
      levels_axis <- sort(unique(parent$Level))
      testable <- parent[!DominantState %in% c("Unassigned", "TiedDominant")]

      for (level_name in levels_axis) {
        for (state_name in state_universe) {
          a <- sum(testable$Level == level_name & testable$DominantState == state_name)
          b <- sum(testable$Level == level_name & testable$DominantState != state_name)
          c_value <- sum(testable$Level != level_name & testable$DominantState == state_name)
          d <- sum(testable$Level != level_name & testable$DominantState != state_name)
          contingency <- matrix(c(a, b, c_value, d), nrow = 2L, byrow = TRUE)
          evaluable <- length(levels_axis) >= 2L && sum(contingency[1L, ]) > 0L && sum(contingency[2L, ]) > 0L && sum(contingency[, 1L]) > 0L && sum(contingency[, 2L]) > 0L
          fisher_result <- if (evaluable) fisher.test(contingency, alternative = "greater") else NULL

          fisher_index <- fisher_index + 1L
          fisher_parts[[fisher_index]] <- data.table(
            MinParentCells = threshold, Scope = scope_name, DatasetScope = dataset_scope, Axis = axis_name, Level = level_name,
            CellType = cell_type, State = state_name, CellState = paste(cell_type, state_name, sep = "_"),
            SamplesTargetLevelState = a, SamplesTargetLevelOtherStates = b, SamplesOtherLevelsState = c_value, SamplesOtherLevelsOtherStates = d,
            OddsRatio = if (evaluable) unname(fisher_result$estimate) else NA_real_, P_one_sided = if (evaluable) fisher_result$p.value else NA_real_,
            EligibleSampleParents = nrow(parent), AssignedDominantSamples = nrow(testable), UnassignedSamples = sum(parent$DominantState == "Unassigned"),
            TiedDominantSamples = sum(parent$DominantState == "TiedDominant"), Donors = uniqueN(parent$DonorID),
            Status = if (evaluable) "EVALUABLE" else "NOT_EVALUABLE_FISHER_MARGIN"
          )
        }
      }
    }
  }
  cat("REPRESENTATION_THRESHOLD_COMPLETE", threshold, "\n")
}

# Step4. 保存原文式 representation 结果
fisher_results <- rbindlist(fisher_parts, use.names = TRUE, fill = TRUE)
fisher_results[, FDR_BH_across_cell_types := p.adjust(P_one_sided, method = "BH"), by = .(MinParentCells, Scope, DatasetScope, Axis, Level)]
fisher_results[, SignificantEnrichment := is.finite(OddsRatio) & OddsRatio > 1 & FDR_BH_across_cell_types < 0.05]
setorder(fisher_results, MinParentCells, Scope, DatasetScope, Axis, Level, CellType, State)

fwrite(fisher_results, file.path(output_dir, "sample_dominant_state_representation_fisher.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(fisher_results[SignificantEnrichment == TRUE], file.path(output_dir, "sample_dominant_state_representation_significant.tsv"), sep = "\t", quote = FALSE, na = "")

completion <- list(
  status = "PASS_IBD_NONIBD_CELL_STATE_REPRESENTATION", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  method = "CELL2021_S2G_S2H_DOMINANT_STATE_FISHER_IBD_NONIBD_EXTENSION", datasets = uniqueN(cells$Dataset),
  input_cells = nrow(cells), biological_samples = uniqueN(paste(cells$Dataset, cells$BiologicalSampleID)),
  thresholds = thresholds, evaluable_fisher_rows = sum(fisher_results$Status == "EVALUABLE"),
  significant_rows = sum(fisher_results$SignificantEnrichment)
)
jsonlite::write_json(completion, file.path(output_dir, "cell_state_representation_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_IBD_NONIBD_CELL_STATE_REPRESENTATION\n")
