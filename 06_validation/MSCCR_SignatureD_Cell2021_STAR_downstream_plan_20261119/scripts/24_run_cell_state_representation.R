#!/usr/bin/env Rscript
# 本脚本参照 Cell 2021 Figure S2G-S2H：先给每个样本确定父类型内的 dominant state，再用 Fisher 检验 representation。
# Disease、Inflammation、Region3 和 dataset 分开输出；当前四队列无 nonIBD，相关层明确记为不可评价。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读入正式逐细胞分配
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 24_run_cell_state_representation.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

axes <- c("Disease", "Inflammation", "Region3")
input_file <- file.path(RUN_ROOT, "output/11_known_phenotype_enrichment_exact/formal_cell_assignments_with_phenotypes.tsv.gz")
output_dir <- file.path(RUN_ROOT, "output/13_cell_state_representation")
if (!file.exists(input_file)) stop("Run script 22 before state representation", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cells <- fread(input_file, colClasses = "character")
required_columns <- c("Dataset", "CellType", "ID", "State", "DonorID", "BiologicalSampleID", axes, "FormalStateUniverse")
if (!all(required_columns %in% names(cells))) stop("Formal cell assignment columns changed", call. = FALSE)

# Step2. 冻结逐细胞 representation 和样本 dominant state
cell_representation <- cells[, .(
  CellCount = .N, BiologicalSamplesWithCells = uniqueN(BiologicalSampleID), DonorsWithCells = uniqueN(DonorID)
), by = .(Dataset, CellType, Disease, Inflammation, Region3, State)]
cell_denominators <- cells[, .(
  InputCells = .N, AssignedCells = sum(State != "Unassigned"), UnassignedCells = sum(State == "Unassigned")
), by = .(Dataset, CellType, Disease, Inflammation, Region3)]
cell_representation <- merge(cell_representation, cell_denominators, by = c("Dataset", "CellType", axes), all.x = TRUE, sort = FALSE)
cell_representation[, FractionOfInputCells := CellCount / InputCells]
cell_representation[, FractionAmongAssignedCells := fifelse(State != "Unassigned" & AssignedCells > 0, CellCount / AssignedCells, NA_real_)]
setorder(cell_representation, Dataset, CellType, Disease, Inflammation, Region3, State)

sample_meta <- unique(cells[, c("Dataset", "CellType", "BiologicalSampleID", "DonorID", axes), with = FALSE])
if (sample_meta[, anyDuplicated(paste(Dataset, CellType, BiologicalSampleID))] != 0L) stop("A biological sample has conflicting phenotype labels", call. = FALSE)
assigned_counts <- cells[State != "Unassigned", .N, by = .(Dataset, CellType, BiologicalSampleID, State)]
max_counts <- assigned_counts[, .(MaxN = max(N)), by = .(Dataset, CellType, BiologicalSampleID)]
dominant_candidates <- merge(assigned_counts, max_counts, by = c("Dataset", "CellType", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)[N == MaxN]
dominant <- dominant_candidates[, .(DominantState = if (.N == 1L) State else "TiedDominant", DominantAssignedCells = MaxN[[1L]], DominantTies = .N), by = .(Dataset, CellType, BiologicalSampleID)]
sample_assignments <- merge(sample_meta, dominant, by = c("Dataset", "CellType", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
sample_assignments[is.na(DominantState), `:=`(DominantState = "Unassigned", DominantAssignedCells = 0L, DominantTies = 0L)]

sample_totals <- cells[, .(InputCells = .N, AssignedCells = sum(State != "Unassigned"), UnassignedCells = sum(State == "Unassigned")), by = .(Dataset, CellType, BiologicalSampleID)]
sample_assignments <- merge(sample_assignments, sample_totals, by = c("Dataset", "CellType", "BiologicalSampleID"), all.x = TRUE, sort = FALSE)
setorder(sample_assignments, Dataset, CellType, BiologicalSampleID)

fwrite(cell_representation, file.path(output_dir, "cell_level_state_representation.tsv.gz"), sep = "\t", quote = FALSE, na = "")
fwrite(sample_assignments, file.path(output_dir, "sample_dominant_state_assignments.tsv"), sep = "\t", quote = FALSE, na = "")

# Step3. 按原文结构进行 dominant-state Fisher 检验
fisher_parts <- list()
fisher_index <- 0L

analysis_specs <- rbindlist(list(
  data.table(Scope = "POOLED_FOUR_DATASETS", DatasetScope = "ALL", Axis = c(axes, "Dataset")),
  CJ(Scope = "WITHIN_DATASET", DatasetScope = sort(unique(sample_assignments$Dataset)), Axis = axes, unique = TRUE)
))

for (spec_index in seq_len(nrow(analysis_specs))) {
  scope_name <- analysis_specs$Scope[[spec_index]]
  dataset_scope <- analysis_specs$DatasetScope[[spec_index]]
  axis_name <- analysis_specs$Axis[[spec_index]]
  scope_samples <- if (scope_name == "WITHIN_DATASET") sample_assignments[Dataset == dataset_scope] else copy(sample_assignments)
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
          Scope = scope_name, DatasetScope = dataset_scope, Axis = axis_name, Level = level_name, CellType = cell_type,
          State = state_name, CellState = paste(cell_type, state_name, sep = "_"), SamplesTargetLevelState = a,
          SamplesTargetLevelOtherStates = b, SamplesOtherLevelsState = c_value, SamplesOtherLevelsOtherStates = d,
          OddsRatio = if (evaluable) unname(fisher_result$estimate) else NA_real_, P_one_sided = if (evaluable) fisher_result$p.value else NA_real_,
          InputSamples = nrow(parent), AssignedDominantSamples = nrow(testable), UnassignedSamples = sum(parent$DominantState == "Unassigned"),
          TiedDominantSamples = sum(parent$DominantState == "TiedDominant"), Donors = uniqueN(parent$DonorID),
          Status = if (evaluable) "EVALUABLE" else "NOT_EVALUABLE_FISHER_MARGIN"
        )
      }
    }
  }
}

fisher_results <- rbindlist(fisher_parts, use.names = TRUE, fill = TRUE)
fisher_results[, FDR_BH_across_cell_types := p.adjust(P_one_sided, method = "BH"), by = .(Scope, DatasetScope, Axis)]
fisher_results[, SignificantEnrichment := is.finite(OddsRatio) & OddsRatio > 1 & FDR_BH_across_cell_types < 0.05]
setorder(fisher_results, Scope, DatasetScope, Axis, Level, CellType, State)

# Step4. 明确记录当前 IBD/nonIBD 不可评价
nonibd_status <- unique(sample_assignments[, .(Dataset, Disease)])
nonibd_status[, NonIBDDetected := grepl("healthy|non[-_ ]?ibd|control|^hc$", Disease, ignore.case = TRUE)]
has_nonibd <- any(nonibd_status$NonIBDDetected, na.rm = TRUE)
ibd_nonibd_gate <- data.table(
  Axis = "IBD_nonIBD", Status = if (has_nonibd) "EVALUABLE" else "NOT_EVALUABLE_NO_NONIBD",
  Reason = if (has_nonibd) "Both IBD and nonIBD detected" else "The four current recovery cohorts contain IBD cells only; defer to G11B",
  Datasets = uniqueN(sample_assignments$Dataset), BiologicalSamples = uniqueN(paste(sample_assignments$Dataset, sample_assignments$BiologicalSampleID)),
  Donors = uniqueN(paste(sample_assignments$Dataset, sample_assignments$DonorID))
)

fwrite(fisher_results, file.path(output_dir, "sample_dominant_state_representation_fisher.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(fisher_results[SignificantEnrichment == TRUE], file.path(output_dir, "sample_dominant_state_representation_significant.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(ibd_nonibd_gate, file.path(output_dir, "IBD_nonIBD_representation_gate.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_CELL_STATE_REPRESENTATION_CURRENT_IBD_SCOPE", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  method = "CELL2021_S2G_S2H_DOMINANT_STATE_FISHER_IBD_EXTENSION", datasets = uniqueN(cells$Dataset),
  input_cells = nrow(cells), biological_samples = uniqueN(paste(sample_assignments$Dataset, sample_assignments$BiologicalSampleID)),
  evaluable_fisher_rows = sum(fisher_results$Status == "EVALUABLE"), significant_rows = sum(fisher_results$SignificantEnrichment),
  ibd_nonibd_status = ibd_nonibd_gate$Status
)
jsonlite::write_json(completion, file.path(output_dir, "cell_state_representation_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_CELL_STATE_REPRESENTATION_CURRENT_IBD_SCOPE\n")
