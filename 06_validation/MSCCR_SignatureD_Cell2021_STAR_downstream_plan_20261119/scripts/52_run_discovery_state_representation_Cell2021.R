#!/usr/bin/env Rscript
# 本脚本按 Cell 2021 Fig. S2G/S2H 方法，在 discovery bulk 样本中计算 cell-state representation。
# 每个父类型直接使用 EcoTyper discovery 的正式 state assignment；样本为统计单位。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定 discovery 模型和样本注释
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 52_run_discovery_state_representation_Cell2021.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

project_root <- normalizePath(file.path(RUN_ROOT, "../.."), mustWork = TRUE)
discovery_name <- "MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101"
model_root <- file.path(project_root, "software/ecotyper/EcoTyper", discovery_name, "Custom")
annotation_file <- file.path(project_root, "05_EcoTyper/input/MSCCR_CORE12/annotation/MSCCR_CORE12_annotation.tsv")
output_dir <- file.path(RUN_ROOT, "output/39_discovery_state_representation_Cell2021")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rank_file <- file.path(model_root, "Analysis/rank_selection/rank_data.txt")
if (!file.exists(rank_file) || !file.exists(annotation_file)) stop("Discovery rank data or annotation is missing", call. = FALSE)
rank_data <- fread(rank_file)
annotation <- fread(annotation_file, colClasses = "character")
required_annotation <- c("ID", "PatientID", "Disease", "Inflammation", "Region3")
if (nrow(rank_data) != 20L || !all(required_annotation %in% names(annotation))) stop("Discovery model or annotation schema changed", call. = FALSE)

# Step2. 读取 20 个父类型的正式 discovery state assignment
assignment_parts <- vector("list", nrow(rank_data))
formal_state_parts <- vector("list", nrow(rank_data))

for (cell_index in seq_len(nrow(rank_data))) {
  cell_type <- rank_data$CellType[[cell_index]]
  chosen_rank <- rank_data$Chosen_Rank[[cell_index]]
  assignment_file <- file.path(
    model_root, "Cell_States/discovery", cell_type, as.character(chosen_rank),
    "state_assignment.txt"
  )
  assignment <- suppressWarnings(fread(assignment_file, colClasses = "character"))
  if ("V1" %in% names(assignment)) {
    row_index <- suppressWarnings(as.integer(assignment$V1))
    if (anyNA(row_index) || anyDuplicated(row_index)) stop("Invalid discovery row index: ", cell_type, call. = FALSE)
    assignment[, V1 := NULL]
  }
  if (!all(c("ID", "State") %in% names(assignment))) stop("Discovery assignment columns changed: ", cell_type, call. = FALSE)
  if (anyDuplicated(assignment$ID)) stop("Duplicated discovery sample assignment: ", cell_type, call. = FALSE)

  assignment <- merge(
    assignment[, .(ID, State)],
    annotation[, .(ID, PatientID, Disease, Inflammation, Region3)],
    by = "ID",
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(assignment$Disease) || any(!assignment$Disease %in% c("CD", "UC"))) stop("Discovery assignment contains missing or non-IBD phenotype: ", cell_type, call. = FALSE)
  assignment[, CellType := cell_type]
  assignment_parts[[cell_index]] <- assignment
  formal_state_parts[[cell_index]] <- data.table(CellType = cell_type, State = sort(unique(assignment$State)))
}

assignments <- rbindlist(assignment_parts, use.names = TRUE)
formal_states <- rbindlist(formal_state_parts)
formal_states[, CellState := paste(CellType, State, sep = "_")]
if (nrow(formal_states) != 85L || anyDuplicated(formal_states$CellState)) stop("Formal discovery state universe is not 85", call. = FALSE)
fwrite(assignments, file.path(output_dir, "discovery_sample_state_assignments.tsv"), sep = "\t", quote = FALSE)
fwrite(formal_states, file.path(output_dir, "formal_state_universe.tsv"), sep = "\t", quote = FALSE)

# Step3. 按原文做 state × phenotype 单侧 Fisher 检验
axes <- c("Disease", "Region3", "Inflammation")
fisher_parts <- list()
fisher_index <- 0L

for (axis_name in axes) {
  for (cell_index in seq_len(nrow(rank_data))) {
    cell_type <- rank_data$CellType[[cell_index]]
    parent <- assignments[CellType == cell_type]
    parent[, Level := get(axis_name)]
    axis_levels <- sort(unique(parent$Level))
    state_levels <- formal_states[CellType == cell_type, State]

    for (level_name in axis_levels) {
      for (state_name in state_levels) {
        target_state <- parent$State == state_name
        target_level <- parent$Level == level_name
        a <- sum(target_level & target_state)
        b <- sum(target_level & !target_state)
        c_value <- sum(!target_level & target_state)
        d <- sum(!target_level & !target_state)
        contingency <- matrix(c(a, b, c_value, d), nrow = 2L, byrow = TRUE)
        evaluable <- sum(contingency[1L, ]) > 0L && sum(contingency[2L, ]) > 0L &&
          sum(contingency[, 1L]) > 0L && sum(contingency[, 2L]) > 0L
        fisher_result <- if (evaluable) fisher.test(contingency, alternative = "greater") else NULL

        fisher_index <- fisher_index + 1L
        fisher_parts[[fisher_index]] <- data.table(
          Axis = axis_name,
          Level = level_name,
          CellType = cell_type,
          State = state_name,
          CellState = paste(cell_type, state_name, sep = "_"),
          SamplesTargetLevelState = a,
          SamplesTargetLevelOtherStates = b,
          SamplesOtherLevelsState = c_value,
          SamplesOtherLevelsOtherStates = d,
          OddsRatio = if (evaluable) unname(fisher_result$estimate) else NA_real_,
          P_one_sided = if (evaluable) fisher_result$p.value else NA_real_,
          AssignedSamples = nrow(parent),
          TargetLevelSamples = sum(target_level),
          Donors = uniqueN(parent$PatientID),
          Status = if (evaluable) "EVALUABLE" else "NOT_EVALUABLE_FISHER_MARGIN"
        )
      }
    }
  }
  cat("DISCOVERY_REPRESENTATION_AXIS_COMPLETE", axis_name, "\n")
}

fisher_results <- rbindlist(fisher_parts, use.names = TRUE)
fisher_results[, FDR_BH_across_cell_states := p.adjust(P_one_sided, method = "BH"), by = .(Axis, Level)]
fisher_results[, EnrichmentScore := fifelse(
  Status == "EVALUABLE" & is.finite(FDR_BH_across_cell_states),
  pmin(-log10(pmax(FDR_BH_across_cell_states, .Machine$double.xmin)), 10),
  NA_real_
)]
fisher_results[, SignificantEnrichment := Status == "EVALUABLE" & OddsRatio > 1 & FDR_BH_across_cell_states < 0.05]
setorder(fisher_results, Axis, Level, CellType, State)

# Step4. 保存正式结果和分母
axis_denominators <- assignments[, .(
  AssignedSamples = .N,
  Donors = uniqueN(PatientID),
  FormalStates = uniqueN(State)
), by = .(CellType)]
phenotype_denominators <- rbindlist(lapply(axes, function(axis_name) {
  assignments[, .(Samples = uniqueN(ID), Donors = uniqueN(PatientID)), by = .(CellType, Level = get(axis_name))][, Axis := axis_name]
}))

fwrite(fisher_results, file.path(output_dir, "discovery_state_representation_fisher.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(fisher_results[SignificantEnrichment == TRUE], file.path(output_dir, "discovery_state_representation_significant.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(axis_denominators, file.path(output_dir, "discovery_cell_type_assignment_denominators.tsv"), sep = "\t", quote = FALSE)
fwrite(phenotype_denominators, file.path(output_dir, "discovery_phenotype_denominators.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_DISCOVERY_STATE_REPRESENTATION_CELL2021",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  discovery = discovery_name,
  parent_cell_types = nrow(rank_data),
  formal_states = nrow(formal_states),
  axes = axes,
  primary_axes = c("Disease", "Region3"),
  secondary_axes = "Inflammation",
  unit = "discovery bulk sample",
  method = "dominant discovery state; one-sided Fisher greater; BH within phenotype level across all cell states",
  non_positive_rows_retained = TRUE
)
jsonlite::write_json(completion, file.path(output_dir, "discovery_state_representation_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(fisher_results[, .(
  Evaluable = sum(Status == "EVALUABLE"),
  Significant = sum(SignificantEnrichment),
  MedianScore = median(EnrichmentScore, na.rm = TRUE),
  MaximumScore = max(EnrichmentScore, na.rm = TRUE)
), by = Axis])
cat("PASS_DISCOVERY_STATE_REPRESENTATION_CELL2021\n")
