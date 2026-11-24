#!/usr/bin/env Rscript
# 本脚本汇总 MSCCR Control 在三个部位的 cell state 和 ecotype 覆盖，仅描述恢复结果，不进行表型显著性检验。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入、模型和输出
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 50_summarize_MSCCR_Control_bulk_recovery_current85.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

project_root <- normalizePath(file.path(RUN_ROOT, "../.."), mustWork = TRUE)
discovery_name <- "MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101"
model_root <- file.path(project_root, "software/ecotyper/EcoTyper", discovery_name, "Custom")
recovery_root <- file.path(RUN_ROOT, "output/36_MSCCR_Control_bulk_recovery_current85_retry1")
input_root <- file.path(
  RUN_ROOT,
  "output/35_MSCCR_Control_CPM_current85_preflight/reconstructed_current_model_gene_union/input"
)
output_dir <- file.path(RUN_ROOT, "output/37_MSCCR_Control_bulk_recovery_summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

regions <- c("Ileum", "Colon", "Rectum")
rank_data <- fread(file.path(model_root, "Analysis/rank_selection/rank_data.txt"))
if (nrow(rank_data) != 20L) stop("Current model is no longer the frozen 20-parent-type model", call. = FALSE)

formal_state_parts <- vector("list", nrow(rank_data))
for (cell_index in seq_len(nrow(rank_data))) {
  discovery_assignment_file <- file.path(
    model_root, "Cell_States/discovery", rank_data$CellType[[cell_index]],
    as.character(rank_data$Chosen_Rank[[cell_index]]), "state_assignment.txt"
  )
  discovery_assignment <- suppressWarnings(fread(discovery_assignment_file, colClasses = "character"))
  if ("V1" %in% names(discovery_assignment) && "ID" %in% names(discovery_assignment)) {
    discovery_row_index <- suppressWarnings(as.integer(discovery_assignment$V1))
    if (anyNA(discovery_row_index) || anyDuplicated(discovery_row_index)) stop("Discovery row index is invalid: ", rank_data$CellType[[cell_index]], call. = FALSE)
    discovery_assignment[, V1 := NULL]
  }
  formal_state_parts[[cell_index]] <- data.table(
    CellType = rank_data$CellType[[cell_index]],
    ChosenRank = rank_data$Chosen_Rank[[cell_index]],
    State = sort(unique(discovery_assignment$State[!is.na(discovery_assignment$State)]))
  )
}
formal_states <- rbindlist(formal_state_parts)
if (nrow(formal_states) != 85L) stop("Current model is no longer the frozen 85-state model after QC", call. = FALSE)

# Step2. 汇总 60 个 cell state assignment 文件
state_assignment_parts <- vector("list", nrow(rank_data) * length(regions))
state_index <- 0L

for (region in regions) {
  annotation_file <- file.path(input_root, paste0("MSCCR_Control_CPM_current85_", region, "_annotation.tsv"))
  annotation <- fread(annotation_file, colClasses = "character")
  expected_ids <- annotation$ID

  for (cell_index in seq_len(nrow(rank_data))) {
    cell_type <- rank_data$CellType[[cell_index]]
    state_file <- file.path(
      recovery_root, region, "RecoveryOutput", paste0("MSCCR_Control_CPM_current85_", region),
      cell_type, "state_assignment.txt"
    )
    if (!file.exists(state_file)) stop("Missing state assignment: ", region, " ", cell_type, call. = FALSE)
    state_assignment <- suppressWarnings(fread(state_file, colClasses = "character"))
    if ("V1" %in% names(state_assignment) && "ID" %in% names(state_assignment)) {
      if (!identical(state_assignment$V1, state_assignment$ID)) stop("Recovery row names do not match ID: ", region, " ", cell_type, call. = FALSE)
      state_assignment[, V1 := NULL]
    }
    if (!all(c("ID", "State") %in% names(state_assignment))) stop("State assignment columns changed: ", region, " ", cell_type, call. = FALSE)
    if (anyDuplicated(state_assignment$ID)) stop("Duplicated state assignment IDs: ", region, " ", cell_type, call. = FALSE)
    if (any(!state_assignment$ID %in% expected_ids)) stop("Unexpected state assignment IDs: ", region, " ", cell_type, call. = FALSE)

    state_index <- state_index + 1L
    state_assignment_parts[[state_index]] <- state_assignment[, .(
      Region3 = region,
      CellType = cell_type,
      ID,
      State
    )]
  }
}

state_assignments <- rbindlist(state_assignment_parts)
state_grid <- CJ(Region3 = regions, CellType = rank_data$CellType, unique = TRUE)
state_grid <- merge(state_grid, formal_states, by = "CellType", allow.cartesian = TRUE)

state_counts <- state_assignments[, .(AssignedSamples = uniqueN(ID)), by = .(Region3, CellType, State)]
state_counts <- merge(state_grid, state_counts, by = c("Region3", "CellType", "State"), all.x = TRUE)
state_counts[is.na(AssignedSamples), AssignedSamples := 0L]

sample_denominators <- rbindlist(lapply(regions, function(region) {
  annotation <- fread(file.path(input_root, paste0("MSCCR_Control_CPM_current85_", region, "_annotation.tsv")), colClasses = "character")
  data.table(Region3 = region, Samples = nrow(annotation), Donors = uniqueN(annotation$PatientID))
}))
state_counts <- merge(state_counts, sample_denominators, by = "Region3")
state_counts[, AssignmentFraction := AssignedSamples / Samples]
state_counts[, RegionOrder := match(Region3, regions)]
state_counts[, CellTypeOrder := match(CellType, rank_data$CellType)]
setorder(state_counts, RegionOrder, CellTypeOrder, State)
state_counts[, c("RegionOrder", "CellTypeOrder") := NULL]

state_parent_summary <- state_counts[, .(
  FormalStates = .N,
  ObservedStates = sum(AssignedSamples > 0L),
  CompletelyUnobservedStates = sum(AssignedSamples == 0L),
  MinimumStateSamples = min(AssignedSamples),
  MaximumStateSamples = max(AssignedSamples)
), by = .(Region3, CellType)]
state_overall_summary <- state_counts[, .(
  FormalStates = .N,
  ObservedStates = sum(AssignedSamples > 0L),
  CompletelyUnobservedStates = sum(AssignedSamples == 0L),
  StatesAtLeastFiveSamples = sum(AssignedSamples >= 5L),
  StatesAtLeastTenPercent = sum(AssignmentFraction >= 0.10)
), by = Region3]

# Step3. 汇总三个部位的 ecotype assignment
ecotype_assignment_parts <- vector("list", length(regions))
ecotype_missing_parts <- vector("list", length(regions))

for (region_index in seq_along(regions)) {
  region <- regions[[region_index]]
  annotation <- fread(file.path(input_root, paste0("MSCCR_Control_CPM_current85_", region, "_annotation.tsv")), colClasses = "character")
  ecotype_file <- file.path(
    model_root, "Ecotypes/recovery", paste0("MSCCR_Control_CPM_current85_", region),
    "initial_ecotype_assignment.txt"
  )
  ecotype_assignment <- suppressWarnings(fread(ecotype_file, colClasses = "character"))
  if ("V1" %in% names(ecotype_assignment) && "ID" %in% names(ecotype_assignment)) {
    ecotype_row_index <- suppressWarnings(as.integer(ecotype_assignment$V1))
    if (anyNA(ecotype_row_index) || anyDuplicated(ecotype_row_index)) stop("Initial ecotype row index is invalid: ", region, call. = FALSE)
    ecotype_assignment[, V1 := NULL]
  }
  required_columns <- c("ID", "Ecotype", "AssignedToEcotypeStates", "AssignmentP", "AssignmentQ")
  if (!all(required_columns %in% names(ecotype_assignment))) stop("Ecotype assignment columns changed: ", region, call. = FALSE)
  if (anyDuplicated(ecotype_assignment$ID)) stop("Duplicated ecotype assignment IDs: ", region, call. = FALSE)
  if (any(!ecotype_assignment$ID %in% annotation$ID)) stop("Unexpected ecotype assignment IDs: ", region, call. = FALSE)

  ecotype_assignment_parts[[region_index]] <- ecotype_assignment[, .(
    Region3 = region,
    ID,
    Ecotype,
    AssignedToEcotypeStates = as.logical(AssignedToEcotypeStates),
    AssignmentP = as.numeric(AssignmentP),
    AssignmentQ = as.numeric(AssignmentQ)
  )]
  ecotype_missing_parts[[region_index]] <- data.table(Region3 = region, ID = setdiff(annotation$ID, ecotype_assignment$ID))
}

ecotype_assignments <- rbindlist(ecotype_assignment_parts)
ecotype_missing <- rbindlist(ecotype_missing_parts)
formal_ecotypes <- paste0("E", seq_len(8L))
ecotype_grid <- CJ(Region3 = regions, Ecotype = formal_ecotypes, unique = TRUE)
ecotype_counts <- ecotype_assignments[!is.na(Ecotype) & Ecotype %in% formal_ecotypes, .(AssignedSamples = uniqueN(ID)), by = .(Region3, Ecotype)]
ecotype_counts <- merge(ecotype_grid, ecotype_counts, by = c("Region3", "Ecotype"), all.x = TRUE)
ecotype_counts[is.na(AssignedSamples), AssignedSamples := 0L]
ecotype_counts <- merge(ecotype_counts, sample_denominators, by = "Region3")
ecotype_counts[, AssignmentFraction := AssignedSamples / Samples]
ecotype_counts[, RegionOrder := match(Region3, regions)]
setorder(ecotype_counts, RegionOrder, Ecotype)
ecotype_counts[, RegionOrder := NULL]

ecotype_region_summary <- ecotype_assignments[, .(
  InitialAssignmentRows = .N,
  AssignedRows = sum(Ecotype %in% formal_ecotypes),
  UnassignedRows = sum(Ecotype == "Unassigned" | is.na(Ecotype)),
  AssignedToEcotypeStatesRows = sum(AssignedToEcotypeStates, na.rm = TRUE),
  SignificantQ025 = sum(AssignmentQ < 0.25, na.rm = TRUE),
  SignificantQ005 = sum(AssignmentQ < 0.05, na.rm = TRUE),
  ObservedEcotypes = uniqueN(Ecotype[!is.na(Ecotype) & Ecotype %in% formal_ecotypes])
), by = Region3]
ecotype_region_summary <- merge(sample_denominators, ecotype_region_summary, by = "Region3", all.x = TRUE)
ecotype_region_summary[, MissingInitialRows := Samples - InitialAssignmentRows]
ecotype_region_summary[, InitialOutputCoverage := InitialAssignmentRows / Samples]
ecotype_region_summary[, RegionOrder := match(Region3, regions)]
setorder(ecotype_region_summary, RegionOrder)
ecotype_region_summary[, RegionOrder := NULL]

# Step4. 保存覆盖统计和用途边界
fwrite(state_assignments, file.path(output_dir, "MSCCR_Control_cell_state_assignments_long.tsv"), sep = "\t", quote = FALSE)
fwrite(state_counts, file.path(output_dir, "MSCCR_Control_cell_state_assignment_counts.tsv"), sep = "\t", quote = FALSE)
fwrite(state_parent_summary, file.path(output_dir, "MSCCR_Control_cell_state_parent_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(state_overall_summary, file.path(output_dir, "MSCCR_Control_cell_state_overall_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(ecotype_assignments, file.path(output_dir, "MSCCR_Control_ecotype_assignments_long.tsv"), sep = "\t", quote = FALSE)
fwrite(ecotype_missing, file.path(output_dir, "MSCCR_Control_ecotype_missing_samples.tsv"), sep = "\t", quote = FALSE)
fwrite(ecotype_counts, file.path(output_dir, "MSCCR_Control_ecotype_assignment_counts.tsv"), sep = "\t", quote = FALSE)
fwrite(ecotype_region_summary, file.path(output_dir, "MSCCR_Control_ecotype_region_summary.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_MSCCR_CONTROL_BULK_RECOVERY_SUMMARY",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  samples = setNames(sample_denominators$Samples, sample_denominators$Region3),
  donors = setNames(sample_denominators$Donors, sample_denominators$Region3),
  formal_states = nrow(formal_states),
  formal_ecotypes = length(formal_ecotypes),
  role = "descriptive internal nonIBD coverage; no phenotype enrichment p-value is calculated from Control-only samples",
  limitation = "MSCCR Control belongs to the discovery cohort source and is not an independent external validation cohort"
)
jsonlite::write_json(completion, file.path(output_dir, "MSCCR_Control_bulk_recovery_summary_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(state_overall_summary)
print(ecotype_region_summary)
cat("PASS_MSCCR_CONTROL_BULK_RECOVERY_SUMMARY\n")
