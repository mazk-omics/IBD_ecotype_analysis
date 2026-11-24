suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript 43_audit_GSE57945_bulk_recovery.R <result_root> <model_root> <expression_file> <annotation_file> <output_dir>")
}

result_root <- normalizePath(args[1], mustWork = TRUE)
model_root <- normalizePath(args[2], mustWork = TRUE)
expression_file <- normalizePath(args[3], mustWork = TRUE)
annotation_file <- normalizePath(args[4], mustWork = TRUE)
output_dir <- args[5]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step1. 冻结输入样本、模型 rank 和 phenotype annotation。
expression_header <- fread(expression_file, nrows = 0, check.names = FALSE)
sample_ids <- names(expression_header)[-1]
annotation <- fread(annotation_file, check.names = FALSE)
rank_data <- fread(file.path(model_root, "Analysis", "rank_selection", "rank_data.txt"))

if (!identical(sample_ids, annotation$ID)) {
  stop("Expression and annotation sample order does not match")
}
if (nrow(rank_data) != 20) {
  stop("Current model does not contain 20 selected cell types")
}

# Step2. 逐父类型核对 state abundance、assignment 和冻结 state 数。
cell_type_audit <- vector("list", nrow(rank_data))
assignment_by_disease <- vector("list", nrow(rank_data))

for (i in seq_len(nrow(rank_data))) {
  cell_type <- rank_data$CellType[i]
  chosen_rank <- rank_data$Chosen_Rank[i]
  abundance_file <- file.path(result_root, cell_type, "state_abundances.txt")
  assignment_file <- file.path(result_root, cell_type, "state_assignment.txt")
  mapping_file <- file.path(
    model_root,
    "Cell_States",
    "discovery",
    cell_type,
    as.character(chosen_rank),
    "mapping_to_initial_states.txt"
  )

  abundance <- read.delim(abundance_file, check.names = FALSE)
  assignment <- fread(assignment_file, check.names = FALSE)
  mapping <- fread(mapping_file, check.names = FALSE)
  expected_states <- unique(as.character(mapping$State))

  sample_match <- identical(names(abundance), sample_ids)
  state_match <- setequal(rownames(abundance), expected_states)
  finite_values <- all(is.finite(as.matrix(abundance)))
  nonnegative_values <- all(as.matrix(abundance) >= 0)
  assignment_ids_valid <- all(assignment$ID %in% sample_ids) && anyDuplicated(assignment$ID) == 0
  assignment_states_valid <- all(assignment$State %in% expected_states)

  cell_type_audit[[i]] <- data.table(
    CellType = cell_type,
    ChosenRank = chosen_rank,
    ExpectedStates = length(expected_states),
    OutputStates = nrow(abundance),
    OutputSamples = ncol(abundance),
    AssignedSamples = nrow(assignment),
    SampleOrderPass = sample_match,
    StateIdentityPass = state_match,
    FinitePass = finite_values,
    NonnegativePass = nonnegative_values,
    AssignmentIDPass = assignment_ids_valid,
    AssignmentStatePass = assignment_states_valid
  )

  assignment_summary <- merge(
    assignment[, .(ID, State)],
    annotation[, .(ID, Disease)],
    by = "ID",
    all.x = TRUE
  )[, .N, by = .(Disease, State)]
  assignment_summary[, CellType := cell_type]
  setcolorder(assignment_summary, c("CellType", "Disease", "State", "N"))
  assignment_by_disease[[i]] <- assignment_summary
}

cell_type_audit <- rbindlist(cell_type_audit)
assignment_by_disease <- rbindlist(assignment_by_disease)

# Step3. 核对 8 个 ecotype 的 abundance 和 assignment。
ecotype_abundance_file <- file.path(result_root, "Ecotypes", "ecotype_abundance.txt")
ecotype_assignment_file <- file.path(result_root, "Ecotypes", "ecotype_assignment.txt")
ecotype_abundance <- read.delim(ecotype_abundance_file, check.names = FALSE)
ecotype_assignment <- as.data.table(read.delim(ecotype_assignment_file, check.names = FALSE))

ecotype_ids_valid <- all(ecotype_assignment$ID %in% sample_ids) && anyDuplicated(ecotype_assignment$ID) == 0
ecotype_audit <- data.table(
  ExpectedEcotypes = 8,
  OutputEcotypes = nrow(ecotype_abundance),
  OutputSamples = ncol(ecotype_abundance),
  AssignedSamples = nrow(ecotype_assignment),
  SampleOrderPass = identical(names(ecotype_abundance), sample_ids),
  FinitePass = all(is.finite(as.matrix(ecotype_abundance))),
  NonnegativePass = all(as.matrix(ecotype_abundance) >= 0),
  AssignmentIDPass = ecotype_ids_valid
)

ecotype_assignment_by_disease <- merge(
  ecotype_assignment[, .(ID, Ecotype)],
  annotation[, .(ID, Disease)],
  by = "ID",
  all.x = TRUE
)[, .N, by = .(Disease, Ecotype)]

# Step4. 写出验收状态和完整分母。
cell_pass_columns <- c(
  "SampleOrderPass",
  "StateIdentityPass",
  "FinitePass",
  "NonnegativePass",
  "AssignmentIDPass",
  "AssignmentStatePass"
)
cell_pass <- all(unlist(cell_type_audit[, ..cell_pass_columns])) &&
  all(cell_type_audit$ExpectedStates == cell_type_audit$OutputStates) &&
  all(cell_type_audit$OutputSamples == length(sample_ids)) &&
  sum(cell_type_audit$OutputStates) == 85
ecotype_pass <- ecotype_audit$OutputEcotypes == 8 &&
  ecotype_audit$OutputSamples == length(sample_ids) &&
  all(unlist(ecotype_audit[, .(SampleOrderPass, FinitePass, NonnegativePass, AssignmentIDPass)]))

status <- if (cell_pass && ecotype_pass) "ACCEPTED" else "FAILED"
detail <- if (status == "ACCEPTED") {
  "Official bulk recovery contains 20 cell types, 85 states, 8 ecotypes and 322 ordered samples"
} else {
  "At least one state/ecotype output identity or dimension check failed"
}

fwrite(cell_type_audit, file.path(output_dir, "cell_type_output_audit.tsv"), sep = "\t")
fwrite(assignment_by_disease, file.path(output_dir, "state_assignment_by_disease.tsv"), sep = "\t")
fwrite(ecotype_audit, file.path(output_dir, "ecotype_output_audit.tsv"), sep = "\t")
fwrite(ecotype_assignment_by_disease, file.path(output_dir, "ecotype_assignment_by_disease.tsv"), sep = "\t")
writeLines(c(paste0("status=", status), paste0("detail=", detail)), file.path(output_dir, "audit_status.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

print(cell_type_audit)
print(ecotype_audit)
cat(paste0("status=", status, "\n"))
if (status != "ACCEPTED") {
  quit(status = 2)
}
