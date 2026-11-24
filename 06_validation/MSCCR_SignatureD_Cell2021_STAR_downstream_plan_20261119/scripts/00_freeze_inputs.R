#!/usr/bin/env Rscript
# 本脚本在服务器只读冻结已验收的正式 CIBERSORTx/EcoTyper 输入、模型和运行身份。
# 它不重算模型、不改源文件；缺失或不一致时停止并保留已写出的审计表。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取固定路径并确认正式状态
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 00_freeze_inputs.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite", "yaml")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)

audit_dir <- file.path(output_root, "00_freeze/formal_model")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  status_file, runtime_config_file, initial_config_file, fractions_file, cell_types_file,
  annotation_file, mixture_tpm_file, rank_data_file,
  file.path(ecotype_discovery_root, "ecotypes.txt"),
  file.path(ecotype_discovery_root, "mapping_all_states.txt"),
  file.path(ecotype_discovery_root, "jaccard_matrix.txt"),
  file.path(ecotype_discovery_root, "ecotype_abundance.txt"),
  file.path(ecotype_discovery_root, "ecotype_assignment.txt")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) stop("Missing formal input: ", paste(missing_files, collapse = "; "), call. = FALSE)

status <- readLines(status_file, warn = FALSE)
if (!"status=ACCEPTED" %in% status) stop("Formal EcoTyper task is not ACCEPTED", call. = FALSE)

runtime_config <- yaml::read_yaml(runtime_config_file)
initial_config <- yaml::read_yaml(initial_config_file)
runtime_skip <- unlist(runtime_config$default$`Pipeline settings`$`Pipeline steps to skip`)
initial_skip <- unlist(initial_config$default$`Pipeline settings`$`Pipeline steps to skip`)
if (!7L %in% as.integer(initial_skip)) stop("Initial formal config no longer records step 7 as skipped", call. = FALSE)

cat("FORMAL_STATUS", paste(status, collapse = "; "), "\n")
cat("CONFIG_SKIP initial=", paste(initial_skip, collapse = ","), " runtime=", paste(runtime_skip, collapse = ","), "\n", sep = "")

# Step2. 核对样本、父类型、Fractions 和 TPM 表头
annotation <- data.table::fread(annotation_file, colClasses = "character", check.names = FALSE)
annotation_required <- c("ID", "PatientID", "Disease", "Inflammation", "Region3")
if (!all(annotation_required %in% names(annotation))) stop("Annotation columns changed", call. = FALSE)
if (nrow(annotation) != expected_samples || anyNA(annotation$ID) || anyDuplicated(annotation$ID)) stop("Annotation identity failed", call. = FALSE)

cell_types <- readLines(cell_types_file, warn = FALSE)
cell_types <- cell_types[nzchar(cell_types)]
if (length(cell_types) != expected_cell_types || anyDuplicated(cell_types)) stop("Cell-type identity failed", call. = FALSE)

fractions <- data.table::fread(fractions_file, check.names = FALSE)
if (!"Mixture" %in% names(fractions)) stop("Fractions file lacks Mixture", call. = FALSE)
if (!identical(as.character(fractions$Mixture), annotation$ID)) stop("Fractions sample order differs from annotation", call. = FALSE)
metric_columns <- intersect(c("P-value", "Correlation", "RMSE", "Absolute score"), names(fractions))
fraction_types <- setdiff(names(fractions), c("Mixture", metric_columns))
if (!setequal(fraction_types, cell_types)) stop("Fractions and frozen cell types differ", call. = FALSE)

fraction_values <- as.matrix(fractions[, ..fraction_types])
storage.mode(fraction_values) <- "double"
if (any(!is.finite(fraction_values)) || any(fraction_values < 0)) stop("Invalid fraction values", call. = FALSE)
fraction_row_sums <- rowSums(fraction_values)

tpm_header <- names(data.table::fread(mixture_tpm_file, nrows = 0L, check.names = FALSE))
if (!identical(tpm_header[-1L], annotation$ID)) stop("TPM sample order differs from annotation", call. = FALSE)

sample_audit <- data.frame(
  Object = c("Annotation", "Fractions", "TPM"),
  Samples = c(nrow(annotation), nrow(fractions), length(tpm_header) - 1L),
  ExactOrder = c(TRUE, identical(as.character(fractions$Mixture), annotation$ID), identical(tpm_header[-1L], annotation$ID))
)
print(sample_audit)

# Step3. 固定 rank 和每个父类型的最终 state 文件
rank_data <- data.table::fread(rank_data_file, colClasses = "character", check.names = FALSE)
if (ncol(rank_data) < 2L) stop("rank_data.txt has fewer than two columns", call. = FALSE)
rank_cell_type <- rank_data[[1L]]
rank_selected <- trimws(rank_data[[2L]])
if (!setequal(rank_cell_type, cell_types) || anyDuplicated(rank_cell_type)) stop("rank_data cell types differ", call. = FALSE)

state_rows <- vector("list", length(cell_types))
state_files <- character()
for (i in seq_along(cell_types)) {
  cell_type <- cell_types[[i]]
  selected_rank <- rank_selected[match(cell_type, rank_cell_type)]
  state_dir <- file.path(state_discovery_root, cell_type, selected_rank)
  files <- file.path(state_dir, c(
    "W.txt", "H.txt", "initial_gene_info.txt", "initial_state_abundances.txt", "initial_state_assignment.txt",
    "first_filter_data.txt", "gene_info.txt", "mapping_to_initial_states.txt", "state_abundances.txt", "state_assignment.txt"
  ))
  if (any(!file.exists(files))) stop("Missing selected state file for ", cell_type, ": ", paste(files[!file.exists(files)], collapse = "; "), call. = FALSE)

  initial_abundance <- as.matrix(read.delim(file.path(state_dir, "initial_state_abundances.txt"), row.names = 1L, check.names = FALSE))
  final_abundance <- as.matrix(read.delim(file.path(state_dir, "state_abundances.txt"), row.names = 1L, check.names = FALSE))
  storage.mode(initial_abundance) <- "double"
  storage.mode(final_abundance) <- "double"
  if (!identical(colnames(initial_abundance), annotation$ID) || !identical(colnames(final_abundance), annotation$ID)) stop("State sample order differs: ", cell_type, call. = FALSE)
  if (any(!is.finite(final_abundance)) || any(final_abundance < 0)) stop("Invalid state abundance: ", cell_type, call. = FALSE)

  initial_assignment <- data.table::fread(file.path(state_dir, "initial_state_assignment.txt"), colClasses = "character", check.names = FALSE)
  final_assignment <- data.table::fread(file.path(state_dir, "state_assignment.txt"), colClasses = "character", check.names = FALSE)
  if (!all(c("ID", "State") %in% names(initial_assignment)) || !all(c("ID", "State") %in% names(final_assignment))) stop("State assignment columns changed: ", cell_type, call. = FALSE)
  if (anyDuplicated(initial_assignment$ID) || anyDuplicated(final_assignment$ID)) stop("Duplicated state assignment IDs: ", cell_type, call. = FALSE)
  if (!all(initial_assignment$ID %in% annotation$ID) || !all(final_assignment$ID %in% annotation$ID)) stop("Unknown state assignment ID: ", cell_type, call. = FALSE)

  mapping <- data.table::fread(file.path(state_dir, "mapping_to_initial_states.txt"), colClasses = "character", check.names = FALSE)
  state_rows[[i]] <- data.frame(
    CellType = cell_type,
    SelectedRank = selected_rank,
    InitialStates = nrow(initial_abundance),
    RetainedStates = nrow(final_abundance),
    InitialAssignedSamples = nrow(initial_assignment),
    RetainedAssignedSamples = nrow(final_assignment),
    MappingRows = nrow(mapping),
    AdvancedQCBackupPresent = file.exists(file.path(state_dir, "state_assignment_first_filter.txt"))
  )
  state_files <- c(state_files, files)
  cat("STATE_FROZEN", cell_type, "rank", selected_rank, "initial", nrow(initial_abundance), "retained", nrow(final_abundance), "\n")
}
state_audit <- data.table::rbindlist(state_rows)

# Step4. 核对 ecotype 成员、abundance、assignment 和并列最大值
ecotypes_file <- file.path(ecotype_discovery_root, "ecotypes.txt")
ecotype_abundance_file <- file.path(ecotype_discovery_root, "ecotype_abundance.txt")
ecotype_assignment_file <- file.path(ecotype_discovery_root, "ecotype_assignment.txt")
members <- data.table::fread(ecotypes_file, colClasses = "character", check.names = FALSE)
if (!all(c("CellType", "State", "ID", "Ecotype") %in% names(members))) stop("ecotypes.txt columns changed", call. = FALSE)
if (anyNA(members$ID) || anyDuplicated(members$ID)) stop("Invalid ecotype member IDs", call. = FALSE)

ecotype_abundance <- as.matrix(read.delim(ecotype_abundance_file, row.names = 1L, check.names = FALSE))
storage.mode(ecotype_abundance) <- "double"
if (!identical(colnames(ecotype_abundance), annotation$ID)) stop("Ecotype sample order differs", call. = FALSE)
if (any(!is.finite(ecotype_abundance)) || any(ecotype_abundance < 0)) stop("Invalid ecotype abundance", call. = FALSE)

ecotype_max <- apply(ecotype_abundance, 2L, max)
ecotype_ties <- colSums(sweep(ecotype_abundance, 2L, ecotype_max, "=="))
tie_audit <- data.frame(ID = colnames(ecotype_abundance), MaxEcotypeAbundance = ecotype_max, MaxTieCount = ecotype_ties)

ecotype_assignment <- data.table::fread(ecotype_assignment_file, colClasses = "character", check.names = FALSE)
if (!all(c("ID", "Ecotype") %in% names(ecotype_assignment))) stop("Ecotype assignment columns changed", call. = FALSE)
if (anyDuplicated(ecotype_assignment$ID) || !all(ecotype_assignment$ID %in% annotation$ID)) stop("Invalid ecotype assignment IDs", call. = FALSE)

# Step5. 记录高级 QC 是否实际存在，不自动补算或修改模型
dropout_files <- c(
  file.path(advanced_qc_root, "dropout_scores.txt"),
  file.path(advanced_qc_root, "dropout_scores_with_outliers.txt")
)
advanced_qc_complete <- all(file.exists(dropout_files)) && all(state_audit$AdvancedQCBackupPresent)
qc_status <- data.frame(
  Component = c("Marker_and_AFI_filter", "Advanced_dropout_QC", "Literal_Cell2021_three_filter_QC"),
  Status = c("PRESENT_IN_SELECTED_OUTPUTS", if (advanced_qc_complete) "COMPLETE" else "SKIPPED_OR_INCOMPLETE", if (advanced_qc_complete) "COMPLETE" else "NOT_COMPLETE"),
  ChangesApplied = c(TRUE, advanced_qc_complete, advanced_qc_complete),
  stringsAsFactors = FALSE
)
print(qc_status)

# Step6. 保存关键文件 SHA256、模型表和阶段状态
ecotype_files <- c(
  file.path(ecotype_discovery_root, c(
    "ecotypes.txt", "mapping_all_states.txt", "jaccard_matrix.txt", "ecotype_abundance.txt", "ecotype_assignment.txt"
  ))
)
input_paths <- unique(c(required_files, state_files, ecotype_files, dropout_files[file.exists(dropout_files)]))
input_paths <- normalizePath(input_paths, mustWork = TRUE)
manifest <- data.frame(
  Path = input_paths,
  Bytes = file.info(input_paths)$size,
  Modified = format(file.info(input_paths)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
  SHA256 = vapply(input_paths, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)

fraction_audit <- data.frame(
  Samples = nrow(fractions),
  CellTypes = length(fraction_types),
  MetricColumns = paste(metric_columns, collapse = ";"),
  RowSumMin = min(fraction_row_sums),
  RowSumMedian = median(fraction_row_sums),
  RowSumMax = max(fraction_row_sums)
)

data.table::fwrite(manifest, file.path(audit_dir, "formal_input_manifest.tsv"), sep = "\t")
data.table::fwrite(sample_audit, file.path(audit_dir, "sample_identity_audit.tsv"), sep = "\t")
data.table::fwrite(fraction_audit, file.path(audit_dir, "fraction_identity_audit.tsv"), sep = "\t")
data.table::fwrite(state_audit, file.path(audit_dir, "state_identity_audit.tsv"), sep = "\t")
data.table::fwrite(members, file.path(audit_dir, "ecotype_members_frozen.tsv"), sep = "\t")
data.table::fwrite(tie_audit, file.path(audit_dir, "ecotype_max_tie_audit.tsv"), sep = "\t")
data.table::fwrite(qc_status, file.path(audit_dir, "state_qc_stage_status.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo.txt"))

summary <- list(
  status = "PASS_FORMAL_IDENTITY",
  model_id = paste(dataset_name, fractions_name, sep = "/"),
  accepted = TRUE,
  samples = nrow(annotation),
  patients = length(unique(annotation$PatientID)),
  cell_types = length(cell_types),
  initial_states = sum(state_audit$InitialStates),
  retained_states = sum(state_audit$RetainedStates),
  ecotype_member_states = nrow(members),
  ecotypes = length(unique(members$Ecotype)),
  ecotype_assigned_samples = nrow(ecotype_assignment),
  dominant_tied_samples = sum(ecotype_ties > 1L),
  advanced_dropout_qc = if (advanced_qc_complete) "COMPLETE" else "SKIPPED_OR_INCOMPLETE",
  literal_cell2021_three_filter_qc = if (advanced_qc_complete) "COMPLETE" else "NOT_COMPLETE",
  downstream_release = if (advanced_qc_complete) "G0_PASS_G2_PENDING" else "G0_PASS_G2_BLOCKED",
  model_changed = FALSE
)
jsonlite::write_json(summary, file.path(audit_dir, "formal_freeze_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat("FORMAL_FREEZE_PASS", summary$samples, "samples", summary$retained_states, "retained states", summary$ecotypes, "ecotypes; advanced QC", summary$advanced_dropout_qc, "\n")
