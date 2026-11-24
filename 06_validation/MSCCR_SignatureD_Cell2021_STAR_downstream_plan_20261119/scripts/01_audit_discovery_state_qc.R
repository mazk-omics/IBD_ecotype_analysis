#!/usr/bin/env Rscript
# 本脚本只读复核正式模型的 marker、AFI 和 dropout 三层 state QC，并估计补做 Step 7 的影响。
# 它不删除 state、不改 assignment、不重聚类 ecotype；所有判断写入独立审计目录。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取 Phase 0 冻结结果和正式路径
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 01_audit_discovery_state_qc.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)

freeze_file <- file.path(output_root, "00_freeze/formal_model/formal_freeze_summary.json")
if (!file.exists(freeze_file)) stop("Run 00_freeze_inputs.R first", call. = FALSE)
freeze <- jsonlite::read_json(freeze_file, simplifyVector = TRUE)
if (!identical(freeze$status, "PASS_FORMAL_IDENTITY")) stop("Formal input freeze did not pass", call. = FALSE)

audit_dir <- file.path(output_root, "01_state_qc_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

rank_data <- data.table::fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_cell_type <- rank_data[[1L]]
rank_selected <- trimws(rank_data[[2L]])
members <- data.table::fread(file.path(ecotype_discovery_root, "ecotypes.txt"), colClasses = "character", check.names = FALSE)

# Step2. 独立复核 marker 数和 AFI 筛选
first_filter_rows <- vector("list", length(rank_cell_type))
first_filter_summary <- vector("list", length(rank_cell_type))

for (i in seq_along(rank_cell_type)) {
  cell_type <- rank_cell_type[[i]]
  selected_rank <- rank_selected[[i]]
  state_dir <- file.path(state_discovery_root, cell_type, selected_rank)
  filter_file <- file.path(state_dir, "first_filter_data.txt")
  mapping_file <- file.path(state_dir, "mapping_to_initial_states.txt")
  if (!file.exists(filter_file) || !file.exists(mapping_file)) stop("Missing first-filter evidence: ", cell_type, call. = FALSE)

  first_filter <- data.table::fread(filter_file, check.names = FALSE)
  mapping <- data.table::fread(mapping_file, colClasses = "character", check.names = FALSE)
  if (!all(c("State", "NumberMarkers", "PosRatio") %in% names(first_filter))) stop("first_filter_data columns changed: ", cell_type, call. = FALSE)
  if (!all(c("State", "InitialState") %in% names(mapping))) stop("mapping columns changed: ", cell_type, call. = FALSE)

  marker_keep <- as.character(first_filter[NumberMarkers > 10, State])
  marker_fallback <- length(marker_keep) < 2L
  if (marker_fallback) marker_keep <- as.character(first_filter$State)

  afi_keep <- as.character(first_filter[PosRatio < 1, State])
  afi_fallback <- length(intersect(marker_keep, afi_keep)) < 2L
  expected_keep <- if (afi_fallback) marker_keep else intersect(marker_keep, afi_keep)
  observed_keep <- as.character(mapping$InitialState)
  filter_match <- setequal(expected_keep, observed_keep)

  first_filter[, CellType := cell_type]
  first_filter[, SelectedRank := selected_rank]
  first_filter[, MarkerPass := NumberMarkers > 10]
  first_filter[, AFIPass := PosRatio < 1]
  first_filter[, ExpectedRetained := State %in% expected_keep]
  first_filter[, ObservedRetained := State %in% observed_keep]
  first_filter_rows[[i]] <- first_filter

  first_filter_summary[[i]] <- data.frame(
    CellType = cell_type,
    SelectedRank = selected_rank,
    InitialStates = nrow(first_filter),
    MarkerPassStates = sum(first_filter$MarkerPass),
    MarkerFallback = marker_fallback,
    AFIPassStates = sum(first_filter$AFIPass),
    AFIFallback = afi_fallback,
    ExpectedRetainedStates = length(expected_keep),
    ObservedRetainedStates = length(observed_keep),
    ExactMatch = filter_match
  )
  cat("FIRST_FILTER_AUDITED", cell_type, "expected", length(expected_keep), "observed", length(observed_keep), "match", filter_match, "\n")
}

first_filter_all <- data.table::rbindlist(first_filter_rows, fill = TRUE)
first_filter_audit <- data.table::rbindlist(first_filter_summary)
if (any(!first_filter_audit$ExactMatch)) stop("Marker/AFI retained states do not reproduce for: ", paste(first_filter_audit[ExactMatch == FALSE, CellType], collapse = ", "), call. = FALSE)

# Step3. 按官方 Step 7 公式重算 dropout 指标
dropout_rows <- vector("list", length(rank_cell_type))

for (i in seq_along(rank_cell_type)) {
  cell_type <- rank_cell_type[[i]]
  selected_rank <- rank_selected[[i]]
  state_dir <- file.path(state_discovery_root, cell_type, selected_rank)
  gene_info_file <- file.path(state_dir, "gene_info.txt")
  assignment_file <- file.path(state_dir, "state_assignment.txt")
  hires_file <- file.path(hires_root, paste0(cell_type, ".txt"))
  if (!file.exists(gene_info_file) || !file.exists(assignment_file) || !file.exists(hires_file)) stop("Missing dropout-QC input: ", cell_type, call. = FALSE)

  gene_info <- data.table::fread(gene_info_file, colClasses = "character", check.names = FALSE)
  assignment <- data.table::fread(assignment_file, colClasses = "character", check.names = FALSE)
  hires <- data.table::fread(hires_file, check.names = FALSE)
  if (!all(c("Gene", "State") %in% names(gene_info))) stop("gene_info columns changed: ", cell_type, call. = FALSE)
  if (!all(c("ID", "State") %in% names(assignment))) stop("state_assignment columns changed: ", cell_type, call. = FALSE)

  hires_genes <- as.character(hires[[1L]])
  hires_samples <- names(hires)[-1L]
  gene_index <- match(gene_info$Gene, hires_genes)
  sample_index <- match(assignment$ID, hires_samples)
  if (anyNA(gene_index)) stop("HiRes lacks marker genes: ", cell_type, call. = FALSE)
  if (anyNA(sample_index)) stop("HiRes lacks assigned samples: ", cell_type, call. = FALSE)

  marker_data <- as.matrix(hires[gene_index, -1L, with = FALSE])[, sample_index, drop = FALSE]
  storage.mode(marker_data) <- "double"
  if (any(!is.finite(marker_data)) || any(marker_data < 0)) stop("Invalid HiRes values for dropout QC: ", cell_type, call. = FALSE)

  expression <- apply(log2(marker_data), 1L, mean)
  dropout_gene <- apply(marker_data, 1L, function(x) max(table(x)) / length(x))
  if (any(!is.finite(expression))) stop("Non-finite log2 HiRes marker expression: ", cell_type, call. = FALSE)

  marker_qc <- data.frame(State = gene_info$State, Dropout = dropout_gene, Expression = expression)
  state_dropout <- aggregate(cbind(Dropout, Expression) ~ State, data = marker_qc, FUN = mean)
  state_dropout$CellType <- cell_type
  state_dropout$CellState <- paste0(cell_type, "_", state_dropout$State)
  state_dropout$Markers <- as.integer(table(gene_info$State)[state_dropout$State])
  state_dropout$AssignedSamples <- as.integer(table(assignment$State)[state_dropout$State])
  dropout_rows[[i]] <- state_dropout
  cat("DROPOUT_COMPONENTS", cell_type, nrow(state_dropout), "states", nrow(marker_data), "markers", ncol(marker_data), "assigned samples\n")
}

dropout <- data.table::rbindlist(dropout_rows, fill = TRUE)
dropout[, DropoutZ := as.numeric(scale(Dropout))]
dropout[, ExpressionZ := as.numeric(scale(Expression))]
dropout[, DropoutIndex := as.numeric(scale(DropoutZ + ExpressionZ))]
dropout[, OutlierRaw := DropoutIndex >= 1.96]

# Step4. 应用官方少于两个 states 时不筛选的回退规则，仅估计影响
dropout[, StatesBefore := .N, by = CellType]
dropout[, StatesAfterRaw := sum(!OutlierRaw), by = CellType]
dropout[, FallbackLessThanTwo := StatesAfterRaw < 2L]
dropout[, WouldRemove := OutlierRaw & !FallbackLessThanTwo]
dropout[, WouldRetain := !WouldRemove]
dropout[, IsEcotypeMember := CellState %in% members$ID]

impact <- dropout[, .(
  StatesBefore = unique(StatesBefore),
  RawOutliers = sum(OutlierRaw),
  FallbackLessThanTwo = unique(FallbackLessThanTwo),
  WouldRemove = sum(WouldRemove),
  StatesAfter = sum(WouldRetain),
  EcotypeMembersBefore = sum(IsEcotypeMember),
  EcotypeMembersWouldRemove = sum(IsEcotypeMember & WouldRemove)
), by = CellType]
impact[, WouldChangeStateModel := WouldRemove > 0L]
impact[, WouldChangeEcotypeInput := EcotypeMembersWouldRemove > 0L]
print(impact)

member_impact <- merge(members, dropout[, .(CellState, DropoutIndex, OutlierRaw, FallbackLessThanTwo, WouldRemove)], by.x = "ID", by.y = "CellState", all.x = TRUE, sort = FALSE)
if (anyNA(member_impact$WouldRemove)) stop("Cannot map dropout audit to ecotype members", call. = FALSE)

# Step5. 保存审计结果；不写入正式模型目录
data.table::fwrite(first_filter_all, file.path(audit_dir, "marker_afi_state_audit.tsv"), sep = "\t")
data.table::fwrite(first_filter_audit, file.path(audit_dir, "marker_afi_celltype_summary.tsv"), sep = "\t")
data.table::fwrite(dropout, file.path(audit_dir, "dropout_scores_recalculated.tsv"), sep = "\t")
data.table::fwrite(impact, file.path(audit_dir, "dropout_qc_impact_by_celltype.tsv"), sep = "\t")
data.table::fwrite(member_impact, file.path(audit_dir, "ecotype_member_dropout_impact.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo.txt"))

summary <- list(
  status = "PASS_READ_ONLY_QC_AUDIT",
  marker_afi_exact_match = all(first_filter_audit$ExactMatch),
  states_before_dropout = nrow(dropout),
  raw_dropout_outliers = sum(dropout$OutlierRaw),
  states_would_remove_after_fallback = sum(dropout$WouldRemove),
  ecotype_members_would_remove = sum(dropout$IsEcotypeMember & dropout$WouldRemove),
  cell_types_affected = as.character(impact[WouldChangeStateModel == TRUE, CellType]),
  ecotype_input_cell_types_affected = as.character(impact[WouldChangeEcotypeInput == TRUE, CellType]),
  model_changed = FALSE,
  automatic_rerun = FALSE,
  release_decision = if (sum(dropout$WouldRemove) == 0L) "STEP7_NO_IMPACT" else "REVIEW_IMPACT_BEFORE_ANY_NEW_MODEL"
)
jsonlite::write_json(summary, file.path(audit_dir, "state_qc_impact_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat("STATE_QC_AUDIT_PASS", nrow(dropout), "states;", sum(dropout$WouldRemove), "would be removed; model unchanged\n")
