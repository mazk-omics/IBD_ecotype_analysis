#!/usr/bin/env Rscript
# 本脚本汇总四队列 Seurat 3.1.3 检验，计算 Cell 2021 n1-n7 和 top-10 markers。
# 每次留出完整数据集后重新做基因归属、七指标排名和预测，禁止 held-out 信息泄漏。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 核对四队列检查点并读取正式 discovery gene assignment
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 11_rank_markers_n1_n7_and_LODO.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite", "Matrix")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(Matrix)

datasets <- c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999")
full_count_root <- file.path(output_root, "09_marker_recovery/full_counts")
de_root <- file.path(output_root, "09_marker_recovery/DE_checkpoints")
result_dir <- file.path(output_root, "09_marker_recovery/ranking_and_LODO")
completed_file <- file.path(result_dir, "ranking_LODO_complete.json")
if (file.exists(completed_file)) stop("Marker ranking and LODO already complete; do not overwrite", call. = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
discovery_parts <- list()
for (cell_type in names(rank_lookup)) {
  gene_info_file <- file.path(state_discovery_root, cell_type, rank_lookup[[cell_type]], "gene_info.txt")
  if (!file.exists(gene_info_file)) next
  gene_info <- as.data.table(read.delim(gene_info_file, check.names = FALSE, stringsAsFactors = FALSE))
  if (!all(c("Gene", "State", "MaxFC") %in% names(gene_info))) stop("Invalid formal gene_info: ", cell_type, call. = FALSE)
  gene_info <- gene_info[State != "-1" & !is.na(State), .(CellType = cell_type, Gene, State, DiscoveryMaxFC = as.numeric(MaxFC))]
  discovery_parts[[cell_type]] <- gene_info
}
discovery <- rbindlist(discovery_parts)
if (anyDuplicated(discovery[, .(CellType, Gene)])) stop("Formal discovery gene assignment is duplicated", call. = FALSE)
discovery[, CellState := paste(CellType, State, sep = "_")]
frozen_states <- unique(discovery[, .(CellState, CellType, State)])
setorder(frozen_states, CellType, State)
if (nrow(frozen_states) != 85L) stop("Formal retained state count is not 85", call. = FALSE)

de_parts <- list()
fold_parts <- list()
expression_parts <- list()
eligibility_parts <- list()
checkpoint_manifest_parts <- list()

for (dataset in datasets) {
  manifest_file <- file.path(full_count_root, dataset, paste0(dataset, "_marker_counts_manifest.tsv"))
  if (!file.exists(manifest_file) && dataset == "IBDverse_E_MTAB_16999") manifest_file <- file.path(full_count_root, dataset, "IBDverse_E_MTAB_16999_marker_counts_manifest.tsv")
  if (!file.exists(manifest_file)) stop("Missing full-count manifest: ", dataset, call. = FALSE)
  manifest <- fread(manifest_file, colClasses = "character")
  checkpoint_manifest_parts[[dataset]] <- manifest

  for (cell_type in manifest$CellType) {
    checkpoint_dir <- file.path(de_root, dataset, cell_type)
    complete <- file.path(checkpoint_dir, "completed.json")
    if (!file.exists(complete)) stop("Missing marker DE checkpoint: ", dataset, " ", cell_type, call. = FALSE)
    complete_info <- jsonlite::read_json(complete, simplifyVector = TRUE)
    if (!identical(complete_info$status, "PASS_MARKER_DE_SEURAT313")) stop("Marker DE checkpoint did not pass: ", dataset, " ", cell_type, call. = FALSE)

    expression <- fread(file.path(checkpoint_dir, "gene_expression.tsv.gz"))
    expression_parts[[paste(dataset, cell_type, sep = "__")]] <- expression
    eligibility <- fread(file.path(checkpoint_dir, "state_eligibility.tsv"))
    eligibility_parts[[paste(dataset, cell_type, sep = "__")]] <- eligibility

    de_files <- list.files(checkpoint_dir, pattern = "__DEG\\.tsv\\.gz$", full.names = TRUE)
    fc_files <- list.files(checkpoint_dir, pattern = "__candidate_full_FC\\.tsv\\.gz$", full.names = TRUE)
    if (length(de_files)) for (path in de_files) de_parts[[length(de_parts) + 1L]] <- fread(path)
    if (length(fc_files)) for (path in fc_files) fold_parts[[length(fold_parts) + 1L]] <- fread(path)
  }
}

de <- rbindlist(de_parts, fill = TRUE)
folds <- rbindlist(fold_parts, fill = TRUE)
expression <- rbindlist(expression_parts, fill = TRUE)
eligibility <- rbindlist(eligibility_parts, fill = TRUE)
if (!nrow(de) || !nrow(folds) || !nrow(expression) || !nrow(eligibility)) stop("Marker checkpoint aggregation is empty", call. = FALSE)
de_keys <- c("Dataset", "CellType", "Scope", "Stratum", "Gene", "State")
if (anyDuplicated(de[, ..de_keys])) stop("Aggregated DE keys are duplicated", call. = FALSE)
if (anyDuplicated(folds[, ..de_keys])) stop("Aggregated candidate full-FC keys are duplicated", call. = FALSE)
if (any(!is.finite(de$p_val) | de$p_val < 0 | de$p_val > 1) || any(!is.finite(de$log2FC))) stop("Aggregated DE contains invalid values", call. = FALSE)
if (any(!is.finite(folds$log2FC))) stop("Aggregated candidate full-FC table contains invalid values", call. = FALSE)
fwrite(de, file.path(result_dir, "external_DEG_by_dataset.tsv.gz"), sep = "\t", quote = FALSE)

# Step2. 固定每个数据集/分层中的外部 gene-state assignment
positive <- folds[log2FC > 0]
positive[, MaxLog2FC := max(log2FC), by = .(Dataset, CellType, Scope, Stratum, Gene)]
assigned <- positive[abs(log2FC - MaxLog2FC) <= 1e-12]
assigned[, AssignmentStates := uniqueN(State), by = .(Dataset, CellType, Scope, Stratum, Gene)]
assigned[, UniqueDatasetAssignment := AssignmentStates == 1L]
fwrite(assigned, file.path(result_dir, "external_gene_state_assignment.tsv.gz"), sep = "\t", quote = FALSE)

training_sets <- list(All = datasets)
for (held_out in datasets) training_sets[[held_out]] <- setdiff(datasets, held_out)
metrics_parts <- list()
ranking_parts <- list()
exclusion_parts <- list()

# Step3. 全数据和每个留出折分别重建 n1-n7
for (held_out in names(training_sets)) {
  training <- training_sets[[held_out]]
  train_de <- de[Dataset %in% training]
  train_folds <- folds[Dataset %in% training]
  train_assign <- assigned[Dataset %in% training]
  train_expression <- expression[Dataset %in% training]

  overall_fc <- train_folds[Scope == "Dataset" & log2FC > 0]
  overall_fc[, GlobalMaxLog2FC := max(log2FC), by = .(CellType, Gene)]
  winners <- overall_fc[abs(log2FC - GlobalMaxLog2FC) <= 1e-12]
  means <- train_folds[Scope == "Dataset", .(MeanAcrossEvaluableDatasets = mean(log2FC)), by = .(CellType, Gene, State)]
  winners <- merge(winners, means, by = c("CellType", "Gene", "State"), all.x = TRUE, sort = FALSE)
  winners[, WinnerMean := max(MeanAcrossEvaluableDatasets), by = .(CellType, Gene)]
  winners <- unique(winners[abs(MeanAcrossEvaluableDatasets - WinnerMean) <= 1e-12, .(CellType, Gene, State)])
  unresolved <- winners[, .N, by = .(CellType, Gene)][N > 1L]
  winners <- winners[!unresolved, on = .(CellType, Gene)]

  candidates <- merge(winners, discovery[, .(CellType, Gene, State_Discovery = State)], by = c("CellType", "Gene"), all = FALSE, sort = FALSE)
  setnames(candidates, "State", "State_External")
  consistent <- candidates[State_External == State_Discovery]
  consistent[, `:=`(State = State_External, CellState = paste(CellType, State_External, sep = "_"))]
  exclusion_parts[[held_out]] <- data.table(
    HeldOut = held_out, TrainingDatasets = paste(training, collapse = ";"), UnresolvedGeneTies = nrow(unresolved), DiscoveryMatchedGenes = nrow(candidates),
    ConsistentCandidates = nrow(consistent), DiscordantExcluded = nrow(candidates) - nrow(consistent)
  )
  if (!nrow(consistent)) stop("No consistent marker candidates in fold: ", held_out, call. = FALSE)

  n1_table <- train_expression[MeanCPM > 0, .(n1 = uniqueN(Dataset)), by = .(CellType, Gene)]
  n2_table <- train_assign[Scope == "Dataset" & UniqueDatasetAssignment == TRUE, .(n2 = uniqueN(Dataset)), by = .(CellType, Gene, State)]
  specific_de <- merge(train_de[Scope == "Dataset"], consistent[, .(CellType, Gene, State)], by = c("CellType", "Gene", "State"), all = FALSE, sort = FALSE)
  n5_meta <- specific_de[, {
    z <- qnorm(pmax(p_val, 1e-300) / 2, lower.tail = FALSE) * sign(log2FC)
    meta_z <- sum(z) / sqrt(.N)
    list(n5 = sum(p_val_adj < 0.05), MetaZUnadjusted = meta_z, MetaP = 2 * pnorm(-abs(meta_z)), TestedDatasets = uniqueN(Dataset))
  }, by = .(CellType, Gene, State)]

  strata_assignment <- train_assign[UniqueDatasetAssignment == TRUE & Scope %in% c("Disease", "Region3"), .(n4 = uniqueN(Stratum)), by = .(CellType, Gene, State, Scope)]
  strata_evaluable <- eligibility[Dataset %in% training & EligibleGE5 == TRUE & EvaluableContrast == TRUE & Scope %in% c("Disease", "Region3"), .(EvaluableStrata = uniqueN(Stratum)), by = .(CellType, State, Scope)]

  metrics <- consistent[, .(HeldOut = held_out, TrainingDatasets = paste(training, collapse = ";"), CellState, CellType, State, Gene)]
  metrics <- merge(metrics, n1_table, by = c("CellType", "Gene"), all.x = TRUE, sort = FALSE)
  metrics <- merge(metrics, n2_table, by = c("CellType", "Gene", "State"), all.x = TRUE, sort = FALSE)
  metrics <- merge(metrics, n5_meta, by = c("CellType", "Gene", "State"), all.x = TRUE, sort = FALSE)
  metrics <- merge(metrics, means, by = c("CellType", "Gene", "State"), all.x = TRUE, sort = FALSE)
  setnames(metrics, "MeanAcrossEvaluableDatasets", "n7")
  metrics[is.na(n1), n1 := 0L]
  metrics[is.na(n2), n2 := 0L]
  metrics[is.na(n5), n5 := 0L]
  metrics[, n3 := fifelse(n1 > 0, n2 / n1, NA_real_)]

  for (scope in c("Disease", "Region3")) {
    scope_assignment <- strata_assignment[Scope == scope, .(CellType, Gene, State, n4)]
    scope_evaluable <- strata_evaluable[Scope == scope, .(CellType, State, EvaluableStrata)]
    metrics <- merge(metrics, scope_assignment, by = c("CellType", "Gene", "State"), all.x = TRUE, sort = FALSE)
    metrics <- merge(metrics, scope_evaluable, by = c("CellType", "State"), all.x = TRUE, sort = FALSE)
    metrics[!is.na(EvaluableStrata) & EvaluableStrata > 0 & is.na(n4), n4 := 0L]
    setnames(metrics, c("n4", "EvaluableStrata"), c(paste0("n4_", scope), paste0("Evaluable_", scope)))
  }

  metrics[, MetaQ := NA_real_]
  metrics[is.finite(MetaP), MetaQ := p.adjust(MetaP, method = "BH"), by = CellState]
  metrics[, n6 := qnorm(pmax(MetaQ, 1e-300) / 2, lower.tail = FALSE) * sign(MetaZUnadjusted)]
  metrics_parts[[held_out]] <- metrics

  for (variant in c("Disease", "Region3")) {
    ranking <- copy(metrics)
    ranking[, n4 := get(paste0("n4_", variant))]
    metric_columns <- c("n1", "n2", "n3", "n4", "n5", "n6", "n7")
    ranking[, AllSevenEvaluable := rowSums(is.na(.SD)) == 0L, .SDcols = metric_columns]
    ranking <- ranking[AllSevenEvaluable == TRUE]
    if (!nrow(ranking)) next
    for (metric in metric_columns) ranking[, paste0(metric, "_rank") := rank(get(metric), ties.method = "average"), by = CellState]
    rank_columns <- paste0(metric_columns, "_rank")
    ranking[, CompositeScore := rowMeans(.SD), .SDcols = rank_columns]
    ranking[, Variant := variant]
    setorder(ranking, CellState, -CompositeScore, Gene)
    ranking[, WithinStateOrder := seq_len(.N), by = CellState]
    ranking_parts[[paste(held_out, variant, sep = "__")]] <- ranking
  }
  cat("MARKER_RANK_FOLD", held_out, length(training), "training datasets", nrow(consistent), "consistent genes\n")
}

metrics_all <- rbindlist(metrics_parts, fill = TRUE)
rankings <- rbindlist(ranking_parts, fill = TRUE)
exclusions <- rbindlist(exclusion_parts)
top_markers <- rankings[WithinStateOrder <= 10L]
fwrite(metrics_all, file.path(result_dir, "marker_n1_n7.tsv.gz"), sep = "\t", quote = FALSE)
fwrite(rankings, file.path(result_dir, "marker_composite_ranking.tsv.gz"), sep = "\t", quote = FALSE)
fwrite(top_markers, file.path(result_dir, "state_top_markers.tsv"), sep = "\t", quote = FALSE)
fwrite(exclusions, file.path(result_dir, "gene_assignment_exclusions.tsv"), sep = "\t", quote = FALSE)

# Step4. 比较 Disease 与 Region3 两版 n4
comparison_parts <- list()
comparison_groups <- unique(rankings[, .(HeldOut, CellState)])
for (group_index in seq_len(nrow(comparison_groups))) {
  held_out <- comparison_groups$HeldOut[[group_index]]
  cell_state <- comparison_groups$CellState[[group_index]]
  disease <- rankings[HeldOut == held_out & CellState == cell_state & Variant == "Disease"]
  region <- rankings[HeldOut == held_out & CellState == cell_state & Variant == "Region3"]
  common <- intersect(disease$Gene, region$Gene)
  top_disease <- disease[WithinStateOrder <= 10L]$Gene
  top_region <- region[WithinStateOrder <= 10L]$Gene
  rho <- if (length(common) >= 3L) cor(disease[match(common, Gene)]$CompositeScore, region[match(common, Gene)]$CompositeScore, method = "spearman") else NA_real_
  union_top <- union(top_disease, top_region)
  comparison_parts[[group_index]] <- data.table(
    HeldOut = held_out, CellState = cell_state, CommonRankableGenes = length(common), RankSpearman = rho, DiseaseTopN = length(top_disease), RegionTopN = length(top_region),
    Top10Intersection = length(intersect(top_disease, top_region)), Top10Jaccard = if (length(union_top)) length(intersect(top_disease, top_region)) / length(union_top) else NA_real_,
    DiseaseOnlyTop = paste(setdiff(top_disease, top_region), collapse = ";"), RegionOnlyTop = paste(setdiff(top_region, top_disease), collapse = ";")
  )
}
comparison <- rbindlist(comparison_parts)
fwrite(comparison, file.path(result_dir, "n4_disease_vs_region_comparison.tsv"), sep = "\t", quote = FALSE)

# Step5. 整数据集留一法评价 top-10 markers
lodo_parts <- list()
lodo_expression_parts <- list()
for (held_out in datasets) {
  train_ranking <- rankings[HeldOut == held_out & WithinStateOrder <= 10L]
  manifest <- checkpoint_manifest_parts[[held_out]]

  for (cell_type in manifest$CellType) {
    type_dir <- file.path(full_count_root, held_out, cell_type)
    raw <- readRDS(file.path(type_dir, "counts.rds"))
    cells <- fread(file.path(type_dir, "cells.tsv.gz"), colClasses = "character")
    cells[, Assigned := Assigned == "TRUE"]
    if (!methods::is(raw, "dgCMatrix") || !identical(colnames(raw), cells$CellID)) stop("LODO count/cell identity mismatch: ", held_out, " ", cell_type, call. = FALSE)
    totals <- Matrix::colSums(raw)
    state_counts <- cells[Assigned == TRUE, .N, by = State]
    allowed <- state_counts[N >= 5L]$State
    candidates <- train_ranking[CellType == cell_type]

    candidate_groups <- unique(candidates[, .(Variant, State)])
    for (candidate_index in seq_len(nrow(candidate_groups))) {
      variant <- candidate_groups$Variant[[candidate_index]]
      state <- candidate_groups$State[[candidate_index]]
      target <- paste(cell_type, state, sep = "_")
      markers <- candidates[Variant == variant & State == state][order(WithinStateOrder)]
      present <- markers$Gene[markers$Gene %in% rownames(raw)]

      if (nrow(markers) != 10L || length(present) != 10L || !state %in% allowed || length(allowed) < 2L) {
        lodo_parts[[length(lodo_parts) + 1L]] <- data.table(
          HeldOut = held_out, Variant = variant, CellState = target, Status = "NOT_EVALUABLE",
          Reason = "Need 10 trained/present markers, >=5 assigned target cells and >=2 eligible held-out states", TrainedMarkers = nrow(markers), PresentMarkers = length(present), Concordant = NA_real_
        )
        next
      }

      values <- as.matrix(raw[present, , drop = FALSE])
      values <- log2(sweep(values, 2, totals, "/") * 1e6 + 1)
      marker_sd <- apply(values, 1, sd)
      if (any(!is.finite(marker_sd) | marker_sd == 0)) {
        lodo_parts[[length(lodo_parts) + 1L]] <- data.table(
          HeldOut = held_out, Variant = variant, CellState = target, Status = "NOT_EVALUABLE", Reason = "Zero marker variance in held-out parent cells",
          TrainedMarkers = 10L, PresentMarkers = 10L, Concordant = NA_real_
        )
        next
      }
      values <- sweep(values, 1, rowMeans(values), "-")
      values <- sweep(values, 1, marker_sd, "/")

      scores <- numeric(length(allowed))
      names(scores) <- allowed
      for (other in allowed) scores[[other]] <- mean(values[, cells$State == other, drop = FALSE])
      winner <- names(scores)[abs(scores - max(scores)) <= 1e-12]
      lodo_expression_parts[[length(lodo_expression_parts) + 1L]] <- data.table(
        HeldOut = held_out, Variant = variant, MarkerState = target, ExpressionState = paste(cell_type, allowed, sep = "_"), MeanUnitVarianceExpression = as.numeric(scores)
      )
      lodo_parts[[length(lodo_parts) + 1L]] <- data.table(
        HeldOut = held_out, Variant = variant, CellState = target, Status = if (length(winner) == 1L) "EVALUABLE" else "TIE", Reason = "",
        PredictedState = paste(paste(cell_type, winner, sep = "_"), collapse = ";"), Concordant = as.numeric(length(winner) == 1L && winner == state),
        TrainedMarkers = 10L, PresentMarkers = 10L, TopScore = max(scores), TargetScore = scores[[state]]
      )
    }
    rm(raw)
    invisible(gc(FALSE))
  }
  cat("MARKER_LODO_HELD_OUT", held_out, "complete\n")
}

lodo_observed <- rbindlist(lodo_parts, fill = TRUE)
full_grid <- CJ(HeldOut = datasets, Variant = c("Disease", "Region3"), CellState = frozen_states$CellState, unique = TRUE)
lodo <- merge(full_grid, lodo_observed, by = c("HeldOut", "Variant", "CellState"), all.x = TRUE, sort = FALSE)
lodo[is.na(Status), `:=`(Status = "NOT_EVALUABLE", Reason = "No trained seven-metric ranking or no input parent lineage")]
lodo[is.na(Reason), Reason := ""]
fwrite(lodo, file.path(result_dir, "marker_LODO_per_dataset.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(lodo_expression_parts, fill = TRUE), file.path(result_dir, "marker_LODO_mean_expression.tsv.gz"), sep = "\t", quote = FALSE)

state_summary <- lodo[, .(TotalFolds = .N, EvaluableFolds = sum(Status == "EVALUABLE"), ConcordantFolds = sum(Concordant, na.rm = TRUE)), by = .(CellState, Variant)]
state_summary[, Concordance := fifelse(EvaluableFolds > 0, ConcordantFolds / EvaluableFolds, NA_real_)]
summary <- lodo[Status == "EVALUABLE", .(EvaluableStates = .N, ConcordantStates = sum(Concordant), Concordance = mean(Concordant)), by = .(HeldOut, Variant)]
fwrite(state_summary, file.path(result_dir, "marker_LODO_per_state.tsv"), sep = "\t", quote = FALSE)
fwrite(summary, file.path(result_dir, "marker_LODO_summary.tsv"), sep = "\t", quote = FALSE)

# Step6. 保存审计身份
final_summary <- list(
  status = "PASS_MARKER_N1_N7_AND_LODO", datasets = datasets, formal_states = nrow(frozen_states), metric_rows = nrow(metrics_all), ranking_rows = nrow(rankings),
  top_marker_rows = nrow(top_markers), no_held_out_DE_used_for_ranking = TRUE, lodo_unit = "whole scRNA-seq dataset",
  n4_variants = c("Disease", "Region3"), disease_n4_method_identity = "IBD_EXTENSION", region3_n4_method_identity = "SUPPLEMENTARY_SENSITIVITY",
  marker_normalization = "held-out raw counts -> log2(CPM+1) -> per-gene unit variance within all parent-type cells",
  lodo_summary = summary
)
jsonlite::write_json(final_summary, completed_file, pretty = TRUE, auto_unbox = TRUE, dataframe = "rows")
writeLines(capture.output(sessionInfo()), file.path(result_dir, "R_sessionInfo.txt"))
print(summary)
cat("MARKER_N1_N7_AND_LODO_PASS", nrow(metrics_all), "metric rows", nrow(top_markers), "top-marker rows\n")
