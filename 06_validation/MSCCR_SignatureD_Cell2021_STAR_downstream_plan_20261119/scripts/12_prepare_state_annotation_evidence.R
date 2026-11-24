#!/usr/bin/env Rscript
# 本脚本整合 85 个正式 cell state 的 discovery marker、四队列 recovery/LODO 和 HiRes 谱系证据。
# 上一轮 53-state 标签仅按同一母细胞类型内的 marker 重叠提供参考，不直接继承状态编号或最终注释。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入并核对正式结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 12_prepare_state_annotation_evidence.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

result_dir <- file.path(output_root, "10_state_annotation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

marker_dir <- file.path(output_root, "09_marker_recovery/ranking_and_LODO")
recovery_dir <- file.path(output_root, "08_scRNA_recovery_summary")
old_annotation_file <- file.path(project_root, "06_validation/paper_aligned_recovery_characterization_20261012/03_state_annotation/53_states_annotation_evidence.tsv")
input_files <- c(
  RankData = rank_data_file,
  Ecotypes = file.path(ecotype_discovery_root, "ecotypes.txt"),
  RecoveryMeta = file.path(recovery_dir, "scRNA_recovery_stouffer_meta.tsv"),
  RecoveryDescriptive = file.path(recovery_dir, "scRNA_recovery_state_across_datasets_descriptive.tsv"),
  TopMarkers = file.path(marker_dir, "state_top_markers.tsv"),
  MarkerLODO = file.path(marker_dir, "marker_LODO_per_state.tsv"),
  HiResContinuous = file.path(scrna_validation_root, "05_posthoc_expression_review/expression_specificity_metrics.tsv"),
  HiResBinary = file.path(scrna_validation_root, "03_results/diagonal_specificity_metrics.tsv"),
  HiResExternalMarkers = file.path(scrna_validation_root, "07_failure_diagnosis/hires_external_marker_recovery.tsv"),
  HiResPartialCorrelation = file.path(scrna_validation_root, "07_failure_diagnosis/hires_own_marker_partial_summary.tsv"),
  PreviousAnnotation = old_annotation_file
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) stop("Missing annotation input: ", paste(names(missing_files), collapse = ", "), call. = FALSE)

manifest <- data.table(Input = names(input_files), Path = unname(input_files))
manifest[, `:=`(Bytes = file.info(Path)$size, SHA256 = vapply(Path, digest::digest, character(1), file = TRUE, algo = "sha256"))]
fwrite(manifest, file.path(result_dir, "state_annotation_input_manifest.tsv"), sep = "\t", quote = FALSE)

# Step2. 汇总正式 discovery marker 和 ecotype 成员关系
rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
discovery_parts <- list()
for (cell_type in names(rank_lookup)) {
  gene_info_file <- file.path(state_discovery_root, cell_type, rank_lookup[[cell_type]], "gene_info.txt")
  if (!file.exists(gene_info_file)) next
  gene_info <- as.data.table(read.delim(gene_info_file, check.names = FALSE, stringsAsFactors = FALSE))
  if (!all(c("Gene", "State", "MaxFC") %in% names(gene_info))) stop("Invalid formal gene_info: ", cell_type, call. = FALSE)
  gene_info <- gene_info[State != "-1" & !is.na(State), .(CellType = cell_type, State, Gene, DiscoveryMaxFC = as.numeric(MaxFC))]
  discovery_parts[[cell_type]] <- gene_info
}
discovery <- rbindlist(discovery_parts)
discovery[, CellState := paste(CellType, State, sep = "_")]
setorder(discovery, CellState, -DiscoveryMaxFC, Gene)
discovery[, DiscoveryOrder := seq_len(.N), by = CellState]
if (uniqueN(discovery$CellState) != 85L) stop("Formal retained state count is not 85", call. = FALSE)

state_summary <- discovery[, .(
  DiscoveryMarkerCount = .N,
  DiscoveryTop20 = paste(Gene[DiscoveryOrder <= 20L], collapse = ";"),
  DiscoveryTop50 = paste(Gene[DiscoveryOrder <= 50L], collapse = ";")
), by = .(CellState, CellType, State)]

ecotypes <- fread(file.path(ecotype_discovery_root, "ecotypes.txt"))
if (!all(c("ID", "Ecotype") %in% names(ecotypes)) || anyDuplicated(ecotypes$ID)) stop("Invalid ecotype membership table", call. = FALSE)
ecotypes <- ecotypes[, .(CellState = ID, Ecotype)]
state_summary <- merge(state_summary, ecotypes, by = "CellState", all.x = TRUE, sort = FALSE)
state_summary[, EcotypeMember := !is.na(Ecotype)]
if (sum(state_summary$EcotypeMember) != 79L || uniqueN(state_summary$Ecotype, na.rm = TRUE) != 8L) stop("Expected 79 states in 8 ecotypes", call. = FALSE)

# Step3. 接入四队列 recovery、marker 排名和留一法结果
recovery_meta <- fread(input_files[["RecoveryMeta"]])
recovery_meta[, CellState := paste(CellType, State, sep = "_")]
recovery_meta <- recovery_meta[, .(
  CellState, RecoveryDatasetsPresent = DatasetsPresent, RecoveryDatasetsEvaluable = DatasetsEvaluable,
  RecoveryDatasetsIndividuallySignificant = DatasetsIndividuallySignificant, RecoveryMetaZ = MetaZ,
  RecoveryMetaP = MetaP_one_sided, RecoveryMetaFDR = MetaFDR_BH, RecoveryMetaSignificant = MetaSignificant_BH
)]

recovery_descriptive <- fread(input_files[["RecoveryDescriptive"]])
recovery_descriptive[, CellState := paste(CellType, State, sep = "_")]
recovery_descriptive <- recovery_descriptive[, .(CellState, RecoverySignificantDatasets = SignificantDatasets, RecoveryAssignedCells = TotalFormalAssignedCells)]

top_markers <- fread(input_files[["TopMarkers"]])
top_markers <- top_markers[HeldOut == "All"]
setorder(top_markers, CellState, Variant, WithinStateOrder, Gene)
top_marker_summary <- top_markers[, .(
  ExternalTopMarkers = paste(Gene, collapse = ";"),
  ExternalTopMarkerN = .N,
  ExternalMetaQSignificantN = sum(is.finite(MetaQ) & MetaQ < 0.05),
  ExternalAllSevenEvaluableN = sum(AllSevenEvaluable == TRUE)
), by = .(CellState, Variant)]
top_marker_wide <- dcast(top_marker_summary, CellState ~ Variant, value.var = c("ExternalTopMarkers", "ExternalTopMarkerN", "ExternalMetaQSignificantN", "ExternalAllSevenEvaluableN"))

lodo <- fread(input_files[["MarkerLODO"]])
lodo_wide <- dcast(lodo, CellState ~ Variant, value.var = c("TotalFolds", "EvaluableFolds", "ConcordantFolds", "Concordance"))

evidence <- merge(state_summary, recovery_meta, by = "CellState", all.x = TRUE, sort = FALSE)
evidence <- merge(evidence, recovery_descriptive, by = "CellState", all.x = TRUE, sort = FALSE)
evidence <- merge(evidence, top_marker_wide, by = "CellState", all.x = TRUE, sort = FALSE)
evidence <- merge(evidence, lodo_wide, by = "CellState", all.x = TRUE, sort = FALSE)

# Step4. 接入 HiRes 母谱系证据；风险等级沿用已完成的表达纯化审计，不重新调阈值
hires_continuous <- fread(input_files[["HiResContinuous"]])
setnames(hires_continuous, c("cell_type", "diagonal_minus_best_offdiagonal", "diagonal_rank", "best_matching_external_type", "auc"), c("CellType", "HiResContinuousMargin", "HiResContinuousRank", "HiResContinuousBestType", "HiResContinuousAUC"))
hires_continuous <- hires_continuous[, .(CellType, HiResContinuousMargin, HiResContinuousRank, HiResContinuousBestType, HiResContinuousAUC)]

hires_binary <- fread(input_files[["HiResBinary"]])
hires_binary <- hires_binary[panel == "COMMON15" & threshold == "any_positive"]
setnames(hires_binary, c("cell_type", "diagonal_jaccard", "diagonal_minus_best_offdiagonal", "diagonal_rank", "best_matching_external_type"), c("CellType", "HiResBinaryJaccard", "HiResBinaryMargin", "HiResBinaryRank", "HiResBinaryBestType"))
hires_binary <- hires_binary[, .(CellType, HiResBinaryJaccard, HiResBinaryMargin, HiResBinaryRank, HiResBinaryBestType)]

hires_external <- fread(input_files[["HiResExternalMarkers"]])
setnames(hires_external, c("hires_cell_type", "markers_in_hires", "marker_recovery_from_signature", "fisher_odds_ratio", "fisher_p_value"), c("CellType", "HiResExternalMarkersRecovered", "HiResExternalMarkerRecovery", "HiResExternalMarkerOR", "HiResExternalMarkerP"))
hires_external <- hires_external[, .(CellType, HiResExternalMarkersRecovered, HiResExternalMarkerRecovery, HiResExternalMarkerOR, HiResExternalMarkerP)]

hires_partial <- fread(input_files[["HiResPartialCorrelation"]])
setnames(hires_partial, c("hires_cell_type", "own_fraction_adjusted_spearman", "strongest_other_fraction", "strongest_other_partial_spearman"), c("CellType", "HiResOwnAdjustedRho", "HiResStrongestOtherFraction", "HiResStrongestOtherAdjustedRho"))
hires_partial <- hires_partial[, .(CellType, HiResOwnAdjustedRho, HiResStrongestOtherFraction, HiResStrongestOtherAdjustedRho)]

evidence <- merge(evidence, hires_continuous, by = "CellType", all.x = TRUE, sort = FALSE)
evidence <- merge(evidence, hires_binary, by = "CellType", all.x = TRUE, sort = FALSE)
evidence <- merge(evidence, hires_external, by = "CellType", all.x = TRUE, sort = FALSE)
evidence <- merge(evidence, hires_partial, by = "CellType", all.x = TRUE, sort = FALSE)

lineage_audit <- data.table(
  CellType = unique(evidence$CellType),
  HiResLineageTier = "SUPPORTED_WITH_CAVEAT",
  HiResLineageCaveat = "母谱系可作条件性解释；cell state 注释仍需 marker 与外部恢复共同支持"
)
lineage_audit[CellType == "ILC", `:=`(HiResLineageTier = "HIGH_RISK", HiResLineageCaveat = "HiRes 母谱系特异性审计失败；仅解释为 ILC 输入槽中的表达程序")]
lineage_audit[CellType %in% c("DC", "Monocyte_Macrophage"), `:=`(HiResLineageTier = "MYELOID_CROSSLINEAGE", HiResLineageCaveat = "DC 与单核/巨噬细胞程序互相渗漏；优先使用髓系程序而非严格亚型名称")]
lineage_audit[CellType %in% c("Enteroendocrine", "Tuft", "Goblet", "Stem"), `:=`(HiResLineageTier = "EPITHELIAL_SHARED_PROGRAM", HiResLineageCaveat = "上皮细胞间共享分化程序明显；优先使用功能程序名称")]
lineage_audit[CellType == "Endothelial", `:=`(HiResLineageTier = "WEAK_LINEAGE_SUPPORT", HiResLineageCaveat = "HiRes 母谱系证据较弱；注释需保守")]
lineage_audit[CellType == "CD8_T", HiResLineageCaveat := "原始 HiRes 矩阵无 NA；EcoTyper 缩放补丁仅将 11 个零方差基因的非有限值置零，未进入最终 gene_info"]
evidence <- merge(evidence, lineage_audit, by = "CellType", all.x = TRUE, sort = FALSE)

# Step5. 仅按同一母细胞类型内 marker 重叠寻找上一轮参考状态
old_annotation <- fread(old_annotation_file, colClasses = "character")
old_annotation[, OldGenes := strsplit(Top20Genes, ",", fixed = TRUE)]
overlap_parts <- list()
for (state_index in seq_len(nrow(state_summary))) {
  current_state <- state_summary$CellState[[state_index]]
  current_type <- state_summary$CellType[[state_index]]
  current_genes <- strsplit(state_summary$DiscoveryTop20[[state_index]], ";", fixed = TRUE)[[1L]]
  candidates <- old_annotation[CellType == current_type]
  if (!nrow(candidates)) {
    overlap_parts[[state_index]] <- data.table(CellState = current_state, PreviousCellState = NA_character_, PreviousProgram = NA_character_, PreviousJudgment = NA_character_, PreviousOverlapN = 0L, PreviousJaccard = 0, PreviousCurrentCoverage = 0, PreviousOldCoverage = 0)
    next
  }
  candidate_scores <- candidates[, {
    old_genes <- OldGenes[[1L]]
    overlap <- intersect(current_genes, old_genes)
    list(PreviousOverlapN = length(overlap), PreviousJaccard = length(overlap) / length(union(current_genes, old_genes)), PreviousCurrentCoverage = length(overlap) / length(current_genes), PreviousOldCoverage = length(overlap) / length(old_genes), PreviousOverlapGenes = paste(overlap, collapse = ";"))
  }, by = .(PreviousCellState = CellState, PreviousProgram = InitialProgram, PreviousJudgment = Judgment)]
  setorder(candidate_scores, -PreviousOverlapN, -PreviousJaccard, PreviousCellState)
  overlap_parts[[state_index]] <- cbind(data.table(CellState = current_state), candidate_scores[1L])
}
previous_overlap <- rbindlist(overlap_parts, fill = TRUE)
previous_overlap[, PreviousMatchStrength := fifelse(PreviousOverlapN >= 5L | PreviousJaccard >= 0.20, "STRONG", fifelse(PreviousOverlapN >= 3L, "MODERATE", fifelse(PreviousOverlapN >= 1L, "WEAK", "NONE")))]
fwrite(previous_overlap, file.path(result_dir, "state_annotation_previous_model_overlap.tsv"), sep = "\t", quote = FALSE)
evidence <- merge(evidence, previous_overlap, by = "CellState", all.x = TRUE, sort = FALSE)

# Step6. 分开评定 recovery、LODO 和母谱系证据，不把单一总分当作验证结论
evidence[, RecoveryEvidence := fifelse(RecoveryMetaSignificant == TRUE & RecoveryDatasetsIndividuallySignificant >= 2L, "HIGH", fifelse(RecoveryMetaSignificant == TRUE, "MODERATE", "LOW"))]
evidence[, LODODiseaseEvidence := fifelse(EvaluableFolds_Disease >= 2L & Concordance_Disease >= 0.75, "HIGH", fifelse(EvaluableFolds_Disease >= 1L & Concordance_Disease >= 0.50, "MODERATE", fifelse(EvaluableFolds_Disease >= 1L, "LOW", "NOT_EVALUABLE")))]
evidence[, LODORegion3Evidence := fifelse(EvaluableFolds_Region3 >= 2L & Concordance_Region3 >= 0.75, "HIGH", fifelse(EvaluableFolds_Region3 >= 1L & Concordance_Region3 >= 0.50, "MODERATE", fifelse(EvaluableFolds_Region3 >= 1L, "LOW", "NOT_EVALUABLE")))]
evidence[, StateEvidenceTier := "EXPLORATORY"]
evidence[RecoveryEvidence == "HIGH" & (LODODiseaseEvidence == "HIGH" | LODORegion3Evidence == "HIGH") & !HiResLineageTier %in% c("HIGH_RISK", "WEAK_LINEAGE_SUPPORT"), StateEvidenceTier := "HIGH"]
evidence[RecoveryMetaSignificant == TRUE & StateEvidenceTier == "EXPLORATORY" & HiResLineageTier != "HIGH_RISK", StateEvidenceTier := "MODERATE"]
evidence[HiResLineageTier %in% c("HIGH_RISK", "MYELOID_CROSSLINEAGE", "EPITHELIAL_SHARED_PROGRAM", "WEAK_LINEAGE_SUPPORT"), AnnotationNamingRule := "PROGRAM_LABEL"]
evidence[is.na(AnnotationNamingRule), AnnotationNamingRule := "LINEAGE_CONDITIONAL_PROGRAM"]
evidence[, PreviousProgramCandidate := fifelse(PreviousMatchStrength %in% c("STRONG", "MODERATE"), PreviousProgram, NA_character_)]
setorder(evidence, CellType, State)

review_template <- evidence[, .(
  CellState, CellType, State, Ecotype, EcotypeMember, DiscoveryTop20,
  ExternalTopMarkers_Disease, ExternalTopMarkers_Region3,
  PreviousCellState, PreviousProgramCandidate, PreviousMatchStrength, PreviousOverlapGenes,
  RecoveryEvidence, LODODiseaseEvidence, LODORegion3Evidence, HiResLineageTier, HiResLineageCaveat,
  StateEvidenceTier, AnnotationNamingRule,
  CandidateProgram = PreviousProgramCandidate, FinalStateLabel = NA_character_, AnnotationConfidence = NA_character_,
  MarkerRationale = NA_character_, LiteratureEvidence = NA_character_, ReviewerDecision = "PENDING"
)]

fwrite(evidence, file.path(result_dir, "state_annotation_evidence.tsv"), sep = "\t", quote = FALSE)
fwrite(review_template, file.path(result_dir, "state_annotation_review_template.tsv"), sep = "\t", quote = FALSE)

summary_table <- rbindlist(list(
  evidence[, .(Section = "StateEvidenceTier", Category = StateEvidenceTier, N = .N), by = StateEvidenceTier][, StateEvidenceTier := NULL],
  evidence[, .(Section = "HiResLineageTier", Category = HiResLineageTier, N = .N), by = HiResLineageTier][, HiResLineageTier := NULL],
  evidence[, .(Section = "PreviousMatchStrength", Category = PreviousMatchStrength, N = .N), by = PreviousMatchStrength][, PreviousMatchStrength := NULL],
  evidence[, .(Section = "EcotypeMember", Category = as.character(EcotypeMember), N = .N), by = EcotypeMember][, EcotypeMember := NULL]
), use.names = TRUE)
fwrite(summary_table, file.path(result_dir, "state_annotation_evidence_summary.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_STATE_ANNOTATION_EVIDENCE_ASSEMBLY",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  states = nrow(evidence), ecotype_members = sum(evidence$EcotypeMember), ecotypes = uniqueN(evidence$Ecotype, na.rm = TRUE),
  high_state_evidence = sum(evidence$StateEvidenceTier == "HIGH"), moderate_state_evidence = sum(evidence$StateEvidenceTier == "MODERATE"), exploratory_state_evidence = sum(evidence$StateEvidenceTier == "EXPLORATORY"),
  previous_strong_or_moderate_matches = sum(evidence$PreviousMatchStrength %in% c("STRONG", "MODERATE"))
)
jsonlite::write_json(completion, file.path(result_dir, "state_annotation_evidence_complete.json"), pretty = TRUE, auto_unbox = TRUE)

print(dim(evidence))
print(table(evidence$StateEvidenceTier, useNA = "ifany"))
print(table(evidence$HiResLineageTier, useNA = "ifany"))
print(table(evidence$PreviousMatchStrength, useNA = "ifany"))
cat("PASS_STATE_ANNOTATION_EVIDENCE_ASSEMBLY\n")
