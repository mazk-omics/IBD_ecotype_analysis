#!/usr/bin/env Rscript
# 本脚本将冻结 marker、通路和 validation 证据整理成一张 85-state 注释总表。
# 注释采用宽泛功能程序；不重新筛选 state，不修改 EcoTyper 模型。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入和输出
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 15_finalize_state_annotation_table.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

annotation_dir <- file.path(output_root, "10_state_annotation")
evidence_file <- file.path(annotation_dir, "state_annotation_evidence.tsv")
pathway_file <- file.path(annotation_dir, "pathway_enrichment/state_marker_pathway_top_for_review.tsv")
label_file <- file.path(RUN_ROOT, "inputs/state_annotation_labels_first_pass.tsv")
if (any(!file.exists(c(evidence_file, pathway_file, label_file)))) stop("State annotation input is missing", call. = FALSE)

# Step2. 读取冻结证据并核对 85 个 state
evidence <- fread(evidence_file, colClasses = "character", check.names = FALSE)
labels <- fread(label_file, colClasses = "character", check.names = FALSE)
pathways <- fread(pathway_file, colClasses = "character", check.names = FALSE)

if (nrow(evidence) != 85L || anyDuplicated(evidence$CellState)) stop("Frozen state evidence identity failed", call. = FALSE)
if (nrow(labels) != 85L || anyDuplicated(labels$CellState) || !setequal(labels$CellState, evidence$CellState)) stop("State label identity failed", call. = FALSE)
if (any(!labels$AnnotationConfidence %in% c("HIGH", "MODERATE", "LOW"))) stop("Unexpected annotation confidence", call. = FALSE)

pathways[, Significant := Significant == "TRUE"]
pathways[, RankWithinCollection := as.integer(RankWithinCollection)]
top_pathways <- pathways[Significant == TRUE][order(CellState, Collection, RankWithinCollection), .SD[1L], by = .(CellState, Collection)]
top_pathways <- dcast(top_pathways, CellState ~ Collection, value.var = "GeneSet")
for (column_name in c("GO_BP", "HALLMARK", "REACTOME")) if (!column_name %in% names(top_pathways)) top_pathways[, (column_name) := NA_character_]

# Step3. 合并标签、marker、通路和 validation
state_table <- merge(evidence, labels, by = "CellState", all.x = TRUE, sort = FALSE)
state_table <- merge(state_table, top_pathways[, .(CellState, GO_BP, HALLMARK, REACTOME)], by = "CellState", all.x = TRUE, sort = FALSE)
if (anyNA(state_table$FinalStateLabel) || nrow(state_table) != 85L) stop("Final annotation merge failed", call. = FALSE)

state_table[, KeyMarkers := vapply(strsplit(DiscoveryTop20, ";", fixed = TRUE), function(marker_set) paste(head(marker_set, 10L), collapse = ";"), character(1L))]
state_table[, ExternalValidatedMarkers := mapply(function(disease_markers, region_markers) {
  marker_set <- unique(c(head(strsplit(disease_markers, ";", fixed = TRUE)[[1L]], 5L), head(strsplit(region_markers, ";", fixed = TRUE)[[1L]], 5L)))
  paste(marker_set[nzchar(marker_set)], collapse = ";")
}, ExternalTopMarkers_Disease, ExternalTopMarkers_Region3, USE.NAMES = FALSE)]

state_table[, KeyPathways := apply(.SD, 1L, function(pathway_set) {
  pathway_set <- pathway_set[!is.na(pathway_set) & nzchar(pathway_set)]
  pathway_set <- gsub("^(HALLMARK_|GOBP_|REACTOME_)", "", pathway_set)
  paste(pathway_set, collapse = " | ")
}), .SDcols = c("HALLMARK", "GO_BP", "REACTOME")]

state_table[, ValidationSummary := paste0(
  "Recovery=", RecoveryEvidence,
  "; LODO-disease=", LODODiseaseEvidence,
  "; LODO-region=", LODORegion3Evidence,
  "; HiRes=", HiResLineageTier
)]

final_table <- state_table[, .(
  CellState, CellType, State, Ecotype, EcotypeMember,
  FinalStateLabel, FinalStateLabel_CN, AnnotationConfidence, StateEvidenceTier,
  KeyMarkers, ExternalValidatedMarkers, KeyPathways, ValidationSummary,
  RecoveryEvidence, LODODiseaseEvidence, LODORegion3Evidence,
  HiResLineageTier, HiResLineageCaveat, AnnotationBasis, AnnotationCaveat
)]
final_table[, EcotypeMember := EcotypeMember == "TRUE"]
final_table[, order_index := match(CellState, evidence$CellState)]
setorder(final_table, order_index)
final_table[, order_index := NULL]

# Step4. 保存主表和摘要
summary_table <- rbindlist(list(
  final_table[, .(Section = "AnnotationConfidence", Category = AnnotationConfidence, N = .N), by = AnnotationConfidence][, AnnotationConfidence := NULL],
  final_table[, .(Section = "StateEvidenceTier", Category = StateEvidenceTier, N = .N), by = StateEvidenceTier][, StateEvidenceTier := NULL],
  final_table[, .(Section = "EcotypeMember", Category = as.character(EcotypeMember), N = .N), by = EcotypeMember][, EcotypeMember := NULL]
), use.names = TRUE)

fwrite(final_table, file.path(annotation_dir, "state_annotation_final_table.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(final_table, file.path(annotation_dir, "state_annotation_final_table.csv"), quote = TRUE, na = "")
fwrite(summary_table, file.path(annotation_dir, "state_annotation_final_summary.tsv"), sep = "\t", quote = FALSE)
print(summary_table)

completion <- list(
  status = "PASS_STATE_ANNOTATION_FINAL_TABLE",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  states = nrow(final_table),
  labeled_states = sum(nzchar(final_table$FinalStateLabel)),
  high_confidence = sum(final_table$AnnotationConfidence == "HIGH"),
  moderate_confidence = sum(final_table$AnnotationConfidence == "MODERATE"),
  low_confidence = sum(final_table$AnnotationConfidence == "LOW")
)
jsonlite::write_json(completion, file.path(annotation_dir, "state_annotation_final_complete.json"), pretty = TRUE, auto_unbox = TRUE)
