#!/usr/bin/env Rscript
# 本脚本合并 marker、recovery、bulk、LODO、表型、文献和 ecotype 证据，生成 85-state 注释审阅优先级。
# 优先级只决定人工审阅顺序，不自动升级标签或修改冻结模型。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取正式注释和验证证据
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 54_build_state_annotation_priority_master.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

annotation_file <- file.path(RUN_ROOT, "output/10_state_annotation/state_annotation_final_table.tsv")
literature_file <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review/state_annotation_manuscript_label_table_v3_cell2021_aligned.tsv")
published_file <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review/published_signature_overlap/current85_published_signature_best_match.tsv")
scrna_file <- file.path(RUN_ROOT, "output/08_scRNA_recovery_summary/scRNA_recovery_stouffer_meta.tsv")
nonibd_file <- file.path(RUN_ROOT, "output/38_integrated_nonIBD_validation_evidence/integrated_nonIBD_state_validation_evidence.tsv")
lodo_file <- file.path(RUN_ROOT, "output/09_marker_recovery/ranking_and_LODO/marker_LODO_per_state.tsv")
pathway_file <- file.path(RUN_ROOT, "output/10_state_annotation/pathway_enrichment/state_marker_pathway_compact.tsv")
representation_file <- file.path(RUN_ROOT, "output/39_discovery_state_representation_Cell2021/discovery_state_representation_fisher.tsv")
output_dir <- file.path(RUN_ROOT, "output/41_state_annotation_priority_master")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(annotation_file, literature_file, published_file, scrna_file, nonibd_file, lodo_file, pathway_file, representation_file)
if (any(!file.exists(input_files))) stop("Missing annotation evidence: ", paste(input_files[!file.exists(input_files)], collapse = "; "), call. = FALSE)

annotation <- fread(annotation_file)
literature <- fread(literature_file)
published <- fread(published_file)
scrna <- fread(scrna_file)
nonibd <- fread(nonibd_file)
lodo <- fread(lodo_file)
pathway <- fread(pathway_file)
representation <- fread(representation_file)

if (nrow(annotation) != 85L || anyDuplicated(annotation$CellState)) stop("Formal annotation is not 85 unique states", call. = FALSE)
if (nrow(literature) != 85L || anyDuplicated(literature$CellState)) stop("Literature review is not 85 unique states", call. = FALSE)

# Step2. 整理逐 state 的外部和表型证据
scrna[, CellState := paste(CellType, State, sep = "_")]
scrna <- scrna[, .(
  CellState, DatasetsEvaluable, DatasetsIndividuallySignificant, MetaZ,
  MetaFDR_BH, MetaSignificant_BH, EvidenceScope
)]

nonibd <- nonibd[, .(
  CellState, BulkZ, BulkSignificant, ScRNADatasetsEvaluable,
  ScRNADatasetsSignificant, ExternalBulkValidated, ExternalScRNAValidated,
  ExternalCrossModalityValidated, MSCCRRegionsObserved, MSCCRAllRegionsObserved,
  EvidenceClass
)]

lodo <- dcast(lodo, CellState ~ Variant, value.var = "Concordance")
if ("Disease" %in% names(lodo)) setnames(lodo, "Disease", "LODODiseaseConcordance")
if ("Region3" %in% names(lodo)) setnames(lodo, "Region3", "LODORegionConcordance")

inflamed <- representation[Axis == "Inflammation" & Level == "Inflamed" & Status == "EVALUABLE", .(
  CellState, InflamedOR = OddsRatio, InflamedFDR = FDR_BH_across_cell_states,
  InflamedScore = EnrichmentScore, InflamedSignificant = SignificantEnrichment
)]
noninflamed <- representation[Axis == "Inflammation" & Level == "Non-inflamed" & Status == "EVALUABLE", .(
  CellState, NoninflamedOR = OddsRatio, NoninflamedFDR = FDR_BH_across_cell_states,
  NoninflamedScore = EnrichmentScore, NoninflamedSignificant = SignificantEnrichment
)]
inflammation <- merge(inflamed, noninflamed, by = "CellState", all = TRUE, sort = FALSE)
inflammation[, InflammationSignificant := InflamedSignificant %in% TRUE | NoninflamedSignificant %in% TRUE]
inflammation[, InflammationDirection := fifelse(
  InflamedSignificant %in% TRUE, "Inflamed",
  fifelse(NoninflamedSignificant %in% TRUE, "Non-inflamed", "Not significant")
)]

literature <- literature[, .(
  CellState, LiteratureInterpretationClass = InterpretationClass,
  LiteratureEvidenceGrade, InterpretableForManuscript, RecommendedUse,
  ManuscriptLabel, ManuscriptLabel_CN, ReviewerDecision, FrozenLabel,
  DatabaseMatchedMarkers, LiteratureMatchedMarkers, DatabaseSourceIDs,
  LiteratureSourceIDs, EvidenceDirection, ConflictSummary, ReviewerNote
)]

published <- published[, .(
  CellState, PublishedMatchPublication = Publication, PublishedMatchDOI = DOI,
  PublishedState, PublishedExpectedParent, PublishedOverlapN = OverlapN,
  PublishedOverlapGenes = OverlapGenes, PublishedJaccard = Jaccard
)]

# Step3. 合并证据并计算人工审阅顺序
evidence <- Reduce(
  function(left, right) merge(left, right, by = "CellState", all.x = TRUE, sort = FALSE),
  list(annotation, literature, published, scrna, nonibd, lodo, pathway, inflammation)
)
if (nrow(evidence) != 85L || anyDuplicated(evidence$CellState)) stop("Evidence merge changed the 85-state universe", call. = FALSE)

evidence[, ReviewPriorityScore :=
  3L * as.integer(ExternalCrossModalityValidated %in% TRUE) +
  2L * as.integer(ExternalBulkValidated %in% TRUE) +
  2L * as.integer(MetaSignificant_BH %in% TRUE) +
  1L * as.integer(MSCCRAllRegionsObserved %in% TRUE) +
  1L * as.integer(fcoalesce(LODODiseaseConcordance >= 0.75, FALSE)) +
  1L * as.integer(fcoalesce(LODORegionConcordance >= 0.75, FALSE)) +
  2L * as.integer(InflammationSignificant %in% TRUE) +
  2L * as.integer(EcotypeMember %in% TRUE)
]

evidence[, ReviewPriorityTier := fifelse(
  LiteratureEvidenceGrade == "D" & EcotypeMember %in% TRUE & ReviewPriorityScore >= 8L,
  "P1_UNRESOLVED_ECOTYPE_MEMBER_HIGH_EVIDENCE",
  fifelse(
    LiteratureEvidenceGrade == "D" & ReviewPriorityScore >= 6L,
    "P2_UNRESOLVED_REPRODUCIBLE",
    fifelse(
      LiteratureEvidenceGrade == "C" & EcotypeMember %in% TRUE & ReviewPriorityScore >= 8L,
      "P3_FUNCTIONAL_PROGRAM_ECOTYPE_REVIEW",
      "P4_RETAIN_CURRENT_LABEL"
    )
  )
)]
evidence[, AutomaticLabelChange := FALSE]
evidence[, ManualReviewRequired := ReviewPriorityTier != "P4_RETAIN_CURRENT_LABEL"]
setorder(evidence, ReviewPriorityTier, -ReviewPriorityScore, CellType, State)

# Step4. 保存主表和审阅清单
priority <- evidence[LiteratureEvidenceGrade == "D"][order(ReviewPriorityTier, -ReviewPriorityScore, CellType, State)]
summary_table <- evidence[, .(
  States = .N,
  EcotypeMembers = sum(EcotypeMember %in% TRUE),
  InflammationAssociated = sum(InflammationSignificant %in% TRUE),
  CrossModalityValidated = sum(ExternalCrossModalityValidated %in% TRUE),
  ManualReviewRequired = sum(ManualReviewRequired)
), by = .(LiteratureEvidenceGrade, ReviewPriorityTier)]

fwrite(evidence, file.path(output_dir, "state_annotation_evidence_master_v4.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(priority, file.path(output_dir, "state_annotation_priority_unresolved.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(summary_table, file.path(output_dir, "state_annotation_priority_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

completion <- list(
  status = "PASS_STATE_ANNOTATION_PRIORITY_MASTER_V4",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  formal_states = nrow(evidence),
  unresolved_states = nrow(priority),
  automatic_label_change = FALSE,
  primary_next_action = "manual literature and marker review of P1 unresolved ecotype-member states"
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_priority_master_v4_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(summary_table)
print(priority[, .(
  CellState, Ecotype, ReviewPriorityScore, ReviewPriorityTier,
  InflammationDirection, ExternalCrossModalityValidated, MetaZ,
  LODODiseaseConcordance, LODORegionConcordance, ConflictSummary
)][1:min(.N, 20L)])
cat("PASS_STATE_ANNOTATION_PRIORITY_MASTER_V4\n")
