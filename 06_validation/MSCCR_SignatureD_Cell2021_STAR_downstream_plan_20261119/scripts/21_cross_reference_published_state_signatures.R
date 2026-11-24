#!/usr/bin/env Rscript
# 本脚本将当前 85 个 state 的可用 marker 与上一轮冻结的 published signatures 做同谱系描述性重叠。
# 结果只提供文献锚定候选，不进行自动改名或证据等级升级。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 21_cross_reference_published_state_signatures.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

annotation_file <- file.path(RUN_ROOT, "output/10_state_annotation/state_annotation_final_table.tsv")
published_file <- normalizePath(file.path(
  RUN_ROOT, "../../analysis_reports/legacy_known_state_correspondence/run_20261013_01/literature_resources/PUBLISHED_SIGNATURE_GENES.tsv"
), mustWork = TRUE)
output_dir <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review/published_signature_overlap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(annotation_file, published_file)
if (any(!file.exists(input_files))) stop("Missing published-signature overlap input", call. = FALSE)

annotation <- fread(annotation_file, colClasses = "character", check.names = FALSE)
published <- fread(published_file, colClasses = "character", check.names = FALSE)
if (nrow(annotation) != 85L || anyDuplicated(annotation$CellState)) stop("Expected 85 unique states", call. = FALSE)

# Step2. 展开当前 marker 和文献正向 signature
marker_parts <- vector("list", nrow(annotation))
for (state_index in seq_len(nrow(annotation))) {
  key_genes <- strsplit(annotation$KeyMarkers[[state_index]], ";", fixed = TRUE)[[1L]]
  external_genes <- strsplit(annotation$ExternalValidatedMarkers[[state_index]], ";", fixed = TRUE)[[1L]]
  genes <- sort(unique(c(key_genes, external_genes)))
  genes <- genes[nzchar(genes)]
  marker_parts[[state_index]] <- data.table(
    CellState = annotation$CellState[[state_index]], CellType = annotation$CellType[[state_index]], Gene = genes,
    IsKeyMarker = genes %in% key_genes, IsExternalMarker = genes %in% external_genes
  )
}
markers <- rbindlist(marker_parts)

published <- published[Direction %in% c("up", "positive_program_list") & nzchar(Gene) & nzchar(PublishedState)]
published <- published[nzchar(ExpectedParent) & ExpectedParent != "Any_shared_program"]
published[, ParentMatchKey := gsub(";", "|", ExpectedParent, fixed = TRUE)]
published <- unique(published[, .(Publication, DOI, PublishedState, ExpectedParent, ParentMatchKey, Gene, SignatureScope)])

# Step3. 计算同谱系 marker 重叠
match_parts <- list()
part_index <- 0L
for (state_index in seq_len(nrow(annotation))) {
  current_state <- annotation$CellState[[state_index]]
  current_type <- annotation$CellType[[state_index]]
  current_genes <- markers[CellState == current_state, Gene]
  current_signatures <- published[grepl(paste0("(^|\\|)", current_type, "($|\\|)"), ParentMatchKey, perl = TRUE)]
  if (!nrow(current_signatures)) next
  signature_groups <- unique(current_signatures[, .(Publication, DOI, PublishedState, ExpectedParent, SignatureScope)])
  for (signature_index in seq_len(nrow(signature_groups))) {
    current_signature <- signature_groups[signature_index]
    signature_genes <- current_signatures[
      Publication == current_signature$Publication & PublishedState == current_signature$PublishedState & ExpectedParent == current_signature$ExpectedParent,
      unique(Gene)
    ]
    overlap <- sort(intersect(current_genes, signature_genes))
    if (!length(overlap)) next
    part_index <- part_index + 1L
    match_parts[[part_index]] <- data.table(
      CellState = current_state, CellType = current_type,
      Publication = current_signature$Publication, DOI = current_signature$DOI,
      PublishedState = current_signature$PublishedState, PublishedExpectedParent = current_signature$ExpectedParent,
      CurrentMarkerN = length(current_genes), PublishedSignatureN = length(signature_genes), OverlapN = length(overlap),
      OverlapGenes = paste(overlap, collapse = ";"),
      CurrentMarkerCoverage = length(overlap) / length(current_genes),
      PublishedSignatureCoverage = length(overlap) / length(signature_genes),
      Jaccard = length(overlap) / length(union(current_genes, signature_genes)),
      SignatureScope = current_signature$SignatureScope
    )
  }
}

matches <- rbindlist(match_parts, fill = TRUE)
setorder(matches, CellState, -OverlapN, -Jaccard, Publication, PublishedState)
fwrite(matches, file.path(output_dir, "current85_published_signature_overlap_all.tsv"), sep = "\t", quote = FALSE, na = "")

best <- matches[, .SD[1L], by = CellState]
state_keys <- annotation[, .(CellState, CellType, AnnotationConfidence)]
best <- merge(state_keys, best, by = c("CellState", "CellType"), all.x = TRUE, sort = FALSE)
best[is.na(OverlapN), `:=`(OverlapN = 0L, OverlapGenes = "", Publication = "", DOI = "", PublishedState = "", PublishedExpectedParent = "", SignatureScope = "")]
setorder(best, CellType, CellState)
fwrite(best, file.path(output_dir, "current85_published_signature_best_match.tsv"), sep = "\t", quote = FALSE, na = "")

# Step4. 输出覆盖审计
coverage <- best[, .(
  States = .N, AnyOverlap = sum(OverlapN >= 1L), AtLeastTwoGenes = sum(OverlapN >= 2L), AtLeastThreeGenes = sum(OverlapN >= 3L)
), by = CellType]
setorder(coverage, CellType)
fwrite(coverage, file.path(output_dir, "current85_published_signature_coverage.tsv"), sep = "\t", quote = FALSE)

manifest <- data.table(
  Input = c("FrozenStateAnnotation", "FrozenPublishedSignatureGenes"),
  Path = input_files, Bytes = file.info(input_files)$size,
  SHA256 = vapply(input_files, digest::digest, character(1L), file = TRUE, algo = "sha256")
)
fwrite(manifest, file.path(output_dir, "published_signature_overlap_input_manifest.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_CURRENT85_PUBLISHED_SIGNATURE_CROSS_REFERENCE",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), states = nrow(annotation),
  published_signatures = uniqueN(published[, .(Publication, PublishedState, ExpectedParent)]),
  states_with_any_overlap = sum(best$OverlapN >= 1L), states_with_at_least_two_genes = sum(best$OverlapN >= 2L),
  automatic_grade_change = FALSE, frozen_annotation_changed = FALSE
)
jsonlite::write_json(completion, file.path(output_dir, "published_signature_overlap_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(coverage)
print(completion)
cat("PASS_CURRENT85_PUBLISHED_SIGNATURE_CROSS_REFERENCE\n")
