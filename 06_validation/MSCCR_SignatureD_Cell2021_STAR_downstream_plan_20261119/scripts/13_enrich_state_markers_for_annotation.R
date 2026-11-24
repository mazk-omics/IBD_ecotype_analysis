#!/usr/bin/env Rscript
# 本脚本用固定 MSigDB Hallmark、GO Biological Process 和 Reactome 为 85 个 state 补充功能富集证据。
# 富集只服务于注释，不参与 recovery、LODO、state 定义或 ecotype 构建。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取正式 discovery marker
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 13_enrich_state_markers_for_annotation.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "digest", "jsonlite", "msigdbr")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

result_dir <- file.path(output_root, "10_state_annotation/pathway_enrichment")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
discovery_parts <- list()
gene_info_manifest <- list()
for (cell_type in names(rank_lookup)) {
  gene_info_file <- file.path(state_discovery_root, cell_type, rank_lookup[[cell_type]], "gene_info.txt")
  if (!file.exists(gene_info_file)) next
  gene_info <- as.data.table(read.delim(gene_info_file, check.names = FALSE, stringsAsFactors = FALSE))
  if (!all(c("Gene", "State", "MaxFC") %in% names(gene_info))) stop("Invalid formal gene_info: ", cell_type, call. = FALSE)
  gene_info <- gene_info[State != "-1" & !is.na(State), .(CellType = cell_type, State, Gene, DiscoveryMaxFC = as.numeric(MaxFC))]
  gene_info[, CellState := paste(CellType, State, sep = "_")]
  discovery_parts[[cell_type]] <- gene_info
  gene_info_manifest[[cell_type]] <- data.table(CellType = cell_type, Path = gene_info_file, SHA256 = digest::digest(gene_info_file, file = TRUE, algo = "sha256"))
}
discovery <- rbindlist(discovery_parts)
manifest <- rbindlist(gene_info_manifest)
if (uniqueN(discovery$CellState) != 85L || anyDuplicated(discovery[, .(CellType, Gene)])) stop("Formal discovery marker identity failed", call. = FALSE)

setorder(discovery, CellState, -DiscoveryMaxFC, Gene)
discovery[, MarkerOrder := seq_len(.N), by = CellState]
top_markers <- discovery[MarkerOrder <= 50L]
state_sizes <- top_markers[, .(TopMarkerN = .N), by = .(CellState, CellType, State)]
print(state_sizes[, summary(TopMarkerN)])

# Step2. 固定三个 MSigDB collection
hallmark <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"))
go_bp <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
reactome <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
gene_sets <- rbindlist(list(
  hallmark[, .(Collection = "HALLMARK", GeneSet = gs_name, Gene = gene_symbol, DatabaseVersion = db_version)],
  go_bp[, .(Collection = "GO_BP", GeneSet = gs_name, Gene = gene_symbol, DatabaseVersion = db_version)],
  reactome[, .(Collection = "REACTOME", GeneSet = gs_name, Gene = gene_symbol, DatabaseVersion = db_version)]
))
gene_sets <- unique(gene_sets[!is.na(Gene) & nzchar(Gene)])
if (uniqueN(gene_sets$DatabaseVersion) != 1L) stop("MSigDB database version is not unique", call. = FALSE)

database_manifest <- data.table(
  Package = "msigdbr", PackageVersion = as.character(packageVersion("msigdbr")), DatabaseVersion = unique(gene_sets$DatabaseVersion),
  Collection = c("HALLMARK", "GO_BP", "REACTOME"), GeneSets = gene_sets[, uniqueN(GeneSet), by = Collection]$V1[match(c("HALLMARK", "GO_BP", "REACTOME"), gene_sets[, uniqueN(GeneSet), by = Collection]$Collection)]
)
fwrite(manifest, file.path(result_dir, "discovery_gene_info_manifest.tsv"), sep = "\t", quote = FALSE)
fwrite(database_manifest, file.path(result_dir, "annotation_database_manifest.tsv"), sep = "\t", quote = FALSE)

# Step3. 在每个母细胞类型自己的 retained-gene universe 内做超几何检验
enrichment_parts <- list()
state_keys <- unique(top_markers[, .(CellState, CellType, State)])
for (state_index in seq_len(nrow(state_keys))) {
  cell_state <- state_keys$CellState[[state_index]]
  cell_type <- state_keys$CellType[[state_index]]
  markers <- unique(top_markers[CellState == cell_state]$Gene)
  universe <- unique(discovery[CellType == cell_type]$Gene)
  collection_sets <- gene_sets[Gene %in% universe, .(SetGenes = list(unique(Gene))), by = .(Collection, GeneSet)]
  collection_sets[, SetSize := lengths(SetGenes)]
  collection_sets <- collection_sets[SetSize >= 5L & SetSize <= 500L]
  collection_sets[, OverlapGenes := lapply(SetGenes, intersect, y = markers)]
  collection_sets[, OverlapN := lengths(OverlapGenes)]
  collection_sets <- collection_sets[OverlapN > 0L]
  collection_sets[, `:=`(
    CellState = cell_state, CellType = cell_type, State = state_keys$State[[state_index]],
    MarkerN = length(markers), UniverseN = length(universe),
    PValue = phyper(OverlapN - 1L, SetSize, length(universe) - SetSize, length(markers), lower.tail = FALSE),
    OverlapGenes = vapply(OverlapGenes, paste, character(1), collapse = ";")
  )]
  collection_sets[, FDR := p.adjust(PValue, method = "BH"), by = Collection]
  setorder(collection_sets, Collection, FDR, PValue, -OverlapN, GeneSet)
  collection_sets[, RankWithinCollection := seq_len(.N), by = Collection]
  enrichment_parts[[state_index]] <- collection_sets[, .(CellState, CellType, State, Collection, GeneSet, MarkerN, UniverseN, SetSize, OverlapN, OverlapGenes, PValue, FDR, RankWithinCollection)]
}
enrichment <- rbindlist(enrichment_parts)
if (uniqueN(enrichment$CellState) != 85L) stop("Some states have no annotation enrichment result", call. = FALSE)
fwrite(enrichment, file.path(result_dir, "state_marker_pathway_enrichment.tsv.gz"), sep = "\t", quote = FALSE)

# Step4. 输出显著结果和每个 collection 的前五项供人工审阅
top_enrichment <- enrichment[FDR < 0.05 | RankWithinCollection <= 5L]
top_enrichment[, Significant := FDR < 0.05]
fwrite(top_enrichment, file.path(result_dir, "state_marker_pathway_top_for_review.tsv"), sep = "\t", quote = FALSE)

compact <- top_enrichment[RankWithinCollection <= 5L, .(
  TopPathways = paste(paste0(GeneSet, " [q=", format(FDR, digits = 3, scientific = TRUE), "; ", OverlapGenes, "]"), collapse = " | ")
), by = .(CellState, Collection)]
compact <- dcast(compact, CellState ~ Collection, value.var = "TopPathways")
fwrite(compact, file.path(result_dir, "state_marker_pathway_compact.tsv"), sep = "\t", quote = FALSE)

completion <- list(
  status = "PASS_STATE_MARKER_PATHWAY_ENRICHMENT",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), states = uniqueN(enrichment$CellState),
  database_version = unique(gene_sets$DatabaseVersion), tested_rows = nrow(enrichment), significant_rows = sum(enrichment$FDR < 0.05),
  states_with_significant_hallmark = uniqueN(enrichment[Collection == "HALLMARK" & FDR < 0.05]$CellState),
  states_with_significant_go_bp = uniqueN(enrichment[Collection == "GO_BP" & FDR < 0.05]$CellState),
  states_with_significant_reactome = uniqueN(enrichment[Collection == "REACTOME" & FDR < 0.05]$CellState)
)
jsonlite::write_json(completion, file.path(result_dir, "state_marker_pathway_enrichment_complete.json"), pretty = TRUE, auto_unbox = TRUE)

print(completion)
cat("PASS_STATE_MARKER_PATHWAY_ENRICHMENT\n")
