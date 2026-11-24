#!/usr/bin/env Rscript
# 本脚本将 85 个 state 的冻结 marker 与 CellMarker 2.0、PanglaoDB 和 HPA v26 交叉核对。
# 数据库命中只形成注释证据，不自动改名，也不参与 state 或 ecotype 计算。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入并检查数据库快照
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 17_cross_reference_state_annotation_databases.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "digest", "jsonlite", "readxl")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

annotation_file <- file.path(RUN_ROOT, "output/10_state_annotation/state_annotation_final_table.tsv")
mapping_file <- file.path(RUN_ROOT, "inputs/state_annotation_database_celltype_map.tsv")
database_dir <- file.path(RUN_ROOT, "inputs/reference_databases/20261123")
cellmarker_file <- file.path(database_dir, "CellMarker2_Human.xlsx")
panglao_file <- file.path(database_dir, "PanglaoDB_markers_27_Mar_2020.tsv.gz")
hpa_file <- file.path(database_dir, "HPA_v26_rna_single_cell_type.tsv.zip")
output_dir <- file.path(RUN_ROOT, "output/10_state_annotation/literature_database_review/database_overlap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(annotation_file, mapping_file, cellmarker_file, panglao_file, hpa_file)
if (any(!file.exists(input_files))) stop("Missing database annotation input", call. = FALSE)
if (any(file.info(input_files)$size == 0)) stop("Database annotation input contains an empty file", call. = FALSE)

input_manifest <- data.table(
  Input = c("FrozenStateAnnotation", "DatabaseCellTypeMap", "CellMarker2", "PanglaoDB", "HPA26"),
  Path = input_files,
  Bytes = file.info(input_files)$size,
  SHA256 = vapply(input_files, digest::digest, character(1L), file = TRUE, algo = "sha256")
)
fwrite(input_manifest, file.path(output_dir, "database_reference_input_manifest.tsv"), sep = "\t", quote = FALSE)
print(input_manifest[, .(Input, Bytes, SHA256)])

# Step2. 展开冻结 marker，保留 discovery 与外部验证来源
annotation <- fread(annotation_file, colClasses = "character", check.names = FALSE)
mapping <- fread(mapping_file, colClasses = "character", check.names = FALSE)
if (nrow(annotation) != 85L || anyDuplicated(annotation$CellState)) stop("Expected 85 unique states", call. = FALSE)
if (nrow(mapping) != 60L || anyDuplicated(mapping[, .(CellType, Database)])) stop("Expected one mapping per cell type and database", call. = FALSE)
if (!setequal(unique(annotation$CellType), unique(mapping$CellType))) stop("Database mapping cell types differ from frozen states", call. = FALSE)

marker_parts <- list()
for (state_index in seq_len(nrow(annotation))) {
  discovery_genes <- unique(strsplit(annotation$KeyMarkers[[state_index]], ";", fixed = TRUE)[[1L]])
  external_genes <- unique(strsplit(annotation$ExternalValidatedMarkers[[state_index]], ";", fixed = TRUE)[[1L]])
  state_genes <- unique(c(discovery_genes, external_genes))
  state_genes <- state_genes[nzchar(state_genes)]
  marker_parts[[state_index]] <- data.table(
    CellState = annotation$CellState[[state_index]], CellType = annotation$CellType[[state_index]], Gene = state_genes,
    IsDiscoveryMarker = state_genes %in% discovery_genes, IsExternalValidatedMarker = state_genes %in% external_genes
  )
}
state_markers <- rbindlist(marker_parts)
if (uniqueN(state_markers$CellState) != 85L) stop("State marker expansion failed", call. = FALSE)
print(state_markers[, .N, by = CellState][, summary(N)])

# Step3. 读取三个冻结数据库
cellmarker <- as.data.table(readxl::read_excel(cellmarker_file, sheet = "human"))
cellmarker <- cellmarker[species == "Human" & !is.na(Symbol) & nzchar(Symbol), .(
  Gene = Symbol, DatabaseCellName = cell_name, DatabaseTissue = tissue_type,
  DatabaseContextDetail = marker_source, ReferenceID = as.character(PMID), ReferenceTitle = Title
)]
cellmarker <- unique(cellmarker)

panglao <- fread(panglao_file, tmpdir = output_dir)
setnames(panglao, c("official gene symbol", "cell type", "canonical marker"), c("Gene", "DatabaseCellName", "CanonicalMarker"))
panglao <- panglao[grepl("Hs", species, fixed = TRUE) & !is.na(Gene) & nzchar(Gene), .(
  Gene, DatabaseCellName, DatabaseTissue = organ, CanonicalMarker,
  SensitivityHuman = sensitivity_human, SpecificityHuman = specificity_human
)]
panglao <- unique(panglao)

hpa_extract_dir <- file.path(output_dir, "hpa_extract_tmp")
if (dir.exists(hpa_extract_dir)) {
  output_check <- tolower(gsub("\\\\", "/", normalizePath(output_dir, mustWork = TRUE)))
  extract_check <- tolower(gsub("\\\\", "/", normalizePath(hpa_extract_dir, mustWork = TRUE)))
  if (!startsWith(extract_check, paste0(output_check, "/"))) stop("Unsafe HPA temporary directory", call. = FALSE)
  unlink(hpa_extract_dir, recursive = TRUE, force = TRUE)
}
dir.create(hpa_extract_dir, recursive = TRUE, showWarnings = FALSE)
unzip(hpa_file, exdir = hpa_extract_dir)
hpa_tsv <- file.path(hpa_extract_dir, "rna_single_cell_type.tsv")
if (!file.exists(hpa_tsv)) stop("HPA archive content changed", call. = FALSE)
hpa <- fread(hpa_tsv, select = c("Gene name", "Cell type", "nCPM"))
setnames(hpa, c("Gene name", "Cell type"), c("Gene", "DatabaseCellName"))
hpa[, nCPM := as.numeric(nCPM)]
if (uniqueN(hpa$DatabaseCellName) != 155L) stop("Unexpected HPA v26 cell-type count", call. = FALSE)
print(c(CellMarkerRows = nrow(cellmarker), PanglaoRows = nrow(panglao), HPARows = nrow(hpa)))

# Step4. 逐母细胞类型匹配数据库，分别保留肠道优先证据和广义证据
evidence_parts <- list()
summary_parts <- list()
part_index <- 0L

for (map_index in seq_len(nrow(mapping))) {
  current_type <- mapping$CellType[[map_index]]
  current_database <- mapping$Database[[map_index]]
  cell_regex <- mapping$CellNameRegex[[map_index]]
  tissue_regex <- mapping$PreferredTissueRegex[[map_index]]
  current_markers <- state_markers[CellType == current_type]

  if (current_database == "CellMarker2") {
    reference <- cellmarker[grepl(cell_regex, DatabaseCellName, ignore.case = TRUE, perl = TRUE)]
    reference[, PreferredContext := nzchar(tissue_regex) & grepl(tissue_regex, DatabaseTissue, ignore.case = TRUE, perl = TRUE)]
    matched <- merge(current_markers, reference, by = "Gene")
    if (nrow(matched)) {
      matched[, `:=`(Database = current_database, HPA_TargetMaxNCPM = NA_real_, HPA_BestOtherNCPM = NA_real_, HPA_TargetRank = NA_integer_, HPA_TargetVsMedianFold = NA_real_)]
      matched[, CanonicalMarker := NA_integer_]
    }
  }

  if (current_database == "PanglaoDB") {
    reference <- panglao[grepl(cell_regex, DatabaseCellName, ignore.case = TRUE, perl = TRUE)]
    reference[, PreferredContext := nzchar(tissue_regex) & grepl(tissue_regex, DatabaseTissue, ignore.case = TRUE, perl = TRUE)]
    matched <- merge(current_markers, reference, by = "Gene")
    if (nrow(matched)) {
      matched[, `:=`(Database = current_database, DatabaseContextDetail = paste0("canonical=", CanonicalMarker, "; sensitivity_human=", SensitivityHuman, "; specificity_human=", SpecificityHuman), ReferenceID = "PanglaoDB-2020-03-27", ReferenceTitle = "")]
      matched[, `:=`(HPA_TargetMaxNCPM = NA_real_, HPA_BestOtherNCPM = NA_real_, HPA_TargetRank = NA_integer_, HPA_TargetVsMedianFold = NA_real_)]
    }
  }

  if (current_database == "HPA26") {
    target_types <- unique(hpa[grepl(cell_regex, DatabaseCellName, ignore.case = TRUE, perl = TRUE), DatabaseCellName])
    gene_expression <- hpa[Gene %in% current_markers$Gene]
    target_expression <- gene_expression[DatabaseCellName %in% target_types, .(HPA_TargetMaxNCPM = max(nCPM, na.rm = TRUE), DatabaseCellName = DatabaseCellName[which.max(nCPM)]), by = Gene]
    other_expression <- gene_expression[!DatabaseCellName %in% target_types, .(HPA_BestOtherNCPM = max(nCPM, na.rm = TRUE), HPA_MedianAllNCPM = median(nCPM, na.rm = TRUE)), by = Gene]
    rank_expression <- gene_expression[, .(HPA_TargetRank = min(rank(-nCPM, ties.method = "min")[DatabaseCellName %in% target_types])), by = Gene]
    reference <- merge(target_expression, other_expression, by = "Gene", all = TRUE)
    reference <- merge(reference, rank_expression, by = "Gene", all = TRUE)
    reference[, HPA_TargetVsMedianFold := HPA_TargetMaxNCPM / (HPA_MedianAllNCPM + 0.1)]
    reference[, `:=`(DatabaseTissue = "HPA aggregated 36 datasets", PreferredContext = HPA_TargetRank <= 10L & HPA_TargetMaxNCPM >= 1, CanonicalMarker = NA_integer_, DatabaseContextDetail = "HPA v26 target expression rank", ReferenceID = "HPA-v26", ReferenceTitle = "Human Protein Atlas Single Cell Type Resource")]
    matched <- merge(current_markers, reference, by = "Gene")
    if (nrow(matched)) matched[, Database := current_database]
  }

  if (!exists("matched") || !nrow(matched)) {
    matched <- data.table()
  } else {
    part_index <- part_index + 1L
    evidence_parts[[part_index]] <- matched[, .(
      CellState, CellType, Gene, IsDiscoveryMarker, IsExternalValidatedMarker, Database,
      DatabaseCellName, DatabaseTissue, PreferredContext, CanonicalMarker, DatabaseContextDetail,
      ReferenceID, ReferenceTitle, HPA_TargetMaxNCPM, HPA_BestOtherNCPM, HPA_TargetRank, HPA_TargetVsMedianFold
    )]
  }

  state_keys <- unique(current_markers[, .(CellState, CellType)])
  if (nrow(matched)) {
    current_summary <- matched[, .(
      DatabaseMatchedMarkerN = uniqueN(Gene),
      PreferredContextMarkerN = uniqueN(Gene[PreferredContext == TRUE]),
      DatabaseMatchedMarkers = paste(sort(unique(Gene)), collapse = ";"),
      PreferredContextMarkers = paste(sort(unique(Gene[PreferredContext == TRUE])), collapse = ";"),
      ReferenceCellTypes = paste(sort(unique(DatabaseCellName)), collapse = ";"),
      ReferenceIDs = paste(sort(unique(ReferenceID[!is.na(ReferenceID) & nzchar(ReferenceID)])), collapse = ";")
    ), by = .(CellState, CellType)]
    current_summary <- merge(state_keys, current_summary, by = c("CellState", "CellType"), all.x = TRUE)
  } else {
    current_summary <- copy(state_keys)
    current_summary[, `:=`(DatabaseMatchedMarkerN = 0L, PreferredContextMarkerN = 0L, DatabaseMatchedMarkers = "", PreferredContextMarkers = "", ReferenceCellTypes = "", ReferenceIDs = "")]
  }
  current_summary[is.na(DatabaseMatchedMarkerN), `:=`(DatabaseMatchedMarkerN = 0L, PreferredContextMarkerN = 0L, DatabaseMatchedMarkers = "", PreferredContextMarkers = "", ReferenceCellTypes = "", ReferenceIDs = "")]
  current_summary[, `:=`(Database = current_database, InterpretationLimit = mapping$InterpretationLimit[[map_index]])]
  summary_parts[[map_index]] <- current_summary
  rm(matched)
}

database_evidence <- rbindlist(evidence_parts, fill = TRUE)
database_summary <- rbindlist(summary_parts, fill = TRUE)
setorder(database_evidence, CellState, Database, -PreferredContext, Gene)
setorder(database_summary, CellState, Database)

if (uniqueN(database_summary$CellState) != 85L || nrow(database_summary) != 255L) stop("Database summary grid is incomplete", call. = FALSE)
fwrite(database_evidence, file.path(output_dir, "state_annotation_database_marker_evidence.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(database_summary, file.path(output_dir, "state_annotation_database_overlap_summary.tsv"), sep = "\t", quote = FALSE, na = "")

# Step5. 输出数据库覆盖审计并清理临时解压文件
coverage <- database_summary[, .(
  States = uniqueN(CellState), StatesWithAnyMatch = uniqueN(CellState[DatabaseMatchedMarkerN > 0L]),
  StatesWithPreferredMatch = uniqueN(CellState[PreferredContextMarkerN > 0L]),
  MatchedStateMarkerPairs = sum(DatabaseMatchedMarkerN), PreferredStateMarkerPairs = sum(PreferredContextMarkerN)
), by = Database]
fwrite(coverage, file.path(output_dir, "state_annotation_database_coverage.tsv"), sep = "\t", quote = FALSE)
print(coverage)

output_check <- tolower(gsub("\\\\", "/", normalizePath(output_dir, mustWork = TRUE)))
extract_check <- tolower(gsub("\\\\", "/", normalizePath(hpa_extract_dir, mustWork = TRUE)))
if (!startsWith(extract_check, paste0(output_check, "/"))) stop("Unsafe HPA temporary cleanup path", call. = FALSE)
unlink(hpa_extract_dir, recursive = TRUE, force = TRUE)

completion <- list(
  status = "PASS_STATE_ANNOTATION_DATABASE_CROSS_REFERENCE",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), states = uniqueN(database_summary$CellState),
  summary_rows = nrow(database_summary), evidence_rows = nrow(database_evidence), databases = unique(database_summary$Database)
)
jsonlite::write_json(completion, file.path(output_dir, "state_annotation_database_cross_reference_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_STATE_ANNOTATION_DATABASE_CROSS_REFERENCE\n")
