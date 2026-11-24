#!/usr/bin/env Rscript
# 本脚本用冻结 Ensembl75 身份表和唯一 HGNC alias 解析外部 scRNA 到正式 W 的基因映射。
# 映射冲突和未解析基因完整保留；不读取或修改表达量，不启动 recovery。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取前序审计、冻结模型和注释依赖
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 03e_resolve_scrna_gene_identity.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "jsonlite", "AnnotationDbi", "org.Hs.eg.db")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

audit_dir <- file.path(output_root, "03_scrna_input_audit")
source_mapping_file <- file.path(audit_dir, "scrna_source_gene_mapping.tsv.gz")
if (!file.exists(source_mapping_file) || !file.exists(ensembl75_identity_file)) stop("Frozen gene mapping inputs missing", call. = FALSE)

source_mapping <- fread(source_mapping_file, colClasses = "character")
source_mapping[, GeneIDNorm := sub("[.].*$", "", GeneID)]
source_mapping[, SyntacticGene := make.names(Gene, unique = FALSE)]
source_mapping[, SourceKey := paste(GeneID, Gene, sep = "|")]

identity <- fread(ensembl75_identity_file, colClasses = "character")
required_identity <- c("msccr_ensembl_id", "ensembl75_gene_name", "mapping_status")
if (!all(required_identity %in% names(identity))) stop("Ensembl75 identity columns changed", call. = FALSE)
identity[, GeneIDNorm := sub("[.].*$", "", msccr_ensembl_id)]
identity[, ModelGene := make.names(ensembl75_gene_name, unique = FALSE)]
identity <- unique(identity[, .(GeneIDNorm, OldSymbol = ensembl75_gene_name, ModelGene, MappingStatus = mapping_status)])

rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_cell_type <- rank_data[[1L]]
rank_selected <- trimws(rank_data[[2L]])
if (anyDuplicated(rank_cell_type) || length(rank_cell_type) != 20L) stop("Frozen rank table changed", call. = FALSE)

model_parts <- vector("list", length(rank_cell_type))
for (i in seq_along(rank_cell_type)) {
  cell_type <- rank_cell_type[[i]]
  selected_rank <- rank_selected[[i]]
  W <- read.delim(file.path(state_discovery_root, cell_type, selected_rank, "W.txt"), check.names = FALSE)
  w_genes <- sub("__pos$", "", rownames(W)[grepl("__pos$", rownames(W))])
  if (anyDuplicated(w_genes)) stop("Duplicated W genes: ", cell_type, call. = FALSE)
  model_parts[[i]] <- data.table(CellType = cell_type, FrozenRank = selected_rank, ModelGene = w_genes)
}
model_genes <- rbindlist(model_parts)

# Step2. 冻结旧 symbol 到当前 symbol 的唯一 alias
model_identity <- unique(merge(model_genes, identity, by = "ModelGene", all.x = TRUE, allow.cartesian = TRUE))
old_symbols <- sort(unique(model_identity[!is.na(OldSymbol) & OldSymbol != "", OldSymbol]))
current_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
alias_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ALIAS")

direct_symbols <- intersect(old_symbols, current_keys)
direct_alias <- data.table(OldSymbol = direct_symbols, CurrentSymbol = direct_symbols, AliasClass = "DIRECT_CURRENT_SYMBOL")

historical_symbols <- intersect(setdiff(old_symbols, current_keys), alias_keys)
alias_raw <- as.data.table(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = historical_symbols,
  keytype = "ALIAS",
  columns = "SYMBOL"
)))
alias_summary <- alias_raw[!is.na(SYMBOL), .(
  CurrentSymbolN = uniqueN(SYMBOL),
  CurrentSymbol = if (uniqueN(SYMBOL) == 1L) unique(SYMBOL) else paste(sort(unique(SYMBOL)), collapse = ";")
), by = .(OldSymbol = ALIAS)]
historical_alias <- alias_summary[CurrentSymbolN == 1L, .(OldSymbol, CurrentSymbol, AliasClass = "HGNC_ALIAS_UNIQUE")]
alias_reference <- unique(rbindlist(list(direct_alias, historical_alias), use.names = TRUE))

alias_audit <- data.table(OldSymbol = old_symbols)
alias_audit <- merge(alias_audit, alias_reference, by = "OldSymbol", all.x = TRUE)
alias_audit <- merge(alias_audit, alias_summary[, .(OldSymbol, AliasTargetN = CurrentSymbolN, AliasTargets = CurrentSymbol)], by = "OldSymbol", all.x = TRUE)
alias_audit[is.na(AliasClass), AliasClass := ifelse(AliasTargetN > 1L, "HGNC_ALIAS_AMBIGUOUS", "UNRESOLVED_ALIAS")]

# Step3. 按 dataset 和父类型建立候选映射
candidate_parts <- list()
datasets <- sort(unique(source_mapping$Dataset))
for (dataset in datasets) {
  source_dataset <- source_mapping[Dataset == dataset]
  source_ensembl <- merge(source_dataset, identity[, .(GeneIDNorm, ModelGene)], by = "GeneIDNorm", allow.cartesian = TRUE)

  for (cell_type in rank_cell_type) {
    w_genes <- model_genes[CellType == cell_type, ModelGene]

    syntactic <- source_dataset[SyntacticGene %in% w_genes, .(
      Dataset,
      CellType = cell_type,
      ModelGene = SyntacticGene,
      SourceKey,
      SourceGeneID = GeneID,
      SourceGene = Gene,
      MappingClass = ifelse(Gene == SyntacticGene, "EXACT_SYMBOL", "SYNTACTIC_SYMBOL"),
      MappingPriority = ifelse(Gene == SyntacticGene, 1L, 2L)
    )]

    ensembl <- source_ensembl[ModelGene %in% w_genes, .(
      Dataset,
      CellType = cell_type,
      ModelGene,
      SourceKey,
      SourceGeneID = GeneID,
      SourceGene = Gene,
      MappingClass = "ENSEMBL75_IDENTITY",
      MappingPriority = 3L
    )]

    alias_model <- unique(model_identity[CellType == cell_type & ModelGene %in% w_genes & !is.na(OldSymbol), .(ModelGene, OldSymbol)])
    alias_model <- merge(alias_model, alias_reference, by = "OldSymbol")
    alias <- merge(source_dataset, alias_model, by.x = "Gene", by.y = "CurrentSymbol", allow.cartesian = TRUE)
    alias <- alias[, .(
      Dataset,
      CellType = cell_type,
      ModelGene,
      SourceKey,
      SourceGeneID = GeneID,
      SourceGene = Gene,
      MappingClass = AliasClass,
      MappingPriority = 4L
    )]

    candidate_parts[[paste(dataset, cell_type, sep = "::")]] <- unique(rbindlist(list(syntactic, ensembl, alias), use.names = TRUE))
  }
}

candidates <- rbindlist(candidate_parts, use.names = TRUE)
best_priority <- candidates[, .(BestPriority = min(MappingPriority)), by = .(Dataset, CellType, SourceKey)]
selected_candidates <- merge(candidates, best_priority, by = c("Dataset", "CellType", "SourceKey"))
selected_candidates <- selected_candidates[MappingPriority == BestPriority]
source_conflicts <- selected_candidates[, .(
  ModelGeneN = uniqueN(ModelGene),
  ModelGenes = paste(sort(unique(ModelGene)), collapse = ";")
), by = .(Dataset, CellType, SourceKey)][ModelGeneN > 1L]

setorder(selected_candidates, Dataset, CellType, ModelGene, SourceKey, MappingPriority)
accepted <- selected_candidates[, .SD[1L], by = .(Dataset, CellType, ModelGene, SourceKey)]
resolved_counts <- accepted[, .(
  SourceRows = uniqueN(SourceKey),
  MappingClasses = paste(sort(unique(MappingClass)), collapse = ";")
), by = .(Dataset, CellType, ModelGene)]

resolution_grid <- merge(CJ(Dataset = datasets, CellType = rank_cell_type, unique = TRUE), model_genes, by = "CellType", allow.cartesian = TRUE)
resolution <- merge(resolution_grid, resolved_counts, by = c("Dataset", "CellType", "ModelGene"), all.x = TRUE)
resolution[is.na(SourceRows), `:=`(SourceRows = 0L, MappingClasses = "UNRESOLVED")]
resolution[, Resolved := SourceRows > 0L]

coverage <- resolution[, .(
  FrozenRank = unique(FrozenRank),
  WGenes = .N,
  WGenesResolved = sum(Resolved),
  WGenesUnresolved = sum(!Resolved),
  WGeneCoverage = mean(Resolved)
), by = .(Dataset, CellType)]

# Step4. 保存映射、冲突和最终 gate
data.table::fwrite(alias_audit, file.path(audit_dir, "scrna_model_symbol_alias_audit.tsv"), sep = "\t")
data.table::fwrite(candidates, file.path(audit_dir, "scrna_gene_mapping_candidates.tsv.gz"), sep = "\t")
data.table::fwrite(selected_candidates, file.path(audit_dir, "scrna_gene_mapping_selected_candidates.tsv.gz"), sep = "\t")
data.table::fwrite(source_conflicts, file.path(audit_dir, "scrna_gene_mapping_source_conflicts.tsv"), sep = "\t")
data.table::fwrite(accepted, file.path(audit_dir, "scrna_resolved_model_gene_mapping.tsv.gz"), sep = "\t")
data.table::fwrite(resolution, file.path(audit_dir, "scrna_model_gene_resolution.tsv.gz"), sep = "\t")
data.table::fwrite(coverage, file.path(audit_dir, "scrna_resolved_model_gene_coverage.tsv"), sep = "\t")

minimum_coverage <- min(coverage$WGeneCoverage)
gate_pass <- nrow(source_conflicts) == 0L && minimum_coverage > 0
low_coverage_groups <- sum(coverage$WGeneCoverage < 0.90)
ensembl75_identity_sha256 <- strsplit(system2("sha256sum", ensembl75_identity_file, stdout = TRUE), " ")[[1L]][[1L]]
resolved_mapping_file <- file.path(audit_dir, "scrna_resolved_model_gene_mapping.tsv.gz")
resolved_mapping_sha256 <- strsplit(system2("sha256sum", resolved_mapping_file, stdout = TRUE), " ")[[1L]][[1L]]
summary <- list(
  status = if (gate_pass && low_coverage_groups == 0L) "PASS_SCRNA_GENE_IDENTITY_RESOLUTION" else if (gate_pass) "PASS_SCRNA_GENE_IDENTITY_RESOLUTION_WITH_LOW_COVERAGE_FLAGS" else "REVIEW_REQUIRED_SCRNA_GENE_IDENTITY_RESOLUTION",
  datasets = datasets,
  model_cell_types = length(rank_cell_type),
  minimum_W_gene_coverage = minimum_coverage,
  median_W_gene_coverage = median(coverage$WGeneCoverage),
  groups_below_90_percent = low_coverage_groups,
  low_coverage_rule = "Descriptive flag only; Cell 2021/EcoTyper zero-fills unresolved W genes and defines no 90 percent exclusion threshold",
  source_conflicts = nrow(source_conflicts),
  ensembl75_identity_sha256 = ensembl75_identity_sha256,
  org_Hs_eg_db_version = as.character(packageVersion("org.Hs.eg.db")),
  annotation_dbi_version = as.character(packageVersion("AnnotationDbi")),
  resolved_mapping_sha256 = resolved_mapping_sha256,
  duplicate_source_rows_for_one_model_gene_will_be_summed_on_raw_counts = TRUE,
  unresolved_W_genes_will_follow_official_zero_fill = TRUE,
  model_changed = FALSE,
  automatic_recovery_started = FALSE,
  release_decision = if (gate_pass && low_coverage_groups == 0L) "G4_SCRNA_INPUTS_PASS_RECOVERY_NOT_STARTED" else if (gate_pass) "G4_PASS_RECOVERY_ALLOWED_WITH_LOW_COVERAGE_FLAGS" else "G4_BLOCKED_GENE_IDENTITY_REVIEW"
)
jsonlite::write_json(summary, file.path(audit_dir, "scrna_gene_identity_resolution.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo_gene_identity.txt"))

if (!gate_pass) stop("scRNA gene identity resolution requires review; outputs saved", call. = FALSE)
cat("SCRNA_GENE_IDENTITY_RESOLUTION_PASS minimum W coverage", minimum_coverage, "\n")
