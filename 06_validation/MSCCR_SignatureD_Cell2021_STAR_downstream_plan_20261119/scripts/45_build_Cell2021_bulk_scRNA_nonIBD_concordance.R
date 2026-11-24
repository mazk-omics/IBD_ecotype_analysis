#!/usr/bin/env Rscript
# 本脚本逐数据集比较 GSE57945 bulk 与四个 scRNA 队列的 nonIBD state enrichment。
# 两侧均沿用 Cell 2021 的 Z > 1.96 判定；一致性按 2x2 表计算，并用 Fisher 精确检验评估关联。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取 bulk 和四个 scRNA 正式结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 45_build_Cell2021_bulk_scRNA_nonIBD_concordance.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

bulk_dir <- file.path(RUN_ROOT, "output/32_GSE57945_bulk_known_phenotype_Cell2021")
bulk_file <- file.path(bulk_dir, "GSE57945_bulk_known_phenotype_enrichment_Cell2021.tsv")
bulk_pass_file <- file.path(bulk_dir, "GSE57945_bulk_known_phenotype_Cell2021_complete.json")
scrna_root <- file.path(RUN_ROOT, "output/23_IBD_nonIBD_known_phenotype_enrichment_Cell2021")
datasets <- c("GSE214695", "IBDverse_E_MTAB_16999", "SCP1884", "SCP259")
scrna_files <- file.path(scrna_root, datasets, "state_known_phenotype_enrichment_Cell2021.tsv")
output_dir <- file.path(RUN_ROOT, "output/33_Cell2021_bulk_scRNA_nonIBD_concordance")

if (!file.exists(bulk_file) || !file.exists(bulk_pass_file) || any(!file.exists(scrna_files))) stop("Bulk or scRNA phenotype results are incomplete", call. = FALSE)
bulk_pass <- jsonlite::read_json(bulk_pass_file, simplifyVector = TRUE)
if (bulk_pass$status != "PASS_GSE57945_BULK_KNOWN_PHENOTYPE_ENRICHMENT_CELL2021") stop("Bulk phenotype enrichment did not pass", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bulk <- fread(bulk_file)
bulk <- bulk[Axis == "BroadPhenotype" & tolower(gsub("-", "", Level)) == "nonibd"]
if (anyDuplicated(bulk$CellState)) stop("Duplicate bulk state rows", call. = FALSE)
bulk <- bulk[, .(
  CellType, State, CellState, BulkZ = Z, BulkP = P_one_sided, BulkStatus = Status,
  BulkSignificant = PositiveZgt1.96, BulkObserved = ObservedAssignedSamples
)]

# Step2. 逐 scRNA 队列冻结可比较 state 和显著性调用
comparison_parts <- vector("list", length(datasets))
summary_parts <- vector("list", length(datasets))

for (dataset_index in seq_along(datasets)) {
  dataset <- datasets[[dataset_index]]
  scrna <- fread(scrna_files[[dataset_index]])
  scrna <- scrna[Axis == "BroadPhenotype" & tolower(gsub("-", "", Level)) == "nonibd"]
  if (anyDuplicated(scrna$CellState)) stop("Duplicate scRNA state rows for ", dataset, call. = FALSE)
  scrna <- scrna[, .(
    CellState, ScRNAZ = Z, ScRNAP = P_one_sided, ScRNAStatus = Status,
    ScRNASignificant = PositiveZgt1.96, ScRNAObserved = ObservedAssignedCells
  )]

  comparison <- merge(bulk, scrna, by = "CellState", all = FALSE, sort = FALSE)
  comparison <- comparison[BulkStatus == "EVALUABLE" & ScRNAStatus == "EVALUABLE" & is.finite(BulkZ) & is.finite(ScRNAZ)]
  if (!nrow(comparison)) stop("No evaluable states for ", dataset, call. = FALSE)
  comparison[, `:=`(
    ScRNADataset = dataset,
    BulkScore = pmin(-log10(pmax(BulkP, 1e-300)), 2),
    ScRNAScore = pmin(-log10(pmax(ScRNAP, 1e-300)), 2),
    ConcordantCall = BulkSignificant == ScRNASignificant,
    CallClass = fifelse(BulkSignificant & ScRNASignificant, "Both significant",
      fifelse(BulkSignificant, "Bulk only", fifelse(ScRNASignificant, "scRNA only", "Both not significant")))
  )]
  setorder(comparison, -BulkScore, -ScRNAScore, CellType, State)
  comparison[, DisplayOrder := seq_len(.N)]

  contingency <- table(
    factor(comparison$BulkSignificant, levels = c(FALSE, TRUE)),
    factor(comparison$ScRNASignificant, levels = c(FALSE, TRUE))
  )
  fisher_two_sided <- fisher.test(contingency, alternative = "two.sided")
  fisher_greater <- fisher.test(contingency, alternative = "greater")
  both_significant <- contingency[2L, 2L]
  bulk_only <- contingency[2L, 1L]
  scrna_only <- contingency[1L, 2L]
  both_not_significant <- contingency[1L, 1L]
  positive_union <- both_significant + bulk_only + scrna_only

  summary_parts[[dataset_index]] <- data.table(
    BulkDataset = "GSE57945", ScRNADataset = dataset, EvaluableStates = nrow(comparison),
    BothSignificant = both_significant, BulkOnly = bulk_only, ScRNAOnly = scrna_only, BothNotSignificant = both_not_significant,
    ConcordantStates = both_significant + both_not_significant,
    ConcordancePercent = 100 * (both_significant + both_not_significant) / nrow(comparison),
    PositiveSetJaccardPercent = fifelse(positive_union > 0, 100 * both_significant / positive_union, NA_real_),
    BulkPositiveReplicationPercent = fifelse(both_significant + bulk_only > 0, 100 * both_significant / (both_significant + bulk_only), NA_real_),
    FisherTwoSidedP = fisher_two_sided$p.value, FisherGreaterP = fisher_greater$p.value,
    FisherOddsRatio = unname(fisher_two_sided$estimate), FisherConfLow = fisher_two_sided$conf.int[[1L]], FisherConfHigh = fisher_two_sided$conf.int[[2L]],
    PrimaryDefinition = "overall agreement of Z>1.96 calls across evaluable states; Fisher exact test on 2x2 calls"
  )
  comparison_parts[[dataset_index]] <- comparison
  cat("CONCORDANCE_READY", dataset, "states", nrow(comparison), "agreement", round(summary_parts[[dataset_index]]$ConcordancePercent, 1), "\n")
}

# Step3. 保存逐 state 证据和 2x2 统计
comparisons <- rbindlist(comparison_parts, use.names = TRUE, fill = TRUE)
summary <- rbindlist(summary_parts, use.names = TRUE, fill = TRUE)
setcolorder(comparisons, c(
  "ScRNADataset", "DisplayOrder", "CellType", "State", "CellState", "BulkZ", "BulkP", "BulkScore", "BulkSignificant", "BulkObserved",
  "ScRNAZ", "ScRNAP", "ScRNAScore", "ScRNASignificant", "ScRNAObserved", "ConcordantCall", "CallClass", "BulkStatus", "ScRNAStatus"
))
fwrite(comparisons, file.path(output_dir, "Cell2021_bulk_scRNA_nonIBD_state_concordance.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(summary, file.path(output_dir, "Cell2021_bulk_scRNA_nonIBD_concordance_summary.tsv"), sep = "\t", quote = FALSE, na = "")

completion <- list(
  status = "PASS_CELL2021_BULK_SCRNA_NONIBD_CONCORDANCE", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  bulk_dataset = "GSE57945", scrna_datasets = datasets, comparisons = length(datasets),
  state_call_threshold = "Z > 1.96", fisher_test = "two-sided primary; greater sensitivity",
  concordance_definition = "overall binary agreement among states evaluable in both datasets"
)
jsonlite::write_json(completion, file.path(output_dir, "Cell2021_bulk_scRNA_nonIBD_concordance_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(summary)
cat("PASS_CELL2021_BULK_SCRNA_NONIBD_CONCORDANCE\n")
