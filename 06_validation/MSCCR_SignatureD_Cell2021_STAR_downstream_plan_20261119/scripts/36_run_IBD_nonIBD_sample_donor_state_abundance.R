#!/usr/bin/env Rscript
# 本脚本以 sample 和 donor 两种统计单位比较 IBD/nonIBD cell-state proportion，不使用逐细胞伪重复检验。
# 20 个父类型细胞是主分析门槛，10 个是敏感性分析；结果写入 output/25_IBD_nonIBD_sample_donor_abundance。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取统一逐细胞分配
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 36_run_IBD_nonIBD_sample_donor_state_abundance.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)

input_dir <- file.path(RUN_ROOT, "output/22_IBD_nonIBD_phenotype_input")
input_file <- file.path(input_dir, "IBD_nonIBD_formal_cell_assignments.tsv.gz")
pass_file <- file.path(input_dir, "IBD_nonIBD_phenotype_input_complete.json")
output_dir <- file.path(RUN_ROOT, "output/25_IBD_nonIBD_sample_donor_abundance")
if (!file.exists(input_file) || !file.exists(pass_file)) stop("Run script 33 before state abundance analysis", call. = FALSE)
pass <- jsonlite::read_json(pass_file, simplifyVector = TRUE)
if (pass$status != "PASS_IBD_NONIBD_PHENOTYPE_INPUT_AUDIT") stop("Unified phenotype input audit did not pass", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cells <- fread(input_file, colClasses = "character")
axes <- c("BroadPhenotype", "Disease", "Inflammation", "Region3")
thresholds <- c(20L, 10L)
min_units <- 3L
required_columns <- c("AnalysisGroup", "Dataset", "CellType", "ID", "State", "DonorID", "BiologicalSampleID", axes, "FormalStateUniverse")
if (!all(required_columns %in% names(cells))) stop("Unified phenotype input columns changed", call. = FALSE)
if (anyDuplicated(cells[, .(Dataset, CellType, ID)])) stop("Duplicate cells in unified phenotype input", call. = FALSE)

state_universe <- unique(cells[, .(CellType, FormalStateUniverse)])
if (state_universe[, anyDuplicated(CellType)] != 0L) stop("Formal state universe differs across datasets", call. = FALSE)
state_universe <- state_universe[, .(State = strsplit(FormalStateUniverse, ";", fixed = TRUE)[[1L]]), by = CellType]

# Step2. 建立 sample-level state proportion
sample_keys <- c("AnalysisGroup", "Dataset", "CellType", "BiologicalSampleID")
sample_meta <- unique(cells[, c(sample_keys, "DonorID", axes), with = FALSE])
if (sample_meta[, anyDuplicated(paste(AnalysisGroup, Dataset, CellType, BiologicalSampleID))] != 0L) stop("A biological sample has conflicting phenotype labels", call. = FALSE)

sample_totals <- cells[, .(
  InputCells = .N, AssignedCells = sum(State != "Unassigned"), UnassignedCells = sum(State == "Unassigned")
), by = sample_keys]
sample_counts <- cells[State != "Unassigned", .(StateCells = .N), by = c(sample_keys, "State")]
sample_grid <- merge(sample_totals, state_universe, by = "CellType", allow.cartesian = TRUE, sort = FALSE)
sample_grid <- merge(sample_grid, sample_counts, by = c(sample_keys, "State"), all.x = TRUE, sort = FALSE)
sample_grid[is.na(StateCells), StateCells := 0L]
sample_grid <- merge(sample_grid, sample_meta, by = sample_keys, all.x = TRUE, sort = FALSE)
sample_grid[, `:=`(StateFractionInput = StateCells / InputCells, AssignmentRate = AssignedCells / InputCells)]
setorder(sample_grid, Dataset, CellType, BiologicalSampleID, State)

cat("Sample-state rows:", nrow(sample_grid), "\n")
print(sample_grid[, .(SampleParents = uniqueN(paste(Dataset, CellType, BiologicalSampleID)), BiologicalSamples = uniqueN(paste(Dataset, BiologicalSampleID)), Donors = uniqueN(DonorID))])

# Step3. 针对每个表型轴建立 donor-level state proportion
donor_parts <- list()

for (axis_name in axes) {
  donor_totals <- cells[!is.na(get(axis_name)) & nzchar(get(axis_name)), .(
    InputCells = .N, AssignedCells = sum(State != "Unassigned"), UnassignedCells = sum(State == "Unassigned"), BiologicalSamples = uniqueN(BiologicalSampleID)
  ), by = .(Dataset, CellType, DonorID, Level = get(axis_name))]
  donor_counts <- cells[State != "Unassigned" & !is.na(get(axis_name)) & nzchar(get(axis_name)), .(StateCells = .N), by = .(Dataset, CellType, DonorID, Level = get(axis_name), State)]
  donor_grid <- merge(donor_totals, state_universe, by = "CellType", allow.cartesian = TRUE, sort = FALSE)
  donor_grid <- merge(donor_grid, donor_counts, by = c("Dataset", "CellType", "DonorID", "Level", "State"), all.x = TRUE, sort = FALSE)
  donor_grid[is.na(StateCells), StateCells := 0L]
  donor_grid[, `:=`(Axis = axis_name, StateFractionInput = StateCells / InputCells, AssignmentRate = AssignedCells / InputCells)]
  donor_parts[[axis_name]] <- donor_grid
}

donor_grid <- rbindlist(donor_parts, use.names = TRUE, fill = TRUE)
setcolorder(donor_grid, c("Dataset", "CellType", "DonorID", "Axis", "Level", "State", "StateCells", "InputCells", "AssignedCells", "UnassignedCells", "BiologicalSamples", "StateFractionInput", "AssignmentRate"))
setorder(donor_grid, Dataset, CellType, DonorID, Axis, Level, State)
cat("Donor-state rows:", nrow(donor_grid), "\n")

# Step4. 固定比较方向并运行 sample/donor 检验
contrast_specs <- rbindlist(list(
  data.table(Axis = "BroadPhenotype", Target = "IBD", Reference = "nonIBD"),
  data.table(Axis = "Disease", Target = c("CD", "UC", "CD"), Reference = c("nonIBD", "nonIBD", "UC")),
  data.table(Axis = "Inflammation", Target = c("Inflamed", "Non-inflamed", "Inflamed"), Reference = c("NonIBD", "NonIBD", "Non-inflamed")),
  data.table(Axis = "Region3", Target = c("Colon", "Rectum", "Rectum"), Reference = c("Ileum", "Colon", "Ileum"))
))
contrast_specs[, Contrast := paste(Target, "vs", Reference)]

test_parts <- list()
summary_parts <- list()
test_index <- 0L
summary_index <- 0L

for (unit_name in c("Sample", "Donor")) {
  for (threshold in thresholds) {
    for (axis_name in axes) {
      if (unit_name == "Sample") {
        analysis_data <- sample_grid[InputCells >= threshold]
        analysis_data[, Level := get(axis_name)]
        unit_id <- "BiologicalSampleID"
      } else {
        analysis_data <- donor_grid[Axis == axis_name & InputCells >= threshold]
        unit_id <- "DonorID"
      }
      analysis_data <- analysis_data[!is.na(Level) & nzchar(Level)]

      group_summary <- analysis_data[, .(
        Units = uniqueN(get(unit_id)), MedianFraction = median(StateFractionInput), Q25Fraction = as.numeric(quantile(StateFractionInput, 0.25)),
        Q75Fraction = as.numeric(quantile(StateFractionInput, 0.75)), MeanFraction = mean(StateFractionInput),
        MedianAssignmentRate = median(AssignmentRate)
      ), by = .(Dataset, CellType, State, Level)]
      group_summary[, `:=`(Unit = unit_name, MinParentCells = threshold, Axis = axis_name)]
      summary_index <- summary_index + 1L
      summary_parts[[summary_index]] <- group_summary

      axis_contrasts <- contrast_specs[Axis == axis_name]
      for (contrast_index in seq_len(nrow(axis_contrasts))) {
        target_level <- axis_contrasts$Target[[contrast_index]]
        reference_level <- axis_contrasts$Reference[[contrast_index]]
        contrast_name <- axis_contrasts$Contrast[[contrast_index]]

        keys <- unique(analysis_data[, .(Dataset, CellType, State)])
        for (key_index in seq_len(nrow(keys))) {
          dataset_name <- keys$Dataset[[key_index]]
          cell_type <- keys$CellType[[key_index]]
          state_name <- keys$State[[key_index]]
          values <- analysis_data[Dataset == dataset_name & CellType == cell_type & State == state_name & Level %in% c(target_level, reference_level)]
          target <- values[Level == target_level, StateFractionInput]
          reference <- values[Level == reference_level, StateFractionInput]
          evaluable <- length(target) >= min_units && length(reference) >= min_units

          if (evaluable) {
            test <- suppressWarnings(wilcox.test(target, reference, alternative = "two.sided", exact = FALSE))
            rank_biserial <- mean(outer(target, reference, FUN = function(x, y) sign(x - y)))
            p_value <- test$p.value
          } else {
            rank_biserial <- NA_real_
            p_value <- NA_real_
          }

          test_index <- test_index + 1L
          test_parts[[test_index]] <- data.table(
            Unit = unit_name, MinParentCells = threshold, Dataset = dataset_name, CellType = cell_type, State = state_name,
            CellState = paste(cell_type, state_name, sep = "_"), Axis = axis_name, Contrast = contrast_name,
            Target = target_level, Reference = reference_level, TargetUnits = length(target), ReferenceUnits = length(reference),
            MedianTarget = if (length(target)) median(target) else NA_real_, MedianReference = if (length(reference)) median(reference) else NA_real_,
            DeltaMedian = if (length(target) && length(reference)) median(target) - median(reference) else NA_real_,
            RankBiserial = rank_biserial, P_two_sided = p_value,
            Status = if (evaluable) "EVALUABLE" else "NOT_EVALUABLE_LT3_UNITS_PER_GROUP"
          )
        }
      }
      cat("ABUNDANCE_TEST_COMPLETE", unit_name, threshold, axis_name, "\n")
    }
  }
}

# Step5. 保存 abundance、效应量和显著性
tests <- rbindlist(test_parts, use.names = TRUE, fill = TRUE)
tests[, FDR_BH_within_dataset_contrast := p.adjust(P_two_sided, method = "BH"), by = .(Unit, MinParentCells, Dataset, Axis, Contrast)]
tests[, SignificantFDR := is.finite(FDR_BH_within_dataset_contrast) & FDR_BH_within_dataset_contrast < 0.05]
tests[, Direction := fifelse(is.na(DeltaMedian), "NotEvaluable", fifelse(DeltaMedian > 0, "TargetHigher", fifelse(DeltaMedian < 0, "ReferenceHigher", "NoMedianDifference")))]
setorder(tests, Unit, MinParentCells, Dataset, Axis, Contrast, CellType, State)

group_summary <- rbindlist(summary_parts, use.names = TRUE, fill = TRUE)
setcolorder(group_summary, c("Unit", "MinParentCells", "Dataset", "CellType", "State", "Axis", "Level", "Units", "MedianFraction", "Q25Fraction", "Q75Fraction", "MeanFraction", "MedianAssignmentRate"))
setorder(group_summary, Unit, MinParentCells, Dataset, Axis, Level, CellType, State)

fwrite(sample_grid, file.path(output_dir, "sample_state_abundance_all.tsv.gz"), sep = "\t", quote = FALSE, na = "")
fwrite(donor_grid, file.path(output_dir, "donor_state_abundance_by_axis_all.tsv.gz"), sep = "\t", quote = FALSE, na = "")
fwrite(group_summary, file.path(output_dir, "sample_donor_state_abundance_group_summary.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(tests, file.path(output_dir, "sample_donor_state_abundance_tests.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(tests[SignificantFDR == TRUE], file.path(output_dir, "sample_donor_state_abundance_significant.tsv"), sep = "\t", quote = FALSE, na = "")

completion <- list(
  status = "PASS_IBD_NONIBD_SAMPLE_DONOR_STATE_ABUNDANCE", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = uniqueN(cells$Dataset), input_cells = nrow(cells), sample_state_rows = nrow(sample_grid), donor_state_rows = nrow(donor_grid),
  thresholds = thresholds, minimum_units_per_group = min_units, evaluable_tests = sum(tests$Status == "EVALUABLE"),
  significant_fdr_tests = sum(tests$SignificantFDR), method = "WITHIN_DATASET_SAMPLE_AND_DONOR_WILCOXON_EFFECT_SIZE"
)
jsonlite::write_json(completion, file.path(output_dir, "sample_donor_state_abundance_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_IBD_NONIBD_SAMPLE_DONOR_STATE_ABUNDANCE\n")
