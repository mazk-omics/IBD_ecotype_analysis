#!/usr/bin/env Rscript
# 本脚本把 IBD/nonIBD cell-state phenotype validation 汇总为一张三面板工作图。
# 只用 R 成图；读取 output/26，写入 output/27_IBD_nonIBD_phenotype_validation_figures。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取汇总表和绘图依赖
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: 38_plot_IBD_nonIBD_phenotype_validation.R RUN_ROOT PANEL_ALIGNMENT_R PANEL_ALIGNMENT_COMPAT_R ALIGNMENT_AUDIT_PY", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)
panel_alignment_file <- normalizePath(args[[2L]], mustWork = TRUE)
panel_alignment_compat_file <- normalizePath(args[[3L]], mustWork = TRUE)
alignment_audit_file <- normalizePath(args[[4L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite", "ggplot2", "patchwork", "svglite", "ragg", "scales")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(ggplot2)
library(patchwork)
source(panel_alignment_file)
source(panel_alignment_compat_file)

input_dir <- file.path(RUN_ROOT, "output/26_IBD_nonIBD_phenotype_validation_summary")
pass_file <- file.path(input_dir, "IBD_nonIBD_phenotype_validation_summary_complete.json")
output_dir <- file.path(RUN_ROOT, "output/27_IBD_nonIBD_phenotype_validation_figures")
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

pass <- jsonlite::read_json(pass_file, simplifyVector = TRUE)
if (pass$status != "PASS_IBD_NONIBD_PHENOTYPE_VALIDATION_SUMMARY") stop("Run script 37 before plotting", call. = FALSE)

method_counts <- fread(file.path(input_dir, "main_broad_phenotype_method_significant_counts.tsv"))
robust <- fread(file.path(input_dir, "sample_donor_robust_state_dataset_pairs.tsv"))
support <- fread(file.path(input_dir, "robust_state_method_support_long.tsv"))
if (nrow(robust) != 4L) stop("The frozen robust state-dataset set changed", call. = FALSE)

# Step2. 固定版式、颜色和显示标签
dataset_levels <- c("GSE214695", "IBDverse_E_MTAB_16999", "SCP1884", "SCP259")
dataset_labels <- c("GSE214695", "IBDverse", "SCP1884", "SCP259")
method_levels <- c("Cell permutation Z>1.96", "Cell permutation BH<0.05", "Dominant-state Fisher", "Sample Wilcoxon", "Donor Wilcoxon", "Sample + donor concordant")
method_labels <- c("Cell permutation\nZ > 1.96", "Cell permutation\nBH < 0.05", "Dominant-state\nFisher", "Sample-level\nWilcoxon", "Donor-level\nWilcoxon", "Sample + donor\nconcordant")
support_method_levels <- c("Cell permutation Z>1.96", "Cell permutation BH<0.05", "Dominant-state Fisher", "Sample Wilcoxon", "Donor Wilcoxon")
support_method_labels <- c("Cell Z", "Cell BH", "Dominant\nFisher", "Sample", "Donor")

palette <- c(IBD = "#B44B45", nonIBD = "#3F78A8", neutral_dark = "#303030", neutral_mid = "#858585", neutral_light = "#E6E6E6", accent = "#D99B2B")

theme_validation <- theme_classic(base_size = 6.4, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(size = 5.8, colour = "black"),
    axis.title = element_text(size = 6.4, colour = "black"),
    legend.title = element_text(size = 5.8),
    legend.text = element_text(size = 5.4),
    plot.title = element_text(size = 7.1, face = "bold", colour = "black", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 5.7, colour = "#454545", margin = margin(b = 4)),
    plot.caption = element_text(size = 5.2, colour = "#555555", hjust = 0),
    plot.tag = element_text(size = 8, face = "bold"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

method_counts[, DatasetDisplay := factor(Dataset, levels = dataset_levels, labels = dataset_labels)]
method_counts[, MethodDisplay := factor(Method, levels = rev(method_levels), labels = rev(method_labels))]
method_counts[, CountLabel := paste0(SignificantStates, "/", EvaluableStates)]

state_label <- gsub("_", " ", robust$CellState, fixed = TRUE)
state_label <- gsub("Monocyte Macrophage", "Monocyte/Macrophage", state_label, fixed = TRUE)
state_label <- gsub("CD8 T", "CD8 T", state_label, fixed = TRUE)
dataset_short <- dataset_labels[match(robust$Dataset, dataset_levels)]
robust[, PairLabel := paste0(dataset_short, " | ", state_label)]
robust[, PairLabel := factor(PairLabel, levels = rev(PairLabel))]

effect_long <- rbindlist(list(
  robust[, .(Dataset, CellState, PairLabel, RobustEnrichedLevel, Unit = "Sample", DeltaPercentagePoints = 100 * DeltaMedian_Sample)],
  robust[, .(Dataset, CellState, PairLabel, RobustEnrichedLevel, Unit = "Donor", DeltaPercentagePoints = 100 * DeltaMedian_Donor)]
))
effect_long[, Unit := factor(Unit, levels = c("Sample", "Donor"))]

support <- merge(support, unique(robust[, .(Dataset, CellState, PairLabel)]), by = c("Dataset", "CellState"), all.x = TRUE, sort = FALSE)
support[, MethodDisplay := factor(Method, levels = support_method_levels, labels = support_method_labels)]
support[, PairLabel := factor(PairLabel, levels = levels(robust$PairLabel))]

fwrite(method_counts[, .(Dataset, DatasetDisplay, Method, SignificantStates, EvaluableStates, SignificantPercent)], file.path(source_dir, "panel_a_method_call_rates.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(effect_long, file.path(source_dir, "panel_b_sample_donor_effects.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(support, file.path(source_dir, "panel_c_method_support.tsv"), sep = "\t", quote = FALSE, na = "")

# Step3. Panel a：不同统计单位下的显著比例
p_a <- ggplot(method_counts, aes(x = DatasetDisplay, y = MethodDisplay, fill = SignificantPercent)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = CountLabel, colour = SignificantPercent >= 28), size = 2.0, show.legend = FALSE) +
  scale_fill_gradientn(colours = c("#F4F4F4", "#B7D2E5", "#5B93BB", "#214E70"), limits = c(0, 50), oob = scales::squish, name = "Significant\n(% evaluable)") +
  scale_colour_manual(values = c(`FALSE` = "#303030", `TRUE` = "white")) +
  labs(x = NULL, y = NULL, title = "Cell-level tests yield more phenotype calls", subtitle = "Labels show significant/evaluable states within each dataset") +
  theme_validation +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "right", legend.key.height = grid::unit(11, "mm"), axis.text.x = element_text(face = "bold"))

# Step4. Panel b：sample 与 donor 均显著的效应量
p_b <- ggplot(effect_long, aes(x = DeltaPercentagePoints, y = PairLabel, colour = RobustEnrichedLevel, shape = Unit)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", colour = palette[["neutral_mid"]]) +
  geom_point(size = 2.5, stroke = 0.55, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = palette[c("IBD", "nonIBD")], name = "Higher in") +
  scale_shape_manual(values = c(Sample = 16, Donor = 15), name = "Analysis unit") +
  scale_x_continuous(breaks = seq(-20, 20, 10), limits = c(-24, 24), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Median state-fraction difference (IBD - nonIBD, percentage points)", y = NULL, title = "Sample- and donor-consistent effects", subtitle = "Four effects pass BH FDR < 0.05 in both units") +
  theme_validation +
  theme(legend.position = "top", legend.box = "horizontal", axis.text.y = element_text(size = 5.6))

# Step5. Panel c：稳健信号的正交方法支持
p_c <- ggplot(support, aes(x = MethodDisplay, y = PairLabel)) +
  geom_point(aes(shape = SupportClass, fill = RobustEnrichedLevel), size = 3.0, stroke = 0.65, colour = palette[["neutral_dark"]]) +
  scale_shape_manual(values = c("Same direction" = 21, "Opposite direction" = 24, "Not significant" = 1), guide = "none") +
  scale_fill_manual(values = palette[c("IBD", "nonIBD")], guide = "none") +
  labs(x = NULL, y = NULL, title = "Orthogonal method support", subtitle = "Cell tests: 4/4; dominant-state Fisher: 2/4") +
  theme_validation +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none", axis.text.x = element_text(size = 5.2, lineheight = 0.9), axis.text.y = element_text(size = 5.6))

# Step6. 组合并执行 panel alignment gate
figure_width_mm <- 183
figure_height_mm <- 155
figure_width_in <- figure_width_mm / 25.4
figure_height_in <- figure_height_mm / 25.4
figure_base <- file.path(output_dir, "IBD_nonIBD_cell_state_phenotype_validation")

figure_design <- "
AA
BC
"
figure <- p_a + p_b + p_c +
  plot_layout(design = figure_design, heights = c(0.95, 1.15)) +
  plot_annotation(tag_levels = "a", caption = "Primary threshold: >=20 parent-type cells per biological sample. Tests are stratified by dataset; shared nonIBD donors in SCP1884/SCP259 are not treated as independent cross-cohort replication.") &
  theme(plot.tag = element_text(size = 8, face = "bold"), plot.caption = element_text(size = 5.2, hjust = 0, colour = "#555555"))

require_patchwork_panel_alignment(
  figure,
  manifest_path = paste0(figure_base, ".alignment-layout.json"),
  report_path = paste0(figure_base, ".alignment.json"),
  overlay_svg = paste0(figure_base, ".alignment.svg"),
  width_in = figure_width_in,
  height_in = figure_height_in,
  panel_ids = c("a", "b", "c"),
  row_groups = list(c("b", "c")),
  audit_script = alignment_audit_file,
  python = Sys.which("python3"),
  tolerance_pt = 1.5,
  gutter_tolerance_pt = 1.5,
  strict = TRUE
)

# Step7. 导出矢量图、印刷图和预览图
svglite::svglite(paste0(figure_base, ".svg"), width = figure_width_in, height = figure_height_in)
print(figure)
dev.off()

grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = figure_width_in, height = figure_height_in, family = "sans")
print(figure)
dev.off()

ragg::agg_tiff(paste0(figure_base, ".tiff"), width = figure_width_in, height = figure_height_in, units = "in", res = 600, compression = "lzw")
print(figure)
dev.off()

ragg::agg_png(paste0(figure_base, ".png"), width = figure_width_in, height = figure_height_in, units = "in", res = 300)
print(figure)
dev.off()

ragg::agg_png(paste0(figure_base, ".preview.png"), width = figure_width_in, height = figure_height_in, units = "in", res = 150)
print(figure)
dev.off()

completion <- list(
  status = "PASS_IBD_NONIBD_PHENOTYPE_VALIDATION_FIGURE",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  backend = "R_GGPLOT2_PATCHWORK", figure_width_mm = figure_width_mm, figure_height_mm = figure_height_mm,
  panels = c("method_call_rates", "sample_donor_effect_sizes", "orthogonal_method_support"), robust_pairs = nrow(robust), source_data_files = 3L
)
jsonlite::write_json(completion, file.path(output_dir, "IBD_nonIBD_phenotype_validation_figure_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_IBD_NONIBD_PHENOTYPE_VALIDATION_FIGURE\n")
