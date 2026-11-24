#!/usr/bin/env Rscript
# 本脚本输出 Cell 2021 式表型富集图，以及 Figure S2G/S2H 对应的 dominant-state representation 图。
# 只使用 output/23 和 output/24 的冻结结果；图件写入 output/28_Cell2021_phenotype_and_S2GH_figures。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取正式结果和绘图门禁
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: 39_plot_Cell2021_phenotype_enrichment_and_S2GH.R RUN_ROOT PANEL_ALIGNMENT_R PANEL_ALIGNMENT_COMPAT_R ALIGNMENT_AUDIT_PY OUTPUT_DIR", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)
panel_alignment_file <- normalizePath(args[[2L]], mustWork = TRUE)
panel_alignment_compat_file <- normalizePath(args[[3L]], mustWork = TRUE)
alignment_audit_file <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- args[[5L]]

required_pkgs <- c("data.table", "jsonlite", "ggplot2", "patchwork", "svglite", "ragg", "scales")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(ggplot2)
library(patchwork)
source(panel_alignment_file)
source(panel_alignment_compat_file)

permutation_dir <- file.path(RUN_ROOT, "output/23_IBD_nonIBD_known_phenotype_enrichment_Cell2021")
representation_dir <- file.path(RUN_ROOT, "output/24_IBD_nonIBD_state_representation")
permutation_files <- list.files(permutation_dir, pattern = "^state_known_phenotype_enrichment_Cell2021\\.tsv$", recursive = TRUE, full.names = TRUE)
representation_file <- file.path(representation_dir, "sample_dominant_state_representation_fisher.tsv")
representation_pass_file <- file.path(representation_dir, "cell_state_representation_complete.json")
if (length(permutation_files) != 4L || !file.exists(representation_file) || !file.exists(representation_pass_file)) stop("Formal phenotype or representation outputs are incomplete", call. = FALSE)

permutation <- rbindlist(lapply(permutation_files, fread), use.names = TRUE, fill = TRUE)
representation <- fread(representation_file)
representation_pass <- jsonlite::read_json(representation_pass_file, simplifyVector = TRUE)
if (representation_pass$status != "PASS_IBD_NONIBD_CELL_STATE_REPRESENTATION") stop("Representation audit did not pass", call. = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

axes <- c("BroadPhenotype", "Disease", "Inflammation", "Region3")
dataset_levels <- c("GSE214695", "IBDverse_E_MTAB_16999", "SCP1884", "SCP259")
dataset_labels <- c("GSE214695", "IBDverse", "SCP1884", "SCP259")
state_order <- unique(permutation[order(CellType, State), CellState])
state_display <- gsub("_", " ", state_order, fixed = TRUE)
state_display <- gsub("Monocyte Macrophage", "Monocyte/Macrophage", state_display, fixed = TRUE)

cat("Permutation rows:", nrow(permutation), "\n")
cat("Representation rows:", nrow(representation), "\n")
cat("Formal states:", length(state_order), "\n")

# Step2. 固定版式和导出函数
palette <- c(blue = "#3F78A8", red = "#B44B45", teal = "#2A8C82", neutral_dark = "#303030", neutral_mid = "#8A8A8A", neutral_light = "#E6E6E6")

theme_heatmap <- theme_minimal(base_size = 6.3, base_family = "sans") +
  theme(
    panel.grid = element_blank(), axis.title = element_text(size = 6.3, colour = "black"),
    axis.text.x = element_text(size = 5.1, colour = "black", angle = 55, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5.0, colour = "black"),
    plot.title = element_text(size = 7.3, face = "bold", colour = "black", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 5.8, colour = "#454545", margin = margin(b = 4)),
    plot.caption = element_text(size = 5.2, colour = "#555555", hjust = 0),
    plot.tag = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 5.8), legend.text = element_text(size = 5.3),
    legend.position = "bottom", legend.key.width = grid::unit(18, "mm"),
    panel.border = element_rect(colour = "#555555", fill = NA, linewidth = 0.25),
    plot.margin = margin(5, 5, 5, 5)
  )

save_figure <- function(plot_object, file_stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(file_stem, ".svg"), width = width_in, height = height_in)
  print(plot_object)
  dev.off()
  grDevices::cairo_pdf(paste0(file_stem, ".pdf"), width = width_in, height = height_in, family = "sans")
  print(plot_object)
  dev.off()
  ragg::agg_tiff(paste0(file_stem, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
  print(plot_object)
  dev.off()
  ragg::agg_png(paste0(file_stem, ".png"), width = width_in, height = height_in, units = "in", res = 300)
  print(plot_object)
  dev.off()
}

# Step3. 输出原文式 known phenotype enrichment Z 热图
permutation[, DatasetDisplay := factor(Dataset, levels = dataset_levels, labels = dataset_labels)]
permutation[, CellStateDisplay := factor(CellState, levels = rev(state_order), labels = rev(state_display))]
permutation[, DisplayZ := pmax(pmin(Z, 4), -4)]
fwrite(permutation, file.path(source_dir, "known_phenotype_enrichment_Cell2021_all.tsv"), sep = "\t", quote = FALSE, na = "")

permutation_summary <- permutation[, .(
  EvaluableStates = uniqueN(CellState[Status == "EVALUABLE"]),
  PositiveZgt1.96States = uniqueN(CellState[Status == "EVALUABLE" & PositiveZgt1.96 == TRUE]),
  SignificantBHStates = uniqueN(CellState[Status == "EVALUABLE" & Significant_BH_sensitivity == TRUE])
), by = .(Dataset, Axis)]
permutation_summary[, `:=`(
  PositivePercent = 100 * PositiveZgt1.96States / EvaluableStates,
  BHPercent = 100 * SignificantBHStates / EvaluableStates
)]
fwrite(permutation_summary, file.path(source_dir, "known_phenotype_enrichment_summary.tsv"), sep = "\t", quote = FALSE, na = "")

for (axis_name in axes) {
  plot_data <- copy(permutation[Axis == axis_name])
  plot_data[, DatasetShort := dataset_labels[match(Dataset, dataset_levels)]]
  column_order <- unique(plot_data[order(match(Dataset, dataset_levels), Level), paste(DatasetShort, Level, sep = " | ")])
  plot_data[, DatasetLevel := factor(paste(DatasetShort, Level, sep = " | "), levels = column_order)]

  phenotype_plot <- ggplot(plot_data, aes(x = DatasetLevel, y = CellStateDisplay, fill = DisplayZ)) +
    geom_tile(colour = "white", linewidth = 0.10) +
    geom_point(data = plot_data[PositiveZgt1.96 == TRUE & Significant_BH_sensitivity != TRUE], shape = 4, size = 0.65, stroke = 0.25, colour = "black") +
    geom_point(data = plot_data[Significant_BH_sensitivity == TRUE], shape = 8, size = 0.82, stroke = 0.20, colour = "black") +
    scale_fill_gradient2(low = "#3F78A8", mid = "#F7F7F7", high = "#B44B45", midpoint = 0, limits = c(-4, 4), na.value = "#D9D9D9", name = "Permutation Z\n(capped at +/-4)") +
    labs(
      title = paste0("Known phenotype enrichment: ", axis_name),
      subtitle = "Cell 2021 state-label permutation within each parent type and dataset",
      x = "Dataset | phenotype level", y = "Cell state",
      caption = "Cross: nominal Z > 1.96; star: BH FDR < 0.05 sensitivity; grey: not evaluable. Permutations = 1,000."
    ) +
    theme_heatmap

  save_figure(phenotype_plot, file.path(output_dir, paste0("known_phenotype_enrichment_", axis_name)), 183, 285)
}

# Step4. 准备 Figure S2G/S2H 对应的 representation 数据
representation_main <- representation[MinParentCells == 20L & Scope == "WITHIN_DATASET" & Axis %in% axes]
representation_main[, RepresentationPercent := 100 * SamplesTargetLevelState / (SamplesTargetLevelState + SamplesTargetLevelOtherStates)]
representation_main[, EnrichmentScore := fifelse(Status == "EVALUABLE" & is.finite(FDR_BH_across_cell_types), pmin(-log10(pmax(FDR_BH_across_cell_types, 1e-300)), 10), NA_real_)]
representation_main[, DatasetShort := dataset_labels[match(DatasetScope, dataset_levels)]]
representation_main[, CellStateDisplay := factor(CellState, levels = rev(state_order), labels = rev(state_display))]
fwrite(representation_main, file.path(source_dir, "dominant_state_representation_within_dataset_min20.tsv"), sep = "\t", quote = FALSE, na = "")

representation_summary <- representation_main[, .(
  EvaluableStates = uniqueN(CellState[Status == "EVALUABLE"]),
  SignificantFDRStates = uniqueN(CellState[Status == "EVALUABLE" & SignificantEnrichment == TRUE])
), by = .(Dataset = DatasetScope, Axis)]
representation_summary[, SignificantPercent := 100 * SignificantFDRStates / EvaluableStates]
fwrite(representation_summary, file.path(source_dir, "dominant_state_representation_summary.tsv"), sep = "\t", quote = FALSE, na = "")

# Step5. 逐表型轴输出 S2G 比例与 S2H Fisher-BH
for (axis_name in axes) {
  plot_data <- copy(representation_main[Axis == axis_name])
  column_order <- unique(plot_data[order(match(DatasetScope, dataset_levels), Level), paste(DatasetShort, Level, sep = " | ")])
  plot_data[, DatasetLevel := factor(paste(DatasetShort, Level, sep = " | "), levels = column_order)]

  p_g <- ggplot(plot_data, aes(x = DatasetLevel, y = CellStateDisplay, fill = RepresentationPercent)) +
    geom_tile(colour = "white", linewidth = 0.10) +
    scale_fill_gradientn(colours = c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"), limits = c(0, 100), na.value = "#D9D9D9", name = "Dominant samples\n(%)") +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(52, "mm"), barheight = unit(2.5, "mm"))) +
    labs(
      tag = "g", title = paste0(axis_name, ": dominant-state representation"),
      subtitle = "Percentage among eligible biological samples within each parent type",
      x = "Dataset | phenotype level", y = "Cell state"
    ) +
    theme_heatmap

  p_h <- ggplot(plot_data, aes(x = DatasetLevel, y = CellStateDisplay, fill = EnrichmentScore)) +
    geom_tile(colour = "white", linewidth = 0.10) +
    geom_point(data = plot_data[SignificantEnrichment == TRUE], shape = 8, size = 0.82, stroke = 0.20, colour = "black") +
    scale_fill_gradientn(colours = c("#F7FCF5", "#C7E9C0", "#74C476", "#238B45", "#005A32"), limits = c(0, 10), na.value = "#D9D9D9", name = "-log10(BH FDR)\n(capped at 10)") +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(52, "mm"), barheight = unit(2.5, "mm"))) +
    labs(
      tag = "h", title = paste0(axis_name, ": over-representation"),
      subtitle = "One-sided Fisher exact test on dominant-state assignments",
      x = "Dataset | phenotype level", y = "Cell state"
    ) +
    theme_heatmap

  figure_stem <- file.path(output_dir, paste0("FigS2GH_state_representation_", axis_name))
  figure_pair <- p_g + p_h + plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(caption = "Primary threshold: >=20 parent-type cells per biological sample. Star: BH FDR < 0.05; grey in panel h: Fisher margin not evaluable.") &
    theme(plot.caption = element_text(size = 5.2, colour = "#555555", hjust = 0))

  require_patchwork_panel_alignment(
    figure_pair,
    manifest_path = paste0(figure_stem, ".alignment-layout.json"),
    report_path = paste0(figure_stem, ".alignment.json"),
    overlay_svg = paste0(figure_stem, ".alignment.svg"),
    width_in = 183 / 25.4,
    height_in = 285 / 25.4,
    panel_ids = c("g", "h"),
    row_groups = list(c("g", "h")),
    audit_script = alignment_audit_file,
    python = Sys.which("python3"),
    tolerance_pt = 1.5,
    gutter_tolerance_pt = 1.5,
    strict = TRUE
  )

  save_figure(figure_pair, figure_stem, 183, 285)
}

# Step6. 输出紧凑结果总览
permutation_raw <- permutation_summary[, .(
  Dataset, Axis, Method = "Permutation Z > 1.96",
  SignificantStates = PositiveZgt1.96States, EvaluableStates,
  SignificantPercent = PositivePercent
)]
permutation_bh <- permutation_summary[, .(
  Dataset, Axis, Method = "Permutation BH FDR < 0.05",
  SignificantStates = SignificantBHStates, EvaluableStates,
  SignificantPercent = BHPercent
)]
representation_bh <- representation_summary[, .(
  Dataset, Axis, Method = "Dominant-state Fisher BH",
  SignificantStates = SignificantFDRStates, EvaluableStates,
  SignificantPercent
)]
summary_data <- rbindlist(list(permutation_raw, permutation_bh, representation_bh), use.names = TRUE)
summary_data[EvaluableStates == 0, SignificantPercent := NA_real_]
summary_data[, CountLabel := fifelse(EvaluableStates == 0, "NE", paste0(SignificantStates, "/", EvaluableStates))]
summary_data[, DatasetDisplay := dataset_labels[match(Dataset, dataset_levels)]]
method_levels <- c("Permutation Z > 1.96", "Permutation BH FDR < 0.05", "Dominant-state Fisher BH")
summary_data[, Method := factor(Method, levels = method_levels)]
summary_data[, Axis := factor(Axis, levels = axes)]
row_levels <- unlist(lapply(dataset_labels, function(dataset_name) paste(dataset_name, method_levels, sep = " | ")))
summary_data[, RowLabel := factor(paste(DatasetDisplay, Method, sep = " | "), levels = rev(row_levels))]
fwrite(summary_data, file.path(source_dir, "phenotype_and_representation_summary.tsv"), sep = "\t", quote = FALSE, na = "")

summary_stem <- file.path(output_dir, "phenotype_and_representation_summary")
summary_figure <- ggplot(summary_data, aes(x = Axis, y = RowLabel, fill = SignificantPercent)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = CountLabel), size = 2.0, colour = "black") +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    limits = c(0, 60), oob = scales::squish, na.value = "#D9D9D9",
    name = "Significant\nstates (%)"
  ) +
  labs(
    title = "Phenotype enrichment and dominant-state representation",
    subtitle = "Labels show significant/evaluable states; NE = not evaluable",
    x = NULL, y = NULL,
    caption = "Permutation: positive Z > 1.96 or within-dataset BH FDR < 0.05. Representation: one-sided Fisher test with within-dataset BH correction."
  ) +
  theme_minimal(base_size = 7.2, base_family = "sans") +
  theme(
    panel.grid = element_blank(), axis.text.x = element_text(size = 6.6, angle = 25, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 5.8, colour = "black"), plot.title = element_text(size = 8.4, face = "bold"),
    plot.subtitle = element_text(size = 6.5, colour = "#454545"), plot.caption = element_text(size = 5.4, colour = "#555555", hjust = 0),
    legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8), legend.position = "right",
    plot.margin = margin(6, 8, 6, 6)
  )
save_figure(summary_figure, summary_stem, 183, 120)

# Step7. 完成记录
completion <- list(
  status = "PASS_CELL2021_PHENOTYPE_AND_S2GH_FIGURES",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  backend = "R_GGPLOT2_PATCHWORK", formal_states = length(state_order), datasets = uniqueN(permutation$Dataset),
  phenotype_axes = axes, phenotype_heatmaps = length(axes), representation_pairs = length(axes),
  main_threshold_parent_cells = 20L, permutations = unique(permutation$Permutations), source_data_files = 5L
)
jsonlite::write_json(completion, file.path(output_dir, "Cell2021_phenotype_and_S2GH_figures_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(permutation_summary[order(Dataset, Axis)])
print(representation_summary[order(Dataset, Axis)])
print(completion)
cat("PASS_CELL2021_PHENOTYPE_AND_S2GH_FIGURES\n")
