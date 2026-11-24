#!/usr/bin/env Rscript
# 本脚本按 Cell 2021 Fig. S2G/S2H 原版式输出 state representation 箱线图。
# 输入为已完成的样本 dominant-state Fisher-BH 结果；不重新计算或筛选表型。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读入正式 Fisher-BH 结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: 40_plot_Cell2021_S2GH_boxplots.R REPRESENTATION_TSV OUTPUT_DIR PANEL_ALIGNMENT_R PANEL_ALIGNMENT_COMPAT_R ALIGNMENT_AUDIT_PY", call. = FALSE)
representation_file <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- args[[2L]]
panel_alignment_file <- normalizePath(args[[3L]], mustWork = TRUE)
panel_alignment_compat_file <- normalizePath(args[[4L]], mustWork = TRUE)
alignment_audit_file <- normalizePath(args[[5L]], mustWork = TRUE)

required_pkgs <- c("data.table", "ggplot2", "patchwork", "jsonlite", "svglite", "ragg")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(ggplot2)
library(patchwork)
source(panel_alignment_file)
source(panel_alignment_compat_file)

representation <- fread(representation_file)
representation <- representation[MinParentCells == 20L & Scope == "WITHIN_DATASET"]
if (!nrow(representation)) stop("No within-dataset min20 representation rows", call. = FALSE)

axes <- c("BroadPhenotype", "Disease", "Inflammation", "Region3")
representation <- representation[Axis %in% axes & Status == "EVALUABLE"]
representation[, OddsRatio := as.numeric(OddsRatio)]
representation[, P_one_sided := as.numeric(P_one_sided)]
representation[, FDR_BH_previous := as.numeric(FDR_BH_across_cell_types)]
representation[, FDR_BH_across_cell_types := p.adjust(P_one_sided, method = "BH"), by = .(MinParentCells, Scope, DatasetScope, Axis, Level)]
representation <- representation[is.finite(OddsRatio) & OddsRatio > 1 & is.finite(FDR_BH_across_cell_types) & FDR_BH_across_cell_types > 0]
representation[, Enrichment := pmin(-log10(FDR_BH_across_cell_types), 10)]
representation[, Context := paste(DatasetShort, Level, sep = " | ")]
if (!nrow(representation)) stop("No positive enrichment rows", call. = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(representation, file.path(source_dir, "positive_state_representation_enrichment.tsv"), sep = "\t", quote = FALSE, na = "")

# Step2. 固定原文颜色和亲本谱系
state_levels <- sort(unique(representation$State))
state_colours <- c(S01 = "#C2185B", S02 = "#4F86A8", S03 = "#5A963E", S04 = "#563B91", S05 = "#D4631B", S06 = "#E6B51C", S07 = "#7A3E1C", S08 = "#CB5BDC", S09 = "#858585")
state_colours <- state_colours[state_levels]

cell_type_levels <- sort(unique(representation$CellType))
cell_type_colours <- c(
  B = "#C8961A", Plasma = "#C02A74", CD8_T = "#5B5AA6", CD4_T = "#A64B16", NK = "#3E855C",
  ILC = "#7867A8", Monocyte_Macrophage = "#B39ACD", DC = "#E5A54C", Mast = "#D57D78", Fibroblast = "#94B6D6",
  Endothelial = "#777777", Glial = "#7AB7A7", Pericyte = "#77A84D", Colonic_absorptive = "#C98EBF",
  Ileal_absorptive = "#D85E9F", Enteroendocrine = "#A46AA6", Goblet = "#5B9DC0", Paneth = "#E8B65C",
  Stem_TA_progenitor = "#8DB85C", Tuft = "#6E87B9"
)
cell_type_colours <- cell_type_colours[cell_type_levels]

lineage_map <- data.table(
  CellType = cell_type_levels,
  Lineage = fifelse(cell_type_levels %in% c("Colonic_absorptive", "Ileal_absorptive", "Enteroendocrine", "Goblet", "Paneth", "Stem_TA_progenitor", "Tuft"), "Epithelial cells",
    fifelse(cell_type_levels %in% c("Fibroblast", "Endothelial", "Glial", "Pericyte"), "Non-immune stroma",
      fifelse(cell_type_levels %in% c("DC", "Mast", "Monocyte_Macrophage"), "Myeloid", "Lymphoid")))
)
lineage_colours <- c("Epithelial cells" = "#A93226", "Non-immune stroma" = "#D7CA2C", Myeloid = "#6B9E45", Lymphoid = "#B87CAF")

display_name <- function(x) {
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("Monocyte Macrophage", "Mono/Mac", x, fixed = TRUE)
  x <- gsub("Stem TA progenitor", "Stem/TA progenitor", x, fixed = TRUE)
  x
}

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

# Step3. 逐表型轴复刻 G/H 箱线图
summary_parts <- list()
for (axis_name in axes) {
  axis_data <- representation[Axis == axis_name]
  if (!nrow(axis_data)) next

  g_data <- axis_data[, .(AverageEnrichment = mean(Enrichment), PositiveContexts = .N), by = .(CellType, State, CellState)]
  g_order <- g_data[, .(MedianEnrichment = median(AverageEnrichment)), by = CellType][order(-MedianEnrichment, CellType)]$CellType
  g_data[, CellTypeDisplay := factor(display_name(CellType), levels = display_name(g_order))]
  g_data[, State := factor(State, levels = state_levels)]
  g_max <- max(g_data$AverageEnrichment, na.rm = TRUE)
  lineage_track <- merge(data.table(CellType = g_order), lineage_map, by = "CellType", sort = FALSE)
  lineage_track[, CellTypeDisplay := factor(display_name(CellType), levels = display_name(g_order))]
  lineage_track[, TrackY := -0.065 * g_max]

  h_data <- axis_data[, .(AverageEnrichment = mean(Enrichment), PositiveStates = .N), by = .(Context, DatasetShort, Level, CellType)]
  h_order <- h_data[, .(MedianEnrichment = median(AverageEnrichment)), by = Context][order(-MedianEnrichment, Context)]$Context
  h_data[, Context := factor(Context, levels = h_order)]
  h_data[, CellType := factor(CellType, levels = cell_type_levels)]
  h_max <- max(h_data$AverageEnrichment, na.rm = TRUE)

  p_g <- ggplot(g_data, aes(x = CellTypeDisplay, y = AverageEnrichment)) +
    geom_boxplot(width = 0.72, outlier.shape = NA, fill = "#D9D9D9", colour = "black", linewidth = 0.38) +
    geom_point(aes(colour = State), position = position_jitter(width = 0.13, height = 0, seed = 20261124), size = 1.75) +
    geom_tile(data = lineage_track, aes(x = CellTypeDisplay, y = TrackY, fill = Lineage), inherit.aes = FALSE, width = 0.98, height = 0.055 * g_max) +
    scale_colour_manual(values = state_colours, drop = FALSE, name = "Cell state") +
    scale_fill_manual(values = lineage_colours, drop = FALSE, name = "Cell lineage") +
    scale_y_continuous(limits = c(-0.11 * g_max, 1.04 * g_max), breaks = pretty(c(0, g_max), n = 5), expand = expansion(mult = c(0, 0.02))) +
    labs(tag = "G", title = "Enrichment of cell states across phenotype contexts", subtitle = "(averaged across datasets and phenotype levels)", x = NULL, y = expression(paste("Average enrichment (", -log[10], " q-value)"))) +
    theme_classic(base_size = 7.4, base_family = "sans") +
    theme(
      axis.text.x = element_text(size = 6.2, angle = 48, hjust = 1, colour = "black"), axis.text.y = element_text(size = 6.5, colour = "black"),
      axis.title.y = element_text(size = 7.3), axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
      plot.title = element_text(size = 9.2, hjust = 0.5, face = "plain"), plot.subtitle = element_text(size = 8.1, hjust = 0.5, colour = "black"),
      plot.tag = element_text(size = 16, face = "bold"), plot.tag.position = c(0.015, 0.985),
      legend.position = "right", legend.title = element_text(size = 7.2), legend.text = element_text(size = 6.6),
      legend.key.height = unit(3.4, "mm"), legend.key.width = unit(3.4, "mm"), plot.margin = margin(8, 4, 5, 8)
    ) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2.5)), fill = guide_legend(order = 2, override.aes = list(shape = 22, size = 3.8)))

  p_h <- ggplot(h_data, aes(x = Context, y = AverageEnrichment)) +
    geom_boxplot(width = 0.72, outlier.shape = NA, fill = "#D9D9D9", colour = "black", linewidth = 0.38) +
    geom_point(aes(colour = CellType), position = position_jitter(width = 0.13, height = 0, seed = 20261124), size = 1.65) +
    scale_colour_manual(values = cell_type_colours, labels = display_name(names(cell_type_colours)), drop = FALSE, name = "Cell type") +
    scale_y_continuous(limits = c(0, 1.04 * h_max), breaks = pretty(c(0, h_max), n = 5), expand = expansion(mult = c(0, 0.02))) +
    labs(tag = "H", title = "Enrichment of cell states in phenotype contexts", subtitle = "(averaged by cell type)", x = NULL, y = expression(paste("Average enrichment (", -log[10], " q-value)"))) +
    theme_classic(base_size = 7.4, base_family = "sans") +
    theme(
      axis.text.x = element_text(size = 6.2, angle = 48, hjust = 1, colour = "black"), axis.text.y = element_text(size = 6.5, colour = "black"),
      axis.title.y = element_text(size = 7.3), axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
      plot.title = element_text(size = 9.2, hjust = 0.5, face = "plain"), plot.subtitle = element_text(size = 8.1, hjust = 0.5, colour = "black"),
      plot.tag = element_text(size = 16, face = "bold"), plot.tag.position = c(0.015, 0.985),
      legend.position = "right", legend.title = element_text(size = 7.2), legend.text = element_text(size = 6.4),
      legend.key.height = unit(3.0, "mm"), legend.key.width = unit(3.0, "mm"), plot.margin = margin(8, 8, 5, 4)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = 1))

  figure_pair <- p_g + p_h + plot_layout(ncol = 2, widths = c(1.08, 1))
  figure_stem <- file.path(output_dir, paste0("Cell2021_FigS2GH_boxplots_", axis_name))
  require_patchwork_panel_alignment(
    figure_pair,
    manifest_path = paste0(figure_stem, ".alignment-layout.json"),
    report_path = paste0(figure_stem, ".alignment.json"),
    overlay_svg = paste0(figure_stem, ".alignment.svg"),
    width_in = 240 / 25.4, height_in = 125 / 25.4,
    panel_ids = c("G", "H"), row_groups = list(c("G", "H")),
    audit_script = alignment_audit_file, python = Sys.which("python"),
    tolerance_pt = 1.5, gutter_tolerance_pt = 1.5, strict = TRUE
  )
  save_figure(figure_pair, figure_stem, 240, 125)

  g_data[, Axis := axis_name]
  h_data[, Axis := axis_name]
  fwrite(g_data, file.path(source_dir, paste0("FigS2G_", axis_name, ".tsv")), sep = "\t", quote = FALSE, na = "")
  fwrite(h_data, file.path(source_dir, paste0("FigS2H_", axis_name, ".tsv")), sep = "\t", quote = FALSE, na = "")
  summary_parts[[axis_name]] <- data.table(Axis = axis_name, PositiveRows = nrow(axis_data), GStates = nrow(g_data), HCellTypeContexts = nrow(h_data))
}

# Step4. 保存完成记录
summary_table <- rbindlist(summary_parts)
fwrite(summary_table, file.path(source_dir, "FigS2GH_boxplot_summary.tsv"), sep = "\t", quote = FALSE)
completion <- list(
  status = "PASS_CELL2021_S2GH_BOXPLOTS", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  formal_axes = summary_table$Axis, min_parent_cells = 20L, positive_enrichment = "OddsRatio>1",
  enrichment_score = "min(-log10(BH_q),10)", bh_family = "DatasetScope x Axis x Level across cell states", layout = "Cell2021_FigS2G_S2H"
)
jsonlite::write_json(completion, file.path(output_dir, "Cell2021_FigS2GH_boxplots_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(summary_table)
print(completion)
cat("PASS_CELL2021_S2GH_BOXPLOTS\n")
