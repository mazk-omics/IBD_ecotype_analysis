#!/usr/bin/env Rscript
# 本脚本用 discovery bulk 的 Fisher-BH 结果复刻 Cell 2021 Fig. S2G/S2H 箱线图。
# G 对每个 state 跨表型水平取平均；H 对每个表型水平在父类型内跨 states 取平均。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取 discovery representation 结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 53_plot_discovery_Cell2021_S2GH_boxplots.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "ggplot2", "patchwork", "jsonlite", "svglite", "ragg")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(ggplot2)
library(patchwork)

input_dir <- file.path(RUN_ROOT, "output/39_discovery_state_representation_Cell2021")
input_file <- file.path(input_dir, "discovery_state_representation_fisher.tsv")
pass_file <- file.path(input_dir, "discovery_state_representation_complete.json")
output_dir <- file.path(RUN_ROOT, "output/40_discovery_Cell2021_S2GH_boxplots")
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

stale_files <- c(
  Sys.glob(file.path(output_dir, "Cell2021_FigS2GH_discovery_DiseaseRegion3*")),
  Sys.glob(file.path(source_dir, "*DiseaseRegion3*"))
)
if (length(stale_files)) file.remove(stale_files)

pass <- jsonlite::read_json(pass_file, simplifyVector = TRUE)
if (pass$status != "PASS_DISCOVERY_STATE_REPRESENTATION_CELL2021") stop("Discovery representation audit did not pass", call. = FALSE)
representation <- fread(input_file)
representation <- representation[Status == "EVALUABLE"]
if (!nrow(representation) || anyNA(representation$EnrichmentScore)) stop("Discovery representation contains missing evaluable scores", call. = FALSE)

panel_alignment_file <- file.path(RUN_ROOT, "scripts/figure_qa/panel_alignment.R")
panel_alignment_compat_file <- file.path(RUN_ROOT, "qa_tools/panel_alignment_patchwork_compat.R")
alignment_audit_file <- file.path(RUN_ROOT, "scripts/figure_qa/audit_panel_alignment.py")
source(panel_alignment_file)
source(panel_alignment_compat_file)
python_bin <- Sys.which("python")
if (!nzchar(python_bin)) python_bin <- Sys.which("python3")
if (!nzchar(python_bin)) stop("Python is required for figure alignment QA", call. = FALSE)

# Step2. 固定原文颜色、谱系和显示顺序
axes <- c("Disease", "Region3", "Inflammation")
axis_titles <- c(
  Disease = "disease groups",
  Region3 = "anatomical regions",
  Inflammation = "inflammation groups"
)
state_colours <- c(
  S01 = "#C2185B", S02 = "#4F86A8", S03 = "#5A963E", S04 = "#563B91",
  S05 = "#D4631B", S06 = "#E6B51C", S07 = "#7A3E1C", S08 = "#CB5BDC",
  S09 = "#858585", S10 = "#2A9D8F", S11 = "#264653"
)
cell_type_colours <- c(
  B = "#C8961A", Plasma = "#C02A74", CD8_T = "#5B5AA6", CD4_T = "#A64B16", NK = "#3E855C",
  ILC = "#7867A8", Monocyte_Macrophage = "#B39ACD", DC = "#E5A54C", Mast = "#D57D78", Fibroblast = "#94B6D6",
  Endothelial = "#777777", Glial = "#7AB7A7", Pericyte = "#77A84D", Colonic_absorptive = "#C98EBF",
  Ileal_absorptive = "#D85E9F", Enteroendocrine = "#A46AA6", Goblet = "#5B9DC0", Paneth = "#E8B65C",
  Stem_TA_progenitor = "#8DB85C", Tuft = "#6E87B9"
)
lineage_colours <- c("Epithelial cells" = "#A93226", "Non-immune stroma" = "#D7CA2C", Myeloid = "#6B9E45", Lymphoid = "#B87CAF")

cell_type_levels <- sort(unique(representation$CellType))
lineage_map <- data.table(
  CellType = cell_type_levels,
  Lineage = fifelse(
    cell_type_levels %in% c("Colonic_absorptive", "Ileal_absorptive", "Enteroendocrine", "Goblet", "Paneth", "Stem_TA_progenitor", "Tuft"),
    "Epithelial cells",
    fifelse(
      cell_type_levels %in% c("Fibroblast", "Endothelial", "Glial", "Pericyte"),
      "Non-immune stroma",
      fifelse(cell_type_levels %in% c("DC", "Mast", "Monocyte_Macrophage"), "Myeloid", "Lymphoid")
    )
  )
)

display_names <- setNames(gsub("_", " ", cell_type_levels, fixed = TRUE), cell_type_levels)
display_names["Monocyte_Macrophage"] <- "Mono/Mac"
display_names["Stem_TA_progenitor"] <- "Stem/TA progenitor"

# Step3. 逐表型轴计算原文 G/H 汇总量
summary_parts <- list()

for (axis_name in axes) {
  axis_data <- representation[Axis == axis_name]
  if (!nrow(axis_data)) next

  expected_levels <- sort(unique(axis_data$Level))
  expected_states <- uniqueN(axis_data$CellState)
  g_data <- axis_data[, .(
    AverageEnrichment = mean(EnrichmentScore),
    PhenotypeLevels = uniqueN(Level),
    MinimumEnrichment = min(EnrichmentScore),
    MaximumEnrichment = max(EnrichmentScore)
  ), by = .(CellType, State, CellState)]
  if (nrow(g_data) != expected_states || any(g_data$PhenotypeLevels != length(expected_levels))) stop("Fig. S2G aggregation lost discovery contexts: ", axis_name, call. = FALSE)

  h_data <- axis_data[, .(
    AverageEnrichment = mean(EnrichmentScore),
    StatesAveraged = .N
  ), by = .(Level, CellType)]
  if (uniqueN(h_data$Level) != length(expected_levels) || uniqueN(h_data$CellType) != 20L) stop("Fig. S2H aggregation lost discovery groups: ", axis_name, call. = FALSE)

  g_order_table <- g_data[, .(MedianEnrichment = median(AverageEnrichment)), by = CellType]
  g_order_table <- merge(g_order_table, lineage_map, by = "CellType", sort = FALSE)
  g_order_table[, Lineage := factor(Lineage, levels = names(lineage_colours))]
  setorder(g_order_table, Lineage, -MedianEnrichment, CellType)
  g_order <- g_order_table$CellType
  h_order <- h_data[, .(MedianEnrichment = median(AverageEnrichment)), by = Level][order(-MedianEnrichment, Level)]$Level
  g_data[, CellTypeDisplay := factor(display_names[CellType], levels = display_names[g_order])]
  observed_states <- names(state_colours)[names(state_colours) %in% sort(unique(g_data$State))]
  g_data[, State := factor(State, levels = observed_states)]
  h_data[, LevelDisplay := factor(gsub("[|]", "·", Level), levels = gsub("[|]", "·", h_order))]
  h_data[, CellType := factor(CellType, levels = cell_type_levels)]

  g_max <- max(g_data$AverageEnrichment)
  h_max <- max(h_data$AverageEnrichment)
  g_plot_max <- max(g_max, 0.1)
  h_plot_max <- max(h_max, 0.1)
  lineage_track <- merge(data.table(CellType = g_order), lineage_map, by = "CellType", sort = FALSE)
  lineage_track[, CellTypeDisplay := factor(display_names[CellType], levels = display_names[g_order])]
  lineage_track[, TrackY := -0.065 * g_plot_max]

  p_g <- ggplot(g_data, aes(x = CellTypeDisplay, y = AverageEnrichment)) +
    geom_boxplot(width = 0.72, outlier.shape = NA, fill = "#D9D9D9", colour = "black", linewidth = 0.38) +
    geom_point(aes(colour = State), position = position_jitter(width = 0.13, height = 0, seed = 20261124), size = 1.65) +
    geom_tile(data = lineage_track, aes(x = CellTypeDisplay, y = TrackY, fill = Lineage), inherit.aes = FALSE, width = 0.98, height = 0.055 * g_plot_max) +
    scale_colour_manual(values = state_colours[observed_states], drop = FALSE, name = "Cell state") +
    scale_fill_manual(values = lineage_colours, drop = FALSE, name = "Cell lineage") +
    scale_y_continuous(limits = c(-0.11 * g_plot_max, 1.04 * g_plot_max), breaks = pretty(c(0, g_plot_max), n = 5), expand = expansion(mult = c(0, 0.02))) +
    labs(
      tag = "G",
      title = "Enrichment of cell states in discovery phenotypes",
      subtitle = paste0("(averaged across ", length(expected_levels), " ", axis_titles[[axis_name]], ")"),
      x = NULL,
      y = expression(paste("Average enrichment (", -log[10], " q-value)"))
    ) +
    theme_classic(base_size = 7.4, base_family = "sans") +
    theme(
      axis.text.x = element_text(size = 6.1, angle = 48, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 6.5, colour = "black"),
      axis.title.y = element_text(size = 7.3),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      plot.title = element_text(size = 9.2, hjust = 0.5),
      plot.subtitle = element_text(size = 8.0, hjust = 0.5, colour = "black"),
      plot.tag = element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.015, 0.985),
      legend.position = "right",
      legend.title = element_text(size = 7.0),
      legend.text = element_text(size = 6.3),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(3.2, "mm"),
      plot.margin = margin(8, 4, 5, 8)
    ) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2.4)), fill = guide_legend(order = 2, override.aes = list(shape = 22, size = 3.6)))

  p_h <- ggplot(h_data, aes(x = LevelDisplay, y = AverageEnrichment)) +
    geom_boxplot(width = 0.72, outlier.shape = NA, fill = "#D9D9D9", colour = "black", linewidth = 0.38) +
    geom_point(aes(colour = CellType), position = position_jitter(width = 0.13, height = 0, seed = 20261124), size = 1.55) +
    scale_colour_manual(values = cell_type_colours, labels = display_names[names(cell_type_colours)], drop = FALSE, name = "Cell type") +
    scale_y_continuous(limits = c(0, 1.04 * h_plot_max), breaks = pretty(c(0, h_plot_max), n = 5), expand = expansion(mult = c(0, 0.02))) +
    labs(
      tag = "H",
      title = "Enrichment of cell states in discovery phenotypes",
      subtitle = "(averaged by cell type)",
      x = NULL,
      y = expression(paste("Average enrichment (", -log[10], " q-value)"))
    ) +
    theme_classic(base_size = 7.4, base_family = "sans") +
    theme(
      axis.text.x = element_text(size = 6.1, angle = 48, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 6.5, colour = "black"),
      axis.title.y = element_text(size = 7.3),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      plot.title = element_text(size = 9.2, hjust = 0.5),
      plot.subtitle = element_text(size = 8.0, hjust = 0.5, colour = "black"),
      plot.tag = element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.015, 0.985),
      legend.position = "right",
      legend.title = element_text(size = 7.0),
      legend.text = element_text(size = 6.1),
      legend.key.height = unit(2.9, "mm"),
      legend.key.width = unit(2.9, "mm"),
      plot.margin = margin(8, 8, 5, 4)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.4), ncol = 1))

  figure_pair <- p_g + p_h + plot_layout(ncol = 2, widths = c(1.12, 1))
  figure_stem <- file.path(output_dir, paste0("Cell2021_FigS2GH_discovery_", axis_name))
  require_patchwork_panel_alignment(
    figure_pair,
    manifest_path = paste0(figure_stem, ".alignment-layout.json"),
    report_path = paste0(figure_stem, ".alignment.json"),
    overlay_svg = paste0(figure_stem, ".alignment.svg"),
    width_in = 260 / 25.4,
    height_in = 125 / 25.4,
    panel_ids = c("G", "H"),
    row_groups = list(c("G", "H")),
    audit_script = alignment_audit_file,
    python = python_bin,
    tolerance_pt = 1.5,
    gutter_tolerance_pt = 1.5,
    strict = TRUE
  )

  svglite::svglite(paste0(figure_stem, ".svg"), width = 260 / 25.4, height = 125 / 25.4)
  print(figure_pair)
  dev.off()
  grDevices::cairo_pdf(paste0(figure_stem, ".pdf"), width = 260 / 25.4, height = 125 / 25.4, family = "sans")
  print(figure_pair)
  dev.off()
  ragg::agg_tiff(paste0(figure_stem, ".tiff"), width = 260 / 25.4, height = 125 / 25.4, units = "in", res = 600, compression = "lzw")
  print(figure_pair)
  dev.off()
  ragg::agg_png(paste0(figure_stem, ".png"), width = 260 / 25.4, height = 125 / 25.4, units = "in", res = 300)
  print(figure_pair)
  dev.off()

  g_data[, Axis := axis_name]
  h_data[, Axis := axis_name]
  fwrite(g_data, file.path(source_dir, paste0("FigS2G_discovery_", axis_name, ".tsv")), sep = "\t", quote = FALSE)
  fwrite(h_data, file.path(source_dir, paste0("FigS2H_discovery_", axis_name, ".tsv")), sep = "\t", quote = FALSE)
  summary_parts[[axis_name]] <- data.table(
    Axis = axis_name,
    PhenotypeLevels = length(expected_levels),
    CellStates = nrow(g_data),
    CellTypeContexts = nrow(h_data),
    GMaximum = g_max,
    HMaximum = h_max
  )
}

# Step4. 保存完成记录
summary_table <- rbindlist(summary_parts)
fwrite(summary_table, file.path(source_dir, "discovery_FigS2GH_boxplot_summary.tsv"), sep = "\t", quote = FALSE)
completion <- list(
  status = "PASS_DISCOVERY_CELL2021_S2GH_BOXPLOTS",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  primary_axes = c("Disease", "Region3"),
  secondary_axes = "Inflammation",
  source = "discovery bulk state assignments",
  non_positive_rows_retained = TRUE,
  score = "-log10(BH q-value), capped at 10"
)
jsonlite::write_json(completion, file.path(output_dir, "discovery_Cell2021_S2GH_boxplots_complete.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
print(summary_table)
cat("PASS_DISCOVERY_CELL2021_S2GH_BOXPLOTS\n")
