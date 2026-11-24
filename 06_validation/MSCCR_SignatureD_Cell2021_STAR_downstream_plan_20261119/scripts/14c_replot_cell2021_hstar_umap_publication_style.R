#!/usr/bin/env Rscript
# 本脚本仅重绘严格 Cell 2021 H* UMAP 坐标，不重新计算 UMAP 或修改 state。
# 版式继承原文 Figure 2B 的点图密度、配色和样本数标注。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定输入、输出和绘图环境
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 14c_replot_cell2021_hstar_umap_publication_style.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("data.table", "ggplot2", "patchwork", "jsonlite", "svglite", "ragg")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(ggplot2)
library(patchwork)

source(file.path(RUN_ROOT, "qa_tools/panel_alignment.R"))
source(file.path(RUN_ROOT, "qa_tools/panel_alignment_patchwork_compat.R"))
result_dir <- file.path(output_root, "10_state_annotation/hstar_umap_cell2021_exact")
coordinate_file <- file.path(result_dir, "hstar_umap_coordinates.tsv.gz")
manifest_file <- file.path(result_dir, "hstar_umap_manifest.tsv")
if (!file.exists(coordinate_file) || !file.exists(manifest_file)) stop("Strict Cell 2021 UMAP output is missing", call. = FALSE)

coordinates <- fread(coordinate_file)
manifest <- fread(manifest_file)
if (nrow(coordinates) != sum(manifest$UmapSamples) || any(manifest$UmapSamples != manifest$FormalAssigned)) stop("Strict UMAP identity failed", call. = FALSE)

cell_order <- manifest$CellType
cell_labels <- c(
  B = "B cells", CD4_T = "CD4 T cells", CD8_T = "CD8 T cells", Colonic_absorptive = "Colonic absorptive",
  DC = "Dendritic cells", Endothelial = "Endothelial cells", Enteroendocrine = "Enteroendocrine", Fibroblast = "Fibroblasts",
  Glial = "Glial cells", Goblet = "Goblet cells", ILC = "ILCs", Ileal_absorptive = "Ileal absorptive",
  Mast = "Mast cells", Monocyte_Macrophage = "Monocytes/macrophages", NK = "NK cells", Paneth = "Paneth cells",
  Pericyte = "Pericytes", Plasma = "Plasma cells", Stem_TA_progenitor = "Stem/TA progenitors", Tuft = "Tuft cells"
)
if (!setequal(cell_order, names(cell_labels))) stop("Cell type label mapping failed", call. = FALSE)

state_palette <- c(S01 = "#E7298A", S02 = "#40A4C1", S03 = "#52B447", S04 = "#7E4AAA", S05 = "#FF7F00", S06 = "#FDC000")
coordinates[, CellType := factor(CellType, levels = cell_order)]
coordinates[, DisplayState := factor(as.character(DisplayState), levels = names(state_palette))]

# Step2. 保留严格坐标并按原文风格绘制 20 个等宽面板
plots <- vector("list", length(cell_order))
for (cell_index in seq_along(cell_order)) {
  cell_type <- cell_order[[cell_index]]
  plot_data <- coordinates[CellType == cell_type]
  sample_n <- nrow(plot_data)
  if (sample_n != manifest[CellType == cell_type, UmapSamples]) stop("Panel sample count failed: ", cell_type, call. = FALSE)

  state_n <- manifest[CellType == cell_type, States]
  row_index <- ceiling(cell_index / 4)
  column_index <- ((cell_index - 1L) %% 4L) + 1L
  x_mid <- mean(range(plot_data$UMAP1))
  y_mid <- mean(range(plot_data$UMAP2))
  display_span <- max(diff(range(plot_data$UMAP1)), diff(range(plot_data$UMAP2)))
  if (!is.finite(display_span) || display_span <= 0) stop("Panel display span failed: ", cell_type, call. = FALSE)
  plot_data[, `:=`(
    UMAP1Display = (UMAP1 - x_mid) / display_span,
    UMAP2Display = (UMAP2 - y_mid) / display_span
  )]

  plots[[cell_index]] <- ggplot(plot_data, aes(UMAP1Display, UMAP2Display, colour = DisplayState)) +
    geom_point(size = 0.36, alpha = 0.90, shape = 16, stroke = 0) +
    scale_colour_manual(values = state_palette, drop = FALSE) +
    coord_equal(xlim = c(-0.55, 0.55), ylim = c(-0.55, 0.55), expand = FALSE, clip = "off") +
    labs(
      title = paste0(cell_labels[[cell_type]], "  (n = ", state_n, ")"),
      x = if (row_index == 5L) "UMAP1" else NULL,
      y = if (column_index == 1L) "UMAP2" else NULL
    ) +
    theme_classic(base_size = 6.5, base_family = "sans") +
    theme(
      plot.title = element_text(size = 6.5, face = "plain", hjust = 0.5, margin = margin(b = 1.5)),
      axis.title = element_text(size = 6.2), axis.text = element_blank(), axis.ticks = element_blank(),
      axis.line = element_line(linewidth = 0.30, colour = "black"), legend.position = "none",
      plot.margin = margin(2.5, 2.5, 2.5, 2.5)
    )
}

figure <- wrap_plots(plots, ncol = 4) +
  plot_annotation(
    title = "Landscape of cell states across the IBD discovery cohort",
    subtitle = "Assigned biopsies only; precomputed Euclidean distances on H*; R umap 0.2.0.0 default parameters",
    theme = theme(
      plot.title = element_text(size = 9, face = "bold", family = "sans", hjust = 0),
      plot.subtitle = element_text(size = 6.5, family = "sans", hjust = 0, margin = margin(b = 4)),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

# Step3. 审计面板对齐并导出可编辑矢量和高分辨率图
width_mm <- 183
height_mm <- 220
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
panel_ids <- sprintf("p%02d", seq_along(cell_order))
row_groups <- split(panel_ids, rep(seq_len(5), each = 4))
column_groups <- split(panel_ids, rep(seq_len(4), times = 5))

require_patchwork_panel_alignment(
  plot = figure,
  manifest_path = file.path(result_dir, "hstar_umap_cell2021_publication.alignment-layout.json"),
  report_path = file.path(result_dir, "hstar_umap_cell2021_publication.alignment.json"),
  overlay_svg = file.path(result_dir, "hstar_umap_cell2021_publication.alignment.svg"),
  width_in = width_in,
  height_in = height_in,
  panel_ids = panel_ids,
  row_groups = row_groups,
  column_groups = column_groups,
  audit_script = file.path(RUN_ROOT, "qa_tools/audit_panel_alignment.py"),
  python = "/usr/bin/python3",
  tolerance_pt = 1.5,
  gutter_tolerance_pt = 1.5,
  strict = TRUE
)

svglite::svglite(file.path(result_dir, "hstar_umap_cell2021_publication.svg"), width = width_in, height = height_in)
print(figure)
dev.off()
grDevices::cairo_pdf(file.path(result_dir, "hstar_umap_cell2021_publication.pdf"), width = width_in, height = height_in, family = "sans")
print(figure)
dev.off()
ragg::agg_png(file.path(result_dir, "hstar_umap_cell2021_publication.png"), width = width_in, height = height_in, units = "in", res = 600, background = "white")
print(figure)
dev.off()
ragg::agg_tiff(file.path(result_dir, "hstar_umap_cell2021_publication.tiff"), width = width_in, height = height_in, units = "in", res = 600, background = "white", compression = "lzw")
print(figure)
dev.off()

completion <- list(
  status = "PASS_CELL2021_HSTAR_UMAP_PUBLICATION_REPLOT",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  coordinate_rows = nrow(coordinates), cell_types = length(cell_order),
  coordinate_source = "strict Cell 2021 H* Euclidean-distance UMAP",
  geometry_changed = FALSE, style_only = TRUE,
  display_transform = "per-panel translation and isotropic scaling; pairwise geometry preserved up to one constant",
  width_mm = width_mm, height_mm = height_mm, dpi = 600
)
jsonlite::write_json(completion, file.path(result_dir, "hstar_umap_cell2021_publication_complete.json"), pretty = TRUE, auto_unbox = TRUE)
cat("PASS_CELL2021_HSTAR_UMAP_PUBLICATION_REPLOT\n")
