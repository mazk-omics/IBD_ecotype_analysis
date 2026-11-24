#!/usr/bin/env Rscript
# 本脚本按 Cell 2021 Figure 2D 的单面板版式，分别绘制 GSE57945 bulk 与四个 scRNA 队列的 nonIBD 富集一致性。
# 热图显示单侧 P 的 -log10 值并截断在 2；下方报告 Z > 1.96 调用的一致率和 Fisher 精确检验。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结的 concordance source data
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 46_plot_Cell2021_Figure2D_nonIBD_concordance.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(grid)

input_dir <- file.path(RUN_ROOT, "output/33_Cell2021_bulk_scRNA_nonIBD_concordance")
state_file <- file.path(input_dir, "Cell2021_bulk_scRNA_nonIBD_state_concordance.tsv")
summary_file <- file.path(input_dir, "Cell2021_bulk_scRNA_nonIBD_concordance_summary.tsv")
pass_file <- file.path(input_dir, "Cell2021_bulk_scRNA_nonIBD_concordance_complete.json")
output_dir <- file.path(RUN_ROOT, "output/34_Cell2021_Figure2D_nonIBD_concordance")
source_dir <- file.path(output_dir, "source_data")

if (!file.exists(state_file) || !file.exists(summary_file) || !file.exists(pass_file)) stop("Concordance outputs are incomplete", call. = FALSE)
pass <- jsonlite::read_json(pass_file, simplifyVector = TRUE)
if (pass$status != "PASS_CELL2021_BULK_SCRNA_NONIBD_CONCORDANCE") stop("Concordance audit did not pass", call. = FALSE)
states <- fread(state_file)
summary <- fread(summary_file)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(states, file.path(source_dir, "Figure2D_nonIBD_state_source_data.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(summary, file.path(source_dir, "Figure2D_nonIBD_concordance_statistics.tsv"), sep = "\t", quote = FALSE, na = "")

# Step2. 固定原文式画布和配色
datasets <- c("GSE214695", "IBDverse_E_MTAB_16999", "SCP1884", "SCP259")
figure_width_in <- 2.60
figure_height_in <- 4.15
heatmap_left <- 0.32
heatmap_right <- 0.52
heatmap_bottom <- 0.26
heatmap_top <- 0.59
heatmap_colours <- colorRampPalette(c("#FFFFFF", "#D8DFF0", "#8FA4D0", "#3F57A6"))(201)

format_p_expression <- function(p_value) {
  if (!is.finite(p_value)) return(expression(italic(P) == NA))
  if (p_value < 0.001) {
    exponent <- floor(log10(p_value))
    mantissa <- signif(p_value / 10^exponent, 2)
    return(bquote(italic(P) == .(mantissa) %*% 10^.(exponent)))
  }
  bquote(italic(P) == .(formatC(p_value, format = "f", digits = 3)))
}

draw_figure <- function(plot_data, plot_summary) {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text("D", x = 0.055, y = 0.965, just = c("left", "top"), gp = gpar(fontfamily = "sans", fontsize = 18, fontface = "bold"))
  grid.text("Validation of cell states", x = 0.57, y = 0.955, just = c("center", "top"), gp = gpar(fontfamily = "sans", fontsize = 10.5))
  grid.text("enriched in non-IBD", x = 0.57, y = 0.910, just = c("center", "top"), gp = gpar(fontfamily = "sans", fontsize = 10.5))
  grid.text("tissues", x = 0.57, y = 0.865, just = c("center", "top"), gp = gpar(fontfamily = "sans", fontsize = 10.5))

  column_centres <- seq(heatmap_left, heatmap_right, length.out = 3L)[1:2] + (heatmap_right - heatmap_left) / 4
  grid.text(c("Bulk RNA-seq", "scRNA-seq"), x = column_centres, y = 0.625, rot = 45, just = c("left", "center"), gp = gpar(fontfamily = "sans", fontsize = 7.8))
  grid.text(paste0(nrow(plot_data), " cell states"), x = 0.245, y = (heatmap_bottom + heatmap_top) / 2, rot = 90, gp = gpar(fontfamily = "sans", fontsize = 9.5))

  row_height <- (heatmap_top - heatmap_bottom) / nrow(plot_data)
  col_width <- (heatmap_right - heatmap_left) / 2
  for (row_index in seq_len(nrow(plot_data))) {
    y_center <- heatmap_top - (row_index - 0.5) * row_height
    values <- c(plot_data$BulkScore[[row_index]], plot_data$ScRNAScore[[row_index]])
    colour_index <- pmax(1L, pmin(201L, as.integer(round(values / 2 * 200)) + 1L))
    grid.rect(x = column_centres, y = y_center, width = col_width, height = row_height, gp = gpar(fill = heatmap_colours[colour_index], col = NA))
  }
  grid.rect(x = (heatmap_left + heatmap_right) / 2, y = (heatmap_bottom + heatmap_top) / 2,
            width = heatmap_right - heatmap_left, height = heatmap_top - heatmap_bottom, gp = gpar(fill = NA, col = "black", lwd = 0.8))
  grid.lines(x = unit(rep((heatmap_left + heatmap_right) / 2, 2), "npc"), y = unit(c(heatmap_bottom, heatmap_top), "npc"), gp = gpar(col = "#B5B5B5", lwd = 0.25))

  grid.text("Non-IBD", x = 0.73, y = 0.425, gp = gpar(fontfamily = "sans", fontsize = 8.5))
  grid.text("enrichment", x = 0.73, y = 0.395, gp = gpar(fontfamily = "sans", fontsize = 8.5))
  grid.text("(-log10 p-value)", x = 0.73, y = 0.365, gp = gpar(fontfamily = "sans", fontsize = 8.5))
  legend_x <- seq(0.63, 0.86, length.out = 201L)
  for (colour_index in seq_along(heatmap_colours)) {
    grid.rect(x = legend_x[[colour_index]], y = 0.315, width = 0.23 / 201, height = 0.024, gp = gpar(fill = heatmap_colours[[colour_index]], col = NA))
  }
  grid.rect(x = 0.745, y = 0.315, width = 0.23, height = 0.024, gp = gpar(fill = NA, col = "black", lwd = 0.7))
  grid.text("0", x = 0.61, y = 0.315, gp = gpar(fontfamily = "sans", fontsize = 8.5))
  grid.text("2", x = 0.88, y = 0.315, gp = gpar(fontfamily = "sans", fontsize = 8.5))

  grid.text(paste0("Concordance: ", round(plot_summary$ConcordancePercent), "%"), x = 0.42, y = 0.165, gp = gpar(fontfamily = "sans", fontsize = 10))
  grid.text(format_p_expression(plot_summary$FisherTwoSidedP), x = 0.42, y = 0.120, gp = gpar(fontfamily = "sans", fontsize = 10))
}

# Step3. 为四个 scRNA validation 队列分别导出同版式图
for (dataset in datasets) {
  plot_data <- states[ScRNADataset == dataset][order(DisplayOrder)]
  plot_summary <- summary[ScRNADataset == dataset]
  if (!nrow(plot_data) || nrow(plot_summary) != 1L) stop("Missing plot data for ", dataset, call. = FALSE)
  file_stem <- file.path(output_dir, paste0("Figure2D_nonIBD_GSE57945_vs_", dataset))

  grDevices::cairo_pdf(paste0(file_stem, ".pdf"), width = figure_width_in, height = figure_height_in, family = "sans")
  draw_figure(plot_data, plot_summary)
  dev.off()

  grDevices::svg(paste0(file_stem, ".svg"), width = figure_width_in, height = figure_height_in, family = "sans", onefile = TRUE)
  draw_figure(plot_data, plot_summary)
  dev.off()

  grDevices::tiff(paste0(file_stem, ".tiff"), width = figure_width_in, height = figure_height_in, units = "in", res = 600, compression = "lzw", type = "cairo")
  draw_figure(plot_data, plot_summary)
  dev.off()

  grDevices::png(paste0(file_stem, ".png"), width = figure_width_in, height = figure_height_in, units = "in", res = 300, type = "cairo-png")
  draw_figure(plot_data, plot_summary)
  dev.off()
  cat("FIGURE2D_EXPORTED", dataset, "\n")
}

# Step4. 写出图件完成记录
completion <- list(
  status = "PASS_CELL2021_FIGURE2D_NONIBD_CONCORDANCE", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  backend = "R_GRID", bulk_dataset = "GSE57945", scrna_datasets = datasets, figures = length(datasets),
  score = "one-sided -log10(P) capped at 2", call_threshold = "Z > 1.96", concordance = "overall binary agreement",
  fisher_test = "two-sided Fisher exact test", exports = c("PDF", "SVG", "TIFF_600dpi", "PNG_300dpi")
)
jsonlite::write_json(completion, file.path(output_dir, "Cell2021_Figure2D_nonIBD_concordance_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(summary)
print(completion)
cat("PASS_CELL2021_FIGURE2D_NONIBD_CONCORDANCE\n")
