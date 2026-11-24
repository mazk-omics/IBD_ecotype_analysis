#!/usr/bin/env Rscript
# 本脚本将 G7A dominant-state Fisher 结果画成 Cell 2021 Figure S2G-S2H 对应的 IBD representation 热图。
# pooled 与 within-dataset 结果分开，灰色保留不可评价格子，IBD/nonIBD 不在无健康对照时作图。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读入 representation 结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 26_plot_cell_state_representation.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "ggplot2", "jsonlite", "svglite", "ragg")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(ggplot2)

input_dir <- file.path(RUN_ROOT, "output/13_cell_state_representation")
output_dir <- file.path(RUN_ROOT, "output/15_cell_state_representation_figures")
input_file <- file.path(input_dir, "sample_dominant_state_representation_fisher.tsv")
gate_file <- file.path(input_dir, "IBD_nonIBD_representation_gate.tsv")
if (!file.exists(input_file) || !file.exists(gate_file)) stop("Run script 24 before representation plotting", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

representation <- fread(input_file)
nonibd_gate <- fread(gate_file, colClasses = "character")
required_columns <- c("Scope", "DatasetScope", "Axis", "Level", "CellType", "State", "CellState", "OddsRatio", "FDR_BH_across_cell_types", "SignificantEnrichment", "Status")
if (!all(required_columns %in% names(representation))) stop("Representation columns changed", call. = FALSE)

representation[, OddsRatio := as.numeric(OddsRatio)]
representation[, FDR_BH_across_cell_types := as.numeric(FDR_BH_across_cell_types)]
representation[, EnrichmentScore := fifelse(Status == "EVALUABLE" & is.finite(FDR_BH_across_cell_types), pmin(-log10(pmax(FDR_BH_across_cell_types, 1e-300)), 10), NA_real_)]
representation[, SignificantEnrichment := as.logical(SignificantEnrichment)]
state_order <- unique(representation[order(CellType, State), CellState])
representation[, CellState := factor(CellState, levels = rev(state_order))]

theme_representation <- theme_minimal(base_size = 6.5, base_family = "sans") +
  theme(
    panel.grid = element_blank(), axis.title = element_text(size = 6.5), axis.text = element_text(size = 5.5),
    plot.title = element_text(size = 7.5, face = "bold"), plot.subtitle = element_text(size = 6),
    legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.5),
    panel.border = element_rect(colour = "#555555", fill = NA, linewidth = 0.25),
    plot.margin = margin(4, 5, 4, 5, unit = "pt")
  )
theme_set(theme_representation)

save_plot <- function(plot_object, file_stem, width_mm, height_mm) {
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
}

# Step2. 输出 pooled four-dataset representation
pooled <- representation[Scope == "POOLED_FOUR_DATASETS"]
fwrite(pooled, file.path(output_dir, "source_representation_pooled.tsv"), sep = "\t", quote = FALSE, na = "")
for (axis_name in c("Disease", "Inflammation", "Region3", "Dataset")) {
  plot_data <- pooled[Axis == axis_name]
  pooled_plot <- ggplot(plot_data, aes(x = Level, y = CellState, fill = EnrichmentScore)) +
    geom_tile(colour = "white", linewidth = 0.1) +
    geom_point(data = plot_data[SignificantEnrichment == TRUE], shape = 8, size = 0.7, colour = "black", stroke = 0.12) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", na.value = "#D9D9D9", limits = c(0, 10), name = "-log10(BH FDR)\n(capped at 10)") +
    labs(title = paste0("Cell-state representation by ", axis_name), subtitle = "Dominant state per biological sample; asterisks mark one-sided Fisher BH FDR < 0.05", x = axis_name, y = "Cell state") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(pooled_plot, file.path(output_dir, paste0("representation_pooled_", axis_name)), 150, 285)
}

# Step3. 输出逐数据集 phenotype representation
within_dataset <- representation[Scope == "WITHIN_DATASET"]
within_dataset[, DatasetLevel := paste(DatasetScope, Level, sep = " | ")]
fwrite(within_dataset, file.path(output_dir, "source_representation_within_dataset.tsv"), sep = "\t", quote = FALSE, na = "")
for (axis_name in c("Disease", "Inflammation", "Region3")) {
  plot_data <- within_dataset[Axis == axis_name]
  within_plot <- ggplot(plot_data, aes(x = DatasetLevel, y = CellState, fill = EnrichmentScore)) +
    geom_tile(colour = "white", linewidth = 0.1) +
    geom_point(data = plot_data[SignificantEnrichment == TRUE], shape = 8, size = 0.65, colour = "black", stroke = 0.12) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", na.value = "#D9D9D9", limits = c(0, 10), name = "-log10(BH FDR)\n(capped at 10)") +
    labs(title = paste0(axis_name, " representation within each dataset"), subtitle = "Dominant state per biological sample; grey = Fisher margin not evaluable", x = "Dataset | phenotype level", y = "Cell state") +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
  save_plot(within_plot, file.path(output_dir, paste0("representation_within_dataset_", axis_name)), 183, 285)
}

# Step4. 保存图件边界和完成记录
fwrite(nonibd_gate, file.path(output_dir, "IBD_nonIBD_representation_gate.tsv"), sep = "\t", quote = FALSE)
figure_contract <- c(
  "Core conclusion: state representation varies across IBD phenotypes and validation datasets, while unassigned and non-evaluable strata remain explicit.",
  "Results-level question: which frozen Signature-D states are over-represented in each phenotype or dataset at the biological-sample level?",
  "Figure archetype: quantitative heatmap grid.",
  "Backend: R only.",
  "Statistics: dominant state per biological sample; one-sided Fisher exact test; BH correction across cell types within each scope and axis.",
  "IBD/nonIBD: not plotted when the current four cohorts contain IBD cells only; deferred to G11B.",
  "Missingness: grey indicates a non-evaluable Fisher margin and is never converted to zero."
)
writeLines(figure_contract, file.path(output_dir, "FIGURE_CONTRACT.txt"))

completion <- list(
  status = "PASS_CELL_STATE_REPRESENTATION_FIGURES", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  backend = "R", formal_states = length(state_order), pooled_axes = 4L, within_dataset_axes = 3L,
  ibd_nonibd_status = nonibd_gate$Status[[1L]]
)
jsonlite::write_json(completion, file.path(output_dir, "cell_state_representation_figures_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_CELL_STATE_REPRESENTATION_FIGURES\n")
