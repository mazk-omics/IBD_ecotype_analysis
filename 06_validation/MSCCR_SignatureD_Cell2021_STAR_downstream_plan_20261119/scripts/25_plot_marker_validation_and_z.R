#!/usr/bin/env Rscript
# 本脚本用 R 输出本轮正式 marker validation 热图，以及 recovery 和 phenotype enrichment 的 Z 值图。
# 图中保留不可评价格子；颜色变换只用于显示，原始数值完整写入 source-data 表。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读入本轮正式结果
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 25_plot_marker_validation_and_z.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1L]], mustWork = TRUE)

required_pkgs <- c("data.table", "ggplot2", "jsonlite", "svglite", "ragg")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
library(data.table)
library(ggplot2)

marker_dir <- file.path(RUN_ROOT, "output/09_marker_recovery/ranking_and_LODO")
recovery_dir <- file.path(RUN_ROOT, "output/08_scRNA_recovery_summary")
phenotype_exact_dir <- file.path(RUN_ROOT, "output/11_known_phenotype_enrichment_exact")
phenotype_patient_dir <- file.path(RUN_ROOT, "output/12_known_phenotype_enrichment_patient_aware")
output_dir <- file.path(RUN_ROOT, "output/14_validation_figures")
detailed_dir <- file.path(output_dir, "marker_LODO_detailed_by_parent")
dir.create(detailed_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  file.path(marker_dir, "marker_LODO_per_dataset.tsv"),
  file.path(marker_dir, "marker_LODO_mean_expression.tsv.gz"),
  file.path(recovery_dir, "scRNA_recovery_state_by_dataset.tsv"),
  file.path(recovery_dir, "scRNA_recovery_stouffer_meta.tsv"),
  file.path(phenotype_exact_dir, "state_known_phenotype_enrichment_exact.tsv"),
  file.path(phenotype_patient_dir, "state_known_phenotype_enrichment_patient_aware.tsv")
)
if (any(!file.exists(input_files))) stop("Missing marker, recovery or phenotype results required for plotting", call. = FALSE)

lodo <- fread(input_files[[1L]], colClasses = "character")
lodo_expression <- fread(input_files[[2L]])
recovery <- fread(input_files[[3L]])
recovery_meta <- fread(input_files[[4L]])
phenotype_exact <- fread(input_files[[5L]])
phenotype_patient <- fread(input_files[[6L]])

lodo[, Concordant := as.numeric(Concordant)]
recovery[, Z := as.numeric(Z)]
recovery[, CellState := paste(CellType, State, sep = "_")]
recovery_meta[, MetaZ := as.numeric(MetaZ)]
phenotype_exact[, Z := as.numeric(Z)]
phenotype_patient[, Z := as.numeric(Z)]
lodo_expression[, MeanUnitVarianceExpression := as.numeric(MeanUnitVarianceExpression)]
lodo_expression[, CellType := sub("_S[0-9]+$", "", MarkerState)]

state_order <- unique(recovery[order(CellType, State), CellState])
dataset_order <- c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999")
recovery[, CellState := factor(CellState, levels = rev(state_order))]
lodo[, CellState := factor(CellState, levels = rev(state_order))]
recovery[, Dataset := factor(Dataset, levels = dataset_order)]
lodo[, HeldOut := factor(HeldOut, levels = dataset_order)]

theme_validation <- theme_minimal(base_size = 6.5, base_family = "sans") +
  theme(
    panel.grid = element_blank(), axis.title = element_text(size = 6.5), axis.text = element_text(size = 5.5),
    plot.title = element_text(size = 7.5, face = "bold"), plot.subtitle = element_text(size = 6),
    legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.5),
    strip.text = element_text(size = 6, face = "bold"), strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    plot.margin = margin(4, 5, 4, 5, unit = "pt")
  )
theme_set(theme_validation)

save_plot <- function(plot_object, file_stem, width_mm, height_mm, tiff = TRUE) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(file_stem, ".svg"), width = width_in, height = height_in)
  print(plot_object)
  dev.off()
  grDevices::cairo_pdf(paste0(file_stem, ".pdf"), width = width_in, height = height_in, family = "sans")
  print(plot_object)
  dev.off()
  if (tiff) {
    ragg::agg_tiff(paste0(file_stem, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
    print(plot_object)
    dev.off()
  }
}

# Step2. 输出 marker LODO 汇总热图
fwrite(lodo, file.path(output_dir, "source_marker_LODO_concordance.tsv"), sep = "\t", quote = FALSE, na = "")
for (variant_name in c("Disease", "Region3")) {
  plot_data <- lodo[Variant == variant_name]
  plot_data[, DisplayValue := fifelse(Status == "EVALUABLE", Concordant, NA_real_)]
  marker_summary_plot <- ggplot(plot_data, aes(x = HeldOut, y = CellState, fill = DisplayValue)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    scale_fill_gradient2(low = "#B5524C", mid = "#F7F7F7", high = "#2C7FB8", midpoint = 0.5, limits = c(0, 1), breaks = c(0, 1), labels = c("Discordant", "Concordant"), na.value = "#D9D9D9", name = "LODO") +
    labs(title = paste0("Marker validation across held-out datasets (", variant_name, ")"), subtitle = "Top-10 markers selected without the held-out dataset; grey = not evaluable", x = "Held-out scRNA-seq dataset", y = "Cell state") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.border = element_rect(colour = "#555555", fill = NA, linewidth = 0.25))
  save_plot(marker_summary_plot, file.path(output_dir, paste0("marker_LODO_concordance_", variant_name)), 120, 285)
}

# Step3. 输出各父类型 marker-expression 对角热图
parent_types <- sort(unique(lodo_expression$CellType))
for (cell_type in parent_types) {
  plot_data <- lodo_expression[CellType == cell_type]
  plot_data[, MarkerState := factor(MarkerState, levels = rev(sort(unique(MarkerState))))]
  plot_data[, ExpressionState := factor(ExpressionState, levels = sort(unique(ExpressionState)))]
  marker_expression_plot <- ggplot(plot_data, aes(x = ExpressionState, y = MarkerState, fill = MeanUnitVarianceExpression)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    facet_grid(Variant ~ HeldOut) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426", midpoint = 0, na.value = "#D9D9D9", name = "Mean z-expression") +
    labs(title = paste0(cell_type, ": held-out marker specificity"), subtitle = "Rows: marker state; columns: expression state; expected validation signal lies on the diagonal", x = "Held-out expression state", y = "Trained marker state") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(2, "pt"), panel.border = element_rect(colour = "#555555", fill = NA, linewidth = 0.25))
  save_plot(marker_expression_plot, file.path(detailed_dir, paste0("marker_LODO_expression_", gsub("[^A-Za-z0-9]+", "_", cell_type))), 183, 110, tiff = FALSE)
}

# Step4. 输出 state recovery Z 热图和 Meta-Z 图
recovery[, DisplayZ := sign(Z) * log10(1 + abs(Z))]
fwrite(recovery, file.path(output_dir, "source_state_recovery_Z.tsv"), sep = "\t", quote = FALSE, na = "")
recovery_z_plot <- ggplot(recovery, aes(x = Dataset, y = CellState, fill = DisplayZ)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  geom_point(data = recovery[Significant_Zgt1.96 == TRUE], shape = 8, size = 0.8, colour = "black", stroke = 0.15) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426", midpoint = 0, na.value = "#D9D9D9", name = "signed log10(1+|Z|)") +
  labs(title = "Cell-state recovery across four scRNA-seq datasets", subtitle = "Asterisks mark the Cell 2021 criterion Z > 1.96; grey = not evaluable", x = "Validation dataset", y = "Cell state") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.border = element_rect(colour = "#555555", fill = NA, linewidth = 0.25))
save_plot(recovery_z_plot, file.path(output_dir, "state_recovery_Z_heatmap"), 125, 285)

recovery_meta[, CellState := paste(CellType, State, sep = "_")]
recovery_meta[, CellState := factor(CellState, levels = state_order)]
fwrite(recovery_meta, file.path(output_dir, "source_state_recovery_MetaZ.tsv"), sep = "\t", quote = FALSE, na = "")
meta_z_plot <- ggplot(recovery_meta, aes(x = CellState, y = MetaZ, colour = MetaSignificant_Zgt1.96)) +
  geom_hline(yintercept = 1.96, linewidth = 0.3, linetype = 2, colour = "#555555") +
  geom_point(size = 1.1, alpha = 0.9) +
  scale_colour_manual(values = c(`TRUE` = "#B40426", `FALSE` = "#777777", `NA` = "#D9D9D9"), name = "Meta Z > 1.96") +
  labs(title = "Equal-weight Stouffer recovery evidence", subtitle = "Meta Z summarizes only datasets in which each state was formally evaluable", x = "Cell state", y = "Stouffer Meta Z") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.25), legend.position = "top")
save_plot(meta_z_plot, file.path(output_dir, "state_recovery_Stouffer_MetaZ"), 183, 95)

# Step5. 输出 exact 和患者感知 phenotype Z 热图
phenotype_sources <- list(Cell2021_exact = phenotype_exact, Patient_aware = phenotype_patient)
for (method_name in names(phenotype_sources)) {
  method_data <- phenotype_sources[[method_name]]
  method_data[, CellState := paste(CellType, State, sep = "_")]
  method_data[, CellState := factor(CellState, levels = rev(state_order))]
  method_data[, DisplayZ := sign(Z) * log10(1 + abs(Z))]
  fwrite(method_data, file.path(output_dir, paste0("source_phenotype_Z_", method_name, ".tsv")), sep = "\t", quote = FALSE, na = "")

  for (axis_name in c("Disease", "Inflammation", "Region3")) {
    plot_data <- method_data[Axis == axis_name]
    plot_data[, LevelDisplay := fifelse(is.na(Level) | !nzchar(Level), "Not evaluable", Level)]
    plot_data[, DatasetLevel := paste(Dataset, LevelDisplay, sep = " | ")]
    phenotype_z_plot <- ggplot(plot_data, aes(x = DatasetLevel, y = CellState, fill = DisplayZ)) +
      geom_tile(colour = "white", linewidth = 0.1) +
      geom_point(data = plot_data[PositiveZgt1.96 == TRUE], shape = 8, size = 0.65, colour = "black", stroke = 0.12) +
      scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426", midpoint = 0, na.value = "#D9D9D9", name = "signed log10(1+|Z|)") +
      labs(title = paste0(axis_name, " enrichment Z (", gsub("_", " ", method_name), ")"), subtitle = "Asterisks mark Z > 1.96; exact and patient-aware tracks are never pooled", x = "Dataset | phenotype level", y = "Cell state") +
      theme(axis.text.x = element_text(angle = 55, hjust = 1), panel.border = element_rect(colour = "#555555", fill = NA, linewidth = 0.25))
    save_plot(phenotype_z_plot, file.path(output_dir, paste0("phenotype_Z_", method_name, "_", axis_name)), 183, 285)
  }
}

# Step6. 保存图件方法和完成记录
figure_contract <- c(
  "Core conclusion: marker specificity and recovery evidence must remain visible across all four validation datasets, including non-evaluable states.",
  "Results-level question: which Signature-D cell states retain state-specific marker expression and significant recovery in independent scRNA-seq cohorts?",
  "Figure archetype: quantitative heatmap grid.",
  "Backend: R only.",
  "Statistics: whole-dataset LODO concordance; official recovery Z; Cell 2021 exact phenotype Z; separate patient-aware IBD extension.",
  "Display transform: signed log10(1+abs(Z)) for heatmap color only; raw Z retained in source-data tables.",
  "Missingness: grey indicates not evaluable and is never converted to zero."
)
writeLines(figure_contract, file.path(output_dir, "FIGURE_CONTRACT.txt"))

completion <- list(
  status = "PASS_MARKER_AND_Z_FIGURES", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  backend = "R", formal_states = length(state_order), parent_marker_heatmaps = length(parent_types),
  recovery_datasets = uniqueN(recovery$Dataset), phenotype_tracks = names(phenotype_sources),
  display_transform = "sign(Z)*log10(1+abs(Z))"
)
jsonlite::write_json(completion, file.path(output_dir, "validation_figures_complete.json"), pretty = TRUE, auto_unbox = TRUE)
print(completion)
cat("PASS_MARKER_AND_Z_FIGURES\n")
