#!/usr/bin/env Rscript
# 本脚本把四个独立 scRNA-seq 队列的 marker 表达整理成 Cell 2021 Figure 2C 式热图。
# 读取既有 top-marker 排名和原始 counts 检查点；不重跑 marker recovery，也不改正式 state。
options(stringsAsFactors=FALSE, warn=1)

library(data.table)
library(Matrix)
library(ggplot2)
library(patchwork)
library(ragg)
library(svglite)
library(jsonlite)

# Step1. 定位输入并固定展示顺序
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: 59_plot_marker_validation_Cell2021_style.R RUN_ROOT", call.=FALSE)
run_root <- normalizePath(args[[1L]], winslash="/", mustWork=TRUE)
marker_root <- file.path(run_root, "output/09_marker_recovery")
full_count_root <- file.path(marker_root, "full_counts")
ranking_file <- file.path(marker_root, "ranking_and_LODO/state_top_markers.tsv")
recovery_file <- file.path(run_root, "output/08_scRNA_recovery_summary/scRNA_recovery_state_by_dataset.tsv")
out <- file.path(run_root, "output/46_marker_validation_Cell2021_figure")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
source(file.path(run_root, "scripts/figure_qa/panel_alignment.R"))

dataset_order <- c("GSE214695", "SCP1884", "SCP259", "IBDverse_E_MTAB_16999")
cell_type_order <- c("B", "Plasma", "CD8_T", "CD4_T", "NK", "Monocyte_Macrophage", "DC", "Mast", "ILC", "Fibroblast", "Endothelial", "Pericyte", "Glial", "Colonic_absorptive", "Ileal_absorptive", "Goblet", "Stem_TA_progenitor", "Paneth", "Enteroendocrine", "Tuft")
cell_type_title <- c(B="B cells", Plasma="Plasma cells", CD8_T="CD8 T cells", CD4_T="CD4 T cells", NK="NK cells", Monocyte_Macrophage="Mono/Mac", DC="Dendritic cells", Mast="Mast cells", ILC="ILC", Fibroblast="Fibroblasts", Endothelial="Endothelial cells", Pericyte="Pericytes", Glial="Glial cells", Colonic_absorptive="Colonic absorptive", Ileal_absorptive="Ileal absorptive", Goblet="Goblet cells", Stem_TA_progenitor="Stem/TA", Paneth="Paneth cells", Enteroendocrine="Enteroendocrine", Tuft="Tuft cells")
state_colours <- c(S01="#C53E73", S02="#5796B7", S03="#6FA64A", S04="#6C50A3", S05="#D46A1D", S06="#E3AE1A")
dataset_colours <- c(GSE214695="#A65E2E", SCP1884="#7FBF6B", SCP259="#9ECAE1", IBDverse_E_MTAB_16999="#B21872")
font_family <- "Helvetica"

stopifnot(file.exists(ranking_file), file.exists(recovery_file))
top_markers <- fread(ranking_file)
top_markers <- top_markers[HeldOut=="All" & Variant=="Disease" & WithinStateOrder<=2L]
top_markers[, StateNumber := as.integer(sub("^S0?", "", State))]
top_markers[, CellTypeOrder := match(CellType, cell_type_order)]
setorder(top_markers, CellTypeOrder, StateNumber, WithinStateOrder, Gene)
stopifnot(uniqueN(top_markers$CellState)==85L, setequal(unique(top_markers$CellType), cell_type_order))
print(top_markers[, .(States=uniqueN(State), Markers=.N), by=CellType][order(factor(CellType, levels=cell_type_order))])

recovery <- fread(recovery_file)
stopifnot(all(c("Dataset", "CellType", "State", "Significant_Zgt1.96") %in% names(recovery)))
recovery_pass <- unique(recovery[Significant_Zgt1.96==TRUE, .(Dataset, CellType, State)])
recovery_pass[, `:=`(DatasetOrder=match(Dataset, dataset_order), StateNumber=as.integer(sub("^S0?", "", State)))]
setorder(recovery_pass, CellType, StateNumber, DatasetOrder)
recovery_state_audit <- recovery[, .(TestedDatasets=.N, RecoveredDatasets=sum(Significant_Zgt1.96==TRUE), RecoveredDatasetNames=paste(Dataset[Significant_Zgt1.96==TRUE], collapse=",")), by=.(CellType, State)]
fwrite(recovery_state_audit, file.path(out, "Figure2C_recovery_state_display_audit.tsv"), sep="\t", quote=FALSE, na="")
print(recovery_state_audit[, .N, by=RecoveredDatasets][order(RecoveredDatasets)])

# Step2. 在每个数据集内部重算 marker 的 state 平均标准化表达
expression_parts <- list()
coverage_parts <- list()

for (dataset in dataset_order) {
    for (cell_type in cell_type_order) {
        count_file <- file.path(full_count_root, dataset, cell_type, "counts.rds")
        cell_file <- file.path(full_count_root, dataset, cell_type, "cells.tsv.gz")
        marker_one <- top_markers[CellType==cell_type]
        formal_states <- sort(unique(marker_one$State))
        if (!file.exists(count_file) || !file.exists(cell_file)) {
            coverage_parts[[length(coverage_parts)+1L]] <- data.table(Dataset=dataset, CellType=cell_type, Status="NO_PARENT_INPUT", Cells=0L, AssignedCells=0L, PresentMarkers=0L)
            next
        }

        counts <- readRDS(count_file)
        cells <- fread(cell_file, colClasses="character")
        cells[, Assigned := Assigned=="TRUE"]
        if (!methods::is(counts, "dgCMatrix") || !identical(colnames(counts), cells$CellID)) stop("Count/cell identity mismatch: ", dataset, " ", cell_type, call.=FALSE)
        totals <- Matrix::colSums(counts)
        if (any(!is.finite(totals) | totals<=0)) stop("Invalid library size: ", dataset, " ", cell_type, call.=FALSE)

        present_genes <- marker_one$Gene[marker_one$Gene %in% rownames(counts)]
        coverage_parts[[length(coverage_parts)+1L]] <- data.table(Dataset=dataset, CellType=cell_type, Status="INPUT_FOUND", Cells=ncol(counts), AssignedCells=sum(cells$Assigned), PresentMarkers=length(present_genes))
        if (!length(present_genes)) {
            rm(counts)
            invisible(gc(FALSE))
            next
        }

        values <- as.matrix(counts[present_genes, , drop=FALSE])
        values <- log2(sweep(values, 2, totals, "/")*1e6+1)
        gene_sd <- apply(values, 1, sd)
        valid_genes <- names(gene_sd)[is.finite(gene_sd) & gene_sd>0]
        if (!length(valid_genes)) {
            rm(counts, values)
            invisible(gc(FALSE))
            next
        }
        values <- values[valid_genes, , drop=FALSE]
        values <- sweep(values, 1, rowMeans(values), "-")
        values <- sweep(values, 1, gene_sd[valid_genes], "/")

        for (state in formal_states) {
            state_cells <- cells$Assigned & cells$State==state
            if (sum(state_cells)<5L) next
            state_mean <- rowMeans(values[, state_cells, drop=FALSE])
            expression_parts[[length(expression_parts)+1L]] <- data.table(Dataset=dataset, CellType=cell_type, Gene=names(state_mean), ExpressionState=state, MeanUnitVarianceExpression=as.numeric(state_mean), AssignedStateCells=sum(state_cells))
        }
        rm(counts, values)
        invisible(gc(FALSE))
    }
    cat("MARKER_FIGURE_DATASET", dataset, "complete\n")
}

expression_observed <- rbindlist(expression_parts, fill=TRUE)
coverage <- rbindlist(coverage_parts, fill=TRUE)
if (!nrow(expression_observed)) stop("Marker-expression aggregation is empty", call.=FALSE)

full_parts <- list()
for (cell_type in cell_type_order) {
    marker_one <- top_markers[CellType==cell_type, .(CellType, MarkerState=State, Gene, WithinStateOrder)]
    state_one <- sort(unique(marker_one$MarkerState))
    grid_one <- CJ(Dataset=dataset_order, Gene=marker_one$Gene, ExpressionState=state_one, unique=TRUE)
    grid_one <- merge(grid_one, marker_one, by="Gene", all.x=TRUE, sort=FALSE)
    full_parts[[cell_type]] <- grid_one
}
expression_grid <- rbindlist(full_parts)
expression_grid <- merge(expression_grid, expression_observed, by=c("Dataset", "CellType", "Gene", "ExpressionState"), all.x=TRUE, sort=FALSE)
recovery_display_flag <- recovery_pass[, .(Dataset, CellType, ExpressionState=State, DisplayEligible=TRUE)]
expression_grid <- merge(expression_grid, recovery_display_flag, by=c("Dataset", "CellType", "ExpressionState"), all.x=TRUE, sort=FALSE)
expression_grid[is.na(DisplayEligible), DisplayEligible := FALSE]
expression_grid[, DisplayExpression := pmax(-0.2, pmin(0.2, MeanUnitVarianceExpression))]
expression_grid[, `:=`(CellTypeOrder=match(CellType, cell_type_order), DatasetOrder=match(Dataset, dataset_order))]
setorder(expression_grid, CellTypeOrder, MarkerState, WithinStateOrder, DatasetOrder, ExpressionState)
fwrite(expression_grid, file.path(out, "Figure2C_marker_expression_source_data.tsv.gz"), sep="\t", quote=FALSE, na="")
fwrite(coverage, file.path(out, "Figure2C_marker_expression_coverage.tsv"), sep="\t", quote=FALSE, na="")
print(coverage[, .(ParentInputs=sum(Status=="INPUT_FOUND"), MedianPresentMarkers=median(PresentMarkers)), by=Dataset])

# Step3. 按原文结构绘制 20 个父类型热图
panel_list <- list()
panel_audit <- list()

for (cell_type in cell_type_order) {
    recovered_state_one <- unique(recovery_pass[CellType==cell_type, State])
    marker_one <- unique(top_markers[CellType==cell_type & State %in% recovered_state_one, .(MarkerState=State, Gene, WithinStateOrder)])
    marker_one[, StateNumber := as.integer(sub("^S0?", "", MarkerState))]
    setorder(marker_one, StateNumber, WithinStateOrder, Gene)
    marker_one[, RowOrder := seq_len(.N)]
    marker_one[, Y := .N-RowOrder+1L]
    state_one <- sort(unique(marker_one$MarkerState))

    column_parts <- list()
    for (state_index in seq_along(state_one)) {
        state <- state_one[[state_index]]
        dataset_one <- recovery_pass[CellType==cell_type & State==state, Dataset]
        dataset_one <- dataset_order[dataset_order %in% dataset_one]
        column_parts[[state]] <- data.table(Dataset=dataset_one, ExpressionState=state)
    }
    columns <- rbindlist(column_parts)
    columns[, X := seq_len(.N)]

    plot_data <- expression_grid[CellType==cell_type & MarkerState %in% state_one]
    plot_data <- merge(plot_data, columns, by=c("Dataset", "ExpressionState"), all=FALSE, sort=FALSE)
    plot_data <- merge(plot_data, marker_one[, .(Gene, MarkerState, Y)], by=c("Gene", "MarkerState"), all.x=TRUE, sort=FALSE)

    row_annotation <- unique(marker_one[, .(MarkerState, Y)])
    row_annotation[, Colour := state_colours[MarkerState]]
    column_annotation <- copy(columns)
    column_annotation[, Y := nrow(marker_one)+0.82]
    column_annotation[, Colour := dataset_colours[Dataset]]
    label_data <- marker_one[, .(Gene, Y)]

    row_blocks <- marker_one[, .(YMin=min(Y)-0.5, YMax=max(Y)+0.5), by=MarkerState]
    diagonal_boxes <- columns[, .(XMin=min(X)-0.5, XMax=max(X)+0.5), by=ExpressionState]
    diagonal_boxes <- merge(diagonal_boxes, row_blocks, by.x="ExpressionState", by.y="MarkerState", all=FALSE, sort=FALSE)

    panel_list[[cell_type]] <- ggplot(plot_data, aes(x=X, y=Y, fill=DisplayExpression)) +
        geom_tile(width=0.98, height=0.98) +
        geom_tile(data=row_annotation, aes(x=0.30, y=Y), inherit.aes=FALSE, width=0.62, height=0.98, fill=row_annotation$Colour) +
        geom_tile(data=column_annotation, aes(x=X, y=Y), inherit.aes=FALSE, width=0.98, height=0.38, fill=column_annotation$Colour) +
        geom_rect(data=diagonal_boxes, aes(xmin=XMin, xmax=XMax, ymin=YMin, ymax=YMax), inherit.aes=FALSE, fill=NA, colour="white", linewidth=0.32) +
        geom_text(data=label_data, aes(x=-0.22, y=Y, label=Gene), inherit.aes=FALSE, family=font_family, fontface="italic", size=1.80, hjust=1) +
        scale_fill_gradientn(colours=c("#352A87", "#2448B8", "#2C7FB8", "#49A65A", "#FDE725"), values=c(0, 0.25, 0.5, 0.75, 1), limits=c(-0.2, 0.2), na.value="#D9D9D9", guide="none") +
        coord_cartesian(xlim=c(-4.9, max(columns$X)+0.55), ylim=c(0.45, nrow(marker_one)+1.15), expand=FALSE, clip="off") +
        labs(title=unname(cell_type_title[cell_type])) +
        theme_void(base_family=font_family, base_size=7) +
        theme(plot.title=element_text(size=8.0, face="plain", hjust=0.58, margin=margin(b=1.5)), plot.margin=margin(1.0, 1.0, 1.0, 1.0))

    panel_audit[[cell_type]] <- data.table(CellType=cell_type, DisplayMarkers=nrow(marker_one), States=length(state_one), HeatmapColumns=nrow(columns), MissingCells=sum(is.na(plot_data$MeanUnitVarianceExpression)), TotalCells=nrow(plot_data))
}

panel_audit <- rbindlist(panel_audit)
fwrite(panel_audit, file.path(out, "Figure2C_panel_audit.tsv"), sep="\t", quote=FALSE)
print(panel_audit)

# Step4. 构建原文式标题和共享图例
title_plot <- ggplot() +
    annotate("text", x=0.02, y=0.5, label="C", family=font_family, fontface="bold", size=4.2, hjust=0) +
    annotate("text", x=0.105, y=0.5, label="Marker expression across scRNA-seq validation datasets", family=font_family, size=3.45, hjust=0) +
    coord_cartesian(xlim=c(0, 1), ylim=c(0, 1), expand=FALSE, clip="off") +
    theme_void()

legend_plot <- ggplot() +
    annotate("text", x=0.01, y=0.74, label="Cell state", family=font_family, size=2.65, hjust=0) +
    annotate("text", x=0.01, y=0.28, label="scRNA\nsequencing dataset", family=font_family, size=2.65, hjust=0, lineheight=0.9) +
    coord_cartesian(xlim=c(0, 1), ylim=c(0, 1), expand=FALSE, clip="off") +
    theme_void()

state_legend <- data.table(X=seq(0.15, 0.52, length.out=6L), State=names(state_colours), Colour=unname(state_colours), Label=sub("^S0", "S", names(state_colours)))
dataset_legend <- data.table(X=c(0.20, 0.36, 0.51, 0.65), Dataset=dataset_order, Colour=unname(dataset_colours), Label=c("GSE214695", "SCP1884", "SCP259", "IBDverse"))
gradient_legend <- data.table(X=seq(0.80, 0.97, length.out=101L), Colour=colorRampPalette(c("#352A87", "#2448B8", "#2C7FB8", "#49A65A", "#FDE725"))(101L))

legend_plot <- legend_plot +
    geom_point(data=state_legend, aes(x=X, y=0.74), inherit.aes=FALSE, shape=22, size=3.1, stroke=0.2, colour="#777777", fill=state_legend$Colour) +
    geom_text(data=state_legend, aes(x=X+0.017, y=0.74, label=Label), inherit.aes=FALSE, family=font_family, size=2.35, hjust=0) +
    geom_point(data=dataset_legend, aes(x=X, y=0.28), inherit.aes=FALSE, shape=22, size=3.1, stroke=0, fill=dataset_legend$Colour) +
    geom_text(data=dataset_legend, aes(x=X+0.017, y=0.28, label=Label), inherit.aes=FALSE, family=font_family, size=2.20, hjust=0) +
    geom_tile(data=gradient_legend, aes(x=X, y=0.32), inherit.aes=FALSE, width=0.0022, height=0.18, fill=gradient_legend$Colour) +
    annotate("text", x=0.885, y=0.76, label="Relative log2\nexpression", family=font_family, size=2.40, hjust=0.5, lineheight=0.9) +
    annotate("text", x=0.80, y=0.11, label="-0.2", family=font_family, size=2.15, hjust=0.5) +
    annotate("text", x=0.97, y=0.11, label="0.2", family=font_family, size=2.15, hjust=0.5)

panel_grid <- patchwork::wrap_plots(panel_list, ncol=4)
figure <- title_plot / panel_grid / legend_plot + patchwork::plot_layout(heights=c(0.045, 1, 0.075))

# Step5. 审计面板对齐并导出多格式终稿
width_in <- 183/25.4
height_in <- 255/25.4
require_patchwork_panel_alignment(
    panel_grid,
    manifest_path=file.path(out, "Figure2C_marker_validation_Cell2021_style.alignment-layout.json"),
    report_path=file.path(out, "Figure2C_marker_validation_Cell2021_style.alignment.json"),
    overlay_svg=file.path(out, "Figure2C_marker_validation_Cell2021_style.alignment.svg"),
    width_in=width_in,
    height_in=height_in/1.12,
    panel_ids=cell_type_order,
    row_groups=unname(split(cell_type_order, rep(seq_len(5), each=4))),
    column_groups=unname(lapply(seq_len(4), function(index) cell_type_order[seq(index, 20, by=4)])),
    audit_script=file.path(run_root, "scripts/figure_qa/audit_panel_alignment.py"),
    python=Sys.which("python3"),
    tolerance_pt=1.5,
    gutter_tolerance_pt=1.5,
    strict=TRUE
)

ragg::agg_png(file.path(out, "Figure2C_marker_validation_Cell2021_style.png"), width=width_in, height=height_in, units="in", res=600, background="white")
print(figure)
dev.off()
ragg::agg_png(file.path(out, "Figure2C_marker_validation_Cell2021_style_preview.png"), width=width_in, height=height_in, units="in", res=300, background="white")
print(figure)
dev.off()
grDevices::pdf(file.path(out, "Figure2C_marker_validation_Cell2021_style.pdf"), width=width_in, height=height_in, family=font_family, useDingbats=FALSE)
print(figure)
dev.off()
grDevices::cairo_pdf(file.path(out, "Figure2C_marker_validation_Cell2021_style_editable.pdf"), width=width_in, height=height_in, family=font_family)
print(figure)
dev.off()
svglite::svglite(file.path(out, "Figure2C_marker_validation_Cell2021_style.svg"), width=width_in, height=height_in, bg="white")
print(figure)
dev.off()
ragg::agg_tiff(file.path(out, "Figure2C_marker_validation_Cell2021_style.tiff"), width=width_in, height=height_in, units="in", res=600, compression="lzw", background="white")
print(figure)
dev.off()

# Step6. 保存方法边界和完成记录
excluded_states <- recovery_state_audit[RecoveredDatasets==0L, paste(CellType, State, sep="_")]
displayed_state_ids <- unique(recovery_pass[, paste(CellType, State, sep="_")])
method_notes <- c(
    "Figure role: descriptive Cell 2021-style display of cross-dataset marker expression; leakage-free quantitative evidence remains the whole-dataset LODO concordance.",
    "Markers: first two HeldOut=All, Disease composite-ranked markers per displayed formal state.",
    "Display eligibility: a dataset-state column is shown only when official Cell 2021 recovery significance is Z > 1.96 in that dataset; columns are ordered state then dataset to form the original staircase.",
    paste0("States with no significant recovery in any validation dataset are omitted from this marker-validation display: ", paste(excluded_states, collapse=", "), "."),
    "Expression: raw counts -> log2(CPM+1) -> per-gene centering and unit-variance scaling within each held-out parent-cell dataset -> assigned-state mean.",
    "Expression eligibility: displayed recovered state means require at least five assigned cells; unavailable genes remain NA and are shown in grey.",
    "Display: values are clipped to [-0.2, 0.2] only for the color scale; unclipped means are retained in source data.",
    "White rectangles: marker-state by matching expression-state diagonal blocks."
)
writeLines(method_notes, file.path(out, "Figure2C_marker_validation_method_notes.txt"))
completion <- list(status="PASS_FIGURE2C_MARKER_VALIDATION", generated_at=format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), backend="R", datasets=dataset_order, parent_cell_types=length(cell_type_order), formal_states=uniqueN(top_markers$CellState), displayed_states=length(displayed_state_ids), displayed_dataset_state_columns=nrow(recovery_pass), displayed_markers=sum(top_markers$CellState %in% displayed_state_ids), state_mean_min_cells=5L)
write_json(completion, file.path(out, "Figure2C_marker_validation_complete.json"), pretty=TRUE, auto_unbox=TRUE)
capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
print(completion)
cat("PASS_FIGURE2C_MARKER_VALIDATION\n")
