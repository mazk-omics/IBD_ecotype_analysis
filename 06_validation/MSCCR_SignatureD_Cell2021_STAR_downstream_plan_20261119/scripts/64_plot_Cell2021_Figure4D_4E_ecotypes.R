# 按 Cell 2021 Figure 4D-4E 版式绘制 discovery bulk、pooled scRNA 共现与代表样本组成。
# 在服务器运行；读取官方 discovery Jaccard、61 号 pooled C/F1 和冻结成员，不重新训练模型。
# Figure 4E 按原文图注展示各样本中 abundance 最高的 EcoType，按 pooled scRNA Z 预先排序展示前三个显著 EcoType。
library(data.table)
library(ggplot2)
library(patchwork)
library(ragg)
library(svglite)
font_family <- "Helvetica"

# Step1. 读入冻结成员、bulk Jaccard、pooled scRNA C/F1 和显著性
args <- commandArgs(trailingOnly=TRUE)
run_root <- if (length(args)>=1) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
input <- file.path(run_root, "output/48_pooled_scRNA_ecotype_cooccurrence_Cell2021")
out <- file.path(run_root, "output/52_Cell2021_Figure4D_4E_ecotypes")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
source(file.path(run_root, "scripts/figure_qa/panel_alignment.R"))

members <- fread(file.path(run_root, "output/00_freeze/formal_model/ecotype_members_frozen.tsv"), data.table=FALSE)
results <- fread(file.path(input, "pooled_four_atlas_ecotype_cooccurrence_8tests.tsv"), data.table=FALSE)
bulk_path <- "/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101/Custom/Ecotypes/discovery/jaccard_matrix.txt"
bulk_C <- as.matrix(read.delim(bulk_path, row.names=1, check.names=FALSE))
scrna_C <- as.matrix(read.delim(file.path(input, "pooled_four_atlas_C_four_matrix_mean.tsv"), row.names=1, check.names=FALSE))
scrna_diagonal <- as.matrix(read.delim(file.path(input, "pooled_four_atlas_C_display_diagonal.tsv"), row.names=1, check.names=FALSE))
F1 <- as.matrix(read.delim(file.path(input, "pooled_four_atlas_F1_assigned_denominator.tsv"), row.names=1, check.names=FALSE))
ecotypes <- paste0("E", 1:8)
ecotype_colours <- c(E1="#A62B49", E2="#D8A51F", E3="#A8C93D", E4="#C0639B", E5="#8064A2", E6="#8B8B8B", E7="#9CB8CC", E8="#CA6F2D")
stopifnot(nrow(members)==79, setequal(members$ID, rownames(scrna_C)), all(members$ID %in% rownames(bulk_C)), setequal(rownames(F1), members$ID), nrow(results)==8)

members$Ecotype <- factor(members$Ecotype, levels=ecotypes)
members <- members[order(members$Ecotype, members$CellType, members$State), ]
bulk_C <- bulk_C[members$ID, members$ID, drop=FALSE]
scrna_C <- scrna_C[members$ID, members$ID, drop=FALSE]
scrna_diagonal <- scrna_diagonal[members$ID, members$ID, drop=FALSE]
F1 <- F1[members$ID, , drop=FALSE]
diag(scrna_C) <- diag(scrna_diagonal)
stopifnot(all(abs(diag(bulk_C)-1)<1e-12), all(abs(diag(scrna_C)-1)<1e-12))

# Step2. 组装 Figure 4D 的两个同序矩阵
counts <- table(factor(members$Ecotype, levels=ecotypes))
starts <- cumsum(c(1, head(as.integer(counts), -1)))
blocks <- data.frame(Ecotype=ecotypes, Start=starts, End=starts+as.integer(counts)-1, Count=as.integer(counts))
blocks$StartNorm <- (blocks$Start-1)/nrow(members)
blocks$EndNorm <- blocks$End/nrow(members)
blocks$CenterNorm <- ((blocks$Start+blocks$End)/2-0.5)/nrow(members)
blocks$Colour <- unname(ecotype_colours[blocks$Ecotype])

bulk_display <- pmin(pmax(bulk_C, 0), 1)
scrna_display <- pmin(pmax(scrna_C, 0), 1)
heatmap_tables <- list()

for (panel in c("Discovery bulk RNA-seq\n(n = 2,025 samples)", "scRNA-seq\n(n = 299 samples)")) {
    display <- if (grepl("Discovery", panel)) bulk_display else scrna_display
    long <- as.data.frame(as.table(display), stringsAsFactors=FALSE)
    names(long) <- c("RowState", "ColState", "Cooccurrence")
    long$Row <- match(long$RowState, members$ID)
    long$Col <- match(long$ColState, members$ID)
    long$RowNorm <- (long$Row-0.5)/nrow(members)
    long$ColNorm <- (long$Col-0.5)/nrow(members)
    long$CellWidth <- 1/nrow(members)
    long$DisplayValue <- sqrt(long$Cooccurrence)
    long$Panel <- panel
    heatmap_tables[[panel]] <- long
}

heatmap_table <- rbindlist(heatmap_tables)
panel_levels <- names(heatmap_tables)
heatmap_table$Panel <- factor(heatmap_table$Panel, levels=panel_levels)
block_table <- rbindlist(lapply(panel_levels, function(panel) cbind(blocks, Panel=panel)))
block_table$Panel <- factor(block_table$Panel, levels=panel_levels)

results <- results[match(ecotypes, results$Ecotype), ]
label_table <- copy(block_table)
label_table$Symbol <- ""
label_table$Symbol[label_table$Panel==panel_levels[2]] <- ifelse(results$P_normal_greater<0.05 & results$AvgDif>results$NullMean, "*", "")
label_table$Label <- paste0(label_table$Ecotype, label_table$Symbol)
bottom_labels <- data.frame(Panel=factor(panel_levels, levels=panel_levels), X=0.5, Y=1.075, Label="79 cell states")

fwrite(heatmap_table, file.path(out, "Figure4D_bulk_and_pooled_scRNA_source.tsv.gz"), sep="\t")
fwrite(label_table, file.path(out, "Figure4D_ecotype_annotations.tsv"), sep="\t")

figure_d <- ggplot(heatmap_table, aes(x=ColNorm, y=RowNorm, fill=DisplayValue)) +
    geom_tile(aes(width=CellWidth, height=CellWidth)) +
    geom_rect(data=block_table, aes(xmin=StartNorm, xmax=EndNorm, ymin=StartNorm, ymax=EndNorm), inherit.aes=FALSE, fill=NA, colour="black", linewidth=0.35) +
    geom_rect(data=block_table, aes(xmin=StartNorm, xmax=EndNorm, ymin=-0.060, ymax=0), inherit.aes=FALSE, fill=block_table$Colour, colour=NA) +
    geom_text(data=label_table, aes(x=CenterNorm, y=-0.105, label=Label), inherit.aes=FALSE, size=1.75, family=font_family) +
    geom_text(data=bottom_labels, aes(x=X, y=Y, label=Label), inherit.aes=FALSE, size=2.25, family=font_family) +
    scale_fill_gradientn(colours=c("#FFFFD9", "#D9F0A3", "#78C679", "#2C7FB8", "#253494"), limits=c(0, 1), na.value="#D9D9D9", breaks=c(0, 1), labels=c("Low", "High"), name="Co-occurrence index") +
    facet_wrap(~Panel, nrow=1) +
    scale_y_reverse(limits=c(1.135, -0.150), expand=c(0, 0)) +
    coord_equal(xlim=c(0, 1), expand=FALSE, clip="off") +
    labs(tag="D", title="Cross-correlation of cell state abundances", x=NULL, y="79 cell states", caption="* P < 0.05") +
    theme_classic(base_family=font_family, base_size=6.5) +
    theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), strip.background=element_blank(), strip.text=element_text(size=7.2), axis.title=element_text(size=6.8), legend.position="bottom", legend.title=element_text(size=6.0), legend.text=element_text(size=5.8), plot.tag=element_text(size=10, face="bold"), plot.tag.position=c(0.025, 0.965), plot.title=element_text(size=8.0, face="plain", hjust=0.54, margin=margin(b=2)), plot.title.position="plot", plot.caption=element_text(size=6.0, hjust=1), plot.caption.position="plot", panel.spacing=grid::unit(3.2, "mm"), plot.margin=margin(3, 2, 2, 2)) +
    guides(fill=guide_colourbar(title.position="top", title.hjust=0.5, barwidth=grid::unit(17, "mm"), barheight=grid::unit(2.2, "mm"), ticks=FALSE))

# Step3. 按 STAR Methods 计算 pooled scRNA 的 EcoType abundance 和正式 assignment
membership <- lapply(ecotypes, function(ecotype) which(members$Ecotype==ecotype))
names(membership) <- ecotypes
eco_raw <- matrix(NA_real_, nrow=length(ecotypes), ncol=ncol(F1), dimnames=list(ecotypes, colnames(F1)))
eco_p <- eco_raw

for (sample_index in seq_len(ncol(F1))) {
    sample_values <- F1[, sample_index]
    for (ecotype in ecotypes) {
        member_index <- membership[[ecotype]]
        inside <- sample_values[member_index]
        outside <- sample_values[-member_index]
        eco_raw[ecotype, sample_index] <- mean(inside[is.finite(inside)])
        eco_p[ecotype, sample_index] <- tryCatch(t.test(inside, outside, var.equal=FALSE)$p.value, error=function(error) 1)
    }
}

eco_abundance <- sweep(eco_raw, 2, colSums(eco_raw, na.rm=TRUE), "/")
eco_q <- apply(eco_p, 2, p.adjust, method="BH")
if (is.vector(eco_q)) eco_q <- matrix(eco_q, nrow=length(ecotypes), dimnames=dimnames(eco_p))
rownames(eco_q) <- ecotypes
colnames(eco_q) <- colnames(F1)
max_ecotype <- rep(NA_character_, ncol(eco_abundance))

for (sample_index in seq_len(ncol(eco_abundance))) {
    sample_abundance <- eco_abundance[, sample_index]
    if (any(is.finite(sample_abundance))) max_ecotype[sample_index] <- ecotypes[which.max(replace(sample_abundance, !is.finite(sample_abundance), -Inf))]
}

max_q <- rep(NA_real_, ncol(eco_q))
defined_sample <- !is.na(max_ecotype)
max_q[defined_sample] <- eco_q[cbind(match(max_ecotype[defined_sample], ecotypes), which(defined_sample))]
has_member_state <- colSums(F1>0, na.rm=TRUE)>0
official_ecotype <- ifelse(max_q<0.25 & has_member_state, max_ecotype, NA_character_)

assignment <- data.frame(Sample=colnames(F1), Dataset=sub("::.*$", "", colnames(F1)), MaxEcotype=max_ecotype, MaxEcotypeAbundance=eco_abundance[cbind(match(max_ecotype, ecotypes), seq_len(ncol(eco_abundance)))], AssignmentQ=max_q, HasMemberState=has_member_state, OfficialEcotype=official_ecotype)
fwrite(assignment, file.path(out, "Figure4E_pooled_scRNA_ecotype_assignment.tsv"), sep="\t")
write.table(eco_abundance, file.path(out, "Figure4E_pooled_scRNA_ecotype_abundance.tsv"), sep="\t", quote=FALSE, col.names=NA)
assignment_summary <- data.frame(
    Ecotype=ecotypes,
    PooledZ=results$Z,
    Cell2021Nominal=results$Cell2021Nominal,
    MostAbundantSamples=vapply(ecotypes, function(ecotype) sum(assignment$MaxEcotype==ecotype, na.rm=TRUE), integer(1)),
    FormalAssignedSamples=vapply(ecotypes, function(ecotype) sum(assignment$OfficialEcotype==ecotype, na.rm=TRUE), integer(1))
)
fwrite(assignment_summary, file.path(out, "Figure4E_ecotype_selection_audit.tsv"), sep="\t")

# Step4. 客观选择前三个显著 EcoType，并按原文图注展示各组 abundance 最高的 CE-high 样本
selected_ecotypes <- results$Ecotype[results$Cell2021Nominal][order(results$Z[results$Cell2021Nominal], decreasing=TRUE)][1:3]
selected_samples <- list()

for (ecotype in selected_ecotypes) {
    candidates <- assignment[assignment$MaxEcotype==ecotype & is.finite(assignment$MaxEcotypeAbundance), ]
    candidates <- candidates[order(candidates$MaxEcotypeAbundance, decreasing=TRUE, candidates$Sample), ]
    selected_samples[[ecotype]] <- head(candidates, min(10, nrow(candidates)))
}

selected_sample_table <- rbindlist(lapply(selected_ecotypes, function(ecotype) {
    table <- selected_samples[[ecotype]]
    table$DisplayGroup <- ecotype
    table
}))
stopifnot(nrow(selected_sample_table)>0, all(table(selected_sample_table$DisplayGroup)>0))
selected_sample_table$X <- seq_len(nrow(selected_sample_table))
fwrite(selected_sample_table, file.path(out, "Figure4E_selected_samples.tsv"), sep="\t")

# Step5. 计算选定样本中三个 EcoType 的成员组成，并按每个 EcoType 的最大总量缩放
composition_tables <- list()

for (target_ecotype in selected_ecotypes) {
    member_index <- membership[[target_ecotype]]
    target_members <- members[member_index, ]
    target_values <- F1[member_index, selected_sample_table$Sample, drop=FALSE]/length(member_index)
    target_values[!is.finite(target_values)] <- 0
    row_max <- max(colSums(target_values), na.rm=TRUE)
    if (!is.finite(row_max) || row_max<=0) row_max <- 1
    target_values <- target_values/row_max
    long <- as.data.frame(as.table(target_values), stringsAsFactors=FALSE)
    names(long) <- c("CellState", "Sample", "Contribution")
    long$CellType <- target_members$CellType[match(long$CellState, target_members$ID)]
    long$Ecotype <- target_ecotype
    long$X <- selected_sample_table$X[match(long$Sample, selected_sample_table$Sample)]
    composition_tables[[target_ecotype]] <- long
}

composition <- rbindlist(composition_tables)
composition$Ecotype <- factor(composition$Ecotype, levels=selected_ecotypes)
composition$CellType <- factor(composition$CellType, levels=c("B", "Plasma", "CD8_T", "CD4_T", "NK", "Monocyte_Macrophage", "DC", "Mast", "ILC", "Fibroblast", "Endothelial", "Pericyte", "Glial", "Colonic_absorptive", "Ileal_absorptive", "Goblet", "Paneth", "Tuft", "Enteroendocrine", "Stem_TA_progenitor"))
fwrite(composition, file.path(out, "Figure4E_selected_ecotype_composition_source.tsv"), sep="\t")

celltype_colours <- c(B="#C99A20", Plasma="#B51F6D", CD8_T="#6B5B95", CD4_T="#B55D20", NK="#4C8B65", Monocyte_Macrophage="#9B86B5", DC="#E0B35A", Mast="#D98783", ILC="#65A6A1", Fibroblast="#9ABBD5", Endothelial="#7F8788", Pericyte="#B7A6D1", Glial="#6C8EAD", Colonic_absorptive="#D7659E", Ileal_absorptive="#F28E2B", Goblet="#7FBF7B", Paneth="#D6C04C", Tuft="#5AAEAE", Enteroendocrine="#8C6BB1", Stem_TA_progenitor="#E58B6A")
dataset_colours <- c(GSE214695="#1B1B1B", SCP1884="#6B6B6B", SCP259="#A0A0A0", IBDverse_E_MTAB_16999="#D0D0D0")
celltype_labels <- c(B="B cells", Plasma="PCs", CD8_T="CD8 T", CD4_T="CD4 T", NK="NKs", Monocyte_Macrophage="Mono/Mac", DC="DCs", Mast="MCs", ILC="ILCs", Fibroblast="Fibroblasts", Endothelial="Endothelial", Pericyte="Pericytes", Glial="Glial", Colonic_absorptive="Colonic absorptive", Ileal_absorptive="Ileal absorptive", Goblet="Goblet", Paneth="Paneth", Tuft="Tuft", Enteroendocrine="EECs", Stem_TA_progenitor="Stem/TA")
dataset_labels <- c(GSE214695="GSE214695", SCP1884="SCP1884", SCP259="SCP259", IBDverse_E_MTAB_16999="IBDverse")

group_summary <- selected_sample_table[, .(Start=min(X), End=max(X), Center=mean(range(X)), N=.N), by=DisplayGroup]
group_summary$Label <- paste0(group_summary$DisplayGroup, "\nhigh")
separators <- group_summary$End[-nrow(group_summary)]+0.5
dataset_strip <- selected_sample_table[, c("Sample", "Dataset", "X")]
dataset_strip$Ecotype <- factor(selected_ecotypes[length(selected_ecotypes)], levels=selected_ecotypes)

figure_e <- ggplot(composition, aes(x=X, y=Contribution, fill=CellType)) +
    geom_col(width=0.86, colour="black", linewidth=0.12) +
    geom_vline(xintercept=separators, linewidth=0.35) +
    geom_segment(data=dataset_strip, aes(x=X-0.50, xend=X+0.50, y=-0.035, yend=-0.035, colour=Dataset), inherit.aes=FALSE, linewidth=3, lineend="butt") +
    facet_grid(Ecotype~., switch="y") +
    scale_fill_manual(values=celltype_colours, labels=celltype_labels, drop=TRUE, name="Cell type") +
    scale_colour_manual(values=dataset_colours, labels=dataset_labels, drop=TRUE, name="scRNA-seq dataset") +
    scale_x_continuous(limits=c(0.5, nrow(selected_sample_table)+0.5), expand=c(0, 0), sec.axis=dup_axis(breaks=group_summary$Center, labels=group_summary$Label, name=NULL)) +
    scale_y_continuous(limits=c(-0.08, 1.02), expand=c(0, 0)) +
    labs(tag="E", title="EcoType composition in selected\nsamples profiled by scRNA-seq", x=paste0(nrow(selected_sample_table), " selected scRNA-seq samples"), y=NULL) +
    theme_classic(base_family=font_family, base_size=6.2) +
    theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text.x.bottom=element_blank(), axis.text.y=element_blank(), axis.text.x.top=element_text(size=6.7, lineheight=0.88, margin=margin(b=1.5)), strip.placement="outside", strip.background=element_blank(), strip.text.y.left=element_text(angle=0, size=7.0), panel.spacing=grid::unit(1.1, "mm"), legend.position="right", legend.box="vertical", legend.spacing.y=grid::unit(0.4, "mm"), legend.margin=margin(0, 0, 0, 0), legend.title=element_text(size=5.5), legend.text=element_text(size=5.0), legend.key.height=grid::unit(1.7, "mm"), legend.key.spacing.x=grid::unit(0.7, "mm"), plot.tag=element_text(size=10, face="bold"), plot.tag.position=c(0.025, 0.965), plot.title=element_text(size=7.7, face="plain", hjust=0.58, lineheight=0.92, margin=margin(b=2)), plot.title.position="plot", axis.title.x=element_text(size=6.2, margin=margin(t=4)), plot.margin=margin(3, 2, 2, 2)) +
    guides(fill=guide_legend(order=1, ncol=1, keyheight=grid::unit(1.7, "mm"), keywidth=grid::unit(1.9, "mm")), colour=guide_legend(order=2, override.aes=list(linewidth=3), keyheight=grid::unit(1.7, "mm")))

# Step6. 分别导出 D、E，并组合为原文式 D-E 页面
width_d <- 96/25.4
width_e <- 87/25.4
height_panel <- 82/25.4

ragg::agg_png(file.path(out, "Figure4D_bulk_and_pooled_scRNA.png"), width=width_d, height=height_panel, units="in", res=600, background="white")
print(figure_d)
dev.off()
grDevices::pdf(file.path(out, "Figure4D_bulk_and_pooled_scRNA.pdf"), width=width_d, height=height_panel, family=font_family, useDingbats=FALSE)
print(figure_d)
dev.off()
svglite::svglite(file.path(out, "Figure4D_bulk_and_pooled_scRNA.svg"), width=width_d, height=height_panel, bg="white")
print(figure_d)
dev.off()

ragg::agg_png(file.path(out, "Figure4E_selected_scRNA_ecotype_composition.png"), width=width_e, height=height_panel, units="in", res=600, background="white")
print(figure_e)
dev.off()
grDevices::pdf(file.path(out, "Figure4E_selected_scRNA_ecotype_composition.pdf"), width=width_e, height=height_panel, family=font_family, useDingbats=FALSE)
print(figure_e)
dev.off()
svglite::svglite(file.path(out, "Figure4E_selected_scRNA_ecotype_composition.svg"), width=width_e, height=height_panel, bg="white")
print(figure_e)
dev.off()

combined <- figure_d + figure_e + patchwork::plot_layout(widths=c(1.04, 1))
width_combined <- 183/25.4
height_combined <- 82/25.4
require_patchwork_panel_alignment(
    combined,
    manifest_path=file.path(out, "Figure4D_4E_Cell2021_style.alignment-layout.json"),
    report_path=file.path(out, "Figure4D_4E_Cell2021_style.alignment.json"),
    overlay_svg=file.path(out, "Figure4D_4E_Cell2021_style.alignment.svg"),
    width_in=width_combined,
    height_in=height_combined,
    panel_ids=c("D", "E"),
    row_groups=list(c("D", "E")),
    exemptions=list(list(panels=c("D", "E"), checks="panel-width", reason="Figure D is a two-heatmap hero panel and Figure E is a narrower composition panel.")),
    audit_script=file.path(run_root, "scripts/figure_qa/audit_panel_alignment.py"),
    python=Sys.which("python3"),
    tolerance_pt=1.5,
    gutter_tolerance_pt=1.5,
    strict=TRUE
)
ragg::agg_png(file.path(out, "Figure4D_4E_Cell2021_style.png"), width=width_combined, height=height_combined, units="in", res=600, background="white")
print(combined)
dev.off()
grDevices::pdf(file.path(out, "Figure4D_4E_Cell2021_style.pdf"), width=width_combined, height=height_combined, family=font_family, useDingbats=FALSE)
print(combined)
dev.off()
svglite::svglite(file.path(out, "Figure4D_4E_Cell2021_style.svg"), width=width_combined, height=height_combined, bg="white")
print(combined)
dev.off()
ragg::agg_tiff(file.path(out, "Figure4D_4E_Cell2021_style.tiff"), width=width_combined, height=height_combined, units="in", res=600, compression="lzw", background="white")
print(combined)
dev.off()
grDevices::cairo_pdf(file.path(out, "Figure4D_4E_Cell2021_style.cairo.pdf"), width=width_combined, height=height_combined, family=font_family)
print(combined)
dev.off()

# Step7. 保存方法边界和环境
writeLines(c(
    "Figure 4D bulk: the official 85-state EcoTyper discovery Jaccard matrix is restricted to the 79 frozen EcoType members and ordered identically to pooled scRNA.",
    "Figure 4D scRNA: the four prespecified atlases are pooled before F1/F2 and the four-matrix co-occurrence index are calculated; no dataset-specific heat maps are drawn.",
    "Figure 4D colour mapping uses a square-root display transform to separate low and intermediate co-occurrence values; exported source tables retain the untransformed Cooccurrence column and no statistics are changed.",
    "Figure 4E audit assignment: EcoType abundance is the mean member-state F1 abundance normalized to sum to one; Welch tests and within-sample BH are used, with q<0.25 and at least one positive member state required for the separately exported formal assignment.",
    paste0("Figure 4E display: ", paste(selected_ecotypes, collapse=", "), " are the three significant EcoTypes with the largest prespecified pooled-scRNA Z values; following the Figure 4E caption, up to 10 samples in which each displayed EcoType is most abundant are shown, without adding the formal q<0.25 filter."),
    "Figure 4E bar segments are member-state F1 contributions coloured by parent cell type and scaled within each displayed EcoType to preserve between-sample height differences.",
    "Selection affects display only; all 8 EcoType abundances and assignments are exported and the formal recovery tests are unchanged."
), file.path(out, "Figure4D_4E_method_boundary.txt"))
capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
