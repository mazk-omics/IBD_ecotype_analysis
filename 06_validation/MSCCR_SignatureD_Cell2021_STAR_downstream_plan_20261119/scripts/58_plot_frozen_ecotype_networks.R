# 绘制本轮冻结 79-state/8-EcoType 的成员网络，版式参照 Cell 2021 Figure 4B。
# 在服务器运行；只读正式成员与发现队列 Jaccard，不修改模型或重新聚类。
# 边表示 bulk 样本共现 Jaccard，不表示配体受体通讯或空间接触。
library(data.table)
library(igraph)
library(ggplot2)
library(patchwork)
library(ragg)
library(svglite)
font_family <- "Helvetica"

# Step1. 定位正式输入和新输出
args <- commandArgs(trailingOnly=TRUE)
run_root <- if (length(args)>=1) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
out <- file.path(run_root, "output/45_frozen_ecotype_networks")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
source(file.path(run_root, "scripts/figure_qa/panel_alignment.R"))

members_path <- file.path(run_root, "output/00_freeze/formal_model/ecotype_members_frozen.tsv")
jaccard_path <- "/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101/Custom/Ecotypes/discovery/jaccard_matrix.txt"
members <- fread(members_path, data.table=FALSE)
jaccard <- as.matrix(read.delim(jaccard_path, row.names=1, check.names=FALSE))
ecotypes <- paste0("E", 1:8)

stopifnot(nrow(members)==79, setequal(unique(members$Ecotype), ecotypes))
stopifnot(all(members$ID %in% rownames(jaccard)), all(members$ID %in% colnames(jaccard)))
stopifnot(max(abs(jaccard- t(jaccard)), na.rm=TRUE)<1e-12)
print(table(members$Ecotype))

# Step2. 固定父类型缩写和 state 颜色
parent_label <- c(B="B", CD4_T="CD4", CD8_T="CD8", Colonic_absorptive="Colon", DC="DC", Endothelial="Endo", Enteroendocrine="EEC", Fibroblast="Fib", Glial="Glia", Goblet="Gob", ILC="ILC", Ileal_absorptive="Ileum", Mast="Mast", Monocyte_Macrophage="Mac", NK="NK", Paneth="Pan", Pericyte="Peri", Plasma="PC", Stem_TA_progenitor="Stem", Tuft="Tuft")
state_colours <- c(S01="#C53E73", S02="#5796B7", S03="#6FA64A", S04="#6C50A3", S05="#D46A1D", S06="#E3AE1A", S07="#8A4B27", S08="#C764D9", S09="#808080", S10="#2C9A8A", S11="#B24C3D")
stopifnot(all(members$CellType %in% names(parent_label)), all(members$State %in% names(state_colours)))
used_states <- names(state_colours)[names(state_colours) %in% sort(unique(members$State))]
print(table(members$State))

# Step3. 每个 EcoType 计算固定 FR 坐标和真实边
set.seed(20261125)
node_tables <- list()
edge_tables <- list()
edge_audit <- list()

for (ecotype in ecotypes) {
    own_members <- members[members$Ecotype==ecotype, ]
    weights <- jaccard[own_members$ID, own_members$ID, drop=FALSE]
    diag(weights) <- 0
    graph <- graph_from_adjacency_matrix(weights, mode="undirected", weighted=TRUE, diag=FALSE)
    set.seed(20261125)
    coordinates <- layout_with_fr(graph, weights=E(graph)$weight)
    nodes <- data.frame(Ecotype=ecotype, ID=V(graph)$name, X=coordinates[,1], Y=coordinates[,2])
    nodes <- merge(nodes, own_members[, c("ID", "CellType", "State")], by="ID", sort=FALSE)
    nodes$X <- nodes$X-mean(nodes$X)
    nodes$Y <- nodes$Y-mean(nodes$Y)
    span <- max(diff(range(nodes$X)), diff(range(nodes$Y)))
    nodes$X <- nodes$X/span*2
    nodes$Y <- nodes$Y/span*2
    nodes$Label <- unname(parent_label[nodes$CellType])

    pairs <- which(upper.tri(weights), arr.ind=TRUE)
    edges <- data.frame(Ecotype=ecotype, From=rownames(weights)[pairs[,1]], To=colnames(weights)[pairs[,2]], Jaccard=weights[pairs])
    edges <- edges[edges$Jaccard>0, ]
    edges$X1 <- nodes$X[match(edges$From, nodes$ID)]
    edges$Y1 <- nodes$Y[match(edges$From, nodes$ID)]
    edges$X2 <- nodes$X[match(edges$To, nodes$ID)]
    edges$Y2 <- nodes$Y[match(edges$To, nodes$ID)]

    node_tables[[ecotype]] <- nodes
    edge_tables[[ecotype]] <- edges
    edge_audit[[ecotype]] <- data.frame(Ecotype=ecotype, Members=nrow(nodes), PossibleEdges=nrow(pairs), PositiveEdges=nrow(edges), MedianPositiveJaccard=median(edges$Jaccard), MaxJaccard=max(edges$Jaccard))

}

nodes_all <- rbindlist(node_tables)
edges_all <- rbindlist(edge_tables)
audit_all <- rbindlist(edge_audit)
fwrite(nodes_all, file.path(out, "frozen_ecotype_network_nodes.tsv"), sep="\t")
fwrite(edges_all, file.path(out, "frozen_ecotype_network_positive_edges.tsv"), sep="\t")
fwrite(audit_all, file.path(out, "frozen_ecotype_network_audit.tsv"), sep="\t")
print(audit_all)

# Step4. 按原文版式输出 4x2 独立网络面板
nodes_all$Ecotype <- factor(nodes_all$Ecotype, levels=ecotypes)
edges_all$Ecotype <- factor(edges_all$Ecotype, levels=ecotypes)
nodes_all$State <- factor(nodes_all$State, levels=used_states)
state_labels <- sub("^S0", "S", used_states)
network_panels <- list()

for (ecotype in ecotypes) {
    nodes_one <- nodes_all[nodes_all$Ecotype==ecotype, ]
    edges_one <- edges_all[edges_all$Ecotype==ecotype, ]
    network_panels[[ecotype]] <- ggplot() +
        geom_segment(data=edges_one, aes(x=X1, y=Y1, xend=X2, yend=Y2, linewidth=Jaccard), colour="black", alpha=0.78, lineend="round") +
        scale_linewidth_continuous(range=c(0.08, 0.72), guide="none") +
        geom_point(data=nodes_one, aes(x=X, y=Y, fill=State), shape=21, size=4.4, stroke=0.25, colour="#555555") +
        geom_text(data=nodes_one, aes(x=X, y=Y, label=Label), size=1.8, family=font_family, lineheight=0.9) +
        scale_fill_manual(values=state_colours[used_states], breaks=used_states, labels=state_labels, drop=FALSE, guide="none") +
        coord_equal(xlim=c(-1.32, 1.32), ylim=c(-1.26, 1.26), expand=FALSE, clip="off") +
        labs(title=ecotype) +
        theme_void(base_family=font_family, base_size=7) +
        theme(plot.title=element_text(size=8.0, face="plain", hjust=0.5, margin=margin(b=1.5)), plot.margin=margin(0.5, 0.5, 0.5, 0.5))
}

title_plot <- ggplot() +
    annotate("text", x=0.02, y=0.5, label="B", family=font_family, fontface="bold", size=4.0, hjust=0) +
    annotate("text", x=0.08, y=0.5, label="Co-occurrence of cell states in IBD ecotypes", family=font_family, size=3.1, hjust=0) +
    coord_cartesian(xlim=c(0, 1), ylim=c(0, 1), expand=FALSE, clip="off") +
    theme_void()

legend_data <- data.frame(X=seq(0.32, 0.78, length.out=length(used_states)), State=factor(used_states, levels=used_states), Label=state_labels)
legend_plot <- ggplot(legend_data, aes(x=X, y=0.5)) +
    annotate("text", x=0.01, y=0.5, label="Cell state", family=font_family, size=2.7, hjust=0) +
    geom_point(aes(fill=State), shape=22, size=3.1, stroke=0.25, colour="#777777") +
    geom_text(aes(x=X+0.017, label=Label), family=font_family, size=2.35, hjust=0) +
    scale_fill_manual(values=state_colours[used_states], guide="none") +
    coord_cartesian(xlim=c(0, 1), ylim=c(0, 1), expand=FALSE, clip="off") +
    theme_void()

network_grid <- patchwork::wrap_plots(network_panels, ncol=4)
figure_layout <- title_plot / network_grid / legend_plot + patchwork::plot_layout(heights=c(0.08, 1, 0.085))

width_in <- 183/25.4
height_in <- 119/25.4
require_patchwork_panel_alignment(
    network_grid,
    manifest_path=file.path(out, "frozen_ecotype_networks_Cell2021_style.alignment-layout.json"),
    report_path=file.path(out, "frozen_ecotype_networks_Cell2021_style.alignment.json"),
    overlay_svg=file.path(out, "frozen_ecotype_networks_Cell2021_style.alignment.svg"),
    width_in=width_in,
    height_in=height_in/1.165,
    panel_ids=ecotypes,
    row_groups=list(ecotypes[1:4], ecotypes[5:8]),
    column_groups=lapply(seq_len(4), function(index) c(ecotypes[index], ecotypes[index+4])),
    audit_script=file.path(run_root, "scripts/figure_qa/audit_panel_alignment.py"),
    python=Sys.which("python3"),
    tolerance_pt=1.5,
    gutter_tolerance_pt=1.5,
    strict=TRUE
)
ragg::agg_png(file.path(out, "frozen_ecotype_networks_Cell2021_style.png"), width=width_in, height=height_in, units="in", res=600, background="white")
print(figure_layout)
dev.off()
ragg::agg_png(file.path(out, "frozen_ecotype_networks_Cell2021_style_preview.png"), width=width_in, height=height_in, units="in", res=300, background="white")
print(figure_layout)
dev.off()
grDevices::cairo_pdf(file.path(out, "frozen_ecotype_networks_Cell2021_style.pdf"), width=width_in, height=height_in, family=font_family)
print(figure_layout)
dev.off()
svglite::svglite(file.path(out, "frozen_ecotype_networks_Cell2021_style.svg"), width=width_in, height=height_in, bg="white")
print(figure_layout)
dev.off()
ragg::agg_tiff(file.path(out, "frozen_ecotype_networks_Cell2021_style.tiff"), width=width_in, height=height_in, units="in", res=600, compression="lzw", background="white")
print(figure_layout)
dev.off()

# Step5. 保存方法边界和环境
writeLines(c("Model: frozen 79 member states in 8 EcoTypes; no reclustering.", "Layout: deterministic Fruchterman-Reingold within each EcoType; display only.", "Edge: discovery bulk sample co-occurrence Jaccard > 0.", "Boundary: the network is not ligand-receptor communication or spatial contact."), file.path(out, "network_method_boundary.txt"))
capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
