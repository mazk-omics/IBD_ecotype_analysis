# 绘制本轮冻结 79-state/8-EcoType 的成员网络，版式参照 Cell 2021 Figure 4B。
# 在服务器运行；只读正式成员与发现队列 Jaccard，不修改模型或重新聚类。
# 边表示 bulk 样本共现 Jaccard，不表示配体受体通讯或空间接触。
library(data.table)
library(igraph)
library(ggplot2)
library(ragg)
library(svglite)

# Step1. 定位正式输入和新输出
args <- commandArgs(trailingOnly=TRUE)
run_root <- if (length(args)>=1) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
out <- file.path(run_root, "output/45_frozen_ecotype_networks")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

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
parent_label <- c(B="B", CD4_T="CD4", CD8_T="CD8", Colonic_absorptive="Colon", DC="DC", Endothelial="Endo", Enteroendocrine="EEC", Fibroblast="Fib", Glial="Glia", Goblet="Gob", ILC="ILC", Ileal_absorptive="Ileum", Mast="Mast", Monocyte_Macrophage="Mac", NK="NK", Paneth="Pan", Pericyte="Peri", Plasma="Plasma", Stem_TA_progenitor="Stem", Tuft="Tuft")
state_colours <- c(S01="#C53E73", S02="#5796B7", S03="#6FA64A", S04="#6C50A3", S05="#D46A1D", S06="#E3AE1A", S07="#8A4B27", S08="#C764D9", S09="#808080", S10="#2C9A8A", S11="#B24C3D")
stopifnot(all(members$CellType %in% names(parent_label)), all(members$State %in% names(state_colours)))

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

# Step4. 用单一 facet 画布输出 4x2 网络，保证各面板坐标和尺寸一致
nodes_all$Ecotype <- factor(nodes_all$Ecotype, levels=ecotypes)
edges_all$Ecotype <- factor(edges_all$Ecotype, levels=ecotypes)
nodes_all$State <- factor(nodes_all$State, levels=names(state_colours))
figure <- ggplot() +
    geom_segment(data=edges_all, aes(x=X1, y=Y1, xend=X2, yend=Y2, linewidth=Jaccard), colour="#777777", alpha=0.72, lineend="round") +
    scale_linewidth_continuous(range=c(0.15, 1.2), guide="none") +
    geom_point(data=nodes_all, aes(x=X, y=Y, fill=State), shape=21, size=6.5, stroke=0.35, colour="white") +
    geom_text(data=nodes_all, aes(x=X, y=Y, label=Label), size=1.95, family="sans", lineheight=0.9) +
    scale_fill_manual(values=state_colours, drop=FALSE, name="Cell state", guide=guide_legend(nrow=1, byrow=TRUE, override.aes=list(shape=22, size=3.2, colour="#666666"))) +
    facet_wrap(~Ecotype, nrow=2) +
    coord_equal(xlim=c(-1.35, 1.35), ylim=c(-1.30, 1.30), expand=FALSE, clip="off") +
    labs(title="Co-occurrence of cell states in IBD EcoTypes") +
    theme_void(base_family="sans", base_size=7) +
    theme(strip.text=element_text(size=8.5), plot.title=element_text(size=11, face="bold", hjust=0.02), legend.position="bottom", legend.title=element_text(size=6.5), legend.text=element_text(size=6), legend.spacing.x=grid::unit(1.2, "pt"), panel.spacing=grid::unit(2.5, "mm"), plot.margin=margin(3, 3, 2, 3))

width_in <- 183/25.4
height_in <- 128/25.4
ragg::agg_png(file.path(out, "frozen_ecotype_networks_Cell2021_style.png"), width=width_in, height=height_in, units="in", res=600, background="white")
print(figure)
dev.off()
grDevices::cairo_pdf(file.path(out, "frozen_ecotype_networks_Cell2021_style.pdf"), width=width_in, height=height_in, family="sans")
print(figure)
dev.off()
svglite::svglite(file.path(out, "frozen_ecotype_networks_Cell2021_style.svg"), width=width_in, height=height_in, bg="white")
print(figure)
dev.off()

# Step5. 保存方法边界和环境
writeLines(c("Model: frozen 79 member states in 8 EcoTypes; no reclustering.", "Layout: deterministic Fruchterman-Reingold within each EcoType; display only.", "Edge: discovery bulk sample co-occurrence Jaccard > 0.", "Boundary: the network is not ligand-receptor communication or spatial contact."), file.path(out, "network_method_boundary.txt"))
capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
