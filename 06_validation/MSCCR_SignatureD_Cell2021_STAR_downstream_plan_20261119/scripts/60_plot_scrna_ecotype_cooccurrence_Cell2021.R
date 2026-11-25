# 将 pooled 单细胞 EcoType 共现结果绘制为一张 Cell 2021 Figure 4D 风格热图。
# 在服务器运行；只读取 61 号脚本的 pooled C 矩阵、显著性和冻结成员，不重新计算检验。
# 热图显示非负共现强度；统计始终使用带符号且对角线为 NA 的原始 C。
library(data.table)
library(ggplot2)
library(ragg)
library(svglite)

# Step1. 读入 pooled 结果和冻结成员
args <- commandArgs(trailingOnly=TRUE)
run_root <- if (length(args)>=1) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
input <- file.path(run_root, "output/48_pooled_scRNA_ecotype_cooccurrence_Cell2021")
out <- file.path(run_root, "output/49_pooled_scRNA_ecotype_cooccurrence_figure")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

members <- fread(file.path(run_root, "output/00_freeze/formal_model/ecotype_members_frozen.tsv"), data.table=FALSE)
results <- fread(file.path(input, "pooled_four_atlas_ecotype_cooccurrence_8tests.tsv"), data.table=FALSE)
C <- as.matrix(read.delim(file.path(input, "pooled_four_atlas_C_four_matrix_mean.tsv"), row.names=1, check.names=FALSE))
C_diagonal <- as.matrix(read.delim(file.path(input, "pooled_four_atlas_C_display_diagonal.tsv"), row.names=1, check.names=FALSE))
ecotypes <- paste0("E", 1:8)
ecotype_colours <- c(E1="#9E2946", E2="#D49B20", E3="#62A066", E4="#3E9A9A", E5="#4E78A8", E6="#775DA6", E7="#C36B32", E8="#C15D9A")
stopifnot(nrow(members)==79, nrow(results)==8, setequal(rownames(C), members$ID), setequal(results$Ecotype, ecotypes))

# Step2. 按冻结 EcoType、父类型和 state 排列单一 pooled C 矩阵
members$Ecotype <- factor(members$Ecotype, levels=ecotypes)
members <- members[order(members$Ecotype, members$CellType, members$State), ]
C <- C[members$ID, members$ID, drop=FALSE]
C_diagonal <- C_diagonal[members$ID, members$ID, drop=FALSE]
diag(C) <- diag(C_diagonal)
display <- pmin(pmax(C, 0), 1)

heatmap_table <- as.data.frame(as.table(display), stringsAsFactors=FALSE)
names(heatmap_table) <- c("RowState", "ColState", "Cooccurrence")
heatmap_table$Row <- match(heatmap_table$RowState, members$ID)
heatmap_table$Col <- match(heatmap_table$ColState, members$ID)

counts <- table(factor(members$Ecotype, levels=ecotypes))
starts <- cumsum(c(1, head(as.integer(counts), -1)))
blocks <- data.frame(Ecotype=ecotypes, Start=starts, End=starts+as.integer(counts)-1, Count=as.integer(counts))
blocks$Center <- (blocks$Start+blocks$End)/2
blocks$Colour <- unname(ecotype_colours[blocks$Ecotype])

results <- results[match(ecotypes, results$Ecotype), ]
blocks$P_normal_greater <- results$P_normal_greater
blocks$P_empirical_greater <- results$P_empirical_greater
blocks$Z <- results$Z
blocks$Symbol <- ifelse(results$P_normal_greater<0.05 & results$AvgDif>results$NullMean, "*", "")

fwrite(heatmap_table, file.path(out, "Figure4D_pooled_scRNA_display_matrix.tsv.gz"), sep="\t")
fwrite(blocks, file.path(out, "Figure4D_pooled_scRNA_ecotype_annotations.tsv"), sep="\t")

# Step3. 绘制一张 pooled Figure 4D 风格热图
figure_4d <- ggplot(heatmap_table, aes(x=Col, y=Row, fill=Cooccurrence)) +
    geom_tile() +
    geom_rect(data=blocks, aes(xmin=Start-0.5, xmax=End+0.5, ymin=Start-0.5, ymax=End+0.5), inherit.aes=FALSE, fill=NA, colour="black", linewidth=0.35) +
    geom_rect(data=blocks, aes(xmin=Start-0.5, xmax=End+0.5, ymin=0.0, ymax=0.9), inherit.aes=FALSE, fill=blocks$Colour, colour=NA) +
    geom_text(data=blocks, aes(x=Center, y=-2.3, label=paste0(Ecotype, Symbol)), inherit.aes=FALSE, size=2.4, family="sans") +
    scale_fill_gradientn(colours=c("#FFF7EC", "#FDD49E", "#FC8D59", "#7F0000"), limits=c(0, 1), na.value="#D9D9D9", name="Co-occurrence\nindex") +
    scale_x_continuous(limits=c(0.5, nrow(C)+0.5), expand=c(0, 0)) +
    scale_y_reverse(limits=c(nrow(C)+0.5, -4.0), expand=c(0, 0)) +
    coord_equal(clip="off") +
    labs(title="EcoType recovery in pooled scRNA-seq atlases", subtitle="Four independent IBD atlases; * Cell 2021 nominal one-sided P < 0.05", x="Cell states", y="Cell states") +
    theme_classic(base_family="sans", base_size=7) +
    theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_text(size=7), legend.position="right", legend.title=element_text(size=6.5), legend.text=element_text(size=6), plot.title=element_text(size=9.5, face="bold"), plot.subtitle=element_text(size=6.5), plot.margin=margin(4, 4, 4, 4))

width_in <- 105/25.4
height_in <- 98/25.4
ragg::agg_png(file.path(out, "Figure4D_pooled_scRNA_ecotype_cooccurrence.png"), width=width_in, height=height_in, units="in", res=600, background="white")
print(figure_4d)
dev.off()
grDevices::pdf(file.path(out, "Figure4D_pooled_scRNA_ecotype_cooccurrence.pdf"), width=width_in, height=height_in, family="sans", useDingbats=FALSE)
print(figure_4d)
dev.off()
svglite::svglite(file.path(out, "Figure4D_pooled_scRNA_ecotype_cooccurrence.svg"), width=width_in, height=height_in, bg="white")
print(figure_4d)
dev.off()

# Step4. 保存图注口径和环境
writeLines(c(
    "One heat map is shown: all four prespecified scRNA-seq atlases were pooled before F1/F2 and C were calculated.",
    "Rows and columns are the 79 frozen EcoType-member states, ordered by E1-E8; black diagonal blocks mark frozen EcoTypes.",
    "The display truncates negative C values to the low end and restores diagonal self-correlation only for visualization.",
    "The statistical test uses the signed C matrix with diagonal NA; * denotes the Cell 2021 nominal one-sided normal-tail P < 0.05 without BH.",
    "Dataset-specific matrices and BH results remain sensitivity analyses and are not drawn as four main heat maps."
), file.path(out, "Figure4D_legend_and_method_boundary.txt"))
capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
