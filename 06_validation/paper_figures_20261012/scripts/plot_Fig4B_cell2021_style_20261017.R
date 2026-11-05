# 使用冻结的53个状态、7个Ecotype和样本共现Jaccard，绘制Cell 2021 Figure 4B风格网络。
# 本脚本只改变可视化；不修改冻结表、不重跑EcoTyper，也不把边解释为通讯或空间接触。

# Step1. 定位项目和输出，载入绘图依赖
args = commandArgs(trailingOnly=TRUE)
project = if (length(args) >= 1) normalizePath(args[1], winslash='/', mustWork=TRUE) else normalizePath(getwd(), winslash='/', mustWork=TRUE)
paper = file.path(project,'06_validation/paper_figures_20261012')
source_dir = file.path(paper,'source_data')
figure_dir = file.path(paper,'figures')
report_dir = file.path(paper,'reports')
log_dir = file.path(paper,'logs')
qa_dir = if (length(args) >= 2) normalizePath(args[2], winslash='/', mustWork=TRUE) else file.path(paper,'scripts','nature_figure_qa')
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)

library(ggplot2)
library(patchwork)
library(svglite)
library(ragg)
library(jsonlite)
source(file.path(qa_dir,'panel_alignment.R'))
python_bin = Sys.which('python')
if (!nzchar(python_bin)) python_bin = Sys.which('python3')

prefix = 'Fig4B_ecotype_member_networks_Cell2021_style_20261017'
width_mm = 183
height_mm = 126
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Step2. 读入正式冻结成员、既有FR坐标和全部生态型内Jaccard边
official_members_path = file.path(project,'ctspec_region_ibd/Ecotypes/ecotypes.txt')
if (!file.exists(official_members_path)) official_members_path = file.path(project,'software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes/Ecotypes/discovery/ecotypes.txt')
stopifnot(file.exists(official_members_path))

official_members = read.delim(official_members_path, check.names=FALSE)
official_members$CellState = official_members$ID
frozen_states = read.delim(file.path(source_dir,'frozen_states.tsv'), check.names=FALSE)
edges = read.delim(file.path(source_dir,'Fig4B_all_within_ecotype_edges.tsv'), check.names=FALSE)
ecotypes = paste0('E',1:7)
expected_counts = c(E1=10,E2=9,E3=12,E4=3,E5=4,E6=9,E7=6)

stopifnot(nrow(official_members) == 53, nrow(frozen_states) == 53)
stopifnot(setequal(official_members$CellState, frozen_states$CellState))
stopifnot(all(table(official_members$Ecotype)[ecotypes] == expected_counts))
stopifnot(!any(official_members$Ecotype == 'E1' & official_members$CellType %in% c('Fibroblast','Goblet')))
stopifnot(all(edges$Ecotype %in% ecotypes), all(edges$Jaccard >= 0), all(edges$Jaccard <= 1))

parent_label = c(B='B', CD4_T='CD4', CD8_T='CD8', Endothelial='Endo', Enteroendocrine='EEC', Fibroblast='Fib', Glial='Glia', Goblet='Gob', ILC='ILC', Mast='Mast', Monocyte_Macrophage='Mac', Pericyte='Peri')
state_colours = c(S01='#D96AA5', S02='#9FC6D8', S03='#A8CF83', S04='#9B84BC', S05='#E9A15F', S06='#F2D15B', S07='#B79A86')
stopifnot(all(frozen_states$CellType %in% names(parent_label)))
stopifnot(all(frozen_states$State %in% names(state_colours)))

# Step3. 保留论文版本FR坐标；每个Ecotype仅平移缩放以统一面板占用范围
node_tables = list()
for (ecotype in ecotypes) {
    nodes = read.delim(file.path(source_dir,paste0(ecotype,'_Fig4B_FR_nodes.tsv')), check.names=FALSE)
    members = frozen_states[frozen_states$Ecotype == ecotype,]
    stopifnot(nrow(nodes) == expected_counts[ecotype], setequal(nodes$CellState,members$CellState))
    nodes = merge(nodes,members[,c('CellState','CellType','State')],by='CellState',all.x=TRUE,sort=FALSE,suffixes=c('.coord','.frozen'))
    if ('CellType.frozen' %in% names(nodes)) nodes$CellType = nodes$CellType.frozen
    if ('State.frozen' %in% names(nodes)) nodes$State = nodes$State.frozen
    nodes$X = nodes$X - mean(nodes$X)
    nodes$Y = nodes$Y - mean(nodes$Y)
    span = max(diff(range(nodes$X)),diff(range(nodes$Y)))
    nodes$PlotX = nodes$X / span * 2
    nodes$PlotY = nodes$Y / span * 2
    nodes$Label = unname(parent_label[nodes$CellType])
    node_tables[[ecotype]] = nodes
}

# Step4. 为每个Ecotype构造真实边坐标；仅不画Jaccard=0，线宽为固定线性映射
plots = list()
edge_audit = data.frame()
for (ecotype in ecotypes) {
    nodes = node_tables[[ecotype]]
    own_edges = edges[edges$Ecotype == ecotype,]
    positive_edges = own_edges[own_edges$Jaccard > 0,]
    positive_edges = merge(positive_edges,nodes[,c('CellState','PlotX','PlotY')],by.x='From',by.y='CellState',all.x=TRUE,sort=FALSE)
    names(positive_edges)[names(positive_edges) %in% c('PlotX','PlotY')] = c('X1','Y1')
    positive_edges = merge(positive_edges,nodes[,c('CellState','PlotX','PlotY')],by.x='To',by.y='CellState',all.x=TRUE,sort=FALSE)
    names(positive_edges)[names(positive_edges) %in% c('PlotX','PlotY')] = c('X2','Y2')
    stopifnot(!anyNA(positive_edges[,c('X1','Y1','X2','Y2')]))
    positive_edges$LineWidth = 0.08 + 0.48 * positive_edges$Jaccard

    # 线段在所有节点圆盘内的部分不绘制，避免视觉上让边穿过节点文字
    edge_segments = data.frame()
    node_radius = 0.17
    for (edge_index in seq_len(nrow(positive_edges))) {
        edge = positive_edges[edge_index,]
        vx = edge$X2 - edge$X1
        vy = edge$Y2 - edge$Y1
        a = vx^2 + vy^2
        blocked = data.frame(Start=numeric(),End=numeric())
        for (node_index in seq_len(nrow(nodes))) {
            fx = edge$X1 - nodes$PlotX[node_index]
            fy = edge$Y1 - nodes$PlotY[node_index]
            b = 2 * (fx * vx + fy * vy)
            c = fx^2 + fy^2 - node_radius^2
            discriminant = b^2 - 4 * a * c
            if (discriminant >= 0) {
                roots = sort(c((-b - sqrt(discriminant)) / (2 * a),(-b + sqrt(discriminant)) / (2 * a)))
                interval = c(max(0,roots[1]),min(1,roots[2]))
                if (interval[2] > interval[1]) blocked = rbind(blocked,data.frame(Start=interval[1],End=interval[2]))
            }
        }
        blocked = blocked[order(blocked$Start),]
        merged = data.frame(Start=numeric(),End=numeric())
        for (interval_index in seq_len(nrow(blocked))) {
            if (!nrow(merged) || blocked$Start[interval_index] > merged$End[nrow(merged)]) {
                merged = rbind(merged,blocked[interval_index,])
            } else {
                merged$End[nrow(merged)] = max(merged$End[nrow(merged)],blocked$End[interval_index])
            }
        }
        visible = data.frame(Start=numeric(),End=numeric())
        cursor = 0
        for (interval_index in seq_len(nrow(merged))) {
            if (merged$Start[interval_index] > cursor) visible = rbind(visible,data.frame(Start=cursor,End=merged$Start[interval_index]))
            cursor = max(cursor,merged$End[interval_index])
        }
        if (cursor < 1) visible = rbind(visible,data.frame(Start=cursor,End=1))
        if (nrow(visible)) {
            visible$X1 = edge$X1 + visible$Start * vx
            visible$Y1 = edge$Y1 + visible$Start * vy
            visible$X2 = edge$X1 + visible$End * vx
            visible$Y2 = edge$Y1 + visible$End * vy
            visible$LineWidth = edge$LineWidth
            edge_segments = rbind(edge_segments,visible[,c('X1','Y1','X2','Y2','LineWidth')])
        }
    }

    edge_audit = rbind(edge_audit,data.frame(Ecotype=ecotype,States=nrow(nodes),AllPairs=nrow(own_edges),PositiveEdges=nrow(positive_edges),ZeroEdges=sum(own_edges$Jaccard == 0),MinPositiveJaccard=min(positive_edges$Jaccard),MaxJaccard=max(positive_edges$Jaccard)))

    plots[[ecotype]] = ggplot() +
        geom_segment(data=edge_segments,aes(x=X1,y=Y1,xend=X2,yend=Y2,linewidth=LineWidth),colour='black',lineend='round',alpha=0.92) +
        scale_linewidth_identity() +
        geom_point(data=nodes,aes(x=PlotX,y=PlotY,fill=State),shape=21,size=7.2,stroke=0.01,colour='#00000000') +
        geom_text(data=nodes,aes(x=PlotX,y=PlotY,label=Label),family='sans',size=2.25,colour='black',lineheight=0.92) +
        scale_fill_manual(values=state_colours,breaks=names(state_colours),drop=FALSE,name='Cell state',labels=names(state_colours)) +
        coord_equal(xlim=c(-1.25,1.25),ylim=c(-1.22,1.22),expand=FALSE,clip='off') +
        labs(title=ecotype) +
        theme_void(base_family='sans',base_size=8) +
        theme(plot.title=element_text(size=10.2,face='plain',hjust=0.12,margin=margin(b=1.2)),plot.margin=margin(1.5,1.5,1.5,1.5),legend.position='none')
}

# Step5. 按原图上4下3排版并统一S01-S07图例，不显示父类型颜色图例
legend_data = data.frame(State=factor(names(state_colours),levels=names(state_colours)),X=seq(1.35,7.95,length.out=7),Y=0.5)
legend_plot = ggplot(legend_data,aes(x=X,y=Y)) +
    annotate('text',x=0.18,y=0.5,label='Cell state',family='sans',size=2.75,hjust=0) +
    geom_point(aes(fill=State),shape=22,size=4.5,stroke=0.28,colour='#666666') +
    geom_text(aes(x=X+0.22,label=State),family='sans',size=2.5,hjust=0) +
    scale_fill_manual(values=state_colours,drop=FALSE,guide='none') +
    coord_cartesian(xlim=c(0,8.75),ylim=c(0,1),clip='off',expand=FALSE) +
    theme_void(base_family='sans') +
    theme(plot.margin=margin(0,0,0,0))

design = '\nABCD\nEFG#\nHHHH\n'
figure = wrap_plots(c(plots[ecotypes],list(legend_plot)),design=design,heights=c(1,1,0.16)) +
    plot_annotation(title=expression(bold(B)~~'Co-occurrence of cell states in IBD ecotypes'),theme=theme(plot.title=element_text(family='sans',size=13,face='plain',hjust=0.015,margin=margin(b=2))))

# Step6. 在最终物理尺寸运行面板对齐门禁，再导出PNG、PDF和SVG
require_patchwork_panel_alignment(
    figure,
    manifest_path=file.path(report_dir,paste0(prefix,'.alignment-layout.json')),
    report_path=file.path(report_dir,paste0(prefix,'.alignment.json')),
    overlay_svg=file.path(report_dir,paste0(prefix,'.alignment.svg')),
    width_in=width_in,
    height_in=height_in,
    panel_ids=c(ecotypes,'Legend'),
    row_groups=list(c('E1','E2','E3','E4'),c('E5','E6','E7')),
    column_groups=list(c('E1','E5'),c('E2','E6'),c('E3','E7')),
    audit_script=file.path(qa_dir,'audit_panel_alignment.py'),
    python=python_bin,
    tolerance_pt=1.5,
    gutter_tolerance_pt=1.5,
    strict=TRUE
)

png_path = file.path(figure_dir,paste0(prefix,'.png'))
pdf_path = file.path(figure_dir,paste0(prefix,'.pdf'))
svg_path = file.path(figure_dir,paste0(prefix,'.svg'))
ragg::agg_png(png_path,width=width_in,height=height_in,units='in',res=600,background='white')
print(figure)
dev.off()
grDevices::pdf(pdf_path,width=width_in,height=height_in,family='sans',onefile=TRUE,useDingbats=FALSE)
print(figure)
dev.off()
svglite::svglite(svg_path,width=width_in,height=height_in,bg='white')
print(figure)
dev.off()

# Step7. 保存成员/边核查和简明复现日志
write.table(edge_audit,file.path(report_dir,paste0(prefix,'.data-audit.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
member_audit = data.frame(Ecotype=ecotypes,States=as.integer(table(frozen_states$Ecotype)[ecotypes]),ParentTypes=vapply(ecotypes,function(x) paste(frozen_states$CellType[frozen_states$Ecotype == x],collapse=';'),character(1)))
write.table(member_audit,file.path(report_dir,paste0(prefix,'.member-audit.tsv')),sep='\t',row.names=FALSE,quote=FALSE)

log_lines = c(
    paste('Completed:',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')),
    paste('R:',R.version.string),
    paste('ggplot2:',as.character(packageVersion('ggplot2'))),
    paste('patchwork:',as.character(packageVersion('patchwork'))),
    paste('Frozen members:',official_members_path),
    paste('Frozen states:',file.path(source_dir,'frozen_states.tsv')),
    paste('Frozen edges:',file.path(source_dir,'Fig4B_all_within_ecotype_edges.tsv')),
    paste('Members:',paste(paste0(ecotypes,'=',expected_counts),collapse=', ')),
    paste('Positive edges:',sum(edge_audit$PositiveEdges),'of',sum(edge_audit$AllPairs),'within-ecotype pairs'),
    'Edge rule: draw Jaccard > 0; linewidth_mm = 0.08 + 0.48 * Jaccard; no inferred edges added.',
    'Node rule: label = parent cell-type abbreviation; fill = state number S01-S07 with one global legend.',
    'Interpretation boundary: edges are sample co-occurrence Jaccard, not communication or spatial contact.',
    paste('PNG:',png_path),
    paste('PDF:',pdf_path),
    paste('SVG:',svg_path)
)
writeLines(log_lines,file.path(log_dir,paste0(prefix,'.log')))
cat(paste(log_lines,collapse='\n'),'\n')
