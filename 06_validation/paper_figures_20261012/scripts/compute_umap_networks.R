# 按Cell 2021 Fig2B/4B方法计算坐标；本脚本不画图、不重拟合状态。
# 在ecotyper_interpretation运行；读取冻结H*欧氏距离和Jaccard，保存UMAP/FR坐标和配置。
# Step1. 载入论文版本umap和已有igraph，记录版本边界
out = '/home/mazekai/IBD_EcoTyper/06_validation/paper_figures_20261012'
.libPaths(c(file.path(out,'R_library'), .libPaths()))
library(umap)
library(igraph)
stopifnot(as.character(packageVersion('umap')) == '0.2.0.0')
states = read.delim(file.path(out,'source_data/frozen_states.tsv'), check.names=FALSE)
set.seed(20261012) # 本项目固定复现种子；论文没有报告种子

# Step2. 每父类型仅使用通过冻结QC分配的活检，完整H*不做重新归一化
for (ct in unique(states$CellType)) {
    Hstar = read.delim(file.path(out,'source_data',paste0(ct,'_Fig2B_Hstar_assigned.tsv')), row.names=1, check.names=FALSE)
    annotation = read.delim(file.path(out,'source_data',paste0(ct,'_Fig2A_sample_order.tsv')), row.names=1, check.names=FALSE)
    stopifnot(identical(rownames(Hstar),rownames(annotation)))
    distance = as.matrix(dist(Hstar, method='euclidean'))
    embedding = umap(distance, input='dist', random_state=20261012)
    coordinates = data.frame(ID=rownames(Hstar), UMAP1=embedding$layout[,1], UMAP2=embedding$layout[,2], State=annotation$State)
    write.table(coordinates, file.path(out,'source_data',paste0(ct,'_Fig2B_UMAP.tsv')), sep='\t', row.names=FALSE, quote=FALSE)
    saveRDS(embedding$config, file.path(out,'source_data',paste0(ct,'_UMAP_config.rds')))
    cat('UMAP complete:',ct,nrow(coordinates),'\n')
}

# Step3. 每生态型全部冻结节点保留，边权使用既有Jaccard；0边不画，FR不表示通讯或因果
J = as.matrix(read.delim(file.path(out,'source_data/Fig4B_frozen_jaccard.tsv'), row.names=1, check.names=FALSE))
for (ecotype in unique(states$Ecotype)) {
    member = states[states$Ecotype==ecotype,]
    weights = J[member$CellState,member$CellState]
    diag(weights) = 0
    graph = graph_from_adjacency_matrix(weights, mode='undirected', weighted=TRUE, diag=FALSE)
    set.seed(20261012)
    position = layout_with_fr(graph, weights=E(graph)$weight)
    nodes = data.frame(Ecotype=ecotype, CellState=V(graph)$name, X=position[,1], Y=position[,2])
    nodes = merge(nodes, member[,c('CellState','CellType','State')],by='CellState',sort=FALSE)
    write.table(nodes,file.path(out,'source_data',paste0(ecotype,'_Fig4B_FR_nodes.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
    cat('FR complete:',ecotype,vcount(graph),ecount(graph),'\n')
}
writeLines(capture.output(sessionInfo()),file.path(out,'logs/R_numerical_sessionInfo.txt'))
writeLines(capture.output(umap.defaults),file.path(out,'reports/umap_0.2.0.0_defaults.txt'))
