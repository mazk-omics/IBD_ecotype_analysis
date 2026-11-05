# 保存真正使用的UMAP配置，避免自定义print绕过capture.output导致空白文本。
out = '/home/mazekai/IBD_EcoTyper/06_validation/paper_figures_20261012'
.libPaths(c(file.path(out,'R_library'),.libPaths()))
library(umap)
config = readRDS(file.path(out,'source_data/B_UMAP_config.rds'))
dput(config,file=file.path(out,'reports/umap_0.2.0.0_used_config.txt'))
