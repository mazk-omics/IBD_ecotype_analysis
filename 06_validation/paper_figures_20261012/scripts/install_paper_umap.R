# 安装论文使用的umap 0.2.0.0到本批次独立R库，不覆盖现有环境。
# 在ecotyper_interpretation环境运行；仅补绘图分析包，原始数据与EcoTyper代码保持只读。
# Step1. 明确版本与安装位置
lib = '/home/mazekai/IBD_EcoTyper/06_validation/paper_figures_20261012/R_library'
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
install.packages('https://cran.r-project.org/src/contrib/Archive/umap/umap_0.2.0.0.tar.gz', repos=NULL, type='source', lib=lib)
stopifnot(as.character(packageVersion('umap')) == '0.2.0.0')
print(packageVersion('umap'))
writeLines(capture.output(sessionInfo()), '/home/mazekai/IBD_EcoTyper/06_validation/paper_figures_20261012/logs/R_sessionInfo.txt')
