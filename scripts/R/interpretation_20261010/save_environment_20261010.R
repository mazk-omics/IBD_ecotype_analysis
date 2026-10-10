# 保存解释环境真实R包信息，避免源码安装与conda登记信息不同造成误判。
# 在新解释环境运行；只写环境记录，不读取生物学数据。
library(Seurat)
library(fgsea)
library(data.table)
library(igraph)
library(metap)
library(jsonlite)
out <- '/home/mazekai/IBD_EcoTyper/06_validation/ecotyper_original_method_annotation_20261010/environment'
p <- installed.packages()
write.table(p[,c('Package','Version','Built','LibPath')], file.path(out,'R_packages.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
stopifnot(as.character(packageVersion('Seurat'))=='3.1.3', as.character(packageVersion('fgsea'))=='1.14.0')
capture.output(sessionInfo(), file=file.path(out,'R_sessionInfo.txt'))
print('PAPER_PACKAGE_LOAD_OK')
