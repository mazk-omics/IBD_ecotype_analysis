# 对外部真实count-based混合物做BayesPrism验证；所有数据集使用同一冻结Joint20参考。
# 服务器bayesprism环境运行；smoke为3例，formal为100例；保留完整20类输出和initial/final theta。
suppressPackageStartupMessages(library(BayesPrism))
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
parent <- args[1]
dataset <- args[2]
mode <- args[3]
folder <- file.path(parent, '08_external_pseudobulk', dataset)
input <- file.path(folder, if (mode == 'smoke') '01_smoke' else '02_simulation')
output <- if (mode == 'smoke') input else file.path(folder, '03_bayesprism')
started <- Sys.time()

# Step1. 同一参考及输入检查，不删类型，不用优化后的signature作为reference
reference <- readRDS(file.path(parent, '01_reference/clean_reference.rds'))
types <- reference$types
bulk <- read.delim(gzfile(file.path(input, 'pseudo_bulk_raw_counts.tsv.gz')), row.names = 1, check.names = FALSE)
mixture <- t(as.matrix(bulk))
stopifnot(length(types) == 20, setequal(unique(reference$labels), types), nrow(mixture) == if (mode == 'smoke') 3 else 100, all(is.finite(mixture)), all(mixture >= 0), all(abs(mixture - round(mixture)) < 1e-8), all(rowSums(mixture) > 0))
prism <- new.prism(reference = reference$reference, input.type = 'count.matrix', cell.type.labels = reference$labels, cell.state.labels = reference$labels, key = NULL, mixture = mixture, outlier.cut = 0.01, outlier.fraction = 0.1, pseudo.min = 1e-8)
stopifnot(setequal(rownames(prism@phi_cellType@phi), types))
writeLines(colnames(prism@mixture), file.path(output, 'retained_genes.txt'))
write.table(data.frame(gene = colnames(mixture), in_frozen_clean_reference = colnames(mixture) %in% colnames(reference$reference), retained_in_prism = colnames(mixture) %in% colnames(prism@mixture)), file.path(output, 'gene_filter_audit.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
saveRDS(prism, file.path(output, 'prism_input.rds'))
cat(dataset, mode, 'prism:', dim(prism@mixture), '\n')

# Step2. 默认迭代，smoke通过后才允许生成正式模拟
ncores <- if (mode == 'smoke') 3 else 12
set.seed(20261013)
bp <- run.prism(prism = prism, n.cores = ncores, update.gibbs = TRUE)
initial <- get.fraction(bp, which.theta = 'first', state.or.type = 'type')[, types, drop = FALSE]
final <- get.fraction(bp, which.theta = 'final', state.or.type = 'type')[, types, drop = FALSE]
stopifnot(identical(rownames(final), rownames(mixture)), ncol(final) == 20, all(is.finite(final)), all(final >= 0), max(abs(rowSums(final) - 1)) < 1e-6)
write.table(initial, file.path(output, 'theta_initial.tsv'), sep = '\t', col.names = NA, quote = FALSE)
write.table(final, file.path(output, 'theta_final.tsv'), sep = '\t', col.names = NA, quote = FALSE)
write.table(bp@posterior.initial.cellState@theta.cv[, types, drop = FALSE], file.path(output, 'theta_initial_cv.tsv'), sep = '\t', col.names = NA, quote = FALSE)
write.table(bp@posterior.theta_f@theta.cv[, types, drop = FALSE], file.path(output, 'theta_final_cv.tsv'), sep = '\t', col.names = NA, quote = FALSE)
saveRDS(bp, file.path(output, 'BayesPrism_result.rds'))
writeLines(capture.output(bp@control_param), file.path(output, 'installed_run_control.txt'))
writeLines(capture.output(sessionInfo()), file.path(output, 'R_sessionInfo.txt'))
writeLines(c(paste('elapsed_seconds', as.numeric(difftime(Sys.time(), started, units = 'secs'))), paste('n.cores', ncores), 'update.gibbs TRUE', 'seed 20261013', 'input.type count.matrix', 'same frozen clean_reference.rds', 'outlier.cut 0.01', 'outlier.fraction 0.1', 'pseudo.min 1e-8'), file.path(output, 'run_parameters.txt'))
if (mode == 'smoke') writeLines('COMPLETE: raw counts; full20 valid final theta; author-derived known truth', file.path(input, 'smoke_complete.txt'))
cat(dataset, mode, 'BAYESPRISM_COMPLETE\n')
