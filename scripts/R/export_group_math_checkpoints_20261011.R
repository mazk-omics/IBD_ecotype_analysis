# 导出Group抽样索引和系数供独立SciPy核对，同时比较完整残差与旧片段短向量残差。
# Rscript 本脚本 Group研究结果目录 原始mixture 固定fraction 新诊断目录。
# 仅导出24个预先均匀选定检查点，不执行缓存供应商函数，不修改研究表达。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
state <- readRDS(file.path(args[1], 'group_unfiltered_checkpoint.rds'))
mixture <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
fractions <- as.matrix(read.delim(args[3], row.names = 1, check.names = FALSE))
outdir <- args[4]
stopifnot(!dir.exists(outdir), identical(colnames(mixture), rownames(fractions)))
dir.create(outdir, recursive = TRUE)
coefficient_tables <- index_tables <- records <- vector('list', length(state$checkpoints))

# Step1. 每个基因保留完整30次抽样，Python可从同一输入直接再拟合
for (i in seq_along(state$checkpoints)) {
    checkpoint <- state$checkpoints[[i]]
    coefficient_tables[[i]] <- data.frame(GeneSymbol = checkpoint$gene, Bootstrap = seq_len(nrow(checkpoint$coefficients)), checkpoint$coefficients, check.names = FALSE)
    colnames(coefficient_tables[[i]])[-c(1, 2)] <- colnames(fractions)
    index_tables[[i]] <- data.frame(GeneSymbol = checkpoint$gene, Bootstrap = seq_len(nrow(checkpoint$samples)), checkpoint$samples, check.names = FALSE)
    gene <- checkpoint$gene
    prediction <- as.vector(fractions %*% state$beta[gene, ])
    short_y <- mixture[gene, checkpoint$samples[nrow(checkpoint$samples), ]]
    legacy_y <- rep(short_y, length.out = nrow(fractions))  # 明确展示旧片段发生的回收，而非科学残差
    records[[i]] <- data.frame(GeneSymbol = gene, full_rows = length(prediction), legacy_rows = length(short_y), full_residual_RSS = sum((mixture[gene, ] - prediction)^2), recycled_residual_RSS = sum((legacy_y - prediction)^2))
}
write.table(do.call(rbind, coefficient_tables), file.path(outdir, 'bootstrap_coefficients.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, index_tables), file.path(outdir, 'bootstrap_indices.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, records), file.path(outdir, 'legacy_residual_dimension_diagnostic.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
print(c(checked_genes = length(state$checkpoints), NNLS_checkpoints = sum(vapply(state$checkpoints, function(x) nrow(x$samples), integer(1)))))
