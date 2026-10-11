# 诊断旧Group R片段中的短表达向量回收，会怎样改变p/q值；不是正式估计方法。
# Rscript 本脚本 Group研究目录 原始mixture 固定fraction 新诊断目录。
# 复用已保存系数，只重建最后一次抽样索引，不再拟合，也不替换主流程结果。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
state <- readRDS(file.path(args[1], 'group_unfiltered_checkpoint.rds'))
mixture <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
fractions <- as.matrix(read.delim(args[3], row.names = 1, check.names = FALSE))
outdir <- args[4]
stopifnot(!dir.exists(outdir), identical(colnames(mixture), rownames(fractions)))
dir.create(outdir, recursive = TRUE)

# Step1. 显式回收90个随机顺序值到120行，不让R隐式警告淹没诊断
legacy_pvalues <- state$raw_pvalues
covariance_basis <- chol2inv(chol(crossprod(fractions)))
for (g in seq_len(nrow(mixture))) {
    set.seed(20261015L + g, kind = "L'Ecuyer-CMRG")
    for (b in seq_len(30L)) chosen <- sample.int(nrow(fractions), as.integer(max(5, ncol(fractions), nrow(fractions) * .75)), replace = FALSE)
    y_recycled <- rep(mixture[g, chosen], length.out = nrow(fractions))
    prediction <- as.vector(fractions %*% state$beta[g, ])
    variance <- sum((y_recycled - prediction)^2) / (nrow(fractions) - ncol(fractions))
    se <- sqrt(variance * diag(covariance_basis))
    legacy_pvalues[g, ] <- 2 * pt(state$beta[g, ] / se, df = nrow(fractions) - ncol(fractions), lower.tail = FALSE)
}
legacy_qvalues <- legacy_pvalues
for (j in seq_len(ncol(fractions))) legacy_qvalues[, j] <- p.adjust(legacy_pvalues[, j], method = 'BH')
legacy_qvalues[is.na(legacy_qvalues)] <- 1
write.table(legacy_pvalues, file.path(outdir, 'recycled_residual_Pvals.tsv'), sep = '\t', quote = FALSE, col.names = NA)
write.table(legacy_qvalues, file.path(outdir, 'recycled_residual_Qvals.tsv'), sep = '\t', quote = FALSE, col.names = NA)
print(c(recycled_residual_strong_genes = sum(apply(legacy_qvalues, 1, min) < 1e-5), full_residual_strong_genes = sum(apply(read.delim(file.path(args[1], 'GEPs_Qvals.tsv'), row.names = 1), 1, min) < 1e-5)))
writeLines('Diagnostic only: reproducing vector recycling does not establish current web implementation or make this valid inference. Main Group outputs remain unchanged.', file.path(outdir, 'LIMITATION.txt'))
