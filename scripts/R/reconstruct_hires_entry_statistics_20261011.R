# 独立重建旧版 HiRes 进入逐基因 worker 前的 bootstrap、解析 t 检验及全基因 BH。
# 在含 nnls 的 R 环境运行：Rscript 本脚本 原始mixture文件 固定fraction文件 新输出目录。
# 仅重建一个数值阶段；使用独立 R 随机流，不能当作完整官方 Group-mode/HiRes 软件。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
library(nnls)
mixture <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
fractions <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
outdir <- args[3]
stopifnot(!dir.exists(outdir), identical(colnames(mixture), rownames(fractions)), all(is.finite(mixture)), all(is.finite(fractions)), all(mixture >= 0), all(fractions >= 0))
dir.create(outdir, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))

# Step1. 固定整张表达表的统计范围，不能对子集24基因重新BH
gene_count <- nrow(mixture)
sample_count <- nrow(fractions)
type_count <- ncol(fractions)
stopifnot(sample_count > type_count, qr(fractions)$rank == type_count)
covariance_basis <- chol2inv(chol(crossprod(fractions)))
bootstrap_count <- 30L  # 固定旧版 HiRes 主调度，非 worker 的100次有放回 bootstrap
subsample_count <- as.integer(floor(sample_count * (2 / 3)))  # 旧 native doBootstrap(replace=FALSE)
stopifnot(subsample_count > type_count)
global_medians <- matrix(0, gene_count, type_count, dimnames = list(rownames(mixture), colnames(fractions)))
raw_pvalues <- matrix(NA_real_, gene_count, type_count, dimnames = dimnames(global_medians))
bootstrap_sds <- global_medians
empty_genes <- rowSums(mixture) == 0
solver_count <- 0L
started <- proc.time()[['elapsed']]
print(c(genes = gene_count, samples = sample_count, types = type_count, bootstrap = bootstrap_count, subsampled_rows = subsample_count, all_zero_genes = sum(empty_genes)))

# Step2. 各基因的系数中位数与全样本残差方差；空表达产生不可检验的p值
for (g in seq_len(gene_count)) {
    set.seed(20261011L + g)  # 独立研究约定；不宣称复制官方跨线程随机流
    if (!empty_genes[g]) {
        y <- mixture[g, ]
        coefficients <- matrix(NA_real_, bootstrap_count, type_count)
        for (b in seq_len(bootstrap_count)) {
            sampled_rows <- sample.int(sample_count, subsample_count, replace = FALSE)
            fitted <- nnls(fractions[sampled_rows, , drop = FALSE], y[sampled_rows])
            stopifnot(fitted$mode == 1)
            coefficients[b, ] <- coef(fitted)
        }
        beta <- apply(coefficients, 2, median)
        residuals <- y - as.vector(fractions %*% beta)
        residual_variance <- sum(residuals^2) / (sample_count - type_count)
        analytic_se <- sqrt(residual_variance * diag(covariance_basis))
        raw_pvalues[g, ] <- 2 * pt(beta / analytic_se, df = sample_count - type_count, lower.tail = FALSE)
        global_medians[g, ] <- beta
        bootstrap_sds[g, ] <- apply(coefficients, 2, sd)
        solver_count <- solver_count + bootstrap_count
    }
    if (g %% 2000L == 0L) print(c(completed_genes = g, elapsed_seconds = proc.time()[['elapsed']] - started))
}
entry_seconds <- proc.time()[['elapsed']] - started
write.table(global_medians, file.path(outdir, 'entry_bootstrap_medians.tsv'), sep = '\t', quote = FALSE, col.names = NA)
write.table(raw_pvalues, file.path(outdir, 'entry_raw_pvalues.tsv'), sep = '\t', quote = FALSE, col.names = NA, na = 'NA')
write.table(bootstrap_sds, file.path(outdir, 'entry_bootstrap_sds.tsv'), sep = '\t', quote = FALSE, col.names = NA)

# Step3. 每个类型在完整基因范围BH，再处理强marker排斥其他类型的规则
bh_qvalues <- raw_pvalues
for (j in seq_len(type_count)) {
    tested <- is.finite(raw_pvalues[, j])
    bh_qvalues[tested, j] <- p.adjust(raw_pvalues[tested, j], method = 'BH')
}
worker_qvalues <- bh_qvalues
strong_gene <- apply(worker_qvalues, 1, function(x) any(x < 1e-5, na.rm = TRUE))
worker_qvalues[strong_gene, ] <- ifelse(worker_qvalues[strong_gene, , drop = FALSE] > 0.01, 1, worker_qvalues[strong_gene, , drop = FALSE])
stopifnot(all(raw_pvalues[is.finite(raw_pvalues)] >= 0), all(bh_qvalues[is.finite(bh_qvalues)] <= 1))
write.table(bh_qvalues, file.path(outdir, 'entry_bh_qvalues.tsv'), sep = '\t', quote = FALSE, col.names = NA, na = 'NA')
write.table(worker_qvalues, file.path(outdir, 'entry_worker_qvalues.tsv'), sep = '\t', quote = FALSE, col.names = NA, na = 'NA')
summary <- data.frame(genes = gene_count, samples = sample_count, types = type_count, bootstrap_count = bootstrap_count, subsample_count = subsample_count, solver_count = solver_count, all_zero_genes = sum(empty_genes), strong_marker_genes = sum(strong_gene), entry_seconds = entry_seconds)
write.table(summary, file.path(outdir, 'entry_stage_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
saveRDS(list(global_medians = global_medians, raw_pvalues = raw_pvalues, bh_qvalues = bh_qvalues, worker_qvalues = worker_qvalues, samples = colnames(mixture), types = colnames(fractions), stage = 'Independent old HiRes entry statistics; not official Group output'), file.path(outdir, 'entry_statistics.rds'))
writeLines(c(paste('Original mixture:', args[1]), paste('Fixed fractions:', args[2]), 'Full gene-space BH; independent R sampling, no token or vendor program used.', 'This is the old HiRes entry statistic, not the official Group-mode Qvals output.', 'No worker refinement, S-mode, final NA/floor/auxiliary fallback or official RNG parity implemented here.'), file.path(outdir, 'METHOD_LIMITATIONS.txt'))
print(summary)
