# 独立研究实现Group表达估计和自适应CV过滤，读取线性mixture和已固定的fraction。
# 含nnls的R环境运行：Rscript 本脚本 原始mixture 固定fraction 新输出目录。
# 采用完整表达计算残差；旧缓存R片段的短向量残差存在维数疑点，本实现不宣称版本等价。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
library(nnls)
mixture <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
fractions <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
outdir <- args[3]
stopifnot(!dir.exists(outdir), identical(colnames(mixture), rownames(fractions)), !anyDuplicated(rownames(mixture)), all(is.finite(mixture)), all(mixture >= 0), all(is.finite(fractions)), all(fractions >= 0))
dir.create(outdir, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))

# Step1. 完整基因范围决定BH和CV阈值，不能只在1249个目标基因中计算过滤
gene_count <- nrow(mixture)
sample_count <- nrow(fractions)
type_count <- ncol(fractions)
stopifnot(sample_count > type_count, qr(fractions)$rank == type_count)
bootstrap_count <- 30L  # 缓存Group函数默认，与HiRes worker的100次抽样不同
subsample_count <- as.integer(max(5, type_count, sample_count * 0.75))  # Group无放回75%样本
covariance_basis <- chol2inv(chol(crossprod(fractions)))
beta <- matrix(0, gene_count, type_count, dimnames = list(rownames(mixture), colnames(fractions)))
geometric_cv <- bootstrap_sd <- beta
raw_pvalues <- matrix(NA_real_, gene_count, type_count, dimnames = dimnames(beta))
analytic_se <- raw_pvalues
checkpoint_rows <- unique(as.integer(round(seq(1, gene_count, length.out = 24))))
checkpoints <- vector('list', length(checkpoint_rows))
started <- proc.time()[['elapsed']]
print(c(genes = gene_count, samples = sample_count, types = type_count, bootstrap = bootstrap_count, subsample_rows = subsample_count))

# Step2. 非负回归系数在log空间取中位数；StdErr输出是bootstrap SD，解析SE另存
for (g in seq_len(gene_count)) {
    set.seed(20261015L + g, kind = "L'Ecuyer-CMRG")  # 独立逐基因随机流，不复制官方线程/种子
    y <- mixture[g, ]
    coefficients <- matrix(0, bootstrap_count, type_count)
    sampled_indices <- matrix(NA_integer_, bootstrap_count, subsample_count)
    for (b in seq_len(bootstrap_count)) {
        sampled_rows <- sample.int(sample_count, subsample_count, replace = FALSE)
        sampled_indices[b, ] <- sampled_rows
        fitted <- nnls(fractions[sampled_rows, , drop = FALSE], y[sampled_rows])
        stopifnot(fitted$mode == 1)
        coefficients[b, ] <- coef(fitted)
    }
    beta[g, ] <- 2^apply(log2(coefficients + 1), 2, median) - 1
    geometric_cv[g, ] <- sqrt(exp(apply(log1p(coefficients), 2, sd)^2) - 1)
    bootstrap_sd[g, ] <- apply(coefficients, 2, sd)
    residuals <- y - as.vector(fractions %*% beta[g, ])  # 完整120表达与120行设计逐样本匹配
    residual_variance <- sum(residuals^2) / (sample_count - type_count)
    analytic_se[g, ] <- sqrt(residual_variance * diag(covariance_basis))
    raw_pvalues[g, ] <- 2 * pt(beta[g, ] / analytic_se[g, ], df = sample_count - type_count, lower.tail = FALSE)
    if (g %in% checkpoint_rows) checkpoints[[match(g, checkpoint_rows)]] <- list(gene = rownames(mixture)[g], row = g, samples = sampled_indices, coefficients = coefficients)
    if (g %% 2000L == 0L) print(c(completed_genes = g, elapsed_seconds = proc.time()[['elapsed']] - started))
}
saveRDS(list(beta = beta, geometric_cv = geometric_cv, bootstrap_sd = bootstrap_sd, analytic_se = analytic_se, raw_pvalues = raw_pvalues, checkpoints = checkpoints), file.path(outdir, 'group_unfiltered_checkpoint.rds'))

# Step3. 强marker排斥规则使用q>0.25，零表达保留为0；HiRes入口的0.01规则不能借用
qvalues <- raw_pvalues
for (j in seq_len(type_count)) qvalues[, j] <- p.adjust(raw_pvalues[, j], method = 'BH')
qvalues[is.na(qvalues)] <- 1
expression <- beta
strong_gene <- apply(qvalues, 1, min) < 1e-5
expression[strong_gene, ] <- ifelse(qvalues[strong_gene, , drop = FALSE] > 0.25, 0, expression[strong_gene, , drop = FALSE])
# 论文写置0，缓存旧R写置1；当前Job24未过滤表没有1。保持可解释的论文分支，不暗中拟合版本。

# Step4. 使用缓存实现的CV曲线几何规则，边界CV等于阈值也置NA
filtered <- expression
thresholds <- vector('list', type_count)
for (j in seq_len(type_count)) {
    sorted_cv <- sort(geometric_cv[, j])
    ycut <- max(1, quantile(geometric_cv[geometric_cv[, j] < 999999, j], 0.75))
    endpoint <- which.min(abs(sorted_cv - ycut))
    normalized_rank <- seq_len(endpoint) / endpoint
    squared_distance <- (normalized_rank - ycut)^2 + sorted_cv[seq_len(endpoint)]^2
    threshold <- sorted_cv[which.min(squared_distance)]
    excluded <- geometric_cv[, j] >= threshold
    filtered[excluded, j] <- NA_real_
    thresholds[[j]] <- data.frame(CellType = colnames(fractions)[j], ycut = ycut, endpoint_rank = endpoint, cv_threshold = threshold, excluded_pairs = sum(excluded))
}
write.table(do.call(rbind, thresholds), file.path(outdir, 'group_cv_thresholds.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

# Step5. 保存可审阅中间量和外层回填入口；这些研究表不用于正式discovery/recovery
tables <- list(GEPs = expression, GEPs_Filtered = filtered, GEPs_CVs = geometric_cv, GEPs_Pvals = raw_pvalues, GEPs_Qvals = qvalues, GEPs_StdErrs = bootstrap_sd, analytic_SE = analytic_se, bootstrap_beta_before_q_guard = beta)
for (name in names(tables)) write.table(tables[[name]], file.path(outdir, paste0(name, '.tsv')), sep = '\t', quote = FALSE, col.names = NA, na = 'NA')
summary <- data.frame(genes = gene_count, samples = sample_count, types = type_count, bootstrap_count = bootstrap_count, subsample_count = subsample_count, solver_count = gene_count * bootstrap_count, strong_marker_genes = sum(strong_gene), filtered_NA_pairs = sum(is.na(filtered)), seconds = proc.time()[['elapsed']] - started)
write.table(summary, file.path(outdir, 'group_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(c('Independent research Group calculation; upstream fractions remain fixed official input.', 'Original linear expression, no S-mode or signature estimation implemented here.', 'Full expression residual replaces the dimensionally inconsistent short-vector residual in cached legacy R.', 'Paper zero guard retained; cached legacy R uses one. No exact web-version or RNG equivalence claimed.', 'CV geometry matches cached R code, conditional official-CV validation must be reported separately.', 'Bootstrap SD is not the analytical regression SE; both tables are saved.'), file.path(outdir, 'METHOD_LIMITATIONS.txt'))
print(summary)
