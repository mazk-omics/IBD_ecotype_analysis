# 独立实现 NNLS、全局重抽样和 HiRes 滑窗核心；固定 Job24 比例，显式选择原始或调整后表达输入。
# 在安装 nnls 的 R 环境运行：Rscript 本脚本 输入目录 新输出目录；可附基因起止行及第五参数 mixture 路径。
# 这不是完整 HiRes：尚不含显著性筛选、两段收缩、OLS 校准、S-mode 或官方随机流，禁止用于正式 EcoTyper。
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
outdir <- args[2]
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
library(nnls)
stopifnot(length(args) >= 2, !dir.exists(outdir))  # 保留每次独立结果，避免覆盖原件或先前检查点
dir.create(outdir, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))

# Step1. 确认基因、样本和类型顺序；全部样本固定，仅拆目标基因
mixture_path <- if (length(args) >= 5) args[5] else file.path(input_dir, 'mixture_adjusted_24genes.tsv')
mixture <- as.matrix(read.delim(mixture_path, row.names = 1, check.names = FALSE))
fractions <- as.matrix(read.delim(file.path(input_dir, 'fractions_fixed.tsv'), row.names = 1, check.names = FALSE))
official_group <- as.matrix(read.delim(file.path(input_dir, 'official_group_24genes.tsv'), row.names = 1, check.names = FALSE))
selection <- read.delim(file.path(input_dir, 'gene_selection.tsv'), check.names = FALSE)
stopifnot(identical(colnames(mixture), rownames(fractions)), identical(colnames(fractions), colnames(official_group)), identical(rownames(mixture), rownames(official_group)))
stopifnot(all(is.finite(mixture)), all(is.finite(fractions)), all(fractions >= 0), max(abs(rowSums(fractions) - 1)) < 1e-8)
gene_start <- if (length(args) >= 3) as.integer(args[3]) else 1L
gene_end <- if (length(args) >= 4) as.integer(args[4]) else nrow(mixture)
stopifnot(gene_start >= 1, gene_end <= nrow(mixture), gene_start <= gene_end)
gene_rows <- seq.int(gene_start, gene_end)
sample_count <- nrow(fractions)
type_count <- ncol(fractions)
window_size <- ceiling(max(type_count + 1, min(sample_count / 2, 4 * type_count)))  # 旧数值版本的默认公式
window_count <- sample_count - window_size + 1
print(c(genes = length(gene_rows), samples = sample_count, types = type_count, window = window_size))
fraction_diagnostics <- data.frame(samples = sample_count, types = type_count, rank = qr(fractions)$rank, kappa = kappa(fractions, exact = TRUE))
write.table(fraction_diagnostics, file.path(outdir, 'fraction_diagnostics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

# Step2. 用已知生成值及边界解检查求解器，而不是把实现自己当作答案
known_expression <- seq_len(type_count) * 3
known_mixture <- as.vector(fractions %*% known_expression)
known_fit <- nnls(fractions, known_mixture)
boundary_fit <- nnls(diag(3), c(2, -1, 4))
known_coefficient_error <- max(abs(coef(known_fit) - known_expression))
known_residual_error <- max(abs(known_fit$residuals))
boundary_error <- max(abs(coef(boundary_fit) - c(2, 0, 4)))
checks <- data.frame(check = c('known_coefficients', 'known_mixture_reconstruction', 'nonnegative_boundary'), error = c(known_coefficient_error, known_residual_error, boundary_error))
checks$passed <- checks$error < 1e-8
write.table(checks, file.path(outdir, 'solver_checks.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
stopifnot(all(checks$passed))
print(checks)

# Step3. 每基因保存检查点；按原始基因行固定独立种子，便于以后分工，不冒充官方 MT 随机流
global_nnls <- matrix(NA_real_, length(gene_rows), type_count, dimnames = list(rownames(mixture)[gene_rows], colnames(fractions)))
global_median <- global_nnls
core_smoothed <- array(NA_real_, c(length(gene_rows), sample_count, type_count), dimnames = list(rownames(mixture)[gene_rows], colnames(mixture), colnames(fractions)))
timings <- vector('list', length(gene_rows))
kkt_records <- vector('list', length(gene_rows))
bootstrap_global <- 100L  # 论文全局 bootstrap 默认；这里只计算其系数中位数，未执行统计过滤
bootstrap_window <- 5L  # 旧数值版本的每窗重抽样次数
for (output_row in seq_along(gene_rows)) {
    gene_row <- gene_rows[output_row]
    gene <- rownames(mixture)[gene_row]
    set.seed(20261011L + selection$full_mixture_row_1based[gene_row])
    y <- mixture[gene_row, ]
    started <- proc.time()[['elapsed']]
    fit <- nnls(fractions, y)
    stopifnot(fit$mode == 1)
    global_nnls[output_row, ] <- coef(fit)
    gradient <- as.vector(crossprod(fractions, as.vector(fractions %*% coef(fit)) - y))
    gradient_scale <- max(1, max(abs(crossprod(fractions, y))))
    kkt_records[[output_row]] <- data.frame(GeneSymbol = gene, nonnegative_violation = max(0, -min(coef(fit))), dual_violation_scaled = max(0, -min(gradient)) / gradient_scale, complementarity_scaled = max(abs(coef(fit) * gradient)) / (gradient_scale * max(1, max(coef(fit)))))
    global_seconds <- proc.time()[['elapsed']] - started

    # 全局有放回抽样；不把系数归一化为比例，它们表示表达强度
    started <- proc.time()[['elapsed']]
    bootstrap_coefficients <- matrix(NA_real_, bootstrap_global, type_count)
    for (b in seq_len(bootstrap_global)) {
        sampled_rows <- sample.int(sample_count, sample_count, replace = TRUE)
        fit_bootstrap <- nnls(fractions[sampled_rows, , drop = FALSE], y[sampled_rows])
        stopifnot(fit_bootstrap$mode == 1)
        bootstrap_coefficients[b, ] <- coef(fit_bootstrap)
    }
    global_median[output_row, ] <- apply(bootstrap_coefficients, 2, median)
    bootstrap_seconds <- proc.time()[['elapsed']] - started

    # 对 bulk 排序并同步排列比例；排序窗口不使用疾病、炎症或区域标签
    started <- proc.time()[['elapsed']]
    sample_order <- order(y, seq_along(y))
    sorted_y <- y[sample_order]
    sorted_fractions <- fractions[sample_order, , drop = FALSE]
    windows <- matrix(NA_real_, window_count, type_count)
    window_ranks <- integer(window_count)
    for (w in seq_len(window_count)) {
        window_rows <- seq.int(w, w + window_size - 1L)
        window_ranks[w] <- qr(sorted_fractions[window_rows, , drop = FALSE])$rank
        window_coefficients <- matrix(NA_real_, bootstrap_window, type_count)
        for (b in seq_len(bootstrap_window)) {
            sampled_rows <- sample(window_rows, window_size, replace = TRUE)
            fit_window <- nnls(sorted_fractions[sampled_rows, , drop = FALSE], sorted_y[sampled_rows])
            stopifnot(fit_window$mode == 1)
            window_coefficients[b, ] <- coef(fit_window)
        }
        windows[w, ] <- colMeans(window_coefficients)
    }
    forward <- rbind(windows, windows[rep(window_count, window_size - 1L), , drop = FALSE])
    backward <- rbind(windows[rep(1L, window_size - 1L), , drop = FALSE], windows)
    smoothed_sorted <- (forward + backward) / 2
    core_smoothed[output_row, , ] <- smoothed_sorted[order(sample_order), , drop = FALSE]
    stopifnot(nrow(smoothed_sorted) == sample_count, all(is.finite(smoothed_sorted)), all(smoothed_sorted >= 0))
    window_seconds <- proc.time()[['elapsed']] - started
    timings[[output_row]] <- data.frame(GeneSymbol = gene, global_nnls_seconds = global_seconds, bootstrap_seconds = bootstrap_seconds, window_seconds = window_seconds, windows = window_count, min_window_rank = min(window_ranks), nnls_fits = 1L + bootstrap_global + window_count * bootstrap_window)
    saveRDS(list(method = 'Independent uncalibrated NNLS window core', gene = gene, samples = colnames(mixture), types = colnames(fractions), global_nnls = global_nnls[output_row, ], global_bootstrap_median = global_median[output_row, ], smoothed = core_smoothed[output_row, , ], window_ranks = window_ranks, seed = 20261011L + selection$full_mixture_row_1based[gene_row]), file.path(outdir, sprintf('gene_%05d.rds', selection$full_mixture_row_1based[gene_row])))
    print(timings[[output_row]])
}
write.table(global_nnls, file.path(outdir, 'core_global_nnls.tsv'), sep = '\t', quote = FALSE, col.names = NA)
write.table(global_median, file.path(outdir, 'core_global_bootstrap_median.tsv'), sep = '\t', quote = FALSE, col.names = NA)
timings <- do.call(rbind, timings)
kkt_records <- do.call(rbind, kkt_records)
write.table(timings, file.path(outdir, 'timings.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(kkt_records, file.path(outdir, 'real_gene_kkt_checks.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
stopifnot(max(kkt_records$nonnegative_violation) == 0, max(kkt_records$dual_violation_scaled) < 1e-8, max(kkt_records$complementarity_scaled) < 1e-8)

# Step4. 比较官方输出；没有校准与过滤的滑窗差异不能解释成完整重构成功或失败
comparison <- vector('list', type_count)
for (j in seq_len(type_count)) {
    cell_type <- colnames(fractions)[j]
    predicted <- core_smoothed[, , j, drop = FALSE]
    dim(predicted) <- c(length(gene_rows), sample_count)
    dimnames(predicted) <- list(rownames(global_nnls), colnames(mixture))
    write.table(predicted, file.path(outdir, paste0('CORE_ONLY_', cell_type, '.tsv')), sep = '\t', quote = FALSE, col.names = NA)
    official <- as.matrix(read.delim(file.path(input_dir, paste0('official_hires_', cell_type, '.tsv')), row.names = 1, check.names = FALSE))
    official <- official[rownames(predicted), , drop = FALSE]
    stopifnot(identical(colnames(official), colnames(predicted)))
    complete <- rowSums(is.na(official)) == 0
    varying <- complete & apply(official, 1, function(x) if (anyNA(x)) FALSE else diff(range(x)) > 1e-8)
    correlations <- rep(NA_real_, nrow(official))
    for (g in which(varying)) {
        if (sd(predicted[g, ]) > 0) correlations[g] <- cor(log2(predicted[g, ] + 1), log2(official[g, ] + 1))
    }
    errors <- abs(log2(predicted + 1) - log2(official + 1))
    comparison[[j]] <- data.frame(CellType = cell_type, n_official_complete = sum(complete), n_official_varying = sum(varying), n_finite_correlations = sum(is.finite(correlations)), median_variable_gene_log2_correlation = if (any(is.finite(correlations))) median(correlations, na.rm = TRUE) else NA_real_, median_log2_absolute_error_complete = if (any(complete)) median(errors[complete, , drop = FALSE]) else NA_real_)
}
comparison <- do.call(rbind, comparison)
write.table(comparison, file.path(outdir, 'core_vs_official_hires.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
group_comparison <- data.frame(CellType = colnames(fractions), direct_nnls_log2_mae = colMeans(abs(log2(global_nnls + 1) - log2(official_group[gene_rows, , drop = FALSE] + 1))), bootstrap_median_log2_mae = colMeans(abs(log2(global_median + 1) - log2(official_group[gene_rows, , drop = FALSE] + 1))))
write.table(group_comparison, file.path(outdir, 'core_vs_official_group.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(c('Independent research prototype: numerical core only; NOT official CIBERSORTx or a validated HiRes replacement.', paste('Mixture input:', mixture_path), 'Official fixed weights reused; upstream fractions and S-mode are not implemented.', 'No global gene-wide BH, detectability filter, two-segment refinement, OLS calibration or official output filter.', 'Independent R seeds are fixed per original gene row; they do not reproduce the vendor random stream.', 'CORE_ONLY files must not be passed into formal EcoTyper discovery/recovery.', 'Gene checkpoints preserve all samples. No sample batches, automatic cache reuse or SSH operations.'), file.path(outdir, 'METHOD_LIMITATIONS.txt'))
print(comparison)
print(colSums(timings[, c('global_nnls_seconds', 'bootstrap_seconds', 'window_seconds', 'nnls_fits')]))
