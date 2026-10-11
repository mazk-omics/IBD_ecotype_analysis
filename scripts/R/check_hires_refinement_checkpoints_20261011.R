# 检查HiRes收缩与校准检查点，导出独立SciPy复核所需的小输入，不覆盖表达结果。
# 在本机R环境运行：Rscript 本脚本 基准目录 原始24基因mixture 收缩结果目录 新检查目录。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
stopifnot(length(args) == 4, !dir.exists(args[4]))
outdir <- args[4]
dir.create(outdir, recursive = TRUE)
mixture <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
fractions <- as.matrix(read.delim(file.path(args[1], 'fractions_fixed.tsv'), row.names = 1, check.names = FALSE))
type_statistics <- read.delim(file.path(args[3], 'refinement_type_statistics.tsv'), check.names = FALSE)

# Step1. 小无并列样例供Python穷举置换，真实残差统计另与R公开cor.test比较
toy_y <- list(c(1, 3, 2, 4, 6, 5, 7, 8), c(8, 6, 7, 5, 3, 4, 2, 1), c(1, 2, 3, 4, 5, 6, 7, 8))
toy_records <- vector('list', length(toy_y))
for (i in seq_along(toy_y)) {
    y <- toy_y[[i]]
    rho <- cor(seq_along(y), rank(y))
    n <- length(y)
    s <- (1 - rho) * (n^3 - n) / 6
    lower_tail <- s <= (n^3 - n) / 6
    p <- min(2 * .Call(stats:::C_pRho, round(s) + 2 * lower_tail, as.integer(n), lower_tail), 1)
    toy_records[[i]] <- data.frame(example = i, y_order = paste(y, collapse = ','), rho = rho, s = s, two_sided_p = p)
}
write.table(do.call(rbind, toy_records), file.path(outdir, 'spearman_toys_for_enumeration.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

# Step2. 验证真实检查点，保留分界、弱类型与OLS输入供另一语言求解
correlation_records <- list()
shrink_records <- list()
refit_inputs <- list()
weak_inputs <- list()
weak_coefficients <- list()
calibration_inputs <- list()
for (g in seq_len(nrow(mixture))) {
    checkpoint <- readRDS(file.path(args[3], sprintf('gene_%04d_refinement.rds', g)))
    gene <- rownames(mixture)[g]
    stopifnot(identical(checkpoint$gene, gene))
    if (!checkpoint$gene_detectable) next
    info <- type_statistics[type_statistics$GeneSymbol == gene, ]
    stopifnot(identical(info$CellType, colnames(fractions)))
    sorted_fractions <- fractions[checkpoint$sample_order, , drop = FALSE]
    sorted_y <- as.vector(mixture[g, checkpoint$sample_order])
    for (j in seq_len(ncol(fractions))) {
        others <- setdiff(seq_len(ncol(fractions)), j)
        residual <- sorted_y - as.vector(sorted_fractions[, others, drop = FALSE] %*% info$global_initial[others])
        x <- sorted_fractions[, j]
        expected_p <- 1
        if (sd(rank(x)) > 0 && sd(rank(residual)) > 0) {
            no_ties <- !anyDuplicated(x) && !anyDuplicated(residual)
            tested <- cor.test(x, residual, method = 'spearman', exact = nrow(fractions) <= 1290 && no_ties)
            expected_p <- if (tested$p.value == 0) 1e-100 else tested$p.value
        }
        # 输出表中缺少检测证据的类型p已设为1，比较同一阶段语义
        if (checkpoint$global_p[j] > 0.05) expected_p <- 1
        error <- abs(expected_p - info$residual_p[j])
        correlation_records[[length(correlation_records) + 1L]] <- data.frame(GeneSymbol = gene, CellType = colnames(fractions)[j], branch = info$residual_branch[j], error = error, passed = error < 1e-12)
    }
    midpoint_error <- max(abs((checkpoint$adjusted1_before_weak + checkpoint$adjusted2_before_weak) / 2 - (checkpoint$segment1_before_shrink + checkpoint$segment2_before_shrink) / 2))
    difference_error <- max(abs(abs(checkpoint$adjusted1_before_weak - checkpoint$adjusted2_before_weak) - abs(checkpoint$segment1_before_shrink - checkpoint$segment2_before_shrink) * checkpoint$shrink_weights^2))
    scale <- max(1, max(abs(checkpoint$segment1_before_shrink)), max(abs(checkpoint$segment2_before_shrink)))
    shrink_records[[length(shrink_records) + 1L]] <- data.frame(GeneSymbol = gene, midpoint_scaled_error = midpoint_error / scale, squared_weight_difference_scaled_error = difference_error / scale, passed = max(midpoint_error, difference_error) / scale < 1e-12)
    refit_inputs[[length(refit_inputs) + 1L]] <- data.frame(GeneSymbol = gene, sorted_row = seq_len(nrow(fractions)), mixture = checkpoint$variable_y, best_cut = checkpoint$candidate_cuts[which.min(checkpoint$refit_rss)], checkpoint$variable_fractions, check.names = FALSE)
    weak_inputs[[length(weak_inputs) + 1L]] <- data.frame(GeneSymbol = gene, sorted_row = seq_len(nrow(fractions)), residual = checkpoint$weak_residual, checkpoint$weak_fractions, check.names = FALSE)
    weak_coefficients[[length(weak_coefficients) + 1L]] <- data.frame(GeneSymbol = gene, CellType = colnames(fractions), coefficient = checkpoint$weak_fit_coefficients)
    for (j in checkpoint$calibrated_types) {
        calibration_inputs[[length(calibration_inputs) + 1L]] <- data.frame(GeneSymbol = gene, CellType = colnames(fractions)[j], sorted_row = seq_len(nrow(fractions)), smooth_x = checkpoint$smoothed[, j], target_y = checkpoint$anchors[, j], r_slope = checkpoint$slope[j], r_intercept = checkpoint$intercept[j])
    }
}
correlation_checks <- do.call(rbind, correlation_records)
shrink_checks <- do.call(rbind, shrink_records)
stopifnot(all(correlation_checks$passed), all(shrink_checks$passed))
write.table(correlation_checks, file.path(outdir, 'residual_spearman_vs_R_public.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(shrink_checks, file.path(outdir, 'shrink_invariants.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, refit_inputs), file.path(outdir, 'refit_inputs.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, weak_inputs), file.path(outdir, 'weak_inputs.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, weak_coefficients), file.path(outdir, 'weak_coefficients.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, calibration_inputs), file.path(outdir, 'calibration_inputs.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
print(c(residual_correlation_checks = nrow(correlation_checks), shrink_gene_checks = nrow(shrink_checks), calibration_pairs = length(calibration_inputs)))
