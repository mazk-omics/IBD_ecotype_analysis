# 接续固定基因集合的HiRes统计，重拟合、收缩、残差重估、滑窗与OLS；保存每步可核对对象。
# 在 nnls 的R环境运行：Rscript 本脚本 基准目录 原始mixture子集 统计目录 新输出目录 [阶段种子]。
# 这是旧数值路径的独立研究原型；最终网页NA掩码、辅助回填和官方随机流尚未复现。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
library(nnls)
stopifnot(length(args) %in% c(4L, 5L), !dir.exists(args[4]))
input_dir <- args[1]
statistics_dir <- args[3]
outdir <- args[4]
dir.create(outdir, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))

# Step1. 固定统计与表达的行列；这里不重新估计上游fraction
mixture <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
fractions <- as.matrix(read.delim(file.path(input_dir, 'fractions_fixed.tsv'), row.names = 1, check.names = FALSE))
selection <- read.delim(file.path(input_dir, 'gene_selection.tsv'), check.names = FALSE)
statistics <- read.delim(file.path(statistics_dir, 'worker_type_statistics.tsv'), check.names = FALSE)
gene_statistics <- read.delim(file.path(statistics_dir, 'worker_gene_statistics.tsv'), check.names = FALSE)
stopifnot(identical(rownames(mixture), selection$GeneSymbol), identical(rownames(mixture), gene_statistics$GeneSymbol), identical(colnames(mixture), rownames(fractions)))
sample_count <- nrow(fractions)
type_count <- ncol(fractions)
window_size <- ceiling(max(type_count + 1, min(sample_count / 2, 4 * type_count)))
window_count <- sample_count - window_size + 1L
split_boundary <- if (sample_count - 2 * window_size < 10) window_size / 2 else window_size
candidate_cuts <- seq.int(ceiling(min(split_boundary, sample_count - split_boundary)), floor(max(split_boundary, sample_count - split_boundary)))
alpha <- 0.05  # 固定旧数值路径，不为提高官方相关性调阈值
minimum_score <- -log10(alpha)
bootstrap_refit <- 100L
bootstrap_window <- 5L
stage_seed <- if (length(args) == 5L) as.integer(args[5]) else 20261012L  # 独立重抽样检查用，不改变统计门槛
fraction_means <- colMeans(fractions)
rarest_types <- which(fraction_means == min(fraction_means))
refined <- array(NA_real_, c(nrow(mixture), sample_count, type_count), dimnames = list(rownames(mixture), colnames(mixture), colnames(fractions)))
type_records <- vector('list', nrow(mixture))
gene_records <- vector('list', nrow(mixture))
solver_count <- 0L
started <- proc.time()[['elapsed']]
print(c(genes = nrow(mixture), samples = sample_count, types = type_count, window = window_size))

# Step2. 逐基因保留原始统计，在同一排序中完成后续计算
for (g in seq_len(nrow(mixture))) {
    gene <- rownames(mixture)[g]
    checkpoint <- readRDS(file.path(statistics_dir, sprintf('gene_%04d_statistics.rds', g)))
    stopifnot(identical(checkpoint$gene, gene))
    gene_types <- statistics[statistics$GeneSymbol == gene, ]
    stopifnot(identical(gene_types$CellType, colnames(fractions)))
    global_p <- gene_types$global_p
    differential_p <- gene_types$differential_p
    global_initial <- gene_types$global_beta
    sample_order <- checkpoint$sample_order
    sorted_y <- as.vector(mixture[g, sample_order])
    sorted_fractions <- fractions[sample_order, , drop = FALSE]
    set.seed(stage_seed + selection$full_mixture_row_1based[g])  # 独立阶段种子；不声称延续官方MT线程状态
    detectable <- min(global_p) < alpha
    global_refit <- global_initial
    segment1 <- segment2 <- global_initial
    shrink_weights <- rep(NA_real_, type_count)
    shrink_score <- shrink_weights
    residual_p <- rep(1, type_count)
    residual_rho <- rep(NA_real_, type_count)
    residual_branch <- rep('not_run_gene_undetectable', type_count)
    weak_types <- integer(0)
    rare_override <- integer(0)
    fallback_used <- FALSE
    refined_cut <- NA_integer_
    slope <- rep(NA_real_, type_count)
    intercept <- slope
    calibrated_types <- integer(0)
    final_sorted <- matrix(1, sample_count, type_count)
    checkpoint_refined <- list(gene = gene, sample_order = sample_order, global_p = global_p, differential_p = differential_p, gene_detectable = detectable, stage_seed = stage_seed)

    if (detectable) {
        # 残差相关使用最初全局系数，减去其他类型贡献；秩与并列值决定统计分支
        for (j in seq_len(type_count)) {
            other_types <- setdiff(seq_len(type_count), j)
            residual_y <- sorted_y - as.vector(sorted_fractions[, other_types, drop = FALSE] %*% global_initial[other_types])
            x_rank <- rank(sorted_fractions[, j], ties.method = 'average')
            y_rank <- rank(residual_y, ties.method = 'average')
            ties_present <- anyDuplicated(sorted_fractions[, j]) > 0 || anyDuplicated(residual_y) > 0
            if (sd(x_rank) == 0 || sd(y_rank) == 0) {
                residual_branch[j] <- 'constant_rank_p1'
            } else {
                rho <- cor(x_rank, y_rank)
                residual_rho[j] <- rho
                s_statistic <- (1 - rho) * (sample_count^3 - sample_count) / 6
                if (sample_count <= 1290L && !ties_present) {
                    lower_tail <- s_statistic <= (sample_count^3 - sample_count) / 6
                    # R的公开AS89例程：n<=1290且无ties；不调用供应商二进制
                    one_tail <- .Call(stats:::C_pRho, round(s_statistic) + 2 * lower_tail, as.integer(sample_count), lower_tail)
                    residual_p[j] <- min(2 * one_tail, 1)
                    residual_branch[j] <- 'R_AS89_no_ties'
                } else {
                    residual_p[j] <- 2 * pt(-abs(rho / sqrt((1 - rho^2) / (sample_count - 2))), df = sample_count - 2)
                    residual_branch[j] <- 'rank_t_with_ties_or_large_n'
                }
                if (!is.finite(residual_p[j])) residual_p[j] <- 1
                if (residual_p[j] == 0) residual_p[j] <- 1e-100
            }
        }

        # 屏蔽缺少检测证据的比例列，重估全局表达；最初系数为零的类型维持零
        active_types <- which(global_p < alpha)
        detectable_fractions <- sorted_fractions
        detectable_fractions[, global_p > alpha] <- 0
        refit_coefficients <- matrix(NA_real_, bootstrap_refit, type_count)
        for (b in seq_len(bootstrap_refit)) {
            sampled <- sample.int(sample_count, sample_count, replace = TRUE)
            fitted <- nnls(detectable_fractions[sampled, , drop = FALSE], sorted_y[sampled])
            stopifnot(fitted$mode == 1)
            refit_coefficients[b, ] <- coef(fitted)
        }
        solver_count <- solver_count + bootstrap_refit
        global_refit <- apply(refit_coefficients, 2, median)
        global_refit[global_initial == 0] <- 0
        detectable_fractions[, global_refit == 0] <- 0
        variable_fractions <- sorted_fractions
        variable_fractions[, setdiff(seq_len(type_count), active_types)] <- 0
        variable_fractions[, global_refit == 0] <- 0
        other_types <- setdiff(seq_len(type_count), active_types)
        variable_y <- sorted_y - as.vector(detectable_fractions[, other_types, drop = FALSE] %*% global_refit[other_types])

        # 先找新分界再重抽样，不能拿未经屏蔽的旧滑窗直接校准
        refit_rss <- numeric(length(candidate_cuts))
        for (i in seq_along(candidate_cuts)) {
            cut <- candidate_cuts[i]
            fit1 <- nnls(variable_fractions[seq_len(cut), , drop = FALSE], variable_y[seq_len(cut)])
            fit2 <- nnls(variable_fractions[seq.int(cut + 1L, sample_count), , drop = FALSE], variable_y[seq.int(cut + 1L, sample_count)])
            stopifnot(fit1$mode == 1, fit2$mode == 1)
            refit_rss[i] <- sum(fit1$residuals^2) + sum(fit2$residuals^2)
        }
        solver_count <- solver_count + 2L * length(candidate_cuts)
        refined_cut <- ceiling(max(split_boundary, candidate_cuts[which.min(refit_rss)]))
        segment_rows <- list(seq_len(refined_cut), seq.int(refined_cut + 1L, sample_count))
        segment_refits <- array(NA_real_, c(bootstrap_refit, type_count, 2))
        for (segment in 1:2) {
            rows <- segment_rows[[segment]]
            for (b in seq_len(bootstrap_refit)) {
                sampled <- rows[sample.int(length(rows), length(rows), replace = TRUE)]
                fitted <- nnls(variable_fractions[sampled, , drop = FALSE], variable_y[sampled])
                stopifnot(fitted$mode == 1)
                segment_refits[b, , segment] <- coef(fitted)
            }
        }
        solver_count <- solver_count + 2L * bootstrap_refit
        segment_means <- apply(segment_refits, c(2, 3), mean)
        segment1 <- segment_means[, 1]
        segment2 <- segment_means[, 2]
        segment1[other_types] <- global_refit[other_types]
        segment2[other_types] <- global_refit[other_types]
        segment1[global_p > alpha] <- 0
        segment2[global_p > alpha] <- 0
        checkpoint_refined$global_refit_coefficients <- refit_coefficients
        checkpoint_refined$candidate_cuts <- candidate_cuts
        checkpoint_refined$refit_rss <- refit_rss
        checkpoint_refined$segment_refits <- segment_refits
        checkpoint_refined$segment1_before_shrink <- segment1
        checkpoint_refined$segment2_before_shrink <- segment2

        # 旧路径先将不显著DE分数设为1；全部为1时采用残差Spearman证据
        shrink_score <- -log10(differential_p)
        shrink_score[global_p > alpha] <- 1
        residual_p[global_p > alpha] <- 1
        shrink_score[shrink_score < minimum_score] <- 1
        if (max(shrink_score) == 1) {
            shrink_score <- -log10(residual_p)
            fallback_used <- TRUE
        }
        shrink_weights <- rep(0, type_count)
        decreasing_types <- which(segment1 > segment2)
        increasing_types <- which(segment1 < segment2)
        decreasing_evidence <- decreasing_types[residual_p[decreasing_types] < alpha]
        increasing_evidence <- increasing_types[residual_p[increasing_types] < alpha]
        decreasing_max <- max(c(1, shrink_score[decreasing_evidence]))
        increasing_max <- max(c(1, shrink_score[increasing_evidence]))
        shrink_weights[decreasing_types] <- shrink_score[decreasing_types] / decreasing_max
        shrink_weights[increasing_types] <- shrink_score[increasing_types] / increasing_max
        shrink_weights <- pmin(shrink_weights, 1)
        stopifnot(all(is.finite(shrink_weights)), all(shrink_weights >= 0))
        delta <- (1 - shrink_weights^2) * abs(segment1 - segment2) / 2
        direction <- sign(segment2 - segment1)
        adjusted1 <- segment1 + direction * delta
        adjusted2 <- segment2 - direction * delta
        crossed <- (segment1 > segment2 & adjusted1 < adjusted2) | (segment1 < segment2 & adjusted1 > adjusted2)
        adjusted1[crossed] <- global_refit[crossed]
        adjusted2[crossed] <- global_refit[crossed]
        fold <- (pmax(adjusted1, adjusted2) + 1) / (pmin(adjusted1, adjusted2) + 1)
        weak_types <- which(fold < 1.1 & shrink_weights != 1 & global_refit > 0)
        checkpoint_refined$shrink_weights <- shrink_weights
        checkpoint_refined$shrink_delta <- delta
        checkpoint_refined$adjusted1_before_weak <- adjusted1
        checkpoint_refined$adjusted2_before_weak <- adjusted2
        adjusted1[weak_types] <- global_refit[weak_types]
        adjusted2[weak_types] <- global_refit[weak_types]

        # 稀有类型覆盖先记录选择，实际在OLS之后执行；严格按最小平均比例集合判定
        rare_direction <- sign(adjusted2[rarest_types[1]] - adjusted1[rarest_types[1]])
        higher_types <- which(fraction_means > min(fraction_means))
        if (any(rare_direction * sign(adjusted2[higher_types] - adjusted1[higher_types]) < 0)) rare_override <- rarest_types

        # 弱变化类型通过残差NNLS重新估计，不直接使用刚才的bootstrap中位数
        piecewise <- rbind(matrix(adjusted1, refined_cut, type_count, byrow = TRUE), matrix(adjusted2, sample_count - refined_cut, type_count, byrow = TRUE))
        weak_fractions <- sorted_fractions
        weak_fractions[, setdiff(seq_len(type_count), weak_types)] <- 0
        fixed_fractions <- sorted_fractions
        fixed_fractions[, weak_types] <- 0
        weak_residual <- sorted_y - rowSums(fixed_fractions * piecewise)
        weak_fit <- nnls(weak_fractions, weak_residual)
        stopifnot(weak_fit$mode == 1)
        solver_count <- solver_count + 1L
        global_refit[weak_types] <- coef(weak_fit)[weak_types]
        adjusted1[weak_types] <- global_refit[weak_types]
        adjusted2[weak_types] <- global_refit[weak_types]
        piecewise <- rbind(matrix(adjusted1, refined_cut, type_count, byrow = TRUE), matrix(adjusted2, sample_count - refined_cut, type_count, byrow = TRUE))

        # 滑窗重新使用可检测类型的比例；固定全部样本，窗口系数不归一化成比例
        window_fractions <- sorted_fractions
        window_fractions[, setdiff(seq_len(type_count), active_types)] <- 0
        windows <- matrix(NA_real_, window_count, type_count)
        for (w in seq_len(window_count)) {
            rows <- seq.int(w, w + window_size - 1L)
            window_coefficients <- matrix(NA_real_, bootstrap_window, type_count)
            for (b in seq_len(bootstrap_window)) {
                sampled <- rows[sample.int(length(rows), length(rows), replace = TRUE)]
                fitted <- nnls(window_fractions[sampled, , drop = FALSE], variable_y[sampled])
                stopifnot(fitted$mode == 1)
                window_coefficients[b, ] <- coef(fitted)
            }
            windows[w, ] <- colMeans(window_coefficients)
        }
        solver_count <- solver_count + window_count * bootstrap_window
        forward <- rbind(windows, windows[rep(window_count, window_size - 1L), , drop = FALSE])
        backward <- rbind(windows[rep(1L, window_size - 1L), , drop = FALSE], windows)
        smoothed <- (forward + backward) / 2

        # 用两段表达作y、平滑表达作x校准幅度；常数x显式退化为锚点均值
        anchors <- piecewise + matrix(rnorm(sample_count * type_count, mean = 0, sd = 0.001), sample_count, type_count, byrow = TRUE)
        final_sorted <- piecewise
        calibrated_types <- which(abs(adjusted1 - adjusted2) > 1e-5)
        for (j in calibrated_types) {
            x <- smoothed[, j]
            target <- anchors[, j]
            slope[j] <- if (sd(x) == 0 || sd(target) == 0) 0 else cor(x, target) * sd(target) / sd(x)
            intercept[j] <- mean(target) - slope[j] * mean(x)
            final_sorted[, j] <- intercept[j] + slope[j] * x
        }
        final_sorted <- pmax(final_sorted, 1)
        if (length(rare_override)) final_sorted[, rare_override] <- matrix(global_refit[rare_override], sample_count, length(rare_override), byrow = TRUE)
        checkpoint_refined$variable_fractions <- variable_fractions
        checkpoint_refined$variable_y <- variable_y
        checkpoint_refined$active_types <- active_types
        checkpoint_refined$weak_types <- weak_types
        checkpoint_refined$weak_fractions <- weak_fractions
        checkpoint_refined$weak_residual <- weak_residual
        checkpoint_refined$weak_fit_coefficients <- coef(weak_fit)
        checkpoint_refined$piecewise <- piecewise
        checkpoint_refined$windows <- windows
        checkpoint_refined$smoothed <- smoothed
        checkpoint_refined$anchors <- anchors
        checkpoint_refined$calibrated_types <- calibrated_types
        checkpoint_refined$slope <- slope
        checkpoint_refined$intercept <- intercept
        checkpoint_refined$rare_override <- rare_override
    }

    # 外层表达下限与原始样本顺序；这里未套用未知的网页NA输出掩码
    final_original_order <- pmax(final_sorted[order(sample_order), , drop = FALSE], 1)
    stopifnot(all(is.finite(final_original_order)), all(final_original_order >= 1))
    refined[g, , ] <- final_original_order
    checkpoint_refined$final <- final_original_order
    checkpoint_refined$global_refit <- global_refit
    saveRDS(checkpoint_refined, file.path(outdir, sprintf('gene_%04d_refinement.rds', g)))
    type_records[[g]] <- data.frame(GeneSymbol = gene, CellType = colnames(fractions), global_initial = global_initial, global_refit = global_refit, residual_rho = residual_rho, residual_p = residual_p, residual_branch = residual_branch, shrink_score = shrink_score, shrink_weight = shrink_weights, segment1_final = if (detectable) adjusted1 else segment1, segment2_final = if (detectable) adjusted2 else segment2, weak_type_refitted = seq_len(type_count) %in% weak_types, rare_type_override = seq_len(type_count) %in% rare_override, slope = slope, intercept = intercept, row.names = NULL)
    gene_records[[g]] <- data.frame(GeneSymbol = gene, detectable = detectable, initial_cut = checkpoint$best_cut, refined_cut = refined_cut, spearman_fallback_used = fallback_used, weak_types_refitted = length(weak_types), calibrated_types = length(calibrated_types), rare_types_overridden = length(rare_override))
    if (g == 1L || g %% 50L == 0L || g == nrow(mixture)) print(c(completed_genes = g, elapsed_seconds = proc.time()[['elapsed']] - started))
}

# Step3. 新文件名明确标记研究原型，避免误当成官方HiRes或正式G
for (j in seq_len(type_count)) {
    matrix_out <- refined[, , j]
    dimnames(matrix_out) <- list(rownames(mixture), colnames(mixture))
    write.table(matrix_out, file.path(outdir, paste0('REFINEMENT_ONLY_', colnames(fractions)[j], '.tsv')), sep = '\t', quote = FALSE, col.names = NA)
}
write.table(do.call(rbind, type_records), file.path(outdir, 'refinement_type_statistics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
write.table(do.call(rbind, gene_records), file.path(outdir, 'refinement_gene_statistics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
summary <- data.frame(genes = nrow(mixture), samples = sample_count, types = type_count, nnls_fits = solver_count, stage_seed = stage_seed, elapsed_seconds = proc.time()[['elapsed']] - started)
write.table(summary, file.path(outdir, 'refinement_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(c('Independent old HiRes numerical-path prototype, not current official software equivalence.', 'Reads official fixed fractions and independently reconstructed entry/worker statistics.', 'R AS89 or rank-t residual Spearman branches; vendor prho executable parity not tested.', 'Stage-specific R seeds, different from official threaded MT or MATCH_R streams.', 'No webpage NA/reliability mask, optional range filter or auxiliary Group expression refill.', 'No independent signature, fractions or S-mode. Do not use in formal EcoTyper discovery/recovery.'), file.path(outdir, 'METHOD_LIMITATIONS.txt'))
print(summary)
