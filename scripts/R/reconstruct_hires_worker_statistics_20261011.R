# 重建旧 HiRes 的可检测性、两段拟合和差异统计，在显式固定基因集合上保存逐步结果。
# 在 nnls 的 R 环境运行：Rscript 本脚本 基准目录 原始mixture子集 全基因统计RDS 新输出目录。
# 输入比例来自官方；这里只实现统计阶段，不产生可用于正式 EcoTyper 的 HiRes 表达矩阵。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
library(nnls)
stopifnot(length(args) == 4, !dir.exists(args[4]))
input_dir <- args[1]
outdir <- args[4]
dir.create(outdir, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))

# Step1. 固定样本和细胞类型，使用完整基因范围已计算的入口q值
mixture <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
fractions <- as.matrix(read.delim(file.path(input_dir, 'fractions_fixed.tsv'), row.names = 1, check.names = FALSE))
selection <- read.delim(file.path(input_dir, 'gene_selection.tsv'), check.names = FALSE)
entry <- readRDS(args[3])
stopifnot(identical(colnames(mixture), rownames(fractions)), identical(entry$samples, colnames(mixture)), identical(entry$types, colnames(fractions)))
stopifnot(all(rownames(mixture) %in% rownames(entry$worker_qvalues)), all(is.finite(mixture)), all(mixture >= 0))
stopifnot(identical(selection$GeneSymbol, rownames(mixture)))  # 防止基因重排后使用别的基因的种子
qvalues <- entry$worker_qvalues[rownames(mixture), colnames(fractions), drop = FALSE]
sample_count <- nrow(fractions)
type_count <- ncol(fractions)
window_size <- ceiling(max(type_count + 1, min(sample_count / 2, 4 * type_count)))
split_boundary <- if (sample_count - 2 * window_size < 10) window_size / 2 else window_size
candidate_cuts <- seq.int(ceiling(min(split_boundary, sample_count - split_boundary)), floor(max(split_boundary, sample_count - split_boundary)))
print(c(genes = nrow(mixture), samples = sample_count, types = type_count, window = window_size, split_min = min(candidate_cuts), split_max = max(candidate_cuts)))

# Step2. 用已知两段生成值核对分界拟合；不靠官方输出来挑选分界
known_cut <- sample_count %/% 2L
known_left <- seq_len(type_count) * 3
known_right <- rev(known_left) * 2
known_y <- c(as.vector(fractions[seq_len(known_cut), , drop = FALSE] %*% known_left), as.vector(fractions[seq.int(known_cut + 1L, sample_count), , drop = FALSE] %*% known_right))
known_rss <- numeric(length(candidate_cuts))
for (i in seq_along(candidate_cuts)) {
    cut <- candidate_cuts[i]
    fit_left <- nnls(fractions[seq_len(cut), , drop = FALSE], known_y[seq_len(cut)])
    fit_right <- nnls(fractions[seq.int(cut + 1L, sample_count), , drop = FALSE], known_y[seq.int(cut + 1L, sample_count)])
    known_rss[i] <- sum(fit_left$residuals^2) + sum(fit_right$residuals^2)
}
known_best_cut <- candidate_cuts[which.min(known_rss)]
stopifnot(known_best_cut == known_cut, min(known_rss) < 1e-16)
write.table(data.frame(expected_cut = known_cut, recovered_cut = known_best_cut, residual_sum_squares = min(known_rss), passed = TRUE), file.path(outdir, 'known_two_segment_check.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

# Step3. 每基因先完成全局和分界估计，保留抽样系数以便独立检查统计量
bootstrap_count <- 100L  # 旧 worker 全局及每段有放回重抽样次数
random_count <- 50L  # 旧 worker 的列置换次数
random_bootstrap_count <- 5L  # 每个置换的左右段各5次有放回抽样
type_records <- vector('list', nrow(mixture))
gene_records <- vector('list', nrow(mixture))
solver_count <- 0L
started <- proc.time()[['elapsed']]
for (g in seq_len(nrow(mixture))) {
    gene <- rownames(mixture)[g]
    set.seed(20261011L + selection$full_mixture_row_1based[g])  # 项目独立随机流，不冒充官方线程随机流
    y <- mixture[g, ]
    global_coefficients <- matrix(NA_real_, bootstrap_count, type_count)
    for (b in seq_len(bootstrap_count)) {
        sampled <- sample.int(sample_count, sample_count, replace = TRUE)
        fitted <- nnls(fractions[sampled, , drop = FALSE], y[sampled])
        stopifnot(fitted$mode == 1)
        global_coefficients[b, ] <- coef(fitted)
    }
    solver_count <- solver_count + bootstrap_count
    global_beta <- apply(global_coefficients, 2, median)
    probability_nonpositive <- colMeans(global_coefficients <= 0)
    probability_nonpositive[probability_nonpositive == 0] <- 1 / bootstrap_count
    entry_q <- qvalues[g, ]
    entry_q_for_fisher <- entry_q
    entry_q_for_fisher[!is.na(entry_q_for_fisher) & entry_q_for_fisher == 0] <- 1e-100
    global_p <- pchisq(-2 * (log(probability_nonpositive) + log(entry_q_for_fisher)), df = 4, lower.tail = FALSE)
    global_p[is.na(entry_q) | entry_q == 1] <- 1

    sample_order <- order(y, seq_along(y))
    sorted_y <- y[sample_order]
    sorted_fractions <- fractions[sample_order, , drop = FALSE]
    split_rss <- numeric(length(candidate_cuts))
    for (i in seq_along(candidate_cuts)) {
        cut <- candidate_cuts[i]
        fit_left <- nnls(sorted_fractions[seq_len(cut), , drop = FALSE], sorted_y[seq_len(cut)])
        fit_right <- nnls(sorted_fractions[seq.int(cut + 1L, sample_count), , drop = FALSE], sorted_y[seq.int(cut + 1L, sample_count)])
        stopifnot(fit_left$mode == 1, fit_right$mode == 1)
        split_rss[i] <- sum(fit_left$residuals^2) + sum(fit_right$residuals^2)
    }
    solver_count <- solver_count + 2L * length(candidate_cuts)
    best_cut <- ceiling(max(split_boundary, candidate_cuts[which.min(split_rss)]))
    segment_coefficients <- array(NA_real_, c(bootstrap_count, type_count, 2))
    segment_rows <- list(seq_len(best_cut), seq.int(best_cut + 1L, sample_count))
    for (segment in 1:2) {
        rows <- segment_rows[[segment]]
        for (b in seq_len(bootstrap_count)) {
            sampled <- rows[sample.int(length(rows), length(rows), replace = TRUE)]
            fitted <- nnls(sorted_fractions[sampled, , drop = FALSE], sorted_y[sampled])
            stopifnot(fitted$mode == 1)
            segment_coefficients[b, , segment] <- coef(fitted)
        }
    }
    solver_count <- solver_count + 2L * bootstrap_count
    segment_means <- apply(segment_coefficients, c(2, 3), mean)
    segment_sds <- apply(segment_coefficients, c(2, 3), sd)
    difference <- segment_means[, 2] - segment_means[, 1]
    difference_z <- difference / sqrt(segment_sds[, 1]^2 + segment_sds[, 2]^2)
    difference_z[is.nan(difference_z)] <- 0
    z_p <- 2 * pnorm(-abs(difference_z))
    z_p[!is.finite(z_p)] <- 1
    z_p[z_p == 0] <- 1e-100
    detectable <- min(global_p) < 0.05
    random_differences <- matrix(NA_real_, random_count, type_count)
    empirical_p <- rep(NA_real_, type_count)
    differential_p <- empirical_p

    # 原数值实现逐列打乱比例；不使用同一整行置换，不重新归一化比例行
    if (detectable) {
        for (r in seq_len(random_count)) {
            random_left <- matrix(NA_real_, random_bootstrap_count, type_count)
            rows <- segment_rows[[1]]
            for (b in seq_len(random_bootstrap_count)) {
                sampled <- rows[sample.int(length(rows), length(rows), replace = TRUE)]
                fitted <- nnls(sorted_fractions[sampled, , drop = FALSE], sorted_y[sampled])
                stopifnot(fitted$mode == 1)
                random_left[b, ] <- coef(fitted)
            }
            shuffled_fractions <- sorted_fractions
            for (j in seq_len(type_count)) shuffled_fractions[, j] <- sample(sorted_fractions[, j], sample_count, replace = FALSE)
            random_right <- random_left
            rows <- segment_rows[[2]]
            for (b in seq_len(random_bootstrap_count)) {
                sampled <- rows[sample.int(length(rows), length(rows), replace = TRUE)]
                fitted <- nnls(shuffled_fractions[sampled, , drop = FALSE], sorted_y[sampled])
                stopifnot(fitted$mode == 1)
                random_right[b, ] <- coef(fitted)
            }
            random_differences[r, ] <- colMeans(random_right) - colMeans(random_left)
        }
        solver_count <- solver_count + random_count * 2L * random_bootstrap_count
        for (j in seq_len(type_count)) {
            ordered_null <- sort(random_differences[, j], decreasing = difference[j] < 0)
            closest_index <- which.min(abs(ordered_null - difference[j]))
            empirical_p[j] <- 1 - closest_index / random_count
        }
        empirical_p[empirical_p == 0] <- 1 / random_count
        differential_p <- pchisq(-2 * (log(empirical_p) + log(z_p)), df = 4, lower.tail = FALSE)
    }

    type_records[[g]] <- data.frame(GeneSymbol = gene, CellType = colnames(fractions), global_beta = global_beta, entry_q = entry_q, probability_nonpositive = probability_nonpositive, global_p = global_p, segment1_mean = segment_means[, 1], segment2_mean = segment_means[, 2], segment1_sd = segment_sds[, 1], segment2_sd = segment_sds[, 2], difference = difference, difference_z = difference_z, z_p = z_p, empirical_p = empirical_p, differential_p = differential_p, row.names = NULL)
    gene_records[[g]] <- data.frame(GeneSymbol = gene, best_cut = best_cut, minimum_split_rss = min(split_rss), gene_detectable = detectable, globally_detectable_types = sum(global_p < 0.05), differentially_supported_types = sum(differential_p < 0.05, na.rm = TRUE))
    saveRDS(list(gene = gene, sample_order = sample_order, global_coefficients = global_coefficients, candidate_cuts = candidate_cuts, split_rss = split_rss, best_cut = best_cut, segment_coefficients = segment_coefficients, random_differences = random_differences, entry_q = entry_q, global_p = global_p, differential_p = differential_p, stage = 'Independent worker statistics only'), file.path(outdir, sprintf('gene_%04d_statistics.rds', g)))
    print(c(completed_genes = g, elapsed_seconds = proc.time()[['elapsed']] - started))
}

# Step4. 导出可检查的统计阶段；无最终表达、NA掩码或正式G矩阵
type_statistics <- do.call(rbind, type_records)
gene_statistics <- do.call(rbind, gene_records)
write.table(type_statistics, file.path(outdir, 'worker_type_statistics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
write.table(gene_statistics, file.path(outdir, 'worker_gene_statistics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
summary <- data.frame(genes = nrow(mixture), samples = sample_count, types = type_count, detectable_genes = sum(gene_statistics$gene_detectable), nnls_fits = solver_count, elapsed_seconds = proc.time()[['elapsed']] - started)
write.table(summary, file.path(outdir, 'worker_statistics_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(c('Independent reconstruction of old numerical worker statistics; not current website parity.', 'Official fixed fractions, original linear mixture, full gene-space entry Q values.', 'Independent RNG; Fisher combines dependent statistics and is reproduced as an algorithmic rule, not a new inferential justification.', 'Missing residual Spearman fallback, shrinkage, targeted re-fit, sliding calibration, rare-type override and final export masking.', 'These statistics do not constitute a discovery/recovery result or a formal HiRes matrix.'), file.path(outdir, 'METHOD_LIMITATIONS.txt'))
print(summary)
