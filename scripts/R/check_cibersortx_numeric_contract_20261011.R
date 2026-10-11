# 独立核对矩阵方向、列优先展开、窗口公式和非负最小二乘的数学条件。
# 本机基础 R 运行，只用人造数据；输出核对表，不运行供应商代码，也不验证与 CIBERSORTx 数值等价。
args <- commandArgs(TRUE)
outdir <- args[1]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Step1. 一基因在八个样本中的表达 = 样本 × 三类细胞比例矩阵乘以三类表达
fractions <- rbind(c(.7, .2, .1), c(.2, .7, .1), c(.1, .2, .7), c(.4, .4, .2), c(.3, .2, .5), c(.2, .3, .5), c(.6, .1, .3), c(.1, .6, .3))
colnames(fractions) <- c('A', 'B', 'C')
expression_truth <- c(10, 30, 5)
mixture_gene <- drop(fractions %*% expression_truth)
print(dim(fractions))
print(head(mixture_gene))

# Step2. C++ 的 NNLS 封装按列展开；R 默认矩阵向量化也按列展开
packed <- as.vector(fractions)
unpacked <- matrix(packed, nrow(fractions), ncol(fractions))
stopifnot(identical(unname(fractions), unpacked))
subset_indices <- c(1, 1, 3, 5, 7, 8) # 有放回抽样的重复行同时作用于比例和混合表达
fractions_subset <- fractions[subset_indices, , drop = FALSE]
mixture_subset <- mixture_gene[subset_indices]
subset_estimate <- drop(qr.solve(fractions_subset, mixture_subset))
stopifnot(max(abs(subset_estimate - expression_truth)) < 1e-10)

# Step3. 仅三类人造数据，穷举活跃集合检验 NNLS 数学约束；不用它处理真实研究矩阵
mixture_boundary <- drop(fractions %*% c(10, -5, 5))
active_sets <- list(integer(0))
for (active_count in 1:ncol(fractions)) {
    active_sets <- c(active_sets, combn(seq_len(ncol(fractions)), active_count, simplify = FALSE))
}
nnls_estimate <- rep(0, ncol(fractions))
best_rss <- sum(mixture_boundary^2)
for (active in active_sets) { # 此循环是小矩阵 NNLS 的明确数学求解步骤
    candidate <- rep(0, ncol(fractions))
    if (length(active) > 0) candidate[active] <- qr.solve(fractions[, active, drop = FALSE], mixture_boundary)
    if (any(candidate < 0)) next
    candidate_residual <- drop(fractions %*% candidate) - mixture_boundary
    candidate_rss <- sum(candidate_residual^2)
    if (candidate_rss < best_rss) {
        best_rss <- candidate_rss
        nnls_estimate <- candidate
    }
}
nnls_residual <- drop(fractions %*% nnls_estimate) - mixture_boundary
gradient <- drop(crossprod(fractions, nnls_residual))
kkt_error <- max(abs(gradient[nnls_estimate > 1e-10]), pmax(-gradient[nnls_estimate <= 1e-10], 0))
stopifnot(all(nnls_estimate >= 0), kkt_error < 1e-10)
print(data.frame(cell_type = colnames(fractions), nnls_estimate, gradient))

# Step4. 核对固定镜像的自适应窗口与分割点；参数仅用于解释，不能视作网页实测值
cell_types <- 20
samples <- 120
window <- ceiling(max(cell_types + 1, min(samples / 2, 4 * cell_types)))
stopifnot(window == 60)
split_bound <- window
if (samples - 2 * window < 10) split_bound <- window / 2
split_candidates <- seq(min(split_bound, samples - split_bound), max(split_bound, samples - split_bound), by = 1)
stopifnot(split_bound == 30, length(split_candidates) == 61)
print(data.frame(samples, cell_types, window, split_bound, split_count = length(split_candidates)))

# 原 R 中 kappa(BG.comb.xpr, exact <- T) 的第二位置实参确实传给 exact；赋值写法不改变这一点
barcode <- cbind(c(1, 2, 3), c(2, 1, 4))
kappa_positional <- kappa(barcode, exact <- TRUE)
kappa_named <- kappa(barcode, exact = TRUE)
stopifnot(identical(kappa_positional, kappa_named))

# Step5. 保存数学核对范围，明确不能据此宣称反编译已正确或算法可用于正式分析
checks <- data.frame(check = c('column_major_packing', 'subset_row_alignment', 'nnls_kkt', 'default_window_120_20', 'kappa_argument'), passed = TRUE, value = c(max(abs(unpacked - fractions)), max(abs(subset_estimate - expression_truth)), kkt_error, window, kappa_named))
write.table(checks, file.path(outdir, 'independent_numeric_contract_checks.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(data.frame(cell_type = colnames(fractions), nnls_estimate, gradient), file.path(outdir, 'toy_nnls_solution.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, 'numeric_contract_R_sessionInfo.txt'))
print(checks)
