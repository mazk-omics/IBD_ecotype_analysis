# 核对原论文收缩公式与静态实现中 +1 倍数的边界；使用人造值，不拟合研究样本。
# 本机 R 运行，结果保存到源码审计目录；这是数学核对，不是官方程序等价性验证。
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

# Step1. 正反方向都向中点收缩，收缩后两段差应等于原差乘以权重平方
g1 <- c(2, 10, 5, 20)
g2 <- c(10, 2, 5, 4)
weight <- c(0, 0.5, 1, 1)
delta <- (1 - weight^2) * abs(g1 - g2) / 2
g1_adjusted <- g1 + sign(g2 - g1) * delta
g2_adjusted <- g2 + sign(g1 - g2) * delta
shrink_error <- max(abs(abs(g2_adjusted - g1_adjusted) - abs(g2 - g1) * weight^2))
midpoint_error <- max(abs((g1_adjusted + g2_adjusted) / 2 - (g1 + g2) / 2))
print(data.frame(g1, g2, weight, delta, g1_adjusted, g2_adjusted))

# Step2. +1 会改变低表达的替换判定，必须保留精确不等式而非随意近似
low <- c(0, 0.1, 1, 100, 1)
high <- c(0.05, 0.15, 1.15, 115, 1.2)
fold_pseudocount <- (high + 1) / (low + 1)
replace_weak <- fold_pseudocount < 1.1
expected_replace <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
boundary_error <- max(abs(as.integer(replace_weak) - as.integer(expected_replace)))
print(data.frame(low, high, fold_pseudocount, replace_weak))

# Step3. Pearson 与 Spearman 对非线性表达关系可得不同结果，不能按论文名称替换旧实现
x <- seq_len(12)
y <- exp(x / 3)
pearson <- cor.test(x, y, method = 'pearson')$p.value
spearman <- cor.test(x, y, method = 'spearman', exact = FALSE)$p.value
print(c(pearson_p = pearson, spearman_p = spearman))
stopifnot(shrink_error < 1e-12, midpoint_error < 1e-12, boundary_error == 0, abs(pearson - spearman) > 0)
checks <- data.frame(check = c('squared_weight_difference', 'midpoint_preserved', 'pseudocount_strict_boundary', 'pearson_spearman_not_interchangeable'), passed = TRUE, numeric_evidence = c(shrink_error, midpoint_error, boundary_error, abs(pearson - spearman)))
write.table(checks, file.path(outdir, 'independent_paper_refinement_checks.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(data.frame(g1, g2, weight, delta, g1_adjusted, g2_adjusted), file.path(outdir, 'paper_refinement_toy_values.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
print(checks)
