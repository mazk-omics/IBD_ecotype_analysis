# 独立核对 HiRes 的统计辅助、窗口边缘与样本回排规则，不载入或执行供应商代码。
# 本机 R 4.4.1 运行；人造数据不接触研究样本，结果写入源码审计目录，仅证明解释的数学自洽性。
outdir <- "D:/IBD_ecotype/IBD_ecotype_analysis/CIBERSORTx_server_plan_20261011"
checks <- data.frame(Check=character(), Error=numeric(), Passed=logical())

# Step1. Bootstrap 非正表达频率为统计辅助值；零频率在调用者替换为 1/B
bootstrap_expression <- cbind(c(0, 2, 3, 0, 4), c(1, 2, 3, 4, 5), c(0, 0, 0, 0, 0))
bootstrap_probability <- colMeans(bootstrap_expression <= 0)
bootstrap_probability[bootstrap_probability <= 0] <- 1 / nrow(bootstrap_expression)
probability_error <- max(abs(bootstrap_probability - c(0.4, 0.2, 1)))
checks <- rbind(checks, data.frame(Check="bootstrap_nonpositive_frequency", Error=probability_error, Passed=probability_error < 1e-12))
print(bootstrap_probability)

# Step2. 按细胞类型跨基因校正；静态路径先去掉 NaN，用非缺失基因数作为 BH 分母
pvalues <- c(0.001, 0.04, 0.04, 0.8, NA, 0.0005, 0.25)
valid_index <- which(!is.na(pvalues))
valid_pvalues <- pvalues[valid_index]
descending <- order(-valid_pvalues, -seq_along(valid_pvalues))
valid_count <- length(valid_pvalues)
adjusted_descending <- pmin(1, cummin(valid_pvalues[descending] * valid_count / seq(valid_count, 1)))
adjusted <- rep(NA_real_, length(pvalues))
adjusted[valid_index[descending]] <- adjusted_descending
adjustment_error <- max(abs(adjusted[valid_index] - p.adjust(valid_pvalues, method="BH")))
checks <- rbind(checks, data.frame(Check="BH_valid_gene_universe", Error=adjustment_error, Passed=adjustment_error < 1e-12))
print(data.frame(pvalues, adjusted))

# Step3. 两项 Fisher 合并：卡方右尾与 df=4 的闭式表达式独立核对
combined_input <- c(1 / 30, 0.002)
fisher_statistic <- -2 * sum(log(combined_input))
combined_probability <- pchisq(fisher_statistic, df=2 * length(combined_input), lower.tail=FALSE)
half_statistic <- fisher_statistic / 2
closed_probability <- exp(-half_statistic) * (1 + half_statistic)
fisher_error <- abs(combined_probability - closed_probability)
checks <- rbind(checks, data.frame(Check="Fisher_two_statistical_inputs", Error=fisher_error, Passed=fisher_error < 1e-12))

# Step4. W=4、N=9 时六个窗口估计；两端分别延长 W-1 行，再将两种对齐结果平均
sample_count <- 9
window_size <- 4
window_expression <- cbind(seq_len(sample_count - window_size + 1), 10 * seq_len(sample_count - window_size + 1))
forward_expression <- rbind(window_expression, window_expression[rep(nrow(window_expression), window_size - 1), , drop=FALSE])
backward_expression <- rbind(window_expression[rep(1, window_size - 1), , drop=FALSE], window_expression)
smoothed_expression <- (forward_expression + backward_expression) / 2
expected_first_type <- c(1, 1.5, 2, 2.5, 3.5, 4.5, 5, 5.5, 6)
smoothing_error <- max(abs(smoothed_expression[, 1] - expected_first_type))
checks <- rbind(checks, data.frame(Check="two_edge_aligned_window_arrays", Error=smoothing_error, Passed=smoothing_error < 1e-12))
print(dim(smoothed_expression))

# Step5. 基因内部按表达排序后，排序索引的逆排列必须恢复输入样本顺序
original_expression <- c(8, 1, 7, 3, 5, 2, 9, 6, 4)
sample_order <- order(original_expression)
sorted_expression <- original_expression[sample_order]
inverse_order <- order(sample_order)
restored_expression <- sorted_expression[inverse_order]
ordering_error <- max(abs(restored_expression - original_expression))
checks <- rbind(checks, data.frame(Check="inverse_sample_permutation", Error=ordering_error, Passed=ordering_error == 0))

# Step6. Worker 下限之后仍可覆盖为基线；写文件前的最终下限需单独保留
calibrated_expression <- c(-0.2, 0.3, 2, 4)
worker_expression <- pmax(calibrated_expression, 1)
worker_expression[c(1, 2)] <- c(0, 0.25) # 模拟被指定为恒定输出的类型，不代表完整筛选规则
written_expression <- pmax(worker_expression, 1)
floor_error <- max(abs(written_expression - c(1, 1, 2, 4)))
checks <- rbind(checks, data.frame(Check="floor_after_constant_override", Error=floor_error, Passed=floor_error == 0))

write.table(checks, file.path(outdir, "independent_hires_output_contract_checks.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(smoothed_expression, file.path(outdir, "toy_hires_edge_alignment.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
capture.output(sessionInfo(), file=file.path(outdir, "hires_output_contract_R_sessionInfo.txt"))
print(checks)
stopifnot(all(checks$Passed))
