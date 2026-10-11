# 核对逐列置换和整行置换的科学含义；使用固定120例比例，未运行官方随机流。
# 本机 R 运行：Rscript 本脚本 输入目录 输出目录；只保存性质检查，不输出新的研究比例。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
fractions <- as.matrix(read.delim(file.path(args[1], 'fractions_fixed.tsv'), row.names = 1, check.names = FALSE))

# Step1. 旧代码每列分别无放回重排，保留列边际但破坏样本内组成
set.seed(20261011)
column_shuffled <- fractions
for (j in seq_len(ncol(fractions))) column_shuffled[, j] <- sample(fractions[, j], nrow(fractions), replace = FALSE)
column_marginal_error <- max(abs(apply(fractions, 2, sort) - apply(column_shuffled, 2, sort)))
column_mean_error <- max(abs(colMeans(fractions) - colMeans(column_shuffled)))
column_shuffle_row_sum_change <- max(abs(rowSums(column_shuffled) - 1))
row_shuffled <- fractions[sample.int(nrow(fractions)), , drop = FALSE]
row_shuffle_row_sum_change <- max(abs(rowSums(row_shuffled) - 1))
stopifnot(column_marginal_error == 0, column_mean_error < 1e-12, column_shuffle_row_sum_change > 1e-3, row_shuffle_row_sum_change < 1e-8)
checks <- data.frame(check = c('column_values_preserved', 'column_means_preserved', 'column_shuffle_breaks_row_composition', 'whole_row_shuffle_keeps_composition'), evidence = c(column_marginal_error, column_mean_error, column_shuffle_row_sum_change, row_shuffle_row_sum_change), passed = TRUE)
write.table(checks, file.path(args[2], 'column_randomization_contract.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
print(checks)
