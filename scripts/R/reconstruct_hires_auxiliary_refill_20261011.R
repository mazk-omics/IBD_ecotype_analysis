# 独立实现HiRes外层辅助表达回填：仅全1行回填，整类型全1时保留原值，然后设表达下限。
# nnls阶段已完成后运行：Rscript 本脚本 收缩结果目录 辅助Group表达表 新输出目录。
# 辅助表是显式外部依赖，可输入官方表或独立Group研究表；来源须由相应审计单独记录。
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_ALL', '.UTF-8')
stopifnot(length(args) == 3, !dir.exists(args[3]))
outdir <- args[3]
dir.create(outdir, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))

# Step1. 辅助表允许NA，但不允许缺基因、错类型或把0误当成缺失
auxiliary <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
stopifnot(!anyDuplicated(rownames(auxiliary)), !anyDuplicated(colnames(auxiliary)))
auxiliary_names <- gsub('-', '.', rownames(auxiliary), fixed = TRUE)
stopifnot(!anyDuplicated(auxiliary_names))  # 外部表可能使用R规范化名称，只接受无歧义映射
records <- vector('list', ncol(auxiliary))
reference_genes <- reference_samples <- NULL

# Step2. 在完整目标集合中检查全1保护，不逐基因块单独判断
for (j in seq_len(ncol(auxiliary))) {
    celltype <- colnames(auxiliary)[j]
    expression <- as.matrix(read.delim(file.path(args[1], paste0('REFINEMENT_ONLY_', celltype, '.tsv')), row.names = 1, check.names = FALSE))
    auxiliary_rows <- match(gsub('-', '.', rownames(expression), fixed = TRUE), auxiliary_names)
    stopifnot(all(is.finite(expression)), all(expression >= 1), !anyNA(auxiliary_rows), !anyDuplicated(auxiliary_rows))
    if (j == 1) {
        reference_genes <- rownames(expression)
        reference_samples <- colnames(expression)
    }
    stopifnot(identical(rownames(expression), reference_genes), identical(colnames(expression), reference_samples))
    row_maximum <- apply(expression, 1, max)
    type_all_one <- all(row_maximum == 1)
    refill_rows <- if (type_all_one) integer(0) else which(row_maximum == 1)
    auxiliary_values <- auxiliary[auxiliary_rows, j]
    if (length(refill_rows)) expression[refill_rows, ] <- matrix(auxiliary_values[refill_rows], length(refill_rows), ncol(expression))
    expression <- pmax(expression, 1)  # NA保持NA，0或低于1的回填值变为1
    records[[j]] <- data.frame(CellType = celltype, genes = nrow(expression), samples = ncol(expression), full_target_type_all_one = type_all_one, refilled_rows = length(refill_rows), auxiliary_NA_refilled = sum(is.na(auxiliary_values[refill_rows])), output_NA_rows = sum(rowSums(is.na(expression)) == ncol(expression)), variable_rows = sum(apply(expression, 1, function(x) if (all(is.na(x))) FALSE else max(x) - min(x) > 1e-8)))
    write.table(expression, file.path(outdir, paste0('AUX_DEPENDENT_', celltype, '.tsv')), sep = '\t', quote = FALSE, col.names = NA, na = 'NA')
}
write.table(do.call(rbind, records), file.path(outdir, 'auxiliary_refill_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(c('Outer numerical refill research implementation; not full CIBERSORTx equivalence.', paste('Explicit auxiliary dependency:', normalizePath(args[2], winslash = '/')), 'Guard is evaluated over the entire supplied target-gene set, not each separately exported block.', 'No official HiRes NA mask is copied. NA derives only from the supplied auxiliary Group values.', 'Auxiliary estimation is external to this script; audit its method and upstream dependencies separately.', 'The worker already imposes its expression floor before this outer stage; this ordering remains under version review.'), file.path(outdir, 'METHOD_LIMITATIONS.txt'))
print(do.call(rbind, records))
