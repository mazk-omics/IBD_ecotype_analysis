# 在父类型正式冻结后，把本轮候选G-like转成EcoTyper presorted所需每类型TXT。
# 当前没有冻结父类型清单，因此本轮不执行；只导出输入，不做log、Region3或NMF。
options(stringsAsFactors = FALSE)
root <- '/home/mazekai/IBD_EcoTyper/06_validation/joint20_bayesprism_ecotyper_matched_sensitivity_20261014'
parent_list <- file.path(root, '04_ecotyper_preparation/FROZEN_DISCOVERY_PARENT_TYPES.txt')
stopifnot(file.exists(parent_list))
types <- readLines(parent_list)
reference_types <- readLines(file.path(root, '00_audit/frozen_types.txt'))
stopifnot(length(types) > 0, !anyDuplicated(types), all(types %in% reference_types))
annotation <- read.delim(file.path(root, '00_audit/MSCCR_CD_UC_discovery_manifest.tsv'), check.names = FALSE)
target <- file.path(root, '04_ecotyper_preparation/frozen_parent_inputs')
stopifnot(!dir.exists(target))
dir.create(target)
write.table(annotation[, setdiff(names(annotation), 'frozen_order')], file.path(target, 'annotation.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

# Step1. 按冻结类型逐一导出原始CPM-like单位，样本顺序保持主线不变
for (ct in types) {
  g <- readRDS(file.path(root, '03_expression_qc/candidate_glike_rds', paste0(ct, '.rds')))
  stopifnot(identical(colnames(g), annotation$ID), !anyDuplicated(rownames(g)), all(is.finite(g)), all(g >= 0))
  expression <- data.frame(Gene = rownames(g), g, check.names = FALSE)
  write.table(expression, file.path(target, paste0(ct, '.txt')), sep = '\t', row.names = FALSE, quote = FALSE)
  cat('PRESORTED_TEXT_EXPORTED', ct, dim(g), '\n')
  rm(g, expression)
  gc()
}
writeLines(types, file.path(target, 'parent_types_used.txt'))
cat('INPUT_EXPORT_COMPLETE; EcoTyper NMF not launched\n')
