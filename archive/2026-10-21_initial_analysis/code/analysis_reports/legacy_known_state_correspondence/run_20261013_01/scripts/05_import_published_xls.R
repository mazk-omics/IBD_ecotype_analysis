# 用本机已有readxl只读导入公开发表的.xls；原工作簿保持原样。
# 单线程，把各sheet保存为可复用TSV及结构预览。
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
library(readxl)

# Step1. 检查发表表的真实结构
out <- "D:/IBD_ecotype/IBD_ecotype_analysis/analysis_reports/legacy_known_state_correspondence/run_20261013_01"
input <- file.path(out, "literature_resources/41586_2021_3852_MOESM9_ESM.xls")
sheets <- excel_sheets(input)
print(sheets)

# Step2. 导出全部sheet原始单元格，不重建或筛选signature
for (i in seq_along(sheets)) {
  table <- read_excel(input, sheet=sheets[i], col_names=FALSE, .name_repair="minimal")
  target <- file.path(out, paste0("literature_resources/Elmentaite2021_sheet", i, "_raw.tsv"))
  write.table(table, target, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, na="")
  cat("Sheet", i, sheets[i], "rows", nrow(table), "columns", ncol(table), "\n")
  print(head(table, 5))
}
write.table(data.frame(SheetIndex=seq_along(sheets), SheetName=sheets), file.path(out, "audit/Elmentaite2021_SHEET_MAPPING.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
