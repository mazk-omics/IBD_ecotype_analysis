# 本机复核轻量ORA的P和BH，不读取表达矩阵、不改变既有环境。
# 独立目录运行，数值检查仅使用本次2650行结果。
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

# Step1. 读入本次ORA结果
out <- "D:/IBD_ecotype/IBD_ecotype_analysis/analysis_reports/legacy_bayesprism_annotation/run_20261013_01"
ora <- read.delim(file.path(out, "CELL_STATE_FUNCTIONAL_ENRICHMENT.tsv"), check.names=FALSE)

# Step2. 用独立R实现核对上尾和全局BH
p_check <- phyper(ora$OverlapN-1, ora$SetInBackgroundN, ora$BackgroundN-ora$SetInBackgroundN, ora$SelectedN, lower.tail=FALSE)
q_check <- p.adjust(p_check, method="BH")
qc <- data.frame(CellState=ora$CellState, GeneSet=ora$GeneSet, PythonP=ora$P, R_P=p_check, AbsolutePError=abs(ora$P-p_check), PythonGlobalQ=ora$FDR_Global2650, R_GlobalQ=q_check, AbsoluteQError=abs(ora$FDR_Global2650-q_check))
write.table(qc, file.path(out, "audit/ORA_INDEPENDENT_R_QC.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
print(c(Tests=nrow(qc), MaxPError=max(qc$AbsolutePError), MaxQError=max(qc$AbsoluteQError)))
stopifnot(max(qc$AbsolutePError)<1e-10, max(qc$AbsoluteQError)<1e-10)
