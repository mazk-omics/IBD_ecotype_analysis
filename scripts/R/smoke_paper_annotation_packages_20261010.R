# 用人工小矩阵验证解释环境中的Seurat与fgsea接口能够实际计算。
# 服务器ecotyper_interpretation环境运行，只输出测试结果，不接触项目原始数据。
# 人工结果不能作为细胞状态证据；参数保留论文接口以检查兼容性。
library(Seurat)
library(fgsea)
library(metap)

# Step1. 核对真实加载版本，避免conda登记版本掩盖源码覆盖
out <- "/home/mazekai/IBD_EcoTyper/06_validation/ecotyper_original_method_annotation_20261010/environment"
versions <- data.frame(Package=c("Seurat", "fgsea", "metap"), Version=c(as.character(packageVersion("Seurat")), as.character(packageVersion("fgsea")), as.character(packageVersion("metap"))))
write.table(versions, file.path(out, "smoke_actual_versions.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
stopifnot(versions$Version[versions$Package=="Seurat"]=="3.1.3", versions$Version[versions$Package=="fgsea"]=="1.14.0")

# Step2. 人工两组细胞的Wilcoxon marker测试
set.seed(20261010)
counts <- matrix(rpois(200*20, lambda=3), nrow=200, dimnames=list(paste0("G", seq_len(200)), paste0("Cell", seq_len(20))))
counts[1:10, 1:10] <- counts[1:10, 1:10] + 20
object <- CreateSeuratObject(counts=counts, min.cells=0, min.features=0)
object <- NormalizeData(object, normalization.method="LogNormalize", scale.factor=1000000, verbose=FALSE)
Idents(object) <- rep(c("A", "B"), each=10)
de <- FindMarkers(object, ident.1="A", ident.2="B", test.use="wilcox", min.pct=0.1, logfc.threshold=0.05, verbose=FALSE)
write.table(de, file.path(out, "smoke_Seurat_FindMarkers.tsv"), sep="\t", quote=FALSE, col.names=NA)
stopifnot(nrow(de)>0, "p_val_adj" %in% names(de), all(de$p_val_adj>=0 & de$p_val_adj<=1))
print(head(de, 3))

# Step3. 原文10,000次置换的fgsea接口测试，不使用项目基因集
ranked <- setNames(seq(200, 1), paste0("G", seq_len(200)))
pathways <- list(ArtificialPositive=paste0("G", seq_len(20)))
enrichment <- fgsea(pathways=pathways, stats=ranked, nperm=10000, nproc=1)
enrichment$leadingEdge <- vapply(enrichment$leadingEdge, paste, collapse=";", FUN.VALUE=character(1))
write.table(as.data.frame(enrichment), file.path(out, "smoke_fgsea_10000.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
stopifnot(nrow(enrichment)==1, is.finite(enrichment$NES), enrichment$NES>0)
print(enrichment)
print(metap::sumz(c(0.01, 0.02, 0.03)))
capture.output(sessionInfo(), file=file.path(out, "numerical_smoke_sessionInfo.txt"))
writeLines("NUMERICAL_SMOKE_OK", file.path(out, "numerical_smoke_completed.txt"))
