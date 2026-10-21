# 在服务器 ecotyper 环境读取已有 log2 RMA，准备本轮独立 Recovery 输入。
# 不逆变换或更改基因映射，Recovery标准化仅用67份W0，避免随访表达泄漏。
suppressPackageStartupMessages(library(data.table))
root = "/home/mazekai/IBD_EcoTyper"
out = file.path(root, "06_validation/external_clinical_outcomes_20261015")

# Step1. 核对 RDS 的对象结构、样本及数值范围。
input = file.path(root, "06_validation/02_bulk_treatment/GSE73661/prepared/03_canonical_expression/GSE73661_canonical_log2_RMA_18407x178.rds")
x = readRDS(input)
print(class(x))
if (is.list(x)) {
  print(names(x))
  x = x$expression
}
stopifnot(is.matrix(x), identical(dim(x), c(18407L, 178L)), !anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)), all(is.finite(x)))
print(dim(x))
print(range(x))
stopifnot(max(x) <= 50)

# Step2. 先排除所有治疗后样本及对照，再进行官方Recovery的跨样本基因标准化。
metadata = fread(file.path(root, "06_validation/02_bulk_treatment/GSE73661/prepared/04_canonical_metadata/GSE73661_canonical_metadata_178_samples.tsv"), data.table = FALSE)
baseline = metadata[metadata$Timepoint == "W0", ]
stopifnot(nrow(baseline) == 67, !anyDuplicated(baseline$PatientKey), all(baseline$GSM %in% colnames(x)))
x = x[, baseline$GSM, drop = FALSE]
stopifnot(ncol(x) == 67, all(is.finite(x)))
fwrite(data.frame(GENES = rownames(x), x, check.names = FALSE), file.path(out, "03_GSE73661_mucosal_healing/GSE73661_RMA_W0_67_outcomes_20261015.txt"), sep = "\t")
fwrite(data.frame(ID = colnames(x), Cohort = "GSE73661_W0_only"), file.path(out, "03_GSE73661_mucosal_healing/recovery_annotation.txt"), sep = "\t")
fwrite(data.frame(Accession = "GSE73661", Genes = nrow(x), Samples = ncol(x), OriginalDepositedSamples = 178, Min = min(x), Max = max(x), NegativeN = sum(x < 0), Scale = "GEO log2 RMA; as deposited", Transformation = "None; Recovery gene standardization on W0-only 67 biopsies", MD5 = unname(tools::md5sum(input))), file.path(out, "03_GSE73661_mucosal_healing/expression_QC.tsv"), sep = "\t")
