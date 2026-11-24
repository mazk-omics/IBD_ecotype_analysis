suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 42_prepare_GSE57945_bulk_annotation.R <expression_file> <metadata_file> <output_file>")
}

expression_file <- normalizePath(args[1], mustWork = TRUE)
metadata_file <- normalizePath(args[2], mustWork = TRUE)
output_file <- args[3]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Step1. 读取冻结表达样本和规范化 metadata。
expression_header <- fread(expression_file, nrows = 0, check.names = FALSE)
sample_ids <- names(expression_header)[-1]
metadata <- fread(metadata_file, check.names = FALSE)

required_columns <- c(
  "SampleID",
  "GSM",
  "Disease",
  "DiseaseSubtype",
  "Tissue",
  "Histopathology",
  "DeepUlcer",
  "Sex"
)
if (!all(required_columns %in% names(metadata))) {
  stop("Canonical metadata is missing required columns")
}
if (anyDuplicated(metadata$SampleID) > 0) {
  stop("Canonical metadata SampleID is not unique")
}

# Step2. 按表达矩阵顺序生成多列 EcoTyper annotation。
annotation <- metadata[match(sample_ids, SampleID)]
if (any(is.na(annotation$SampleID))) {
  stop("At least one expression sample is absent from canonical metadata")
}

annotation <- annotation[, .(
  ID = SampleID,
  GSM,
  Disease,
  DiseaseSubtype,
  Tissue,
  Histopathology,
  DeepUlcer,
  Sex
)]
if (!identical(annotation$ID, sample_ids)) {
  stop("Prepared annotation order does not match expression columns")
}

# Step3. 保存输入和疾病分母审计。
fwrite(annotation, output_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(
  annotation[, .N, by = .(Disease, DiseaseSubtype, Tissue, Histopathology)],
  paste0(output_file, ".group_counts.tsv"),
  sep = "\t"
)
writeLines(
  c(
    "status=ACCEPTED",
    paste0("samples=", nrow(annotation)),
    paste0("nonIBD=", sum(annotation$Disease == "Non-IBD")),
    paste0("CD=", sum(annotation$Disease == "CD")),
    paste0("UC=", sum(annotation$Disease == "UC"))
  ),
  paste0(output_file, ".status.txt")
)

print(annotation[, .N, by = Disease])
