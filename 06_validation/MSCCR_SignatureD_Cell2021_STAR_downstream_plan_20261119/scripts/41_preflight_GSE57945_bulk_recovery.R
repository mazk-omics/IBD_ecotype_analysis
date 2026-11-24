suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 41_preflight_GSE57945_bulk_recovery.R <model_root> <expression_file> <annotation_file> <output_dir>")
}

model_root <- normalizePath(args[1], mustWork = TRUE)
expression_file <- normalizePath(args[2], mustWork = TRUE)
annotation_file <- normalizePath(args[3], mustWork = TRUE)
output_dir <- args[4]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step1. 读取当前模型选定 rank 和 GSE57945 冻结输入。
rank_file <- file.path(model_root, "Analysis", "rank_selection", "rank_data.txt")
rank_data <- fread(rank_file)
expression <- fread(expression_file, check.names = FALSE)
annotation <- fread(annotation_file, check.names = FALSE)

gene_column <- names(expression)[1]
genes <- as.character(expression[[gene_column]])
sample_ids <- names(expression)[-1]

if (anyDuplicated(genes) > 0) {
  stop("Expression gene identifiers are not unique")
}
if (anyDuplicated(sample_ids) > 0) {
  stop("Expression sample identifiers are not unique")
}
if (!"ID" %in% names(annotation)) {
  stop("Annotation does not contain ID")
}
if (!identical(sample_ids, as.character(annotation$ID))) {
  stop("Expression and annotation sample order does not match")
}

expression_matrix <- as.matrix(expression[, -1])
storage.mode(expression_matrix) <- "double"
if (any(!is.finite(expression_matrix))) {
  stop("Expression matrix contains non-finite values")
}
if (any(expression_matrix < 0)) {
  stop("Expression matrix contains negative values")
}

apply_log2 <- all(expression_matrix >= 0) && max(expression_matrix) > 50
if (apply_log2) {
  expression_matrix <- log2(expression_matrix + 1)
}
rownames(expression_matrix) <- genes

# Step2. 按当前 20 个父类型和选定 rank 审计模型基因覆盖与方差。
coverage <- vector("list", nrow(rank_data))
missing_rows <- vector("list", nrow(rank_data))

for (i in seq_len(nrow(rank_data))) {
  cell_type <- rank_data$CellType[i]
  chosen_rank <- rank_data$Chosen_Rank[i]
  gene_info_file <- file.path(
    model_root,
    "Cell_States",
    "discovery",
    cell_type,
    as.character(chosen_rank),
    "initial_gene_info.txt"
  )
  gene_info <- fread(gene_info_file)
  model_genes <- unique(as.character(gene_info$Gene))
  matched <- match(model_genes, genes)
  present <- !is.na(matched)
  gene_sd <- rep(NA_real_, length(model_genes))
  gene_sd[present] <- apply(expression_matrix[matched[present], , drop = FALSE], 1, sd)
  variable <- present & is.finite(gene_sd) & gene_sd > 0

  coverage[[i]] <- data.table(
    CellType = cell_type,
    ChosenRank = chosen_rank,
    ModelGenes = length(model_genes),
    PresentGenes = sum(present),
    MissingGenes = sum(!present),
    VariableGenes = sum(variable),
    ConstantGenes = sum(present & !variable),
    PresentFraction = sum(present) / length(model_genes),
    UsableFraction = sum(variable) / length(model_genes)
  )
  missing_rows[[i]] <- data.table(
    CellType = cell_type,
    ChosenRank = chosen_rank,
    Gene = model_genes[!present]
  )
}

coverage <- rbindlist(coverage)
missing_rows <- rbindlist(missing_rows)

# Step3. 冻结预检查结论和输入身份。
sample_summary <- data.table(
  Metric = c(
    "Genes",
    "Samples",
    "AnnotationRows",
    "CellTypes",
    "InputMaxBeforeEcoTyper",
    "EcoTyperWillLog2",
    "MinimumPresentFraction",
    "MinimumUsableFraction",
    "DecisionThreshold",
    "WarningThreshold"
  ),
  Value = c(
    length(genes),
    length(sample_ids),
    nrow(annotation),
    nrow(rank_data),
    max(as.matrix(expression[, -1])),
    apply_log2,
    min(coverage$PresentFraction),
    min(coverage$UsableFraction),
    0.75,
    0.80
  )
)

input_identity <- data.table(
  File = c(expression_file, annotation_file, rank_file),
  MD5 = unname(tools::md5sum(c(expression_file, annotation_file, rank_file)))
)

status <- "ACCEPTED"
detail <- "Current 20-cell-type model coverage and sample identity passed"
if (min(coverage$PresentFraction) < 0.75 || min(coverage$UsableFraction) < 0.75) {
  status <- "FAILED"
  detail <- "At least one cell type has less than 75% present or variable model genes"
} else if (min(coverage$PresentFraction) < 0.80 || min(coverage$UsableFraction) < 0.80) {
  status <- "ACCEPTED_WITH_WARNING"
  detail <- "Official recovery has no coverage cutoff; CD4 T is below the prespecified 80% warning line but above the 75% stop line"
}

fwrite(coverage, file.path(output_dir, "GSE57945_current85_gene_coverage.tsv"), sep = "\t")
fwrite(missing_rows, file.path(output_dir, "GSE57945_current85_missing_genes.tsv"), sep = "\t")
fwrite(sample_summary, file.path(output_dir, "GSE57945_current85_input_summary.tsv"), sep = "\t")
fwrite(input_identity, file.path(output_dir, "GSE57945_current85_input_md5.tsv"), sep = "\t")
writeLines(c(paste0("status=", status), paste0("detail=", detail)), file.path(output_dir, "preflight_status.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

print(coverage)
cat(paste0("status=", status, "\n"))
if (status == "FAILED") {
  quit(status = 2)
}
