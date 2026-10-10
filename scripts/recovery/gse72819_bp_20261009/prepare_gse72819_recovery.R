# 在服务器 conda ecotyper 环境运行，核对冻结模型与 GSE72819 输入，生成24例测试和73例完整输入。
# 输入为既有 RPKM 和 discovery 模型；输出写到命令行指定的新目录，不修改源数据或模型。
suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(data.table))
args = commandArgs(TRUE)
run_dir = normalizePath(args[1], mustWork = TRUE)
project_dir = "/home/mazekai/IBD_EcoTyper"
software_dir = file.path(project_dir, "software/ecotyper")
model_dir = file.path(software_dir, "EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes")
source_dir = file.path(project_dir, "06_validation/02_bulk_treatment/GSE72819/prepared/06_ecotyper_recovery_input")

# Step1. 核对冻结输入和样本；缺失的临床响应标签不作为 recovery 排除条件。
source_files = file.path(source_dir, c("GSE72819_bulk_recovery_expression.txt", "GSE72819_bulk_recovery_annotation.txt"))
source_md5 = tools::md5sum(source_files)
stopifnot(identical(unname(source_md5), c("8344cca30ec077c05c6877f151c7e4a4", "f58807ab24a7d3e5350a2f1c4fade877")))
bulk_table = fread(source_files[1], data.table = FALSE, check.names = FALSE)
annotation = read.delim(source_files[2], check.names = FALSE, stringsAsFactors = FALSE)
genes = as.character(bulk_table[[1]])
bulk = as.matrix(bulk_table[, -1])
rownames(bulk) = genes
stopifnot(identical(dim(bulk), c(24656L, 73L)), !anyDuplicated(genes), !anyNA(genes), all(nzchar(genes)))
stopifnot(is.numeric(bulk), all(is.finite(bulk)), all(bulk >= 0), !anyDuplicated(colnames(bulk)))
stopifnot(nrow(annotation) == 73, !anyDuplicated(annotation$ID), setequal(annotation$ID, colnames(bulk)))
annotation = annotation[match(colnames(bulk), annotation$ID), ]
stopifnot(identical(annotation$ID, colnames(bulk)), all(annotation$Disease == "UC"), all(annotation$Timepoint == "Baseline"))
cat("Input:", nrow(bulk), "genes x", ncol(bulk), "samples; RPKM; max", max(bulk), "\n")
print(table(annotation$ResponseWeek10, useNA = "ifany"))
fwrite(data.frame(File = names(source_md5), MD5 = unname(source_md5)), file.path(run_dir, "source_input_md5.tsv"), sep = "\t")

# Step2. 核对12类、53个保留状态和7个生态型；有问题先停，不改变状态或分配阈值。
rank_file = file.path(model_dir, "Analysis/rank_selection/rank_data.txt")
ecotype_file = file.path(model_dir, "Ecotypes/discovery/ecotypes.txt")
ranks = read.delim(rank_file, stringsAsFactors = FALSE)
ecotypes = read.delim(ecotype_file, stringsAsFactors = FALSE)
expected_types = c("B", "CD4_T", "CD8_T", "Endothelial", "Enteroendocrine", "Fibroblast", "Glial", "Goblet", "ILC", "Mast", "Monocyte_Macrophage", "Pericyte")
expected_ranks = c(10, 8, 7, 6, 6, 9, 10, 6, 9, 11, 7, 7)
expected_retained = c(6, 4, 4, 3, 3, 5, 6, 2, 6, 7, 3, 4)
stopifnot(identical(as.character(ranks[[1]]), expected_types), identical(as.numeric(ranks[[2]]), expected_ranks))
stopifnot(nrow(ecotypes) == 53, !anyDuplicated(ecotypes$ID), setequal(ecotypes$Ecotype, paste0("E", 1:7)))
model_files = c(file.path(dirname(model_dir), "config_used.yml"), rank_file, ecotype_file)
coverage_rows = list()
state_support_rows = list()
marker_rows = list()

# 每类的算法输入不同，必要循环只用于逐类核对。覆盖按实际 W 输入计算，不能只看热图 marker。
log_bulk = log2(bulk + 1)
for (i in seq_along(expected_types)) {
  cell_type = expected_types[i]
  state_dir = file.path(model_dir, "Cell_States/discovery", cell_type, expected_ranks[i])
  required_files = file.path(state_dir, c("W.txt", "mapping_to_initial_states.txt", "initial_gene_info.txt", "gene_info.txt", "heatmap_data.txt", "heatmap_top_ann.txt"))
  stopifnot(all(file.exists(required_files)), all(file.info(required_files)$size > 0))
  model_files = c(model_files, required_files[1:4])
  W = as.matrix(read.delim(required_files[1], check.names = FALSE))
  mapping = read.delim(required_files[2], stringsAsFactors = FALSE)
  initial_gene_info = read.delim(required_files[3], stringsAsFactors = FALSE)
  gene_info = read.delim(required_files[4], stringsAsFactors = FALSE)
  stopifnot(all(is.finite(W)), all(W >= 0), ncol(W) == expected_ranks[i], nrow(mapping) == expected_retained[i])
  stopifnot(!anyDuplicated(mapping$State), !anyDuplicated(mapping$InitialState), all(mapping$InitialState %in% colnames(W)))
  stopifnot(setequal(paste0(cell_type, "_", mapping$State), ecotypes$ID[ecotypes$CellType == cell_type]))
  training_genes = sub("__pos$", "", rownames(W)[grepl("__pos$", rownames(W))])
  stopifnot(length(training_genes) > 0, !anyDuplicated(training_genes), all(training_genes %in% initial_gene_info$Gene))
  matched_expression = log_bulk[match(training_genes, rownames(log_bulk)), , drop = FALSE]
  rownames(matched_expression) = training_genes
  matched_expression[is.na(matched_expression)] = 0
  scaled_expression = t(scale(t(matched_expression)))
  scaled_expression[!is.finite(scaled_expression)] = 0
  signed_expression = as.matrix(posneg(scaled_expression))
  # NMF 0.21 的 posneg 不添加后缀；官方预测按位置使用 W，故核对正负数值及模型行顺序。
  expected_signed_names = c(paste0(training_genes, "__pos"), paste0(training_genes, "__neg"))
  expected_signed_values = rbind(pmax(scaled_expression, 0), pmax(-scaled_expression, 0))
  stopifnot(identical(rownames(W), expected_signed_names), isTRUE(all.equal(unname(signed_expression), unname(expected_signed_values), tolerance = 1e-12)))
  variable_signed_rows = apply(signed_expression, 1, var) > 0
  usable_mass = colSums(W[variable_signed_rows, , drop = FALSE]) / colSums(W)
  stopifnot(any(variable_signed_rows), all(is.finite(usable_mass)), all(usable_mass > 0))
  coverage_rows[[i]] = data.frame(CellType = cell_type, InitialRank = ncol(W), RetainedStates = nrow(mapping), TrainingGenes = length(training_genes), PresentGenes = sum(training_genes %in% genes), VariableSignedRows = sum(variable_signed_rows), MinUsableWeightFraction = min(usable_mass))
  state_support_rows[[i]] = data.frame(CellType = cell_type, InitialState = colnames(W), Retained = colnames(W) %in% mapping$InitialState, UsableWeightFraction = usable_mass)
  gene_info = gene_info[order(gene_info$State, -gene_info$MaxFC), ]
  marker_rows[[i]] = data.frame(CellType = cell_type, gene_info, check.names = FALSE)
}
coverage = do.call(rbind, coverage_rows)
print(coverage, row.names = FALSE)
fwrite(coverage, file.path(run_dir, "model_gene_coverage.tsv"), sep = "\t")
fwrite(do.call(rbind, state_support_rows), file.path(run_dir, "initial_state_weight_support.tsv"), sep = "\t")
fwrite(do.call(rbind, marker_rows), file.path(run_dir, "server_discovery_state_genes.tsv"), sep = "\t")
model_md5 = tools::md5sum(model_files)
fwrite(data.frame(File = names(model_md5), MD5 = unname(model_md5)), file.path(run_dir, "frozen_model_md5.tsv"), sep = "\t")
runtime_files = file.path(software_dir, c("EcoTyper_recovery_bulk.R", "pipeline/state_recovery_bulk.R", "pipeline/ecotypes_recovery.R", "pipeline/lib/ecotyper.R", "pipeline/lib/heatmaps.R"))
runtime_md5 = tools::md5sum(runtime_files)
fwrite(data.frame(File = names(runtime_md5), MD5 = unname(runtime_md5)), file.path(run_dir, "runtime_source_md5.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(run_dir, "R_sessionInfo.txt"))

# Step3. 测试均匀选取排序后的24例；保留全73例，输出仍为非log RPKM，不进行临床关联分析。
run_id = basename(run_dir)
sorted_ids = sort(colnames(bulk))
smoke_ids = sorted_ids[unique(round(seq(1, length(sorted_ids), length.out = 24)))]
stopifnot(length(smoke_ids) == 24)
smoke_table = data.frame(GeneSymbol = genes, bulk[, smoke_ids, drop = FALSE], check.names = FALSE)
fwrite(smoke_table, file.path(run_dir, paste0("GSE72819_BPsmoke24_", run_id, ".txt")), sep = "\t", quote = FALSE)
fwrite(annotation[match(smoke_ids, annotation$ID), ], file.path(run_dir, "smoke_annotation.txt"), sep = "\t", quote = FALSE, na = "NA")
full_table = data.frame(GeneSymbol = genes, bulk, check.names = FALSE)
fwrite(full_table, file.path(run_dir, paste0("GSE72819_BP73_", run_id, ".txt")), sep = "\t", quote = FALSE)
fwrite(annotation, file.path(run_dir, "full_annotation.txt"), sep = "\t", quote = FALSE, na = "NA")
fwrite(data.frame(ID = colnames(bulk), Smoke24 = colnames(bulk) %in% smoke_ids), file.path(run_dir, "sample_manifest.tsv"), sep = "\t")
cat("PREFLIGHT_PASSED: 12 cell types, 53 retained states, 7 ecotypes; source and model preserved.\n")
