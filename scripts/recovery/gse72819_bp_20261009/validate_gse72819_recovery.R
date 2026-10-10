# 在 conda ecotyper 环境核对官方 recovery 输出，检查全样本分母和冻结模型是否保持。
# 参数：本轮目录、smoke/full；只保存运行审计，不进行响应或预后关联。
suppressPackageStartupMessages(library(data.table))
args = commandArgs(TRUE)
run_dir = normalizePath(args[1])
stage = args[2]
stopifnot(stage %in% c("smoke", "full"))
software_dir = "/home/mazekai/IBD_EcoTyper/software/ecotyper"
model_dir = file.path(software_dir, "EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes")

# Step1. 用实际输入的样本顺序检查输出，不把未分配样本从分母删除。
prefix = if (stage == "smoke") "GSE72819_BPsmoke24_" else "GSE72819_BP73_"
dataset = paste0(prefix, basename(run_dir))
sample_ids = names(fread(file.path(run_dir, paste0(dataset, ".txt")), nrows = 0))[-1]
expected_n = if (stage == "smoke") 24L else 73L
stopifnot(length(sample_ids) == expected_n, !anyDuplicated(sample_ids))
result_dir = file.path(run_dir, paste0(stage, "_results"), dataset)
ranks = read.delim(file.path(model_dir, "Analysis/rank_selection/rank_data.txt"), stringsAsFactors = FALSE)
audit_rows = list()

# Step2. 保留状态比例可能不满1：被过滤初始状态的质量不能再次归一化隐藏。
for (i in seq_len(nrow(ranks))) {
  cell_type = ranks[i, 1]
  cell_dir = file.path(result_dir, cell_type)
  mapping = read.delim(file.path(model_dir, "Cell_States/discovery", cell_type, ranks[i, 2], "mapping_to_initial_states.txt"), stringsAsFactors = FALSE)
  H = as.matrix(read.delim(file.path(cell_dir, "state_abundances.txt"), check.names = FALSE))
  assignment = read.delim(file.path(cell_dir, "state_assignment.txt"), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(H), sample_ids), identical(rownames(H), mapping$State), all(is.finite(H)), all(H >= 0))
  stopifnot(all(colSums(H) <= 1 + 1e-8), !anyDuplicated(assignment$ID), all(assignment$ID %in% sample_ids), all(assignment$State %in% mapping$State))
  audit_rows[[i]] = data.frame(CellType = cell_type, Samples = ncol(H), RetainedStates = nrow(H), AssignedSamples = nrow(assignment), UnassignedSamples = expected_n - nrow(assignment), MedianRetainedMass = median(colSums(H)), MinRetainedMass = min(colSums(H)))
}
audit = do.call(rbind, audit_rows)
print(audit, row.names = FALSE)
fwrite(audit, file.path(run_dir, paste0(stage, "_state_recovery_audit.tsv")), sep = "\t")

# Step3. 生态型核对包含 Unassigned 的全样本记录，保持官方 Q<0.25 规则。
internal_ecotype_dir = file.path(model_dir, "Ecotypes/recovery", dataset)
initial_assignment = read.delim(file.path(internal_ecotype_dir, "initial_ecotype_assignment.txt"), stringsAsFactors = FALSE)
eco_H = as.matrix(read.delim(file.path(result_dir, "Ecotypes/ecotype_abundance.txt"), check.names = FALSE))
stopifnot(nrow(initial_assignment) == expected_n, !anyDuplicated(initial_assignment$ID), setequal(initial_assignment$ID, sample_ids))
stopifnot(setequal(rownames(eco_H), paste0("E", 1:7)), setequal(colnames(eco_H), sample_ids), all(is.finite(eco_H)), all(eco_H >= 0), all(abs(colSums(eco_H) - 1) < 1e-8))
stopifnot(all(initial_assignment$Ecotype %in% c(paste0("E", 1:7), "Unassigned")))
fwrite(initial_assignment, file.path(run_dir, paste0(stage, "_all_sample_ecotype_assignment.tsv")), sep = "\t", na = "NA")
counts = as.data.frame(table(factor(initial_assignment$Ecotype, levels = c(paste0("E", 1:7), "Unassigned"))))
names(counts) = c("Ecotype", "Samples")
print(counts, row.names = FALSE)
fwrite(counts, file.path(run_dir, paste0(stage, "_ecotype_counts.tsv")), sep = "\t")

# Step4. 验证计算没有改变 discovery 模型；运行结果通过不等于生物学验证通过。
frozen_hashes = fread(file.path(run_dir, "frozen_model_md5.tsv"), data.table = FALSE)
stopifnot(identical(unname(tools::md5sum(frozen_hashes$File)), frozen_hashes$MD5))
writeLines(paste(toupper(stage), "PASSED", sep = "_"), file.path(run_dir, paste0(toupper(stage), "_PASSED")))
cat(toupper(stage), "output validation passed;", expected_n, "samples; model unchanged.\n")
