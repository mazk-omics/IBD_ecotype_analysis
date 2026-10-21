# 使用与 CIBERSORTx 相同训练细胞及20类型，重新推断 MSCCR120 与自然模拟28例。
# 在服务器 bayesprism 环境运行；先做三样本 smoke test，通过后完整推断，保留初始/更新 theta 和 CV。
suppressPackageStartupMessages(library(BayesPrism))
suppressPackageStartupMessages(library(Matrix))
options(stringsAsFactors = FALSE)
out <- commandArgs(trailingOnly = TRUE)[1]
started <- Sys.time()

# Step1. 读取细胞原始counts，核对20类一一对应
types <- readLines(file.path(out, '00_audit/frozen_types.txt'))
cells <- read.delim(file.path(out, '01_reference/reference_cells.tsv'), check.names = FALSE)
genes <- readLines(file.path(out, '01_reference/reference_genes.tsv'))
reference_raw <- as(readMM(file.path(out, '01_reference/reference_counts_cells_by_genes.mtx')), 'dgCMatrix')
rownames(reference_raw) <- cells$cell_id
colnames(reference_raw) <- genes
labels <- cells$coarse_identity
stopifnot(nrow(reference_raw) == 9294, length(types) == 20, setequal(unique(labels), types), !anyDuplicated(genes), all(rowSums(reference_raw) > 0))
cat('BayesPrism version:', as.character(packageVersion('BayesPrism')), '\n')
print(dim(reference_raw))
print(table(labels))
writeLines(capture.output(sessionInfo()), file.path(out, 'logs/R_sessionInfo.txt'))

# Step2. 沿用正式BayesPrism清理参数，不按signature4474基因限制输入
groups_remove <- c('Rb', 'Mrp', 'other_Rb', 'chrM', 'MALAT1', 'chrX', 'chrY')
reference_clean <- cleanup.genes(input = reference_raw, input.type = 'count.matrix', species = 'hs', gene.group = groups_remove, exp.cells = 1)
stopifnot(all(rowSums(reference_clean) > 0), setequal(unique(labels), types))
write.table(data.frame(gene = genes, retained_after_cleanup = genes %in% colnames(reference_clean)), file.path(out, '01_reference/gene_cleanup_audit.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
saveRDS(list(reference = reference_clean, labels = labels, types = types), file.path(out, '01_reference/clean_reference.rds'))
rm(reference_raw)
gc()

# Step3. 完整120例输入先冻结共同基因和outlier过滤，smoke使用相同基因
bulk <- read.delim(gzfile(file.path(out, '02_bayesprism_joint20_MSCCR120/mixture_raw_counts.tsv.gz')), row.names = 1, check.names = FALSE)
mixture <- t(as.matrix(bulk))
annotation <- read.delim(file.path(out, '00_audit/MSCCR120_annotation.tsv'), check.names = FALSE)
stopifnot(nrow(mixture) == 120, identical(rownames(mixture), annotation$ID), all(is.finite(mixture)), all(mixture >= 0), all(rowSums(mixture) > 0), !anyDuplicated(colnames(mixture)))
prism <- new.prism(reference = reference_clean, input.type = 'count.matrix', cell.type.labels = labels, cell.state.labels = labels, key = NULL, mixture = mixture, outlier.cut = 0.01, outlier.fraction = 0.1, pseudo.min = 1e-8)
stopifnot(setequal(rownames(prism@phi_cellType@phi), types))
saveRDS(prism, file.path(out, '02_bayesprism_joint20_MSCCR120/prism_input.rds'))
write.table(data.frame(gene = colnames(mixture), retained_after_cleanup = colnames(mixture) %in% colnames(reference_clean), retained_in_prism = colnames(mixture) %in% colnames(prism@mixture)), file.path(out, '02_bayesprism_joint20_MSCCR120/input_gene_filter_audit.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
writeLines(colnames(prism@mixture), file.path(out, '02_bayesprism_joint20_MSCCR120/input_genes.txt'))
cat('Full prism dimensions:', dim(prism@mixture), '\n')
smoke_indices <- match(vapply(c('Colon', 'Rectum', 'Ileum'), function(region) annotation$ID[which(annotation$Region3 == region)[1]], character(1)), rownames(prism@mixture))
smoke_prism <- prism
smoke_prism@mixture <- prism@mixture[smoke_indices, , drop = FALSE]
set.seed(20261013)
smoke <- run.prism(prism = smoke_prism, n.cores = 3, update.gibbs = TRUE)
smoke_final <- get.fraction(smoke, which.theta = 'final', state.or.type = 'type')
stopifnot(nrow(smoke_final) == 3, setequal(colnames(smoke_final), types), all(is.finite(smoke_final)), all(smoke_final >= 0), max(abs(rowSums(smoke_final) - 1)) < 1e-6)
write.table(smoke_final[, types, drop = FALSE], file.path(out, '02_bayesprism_joint20_MSCCR120/smoke_final_theta.tsv'), sep = '\t', quote = FALSE, col.names = NA)
writeLines('COMPLETE: three regions; full 20 types; valid final theta', file.path(out, 'logs/smoke_status.txt'))
cat('SMOKE_COMPLETE\n')
rm(smoke, smoke_prism, bulk)
gc()

# Step4. 完整默认Gibbs计算，不修改迭代数，8核避免占满共享服务器
ncores <- 8
set.seed(20261013)
formal_start <- Sys.time()
bp <- run.prism(prism = prism, n.cores = ncores, update.gibbs = TRUE)
saveRDS(bp, file.path(out, '02_bayesprism_joint20_MSCCR120/BayesPrism_result.rds'))
initial <- get.fraction(bp, which.theta = 'first', state.or.type = 'type')
final <- get.fraction(bp, which.theta = 'final', state.or.type = 'type')
stopifnot(setequal(colnames(final), types), identical(rownames(final), annotation$ID), all(is.finite(final)), all(final >= 0), max(abs(rowSums(final) - 1)) < 1e-6)
write.table(initial[, types, drop = FALSE], file.path(out, '02_bayesprism_joint20_MSCCR120/theta_initial.tsv'), sep = '\t', quote = FALSE, col.names = NA)
write.table(final[, types, drop = FALSE], file.path(out, '02_bayesprism_joint20_MSCCR120/theta_final.tsv'), sep = '\t', quote = FALSE, col.names = NA)
# 初始state与type一一对应，因此使用同名state CV；更新CV来自thetaf后验。
write.table(bp@posterior.initial.cellState@theta.cv[, types, drop = FALSE], file.path(out, '02_bayesprism_joint20_MSCCR120/theta_initial_cv.tsv'), sep = '\t', quote = FALSE, col.names = NA)
write.table(bp@posterior.theta_f@theta.cv[, types, drop = FALSE], file.path(out, '02_bayesprism_joint20_MSCCR120/theta_final_cv.tsv'), sep = '\t', quote = FALSE, col.names = NA)
writeLines(c(paste('formal_elapsed_seconds', as.numeric(difftime(Sys.time(), formal_start, units = 'secs'))), paste('n.cores', ncores), 'update.gibbs TRUE', 'gibbs.control installed defaults', 'opt.control installed defaults', 'seed 20261013', 'input.type count.matrix', 'key NULL', 'outlier.cut 0.01', 'outlier.fraction 0.1', 'pseudo.min 1e-8', paste('cleanup_groups', paste(groups_remove, collapse = ','))), file.path(out, '02_bayesprism_joint20_MSCCR120/run_parameters.txt'))
cat('MSCCR120_COMPLETE\n')
rm(bp, prism, mixture)
gc()

# Step5. 精确恢复的自然counts模拟；另存本轮保留基因上的UMI真值
sim <- read.delim(gzfile(file.path(out, '05_existing_pseudobulk_benchmark/natural28_raw_counts.tsv.gz')), row.names = 1, check.names = FALSE)
sim_mixture <- t(as.matrix(sim))
sim_prism <- new.prism(reference = reference_clean, input.type = 'count.matrix', cell.type.labels = labels, cell.state.labels = labels, key = NULL, mixture = sim_mixture, outlier.cut = 0.01, outlier.fraction = 0.1, pseudo.min = 1e-8)
writeLines(colnames(sim_prism@mixture), file.path(out, '05_existing_pseudobulk_benchmark/BayesPrism_retained_genes.txt'))
saveRDS(sim_prism, file.path(out, '05_existing_pseudobulk_benchmark/prism_input.rds'))
set.seed(20261013)
sim_bp <- run.prism(prism = sim_prism, n.cores = ncores, update.gibbs = TRUE)
saveRDS(sim_bp, file.path(out, '05_existing_pseudobulk_benchmark/BayesPrism_natural28_result.rds'))
sim_initial <- get.fraction(sim_bp, which.theta = 'first', state.or.type = 'type')
sim_final <- get.fraction(sim_bp, which.theta = 'final', state.or.type = 'type')
stopifnot(nrow(sim_final) == 28, setequal(colnames(sim_final), types), all(is.finite(sim_final)), all(sim_final >= 0), max(abs(rowSums(sim_final) - 1)) < 1e-6)
write.table(sim_initial[, types, drop = FALSE], file.path(out, '05_existing_pseudobulk_benchmark/theta_initial.tsv'), sep = '\t', quote = FALSE, col.names = NA)
write.table(sim_final[, types, drop = FALSE], file.path(out, '05_existing_pseudobulk_benchmark/theta_final.tsv'), sep = '\t', quote = FALSE, col.names = NA)
write.table(sim_bp@posterior.theta_f@theta.cv[, types, drop = FALSE], file.path(out, '05_existing_pseudobulk_benchmark/theta_final_cv.tsv'), sep = '\t', quote = FALSE, col.names = NA)
writeLines(paste('total_elapsed_seconds', as.numeric(difftime(Sys.time(), started, units = 'secs'))), file.path(out, 'logs/computation_complete.txt'))
cat('NATURAL28_COMPLETE\n')
