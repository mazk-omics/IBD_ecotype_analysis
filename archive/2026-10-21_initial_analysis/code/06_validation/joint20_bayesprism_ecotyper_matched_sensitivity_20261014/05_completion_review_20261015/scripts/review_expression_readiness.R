# 复核本轮20类候选表达、posterior信息量及预定marker分布，准备EcoTyper接口审计。
# 服务器bayesprism环境运行；逐类型读取RDS控制内存，不删样本、不筛父类型、不运行NMF。
suppressPackageStartupMessages(library(Matrix))
options(stringsAsFactors = FALSE, warn = 1)
root <- '/home/mazekai/IBD_EcoTyper/06_validation/joint20_bayesprism_ecotyper_matched_sensitivity_20261014'
out <- file.path(root, '05_completion_review_20261015')

# Step1. 读取冻结轴、信息量及原始同源参考
annotation <- read.delim(file.path(root, '00_audit/MSCCR_CD_UC_discovery_manifest.tsv'), check.names = FALSE)
types <- readLines(file.path(root, '00_audit/frozen_types.txt'))
genes <- readLines(file.path(root, '02_bayesprism/retained_genes.txt'))
mass <- read.delim(file.path(root, '03_expression_qc/posterior_mass_by_sample_type.tsv'), check.names = FALSE)
theta_initial <- read.delim(file.path(root, '02_bayesprism/theta_initial.tsv'), row.names = 1, check.names = FALSE)
theta_final <- read.delim(file.path(root, '02_bayesprism/theta_final.tsv'), row.names = 1, check.names = FALSE)
ref <- readRDS(file.path(root, '01_reference/clean_reference.rds'))
stopifnot(identical(rownames(theta_initial), annotation$ID), identical(rownames(theta_final), annotation$ID), setequal(colnames(theta_final), types), length(genes) == length(unique(genes)))
stopifnot(all(is.finite(as.matrix(theta_final))), all(theta_final >= 0), max(abs(rowSums(theta_final) - 1)) < 1e-6)

# Step2. 使用昨天QC计划中的marker，不把共享marker当互斥类型真值
panel <- data.frame(parent = types, genes = c('MS4A1,CD79A', 'MZB1,JCHAIN,SDC1', 'CD3D,CD3E,TRAC', 'CD3D,CD8A,CD8B', 'GNLY,NKG7', 'IL7R,KIT', 'LST1,TYROBP,CTSS', 'CTSS', 'TPSAB1,KIT', 'COL1A1,DCN', 'RGS5,PDGFRB', 'PECAM1,VWF', 'ALPI,FABP1', 'CA1,CA2', 'MUC2,CLCA1', 'OLFM4,MKI67', 'DEFA5,DEFA6', 'CHGA,CHGB', 'S100B,PLP1', 'POU2F3,TRPM5'))
panel$interpretation <- 'descriptive lineage/program markers; not an exclusive identity test'
panel$interpretation[panel$parent %in% c('CD4_T', 'CD8_T', 'NK', 'ILC', 'DC', 'Stem_TA_progenitor')] <- 'shared or mixed-lineage markers; review with reference and low-mass strata'
markers <- unique(unlist(strsplit(panel$genes, ',', fixed = TRUE)))
available_markers <- intersect(markers, genes)
write.table(data.frame(gene = markers, in_retained_genes = markers %in% genes), file.path(out, 'marker_gene_coverage.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(panel, file.path(out, 'prespecified_marker_panels.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
marker_rows <- list()
type_rows <- list()

# Step3. 每类型独立核对矩阵并比较低mass分层；阈值只用于审计，不删样本
for (ct in types) {
  path <- file.path(root, '03_expression_qc/candidate_glike_rds', paste0(ct, '.rds'))
  g <- readRDS(path)
  stopifnot(is.matrix(g), identical(rownames(g), genes), identical(colnames(g), annotation$ID), all(is.finite(g)), all(g >= 0), max(abs(colSums(g) - 1e6)) < 1e-5)
  m <- mass[mass$cell_type == ct, ]
  stopifnot(identical(m$ID, annotation$ID))
  means <- rowMeans(g)
  variances <- rowSums((g - means)^2) / (ncol(g) - 1)
  reference_counts <- Matrix::colSums(ref$reference[ref$labels == ct, genes, drop = FALSE])
  reference_cpm <- reference_counts / sum(reference_counts) * 1e6
  marker_g <- g[available_markers, , drop = FALSE]
  mean_marker <- rowMeans(marker_g)
  mean_marker_high <- rep(NA_real_, length(available_markers))
  if (any(m$posterior_mass >= 1000)) mean_marker_high <- rowMeans(marker_g[, m$posterior_mass >= 1000, drop = FALSE])
  mean_marker_z <- rowMeans(sweep(marker_g, 2, m$posterior_mass, '*') / 1e6)
  marker_rows[[ct]] <- data.frame(cell_type = ct, gene = available_markers, mean_G_like_all = mean_marker, mean_G_like_mass_ge1000 = mean_marker_high, n_mass_ge1000 = sum(m$posterior_mass >= 1000), mean_Z_contribution = mean_marker_z, reference_common_gene_CPM = reference_cpm[available_markers])
  type_rows[[ct]] <- data.frame(cell_type = ct, genes = nrow(g), samples = ncol(g), sample_order_exact = TRUE, genes_unique = TRUE, finite_nonnegative = TRUE, zero_mass_samples = sum(m$zero_mass), mass_min = min(m$posterior_mass), mass_q25 = unname(quantile(m$posterior_mass, 0.25)), mass_median = median(m$posterior_mass), mass_q75 = unname(quantile(m$posterior_mass, 0.75)), mass_max = max(m$posterior_mass), n_mass_lt100 = sum(m$mass_lt100), pct_mass_lt100 = mean(m$mass_lt100) * 100, n_mass_lt1000 = sum(m$mass_lt1000), pct_mass_lt1000 = mean(m$mass_lt1000) * 100, median_initial_fraction = median(theta_initial[, ct]), median_final_fraction = median(theta_final[, ct]), n_final_fraction_lt0_001 = sum(theta_final[, ct] < 0.001), nonzero_genes = sum(means > 0), varying_genes = sum(variances > 1e-12), max_CPM_sum_error = max(abs(colSums(g) - 1e6)), candidate_bytes = file.info(path)$size)
  cat('MATRIX_QC_PASS', ct, dim(g), 'mass_lt100', sum(m$mass_lt100), 'mass_ge1000', sum(m$posterior_mass >= 1000), '\n')
  rm(g, marker_g, variances)
  gc()
}
marker_summary <- do.call(rbind, marker_rows)
type_summary <- do.call(rbind, type_rows)
write.table(marker_summary, file.path(out, 'marker_distribution_by_type.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(type_summary, file.path(out, 'parent_type_information_QC.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

# Step4. 逐面板报告跨类型排序；保留原始read contribution及参考背景帮助识别分母放大
panel_rows <- list()
for (i in seq_len(nrow(panel))) {
  panel_genes <- intersect(strsplit(panel$genes[i], ',', fixed = TRUE)[[1]], available_markers)
  part <- marker_summary[marker_summary$gene %in% panel_genes, ]
  observed <- vapply(types, function(ct) mean(log2(part$mean_G_like_all[part$cell_type == ct] + 1)), numeric(1))
  high_mass <- vapply(types, function(ct) mean(log2(part$mean_G_like_mass_ge1000[part$cell_type == ct] + 1)), numeric(1))
  reference <- vapply(types, function(ct) mean(log2(part$reference_common_gene_CPM[part$cell_type == ct] + 1)), numeric(1))
  contributions <- vapply(types, function(ct) mean(part$mean_Z_contribution[part$cell_type == ct]), numeric(1))
  panel_rows[[i]] <- data.frame(expected_parent = panel$parent[i], observed_parent = types, available_markers = paste(panel_genes, collapse = ','), score_log2_mean_CPM_all = observed, score_log2_mean_CPM_mass_ge1000 = high_mass, reference_score = reference, mean_marker_Z = contributions, observed_rank = rank(-observed, ties.method = 'min'), reference_rank = rank(-reference, ties.method = 'min'), interpretation = panel$interpretation[i])
}
write.table(do.call(rbind, panel_rows), file.path(out, 'cross_parent_marker_panel_review.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
writeLines(c('EXPRESSION_READINESS_QC_COMPLETE', '20 matrices numerically checked; all frozen samples retained.', 'No parent type selected or removed; no NMF or Region3 processing launched.', 'Marker panel ranks are descriptive, shared-marker and low-mass effects require biological interpretation.'), file.path(out, 'expression_review_status.txt'))
writeLines(capture.output(sessionInfo()), file.path(out, 'R_sessionInfo.txt'))
cat('EXPRESSION_READINESS_QC_COMPLETE\n')
