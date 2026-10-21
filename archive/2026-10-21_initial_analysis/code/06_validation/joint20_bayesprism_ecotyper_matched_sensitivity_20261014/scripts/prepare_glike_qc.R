# 从本轮单阶段初始posterior Z导出按细胞类型条件化的CPM-like表达和QC。
# BayesPrism结束后运行，仅保存候选矩阵和标记，不删除样本、不启动NMF。
suppressPackageStartupMessages(library(BayesPrism))
options(stringsAsFactors = FALSE, warn = 1)
out <- commandArgs(trailingOnly = TRUE)[1]

# Step1. 明确Z来源、样本基因轴和完整质量守恒目标
bp <- readRDS(file.path(out, '02_bayesprism/BayesPrism_result.rds'))
annotation <- read.delim(file.path(out, '00_audit/MSCCR_CD_UC_discovery_manifest.tsv'), check.names = FALSE)
types <- readLines(file.path(out, '00_audit/frozen_types.txt'))
mixture <- bp@prism@mixture
theta <- get.fraction(bp, which.theta = 'first', state.or.type = 'type')
stopifnot(identical(rownames(mixture), annotation$ID), identical(rownames(theta), annotation$ID))
reconstructed <- matrix(0, nrow(mixture), ncol(mixture), dimnames = dimnames(mixture))
mass_by_type <- matrix(0, nrow(mixture), length(types), dimnames = list(annotation$ID, types))
qc_records <- list()
dir.create(file.path(out, '03_expression_qc/candidate_glike_rds'), showWarnings = FALSE)

# Step2. 逐类型归一化；零mass保持NA并标记，低mass规则预先定义不依赖NMF
for (ct in types) {
  z <- get.exp(bp, state.or.type = 'type', cell.name = ct)
  stopifnot(identical(dimnames(z), dimnames(mixture)), all(is.finite(z)), all(z >= 0))
  mass <- rowSums(z)
  mass_by_type[, ct] <- mass
  reconstructed <- reconstructed + z
  g <- z
  nonzero <- mass > 0
  g[nonzero, ] <- sweep(z[nonzero, , drop = FALSE], 1, mass[nonzero], '/') * 1e6
  g[!nonzero, ] <- NA_real_
  stopifnot(all(is.finite(g[nonzero, , drop = FALSE])), all(g[nonzero, , drop = FALSE] >= 0), max(abs(rowSums(g[nonzero, , drop = FALSE]) - 1e6), 0) < 1e-5)
  qc <- data.frame(ID = annotation$ID, cell_type = ct, posterior_mass = mass, initial_fraction = theta[, ct], zero_mass = !nonzero, mass_lt100 = mass < 100, mass_lt1000 = mass < 1000, initial_fraction_lt0_001 = theta[, ct] < 0.001, expressed_genes = rowSums(z > 0), normalized_mass = rowSums(g))
  qc_records[[ct]] <- qc
  saveRDS(t(g), file.path(out, '03_expression_qc/candidate_glike_rds', paste0(ct, '.rds')), compress = FALSE)
  cat('CANDIDATE_GLIKE', ct, 'genes x samples', ncol(g), nrow(g), 'zero', sum(!nonzero), 'mass_lt100', sum(mass < 100), '\n')
  rm(z, g)
  gc()
}
write.table(do.call(rbind, qc_records), file.path(out, '03_expression_qc/posterior_mass_by_sample_type.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

# Step3. 对照原始counts与整数化质量，异常只报告并阻止后续release
residual <- reconstructed - mixture
mass_qc <- data.frame(ID = annotation$ID, input_mass = rowSums(mixture), Z_mass = rowSums(reconstructed), residual = rowSums(residual), max_gene_abs_residual = apply(abs(residual), 1, max), fractional_loss_bound = rowSums(mixture - floor(mixture)))
integer_residual <- abs(residual)
integer_residual[mixture != floor(mixture)] <- 0
mass_qc$max_integer_gene_abs_residual <- apply(integer_residual, 1, max)
mass_qc$within_known_integer_sampler_bound <- mass_qc$max_gene_abs_residual <= 1.0001 & abs(mass_qc$residual) <= mass_qc$fractional_loss_bound + 1e-3
write.table(mass_qc, file.path(out, '03_expression_qc/Z_mass_conservation.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
posterior_theta <- mass_by_type / rowSums(mass_by_type)
write.table(data.frame(cell_type = types, max_abs_Z_fraction_minus_initial_theta = apply(abs(posterior_theta - theta[, types]), 2, max)), file.path(out, '03_expression_qc/Z_theta_consistency.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
stopifnot(all(mass_qc$within_known_integer_sampler_bound), all(mass_qc$max_integer_gene_abs_residual <= 1e-6))
writeLines(c('CANDIDATE_GLIKE_QC_COMPLETE', 'No sample or type removed; matrices remain candidates until marker/crosstalk and paired HiRes audit.', 'Zero posterior mass is NA, never imputed or renormalized across types.', 'No EcoTyper discovery was launched.'), file.path(out, 'logs/postprocessing_status.txt'))
cat('CANDIDATE_GLIKE_QC_COMPLETE\n')
