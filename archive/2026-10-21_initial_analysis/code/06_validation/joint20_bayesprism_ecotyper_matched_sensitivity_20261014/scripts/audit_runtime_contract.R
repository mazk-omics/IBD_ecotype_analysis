# 只读检查已安装BayesPrism及上次120例结果，核对本轮输出提取接口。
# 不重新推断，只验证对象slots、默认参数和一个类型的Z→G-like转换。
suppressPackageStartupMessages(library(BayesPrism))
out <- commandArgs(trailingOnly = TRUE)[1]
benchmark <- '/home/mazekai/IBD_EcoTyper/06_validation/deconv_fraction_sensitivity_joint20_20261013'
stopifnot(as.character(packageVersion('BayesPrism')) == '2.2.3')
for (script in c('run_full_joint20.R', 'prepare_glike_qc.R')) invisible(parse(file.path(out, 'scripts', script)))
bp <- readRDS(file.path(benchmark, '02_bayesprism_joint20_MSCCR120/BayesPrism_result.rds'))
print(slotNames(bp))
print(dim(bp@prism@mixture))
print(dim(bp@posterior.initial.cellType@Z))
z <- get.exp(bp, state.or.type = 'type', cell.name = 'B')
stopifnot(identical(dimnames(z), dimnames(bp@prism@mixture)), all(is.finite(z)), all(z >= 0))
mass <- rowSums(z)
g <- sweep(z[mass > 0, , drop = FALSE], 1, mass[mass > 0], '/') * 1e6
stopifnot(all(is.finite(g)), all(g >= 0), max(abs(rowSums(g) - 1e6)) < 1e-5)
reconstructed <- apply(bp@posterior.initial.cellType@Z, c(1, 2), sum)
residual <- reconstructed - bp@prism@mixture
print(range(residual))
print(range(rowSums(residual)))
print(range(rowSums(bp@prism@mixture - floor(bp@prism@mixture))))
print(get('valid.gibbs.control', asNamespace('BayesPrism'))(list()))
print(get('valid.opt.control', asNamespace('BayesPrism'))(list()))
print(get.exp)
writeLines(capture.output(sessionInfo()), file.path(out, '00_audit/runtime_preflight_sessionInfo.txt'))
cat('RUNTIME_CONTRACT_PASS\n')
