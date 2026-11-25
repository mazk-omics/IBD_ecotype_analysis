# 按 Cell 2021 STAR Methods 计算 8 个 EcoType 在 pooled scRNA 中的联合检测概率。
# 在服务器运行；读取 61 号脚本冻结的 pooled C，不重新计算 state、EcoType 或单个检验。
# 每次独立置换 C 的每一行，统计成员内均值达到 reference 的 EcoType 数量，共 1,000,000 次。
library(data.table)
library(parallel)

# Step1. 读入 pooled C、冻结成员和 reference 成员内均值
args <- commandArgs(trailingOnly=TRUE)
run_root <- if (length(args)>=1) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
input <- file.path(run_root, "output/48_pooled_scRNA_ecotype_cooccurrence_Cell2021")
out <- file.path(run_root, "output/51_pooled_scRNA_ecotype_joint_probability")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

members <- fread(file.path(run_root, "output/00_freeze/formal_model/ecotype_members_frozen.tsv"), data.table=FALSE)
single_results <- fread(file.path(input, "pooled_four_atlas_ecotype_cooccurrence_8tests.tsv"), data.table=FALSE)
C <- as.matrix(read.delim(file.path(input, "pooled_four_atlas_C_four_matrix_mean.tsv"), row.names=1, check.names=FALSE))
ecotypes <- paste0("E", 1:8)
stopifnot(nrow(C)==79, ncol(C)==79, setequal(rownames(C), members$ID), setequal(single_results$Ecotype, ecotypes), all(is.na(diag(C))))

members <- members[match(rownames(C), members$ID), ]
membership <- lapply(ecotypes, function(ecotype) which(members$Ecotype==ecotype))
names(membership) <- ecotypes
reference <- single_results$MeanWithin[match(ecotypes, single_results$Ecotype)]
stopifnot(all(is.finite(reference)), all(lengths(membership)>=3))

# Step2. 分 32 个可复现随机流执行一百万次逐行置换
total_permutations <- 1000000L
workers <- 32L
seed <- 20261125L
iterations_by_worker <- rep(total_permutations%/%workers, workers)
iterations_by_worker[seq_len(total_permutations%%workers)] <- iterations_by_worker[seq_len(total_permutations%%workers)]+1L
worker_seeds <- seed+seq_len(workers)*100003L

worker_counts <- mclapply(seq_len(workers), function(worker_index) {
    set.seed(worker_seeds[worker_index])
    counts <- integer(length(ecotypes)+1L)
    iterations <- iterations_by_worker[worker_index]

    for (iteration in seq_len(iterations)) {
        shuffled <- C
        for (row_index in seq_len(nrow(C))) {
            order_index <- sample.int(ncol(C))
            shuffled[row_index, ] <- C[row_index, order_index]
            original_diagonal_position <- which(order_index==row_index)
            diagonal_value <- shuffled[row_index, row_index]
            shuffled[row_index, row_index] <- shuffled[row_index, original_diagonal_position]
            shuffled[row_index, original_diagonal_position] <- diagonal_value
        }

        detected <- logical(length(ecotypes))
        for (ecotype_index in seq_along(ecotypes)) {
            member_index <- membership[[ecotype_index]]
            within <- rowMeans(shuffled[member_index, member_index, drop=FALSE], na.rm=TRUE)
            shuffled_mean <- mean(within[is.finite(within)])
            detected[ecotype_index] <- is.finite(shuffled_mean) && shuffled_mean>=reference[ecotype_index]
        }
        counts[sum(detected)+1L] <- counts[sum(detected)+1L]+1L
    }
    counts
}, mc.cores=workers, mc.preschedule=TRUE)

counts <- Reduce("+", worker_counts)
stopifnot(sum(counts)==total_permutations)
joint_exceedances <- counts[length(ecotypes)+1L]
joint_empirical_p <- (joint_exceedances+1)/(total_permutations+1)

# Step3. 保存联合分布、运行参数和结论边界
distribution <- data.frame(DetectedEcoTypes=0:length(ecotypes), Permutations=counts, Fraction=counts/total_permutations)
summary <- data.frame(EcoTypes=length(ecotypes), TotalPermutations=total_permutations, JointExceedances=joint_exceedances, JointEmpiricalP=joint_empirical_p, MinimumResolution=1/(total_permutations+1), Workers=workers, Seed=seed)
worker_manifest <- data.frame(Worker=seq_len(workers), Iterations=iterations_by_worker, Seed=worker_seeds)
fwrite(distribution, file.path(out, "joint_detected_ecotype_count_distribution.tsv"), sep="\t")
fwrite(summary, file.path(out, "joint_probability_summary.tsv"), sep="\t")
fwrite(worker_manifest, file.path(out, "joint_probability_worker_manifest.tsv"), sep="\t")
print(distribution)
print(summary)

capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
writeLines(c(
    "Global test: each row of the signed pooled C matrix is independently permuted; the original diagonal NA is moved back to the diagonal after each shuffle.",
    "For each EcoType, the statistic is the mean member-within co-occurrence, matching the Cell 2021 global-test description rather than the individual AvgDif statistic.",
    "Each permutation records how many of the 8 shuffled EcoTypes equal or exceed their own observed member-within mean; the joint P is the probability that all 8 do so.",
    "One million permutations are split across 32 deterministic worker seeds. The plus-one empirical estimator has minimum resolution 1/1,000,001.",
    "This global P does not change which individual EcoTypes pass the separate 1,000-permutation test."
), file.path(out, "method_boundary.txt"))
