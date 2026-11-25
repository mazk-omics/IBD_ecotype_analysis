# 按 Cell 2021 STAR Methods 合并四个单细胞 atlas，正式验证冻结的 8 个 EcoType 共现。
# 在服务器运行；读取已有 recovery 分配和原始样本注释，不重新拟合 state 或 EcoType。
# 原文展示口径为 1000 次逐行置换得到 Z 后转单侧 P；同时保存经验 P 和 BH 作为审计。
library(data.table)

# Step1. 定位冻结模型、四个 recovery 任务表和输出
args <- commandArgs(trailingOnly=TRUE)
run_root <- if (length(args)>=1) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
out <- file.path(run_root, "output/48_pooled_scRNA_ecotype_cooccurrence_Cell2021")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

members <- fread(file.path(run_root, "output/00_freeze/formal_model/ecotype_members_frozen.tsv"), data.table=FALSE)
ecotypes <- paste0("E", 1:8)
stopifnot(nrow(members)==79, setequal(unique(members$Ecotype), ecotypes), !anyDuplicated(members$ID))

jobs_paths <- c(
    GSE214695=file.path(run_root, "output/04c_GSE214695_official_recovery/GSE214695_official_recovery_jobs.tsv"),
    SCP1884=file.path(run_root, "output/05c_SCP1884_official_recovery/SCP1884_official_recovery_jobs.tsv"),
    SCP259=file.path(run_root, "output/06c_SCP259_official_recovery/SCP259_official_recovery_jobs.tsv"),
    IBDverse_E_MTAB_16999=file.path(run_root, "output/07c_IBDverse_official_recovery/IBDverse_E_MTAB_16999_official_recovery_jobs.tsv")
)
stopifnot(all(file.exists(jobs_paths)))

# Step2. 合并四个 atlas，给样本和细胞加数据集前缀防止重名
assigned_atlases <- list()
annotation_atlases <- list()
input_audit <- list()

for (dataset in names(jobs_paths)) {
    jobs <- fread(jobs_paths[[dataset]], data.table=FALSE)
    stopifnot(!anyDuplicated(jobs$CellType), all(file.exists(jobs$AnnotationFile)))
    assigned_tables <- list()
    annotation_tables <- list()

    for (job_index in seq_len(nrow(jobs))) {
        assignment_path <- file.path(jobs$ResultDir[job_index], "state_assignment.txt")
        stopifnot(file.exists(assignment_path))
        assignment <- fread(assignment_path, data.table=FALSE)
        annotation <- fread(jobs$AnnotationFile[job_index], data.table=FALSE)
        stopifnot(all(c("ID", "State") %in% names(assignment)), all(c("ID", "BiologicalSampleID") %in% names(annotation)))
        stopifnot(!anyDuplicated(assignment$ID), !anyDuplicated(annotation$ID), all(assignment$ID %in% annotation$ID))

        assignment$CellType <- jobs$CellType[job_index]
        assignment$CellState <- paste(assignment$CellType, assignment$State, sep="_")
        assignment$BiologicalSampleID <- annotation$BiologicalSampleID[match(assignment$ID, annotation$ID)]
        annotation$CellType <- jobs$CellType[job_index]
        assigned_tables[[job_index]] <- assignment[, c("ID", "CellType", "CellState", "BiologicalSampleID")]
        annotation_tables[[job_index]] <- annotation[, c("ID", "CellType", "BiologicalSampleID")]
    }

    assigned <- rbindlist(assigned_tables, use.names=TRUE)
    annotation <- rbindlist(annotation_tables, use.names=TRUE)
    stopifnot(!anyDuplicated(assigned$ID), !anyDuplicated(annotation$ID), all(assigned$ID %in% annotation$ID))
    assigned$ID <- paste(dataset, assigned$ID, sep="::")
    annotation$ID <- paste(dataset, annotation$ID, sep="::")
    assigned$BiologicalSampleID <- paste(dataset, assigned$BiologicalSampleID, sep="::")
    annotation$BiologicalSampleID <- paste(dataset, annotation$BiologicalSampleID, sep="::")
    assigned$Dataset <- dataset
    annotation$Dataset <- dataset
    assigned_atlases[[dataset]] <- assigned
    annotation_atlases[[dataset]] <- annotation
    input_audit[[dataset]] <- data.frame(Dataset=dataset, ParentTypes=nrow(jobs), BiologicalSamples=length(unique(annotation$BiologicalSampleID)), InputCells=nrow(annotation), AssignedCells=nrow(assigned))
}

assigned <- rbindlist(assigned_atlases, use.names=TRUE)
annotation <- rbindlist(annotation_atlases, use.names=TRUE)
stopifnot(!anyDuplicated(assigned$ID), !anyDuplicated(annotation$ID), all(assigned$ID %in% annotation$ID))
fwrite(rbindlist(input_audit), file.path(out, "pooled_four_atlas_input_audit.tsv"), sep="\t")
print(rbindlist(input_audit))

# Step3. 构建 pooled F1 和 F2；F1 分母包含父类型内全部正式分配状态
sample_names <- unique(annotation$BiologicalSampleID)
cell_types <- unique(members$CellType)
numerator <- table(factor(assigned$CellState, levels=members$ID), factor(assigned$BiologicalSampleID, levels=sample_names))
all_denominator <- table(factor(annotation$CellType, levels=cell_types), factor(annotation$BiologicalSampleID, levels=sample_names))
assigned_denominator <- table(factor(assigned$CellType, levels=cell_types), factor(assigned$BiologicalSampleID, levels=sample_names))
D1 <- assigned_denominator[members$CellType, , drop=FALSE]
D2 <- all_denominator[members$CellType, , drop=FALSE]
F1 <- numerator/D1
F2 <- numerator/D2
F1[D1==0] <- NA_real_
F2[D2==0] <- NA_real_
rownames(F1) <- rownames(F2) <- members$ID
stopifnot(all(F1[is.finite(F1)]>=0 & F1[is.finite(F1)]<=1), all(F2[is.finite(F2)]>=0 & F2[is.finite(F2)]<=1))

write.table(F1, file.path(out, "pooled_four_atlas_F1_assigned_denominator.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(F2, file.path(out, "pooled_four_atlas_F2_all_cells_denominator.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(D2, file.path(out, "pooled_four_atlas_parent_type_input_counts.tsv"), sep="\t", quote=FALSE, col.names=NA)

# Step4. 计算四套 Spearman 矩阵并严格平均为 C
F1_nozero <- F1
F2_nozero <- F2
F1_nozero[F1_nozero==0 & !is.na(F1_nozero)] <- NA_real_
F2_nozero[F2_nozero==0 & !is.na(F2_nozero)] <- NA_real_
fractions <- list(F1=F1, F2=F2, F1_nozero=F1_nozero, F2_nozero=F2_nozero)
correlations <- list()

for (version in names(fractions)) {
    correlations[[version]] <- suppressWarnings(cor(t(fractions[[version]]), method="spearman", use="pairwise.complete.obs"))
    pair_n <- tcrossprod(1L*!is.na(fractions[[version]]))
    write.table(correlations[[version]], file.path(out, paste0("pooled_four_atlas_C_", version, ".tsv")), sep="\t", quote=FALSE, col.names=NA)
    write.table(pair_n, file.path(out, paste0("pooled_four_atlas_pair_N_", version, ".tsv")), sep="\t", quote=FALSE, col.names=NA)
}

C_display_diagonal <- (correlations$F1+correlations$F2+correlations$F1_nozero+correlations$F2_nozero)/4
C <- C_display_diagonal
diag(C) <- NA_real_
write.table(C, file.path(out, "pooled_four_atlas_C_four_matrix_mean.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(C_display_diagonal, file.path(out, "pooled_four_atlas_C_display_diagonal.tsv"), sep="\t", quote=FALSE, col.names=NA)

# Step5. 按原文计算每个 EcoType 的 AvgDif，并逐行置换 1000 次
membership <- lapply(ecotypes, function(ecotype) which(rownames(C) %in% members$ID[members$Ecotype==ecotype]))
names(membership) <- ecotypes
reference <- matrix(NA_real_, nrow=length(ecotypes), ncol=4, dimnames=list(ecotypes, c("AvgDif", "DefinedDeltaMembers", "MeanWithin", "DefinedWithinMembers")))

for (ecotype in ecotypes) {
    member_index <- membership[[ecotype]]
    within <- rowMeans(C[member_index, member_index, drop=FALSE], na.rm=TRUE)
    outside <- rowMeans(C[member_index, -member_index, drop=FALSE], na.rm=TRUE)
    delta <- within-outside
    reference[ecotype, ] <- c(mean(delta[is.finite(delta)]), sum(is.finite(delta)), mean(within[is.finite(within)]), sum(is.finite(within)))
}

permutations <- 1000
seed <- 20261125
set.seed(seed)
null <- matrix(NA_real_, nrow=permutations, ncol=length(ecotypes), dimnames=list(NULL, ecotypes))
null_within <- null
null_defined <- null

for (iteration in seq_len(permutations)) {
    shuffled <- C
    for (row_index in seq_len(nrow(C))) {
        order_index <- sample.int(ncol(C))
        shuffled[row_index, ] <- C[row_index, order_index]
        original_diagonal_position <- which(order_index==row_index)
        diagonal_value <- shuffled[row_index, row_index]
        shuffled[row_index, row_index] <- shuffled[row_index, original_diagonal_position]
        shuffled[row_index, original_diagonal_position] <- diagonal_value
    }

    for (ecotype in ecotypes) {
        member_index <- membership[[ecotype]]
        within <- rowMeans(shuffled[member_index, member_index, drop=FALSE], na.rm=TRUE)
        outside <- rowMeans(shuffled[member_index, -member_index, drop=FALSE], na.rm=TRUE)
        delta <- within-outside
        null[iteration, ecotype] <- mean(delta[is.finite(delta)])
        null_within[iteration, ecotype] <- mean(within[is.finite(within)])
        null_defined[iteration, ecotype] <- sum(is.finite(delta))
    }
}

null_mean <- colMeans(null, na.rm=TRUE)
null_sd <- apply(null, 2, sd, na.rm=TRUE)
Z <- (reference[, "AvgDif"]-null_mean)/null_sd
P_normal <- pnorm(Z, lower.tail=FALSE)
P_empirical <- vapply(ecotypes, function(ecotype) (1+sum(null[,ecotype]>=reference[ecotype, "AvgDif"], na.rm=TRUE))/(1+sum(is.finite(null[,ecotype]))), numeric(1))

result <- data.frame(Ecotype=ecotypes, FrozenMembers=as.integer(table(factor(members$Ecotype, levels=ecotypes))), DefinedDeltaMembers=reference[, "DefinedDeltaMembers"], AvgDif=reference[, "AvgDif"], NullMean=null_mean, NullSD=null_sd, Z=Z, P_normal_greater=P_normal, P_empirical_greater=P_empirical, MeanWithin=reference[, "MeanWithin"], DefinedWithinMembers=reference[, "DefinedWithinMembers"], ValidPermutations=colSums(is.finite(null)), MinNullDefinedMembers=apply(null_defined, 2, min, na.rm=TRUE), MaxNullDefinedMembers=apply(null_defined, 2, max, na.rm=TRUE), Permutations=permutations, Seed=seed)
result$Q_BH8_normal <- p.adjust(result$P_normal_greater, method="BH")
result$Q_BH8_empirical <- p.adjust(result$P_empirical_greater, method="BH")
result$Cell2021Nominal <- result$P_normal_greater<0.05 & result$AvgDif>result$NullMean
result$EmpiricalNominal <- result$P_empirical_greater<0.05 & result$AvgDif>result$NullMean

fwrite(result, file.path(out, "pooled_four_atlas_ecotype_cooccurrence_8tests.tsv"), sep="\t")
write.table(null, file.path(out, "pooled_four_atlas_AvgDif_permutation_null.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(null_within, file.path(out, "pooled_four_atlas_within_mean_permutation_null.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Step6. 保存 state 覆盖、方法边界和环境
coverage <- data.frame(CellState=members$ID, CellType=members$CellType, State=members$State, Ecotype=members$Ecotype, AssignedCells=as.integer(rowSums(numerator)), PositiveSamples=as.integer(rowSums(numerator>0)), SamplesWithParentType=as.integer(rowSums(D2>0)), DefinedEdges=as.integer(rowSums(is.finite(C))))
fwrite(coverage, file.path(out, "pooled_four_atlas_ecotype_member_coverage.tsv"), sep="\t")
print(result[, c("Ecotype", "FrozenMembers", "DefinedDeltaMembers", "AvgDif", "Z", "P_normal_greater", "P_empirical_greater", "Cell2021Nominal")])
capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
writeLines(c(
    "Primary Cell 2021 display test: Z=(AvgDif-reference null mean)/null SD; one-sided normal-tail P; no BH used for the figure star.",
    "Sensitivity: one-sided empirical P=(1+count(null>=observed))/(1000+1), with BH8 for both P definitions saved but not substituted for the paper-style star.",
    "The four prespecified scRNA-seq atlases are pooled at the biological-sample level before F1/F2 and C are calculated, matching the STAR Methods structure.",
    "C is the strict mean of F1, F2, F1-nozero and F2-nozero Spearman matrices; an undefined component keeps the edge undefined.",
    "All 79 frozen EcoType-member states enter C; the F1 denominator includes every formally assigned state in each parent cell type.",
    "No sample, cell, state or dataset was selected according to the resulting P values.",
    "The one-million joint probability test is a separate analysis and is not claimed by this script."
), file.path(out, "method_boundary.txt"))
