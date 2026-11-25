# 按 Cell 2021 STAR Methods 在四个单细胞数据中验证冻结的 8 个 EcoType 共现。
# 在服务器运行；读取正式 recovery 分配和原始样本注释，不重新拟合 state 或 EcoType。
# 主检验为每数据集每 EcoType 1000 次逐行置换的单侧经验 P，同时保留 Z 和 BH32。
library(data.table)

# Step1. 定位冻结模型、四个 recovery 任务表和输出
args <- commandArgs(trailingOnly=TRUE)
run_root <- if (length(args)>=1) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
out <- file.path(run_root, "output/46_scRNA_ecotype_cooccurrence_Cell2021")
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

# Step2. 逐数据集合并父类型分配，保留全部输入细胞作为 F2 分母
all_results <- list()
all_coverage <- list()
all_input_audit <- list()
permutations <- 1000
seed <- 20261125

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
    stopifnot(all(grepl("_S[0-9]+$", assigned$CellState)))

    present_members <- members[members$CellType %in% jobs$CellType, ]
    sample_names <- unique(annotation$BiologicalSampleID)
    numerator <- table(factor(assigned$CellState, levels=present_members$ID), factor(assigned$BiologicalSampleID, levels=sample_names))
    all_denominator <- table(factor(annotation$CellType, levels=jobs$CellType), factor(annotation$BiologicalSampleID, levels=sample_names))
    assigned_denominator <- table(factor(assigned$CellType, levels=jobs$CellType), factor(assigned$BiologicalSampleID, levels=sample_names))
    D1 <- assigned_denominator[present_members$CellType, , drop=FALSE]
    D2 <- all_denominator[present_members$CellType, , drop=FALSE]
    F1 <- numerator/D1
    F2 <- numerator/D2
    F1[D1==0] <- NA_real_
    F2[D2==0] <- NA_real_
    rownames(F1) <- rownames(F2) <- present_members$ID
    stopifnot(all(F1[is.finite(F1)]>=0 & F1[is.finite(F1)]<=1), all(F2[is.finite(F2)]>=0 & F2[is.finite(F2)]<=1))

    write.table(F1, file.path(out, paste0(dataset, "_F1_assigned_denominator.tsv")), sep="\t", quote=FALSE, col.names=NA)
    write.table(F2, file.path(out, paste0(dataset, "_F2_all_cells_denominator.tsv")), sep="\t", quote=FALSE, col.names=NA)
    write.table(D2, file.path(out, paste0(dataset, "_parent_type_input_counts.tsv")), sep="\t", quote=FALSE, col.names=NA)

    # Step3. 构造四套 Spearman 相关矩阵并严格平均
    F1_nozero <- F1
    F2_nozero <- F2
    F1_nozero[F1_nozero==0 & !is.na(F1_nozero)] <- NA_real_
    F2_nozero[F2_nozero==0 & !is.na(F2_nozero)] <- NA_real_
    fractions <- list(F1=F1, F2=F2, F1_nozero=F1_nozero, F2_nozero=F2_nozero)
    correlations <- list()

    for (version in names(fractions)) {
        correlations[[version]] <- suppressWarnings(cor(t(fractions[[version]]), method="spearman", use="pairwise.complete.obs"))
        pair_n <- tcrossprod(1L*!is.na(fractions[[version]]))
        write.table(correlations[[version]], file.path(out, paste0(dataset, "_C_", version, ".tsv")), sep="\t", quote=FALSE, col.names=NA)
        write.table(pair_n, file.path(out, paste0(dataset, "_pair_N_", version, ".tsv")), sep="\t", quote=FALSE, col.names=NA)
    }

    C <- (correlations$F1+correlations$F2+correlations$F1_nozero+correlations$F2_nozero)/4
    diag(C) <- NA_real_
    write.table(C, file.path(out, paste0(dataset, "_C_four_matrix_mean.tsv")), sep="\t", quote=FALSE, col.names=NA)

    # Step4. 计算成员内减成员外 AvgDif，并执行 1000 次逐行置换
    membership <- lapply(ecotypes, function(ecotype) which(rownames(C) %in% members$ID[members$Ecotype==ecotype]))
    names(membership) <- ecotypes
    reference <- matrix(NA_real_, nrow=length(ecotypes), ncol=2, dimnames=list(ecotypes, c("Statistic", "DefinedMembers")))

    for (ecotype in ecotypes) {
        member_index <- membership[[ecotype]]
        within <- rowMeans(C[member_index, member_index, drop=FALSE], na.rm=TRUE)
        outside <- rowMeans(C[member_index, -member_index, drop=FALSE], na.rm=TRUE)
        delta <- within-outside
        reference[ecotype, ] <- c(mean(delta[is.finite(delta)]), sum(is.finite(delta)))
    }

    set.seed(seed)
    null <- matrix(NA_real_, nrow=permutations, ncol=length(ecotypes), dimnames=list(NULL, ecotypes))
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
            null_defined[iteration, ecotype] <- sum(is.finite(delta))
        }
    }

    null_mean <- colMeans(null, na.rm=TRUE)
    null_sd <- apply(null, 2, sd, na.rm=TRUE)
    empirical_p <- vapply(ecotypes, function(ecotype) (1+sum(null[,ecotype]>=reference[ecotype, "Statistic"], na.rm=TRUE))/(1+sum(is.finite(null[,ecotype]))), numeric(1))
    result <- data.frame(Dataset=dataset, Ecotype=ecotypes, FrozenMembers=as.integer(table(factor(members$Ecotype, levels=ecotypes))), InputMembers=lengths(membership), DefinedReferenceMembers=reference[,"DefinedMembers"], AvgDif=reference[,"Statistic"], NullMean=null_mean, NullSD=null_sd, Z=(reference[,"Statistic"]-null_mean)/null_sd, P_empirical_greater=empirical_p, P_normal_greater=pnorm((reference[,"Statistic"]-null_mean)/null_sd, lower.tail=FALSE), ValidPermutations=colSums(is.finite(null)), MinNullDefinedMembers=apply(null_defined, 2, min, na.rm=TRUE), MaxNullDefinedMembers=apply(null_defined, 2, max, na.rm=TRUE), Permutations=permutations, Seed=seed)
    all_results[[dataset]] <- result

    coverage <- data.frame(Dataset=dataset, CellState=rownames(C), Ecotype=members$Ecotype[match(rownames(C), members$ID)], AssignedCells=as.integer(rowSums(numerator)), PositiveSamples=as.integer(rowSums(numerator>0)), DefinedEdges=as.integer(rowSums(is.finite(C))))
    all_coverage[[dataset]] <- coverage
    all_input_audit[[dataset]] <- data.frame(Dataset=dataset, ParentTypes=nrow(jobs), BiologicalSamples=length(sample_names), InputCells=nrow(annotation), AssignedCells=nrow(assigned), MemberStates=nrow(C), FiniteCEdges=sum(is.finite(C))/2)
    write.table(null, file.path(out, paste0(dataset, "_permutation_null.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
    print(result[, c("Ecotype", "DefinedReferenceMembers", "AvgDif", "Z", "P_empirical_greater")])
}

# Step5. 固定 32 个检验的 BH 家族并保存总表
results <- rbindlist(all_results)
results$Q_BH32_empirical <- p.adjust(results$P_empirical_greater, method="BH")
results$NominalSignificant <- results$P_empirical_greater<0.05 & results$AvgDif>results$NullMean
results$BH32Significant <- results$Q_BH32_empirical<0.05 & results$AvgDif>results$NullMean
setcolorder(results, c("Dataset", "Ecotype", "FrozenMembers", "InputMembers", "DefinedReferenceMembers", "AvgDif", "NullMean", "NullSD", "Z", "P_empirical_greater", "Q_BH32_empirical", "P_normal_greater", "NominalSignificant", "BH32Significant", "ValidPermutations", "MinNullDefinedMembers", "MaxNullDefinedMembers", "Permutations", "Seed"))
fwrite(results, file.path(out, "four_dataset_ecotype_cooccurrence_32tests.tsv"), sep="\t")
fwrite(rbindlist(all_coverage), file.path(out, "four_dataset_ecotype_member_coverage.tsv"), sep="\t")
fwrite(rbindlist(all_input_audit), file.path(out, "four_dataset_input_audit.tsv"), sep="\t")

summary <- results[, .(DatasetsTested=.N, PositiveAvgDif=sum(AvgDif>NullMean), NominalSignificant=sum(NominalSignificant), BH32Significant=sum(BH32Significant), MedianZ=median(Z, na.rm=TRUE), MaxZ=max(Z, na.rm=TRUE)), by=Ecotype]
fwrite(summary, file.path(out, "ecotype_recovery_summary.tsv"), sep="\t")
print(summary)
capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
writeLines(c("Primary test: one-sided greater-tail empirical P=(1+count(null>=observed))/(1000+1).", "Multiple testing: BH across the prespecified 4 datasets x 8 EcoTypes = 32 tests.", "C uses the strict mean of F1, F2, F1-nozero and F2-nozero Spearman matrices; any undefined component keeps the edge undefined.", "Only the 79 frozen ecotype-member states enter C; the F1 denominator includes all formally assigned states in each parent cell type.", "No sample or cell selection was based on the resulting P values."), file.path(out, "method_boundary.txt"))
