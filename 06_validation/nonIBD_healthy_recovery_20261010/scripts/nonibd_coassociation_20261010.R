# 沿用此前 Cell2021 四套相关矩阵与 1000 次生态型行置换，仅检查健康子集共同变化。
# 在服务器独立 ecotyper_interpretation 环境运行；读取官方恢复的小表，不修改源码或数据。
# 6/10 位供者的相关边覆盖有限；统计不作为正式生态型确认，不进行文献命名。
library(data.table)

# Step1. 读取本轮冻结状态、成员和健康细胞分配，保存环境记录
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
out <- file.path(root, "audit")
states <- fread(file.path(out, "frozen_states.tsv"), data.table=FALSE)
ecotypes <- fread(file.path(out, "frozen_ecotypes.tsv"), data.table=FALSE)
datasets <- c(IBDverse_Healthy_Ileum_20261010="IBDverse", GSE214695_HC_Colon_20261010="GSE214695")
capture.output(sessionInfo(), file=file.path(out, "coassociation_sessionInfo.txt"))
stopifnot(nrow(states)==53, nrow(ecotypes)==7, !anyDuplicated(states$CellState))

# Step2. 下列 F1/F2、缺失值、行置换和统计量沿用上一轮审阅方法
eligibility <- list()
coassociation <- list()
for (dataset in names(datasets)) {
    cells_path <- file.path(out, dataset, "cell_assignments_all.tsv")
    cells_hash <- tools::md5sum(cells_path)
    cells <- fread(cells_path, data.table=FALSE)
    stopifnot(!anyDuplicated(cells$ID))
    print(c(Dataset=dataset, Cells=nrow(cells), Samples=length(unique(cells$Sample)), Donors=length(unique(cells$Donor))))
    cells$CellState <- ifelse(cells$State=="Unassigned", "Unassigned", paste(cells$CellType, cells$State, sep="_"))
    stopifnot(all(cells$CellState %in% c(states$CellState, "Unassigned")))
    sample_donor <- unique(cells[c("Sample", "Donor")])
    stopifnot(!anyDuplicated(sample_donor$Sample))

    # 原文≥5细胞规则仅用于外部marker筛选，不擅自作为生态型相关矩阵的排除阈值。
    assigned <- table(factor(cells$CellState, levels=states$CellState))
    eligible <- data.frame(Dataset=dataset, CellState=states$CellState, CellType=states$CellType, AssignedCells=as.integer(assigned), LineagePresent=states$CellType %in% cells$CellType)
    eligible$PaperMarkerEligible <- eligible$LineagePresent & eligible$AssignedCells>=5
    eligibility[[dataset]] <- eligible
    fwrite(eligible, file.path(out, paste0(dataset, "_marker_eligibility.tsv")), sep="\t")

    # Step4. 原文F1/F2：分母分别为可分配细胞、该父类型全部输入细胞
    present <- states[states$CellType %in% cells$CellType, ]
    sample_names <- unique(cells$Sample)
    numerator <- table(factor(cells$CellState, levels=present$CellState), factor(cells$Sample, levels=sample_names))
    all_denominator <- table(cells$CellType, factor(cells$Sample, levels=sample_names))
    retained <- cells[cells$State!="Unassigned", ]
    assigned_denominator <- table(factor(retained$CellType, levels=rownames(all_denominator)), factor(retained$Sample, levels=sample_names))
    D1 <- assigned_denominator[present$CellType, , drop=FALSE]
    D2 <- all_denominator[present$CellType, , drop=FALSE]
    F1 <- numerator / D1
    F2 <- numerator / D2
    F1[D1==0] <- NA_real_
    F2[D2==0] <- NA_real_
    rownames(F1) <- rownames(F2) <- present$CellState
    stopifnot(all(F1[is.finite(F1)]>=0), all(F1[is.finite(F1)]<=1), all(F2[is.finite(F2)]>=0), all(F2[is.finite(F2)]<=1))
    write.table(F1, file.path(out, paste0(dataset, "_F1_assigned_denominator.tsv")), sep="\t", quote=FALSE, col.names=NA)
    write.table(F2, file.path(out, paste0(dataset, "_F2_all_cells_denominator.tsv")), sep="\t", quote=FALSE, col.names=NA)
    write.table(D2, file.path(out, paste0(dataset, "_parent_type_input_counts.tsv")), sep="\t", quote=FALSE, col.names=NA)

    # Step5. 原文四套Spearman矩阵，零改NA的两套用于减轻欠抽样影响
    F1_nozero <- F1
    F2_nozero <- F2
    F1_nozero[F1_nozero==0 & !is.na(F1_nozero)] <- NA_real_
    F2_nozero[F2_nozero==0 & !is.na(F2_nozero)] <- NA_real_
    fractions <- list(F1=F1, F2=F2, F1_nozero=F1_nozero, F2_nozero=F2_nozero)
    correlations <- list()
    for (version in names(fractions)) {
        correlations[[version]] <- suppressWarnings(cor(t(fractions[[version]]), method="spearman", use="pairwise.complete.obs"))
        pair_n <- tcrossprod(1L * !is.na(fractions[[version]]))
        write.table(correlations[[version]], file.path(out, paste0(dataset, "_C_", version, ".tsv")), sep="\t", quote=FALSE, col.names=NA)
        write.table(pair_n, file.path(out, paste0(dataset, "_pair_N_", version, ".tsv")), sep="\t", quote=FALSE, col.names=NA)
    }
    # 原文未写缺失值合并方式；主分析要求四套均可计算才保留该边，避免悄悄更换平均分母。
    C <- (correlations[[1]] + correlations[[2]] + correlations[[3]] + correlations[[4]]) / 4
    diag(C) <- NA_real_
    write.table(C, file.path(out, paste0(dataset, "_C_four_matrix_mean.tsv")), sep="\t", quote=FALSE, col.names=NA)
    edge_coverage <- data.frame(CellState=rownames(C), DefinedEdges=rowSums(is.finite(C)), PositiveSamples=rowSums(numerator>0))
    fwrite(edge_coverage, file.path(out, paste0(dataset, "_correlation_coverage.tsv")), sep="\t")

    # Step6. 原文生态型统计量：成员内关联减成员外关联，再在成员间取平均
    # 此短函数只用于重复计算同一个置换统计量；原文未明确NA处理，有限边求均值并报告覆盖。
    average_difference <- function(matrix, member_index) {
        within <- rowMeans(matrix[member_index, member_index, drop=FALSE], na.rm=TRUE)
        outside <- rowMeans(matrix[member_index, -member_index, drop=FALSE], na.rm=TRUE)
        delta <- within - outside
        c(Statistic=mean(delta[is.finite(delta)]), DefinedMembers=sum(is.finite(delta)))
    }
    membership <- lapply(ecotypes$Members, function(x) which(rownames(C) %in% strsplit(x, ";", fixed=TRUE)[[1]]))
    names(membership) <- ecotypes$Ecotype
    reference <- t(vapply(membership, function(i) average_difference(C, i), numeric(2)))
    set.seed(20261010)  # 本项目固定种子；论文未指定种子
    permutations <- 1000  # Cell2021单个生态型检验次数，非全局百万次检验
    null <- matrix(NA_real_, nrow=permutations, ncol=nrow(ecotypes), dimnames=list(NULL, ecotypes$Ecotype))
    null_members <- null
    for (iteration in seq_len(permutations)) {
        shuffled <- C
        for (row in seq_len(nrow(C))) {
            order <- sample.int(ncol(C))
            shuffled[row, ] <- C[row, order]
            original_diagonal_position <- which(order==row)
            swapped <- shuffled[row, row]
            shuffled[row, row] <- shuffled[row, original_diagonal_position]
            shuffled[row, original_diagonal_position] <- swapped
        }
        values <- t(vapply(membership, function(i) average_difference(shuffled, i), numeric(2)))
        null[iteration, ] <- values[, "Statistic"]
        null_members[iteration, ] <- values[, "DefinedMembers"]
    }
    null_mean <- colMeans(null, na.rm=TRUE)
    null_sd <- apply(null, 2, sd, na.rm=TRUE)
    result <- data.frame(Dataset=dataset, Ecotype=ecotypes$Ecotype, FrozenMembers=ecotypes$MemberStates, InputLineageMembers=lengths(membership), DefinedReferenceMembers=reference[, "DefinedMembers"], AvgDif=reference[, "Statistic"], NullMean=null_mean, NullSD=null_sd, Z=(reference[, "Statistic"]-null_mean)/null_sd, ValidPermutations=colSums(is.finite(null)), MinNullDefinedMembers=apply(null_members, 2, min), MaxNullDefinedMembers=apply(null_members, 2, max), Permutations=permutations, Seed=20261010)
    # 不推断原文未明确的P值尾部；Z及零分布可检查，缺失成员不等同完整生态型验证。
    coassociation[[dataset]] <- result
    fwrite(result, file.path(out, paste0(dataset, "_ecotype_coassociation.tsv")), sep="\t")
    write.table(null, file.path(out, paste0(dataset, "_permutation_null.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
    print(result[, c("Ecotype", "DefinedReferenceMembers", "AvgDif", "Z")], row.names=FALSE)
    stopifnot(identical(cells_hash, tools::md5sum(cells_path)))
}


# Step7. 单独保存共同变化支持与可计算成员数，不进行状态命名
fwrite(do.call(rbind, eligibility), file.path(out, "healthy_marker_eligibility.tsv"), sep="\t")
fwrite(do.call(rbind, coassociation), file.path(out, "healthy_ecotype_coassociation.tsv"), sep="\t")
writeLines("Completed healthy-subset F1/F2 and 1000 row-shuffle coassociation; few donors, missing members and sparse correlations restrict biological interpretation. No literature naming or clinical associations.", file.path(out, "coassociation_completed.txt"))