# 按Cell2021的文献marker对照与单细胞生态型共关联步骤计算初步解释。
# 在服务器独立ecotyper_interpretation环境运行；读取冻结结果，写入新审阅目录。
# 不修改模型；缺失值规则与原文未明确之处另行记录，不声称数值等价。
library(data.table)

# Step1. 读取冻结定义和文献面板，记录本次环境与输入校验
root <- "/home/mazekai/IBD_EcoTyper"
out <- file.path(root, "06_validation/ecotyper_original_method_annotation_20261010")
input <- file.path(out, "inputs")
states <- fread(file.path(input, "frozen_states.tsv"), data.table=FALSE)
ecotypes <- fread(file.path(input, "frozen_ecotypes.tsv"), data.table=FALSE)
markers <- fread(file.path(input, "discovery_associated_genes.tsv"), data.table=FALSE)
panels <- fread(file.path(input, "literature_marker_panels.tsv"), data.table=FALSE)
stopifnot(nrow(states)==53, nrow(ecotypes)==7, !anyDuplicated(states$CellState))
input_files <- list.files(input, full.names=TRUE)
input_hash <- tools::md5sum(input_files)
capture.output(sessionInfo(), file=file.path(out, "analysis_sessionInfo.txt"))
print(c(States=nrow(states), Ecotypes=nrow(ecotypes), AssociatedGenes=nrow(markers)))

# Step2. 文献程序匹配，保留实际命中基因；发现marker不是全表达矩阵
overlap <- list()
for (state in states$CellState) {
    genes <- markers[markers$CellState==state, ]
    genes <- genes[order(genes$WithinStateRank), ]
    for (program in unique(panels$Program)) {
        panel <- panels[panels$Program==program, ]
        hits <- genes[genes$Gene %in% panel$Gene, ]
        overlap[[length(overlap)+1]] <- data.frame(CellState=state, Program=program, Hits=nrow(hits), PanelSize=nrow(panel), HitGenes=paste(hits$Gene, collapse=";"), HitDiscoveryRanks=paste(hits$WithinStateRank, collapse=";"), Top20Hits=sum(hits$WithinStateRank<=20), ExpectedLineage=panel$ExpectedLineage[1], Publication=panel$Publication[1], URL=panel$URL[1], PanelStatus=panel$PanelStatus[1])
    }
}
overlap <- do.call(rbind, overlap)
fwrite(overlap, file.path(out, "state_literature_marker_correspondence.tsv"), sep="\t")
print(overlap[overlap$CellState %in% c("Fibroblast_S02", "Goblet_S01", "Glial_S05") & overlap$Hits>=2, c("CellState", "Program", "Hits", "HitGenes")], row.names=FALSE)

# Step3. 两个恢复队列分别分析，不将回肠与结肠细胞合并
datasets <- c(Colon_GSE214695="GSE214695_scRNA_recovery_20261009", Ileum_IBDverse="IBDverse_scRNA_recovery_20261009")
eligibility <- list()
coassociation <- list()
for (dataset in names(datasets)) {
    cells_path <- file.path(root, "06_validation", datasets[[dataset]], "audit/cell_assignments_all.tsv")
    cells_hash <- tools::md5sum(cells_path)
    cells <- fread(cells_path, data.table=FALSE)
    stopifnot(!anyDuplicated(cells$ID), all(cells$FormalRecovery))
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

# Step7. 保存可核对的解释入口，不将marker重叠写成统计富集或最终标签
eligible <- do.call(rbind, eligibility)
fwrite(eligible, file.path(out, "marker_eligibility_both_datasets.tsv"), sep="\t")
fwrite(do.call(rbind, coassociation), file.path(out, "ecotype_coassociation_both_datasets.tsv"), sep="\t")
states$InterpretationStatus <- "文献读出对应与恢复支持分开；尚未完成原文七指标外部marker和整转录组GSEA"
states$LiteratureCorrespondence <- NA_character_
for (i in seq_len(nrow(states))) {
    hits <- overlap[overlap$CellState==states$CellState[i] & overlap$Hits>0, ]
    hits <- hits[order(-hits$Top20Hits, -hits$Hits), ]
    states$LiteratureCorrespondence[i] <- paste(paste0(hits$Program, ":", hits$HitGenes, " [", hits$Publication, "]"), collapse=" | ")
}
fwrite(states, file.path(out, "53_states_preliminary_literature_evidence.tsv"), sep="\t")
stopifnot(identical(input_hash, tools::md5sum(input_files)))
write.table(data.frame(File=names(input_hash), MD5=unname(input_hash)), file.path(out, "input_md5_unchanged.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
print(table(eligible$Dataset, eligible$PaperMarkerEligible))
writeLines("Completed: literature marker correspondence; paper marker eligibility; F1/F2 four Spearman matrices; 1000 row-shuffle coassociation statistics. Not completed: seven-metric external marker meta-analysis, LOOCV, full-transcriptome fgsea, ligand-receptor test, million-permutation global CE test.", file.path(out, "analysis_completed.txt"))
