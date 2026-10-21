# 对冻结恢复细胞做全基因状态间 DEG，另按疾病和 Region3 分层，供 n4 两版比较。
# 服务器 ecotyper_interpretation 环境运行；不改状态，不重做 recovery，不合并两个数据集。
# Seurat3.1.3 原生 avg_logFC 是自然对数；原生0.05过滤阈值不变，输出额外换算 log2FC。
library(Seurat)
library(Matrix)
library(data.table)
library(future)
plan(multicore, workers=4)  # 同一Seurat检验分配4个Linux worker；不改基因、细胞、P值算法
stopifnot(as.character(packageVersion("Seurat")) == "3.1.3")
out <- "/home/mazekai/IBD_EcoTyper/06_validation/paper_aligned_recovery_characterization_20261012"
capture.output(sessionInfo(), file=file.path(out, "00_audit/marker_sessionInfo_parallel4.txt"))

# Step1. 先用可手算的计数确认归一化、fold change 和 Bonferroni 默认行为
set.seed(20261012)
toy <- matrix(rpois(120, 4), nrow=10, dimnames=list(paste0("g", 1:10), paste0("c", 1:12)))
toy[1, 1:6] <- 40
toy_log <- NormalizeData(toy, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
toy_expected <- log1p(sweep(toy, 2, colSums(toy), "/") * 10000)
stopifnot(max(abs(toy_log-toy_expected)) < 1e-10)
toy_de <- FindMarkers(toy_log, cells.1=colnames(toy)[1:6], cells.2=colnames(toy)[7:12], min.pct=0.1, logfc.threshold=0.05, verbose=FALSE)
toy_fc <- log(rowMeans(expm1(toy_log[, 1:6]))+1)-log(rowMeans(expm1(toy_log[, 7:12]))+1)
stopifnot(max(abs(toy_de$avg_logFC-toy_fc[rownames(toy_de)])) < 1e-10)
stopifnot(max(abs(toy_de$p_val_adj-pmin(toy_de$p_val*nrow(toy), 1))) < 1e-10)
writeLines("PASS: NormalizeData ln(CP10k+1); native lnFC; native Bonferroni total genes; threshold0.05 unchanged", file.path(out, "00_audit/Seurat313_numerical_smoke_parallel4.txt"))
print("SEURAT313_NUMERICAL_SMOKE_PASS")
if ("--smoke-only" %in% commandArgs(TRUE)) quit(status=0)

# Step2. 每个数据集和父类型独立读取全表达矩阵；保存每个分层的可评价状态及细胞数
discovery_genes <- fread(file.path(out, "00_audit/discovery_gene_assignment.tsv"), data.table=FALSE)
directories <- list.dirs(file.path(out, "02_external_markers/full_counts"), recursive=FALSE, full.names=TRUE)
eligibility <- list()
for (dataset_dir in directories) {
    dataset <- basename(dataset_dir)
    for (type_dir in list.dirs(dataset_dir, recursive=FALSE, full.names=TRUE)) {
        cell_type <- basename(type_dir)
        counts <- as(readMM(gzfile(file.path(type_dir, "counts.mtx.gz"))), "dgCMatrix")
        genes <- readLines(file.path(type_dir, "genes.tsv"))
        cells <- fread(file.path(type_dir, "cells.tsv"), data.table=FALSE)
        stopifnot(nrow(counts)==length(genes), ncol(counts)==nrow(cells), !anyDuplicated(genes), !anyDuplicated(cells$ID), all(colSums(counts)>0))
        rownames(counts) <- genes
        colnames(counts) <- cells$ID
        data <- NormalizeData(counts, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
        cpm_mean <- Matrix::rowMeans(counts %*% Diagonal(x=1e6/Matrix::colSums(counts)))
        output_dir <- file.path(out, "02_external_markers/DE_checkpoints_parallel4", dataset, cell_type)
        dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
        fwrite(data.frame(Gene=genes, MeanCPM=cpm_mean), file.path(output_dir, "gene_expression.tsv"), sep="\t")
        scope_table <- rbind(data.frame(Scope="Dataset", Stratum="All"), data.frame(Scope="DiseaseGroup", Stratum=sort(unique(cells$DiseaseGroup))), data.frame(Scope="Region3", Stratum=sort(unique(cells$Region3))))
        for (j in seq_len(nrow(scope_table))) {
            scope <- scope_table$Scope[j]
            stratum <- scope_table$Stratum[j]
            group_cells <- if (scope=="Dataset") cells else cells[cells[[scope]]==stratum, ]
            state_counts <- table(group_cells$State[group_cells$State!="Unassigned"])
            allowed <- names(state_counts)[state_counts>=5]  # 原文排除<5细胞状态及其细胞
            eligible_cells <- group_cells[group_cells$State %in% allowed, ]
            eligible_rows <- data.frame(Dataset=dataset, CellType=cell_type, Scope=scope, Stratum=stratum, State=names(state_counts), AssignedCells=as.integer(state_counts), Eligible=names(state_counts) %in% allowed, EvaluableContrast=length(allowed)>=2)
            eligibility[[length(eligibility)+1]] <- eligible_rows
            tag <- paste(scope, stratum, sep="_")
            if (length(allowed)<2) {
                writeLines("NOT_EVALUABLE: fewer than two states with >=5 cells", file.path(output_dir, paste0(tag, "_not_evaluable.txt")))
                next
            }
            print(c(Dataset=dataset, CellType=cell_type, Scope=scope, Stratum=stratum, Cells=nrow(eligible_cells), States=length(allowed)))

            # Step3. 保存全基因测试结果，同时为发现模型候选基因保留未通过预筛选时的真实FC；P不补0
            results <- list()
            folds <- list()
            candidate_genes <- intersect(genes, discovery_genes$Gene[discovery_genes$CellType==cell_type])
            for (state in allowed) {
                own <- eligible_cells$ID[eligible_cells$State==state]
                other <- eligible_cells$ID[eligible_cells$State!=state]
                fold <- log(Matrix::rowMeans(expm1(data[candidate_genes, own, drop=FALSE]))+1)-log(Matrix::rowMeans(expm1(data[candidate_genes, other, drop=FALSE]))+1)
                folds[[state]] <- data.frame(Gene=candidate_genes, State=state, log2FC=as.numeric(fold/log(2)))
                result <- tryCatch(FindMarkers(data[, eligible_cells$ID, drop=FALSE], cells.1=own, cells.2=other, min.pct=0.1, logfc.threshold=0.05, test.use="wilcox", verbose=FALSE), error=function(e) e)
                if (inherits(result, "error")) {
                    if (!grepl("No features pass", conditionMessage(result))) stop(result)
                    writeLines(conditionMessage(result), file.path(output_dir, paste0(tag, "_", state, "_empty.txt")))
                    next
                }
                result$Gene <- rownames(result)
                result$State <- state
                result$log2FC <- result$avg_logFC/log(2)
                result$OwnCells <- length(own)
                result$OtherCells <- length(other)
                results[[state]] <- result
            }
            fwrite(do.call(rbind, folds), file.path(output_dir, paste0(tag, "_candidate_full_FC.tsv")), sep="\t")
            if (length(results)>0) fwrite(do.call(rbind, results), file.path(output_dir, paste0(tag, "_DEG.tsv")), sep="\t")
        }
        fwrite(do.call(rbind, eligibility), file.path(out, "02_external_markers/DE_stratum_eligibility_parallel4.tsv"), sep="\t")
        writeLines("PASS", file.path(output_dir, "completed.txt"))
    }
}
writeLines("PASS", file.path(out, "02_external_markers/external_DE_completed_parallel4.txt"))
print("EXTERNAL_DE_COMPLETE")
