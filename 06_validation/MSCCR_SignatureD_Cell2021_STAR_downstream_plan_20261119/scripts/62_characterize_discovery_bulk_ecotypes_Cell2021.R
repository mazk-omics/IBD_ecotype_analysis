# 使用 discovery bulk cohort 表征冻结的 8 个 EcoType，并与 pooled scRNA 外部恢复合并证据。
# 在服务器运行；读取 EcoTyper 正式 discovery 输出和原始注释，不重新训练或修改模型。
# Disease、Inflammation、Region3 分开检验；全部 2025 样本进入连续丰度分析。
library(data.table)
library(ggplot2)
library(ragg)
library(svglite)

# Step1. 读入 discovery abundance、official assignment 和样本注释
args <- commandArgs(trailingOnly=TRUE)
run_root <- if (length(args)>=1) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
out <- file.path(run_root, "output/50_discovery_bulk_ecotype_characterization")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

discovery_dir <- "/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101/Custom/Ecotypes/discovery"
annotation_path <- "/home/mazekai/IBD_EcoTyper/03_reference/MSCCR_SignatureD_IBD2025_TPM_20261101/MSCCR_CD_UC_2025_annotation.tsv"
pooled_scrna_path <- file.path(run_root, "output/48_pooled_scRNA_ecotype_cooccurrence_Cell2021/pooled_four_atlas_ecotype_cooccurrence_8tests.tsv")
stopifnot(file.exists(file.path(discovery_dir, "ecotype_abundance.txt")), file.exists(file.path(discovery_dir, "ecotype_assignment.txt")), file.exists(annotation_path), file.exists(pooled_scrna_path))

abundance_matrix <- as.matrix(read.delim(file.path(discovery_dir, "ecotype_abundance.txt"), row.names=1, check.names=FALSE))
assignment <- read.delim(file.path(discovery_dir, "ecotype_assignment.txt"), row.names=1, check.names=FALSE)
annotation <- fread(annotation_path, data.table=FALSE)
pooled_scrna <- fread(pooled_scrna_path, data.table=FALSE)
ecotypes <- paste0("E", 1:8)
stopifnot(identical(rownames(abundance_matrix), ecotypes), nrow(annotation)==ncol(abundance_matrix), setequal(annotation$ID, colnames(abundance_matrix)), !anyDuplicated(annotation$ID), !anyDuplicated(assignment$ID))

abundance <- as.data.frame(t(abundance_matrix), stringsAsFactors=FALSE)
abundance$ID <- rownames(abundance)
abundance <- merge(annotation, abundance, by="ID", all.y=TRUE, sort=FALSE)
abundance <- abundance[match(colnames(abundance_matrix), abundance$ID), ]
abundance$MaxEcotypeFromAbundance <- ecotypes[max.col(as.matrix(abundance[, ecotypes]), ties.method="first")]
abundance$OfficialAssigned <- abundance$ID %in% assignment$ID
abundance$OfficialEcotype <- assignment$Ecotype[match(abundance$ID, assignment$ID)]
stopifnot(nrow(abundance)==2025, sum(abundance$OfficialAssigned)==nrow(assignment), all(abs(rowSums(abundance[, ecotypes])-1)<1e-6))

assigned_match <- abundance[abundance$OfficialAssigned, c("ID", "MaxEcotypeFromAbundance", "OfficialEcotype")]
stopifnot(all(assigned_match$MaxEcotypeFromAbundance==assigned_match$OfficialEcotype))
fwrite(abundance, file.path(out, "discovery_2025_ecotype_abundance_and_metadata.tsv"), sep="\t")

# Step2. 汇总 official assignment 和全样本最高丰度类别
official_counts <- as.data.frame(table(factor(abundance$OfficialEcotype[!is.na(abundance$OfficialEcotype)], levels=ecotypes)), stringsAsFactors=FALSE)
names(official_counts) <- c("Ecotype", "OfficialAssignedSamples")
official_counts$OfficialAssignedFractionAmongAll <- official_counts$OfficialAssignedSamples/nrow(abundance)

max_counts <- as.data.frame(table(factor(abundance$MaxEcotypeFromAbundance, levels=ecotypes)), stringsAsFactors=FALSE)
names(max_counts) <- c("Ecotype", "HighestAbundanceSamples")
max_counts$HighestAbundanceFraction <- max_counts$HighestAbundanceSamples/nrow(abundance)
assignment_summary <- merge(official_counts, max_counts, by="Ecotype", sort=FALSE)
assignment_summary$TotalSamples <- nrow(abundance)
assignment_summary$OfficialAssignedTotal <- sum(abundance$OfficialAssigned)
assignment_summary$UnassignedTotal <- sum(!abundance$OfficialAssigned)
fwrite(assignment_summary, file.path(out, "discovery_ecotype_assignment_summary.tsv"), sep="\t")
print(assignment_summary)

# Step3. 对最高丰度类别做 phenotype level 一对其余 Fisher 富集
axes <- c("Disease", "Inflammation", "Region3")
assignment_enrichment <- list()

for (axis in axes) {
    levels_axis <- unique(abundance[[axis]])
    axis_results <- list()

    for (level_axis in levels_axis) {
        for (ecotype in ecotypes) {
            in_level <- abundance[[axis]]==level_axis
            is_ecotype <- abundance$MaxEcotypeFromAbundance==ecotype
            contingency <- matrix(c(sum(in_level & is_ecotype), sum(in_level & !is_ecotype), sum(!in_level & is_ecotype), sum(!in_level & !is_ecotype)), nrow=2, byrow=TRUE)
            test <- fisher.test(contingency, alternative="greater")
            axis_results[[length(axis_results)+1]] <- data.frame(Axis=axis, Level=level_axis, Ecotype=ecotype, InLevelEcotype=contingency[1,1], InLevelOther=contingency[1,2], OutsideEcotype=contingency[2,1], OutsideOther=contingency[2,2], OddsRatio=unname(test$estimate), P_FisherGreater=test$p.value)
        }
    }

    axis_results <- rbindlist(axis_results)
    axis_results$Q_BH_axis <- p.adjust(axis_results$P_FisherGreater, method="BH")
    assignment_enrichment[[axis]] <- axis_results
}

assignment_enrichment <- rbindlist(assignment_enrichment)
fwrite(assignment_enrichment, file.path(out, "discovery_highest_ecotype_phenotype_enrichment.tsv"), sep="\t")

# Step4. 对连续 abundance 做 level 对其余的 Wilcoxon，并另做 Region3 全局检验
abundance_long <- melt(as.data.table(abundance), id.vars=c("ID", "PatientID", "Disease", "Inflammation", "Region3", "OfficialAssigned", "OfficialEcotype", "MaxEcotypeFromAbundance"), measure.vars=ecotypes, variable.name="Ecotype", value.name="Abundance")
abundance_enrichment <- list()

for (axis in axes) {
    levels_axis <- unique(abundance[[axis]])
    axis_results <- list()

    for (level_axis in levels_axis) {
        for (ecotype in ecotypes) {
            values <- abundance_long[abundance_long$Ecotype==ecotype, ]
            in_level <- values[[axis]]==level_axis
            test <- wilcox.test(values$Abundance[in_level], values$Abundance[!in_level], alternative="greater", exact=FALSE)
            axis_results[[length(axis_results)+1]] <- data.frame(Axis=axis, Level=level_axis, Ecotype=ecotype, N_level=sum(in_level), N_other=sum(!in_level), MedianLevel=median(values$Abundance[in_level]), MedianOther=median(values$Abundance[!in_level]), MedianDifference=median(values$Abundance[in_level])-median(values$Abundance[!in_level]), P_WilcoxonGreater=test$p.value)
        }
    }

    axis_results <- rbindlist(axis_results)
    axis_results$Q_BH_axis <- p.adjust(axis_results$P_WilcoxonGreater, method="BH")
    abundance_enrichment[[axis]] <- axis_results
}

abundance_enrichment <- rbindlist(abundance_enrichment)
fwrite(abundance_enrichment, file.path(out, "discovery_ecotype_abundance_phenotype_enrichment.tsv"), sep="\t")

region_global <- rbindlist(lapply(ecotypes, function(ecotype) {
    values <- abundance_long[abundance_long$Ecotype==ecotype, ]
    test <- kruskal.test(Abundance~Region3, data=values)
    data.frame(Ecotype=ecotype, P_Kruskal=test$p.value)
}))
region_global$Q_BH8 <- p.adjust(region_global$P_Kruskal, method="BH")
fwrite(region_global, file.path(out, "discovery_ecotype_abundance_Region3_global.tsv"), sep="\t")

# Step5. 做 patient-aware sensitivity：Disease 用患者中位数，Inflammation 用患者内配对
patient_disease <- abundance_long[, .(Abundance=median(Abundance), Disease=unique(Disease)[1]), by=.(PatientID, Ecotype)]
patient_disease_results <- rbindlist(lapply(ecotypes, function(ecotype) {
    values <- patient_disease[patient_disease$Ecotype==ecotype, ]
    test <- wilcox.test(Abundance~Disease, data=values, exact=FALSE)
    data.frame(Ecotype=ecotype, CD_patients=sum(values$Disease=="CD"), UC_patients=sum(values$Disease=="UC"), MedianCD=median(values$Abundance[values$Disease=="CD"]), MedianUC=median(values$Abundance[values$Disease=="UC"]), P_WilcoxonTwoSided=test$p.value)
}))
patient_disease_results$Q_BH8 <- p.adjust(patient_disease_results$P_WilcoxonTwoSided, method="BH")
fwrite(patient_disease_results, file.path(out, "discovery_patient_median_ecotype_abundance_Disease.tsv"), sep="\t")

patient_inflammation <- abundance_long[, .(Abundance=median(Abundance)), by=.(PatientID, Inflammation, Ecotype)]
patient_inflammation_wide <- dcast(patient_inflammation, PatientID+Ecotype~Inflammation, value.var="Abundance")
paired_inflammation_results <- rbindlist(lapply(ecotypes, function(ecotype) {
    values <- patient_inflammation_wide[patient_inflammation_wide$Ecotype==ecotype & is.finite(patient_inflammation_wide$Inflamed) & is.finite(patient_inflammation_wide$`Non-inflamed`), ]
    test <- wilcox.test(values$Inflamed, values$`Non-inflamed`, paired=TRUE, exact=FALSE)
    data.frame(Ecotype=ecotype, PairedPatients=nrow(values), MedianWithinPatientDifference=median(values$Inflamed-values$`Non-inflamed`), P_PairedWilcoxon=test$p.value)
}))
paired_inflammation_results$Q_BH8 <- p.adjust(paired_inflammation_results$P_PairedWilcoxon, method="BH")
fwrite(paired_inflammation_results, file.path(out, "discovery_patient_paired_ecotype_abundance_Inflammation.tsv"), sep="\t")

# Step6. 分别绘制 Disease、Inflammation、Region3 abundance 箱线图
axis_colours <- list(Disease=c(CD="#C65D47", UC="#497EB3"), Inflammation=c(Inflamed="#B13A3A", `Non-inflamed`="#5B8E61"), Region3=c(Colon="#4E79A7", Ileum="#F28E2B", Rectum="#76B7B2"))
axis_widths <- c(Disease=145, Inflammation=145, Region3=165)

for (axis in axes) {
    plot_data <- copy(abundance_long)
    plot_data$Group <- factor(plot_data[[axis]], levels=names(axis_colours[[axis]]))
    plot_data$Ecotype <- factor(plot_data$Ecotype, levels=ecotypes)
    p <- ggplot(plot_data, aes(x=Group, y=Abundance, fill=Group)) +
        geom_boxplot(width=0.62, outlier.shape=NA, linewidth=0.35) +
        facet_wrap(~Ecotype, nrow=2) +
        scale_fill_manual(values=axis_colours[[axis]], drop=FALSE) +
        labs(title=paste0("Discovery bulk EcoType abundance by ", axis), x=NULL, y="EcoType abundance") +
        theme_classic(base_family="sans", base_size=7) +
        theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none", strip.background=element_blank(), strip.text=element_text(face="bold"), plot.title=element_text(size=10, face="bold"), panel.spacing=grid::unit(2.5, "mm"))

    width_in <- axis_widths[[axis]]/25.4
    height_in <- 108/25.4
    ragg::agg_png(file.path(out, paste0("discovery_ecotype_abundance_by_", axis, ".png")), width=width_in, height=height_in, units="in", res=600, background="white")
    print(p)
    dev.off()
    grDevices::cairo_pdf(file.path(out, paste0("discovery_ecotype_abundance_by_", axis, ".pdf")), width=width_in, height=height_in, family="sans")
    print(p)
    dev.off()
    svglite::svglite(file.path(out, paste0("discovery_ecotype_abundance_by_", axis, ".svg")), width=width_in, height=height_in, bg="white")
    print(p)
    dev.off()
}

# Step7. 合并 E1-E8 的 discovery bulk 与 pooled scRNA 主证据
evidence <- merge(assignment_summary, pooled_scrna[, c("Ecotype", "Z", "P_normal_greater", "P_empirical_greater", "Q_BH8_normal", "Q_BH8_empirical", "Cell2021Nominal", "EmpiricalNominal")], by="Ecotype", all=TRUE, sort=FALSE)
evidence <- evidence[match(ecotypes, evidence$Ecotype), ]
fwrite(evidence, file.path(out, "ecotype_discovery_bulk_and_pooled_scRNA_evidence.tsv"), sep="\t")
print(evidence)

capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
writeLines(c(
    "Discovery bulk cohort: 2,025 samples and 919 patients; this is internal reconstruction and biological characterization, not independent validation.",
    "Official assignment counts retain the EcoTyper assignment threshold; abundance-first analyses use all samples and classify each sample by its largest continuous EcoType abundance.",
    "Disease, Inflammation and Region3 are tested and plotted separately. Fisher tests characterize highest-abundance classes; Wilcoxon/Kruskal tests characterize continuous abundance.",
    "Disease patient sensitivity uses one median abundance per patient. Inflammation sensitivity uses paired within-patient medians for patients observed in both phenotypes.",
    "External validation is supplied by the pooled four-atlas scRNA co-occurrence test and is joined only at the final evidence-table stage."
), file.path(out, "method_boundary.txt"))
