# 本机单线程核对文献匹配统计，并对预定7候选做同患者同肠段炎症比较。
# 不改冻结丰度，病例不足10对不作推断；仅研究背景关联，不作为独立临床验证。
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
out <- "D:/IBD_ecotype/IBD_ecotype_analysis/analysis_reports/legacy_known_state_correspondence/run_20261013_01"

# Step1. 用R独立核对超几何P与所有可评价匹配的BH
matches <- read.delim(file.path(out, "audit/ALL_PUBLISHED_SIGNATURE_MATCH_STATISTICS.tsv"), check.names=FALSE)
eligible <- which(matches$Evaluable=="True")
q <- matches[eligible, ]
p_check <- phyper(q$OverlapN-1, q$SignatureInBackgroundN, q$BackgroundN-q$SignatureInBackgroundN, q$QueryInBackgroundN, lower.tail=FALSE)
q_check <- p.adjust(p_check, "BH")
qc <- data.frame(CellState=q$CellState, PublishedState=q$PublishedStateOrProgram, Query=q$Query, PError=abs(p_check-q$P), QError=abs(q_check-q$BH_Q_GlobalEvaluableMatches))
stopifnot(max(qc$PError)<1e-10, max(qc$QError)<1e-10)
write.table(qc, file.path(out, "audit/SIGNATURE_MATCH_INDEPENDENT_R_QC.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
print(c(Evaluable=nrow(qc), MaxPError=max(qc$PError), MaxQError=max(qc$QError)))

# Step2. 以患者×疾病×肠段配对，控制重复活检与主要解剖/诊断混杂
pairs <- read.delim(file.path(out, "audit/MSCCR_WITHIN_PATIENT_INFLAMMATION_PAIRS.tsv"), check.names=FALSE)
keys <- unique(pairs[, c("CellState", "Disease", "Region3")])
results <- list()
for (i in seq_len(nrow(keys))) {
  key <- keys[i, ]
  group <- pairs[pairs$CellState==key$CellState & pairs$Disease==key$Disease & pairs$Region3==key$Region3 & pairs$EvaluablePair=="True", ]
  delta <- group$WithinPatientDifference
  # <1e-10视为求解器数值零，只用于此敏感性统计，不回写原丰度。
  delta_tolerance <- delta
  delta_tolerance[abs(delta_tolerance)<1e-10] <- 0
  p <- NA_real_
  p_original <- NA_real_
  if (length(delta)>=10 && any(delta_tolerance!=0)) {
    p <- wilcox.test(delta_tolerance, mu=0, exact=FALSE, correct=TRUE)$p.value
    p_original <- wilcox.test(delta, mu=0, exact=FALSE, correct=TRUE)$p.value
  }
  results[[i]] <- data.frame(CellState=key$CellState, Disease=key$Disease, Region3=key$Region3, PairedPatients=length(delta), NonzeroPairsNumericTolerance=sum(delta_tolerance!=0), MedianInflamedMinusNonInflamed=if(length(delta)) median(delta) else NA_real_, FractionPatientDifferencesPositive=if(length(delta)) mean(delta_tolerance>0) else NA_real_, P=p, P_OriginalNoTolerance=p_original, Status=if(length(delta)>=10) "EVALUABLE_PATIENT_PAIRED" else "NOT_EVALUABLE_FEWER10_PATIENT_PAIRS", Method="WithinPatient_sameDisease_sameRegion; median repeated biopsies; two-sided Wilcoxon signed rank; numeric zero tolerance1e-10", EvidenceScope="Legacy discovery cohort association; not independent clinical validation; fractions within parent are compositional")
}
results <- do.call(rbind, results)
results$BH_Q_AllEvaluablePatientComparisons <- NA_real_
valid <- which(!is.na(results$P))
results$BH_Q_AllEvaluablePatientComparisons[valid] <- p.adjust(results$P[valid], "BH")
results$BH_Q_OriginalNoTolerance <- NA_real_
valid_original <- which(!is.na(results$P_OriginalNoTolerance))
results$BH_Q_OriginalNoTolerance[valid_original] <- p.adjust(results$P_OriginalNoTolerance[valid_original], "BH")
write.table(results, file.path(out, "audit/CANDIDATE_MSCCR_PATIENT_PAIRED_CONTEXT.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="NA")
print(results[, c("CellState", "Disease", "Region3", "PairedPatients", "MedianInflamedMinusNonInflamed", "BH_Q_AllEvaluablePatientComparisons")])
