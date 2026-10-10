suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE72819 Stage 6
# Build EcoTyper bulk-recovery input files
#
# Frozen inputs:
#   Stage 3 canonical expression:
#       24,656 genes x 73 samples
#       author-provided RPKM
#
#   Stage 4 canonical metadata:
#       73 samples
#
#   Stage 5 final sample QC:
#       73 evaluated
#       73 retained
#       0 technical exclusions
#
# This stage ONLY:
#   - converts frozen canonical objects to EcoTyper text input
#   - verifies exact sample/gene alignment
#   - verifies text round-trip
#
# This stage DOES NOT:
#   - normalize
#   - log-transform
#   - batch-correct
#   - remove samples
#   - perform signature-gene compatibility filtering
#   - run EcoTyper recovery
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

expr_file <- file.path(
  root,
  "prepared/03_canonical_expression",
  "GSE72819_canonical_RPKM_24656x73.rds"
)

meta_file <- file.path(
  root,
  "prepared/04_sample_metadata",
  "GSE72819_canonical_sample_metadata_73.rds"
)

qc_file <- file.path(
  root,
  "prepared/05_sample_expression_qc",
  "GSE72819_sample_level_qc.tsv"
)

forensic_file <- file.path(
  root,
  "prepared/05_sample_expression_qc",
  "GSE72819_stage5b_flagged_sample_forensic_qc.tsv"
)

out_dir <- file.path(
  root,
  "prepared/06_ecotyper_recovery_input"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Helpers
# ============================================================

check <- function(condition, message) {

  if (!isTRUE(condition)) {
    stop(
      paste0("[FAIL] ", message),
      call. = FALSE
    )
  }

  cat("[PASS] ", message, "\n", sep = "")
}


to_logical <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  toupper(
    trimws(
      as.character(x)
    )
  ) == "TRUE"
}


same_chr <- function(x, y) {

  x <- as.character(x)
  y <- as.character(y)

  length(x) == length(y) &&
    all(
      !is.na(x) &
      !is.na(y) &
      x == y
    )
}


cat("============================================================\n")
cat("GSE72819 Stage 6: EcoTyper recovery input preparation\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load frozen canonical expression
# ============================================================

expr <- readRDS(expr_file)

check(
  is.matrix(expr),
  "Canonical expression is a matrix"
)

check(
  nrow(expr) == 24656,
  "Canonical expression contains 24,656 genes"
)

check(
  ncol(expr) == 73,
  "Canonical expression contains 73 samples"
)

check(
  uniqueN(rownames(expr)) == 24656,
  "Canonical gene symbols are unique"
)

check(
  uniqueN(colnames(expr)) == 73,
  "Canonical sample IDs are unique"
)

check(
  !anyNA(expr),
  "Canonical expression contains no NA"
)

check(
  all(is.finite(expr)),
  "Canonical expression contains only finite values"
)

check(
  all(expr >= 0),
  "Canonical expression contains no negative values"
)

sample_ids <- as.character(
  colnames(expr)
)

gene_symbols <- as.character(
  rownames(expr)
)

numeric_payload_before <- as.vector(expr)


# ============================================================
# 2. Validate sample IDs against EcoTyper-safe naming
# ============================================================

check(
  all(
    make.names(
      sample_ids
    ) == sample_ids
  ),
  "All GSM sample IDs are stable under R make.names"
)


# ============================================================
# 3. Load frozen Stage-4 metadata
# ============================================================

meta <- as.data.table(
  readRDS(meta_file)
)

required_meta <- c(
  "GSM",
  "PatientID",
  "Disease",
  "Tissue",
  "TissueDiagnosis",
  "Timepoint",
  "Treatment",
  "DoseArm",
  "ResponseWeek10",
  "ResponseEvaluable",
  "PriorAntiTNF",
  "TNF_IR",
  "Run",
  "Lane",
  "RecoveryEligible",
  "ResponseAnalysisEligible"
)

check(
  all(
    required_meta %in%
      names(meta)
  ),
  "Canonical metadata contains required recovery fields"
)

check(
  nrow(meta) == 73,
  "Canonical metadata contains 73 rows"
)

check(
  same_chr(
    meta$GSM,
    sample_ids
  ),
  "Metadata order exactly matches expression columns"
)

check(
  all(
    to_logical(
      meta$RecoveryEligible
    )
  ),
  "All 73 Stage-4 samples are recovery eligible"
)


# ============================================================
# 4. Revalidate Stage-5 final sample decision
# ============================================================

qc <- fread(qc_file)

required_qc <- c(
  "GSM",
  "DiagnosticFlag",
  "TechnicalExclusion",
  "Retained"
)

check(
  all(
    required_qc %in%
      names(qc)
  ),
  "Stage-5 QC table contains required fields"
)

check(
  nrow(qc) == 73,
  "Stage-5 QC contains 73 samples"
)

qc <- qc[
  match(
    sample_ids,
    GSM
  )
]

check(
  same_chr(
    qc$GSM,
    sample_ids
  ),
  "Stage-5 QC aligns exactly to canonical expression"
)

check(
  sum(
    to_logical(
      qc$TechnicalExclusion
    )
  ) == 0,
  "Stage-5 technical exclusions = 0"
)

check(
  sum(
    to_logical(
      qc$Retained
    )
  ) == 73,
  "Stage-5 retained samples = 73"
)

check(
  sum(
    to_logical(
      qc$DiagnosticFlag
    )
  ) == 3,
  "Stage-5 diagnostic flags = 3"
)


# ============================================================
# 5. Revalidate targeted forensic decision
# ============================================================

forensic <- fread(
  forensic_file
)

check(
  nrow(forensic) == 3,
  "Stage-5B forensic table contains 3 diagnostic samples"
)

check(
  "TechnicalConcern" %in%
    names(forensic),
  "Stage-5B forensic table contains TechnicalConcern"
)

check(
  sum(
    to_logical(
      forensic$TechnicalConcern
    )
  ) == 0,
  "Stage-5B targeted technical concerns = 0"
)


# ============================================================
# 6. Build EcoTyper expression input
#
# Required orientation:
#   first column = gene symbols
#   subsequent columns = samples
#
# NO numeric transformation occurs.
# ============================================================

expr_dt <- as.data.table(
  expr
)

setnames(
  expr_dt,
  sample_ids
)

expr_out_dt <- cbind(
  data.table(
    GENES = gene_symbols
  ),
  expr_dt
)

check(
  nrow(expr_out_dt) == 24656,
  "Recovery expression table contains 24,656 genes"
)

check(
  ncol(expr_out_dt) == 74,
  "Recovery expression table contains GENES + 73 samples"
)

check(
  same_chr(
    names(expr_out_dt)[-1],
    sample_ids
  ),
  "Recovery expression column order matches canonical matrix"
)


# ============================================================
# 7. Build EcoTyper annotation input
#
# EcoTyper requires an ID column matching expression
# sample column names.
# ============================================================

annotation <- data.table(
  ID =
    sample_ids,

  PatientID =
    as.character(
      meta$PatientID
    ),

  Disease =
    as.character(
      meta$Disease
    ),

  Tissue =
    as.character(
      meta$Tissue
    ),

  TissueDiagnosis =
    as.character(
      meta$TissueDiagnosis
    ),

  Timepoint =
    as.character(
      meta$Timepoint
    ),

  Treatment =
    as.character(
      meta$Treatment
    ),

  DoseArm =
    as.character(
      meta$DoseArm
    ),

  ResponseWeek10 =
    as.character(
      meta$ResponseWeek10
    ),

  ResponseWeek10Plot =
    fifelse(
      is.na(
        meta$ResponseWeek10
      ),
      "Not_evaluable",
      as.character(
        meta$ResponseWeek10
      )
    ),

  ResponseEvaluable =
    to_logical(
      meta$ResponseEvaluable
    ),

  PriorAntiTNF =
    as.character(
      meta$PriorAntiTNF
    ),

  TNF_IR =
    as.character(
      meta$TNF_IR
    ),

  Run =
    as.character(
      meta$Run
    ),

  Lane =
    as.character(
      meta$Lane
    ),

  ResponseAnalysisEligible =
    to_logical(
      meta$ResponseAnalysisEligible
    )
)


check(
  nrow(annotation) == 73,
  "Recovery annotation contains 73 rows"
)

check(
  uniqueN(annotation$ID) == 73,
  "Recovery annotation contains 73 unique IDs"
)

check(
  same_chr(
    annotation$ID,
    sample_ids
  ),
  "Annotation ID order exactly matches expression columns"
)


# ============================================================
# 8. Frozen biological counts
# ============================================================

check(
  annotation[
    ResponseWeek10 == "Response",
    .N
  ] == 12,
  "Response = 12"
)

check(
  annotation[
    ResponseWeek10 == "Non-response",
    .N
  ] == 58,
  "Non-response = 58"
)

check(
  sum(
    is.na(
      annotation$ResponseWeek10
    )
  ) == 3,
  "Week-10 outcome missing = 3"
)

check(
  sum(
    annotation$ResponseAnalysisEligible
  ) == 70,
  "Response-analysis eligible = 70"
)

check(
  annotation[
    PriorAntiTNF == "Prior anti-TNF",
    .N
  ] == 46,
  "Prior anti-TNF = 46"
)

check(
  annotation[
    PriorAntiTNF == "Naive",
    .N
  ] == 27,
  "Anti-TNF naive = 27"
)

check(
  annotation[
    Run == "R1510",
    .N
  ] == 31,
  "Run R1510 = 31"
)

check(
  annotation[
    Run == "R1518",
    .N
  ] == 42,
  "Run R1518 = 42"
)


# ============================================================
# 9. Write final EcoTyper input files
# ============================================================

expression_out <- file.path(
  out_dir,
  "GSE72819_bulk_recovery_expression.txt"
)

annotation_out <- file.path(
  out_dir,
  "GSE72819_bulk_recovery_annotation.txt"
)

fwrite(
  expr_out_dt,
  expression_out,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  annotation,
  annotation_out,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

check(
  file.exists(
    expression_out
  ) &&
    file.info(
      expression_out
    )$size > 0,
  "Recovery expression file written"
)

check(
  file.exists(
    annotation_out
  ) &&
    file.info(
      annotation_out
    )$size > 0,
  "Recovery annotation file written"
)


# ============================================================
# 10. Expression text round-trip validation
# ============================================================

expr_reload <- fread(
  expression_out,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  na.strings = "NA"
)

check(
  nrow(expr_reload) == 24656,
  "Expression text round-trip retains 24,656 genes"
)

check(
  ncol(expr_reload) == 74,
  "Expression text round-trip retains GENES + 73 samples"
)

check(
  names(expr_reload)[1] ==
    "GENES",
  "Expression first column is GENES"
)

check(
  same_chr(
    expr_reload$GENES,
    gene_symbols
  ),
  "Expression text gene order is exact"
)

check(
  same_chr(
    names(expr_reload)[-1],
    sample_ids
  ),
  "Expression text sample order is exact"
)


reload_matrix <- as.matrix(
  expr_reload[
    ,
    -1,
    with = FALSE
  ]
)

storage.mode(
  reload_matrix
) <- "double"

rownames(
  reload_matrix
) <- as.character(
  expr_reload$GENES
)

colnames(
  reload_matrix
) <- sample_ids


check(
  all(
    dim(reload_matrix) ==
      dim(expr)
  ),
  "Expression text round-trip dimensions are exact"
)


max_abs_diff <- max(
  abs(
    reload_matrix -
      expr
  )
)

check(
  is.finite(
    max_abs_diff
  ),
  "Expression round-trip numeric difference is finite"
)

check(
  isTRUE(
    all.equal(
      reload_matrix,
      expr,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  ),
  "Expression text round-trip preserves RPKM values"
)


# ============================================================
# 11. Annotation text round-trip validation
# ============================================================

ann_reload <- fread(
  annotation_out,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  na.strings = "NA"
)

check(
  nrow(ann_reload) == 73,
  "Annotation text round-trip retains 73 samples"
)

check(
  "ID" %in%
    names(ann_reload),
  "Annotation text contains required ID column"
)

check(
  same_chr(
    ann_reload$ID,
    sample_ids
  ),
  "Annotation text ID order exactly matches expression columns"
)


# ============================================================
# 12. Verify canonical object was never modified
# ============================================================

numeric_payload_after <- as.vector(
  expr
)

check(
  identical(
    numeric_payload_before,
    numeric_payload_after
  ),
  "Frozen canonical RPKM remained unchanged"
)


# ============================================================
# 13. Input audit summary
# ============================================================

audit <- data.table(
  Metric = c(
    "ExpressionScale",
    "Genes",
    "Samples",
    "ExpressionRows",
    "ExpressionColumns",
    "AnnotationRows",
    "UniqueGeneSymbols",
    "UniqueSampleIDs",
    "MinExpression",
    "MaxExpression",
    "NegativeValues",
    "NAValues",
    "NonfiniteValues",
    "Response",
    "NonResponse",
    "ResponseMissing",
    "ResponseAnalysisEligible",
    "RecoveryEligible",
    "DiagnosticFlagsRetained",
    "TechnicalExclusions",
    "Stage5BTechnicalConcerns",
    "ExpressionTextMaxAbsoluteRoundTripDifference",
    "FinalSignatureCompatibility"
  ),

  Value = c(
    "author-provided RPKM",
    as.character(nrow(expr)),
    as.character(ncol(expr)),
    as.character(nrow(expr_out_dt)),
    as.character(ncol(expr_out_dt)),
    as.character(nrow(annotation)),
    as.character(uniqueN(gene_symbols)),
    as.character(uniqueN(sample_ids)),
    as.character(min(expr)),
    as.character(max(expr)),
    as.character(sum(expr < 0)),
    as.character(sum(is.na(expr))),
    as.character(sum(!is.finite(expr))),
    "12",
    "58",
    "3",
    "70",
    "73",
    "3",
    "0",
    "0",
    format(
      max_abs_diff,
      scientific = TRUE,
      digits = 16
    ),
    "PENDING_MSCCR_DISCOVERY_COMPLETION"
  )
)


cat("\n============================================================\n")
cat("ECOTYPER RECOVERY INPUT AUDIT\n")
cat("============================================================\n")

print(audit)


fwrite(
  audit,
  file.path(
    out_dir,
    "GSE72819_ecotyper_recovery_input_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 14. Sample manifest
# ============================================================

sample_manifest <- data.table(
  ExpressionColumn =
    seq_along(
      sample_ids
    ) + 1L,

  ID =
    sample_ids,

  PatientID =
    annotation$PatientID,

  ResponseWeek10 =
    annotation$ResponseWeek10,

  ResponseAnalysisEligible =
    annotation$ResponseAnalysisEligible,

  PriorAntiTNF =
    annotation$PriorAntiTNF,

  Run =
    annotation$Run,

  Lane =
    annotation$Lane,

  RecoveryEligible =
    TRUE,

  TechnicalExclusion =
    FALSE,

  DiagnosticFlag =
    to_logical(
      qc$DiagnosticFlag
    )
)


fwrite(
  sample_manifest,
  file.path(
    out_dir,
    "GSE72819_recovery_sample_manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


# ============================================================
# 15. Artifact hashes
# ============================================================

md5 <- tools::md5sum(
  c(
    expr_file,
    meta_file,
    qc_file,
    forensic_file,
    expression_out,
    annotation_out
  )
)

hash_table <- data.table(
  File =
    names(md5),

  MD5 =
    unname(md5)
)

fwrite(
  hash_table,
  file.path(
    out_dir,
    "GSE72819_stage6_artifact_md5.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 16. Human-readable status note
# ============================================================

status_lines <- c(
  "GSE72819 EcoTyper bulk recovery input",
  "",
  "Expression source: author-provided RPKM",
  "Canonical genes: 24656",
  "Recovery samples: 73",
  "Technical exclusions: 0",
  "Diagnostic flags retained: 3",
  "Response-analysis evaluable samples: 70",
  "",
  "No normalization performed.",
  "No log transformation performed.",
  "No batch correction performed.",
  "No sample removal performed.",
  "",
  "Final state-signature compatibility audit:",
  "PENDING MSCCR discovery completion."
)

writeLines(
  status_lines,
  file.path(
    out_dir,
    "GSE72819_recovery_input_README.txt"
  )
)


# ============================================================
# 17. Final structural gate
# ============================================================

critical_pass <- all(c(
  nrow(expr) == 24656,
  ncol(expr) == 73,
  uniqueN(gene_symbols) == 24656,
  uniqueN(sample_ids) == 73,

  nrow(annotation) == 73,

  same_chr(
    annotation$ID,
    sample_ids
  ),

  sum(
    to_logical(
      qc$TechnicalExclusion
    )
  ) == 0,

  sum(
    to_logical(
      qc$Retained
    )
  ) == 73,

  sum(
    to_logical(
      forensic$TechnicalConcern
    )
  ) == 0,

  all(is.finite(expr)),
  all(expr >= 0),

  isTRUE(
    all.equal(
      reload_matrix,
      expr,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  ),

  same_chr(
    ann_reload$ID,
    sample_ids
  ),

  identical(
    numeric_payload_before,
    numeric_payload_after
  )
))


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] 24,656 canonical genes exported\n")
  cat("[PASS] 73 recovery samples exported\n")
  cat("[PASS] Expression first column contains gene symbols\n")
  cat("[PASS] Expression sample columns exactly match canonical GSM order\n")
  cat("[PASS] Annotation contains required ID column\n")
  cat("[PASS] Annotation ID values exactly match expression columns\n")
  cat("[PASS] All 73 Stage-5 retained samples included\n")
  cat("[PASS] 3 diagnostic samples retained\n")
  cat("[PASS] Technical exclusions = 0\n")
  cat("[PASS] Expression scale remains author-provided RPKM\n")
  cat("[PASS] No normalization performed\n")
  cat("[PASS] No log transformation performed\n")
  cat("[PASS] No batch correction performed\n")
  cat("[PASS] No sample removal performed\n")
  cat("[PASS] Expression text round-trip validated\n")
  cat("[INFO] Final signature compatibility remains pending MSCCR discovery\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_ECOTYPER_RECOVERY_INPUTS_BUILT\n")

} else {

  cat("[REVIEW] Recovery input structural validation failed\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_ECOTYPER_RECOVERY_INPUTS\n")
}

cat("============================================================\n")
