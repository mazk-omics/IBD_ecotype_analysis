suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# ============================================================
# GSE72819 Stage 1
# Source expression audit + authoritative RPKM matrix
#
# Goals:
#   1. Audit all 73 processed GEO files
#   2. Confirm identical gene universe/order
#   3. Confirm count / width / RPKM validity
#   4. Parse GeneID:<EntrezID>
#   5. Align expression, technical manifest and clinical metadata
#   6. Audit sequencing run vs treatment response
#   7. Build authoritative RPKM matrix WITHOUT transformation
#
# IMPORTANT:
#   - No normalization
#   - No log transform
#   - No batch correction
#   - No gene-ID remapping yet
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

raw_dir <- file.path(
  root,
  "02_source_expression/raw_tar_contents"
)

meta_file <- file.path(
  root,
  "01_source_metadata/original/metadata_geo_GSE72819.xlsx"
)

prov_manifest_file <- file.path(
  root,
  "00_provenance/GSE72819_RAW_processed_file_manifest.tsv"
)

out_dir <- file.path(
  root,
  "prepared/01_authoritative_matrix"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE72819 Stage 1 source expression audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Input files / metadata
# ============================================================

files <- list.files(
  raw_dir,
  pattern = "\\.counts_gene_exonic\\.tab\\.gz$",
  full.names = TRUE
)

files <- sort(files)

stopifnot(
  length(files) == 73
)

file_manifest <- fread(
  prov_manifest_file
)

stopifnot(
  nrow(file_manifest) == 73,
  uniqueN(file_manifest$GSM) == 73
)

meta <- as.data.table(
  read_excel(
    meta_file,
    sheet = "GEO数据整理"
  )
)

stopifnot(
  nrow(meta) == 73
)

# Standardize metadata GSM immediately
meta[
  ,
  GSM_clean := trimws(
    as.character(GSM编号)
  )
]

stopifnot(
  uniqueN(meta$GSM_clean) == 73,
  !anyNA(meta$GSM_clean),
  !anyDuplicated(meta$GSM_clean)
)


# ============================================================
# 2. Reader
# ============================================================

read_expr_file <- function(f) {

  x <- fread(
    f,
    sep = "\t",
    header = TRUE,
    data.table = TRUE
  )

  expected_names <- c(
    "name",
    "count",
    "width",
    "rpkm"
  )

  if (!identical(names(x), expected_names)) {

    stop(
      "Unexpected schema in ",
      basename(f),
      "\nObserved: ",
      paste(names(x), collapse = " | "),
      "\nExpected: ",
      paste(expected_names, collapse = " | ")
    )
  }

  x
}


# ============================================================
# 3. Reference source file
# ============================================================

ref <- read_expr_file(
  files[1]
)

ref_gene <- as.character(
  ref$name
)

ref_width <- ref$width

cat(
  "[INFO] Reference file:",
  basename(files[1]),
  "\n"
)

cat(
  "[INFO] Reference rows:",
  nrow(ref),
  "\n"
)


# ============================================================
# 4. Audit every processed file
# ============================================================

audit_list <- vector(
  "list",
  length(files)
)

count_list <- vector(
  "list",
  length(files)
)

rpkm_list <- vector(
  "list",
  length(files)
)


for (i in seq_along(files)) {

  f <- files[i]

  x <- read_expr_file(f)

  gsm <- sub(
    "^(GSM[0-9]+)_.*$",
    "\\1",
    basename(f)
  )

  same_gene_order <- identical(
    as.character(x$name),
    ref_gene
  )

  same_width <- identical(
    x$width,
    ref_width
  )

  numeric_ok <- (
    all(is.finite(x$count)) &&
    all(is.finite(x$width)) &&
    all(is.finite(x$rpkm))
  )

  audit_list[[i]] <- data.table(
    GSM = gsm,

    File =
      basename(f),

    Rows =
      nrow(x),

    UniqueGeneID =
      uniqueN(x$name),

    DuplicateGeneID_N =
      sum(duplicated(x$name)),

    SameGeneOrderAsReference =
      same_gene_order,

    SameWidthAsReference =
      same_width,

    NumericFinite =
      numeric_ok,

    CountMin =
      min(x$count),

    CountMax =
      max(x$count),

    CountSum =
      sum(x$count),

    CountFractional_N =
      sum(
        abs(
          x$count -
            round(x$count)
        ) > 1e-8
      ),

    WidthMin =
      min(x$width),

    WidthMax =
      max(x$width),

    WidthNonPositive_N =
      sum(x$width <= 0),

    RPKMMin =
      min(x$rpkm),

    RPKMMax =
      max(x$rpkm),

    RPKMZero_N =
      sum(x$rpkm == 0),

    RPKMNegative_N =
      sum(x$rpkm < 0),

    RPKMFractional_N =
      sum(
        abs(
          x$rpkm -
            round(x$rpkm)
        ) > 1e-8
      ),

    DetectedGenes_Count =
      sum(x$count > 0),

    DetectedGenes_RPKM =
      sum(x$rpkm > 0)
  )

  if (!same_gene_order) {

    stop(
      "Gene order mismatch in ",
      basename(f)
    )
  }

  if (!same_width) {

    stop(
      "Gene width mismatch in ",
      basename(f)
    )
  }

  count_list[[i]] <- as.numeric(
    x$count
  )

  rpkm_list[[i]] <- as.numeric(
    x$rpkm
  )
}


audit <- rbindlist(
  audit_list
)


# ============================================================
# 5. Dataset-wide source structure
# ============================================================

cat("\n============================================================\n")
cat("SOURCE FILE STRUCTURE\n")
cat("============================================================\n")

cat(
  "Files:",
  nrow(audit),
  "\n"
)

cat("\nRows per file:\n")

print(
  audit[
    ,
    .N,
    by = Rows
  ]
)

cat("\nUnique gene-ID counts:\n")

print(
  audit[
    ,
    .N,
    by = UniqueGeneID
  ]
)

cat("\nGene-order concordance:\n")

print(
  audit[
    ,
    .N,
    by = SameGeneOrderAsReference
  ]
)

cat("\nWidth concordance:\n")

print(
  audit[
    ,
    .N,
    by = SameWidthAsReference
  ]
)


# ============================================================
# 6. Gene identifier structure
#
# Source format is expected to be:
#   GeneID:1
#   GeneID:10
#   GeneID:100
#
# Strip only the explicit GeneID: prefix.
# No biological remapping occurs at this stage.
# ============================================================

has_geneid_prefix <- grepl(
  "^GeneID:",
  ref_gene
)

entrez_id <- sub(
  "^GeneID:",
  "",
  ref_gene
)

gene_id_valid <- grepl(
  "^[0-9]+$",
  entrez_id
)


gene_summary <- data.table(
  Metric = c(
    "Rows",
    "Unique_source_name",
    "GeneID_prefix_present",
    "Valid_Entrez_after_prefix_strip",
    "Invalid_after_prefix_strip",
    "Duplicated_source_name",
    "Duplicated_Entrez_after_prefix_strip"
  ),

  Value = c(
    length(ref_gene),
    uniqueN(ref_gene),
    sum(has_geneid_prefix),
    sum(gene_id_valid),
    sum(!gene_id_valid),
    sum(duplicated(ref_gene)),
    sum(duplicated(entrez_id))
  )
)


cat("\n============================================================\n")
cat("GENE IDENTIFIER STRUCTURE\n")
cat("============================================================\n")

print(
  gene_summary
)


if (any(!gene_id_valid)) {

  cat(
    "\nExamples of invalid IDs after GeneID: prefix stripping:\n"
  )

  print(
    head(
      ref_gene[
        !gene_id_valid
      ],
      30
    )
  )
}


# ============================================================
# 7. Build source matrices
#
# Values copied exactly from GEO processed files.
# ============================================================

gsm_order <- vapply(
  files,
  function(f) {

    sub(
      "^(GSM[0-9]+)_.*$",
      "\\1",
      basename(f)
    )

  },
  character(1)
)

gsm_order <- trimws(
  gsm_order
)


stopifnot(
  length(gsm_order) == 73,
  uniqueN(gsm_order) == 73,
  !anyDuplicated(gsm_order)
)


count_mat <- do.call(
  cbind,
  count_list
)

rpkm_mat <- do.call(
  cbind,
  rpkm_list
)


rownames(count_mat) <- ref_gene
rownames(rpkm_mat) <- ref_gene

colnames(count_mat) <- gsm_order
colnames(rpkm_mat) <- gsm_order

storage.mode(count_mat) <- "double"
storage.mode(rpkm_mat) <- "double"


stopifnot(
  ncol(count_mat) == 73,
  ncol(rpkm_mat) == 73,

  identical(
    rownames(count_mat),
    rownames(rpkm_mat)
  ),

  identical(
    colnames(count_mat),
    colnames(rpkm_mat)
  ),

  all(is.finite(count_mat)),
  all(is.finite(rpkm_mat))
)


# ============================================================
# 8. Metadata alignment
#
# IMPORTANT:
# Use value-wise character comparison rather than identical(),
# because readxl character attributes/encodings may differ.
# ============================================================

matrix_gsm <- trimws(
  as.character(
    colnames(rpkm_mat)
  )
)


meta_idx <- match(
  matrix_gsm,
  meta$GSM_clean
)


if (anyNA(meta_idx)) {

  missing_gsm <- matrix_gsm[
    is.na(meta_idx)
  ]

  cat(
    "\nMissing expression GSM in metadata:\n"
  )

  print(
    missing_gsm
  )

  stop(
    "Expression-to-metadata GSM matching failed"
  )
}


meta_aligned <- meta[
  meta_idx
]


alignment_ok <- all(
  as.character(
    meta_aligned$GSM_clean
  ) ==
    matrix_gsm
)


if (!alignment_ok) {

  bad <- which(
    as.character(
      meta_aligned$GSM_clean
    ) !=
      matrix_gsm
  )

  cat(
    "\nMetadata/expression GSM mismatch:\n"
  )

  print(
    data.table(
      Position =
        bad,

      Matrix_GSM =
        matrix_gsm[
          bad
        ],

      Metadata_GSM =
        meta_aligned$GSM_clean[
          bad
        ]
    )
  )

  stop(
    "Metadata GSM order does not match expression matrix"
  )
}


cat(
  "\n[PASS] Metadata aligned exactly to all 73 expression samples\n"
)


# ============================================================
# 9. Technical manifest alignment
# ============================================================

file_manifest[
  ,
  GSM := trimws(
    as.character(GSM)
  )
]


tech_idx <- match(
  matrix_gsm,
  file_manifest$GSM
)


if (anyNA(tech_idx)) {

  stop(
    "Expression-to-technical-manifest GSM matching failed"
  )
}


tech <- file_manifest[
  tech_idx
]


stopifnot(
  all(
    as.character(
      tech$GSM
    ) ==
      matrix_gsm
  )
)


cat(
  "[PASS] Technical manifest aligned exactly to expression samples\n"
)


# ============================================================
# 10. Response variables
# ============================================================

response <- fifelse(
  as.character(
    meta_aligned$缓解情况w10
  ) == "Remitter",

  "Response",

  fifelse(
    as.character(
      meta_aligned$缓解情况w10
    ) == "Non-remitter",

    "Non-response",

    NA_character_
  )
)


stopifnot(
  sum(
    response == "Response",
    na.rm = TRUE
  ) == 12,

  sum(
    response == "Non-response",
    na.rm = TRUE
  ) == 58,

  sum(
    is.na(response)
  ) == 3
)


# ============================================================
# 11. Unified sample manifest
# ============================================================

sample_manifest <- data.table(
  GSM =
    matrix_gsm,

  Run =
    as.character(
      tech$Run
    ),

  Lane =
    as.character(
      tech$Lane
    ),

  Library =
    as.character(
      tech$Library
    ),

  SAM =
    as.character(
      tech$SAM
    ),

  NXG =
    as.character(
      tech$NXG
    ),

  PatientID =
    as.character(
      meta_aligned$样本名
    ),

  Disease =
    as.character(
      meta_aligned$疾病
    ),

  Tissue =
    as.character(
      meta_aligned$组织
    ),

  TissueDiagnosis =
    as.character(
      meta_aligned$组织诊断
    ),

  ResponseWeek10 =
    response,

  PriorAntiTNF =
    as.character(
      meta_aligned$先前抗tnf治疗
    ),

  TNF_IR =
    as.character(
      meta_aligned$`tnf ir`
    )
)


stopifnot(
  nrow(sample_manifest) == 73,
  uniqueN(sample_manifest$GSM) == 73,
  uniqueN(sample_manifest$PatientID) == 73
)


# ============================================================
# 12. Run × clinical variables
# ============================================================

cat("\n============================================================\n")
cat("RUN x WEEK-10 RESPONSE\n")
cat("============================================================\n")


run_response <- dcast(
  sample_manifest[
    !is.na(ResponseWeek10)
  ],
  Run ~ ResponseWeek10,
  value.var = "GSM",
  fun.aggregate = length
)

print(
  run_response
)


cat("\n============================================================\n")
cat("RUN x PRIOR ANTI-TNF\n")
cat("============================================================\n")


run_anti_tnf <- dcast(
  sample_manifest,
  Run ~ PriorAntiTNF,
  value.var = "GSM",
  fun.aggregate = length
)

print(
  run_anti_tnf
)


cat("\n============================================================\n")
cat("RUN x TNF IR\n")
cat("============================================================\n")


run_tnf_ir <- dcast(
  sample_manifest,
  Run ~ TNF_IR,
  value.var = "GSM",
  fun.aggregate = length
)

print(
  run_tnf_ir
)


# ============================================================
# 13. Fisher exact test:
#     sequencing run vs Week-10 response
# ============================================================

response_tab <- table(
  sample_manifest$Run[
    !is.na(
      sample_manifest$ResponseWeek10
    )
  ],

  sample_manifest$ResponseWeek10[
    !is.na(
      sample_manifest$ResponseWeek10
    )
  ]
)


cat("\n============================================================\n")
cat("RUN x RESPONSE FISHER TEST\n")
cat("============================================================\n")

print(
  response_tab
)


if (
  nrow(response_tab) >= 2 &&
  ncol(response_tab) >= 2
) {

  fisher_response <- fisher.test(
    response_tab
  )

  print(
    fisher_response
  )

} else {

  fisher_response <- NULL

  cat(
    "Fisher test not applicable: insufficient table dimensions\n"
  )
}


# ============================================================
# 14. Per-sample expression QC
# ============================================================

sample_expr_qc <- data.table(
  GSM =
    matrix_gsm,

  CountSum =
    colSums(
      count_mat
    ),

  DetectedGenes_Count =
    colSums(
      count_mat > 0
    ),

  RPKM_Median =
    apply(
      rpkm_mat,
      2,
      median
    ),

  RPKM_Mean =
    colMeans(
      rpkm_mat
    ),

  RPKM_SD =
    apply(
      rpkm_mat,
      2,
      sd
    ),

  RPKM_IQR =
    apply(
      rpkm_mat,
      2,
      IQR
    ),

  DetectedGenes_RPKM =
    colSums(
      rpkm_mat > 0
    ),

  ZeroFraction_RPKM =
    colMeans(
      rpkm_mat == 0
    )
)


sample_expr_qc <- merge(
  sample_manifest,
  sample_expr_qc,
  by = "GSM",
  all.x = TRUE,
  sort = FALSE
)


# Restore expression order after merge
sample_expr_qc <- sample_expr_qc[
  match(
    matrix_gsm,
    GSM
  )
]


stopifnot(
  all(
    sample_expr_qc$GSM ==
      matrix_gsm
  )
)


# ============================================================
# 15. QC by sequencing run
# ============================================================

cat("\n============================================================\n")
cat("SAMPLE EXPRESSION QC BY RUN\n")
cat("============================================================\n")


run_qc <- sample_expr_qc[
  ,
  .(
    N = .N,

    CountSum_median =
      median(
        CountSum
      ),

    CountSum_min =
      min(
        CountSum
      ),

    CountSum_max =
      max(
        CountSum
      ),

    DetectedGenes_Count_median =
      median(
        DetectedGenes_Count
      ),

    RPKM_Median_median =
      median(
        RPKM_Median
      ),

    RPKM_Mean_median =
      median(
        RPKM_Mean
      ),

    RPKM_SD_median =
      median(
        RPKM_SD
      ),

    DetectedGenes_RPKM_median =
      median(
        DetectedGenes_RPKM
      ),

    ZeroFraction_RPKM_median =
      median(
        ZeroFraction_RPKM
      )
  ),

  by = Run
]

print(
  run_qc
)


# ============================================================
# 16. QC by response
# ============================================================

cat("\n============================================================\n")
cat("SAMPLE EXPRESSION QC BY RESPONSE\n")
cat("============================================================\n")


response_qc <- sample_expr_qc[
  !is.na(
    ResponseWeek10
  ),

  .(
    N = .N,

    CountSum_median =
      median(
        CountSum
      ),

    DetectedGenes_Count_median =
      median(
        DetectedGenes_Count
      ),

    RPKM_Median_median =
      median(
        RPKM_Median
      ),

    RPKM_Mean_median =
      median(
        RPKM_Mean
      ),

    RPKM_SD_median =
      median(
        RPKM_SD
      ),

    DetectedGenes_RPKM_median =
      median(
        DetectedGenes_RPKM
      ),

    ZeroFraction_RPKM_median =
      median(
        ZeroFraction_RPKM
      )
  ),

  by = ResponseWeek10
]

print(
  response_qc
)


# ============================================================
# 17. Source-expression summary
# ============================================================

summary <- data.table(
  Metric = c(
    "Files",
    "Samples",
    "Gene_rows",
    "Unique_source_gene_names",
    "Valid_Entrez_after_prefix_strip",
    "Unique_Entrez_after_prefix_strip",
    "All_files_same_gene_order",
    "All_files_same_width",
    "Any_nonfinite",
    "Any_negative_count",
    "Any_fractional_count",
    "Any_nonpositive_width",
    "Any_negative_RPKM",
    "Response_evaluable",
    "Response",
    "Non_response",
    "Missing_response"
  ),

  Value = c(
    length(files),

    ncol(rpkm_mat),

    nrow(rpkm_mat),

    uniqueN(
      rownames(rpkm_mat)
    ),

    sum(
      gene_id_valid
    ),

    uniqueN(
      entrez_id
    ),

    all(
      audit$SameGeneOrderAsReference
    ),

    all(
      audit$SameWidthAsReference
    ),

    any(
      !audit$NumericFinite
    ),

    any(
      audit$CountMin < 0
    ),

    any(
      audit$CountFractional_N > 0
    ),

    any(
      audit$WidthNonPositive_N > 0
    ),

    any(
      audit$RPKMNegative_N > 0
    ),

    sum(
      !is.na(response)
    ),

    sum(
      response == "Response",
      na.rm = TRUE
    ),

    sum(
      response == "Non-response",
      na.rm = TRUE
    ),

    sum(
      is.na(response)
    )
  )
)


cat("\n============================================================\n")
cat("AUTHORITATIVE SOURCE EXPRESSION SUMMARY\n")
cat("============================================================\n")

print(
  summary
)


# ============================================================
# 18. Save authoritative matrices
# ============================================================

rpkm_rds <- file.path(
  out_dir,
  paste0(
    "GSE72819_authoritative_RPKM_",
    nrow(rpkm_mat),
    "x73.rds"
  )
)

count_rds <- file.path(
  out_dir,
  paste0(
    "GSE72819_source_counts_",
    nrow(count_mat),
    "x73.rds"
  )
)


saveRDS(
  rpkm_mat,
  rpkm_rds,
  compress = FALSE
)

saveRDS(
  count_mat,
  count_rds,
  compress = FALSE
)


# ============================================================
# 19. Save audits
# ============================================================

fwrite(
  audit,
  file.path(
    out_dir,
    "GSE72819_source_file_expression_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_summary,
  file.path(
    out_dir,
    "GSE72819_source_gene_identifier_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  data.table(
    SourceName =
      ref_gene,

    EntrezID =
      entrez_id,

    ValidEntrez =
      gene_id_valid,

    Width =
      ref_width
  ),
  file.path(
    out_dir,
    "GSE72819_source_gene_manifest.tsv"
  ),
  sep = "\t"
)

fwrite(
  sample_manifest,
  file.path(
    out_dir,
    "GSE72819_source_sample_manifest.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  sample_expr_qc,
  file.path(
    out_dir,
    "GSE72819_source_sample_expression_qc.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  run_qc,
  file.path(
    out_dir,
    "GSE72819_expression_qc_by_run.tsv"
  ),
  sep = "\t"
)

fwrite(
  response_qc,
  file.path(
    out_dir,
    "GSE72819_expression_qc_by_response.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary,
  file.path(
    out_dir,
    "GSE72819_stage1_summary.tsv"
  ),
  sep = "\t"
)


if (!is.null(fisher_response)) {

  fisher_table <- data.table(
    Test =
      "Run_vs_ResponseWeek10",

    OddsRatio =
      unname(
        fisher_response$estimate
      ),

    PValue =
      fisher_response$p.value,

    ConfLow =
      fisher_response$conf.int[1],

    ConfHigh =
      fisher_response$conf.int[2]
  )

  fwrite(
    fisher_table,
    file.path(
      out_dir,
      "GSE72819_run_response_fisher_test.tsv"
    ),
    sep = "\t"
  )
}


# ============================================================
# 20. RDS round-trip
# ============================================================

rpkm_reload <- readRDS(
  rpkm_rds
)

count_reload <- readRDS(
  count_rds
)


stopifnot(
  identical(
    rpkm_reload,
    rpkm_mat
  ),

  identical(
    count_reload,
    count_mat
  )
)


cat(
  "\n[PASS] RDS round-trip exact for RPKM and count matrices\n"
)


# ============================================================
# 21. Final critical checks
# ============================================================

critical_pass <- (
  nrow(audit) == 73 &&

  length(
    unique(
      audit$Rows
    )
  ) == 1 &&

  all(
    audit$SameGeneOrderAsReference
  ) &&

  all(
    audit$SameWidthAsReference
  ) &&

  all(
    audit$NumericFinite
  ) &&

  all(
    audit$CountMin >= 0
  ) &&

  all(
    audit$CountFractional_N == 0
  ) &&

  all(
    audit$WidthNonPositive_N == 0
  ) &&

  all(
    audit$RPKMNegative_N == 0
  ) &&

  all(
    gene_id_valid
  ) &&

  uniqueN(
    entrez_id
  ) ==
    length(
      entrez_id
    ) &&

  nrow(rpkm_mat) == 26468 &&

  ncol(rpkm_mat) == 73 &&

  sum(
    !is.na(response)
  ) == 70 &&

  sum(
    response == "Response",
    na.rm = TRUE
  ) == 12 &&

  sum(
    response == "Non-response",
    na.rm = TRUE
  ) == 58
)


# ============================================================
# 22. Final status
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")


if (critical_pass) {

  cat("[PASS] 73 processed files audited\n")
  cat("[PASS] 26,468-gene source universe consistent across all files\n")
  cat("[PASS] Gene order identical across all files\n")
  cat("[PASS] Gene width identical across all files\n")
  cat("[PASS] Source gene names parse uniquely to Entrez GeneIDs\n")
  cat("[PASS] Counts are finite, non-negative integers\n")
  cat("[PASS] Width values are positive\n")
  cat("[PASS] RPKM values are finite and non-negative\n")
  cat("[PASS] Metadata aligned to all 73 expression samples\n")
  cat("[PASS] Technical manifest aligned to all 73 expression samples\n")
  cat("[PASS] Authoritative RPKM matrix built without transformation\n")
  cat("[PASS] Source count matrix retained without transformation\n")
  cat("[PASS] 70 response-evaluable samples retained\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE1_SOURCE_EXPRESSION_AUDIT\n")

} else {

  cat("[REVIEW] One or more source-expression checks require review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE1_SOURCE_EXPRESSION_AUDIT\n")
}

cat("============================================================\n")
