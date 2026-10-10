suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 2D FINAL VALIDATION
#
# Purpose:
#   Validate the Stage-2D FINAL decision table after correcting
#   a purely technical SourceProbeID type-comparison issue.
#
# IMPORTANT:
#   - no mapping decisions are changed
#   - no expression values are changed
#   - no transcript-cluster collapse is performed
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

expr_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE73661_authoritative_log2_RMA_33252x178.rds"
)

final_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_canonical_probe_gene_decision_table_FINAL.tsv"
)

summary_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2d_summary.tsv"
)

out_dir <- file.path(
  root,
  "prepared/02_gene_harmonization"
)


check <- function(condition, message) {

  if (!isTRUE(condition)) {
    stop(
      paste0("[FAIL] ", message),
      call. = FALSE
    )
  }

  cat(
    "[PASS] ",
    message,
    "\n",
    sep = ""
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 2D: FINAL VALIDATION\n")
cat("============================================================\n\n")


# ============================================================
# 1. Inputs
# ============================================================

check(
  file.exists(expr_file),
  "Frozen Stage-1 expression exists"
)

check(
  file.exists(final_file),
  "Stage-2D FINAL decision table exists"
)

check(
  file.exists(summary_file),
  "Stage-2D summary exists"
)


expr <- readRDS(
  expr_file
)

audit <- fread(
  final_file,
  showProgress = FALSE
)

stage_summary <- fread(
  summary_file,
  showProgress = FALSE
)


# ============================================================
# 2. Structural guards
# ============================================================

check(
  identical(
    dim(expr),
    c(
      33252L,
      178L
    )
  ),
  "Authoritative expression = 33,252 x 178"
)


check(
  nrow(audit) == 33252L,
  "FINAL decision table contains 33,252 source features"
)


# Critical correction:
# compare probe IDs as CHARACTER on both sides.
check(
  identical(
    as.character(
      audit$SourceProbeID
    ),
    as.character(
      rownames(expr)
    )
  ),
  "FINAL decision-table source IDs match expression order exactly"
)


# ============================================================
# 3. Required columns
# ============================================================

required_cols <- c(
  "SourceProbeID",
  "PlatformCategory",
  "MalformedAssignments",
  "AdjudicationClass",
  "DisambiguationContainedInChipSet",
  "FinalDecision",
  "IncludeFinal",
  "FinalGeneID",
  "FinalSymbol"
)


check(
  all(
    required_cols %in%
      names(audit)
  ),
  "All required FINAL decision columns are present"
)


# ============================================================
# 4. Reconstruct logical IncludeFinal safely
# ============================================================

if (!is.logical(audit$IncludeFinal)) {

  audit[
    ,
    IncludeFinal :=
      as.logical(
        IncludeFinal
      )
  ]
}


check(
  is.logical(
    audit$IncludeFinal
  ),
  "IncludeFinal is logical"
)


check(
  !anyNA(
    audit$IncludeFinal
  ),
  "IncludeFinal contains no NA"
)


# ============================================================
# 5. Exact expected Stage-2D decision counts
# ============================================================

expected <- data.table(
  FinalDecision = c(
    "KEEP_STRONG_CONCORDANT_CURRENT_GENE",
    "EXCLUDE_NO_CURRENT_GENE_IDENTITY",
    "EXCLUDE_NON_MAIN_DESIGN_FEATURE",
    "EXCLUDE_MULTI_GENE_AMBIGUOUS",
    "EXCLUDE_MULTI_GENE_AMBIGUOUS_EXACT_SET",
    "KEEP_PLATFORM_DISAMBIGUATION_WITHIN_CHIP_SET",
    "EXCLUDE_CURRENT_GENE_CONFLICT",
    "KEEP_CHIPDB_CURRENT_ONLY",
    "EXCLUDE_MALFORMED_PLATFORM_ASSIGNMENT",
    "EXCLUDE_PLATFORM_DISAMBIGUATION_CONFLICT"
  ),

  ExpectedN = c(
    18920L,
    7097L,
    4383L,
    1528L,
    988L,
    296L,
    18L,
    13L,
    5L,
    4L
  )
)


observed <- audit[
  ,
  .(
    ObservedN = .N
  ),
  by = FinalDecision
]


decision_check <- merge(
  expected,
  observed,
  by = "FinalDecision",
  all = TRUE,
  sort = FALSE
)


decision_check[
  is.na(ExpectedN),
  ExpectedN := 0L
]

decision_check[
  is.na(ObservedN),
  ObservedN := 0L
]


decision_check[
  ,
  Match :=
    ExpectedN ==
      ObservedN
]


cat("============================================================\n")
cat("FINAL DECISION COUNT VALIDATION\n")
cat("============================================================\n")

print(
  decision_check
)


check(
  all(
    decision_check$Match
  ),
  "All Stage-2D decision counts reproduce exactly"
)


check(
  sum(
    observed$ObservedN
  ) == 33252L,
  "Final decisions account for all 33,252 source features"
)


# ============================================================
# 6. Platform disambiguation containment
# ============================================================

disambig <- audit[
  AdjudicationClass ==
    "CANDIDATE_PLATFORM_DISAMBIGUATION"
]


check(
  nrow(disambig) == 300L,
  "Exactly 300 platform-disambiguation candidates recovered"
)


n_true <- disambig[
  DisambiguationContainedInChipSet == TRUE,
  .N
]

n_false <- disambig[
  DisambiguationContainedInChipSet == FALSE,
  .N
]


cat("\n============================================================\n")
cat("PLATFORM DISAMBIGUATION VALIDATION\n")
cat("============================================================\n")

cat(
  "Contained in ChipDb candidate set: ",
  n_true,
  "\n",
  sep = ""
)

cat(
  "Conflicting with ChipDb candidate set: ",
  n_false,
  "\n",
  sep = ""
)


check(
  n_true == 296L,
  "296 valid platform disambiguations confirmed"
)

check(
  n_false == 4L,
  "4 platform-disambiguation conflicts confirmed"
)


# ============================================================
# 7. Retained mappings
# ============================================================

retained <- audit[
  IncludeFinal == TRUE
]


check(
  nrow(retained) == 19229L,
  "19,229 source transcript clusters retained"
)


check(
  uniqueN(
    retained$FinalGeneID
  ) == 18422L,
  "18,422 unique FinalGeneIDs retained"
)


check(
  uniqueN(
    retained$FinalSymbol
  ) == 18422L,
  "18,422 unique FinalSymbols retained"
)


check(
  !anyNA(
    retained$FinalGeneID
  ),
  "No retained FinalGeneID is missing"
)


check(
  !anyNA(
    retained$FinalSymbol
  ),
  "No retained FinalSymbol is missing"
)


check(
  all(
    grepl(
      "^[0-9]+$",
      as.character(
        retained$FinalGeneID
      )
    )
  ),
  "Every retained FinalGeneID is one numeric GeneID"
)


check(
  !any(
    grepl(
      ";",
      as.character(
        retained$FinalGeneID
      ),
      fixed = TRUE
    )
  ),
  "No retained feature maps to multiple FinalGeneIDs"
)


# ============================================================
# 8. GeneID-symbol uniqueness
# ============================================================

gene_to_symbol <- retained[
  ,
  .(
    N_Symbols =
      uniqueN(
        FinalSymbol
      )
  ),
  by = FinalGeneID
]


check(
  all(
    gene_to_symbol$N_Symbols ==
      1L
  ),
  "Every FinalGeneID maps to one FinalSymbol"
)


symbol_to_gene <- retained[
  ,
  .(
    N_GeneIDs =
      uniqueN(
        FinalGeneID
      )
  ),
  by = FinalSymbol
]


check(
  all(
    symbol_to_gene$N_GeneIDs ==
      1L
  ),
  "Every FinalSymbol maps to one FinalGeneID"
)


# ============================================================
# 9. Forbidden-retention guards
# ============================================================

check(
  audit[
    PlatformCategory != "main" &
      IncludeFinal == TRUE,
    .N
  ] == 0L,
  "No non-main feature retained"
)


check(
  audit[
    MalformedAssignments > 0L &
      IncludeFinal == TRUE,
    .N
  ] == 0L,
  "No malformed platform feature retained"
)


check(
  audit[
    FinalDecision ==
      "EXCLUDE_CURRENT_GENE_CONFLICT" &
      IncludeFinal == TRUE,
    .N
  ] == 0L,
  "No current-gene conflict retained"
)


check(
  audit[
    FinalDecision ==
      "EXCLUDE_PLATFORM_DISAMBIGUATION_CONFLICT" &
      IncludeFinal == TRUE,
    .N
  ] == 0L,
  "No failed platform disambiguation retained"
)


check(
  audit[
    grepl(
      "^EXCLUDE_MULTI_GENE",
      FinalDecision
    ) &
      IncludeFinal == TRUE,
    .N
  ] == 0L,
  "No multi-gene ambiguous feature retained"
)


# ============================================================
# 10. REVIEW rows
# ============================================================

review_rows <- audit[
  grepl(
    "^REVIEW_",
    FinalDecision
  )
]


check(
  nrow(review_rows) == 0L,
  "No Stage-2D rows remain for review"
)


# ============================================================
# 11. Final gene collision structure
# ============================================================

gene_counts <- retained[
  ,
  .(
    N_SourceTranscriptClusters =
      .N
  ),
  by = .(
    FinalGeneID,
    FinalSymbol
  )
]


n_single <- gene_counts[
  N_SourceTranscriptClusters == 1L,
  .N
]

n_multi <- gene_counts[
  N_SourceTranscriptClusters > 1L,
  .N
]

max_clusters <- max(
  gene_counts$N_SourceTranscriptClusters
)


cat("\n============================================================\n")
cat("FINAL COLLISION VALIDATION\n")
cat("============================================================\n")

cat(
  "Genes with one transcript cluster: ",
  n_single,
  "\n",
  sep = ""
)

cat(
  "Genes with multiple transcript clusters: ",
  n_multi,
  "\n",
  sep = ""
)

cat(
  "Maximum transcript clusters per gene: ",
  max_clusters,
  "\n",
  sep = ""
)


check(
  n_single == 17804L,
  "17,804 genes have one source transcript cluster"
)

check(
  n_multi == 618L,
  "618 genes have multiple source transcript clusters"
)

check(
  max_clusters == 26L,
  "Maximum source transcript clusters per gene = 26"
)


# ============================================================
# 12. Validation artifact
# ============================================================

validation_summary <- data.table(
  Metric = c(
    "Source_features",
    "Retained_source_transcript_clusters",
    "Unique_final_genes",
    "Valid_platform_disambiguations",
    "Invalid_platform_disambiguations",
    "Genes_with_one_source_transcript_cluster",
    "Genes_with_multiple_source_transcript_clusters",
    "Maximum_source_transcript_clusters_per_gene",
    "Remaining_review_rows"
  ),

  N = c(
    33252L,
    nrow(retained),
    uniqueN(retained$FinalGeneID),
    n_true,
    n_false,
    n_single,
    n_multi,
    max_clusters,
    nrow(review_rows)
  )
)


fwrite(
  validation_summary,
  file.path(
    out_dir,
    "GSE73661_stage2d_FINAL_validation_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  decision_check,
  file.path(
    out_dir,
    "GSE73661_stage2d_FINAL_decision_count_validation.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 13. FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] Stage-2D FINAL table contains all 33,252 source features\n")
cat("[PASS] SourceProbeID order matches authoritative expression after explicit character coercion\n")
cat("[PASS] All final decision counts reproduce exactly\n")
cat("[PASS] 296 valid platform disambiguations confirmed\n")
cat("[PASS] 4 conflicting platform disambiguations excluded\n")
cat("[PASS] 19,229 transcript clusters retained\n")
cat("[PASS] 18,422 current genes retained\n")
cat("[PASS] No unresolved review rows remain\n")
cat("[PASS] GeneID-symbol identity is one-to-one\n")
cat("[PASS] 618 multi-transcript-cluster genes identified for Stage 2E\n")
cat("[PASS] No expression values modified\n")
cat("[PASS] No transcript-cluster collapse performed\n")

cat("\nFINAL STATUS:\n")
cat("PASS_GSE73661_STAGE2D_FINAL_PROBE_GENE_DECISION_TABLE\n")
cat("============================================================\n")
