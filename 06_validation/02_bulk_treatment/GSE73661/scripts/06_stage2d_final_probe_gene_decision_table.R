suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 2D v2
# FINAL transcript-cluster -> current-gene decision table
#
# Inputs:
#   Stage 1 authoritative log2-RMA expression
#   Stage 2C current-gene adjudication audit
#
# Principles:
#
# 1. Canonical recovery matrix is restricted to GPL6244
#    "main" design transcript clusters.
#
# 2. Retain only unique current-gene identities supported by:
#      A. concordant ChipDb + GPL6244 current GeneID
#      B. unique current ChipDb mapping without contradictory
#         current GPL6244 identity
#      C. GPL6244 unique disambiguation only when that current
#         GeneID is contained in the ChipDb multi-gene set
#
# 3. Exclude:
#      - non-main features
#      - no current gene identity
#      - unresolved multi-gene features
#      - current-gene conflicts
#      - malformed platform assignments
#      - platform disambiguation conflicts
#
# 4. NO expression collapse occurs in Stage 2D.
# ============================================================


# ============================================================
# Paths
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

expr_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE73661_authoritative_log2_RMA_33252x178.rds"
)

audit_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2c_current_gene_adjudication_audit.tsv"
)

out_dir <- file.path(
  root,
  "prepared/02_gene_harmonization"
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

  cat(
    "[PASS] ",
    message,
    "\n",
    sep = ""
  )
}


clean_chr <- function(x) {

  x <- trimws(
    as.character(x)
  )

  bad <- is.na(x) |
    x %in%
    c(
      "",
      "-",
      "---",
      "NA",
      "N/A",
      "NULL",
      "null"
    )

  x[bad] <- NA_character_

  x
}


split_set <- function(x) {

  x <- clean_chr(x)

  if (
    length(x) == 0L ||
    is.na(x[1])
  ) {
    return(character())
  }

  z <- strsplit(
    x[1],
    ";",
    fixed = TRUE
  )[[1]]

  z <- clean_chr(z)

  sort(
    unique(
      z[
        !is.na(z)
      ]
    )
  )
}


gene_in_set <- function(gene_id, gene_set) {

  gene_id <- clean_chr(
    gene_id
  )

  if (
    length(gene_id) == 0L ||
    is.na(gene_id[1])
  ) {
    return(FALSE)
  }

  gene_id[1] %in%
    split_set(
      gene_set
    )
}


cat("============================================================\n")
cat("GSE73661 Stage 2D v2: FINAL probe-to-current-gene decisions\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load frozen inputs
# ============================================================

check(
  file.exists(expr_file),
  "Frozen Stage-1 authoritative expression exists"
)

check(
  file.exists(audit_file),
  "Frozen Stage-2C adjudication audit exists"
)


expr <- readRDS(
  expr_file
)

audit <- fread(
  audit_file,
  showProgress = FALSE
)


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
  "Stage-2C audit contains 33,252 rows"
)


check(
  identical(
    as.character(
      audit$SourceProbeID
    ),
    as.character(
      rownames(expr)
    )
  ),
  "Stage-2C source order matches expression exactly"
)


# ============================================================
# 2. Verify frozen Stage-2C structure
# ============================================================

check(
  audit[
    PlatformCategory == "main",
    .N
  ] == 28869L,
  "28,869 main-design transcript clusters confirmed"
)


check(
  audit[
    PlatformCategory != "main",
    .N
  ] == 4383L,
  "4,383 non-main design features confirmed"
)


check(
  audit[
    AdjudicationClass ==
      "STRONG_CONCORDANT_CURRENT_GENE",
    .N
  ] == 18920L,
  "18,920 strong concordant mappings confirmed"
)


check(
  audit[
    AdjudicationClass ==
      "CANDIDATE_PLATFORM_DISAMBIGUATION",
    .N
  ] == 300L,
  "300 platform-disambiguation candidates confirmed"
)


check(
  audit[
    AdjudicationClass ==
      "REVIEW_UNIQUE_CURRENT_CONFLICT",
    .N
  ] == 18L,
  "18 unique current-gene conflicts confirmed"
)


check(
  audit[
    MalformedAssignments > 0L,
    .N
  ] == 5L,
  "Five malformed platform assignments confirmed"
)


# ============================================================
# 3. Platform-disambiguation containment audit
#
# A real disambiguation requires:
#
#   GPL6244 unique current GeneID
#              ∈
#   ChipDb current multi-gene candidate set
#
# Otherwise it is annotation conflict, not disambiguation.
# ============================================================

audit[
  ,
  DisambiguationContainedInChipSet :=
    as.logical(NA)
]


disambig_idx <- which(
  audit$AdjudicationClass ==
    "CANDIDATE_PLATFORM_DISAMBIGUATION"
)


check(
  length(disambig_idx) == 300L,
  "Exactly 300 disambiguation candidates selected"
)


audit[
  disambig_idx,
  DisambiguationContainedInChipSet :=
    mapply(
      gene_in_set,
      CandidateCurrentGeneID,
      ChipCurrentGeneIDs,
      USE.NAMES = FALSE
    )
]


disambig_audit <- audit[
  AdjudicationClass ==
    "CANDIDATE_PLATFORM_DISAMBIGUATION",
  .N,
  by = DisambiguationContainedInChipSet
][
  order(
    -N
  )
]


cat("============================================================\n")
cat("PLATFORM DISAMBIGUATION CONTAINMENT AUDIT\n")
cat("============================================================\n")

print(
  disambig_audit
)


check(
  sum(
    disambig_audit$N
  ) == 300L,
  "All 300 disambiguation candidates containment-audited"
)


n_disambig_true <- audit[
  AdjudicationClass ==
    "CANDIDATE_PLATFORM_DISAMBIGUATION" &
    DisambiguationContainedInChipSet == TRUE,
  .N
]


n_disambig_false <- audit[
  AdjudicationClass ==
    "CANDIDATE_PLATFORM_DISAMBIGUATION" &
    DisambiguationContainedInChipSet == FALSE,
  .N
]


cat(
  "[INFO] Valid within-ChipDb disambiguations: ",
  n_disambig_true,
  "\n",
  sep = ""
)

cat(
  "[INFO] Disambiguation conflicts: ",
  n_disambig_false,
  "\n",
  sep = ""
)


# ============================================================
# 4. Final source-feature decision
# ============================================================

audit[
  ,
  FinalDecision :=
    fcase(

      # ------------------------------------------------------
      # Gate platform category first
      # ------------------------------------------------------
      PlatformCategory != "main",
      "EXCLUDE_NON_MAIN_DESIGN_FEATURE",

      # ------------------------------------------------------
      # Official annotation parser anomaly
      # ------------------------------------------------------
      MalformedAssignments > 0L,
      "EXCLUDE_MALFORMED_PLATFORM_ASSIGNMENT",

      # ------------------------------------------------------
      # Strong independent agreement
      # ------------------------------------------------------
      AdjudicationClass ==
        "STRONG_CONCORDANT_CURRENT_GENE",
      "KEEP_STRONG_CONCORDANT_CURRENT_GENE",

      # ------------------------------------------------------
      # Supported unique ChipDb class, if ever present
      # ------------------------------------------------------
      AdjudicationClass ==
        "CHIPDB_CURRENT_SYMBOL_SUPPORTED",
      "KEEP_CHIPDB_CURRENT_SYMBOL_SUPPORTED",

      # ------------------------------------------------------
      # Unique current ChipDb identity without contradictory
      # current GPL6244 identity
      # ------------------------------------------------------
      AdjudicationClass ==
        "CHIPDB_CURRENT_ONLY",
      "KEEP_CHIPDB_CURRENT_ONLY",

      # ------------------------------------------------------
      # Valid platform disambiguation
      # ------------------------------------------------------
      AdjudicationClass ==
        "CANDIDATE_PLATFORM_DISAMBIGUATION" &
        DisambiguationContainedInChipSet == TRUE,
      "KEEP_PLATFORM_DISAMBIGUATION_WITHIN_CHIP_SET",

      # ------------------------------------------------------
      # Invalid "disambiguation" = annotation conflict
      # ------------------------------------------------------
      AdjudicationClass ==
        "CANDIDATE_PLATFORM_DISAMBIGUATION" &
        DisambiguationContainedInChipSet == FALSE,
      "EXCLUDE_PLATFORM_DISAMBIGUATION_CONFLICT",

      # ------------------------------------------------------
      # Unique current-gene conflict
      # ------------------------------------------------------
      AdjudicationClass %in%
        c(
          "REVIEW_UNIQUE_CURRENT_CONFLICT",
          "REVIEW_CHIPDB_SYMBOL_CONFLICT"
        ),
      "EXCLUDE_CURRENT_GENE_CONFLICT",

      # ------------------------------------------------------
      # Multiple current genes remain unresolved
      # ------------------------------------------------------
      AdjudicationClass ==
        "AMBIGUOUS_MULTI_CURRENT_EXACT_SET",
      "EXCLUDE_MULTI_GENE_AMBIGUOUS_EXACT_SET",

      AdjudicationClass ==
        "AMBIGUOUS_MULTI_CURRENT",
      "EXCLUDE_MULTI_GENE_AMBIGUOUS",

      # ------------------------------------------------------
      # No usable current gene identity
      # ------------------------------------------------------
      AdjudicationClass ==
        "NO_CURRENT_GENE_IDENTITY",
      "EXCLUDE_NO_CURRENT_GENE_IDENTITY",

      # ------------------------------------------------------
      # These are deliberately NOT auto-kept.
      # If they unexpectedly occur, keep as review.
      # ------------------------------------------------------
      AdjudicationClass ==
        "CANDIDATE_PLATFORM_ENTREZ_RESCUE",
      "REVIEW_PLATFORM_ENTREZ_RESCUE",

      AdjudicationClass ==
        "CANDIDATE_EXACT_CURRENT_SYMBOL_RESCUE",
      "REVIEW_EXACT_SYMBOL_RESCUE",

      default =
        "REVIEW_UNEXPECTED_CLASS"
    )
]


# ============================================================
# 5. Final inclusion flag
# ============================================================

audit[
  ,
  IncludeFinal :=
    startsWith(
      FinalDecision,
      "KEEP_"
    )
]


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
# 6. Final gene identity
# ============================================================

audit[
  ,
  FinalGeneID :=
    fifelse(
      IncludeFinal == TRUE,
      clean_chr(
        CandidateCurrentGeneID
      ),
      NA_character_
    )
]


audit[
  ,
  FinalSymbol :=
    fifelse(
      IncludeFinal == TRUE,
      clean_chr(
        CandidateCurrentSymbol
      ),
      NA_character_
    )
]


# ============================================================
# 7. Retained-source guards
# ============================================================

check(
  !any(
    audit[
      IncludeFinal == TRUE,
      is.na(
        FinalGeneID
      )
    ]
  ),
  "Every retained transcript cluster has FinalGeneID"
)


check(
  !any(
    audit[
      IncludeFinal == TRUE,
      is.na(
        FinalSymbol
      )
    ]
  ),
  "Every retained transcript cluster has FinalSymbol"
)


check(
  all(
    grepl(
      "^[0-9]+$",
      audit[
        IncludeFinal == TRUE,
        FinalGeneID
      ]
    )
  ),
  "Every retained FinalGeneID is one numeric GeneID"
)


check(
  !any(
    grepl(
      ";",
      audit[
        IncludeFinal == TRUE,
        FinalGeneID
      ],
      fixed = TRUE
    )
  ),
  "No retained feature maps to multiple final GeneIDs"
)


# ============================================================
# 8. Forbidden classes must never be retained
# ============================================================

check(
  audit[
    PlatformCategory != "main" &
      IncludeFinal == TRUE,
    .N
  ] == 0L,
  "No non-main design feature retained"
)


check(
  audit[
    MalformedAssignments > 0L &
      IncludeFinal == TRUE,
    .N
  ] == 0L,
  "No malformed platform assignment retained"
)


check(
  audit[
    AdjudicationClass %in%
      c(
        "AMBIGUOUS_MULTI_CURRENT",
        "AMBIGUOUS_MULTI_CURRENT_EXACT_SET",
        "REVIEW_UNIQUE_CURRENT_CONFLICT",
        "REVIEW_CHIPDB_SYMBOL_CONFLICT",
        "NO_CURRENT_GENE_IDENTITY"
      ) &
      IncludeFinal == TRUE,
    .N
  ] == 0L,
  "No ambiguous/conflicting/no-identity feature retained"
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


# ============================================================
# 9. Final decision summary
# ============================================================

decision_summary <- audit[
  ,
  .N,
  by = FinalDecision
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("FINAL PROBE-TO-GENE DECISION SUMMARY\n")
cat("============================================================\n")

print(
  decision_summary
)


review_rows <- audit[
  startsWith(
    FinalDecision,
    "REVIEW_"
  )
]


cat(
  "\nRows remaining for REVIEW: ",
  nrow(review_rows),
  "\n",
  sep = ""
)


# ============================================================
# 10. Final retained source features
# ============================================================

retained <- audit[
  IncludeFinal == TRUE
]


cat("\n============================================================\n")
cat("FINAL RETAINED SOURCE FEATURES\n")
cat("============================================================\n")

cat(
  "Retained transcript clusters: ",
  nrow(retained),
  "\n",
  sep = ""
)

cat(
  "Unique FinalGeneIDs: ",
  uniqueN(
    retained$FinalGeneID
  ),
  "\n",
  sep = ""
)

cat(
  "Unique FinalSymbols: ",
  uniqueN(
    retained$FinalSymbol
  ),
  "\n",
  sep = ""
)


# ============================================================
# 11. GeneID <-> symbol consistency
# ============================================================

gene_symbol_check <- retained[
  ,
  .(
    N_FinalSymbols =
      uniqueN(
        FinalSymbol
      ),

    Symbols =
      paste(
        sort(
          unique(
            FinalSymbol
          )
        ),
        collapse = ";"
      )
  ),
  by = FinalGeneID
]


check(
  all(
    gene_symbol_check$N_FinalSymbols ==
      1L
  ),
  "Every retained FinalGeneID maps to one FinalSymbol"
)


symbol_gene_check <- retained[
  ,
  .(
    N_FinalGeneIDs =
      uniqueN(
        FinalGeneID
      ),

    GeneIDs =
      paste(
        sort(
          unique(
            FinalGeneID
          )
        ),
        collapse = ";"
      )
  ),
  by = FinalSymbol
]


symbol_collision <- symbol_gene_check[
  N_FinalGeneIDs > 1L
]


check(
  nrow(symbol_collision) == 0L,
  "Every retained FinalSymbol maps to one FinalGeneID"
)


# ============================================================
# 12. Multi-transcript-cluster collision structure
#
# This defines Stage 2E workload.
# Expression is STILL untouched.
# ============================================================

gene_probe_counts <- retained[
  ,
  .(
    N_SourceTranscriptClusters =
      .N,

    SourceProbeIDs =
      paste(
        SourceProbeID,
        collapse = ";"
      ),

    MappingClasses =
      paste(
        sort(
          unique(
            FinalDecision
          )
        ),
        collapse = ";"
      )
  ),
  by = .(
    FinalGeneID,
    FinalSymbol
  )
][
  order(
    -N_SourceTranscriptClusters,
    FinalSymbol
  )
]


collision_summary <- data.table(
  Metric = c(
    "Retained_source_transcript_clusters",
    "Unique_final_genes",
    "Genes_with_one_source_transcript_cluster",
    "Genes_with_multiple_source_transcript_clusters",
    "Maximum_source_transcript_clusters_per_gene"
  ),

  N = c(
    nrow(retained),

    nrow(
      gene_probe_counts
    ),

    sum(
      gene_probe_counts$N_SourceTranscriptClusters ==
        1L
    ),

    sum(
      gene_probe_counts$N_SourceTranscriptClusters >
        1L
    ),

    if (
      nrow(
        gene_probe_counts
      ) > 0L
    ) {
      max(
        gene_probe_counts$N_SourceTranscriptClusters
      )
    } else {
      0L
    }
  )
)


cat("\n============================================================\n")
cat("FINAL MULTI-TRANSCRIPT-CLUSTER COLLISION STRUCTURE\n")
cat("============================================================\n")

print(
  collision_summary
)


cat("\nTranscript-cluster count distribution per final gene:\n")

collision_distribution <- gene_probe_counts[
  ,
  .N,
  by = N_SourceTranscriptClusters
][
  order(
    N_SourceTranscriptClusters
  )
]


print(
  collision_distribution
)


# ============================================================
# 13. Decision accounting guard
# ============================================================

check(
  sum(
    decision_summary$N
  ) == 33252L,
  "Final decisions account for all 33,252 source features"
)


check(
  nrow(retained) ==
    sum(
      audit$IncludeFinal == TRUE
    ),
  "Retained-row accounting is exact"
)


check(
  nrow(review_rows) ==
    sum(
      startsWith(
        audit$FinalDecision,
        "REVIEW_"
      )
    ),
  "Review-row accounting is exact"
)


# ============================================================
# 14. Save FINAL artifacts
# ============================================================

final_file <- file.path(
  out_dir,
  "GSE73661_canonical_probe_gene_decision_table_FINAL.tsv"
)


fwrite(
  audit,
  final_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  decision_summary,
  file.path(
    out_dir,
    "GSE73661_stage2d_final_decision_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  disambig_audit,
  file.path(
    out_dir,
    "GSE73661_stage2d_platform_disambiguation_containment.tsv"
  ),
  sep = "\t",
  na = "NA"
)


fwrite(
  retained[
    ,
    .(
      SourceRow,
      SourceProbeID,
      FinalGeneID,
      FinalSymbol,
      FinalDecision
    )
  ],
  file.path(
    out_dir,
    "GSE73661_stage2d_retained_probe_gene_mapping.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  review_rows,
  file.path(
    out_dir,
    "GSE73661_stage2d_remaining_review_rows.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  gene_probe_counts,
  file.path(
    out_dir,
    "GSE73661_stage2d_final_gene_probe_counts.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  collision_summary,
  file.path(
    out_dir,
    "GSE73661_stage2d_final_collision_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  collision_distribution,
  file.path(
    out_dir,
    "GSE73661_stage2d_final_collision_distribution.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 15. Artifact hashes
# ============================================================

hashes <- tools::md5sum(
  c(
    expr_file,
    audit_file,
    final_file
  )
)


hash_table <- data.table(
  File =
    names(hashes),

  MD5 =
    as.character(
      unname(
        hashes
      )
    )
)


fwrite(
  hash_table,
  file.path(
    out_dir,
    "GSE73661_stage2d_artifact_md5.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 16. Stage summary
# ============================================================

stage_summary <- data.table(
  Metric = c(
    "Source_features",
    "Main_design_features",
    "Non_main_features",
    "Strong_concordant",
    "Platform_disambiguation_candidates",
    "Valid_platform_disambiguations",
    "Invalid_platform_disambiguations",
    "Retained_source_features",
    "Unique_final_genes",
    "Genes_with_multiple_source_features",
    "Remaining_review_rows"
  ),

  N = c(
    nrow(audit),

    audit[
      PlatformCategory == "main",
      .N
    ],

    audit[
      PlatformCategory != "main",
      .N
    ],

    audit[
      AdjudicationClass ==
        "STRONG_CONCORDANT_CURRENT_GENE",
      .N
    ],

    length(
      disambig_idx
    ),

    n_disambig_true,

    n_disambig_false,

    nrow(
      retained
    ),

    nrow(
      gene_probe_counts
    ),

    gene_probe_counts[
      N_SourceTranscriptClusters >
        1L,
      .N
    ],

    nrow(
      review_rows
    )
  )
)


cat("\n============================================================\n")
cat("STAGE 2D SUMMARY\n")
cat("============================================================\n")

print(
  stage_summary
)


fwrite(
  stage_summary,
  file.path(
    out_dir,
    "GSE73661_stage2d_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 17. Final gate
# ============================================================

critical_pass <- all(
  c(
    nrow(audit) ==
      33252L,

    identical(
      audit$SourceProbeID,
      as.character(
        rownames(expr)
      )
    ),

    nrow(review_rows) ==
      0L,

    !anyNA(
      retained$FinalGeneID
    ),

    !anyNA(
      retained$FinalSymbol
    ),

    nrow(symbol_collision) ==
      0L,

    all(
      gene_symbol_check$N_FinalSymbols ==
        1L
    ),

    audit[
      PlatformCategory != "main" &
        IncludeFinal == TRUE,
      .N
    ] ==
      0L,

    audit[
      MalformedAssignments > 0L &
        IncludeFinal == TRUE,
      .N
    ] ==
      0L,

    audit[
      FinalDecision ==
        "EXCLUDE_PLATFORM_DISAMBIGUATION_CONFLICT" &
        IncludeFinal == TRUE,
      .N
    ] ==
      0L,

    sum(
      decision_summary$N
    ) ==
      33252L
  )
)


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] Final decisions assigned to all 33,252 source features\n")
  cat("[PASS] Canonical mapping restricted to GPL6244 main-design transcript clusters\n")
  cat("[PASS] Platform disambiguations containment-checked\n")
  cat("[PASS] Failed platform disambiguations excluded\n")
  cat("[PASS] No multi-gene ambiguous feature retained\n")
  cat("[PASS] No current-gene conflict retained\n")
  cat("[PASS] No malformed platform assignment retained\n")
  cat("[PASS] Every retained feature has one current GeneID and one symbol\n")
  cat("[PASS] Final GeneID-symbol mapping is internally unique\n")
  cat("[PASS] Multi-transcript-cluster collision structure quantified\n")
  cat("[PASS] No expression values modified\n")
  cat("[PASS] No transcript-cluster collapse performed yet\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE73661_STAGE2D_FINAL_PROBE_GENE_DECISION_TABLE\n")

} else {

  cat("[REVIEW] Stage 2D has unresolved rows or failed guards\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE73661_STAGE2D_FINAL_PROBE_GENE_DECISION_TABLE\n")
}

cat("============================================================\n")
