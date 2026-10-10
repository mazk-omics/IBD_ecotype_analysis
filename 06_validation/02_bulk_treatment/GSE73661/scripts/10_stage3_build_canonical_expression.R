suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 3
# FINAL canonical gene-level expression construction
#
# Frozen upstream:
#   Stage 1  authoritative author-provided log2 RMA
#   Stage 2D current-gene identity FINAL/FROZEN
#   Stage 2E expression concordance audit PASS/FROZEN
#   Stage 2F annotation/design representative audit PASS/FROZEN
#   Stage 2G exact-duplicate tie audit PASS/FROZEN
#
# Final representative logic:
#
#   1. Single-cluster gene:
#        retain its only transcript cluster.
#
#   2. Multi-cluster gene resolved uniquely by Stage 2F:
#        retain Stage-2F annotation/design candidate.
#
#   3. Stage-2F tied gene whose tied candidates are EXACTLY
#      numerically identical across all 178 samples:
#        retain Stage-2G deterministic candidate.
#
#   4. Stage-2F tied gene with non-identical expression:
#        exclude from canonical matrix for now.
#
# IMPORTANT:
#   NO averaging
#   NO median collapse
#   NO MaxMean selection
#   NO MaxIQR selection
#   NO correlation-based selection
#   NO normalization
#   NO batch correction
# ============================================================


# ============================================================
# 0. Paths
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

expr_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE73661_authoritative_log2_RMA_33252x178.rds"
)

decision_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_canonical_probe_gene_decision_table_FINAL.tsv"
)

stage2f_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2f_annotation_design_candidates.tsv"
)

stage2g_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2g_exact_duplicate_tie_audit.tsv"
)

out_dir <- file.path(
  root,
  "prepared/03_canonical_expression"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. Helpers
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


cat("============================================================\n")
cat("GSE73661 Stage 3: build FINAL canonical expression\n")
cat("============================================================\n\n")


# ============================================================
# 2. Input existence
# ============================================================

for (f in c(
  expr_file,
  decision_file,
  stage2f_file,
  stage2g_file
)) {

  check(
    file.exists(f),
    paste0(
      "Input exists: ",
      basename(f)
    )
  )
}


# ============================================================
# 3. Load frozen upstream artifacts
# ============================================================

expr <- readRDS(
  expr_file
)

decision <- fread(
  decision_file,
  showProgress = FALSE
)

stage2f <- fread(
  stage2f_file,
  showProgress = FALSE
)

stage2g <- fread(
  stage2g_file,
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


# ============================================================
# 4. Freeze identifier types
# ============================================================

decision[
  ,
  SourceProbeID :=
    as.character(SourceProbeID)
]

decision[
  ,
  FinalGeneID :=
    as.character(FinalGeneID)
]

decision[
  ,
  FinalSymbol :=
    as.character(FinalSymbol)
]


stage2f[
  ,
  FinalGeneID :=
    as.character(FinalGeneID)
]

stage2f[
  ,
  FinalSymbol :=
    as.character(FinalSymbol)
]

stage2f[
  ,
  Candidate_AnnotationDesign :=
    as.character(Candidate_AnnotationDesign)
]


stage2g[
  ,
  FinalGeneID :=
    as.character(FinalGeneID)
]

stage2g[
  ,
  FinalSymbol :=
    as.character(FinalSymbol)
]

stage2g[
  ,
  DeterministicCandidateIfExact :=
    as.character(DeterministicCandidateIfExact)
]


if (!is.logical(decision$IncludeFinal)) {

  decision[
    ,
    IncludeFinal :=
      as.logical(IncludeFinal)
  ]
}


if (!is.logical(stage2g$AllTiedCandidatesExact)) {

  stage2g[
    ,
    AllTiedCandidatesExact :=
      as.logical(AllTiedCandidatesExact)
  ]
}


expr_probe_ids <- as.character(
  rownames(expr)
)


check(
  identical(
    decision$SourceProbeID,
    expr_probe_ids
  ),
  "Stage-2D source IDs match authoritative expression order"
)


# ============================================================
# 5. Recover frozen retained transcript clusters
# ============================================================

retained <- decision[
  IncludeFinal == TRUE
]


check(
  nrow(retained) == 19229L,
  "19,229 retained Stage-2D transcript clusters recovered"
)


check(
  uniqueN(
    retained$FinalGeneID
  ) == 18422L,
  "18,422 Stage-2D current genes recovered"
)


# ============================================================
# 6. Gene multiplicity
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


single <- gene_counts[
  N_SourceTranscriptClusters == 1L
]


multi <- gene_counts[
  N_SourceTranscriptClusters > 1L
]


check(
  nrow(single) == 17804L,
  "17,804 single-cluster genes confirmed"
)


check(
  nrow(multi) == 618L,
  "618 multi-cluster genes confirmed"
)


# ============================================================
# 7. Single-cluster representative mapping
# ============================================================

single_mapping <- retained[
  FinalGeneID %in%
    single$FinalGeneID,
  .(
    FinalGeneID,
    FinalSymbol,
    RepresentativeProbeID =
      SourceProbeID,
    RepresentativeRule =
      "SINGLE_RETAINED_TRANSCRIPT_CLUSTER"
  )
]


check(
  nrow(single_mapping) == 17804L,
  "17,804 single-cluster representatives assigned"
)


# ============================================================
# 8. Stage-2F uniquely resolved multi-cluster genes
# ============================================================

stage2f_unique <- stage2f[
  N_AnnotationDesignCandidates == 1L
]


check(
  nrow(stage2f_unique) == 498L,
  "498 annotation/design-unique multi-cluster genes recovered"
)


check(
  !anyNA(
    stage2f_unique$Candidate_AnnotationDesign
  ),
  "All 498 unique Stage-2F candidates have SourceProbeID"
)


stage2f_mapping <- stage2f_unique[
  ,
  .(
    FinalGeneID,
    FinalSymbol,
    RepresentativeProbeID =
      Candidate_AnnotationDesign,
    RepresentativeRule =
      "ANNOTATION_DESIGN_UNIQUE"
  )
]


# ============================================================
# 9. Stage-2G exact-numeric duplicate ties
# ============================================================

exact_ties <- stage2g[
  AllTiedCandidatesExact == TRUE
]


nonexact_ties <- stage2g[
  AllTiedCandidatesExact == FALSE
]


check(
  nrow(exact_ties) == 105L,
  "105 exact-numeric duplicate ties recovered"
)


check(
  nrow(nonexact_ties) == 15L,
  "15 non-identical unresolved ties recovered"
)


check(
  !anyNA(
    exact_ties$DeterministicCandidateIfExact
  ),
  "All exact ties have deterministic representative SourceProbeID"
)


check(
  all(
    is.na(
      nonexact_ties$DeterministicCandidateIfExact
    )
  ),
  "No unresolved non-exact tie has representative SourceProbeID"
)


exact_mapping <- exact_ties[
  ,
  .(
    FinalGeneID,
    FinalSymbol,
    RepresentativeProbeID =
      DeterministicCandidateIfExact,
    RepresentativeRule =
      "EXACT_NUMERIC_DUPLICATE_DETERMINISTIC"
  )
]


# ============================================================
# 10. Assemble FINAL representative mapping
# ============================================================

representative_mapping <- rbindlist(
  list(
    single_mapping,
    stage2f_mapping,
    exact_mapping
  ),
  use.names = TRUE
)


check(
  nrow(representative_mapping) == 18407L,
  "18,407 final representative genes assembled"
)


check(
  uniqueN(
    representative_mapping$FinalGeneID
  ) == 18407L,
  "Final representative GeneIDs are unique"
)


check(
  uniqueN(
    representative_mapping$FinalSymbol
  ) == 18407L,
  "Final representative symbols are unique"
)


check(
  uniqueN(
    representative_mapping$RepresentativeProbeID
  ) == 18407L,
  "Every final gene has a unique representative transcript cluster"
)


# ============================================================
# 11. Verify representative probes belong to correct gene
# ============================================================

probe_gene_lookup <- setNames(
  retained$FinalGeneID,
  retained$SourceProbeID
)


probe_symbol_lookup <- setNames(
  retained$FinalSymbol,
  retained$SourceProbeID
)


mapped_gene <- unname(
  probe_gene_lookup[
    representative_mapping$RepresentativeProbeID
  ]
)


mapped_symbol <- unname(
  probe_symbol_lookup[
    representative_mapping$RepresentativeProbeID
  ]
)


check(
  !anyNA(mapped_gene),
  "Every representative probe exists in retained Stage-2D mapping"
)


check(
  identical(
    mapped_gene,
    representative_mapping$FinalGeneID
  ),
  "Every representative probe belongs to its assigned FinalGeneID"
)


check(
  identical(
    mapped_symbol,
    representative_mapping$FinalSymbol
  ),
  "Every representative probe belongs to its assigned FinalSymbol"
)


# ============================================================
# 12. Ensure unresolved genes are excluded
# ============================================================

unresolved_gene_ids <- nonexact_ties$FinalGeneID


check(
  length(
    unresolved_gene_ids
  ) == 15L,
  "Exactly 15 unresolved gene IDs defined"
)


check(
  !any(
    representative_mapping$FinalGeneID %in%
      unresolved_gene_ids
  ),
  "No unresolved non-exact gene entered canonical mapping"
)


# ============================================================
# 13. Explicit unresolved-gene table
# ============================================================

unresolved <- nonexact_ties[
  ,
  .(
    FinalGeneID,
    FinalSymbol,
    N_TiedCandidates,
    TiedCandidates,
    MinSpearmanTied,
    MedianSpearmanTied,
    MinPearsonTied,
    MedianPairMeanAbsDiff,
    MaxPairMeanAbsDiff,
    MaxAbsoluteDifference,
    TieResolutionClass,
    FinalCanonicalDecision =
      "EXCLUDE_UNRESOLVED_NONEXACT_MULTICLUSTER"
  )
]


check(
  nrow(unresolved) == 15L,
  "Unresolved exclusion table contains 15 genes"
)


# ============================================================
# 14. Deterministic final gene order
#
# Preserve first occurrence in the frozen Stage-2D retained
# mapping rather than alphabetical reordering.
# ============================================================

first_gene_order <- retained[
  ,
  .(
    FirstSourceRow =
      min(SourceRow)
  ),
  by = .(
    FinalGeneID,
    FinalSymbol
  )
]


representative_mapping[
  first_gene_order,
  on = c(
    "FinalGeneID",
    "FinalSymbol"
  ),
  FirstSourceRow :=
    i.FirstSourceRow
]


check(
  !anyNA(
    representative_mapping$FirstSourceRow
  ),
  "Every canonical gene has original source-order anchor"
)


setorder(
  representative_mapping,
  FirstSourceRow
)


# ============================================================
# 15. Construct canonical expression matrix
# ============================================================

selected_idx <- match(
  representative_mapping$RepresentativeProbeID,
  expr_probe_ids
)


check(
  !anyNA(selected_idx),
  "Every representative transcript cluster found in authoritative expression"
)


canonical_expr <- expr[
  selected_idx,
  ,
  drop = FALSE
]


rownames(
  canonical_expr
) <- representative_mapping$FinalSymbol


check(
  identical(
    dim(canonical_expr),
    c(
      18407L,
      178L
    )
  ),
  "Canonical expression = 18,407 genes x 178 samples"
)


check(
  !anyDuplicated(
    rownames(canonical_expr)
  ),
  "Canonical expression gene symbols are unique"
)


check(
  identical(
    colnames(canonical_expr),
    colnames(expr)
  ),
  "Canonical expression sample order preserved exactly"
)


check(
  all(
    is.finite(
      canonical_expr
    )
  ),
  "Canonical expression contains only finite values"
)


# ============================================================
# 16. Exact-value preservation
# ============================================================

source_selected <- expr[
  selected_idx,
  ,
  drop = FALSE
]


max_abs_difference <- max(
  abs(
    as.numeric(canonical_expr) -
      as.numeric(source_selected)
  )
)


check(
  max_abs_difference == 0,
  "Canonical values are exactly unchanged from authoritative log2 RMA"
)


# ============================================================
# 17. Exact-tie validation in final matrix
#
# Reconfirm that choosing deterministic probe among the 105
# exact ties cannot affect canonical expression values.
# ============================================================

for (i in seq_len(nrow(exact_ties))) {

  ids <- strsplit(
    exact_ties$TiedCandidates[i],
    ";",
    fixed = TRUE
  )[[1]]


  ids <- trimws(ids)


  idx <- match(
    ids,
    expr_probe_ids
  )


  if (anyNA(idx)) {
    stop(
      paste0(
        "[FAIL] Missing exact-tie probe for ",
        exact_ties$FinalGeneID[i]
      ),
      call. = FALSE
    )
  }


  x <- expr[
    idx,
    ,
    drop = FALSE
  ]


  ref <- as.numeric(
    x[
      1L,
      ,
      drop = TRUE
    ]
  )


  equal_check <- vapply(
    seq_len(nrow(x)),
    function(j) {

      identical(
        as.numeric(
          x[
            j,
            ,
            drop = TRUE
          ]
        ),
        ref
      )
    },
    logical(1)
  )


  if (!all(equal_check)) {

    stop(
      paste0(
        "[FAIL] Exact-tie invariant violated for GeneID ",
        exact_ties$FinalGeneID[i]
      ),
      call. = FALSE
    )
  }
}


check(
  TRUE,
  "All 105 exact-tie invariants reconfirmed"
)


# ============================================================
# 18. Representative-rule summary
# ============================================================

rule_summary <- representative_mapping[
  ,
  .N,
  by = RepresentativeRule
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("FINAL REPRESENTATIVE RULE SUMMARY\n")
cat("============================================================\n")

print(
  rule_summary
)


# ============================================================
# 19. Canonical matrix summary
# ============================================================

matrix_summary <- data.table(
  Metric = c(
    "Authoritative_source_features",
    "Authoritative_samples",
    "Stage2D_current_genes",
    "Final_canonical_genes",
    "Single_cluster_genes",
    "Annotation_design_unique_genes",
    "Exact_duplicate_tie_genes",
    "Excluded_unresolved_nonexact_genes",
    "Canonical_samples",
    "NA_values",
    "Nonfinite_values",
    "Maximum_value_difference_vs_selected_source"
  ),

  Value = c(
    nrow(expr),
    ncol(expr),
    18422,
    nrow(canonical_expr),
    nrow(single_mapping),
    nrow(stage2f_mapping),
    nrow(exact_mapping),
    nrow(unresolved),
    ncol(canonical_expr),
    sum(is.na(canonical_expr)),
    sum(!is.finite(canonical_expr)),
    max_abs_difference
  )
)


cat("\n============================================================\n")
cat("CANONICAL EXPRESSION SUMMARY\n")
cat("============================================================\n")

print(
  matrix_summary
)


# ============================================================
# 20. Save artifacts
# ============================================================

canonical_rds <- file.path(
  out_dir,
  "GSE73661_canonical_log2_RMA_18407x178.rds"
)


saveRDS(
  canonical_expr,
  canonical_rds,
  compress = FALSE
)


fwrite(
  representative_mapping[
    ,
    .(
      FinalGeneID,
      FinalSymbol,
      RepresentativeProbeID,
      RepresentativeRule,
      FirstSourceRow
    )
  ],
  file.path(
    out_dir,
    "GSE73661_canonical_representative_mapping_18407.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  unresolved,
  file.path(
    out_dir,
    "GSE73661_unresolved_nonexact_multicluster_genes_15.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  rule_summary,
  file.path(
    out_dir,
    "GSE73661_stage3_representative_rule_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  matrix_summary,
  file.path(
    out_dir,
    "GSE73661_stage3_canonical_expression_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 21. RDS round-trip validation
# ============================================================

roundtrip <- readRDS(
  canonical_rds
)


check(
  identical(
    canonical_expr,
    roundtrip
  ),
  "Canonical RDS round-trip is exact"
)


# ============================================================
# 22. Checksums
# ============================================================

hash_files <- c(
  expr_file,
  decision_file,
  stage2f_file,
  stage2g_file,
  canonical_rds
)


hashes <- tools::md5sum(
  hash_files
)


fwrite(
  data.table(
    File =
      names(hashes),

    MD5 =
      as.character(
        unname(hashes)
      )
  ),
  file.path(
    out_dir,
    "GSE73661_stage3_artifact_md5.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 23. FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] Stage-2D gene identity mapping preserved\n")
cat("[PASS] 17,804 single-cluster genes retained directly\n")
cat("[PASS] 498 annotation/design-unique representatives retained\n")
cat("[PASS] 105 exact-numeric duplicate ties deterministically resolved\n")
cat("[PASS] 15 non-exact unresolved genes conservatively excluded\n")
cat("[PASS] Final canonical matrix contains 18,407 unique genes\n")
cat("[PASS] All 178 samples preserved in original order\n")
cat("[PASS] Author-provided log2-RMA values preserved exactly\n")
cat("[PASS] No averaging or median collapse performed\n")
cat("[PASS] No expression-based representative selection performed\n")
cat("[PASS] No normalization or batch correction performed\n")
cat("[PASS] Canonical RDS round-trip exact\n")

cat("\nFINAL STATUS:\n")
cat("PASS_GSE73661_STAGE3_CANONICAL_LOG2_RMA_EXPRESSION\n")
cat("============================================================\n")
