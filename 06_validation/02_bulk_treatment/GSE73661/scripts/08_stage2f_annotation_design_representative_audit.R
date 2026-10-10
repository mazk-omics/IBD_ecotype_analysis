suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 2F
# Annotation/design-based representative transcript-cluster audit
#
# Frozen inputs:
#   Stage 1  authoritative expression
#   Stage 2D FINAL probe -> current-gene mapping
#   Stage 2E expression-concordance audit
#   GPL6244 official platform table
#
# Goal:
#   Determine whether the 618 multi-cluster genes can be reduced
#   to ONE representative transcript cluster using annotation /
#   array-design evidence WITHOUT using expression as a primary
#   selection criterion.
#
# NO expression modification
# NO representative cluster frozen
# NO gene-level matrix generated
# ============================================================


# ============================================================
# Paths
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

soft_file <- file.path(
  root,
  "00_provenance/geo_downloads",
  "GSE73661_family.soft.gz"
)

decision_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_canonical_probe_gene_decision_table_FINAL.tsv"
)

cluster_stats_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2e_transcript_cluster_expression_stats.tsv"
)

gene_summary_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2e_multicluster_gene_summary.tsv"
)

candidate_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2e_candidate_representative_clusters.tsv"
)

out_dir <- file.path(
  root,
  "prepared/02_gene_harmonization"
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


as_num_safe <- function(x) {

  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 2F: annotation/design representative audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Input checks
# ============================================================

for (f in c(
  soft_file,
  decision_file,
  cluster_stats_file,
  gene_summary_file,
  candidate_file
)) {

  check(
    file.exists(f),
    paste0(
      "Input exists: ",
      basename(f)
    )
  )
}


decision <- fread(
  decision_file,
  showProgress = FALSE
)

cluster_stats <- fread(
  cluster_stats_file,
  showProgress = FALSE
)

gene_summary <- fread(
  gene_summary_file,
  showProgress = FALSE
)

candidate_expr <- fread(
  candidate_file,
  showProgress = FALSE
)


# ============================================================
# 2. Freeze identifier types
# ============================================================

for (dt in list(
  decision,
  cluster_stats,
  gene_summary,
  candidate_expr
)) {

  invisible(dt)
}


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


cluster_stats[
  ,
  SourceProbeID :=
    as.character(SourceProbeID)
]

cluster_stats[
  ,
  FinalGeneID :=
    as.character(FinalGeneID)
]

cluster_stats[
  ,
  FinalSymbol :=
    as.character(FinalSymbol)
]


gene_summary[
  ,
  FinalGeneID :=
    as.character(FinalGeneID)
]

gene_summary[
  ,
  FinalSymbol :=
    as.character(FinalSymbol)
]


candidate_expr[
  ,
  FinalGeneID :=
    as.character(FinalGeneID)
]

candidate_expr[
  ,
  FinalSymbol :=
    as.character(FinalSymbol)
]

candidate_expr[
  ,
  Candidate_MaxMean :=
    as.character(Candidate_MaxMean)
]

candidate_expr[
  ,
  Candidate_MaxIQR :=
    as.character(Candidate_MaxIQR)
]

candidate_expr[
  ,
  Candidate_SpearmanMedoid :=
    as.character(Candidate_SpearmanMedoid)
]


if (!is.logical(decision$IncludeFinal)) {

  decision[
    ,
    IncludeFinal :=
      as.logical(IncludeFinal)
  ]
}


# ============================================================
# 3. Frozen Stage-2D / Stage-2E structure
# ============================================================

retained <- decision[
  IncludeFinal == TRUE
]


check(
  nrow(retained) == 19229L,
  "19,229 retained transcript clusters recovered"
)

check(
  uniqueN(retained$FinalGeneID) == 18422L,
  "18,422 final genes recovered"
)


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


multi_genes <- gene_counts[
  N_SourceTranscriptClusters > 1L
]


check(
  nrow(multi_genes) == 618L,
  "618 multi-cluster genes recovered"
)


check(
  nrow(gene_summary) == 618L,
  "Stage-2E gene summary contains 618 genes"
)

check(
  nrow(candidate_expr) == 618L,
  "Stage-2E candidate table contains 618 genes"
)


# ============================================================
# 4. Extract GPL6244 official platform table
# ============================================================

awk_file <- tempfile(
  fileext = ".awk"
)

platform_tmp <- tempfile(
  fileext = ".tsv"
)


writeLines(
  c(
    'tolower($0) ~ /^!platform_table_begin/ { in_table=1; next }',
    'tolower($0) ~ /^!platform_table_end/   { in_table=0; next }',
    'in_table { print }'
  ),
  awk_file
)


status <- system(
  paste(
    "zcat",
    shQuote(soft_file),
    "| awk -f",
    shQuote(awk_file),
    ">",
    shQuote(platform_tmp)
  )
)


check(
  status == 0,
  "GPL6244 platform table extracted"
)


platform <- fread(
  platform_tmp,
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = FALSE,
  showProgress = FALSE
)


unlink(
  c(
    awk_file,
    platform_tmp
  )
)


check(
  nrow(platform) == 33297L,
  "GPL6244 contains 33,297 platform rows"
)

check(
  all(
    c(
      "ID",
      "total_probes",
      "category"
    ) %in%
      names(platform)
  ),
  "GPL6244 provides ID / total_probes / category"
)


platform[
  ,
  ID :=
    as.character(ID)
]


platform[
  ,
  TotalProbes :=
    as_num_safe(
      total_probes
    )
]


check(
  !anyNA(
    platform$ID
  ),
  "Platform IDs contain no NA"
)


# ============================================================
# 5. Attach GPL6244 design support
# ============================================================

platform_idx <- match(
  retained$SourceProbeID,
  platform$ID
)


check(
  !anyNA(platform_idx),
  "Every retained transcript cluster found in GPL6244"
)


retained[
  ,
  TotalProbes :=
    platform$TotalProbes[
      platform_idx
    ]
]


retained[
  ,
  PlatformCategoryCheck :=
    as.character(
      platform$category[
        platform_idx
      ]
    )
]


check(
  all(
    retained$PlatformCategoryCheck ==
      "main"
  ),
  "All retained transcript clusters are GPL6244 main-design features"
)


# ============================================================
# 6. Attach Stage-2E expression diagnostics
#
# These are NOT used for primary annotation/design ranking.
# ============================================================

cs_idx <- match(
  retained$SourceProbeID,
  cluster_stats$SourceProbeID
)


check(
  !anyNA(cs_idx),
  "Every retained cluster found in Stage-2E expression statistics"
)


retained[
  ,
  MeanExpression :=
    cluster_stats$MeanExpression[
      cs_idx
    ]
]


retained[
  ,
  IQRExpression :=
    cluster_stats$IQRExpression[
      cs_idx
    ]
]


# ============================================================
# 7. Restrict to the 618 multi-cluster genes
# ============================================================

multi <- retained[
  FinalGeneID %in%
    multi_genes$FinalGeneID
]


check(
  nrow(multi) == 1425L,
  "1,425 multi-cluster source features recovered"
)


# ============================================================
# 8. Annotation evidence rank
#
# Lower = stronger evidence.
#
# This is deliberately lexicographic rather than a summed score.
# ============================================================

multi[
  ,
  MappingEvidenceRank :=
    fcase(

      FinalDecision ==
        "KEEP_STRONG_CONCORDANT_CURRENT_GENE",
      1L,

      FinalDecision ==
        "KEEP_PLATFORM_DISAMBIGUATION_WITHIN_CHIP_SET",
      2L,

      FinalDecision ==
        "KEEP_CHIPDB_CURRENT_ONLY",
      3L,

      default =
        99L
    )
]


check(
  !any(
    multi$MappingEvidenceRank ==
      99L
  ),
  "Every multi-cluster feature belongs to an accepted mapping class"
)


# ============================================================
# 9. Platform-current support
#
# Prefer direct GPL6244 evidence resolving to exactly the same
# current GeneID.
# ============================================================

multi[
  ,
  PlatformCurrentSupportRank :=
    fifelse(
      N_PlatformCurrent == 1L &
        !is.na(
          PlatformCurrentGeneIDs
        ) &
        as.character(
          PlatformCurrentGeneIDs
        ) ==
          FinalGeneID,
      0L,
      1L
    )
]


# ============================================================
# 10. Annotation burden
#
# Fewer unresolved platform Entrez IDs is preferred.
# Fewer platform symbols is preferred only after current-GeneID
# support has already been considered.
# ============================================================

multi[
  ,
  UnresolvedPlatformEntrezRank :=
    fifelse(
      is.na(
        N_PlatformUnresolvedEntrez
      ),
      999L,
      as.integer(
        N_PlatformUnresolvedEntrez
      )
    )
]


multi[
  ,
  PlatformSymbolBurdenRank :=
    fifelse(
      is.na(
        N_PlatformSymbols
      ),
      999L,
      as.integer(
        N_PlatformSymbols
      )
    )
]


multi[
  ,
  TotalProbesRank :=
    fifelse(
      is.na(
        TotalProbes
      ),
      -Inf,
      TotalProbes
    )
]


# ============================================================
# 11. Lexicographic annotation/design ranking
#
# Priority:
#   1. mapping evidence
#   2. platform-current GeneID support
#   3. fewer unresolved platform Entrez IDs
#   4. lower platform symbol burden
#   5. more array probes
#
# Expression is NOT used.
# ============================================================

setorder(
  multi,
  FinalGeneID,
  MappingEvidenceRank,
  PlatformCurrentSupportRank,
  UnresolvedPlatformEntrezRank,
  PlatformSymbolBurdenRank,
  -TotalProbesRank,
  SourceProbeID
)


# ============================================================
# 12. Determine annotation/design ties
# ============================================================

multi[
  ,
  BestMappingRank :=
    min(
      MappingEvidenceRank
    ),
  by = FinalGeneID
]


multi[
  ,
  BestPlatformSupportRank :=
    min(
      PlatformCurrentSupportRank[
        MappingEvidenceRank ==
          BestMappingRank[1L]
      ]
    ),
  by = FinalGeneID
]


candidate_stage1 <- multi[
  MappingEvidenceRank ==
    BestMappingRank &
  PlatformCurrentSupportRank ==
    BestPlatformSupportRank
]


candidate_stage1[
  ,
  BestUnresolvedRank :=
    min(
      UnresolvedPlatformEntrezRank
    ),
  by = FinalGeneID
]


candidate_stage2 <- candidate_stage1[
  UnresolvedPlatformEntrezRank ==
    BestUnresolvedRank
]


candidate_stage2[
  ,
  BestSymbolBurden :=
    min(
      PlatformSymbolBurdenRank
    ),
  by = FinalGeneID
]


candidate_stage3 <- candidate_stage2[
  PlatformSymbolBurdenRank ==
    BestSymbolBurden
]


candidate_stage3[
  ,
  BestTotalProbes :=
    max(
      TotalProbesRank
    ),
  by = FinalGeneID
]


candidate_final <- candidate_stage3[
  TotalProbesRank ==
    BestTotalProbes
]


# ============================================================
# 13. Gene-level annotation/design candidate summary
# ============================================================

annotation_candidates <- candidate_final[
  ,
  .(
    N_AnnotationDesignCandidates =
      .N,

    Candidate_AnnotationDesign =
      if (
        .N == 1L
      ) {
        SourceProbeID[1L]
      } else {
        NA_character_
      },

    TiedCandidates =
      paste(
        SourceProbeID,
        collapse = ";"
      )
  ),
  by = .(
    FinalGeneID,
    FinalSymbol
  )
]


check(
  nrow(annotation_candidates) == 618L,
  "Annotation/design candidate summary contains 618 genes"
)


annotation_candidates[
  gene_counts,
  on = c(
    "FinalGeneID",
    "FinalSymbol"
  ),
  N_SourceTranscriptClusters :=
    i.N_SourceTranscriptClusters
]


annotation_candidates[
  gene_summary,
  on = c(
    "FinalGeneID",
    "FinalSymbol"
  ),
  `:=`(
    MedianSpearman =
      i.MedianSpearman,

    MinSpearman =
      i.MinSpearman
  )
]


annotation_candidates[
  candidate_expr,
  on = c(
    "FinalGeneID",
    "FinalSymbol"
  ),
  `:=`(
    Candidate_MaxMean =
      i.Candidate_MaxMean,

    Candidate_MaxIQR =
      i.Candidate_MaxIQR,

    Candidate_SpearmanMedoid =
      i.Candidate_SpearmanMedoid
  )
]


# ============================================================
# 14. Diagnostic comparison to expression-based candidates
#
# NOT used to choose the annotation/design candidate.
# ============================================================

annotation_candidates[
  ,
  AnnotationUnique :=
    N_AnnotationDesignCandidates ==
      1L
]


annotation_candidates[
  ,
  Annotation_equals_MaxMean :=
    AnnotationUnique &
    Candidate_AnnotationDesign ==
      Candidate_MaxMean
]


annotation_candidates[
  ,
  Annotation_equals_MaxIQR :=
    AnnotationUnique &
    Candidate_AnnotationDesign ==
      Candidate_MaxIQR
]


annotation_candidates[
  ,
  Annotation_equals_Medoid :=
    AnnotationUnique &
    Candidate_AnnotationDesign ==
      Candidate_SpearmanMedoid
]


# ============================================================
# 15. Overall summary
# ============================================================

n_unique <- annotation_candidates[
  AnnotationUnique == TRUE,
  .N
]


n_tied <- annotation_candidates[
  AnnotationUnique == FALSE,
  .N
]


summary <- data.table(
  Metric = c(
    "Multi_cluster_genes",
    "Annotation_design_unique_candidate",
    "Annotation_design_remaining_ties",
    "Unique_candidate_equals_MaxMean",
    "Unique_candidate_equals_MaxIQR",
    "Unique_candidate_equals_SpearmanMedoid"
  ),

  N = c(
    618L,

    n_unique,

    n_tied,

    annotation_candidates[
      Annotation_equals_MaxMean == TRUE,
      .N
    ],

    annotation_candidates[
      Annotation_equals_MaxIQR == TRUE,
      .N
    ],

    annotation_candidates[
      Annotation_equals_Medoid == TRUE,
      .N
    ]
  )
)


cat("\n============================================================\n")
cat("ANNOTATION/DESIGN REPRESENTATIVE AUDIT SUMMARY\n")
cat("============================================================\n")

print(summary)


# ============================================================
# 16. Stratify by transcript-cluster multiplicity
# ============================================================

multiplicity_summary <- annotation_candidates[
  ,
  .(
    Genes =
      .N,

    UniqueAnnotationCandidate =
      sum(
        AnnotationUnique
      ),

    RemainingTies =
      sum(
        !AnnotationUnique
      ),

    EqualsMaxMean =
      sum(
        Annotation_equals_MaxMean
      ),

    EqualsMaxIQR =
      sum(
        Annotation_equals_MaxIQR
      ),

    EqualsMedoid =
      sum(
        Annotation_equals_Medoid
      )
  ),
  by = N_SourceTranscriptClusters
][
  order(
    N_SourceTranscriptClusters
  )
]


cat("\n============================================================\n")
cat("AUDIT BY SOURCE-CLUSTER MULTIPLICITY\n")
cat("============================================================\n")

print(multiplicity_summary)


# ============================================================
# 17. Concordance of unique annotation candidates
# ============================================================

concordance_summary <- annotation_candidates[
  ,
  .(
    Genes =
      .N,

    UniqueCandidate =
      sum(
        AnnotationUnique
      ),

    RemainingTies =
      sum(
        !AnnotationUnique
      )
  ),
  by = .(
    ConcordanceClass =
      fcase(

        MedianSpearman >= 0.90,
        "MedianSpearman>=0.90",

        MedianSpearman >= 0.75,
        "0.75-0.90",

        MedianSpearman >= 0.50,
        "0.50-0.75",

        MedianSpearman >= 0,
        "0-0.50",

        MedianSpearman < 0,
        "<0",

        default =
          "NA"
      )
  )
]


cat("\n============================================================\n")
cat("ANNOTATION CANDIDATES BY EXPRESSION CONCORDANCE\n")
cat("============================================================\n")

print(concordance_summary)


# ============================================================
# 18. Remaining ties
# ============================================================

remaining_ties <- annotation_candidates[
  AnnotationUnique == FALSE
]


cat("\n============================================================\n")
cat("REMAINING ANNOTATION/DESIGN TIES\n")
cat("============================================================\n")

cat(
  "Genes still tied after non-expression ranking: ",
  nrow(remaining_ties),
  "\n",
  sep = ""
)


if (nrow(remaining_ties) > 0L) {

  print(
    head(
      remaining_ties[
        ,
        .(
          FinalGeneID,
          FinalSymbol,
          N_SourceTranscriptClusters,
          N_AnnotationDesignCandidates,
          TiedCandidates,
          MedianSpearman,
          MinSpearman,
          Candidate_MaxMean,
          Candidate_MaxIQR,
          Candidate_SpearmanMedoid
        )
      ],
      30L
    )
  )
}


# ============================================================
# 19. Save outputs
# ============================================================

fwrite(
  multi,
  file.path(
    out_dir,
    "GSE73661_stage2f_multicluster_annotation_design_features.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  annotation_candidates,
  file.path(
    out_dir,
    "GSE73661_stage2f_annotation_design_candidates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  summary,
  file.path(
    out_dir,
    "GSE73661_stage2f_annotation_design_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  multiplicity_summary,
  file.path(
    out_dir,
    "GSE73661_stage2f_multiplicity_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  concordance_summary,
  file.path(
    out_dir,
    "GSE73661_stage2f_concordance_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  remaining_ties,
  file.path(
    out_dir,
    "GSE73661_stage2f_remaining_annotation_design_ties.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


# ============================================================
# 20. Final structural checks
# ============================================================

check(
  nrow(annotation_candidates) == 618L,
  "All 618 multi-cluster genes evaluated"
)


check(
  n_unique + n_tied == 618L,
  "Unique + tied annotation candidates account for all 618 genes"
)


check(
  nrow(multi) == 1425L,
  "All 1,425 multi-cluster source features audited"
)


# ============================================================
# 21. FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] 618 multi-cluster genes annotation/design-audited\n")
cat("[PASS] Mapping evidence incorporated\n")
cat("[PASS] Current GPL6244 GeneID support incorporated\n")
cat("[PASS] Unresolved annotation burden incorporated\n")
cat("[PASS] Platform symbol burden incorporated\n")
cat("[PASS] GPL6244 total-probe design support incorporated\n")
cat("[PASS] Expression was not used for primary selection\n")
cat("[PASS] Expression-based candidates retained only for diagnostics\n")
cat("[PASS] No transcript cluster frozen as representative\n")
cat("[PASS] No expression values modified\n")
cat("[PASS] No gene-level matrix generated\n")

cat("\nFINAL STATUS:\n")
cat("PASS_GSE73661_STAGE2F_ANNOTATION_DESIGN_REPRESENTATIVE_AUDIT\n")
cat("============================================================\n")
