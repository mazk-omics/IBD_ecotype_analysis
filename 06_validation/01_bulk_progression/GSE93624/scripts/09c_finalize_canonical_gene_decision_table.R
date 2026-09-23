suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 2F FINAL
#
# Final canonical gene reconciliation.
#
# Key corrections:
#   1. retired GeneIDs are allowed to redirect to current GeneIDs
#   2. grouping is by FINAL CURRENT GeneID only
#   3. stable/non-redirect source rows are preferred over
#      historical redirected rows
#   4. if multiple stable rows remain, current-primary source
#      symbol is preferred when uniquely present
#   5. expression rows are NEVER summed/averaged/collapsed
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

harm_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

decision_file <- file.path(
  harm_dir,
  "GSE93624_canonical_gene_decision_table.tsv"
)

expr_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE93624_authoritative_expression_245.rds"
)

gene_info_file <- paste0(
  "/home/mazekai/IBD_EcoTyper/06_validation/99_common/",
  "gene_annotation/NCBI_Gene_2026-09-21/",
  "Homo_sapiens.gene_info.gz"
)

ipw_redirect_file <- file.path(
  harm_dir,
  "final_retired_geneid_audit",
  "GeneID_3653.tsv"
)

cat("============================================================\n")
cat("GSE93624 FINAL current-GeneID canonical reconciliation\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load Stage 09 decision table
# ============================================================

d <- fread(
  decision_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  nrow(d) == 13769,
  !anyDuplicated(d$SourceSymbol)
)

cat("[PASS] Stage 09 decision table loaded: 13,769 features\n")


# ============================================================
# 2. Candidate vs locked exclusions
# ============================================================

candidate_classes <- c(
  "KEEP_UNIQUE_CURRENT_IDENTITY",
  "KEEP_CURRENT_PRIMARY_COLLISION_GROUP",
  "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",
  "REVIEW_MANY_TO_ONE_NO_PRIMARY",
  "REVIEW_UNCLASSIFIED",
  "REVIEW_FINAL_SYMBOL_COLLISION"
)

cand <- copy(
  d[
    Decision %chin% candidate_classes
  ]
)

locked_exclusion <- copy(
  d[
    !Decision %chin% candidate_classes
  ]
)

stopifnot(
  nrow(cand) == 13334,
  nrow(locked_exclusion) == 435,
  nrow(cand) + nrow(locked_exclusion) == 13769
)

cat("[PASS] Candidate identity rows: 13,334\n")
cat("[PASS] Locked non-candidate exclusions: 435\n")


# ============================================================
# 3. Load authoritative expression
#    Used only for diagnostic collision correlations.
# ============================================================

obj <- readRDS(expr_file)
expr <- obj$expression

stopifnot(
  nrow(expr) == 13769,
  ncol(expr) == 245,
  all(cand$SourceSymbol %chin% rownames(expr))
)


# ============================================================
# 4. Current NCBI gene_info
# ============================================================

gi <- fread(
  gene_info_file,
  header = TRUE,
  sep = "\t",
  select = c(
    "GeneID",
    "Symbol",
    "description",
    "type_of_gene"
  ),
  colClasses = "character",
  na.strings = c("-", "")
)

setnames(
  gi,
  c(
    "GeneID",
    "Symbol",
    "description",
    "type_of_gene"
  ),
  c(
    "NCBI_GeneID",
    "NCBI_CurrentSymbol",
    "NCBI_Description",
    "NCBI_GeneType"
  )
)

stopifnot(
  !anyDuplicated(gi$NCBI_GeneID)
)

cat(
  "[PASS] Current NCBI records:",
  nrow(gi),
  "\n"
)


# ============================================================
# 5. Identify candidate GeneIDs absent from current gene_info
# ============================================================

cand[
  ,
  OriginalCandidateGeneID :=
    CandidateGeneID
]

cand[
  ,
  FinalResolvedGeneID :=
    CandidateGeneID
]

initial_idx <- match(
  cand$FinalResolvedGeneID,
  gi$NCBI_GeneID
)

missing_ids <- unique(
  cand$FinalResolvedGeneID[
    is.na(initial_idx)
  ]
)

cat("\n============================================================\n")
cat("RETIRED / ABSENT CANDIDATE GeneIDs\n")
cat("============================================================\n")

cat(
  "Unique IDs absent from current gene_info:",
  length(missing_ids),
  "\n"
)

print(missing_ids)

# At the current audit state this should be exactly GeneID 3653.
stopifnot(
  identical(
    sort(missing_ids),
    "3653"
  )
)


# ============================================================
# 6. Parse the independently queried IPW redirect
# ============================================================

if (!file.exists(ipw_redirect_file)) {
  stop(
    "Missing independent NCBI redirect evidence: ",
    ipw_redirect_file
  )
}

ipw <- fread(
  ipw_redirect_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = TRUE,
  colClasses = "character",
  na.strings = c("", "NA")
)

if (nrow(ipw) != 1) {
  stop(
    "Expected exactly one returned current Gene record for GeneID 3653."
  )
}

if (ncol(ipw) != 5) {
  stop(
    "Expected five dataformat columns in GeneID_3653.tsv."
  )
}

setnames(
  ipw,
  c(
    "ResolvedGeneID",
    "ResolvedSymbol",
    "GeneType",
    "Description",
    "ReplacedGeneID"
  )
)

cat("\n============================================================\n")
cat("IPW RETIRED GeneID RECONCILIATION\n")
cat("============================================================\n")

print(ipw)

stopifnot(
  ipw$ResolvedGeneID[1] != "3653",
  ipw$ResolvedGeneID[1] %chin% gi$NCBI_GeneID
)


# ============================================================
# 7. Apply retired-ID redirect
# ============================================================

cand[
  OriginalCandidateGeneID == "3653",
  FinalResolvedGeneID :=
    ipw$ResolvedGeneID[1]
]

cand[
  ,
  DynamicGeneIDRedirect :=
    OriginalCandidateGeneID !=
      FinalResolvedGeneID
]


# ============================================================
# 8. Resolve every candidate to current NCBI symbol
# ============================================================

idx <- match(
  cand$FinalResolvedGeneID,
  gi$NCBI_GeneID
)

stopifnot(
  !anyNA(idx)
)

cand[
  ,
  FinalResolvedSymbol :=
    gi$NCBI_CurrentSymbol[idx]
]

cand[
  ,
  FinalGeneType :=
    gi$NCBI_GeneType[idx]
]

cand[
  ,
  FinalDescription :=
    gi$NCBI_Description[idx]
]


cat("\n============================================================\n")
cat("FINAL CURRENT NCBI COVERAGE\n")
cat("============================================================\n")

cat(
  "Candidate rows with current GeneID:",
  sum(!is.na(idx)),
  "/",
  nrow(cand),
  "\n"
)

stopifnot(
  all(!is.na(cand$FinalResolvedGeneID)),
  all(!is.na(cand$FinalResolvedSymbol))
)


# ============================================================
# 9. Define historical redirect evidence
#
# LOC_GENEID_REDIRECT:
#   old LOC GeneID was externally redirected before Stage 09.
#
# DynamicGeneIDRedirect:
#   IPW 3653 -> current GeneID.
#
# These are weaker row-level representatives than an already
# deposited row whose stable GeneID is the current GeneID.
# ============================================================

cand[
  ,
  HistoricalGeneIDRedirect :=
    EvidenceClass ==
      "LOC_GENEID_REDIRECT" |
    DynamicGeneIDRedirect
]


# ============================================================
# 10. Group by FINAL GeneID ONLY
# ============================================================

groups <- cand[
  ,
  .(
    SourceFeatureN = .N,

    SourceSymbols =
      paste(
        SourceSymbol,
        collapse = "|"
      ),

    NonRedirectSourceN =
      sum(
        !HistoricalGeneIDRedirect
      ),

    CurrentPrimarySourceN =
      sum(
        SourceSymbol ==
          FinalResolvedSymbol
      ),

    NonRedirectPrimarySourceN =
      sum(
        !HistoricalGeneIDRedirect &
        SourceSymbol ==
          FinalResolvedSymbol
      )
  ),
  by = FinalResolvedGeneID
]


gidx <- match(
  cand$FinalResolvedGeneID,
  groups$FinalResolvedGeneID
)

cand[
  ,
  GeneIDSourceN :=
    groups$SourceFeatureN[gidx]
]

cand[
  ,
  NonRedirectSourceN :=
    groups$NonRedirectSourceN[gidx]
]

cand[
  ,
  CurrentPrimarySourceN :=
    groups$CurrentPrimarySourceN[gidx]
]

cand[
  ,
  NonRedirectPrimarySourceN :=
    groups$NonRedirectPrimarySourceN[gidx]
]

cand[
  ,
  CollisionSourceSymbols :=
    groups$SourceSymbols[gidx]
]


# ============================================================
# 11. Final collision decision hierarchy
# ============================================================

cand[
  ,
  FinalDecision := fcase(

    # ---------------------------------------
    # A. Unique current GeneID representation
    # ---------------------------------------
    GeneIDSourceN == 1,
    "KEEP_UNIQUE_CURRENT_IDENTITY",


    # ---------------------------------------
    # B. Exactly one non-redirect source row
    #
    # Stable/current GeneID representation wins over historical
    # retired-ID redirect.
    # ---------------------------------------
    GeneIDSourceN > 1 &
      NonRedirectSourceN == 1 &
      !HistoricalGeneIDRedirect,
    "KEEP_STABLE_GENEID_SOURCE",

    GeneIDSourceN > 1 &
      NonRedirectSourceN == 1 &
      HistoricalGeneIDRedirect,
    "EXCLUDE_HISTORICAL_REDIRECT_SOURCE",


    # ---------------------------------------
    # C. Multiple non-redirect rows remain,
    # but exactly one uses current primary symbol
    # ---------------------------------------
    GeneIDSourceN > 1 &
      NonRedirectSourceN > 1 &
      CurrentPrimarySourceN == 1 &
      SourceSymbol ==
        FinalResolvedSymbol,
    "KEEP_CURRENT_PRIMARY_COLLISION_GROUP",

    GeneIDSourceN > 1 &
      NonRedirectSourceN > 1 &
      CurrentPrimarySourceN == 1 &
      SourceSymbol !=
        FinalResolvedSymbol,
    "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",


    # ---------------------------------------
    # D. Multiple rows, no defensible winner
    # ---------------------------------------
    GeneIDSourceN > 1,
    "REVIEW_UNRESOLVED_GENEID_COLLISION",

    default =
      "REVIEW_UNCLASSIFIED"
  )
]


cand[
  ,
  FinalDecisionReason := fcase(

    FinalDecision ==
      "KEEP_UNIQUE_CURRENT_IDENTITY",
    "Only one deposited source row resolves to this current GeneID",

    FinalDecision ==
      "KEEP_STABLE_GENEID_SOURCE",
    paste0(
      "Multiple source rows resolve to the same current GeneID; ",
      "this row already represents the stable/current GeneID, ",
      "whereas competing row(s) entered through historical GeneID redirects"
    ),

    FinalDecision ==
      "EXCLUDE_HISTORICAL_REDIRECT_SOURCE",
    paste0(
      "Historical/retired GeneID redirects to a current GeneID ",
      "already represented by a stable deposited source row; ",
      "no expression aggregation performed"
    ),

    FinalDecision ==
      "KEEP_CURRENT_PRIMARY_COLLISION_GROUP",
    paste0(
      "Multiple non-redirect source rows resolve to the same GeneID; ",
      "this row uniquely uses the current NCBI primary symbol"
    ),

    FinalDecision ==
      "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",
    paste0(
      "Same current GeneID is represented by another deposited row ",
      "using the unique current-primary symbol; no collapse performed"
    ),

    FinalDecision ==
      "REVIEW_UNRESOLVED_GENEID_COLLISION",
    paste0(
      "Multiple deposited expression rows resolve to the same current GeneID ",
      "without a unique stable/current-primary representative"
    ),

    default =
      "Unclassified"
  )
]


# ============================================================
# 12. Collision report
# ============================================================

collision <- cand[
  GeneIDSourceN > 1
][
  order(
    FinalResolvedGeneID,
    SourceSymbol
  )
]


cat("\n============================================================\n")
cat("FINAL GeneID-ONLY COLLISION AUDIT\n")
cat("============================================================\n")

cat(
  "Collided current GeneIDs:",
  uniqueN(
    collision$FinalResolvedGeneID
  ),
  "\n"
)

cat(
  "Source rows involved:",
  nrow(collision),
  "\n"
)

print(
  collision[
    ,
    .(
      SourceSymbol,
      EvidenceClass,
      OriginalCandidateGeneID,
      FinalResolvedGeneID,
      FinalResolvedSymbol,
      HistoricalGeneIDRedirect,
      GeneIDSourceN,
      NonRedirectSourceN,
      CurrentPrimarySourceN,
      CollisionSourceSymbols,
      FinalDecision
    )
  ]
)


# ============================================================
# 13. Diagnostic expression collision audit
#
# Never used to automatically collapse rows.
# ============================================================

collision_ids <- unique(
  collision$FinalResolvedGeneID
)

pair_list <- list()
k <- 0L

for (gid in collision_ids) {

  ss <- cand[
    FinalResolvedGeneID == gid,
    SourceSymbol
  ]

  if (length(ss) < 2)
    next

  pp <- combn(
    ss,
    2,
    simplify = FALSE
  )

  for (p in pp) {

    k <- k + 1L

    x <- as.numeric(
      expr[p[1], ]
    )

    y <- as.numeric(
      expr[p[2], ]
    )

    pair_list[[k]] <- data.table(
      FinalGeneID =
        gid,

      FinalSymbol =
        cand[
          FinalResolvedGeneID == gid,
          FinalResolvedSymbol
        ][1],

      Source1 =
        p[1],

      Source2 =
        p[2],

      ExactDuplicate =
        identical(x, y),

      Pearson =
        suppressWarnings(
          cor(
            x,
            y,
            method = "pearson"
          )
        ),

      Spearman =
        suppressWarnings(
          cor(
            x,
            y,
            method = "spearman"
          )
        ),

      MeanAbsDiff =
        mean(
          abs(x - y)
        ),

      MaxAbsDiff =
        max(
          abs(x - y)
        )
    )
  }
}

collision_expr <- rbindlist(
  pair_list,
  use.names = TRUE,
  fill = TRUE
)


cat("\n============================================================\n")
cat("COLLISION EXPRESSION DIAGNOSTICS\n")
cat("============================================================\n")

print(
  collision_expr
)


# ============================================================
# 14. Review rows
# ============================================================

review <- cand[
  grepl(
    "^REVIEW_",
    FinalDecision
  )
]


cat("\n============================================================\n")
cat("REVIEW-REQUIRED CANDIDATE ROWS\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(review),
  "\n"
)

if (nrow(review) > 0) {

  print(
    review[
      ,
      .(
        SourceSymbol,
        EvidenceClass,
        FinalResolvedGeneID,
        FinalResolvedSymbol,
        HistoricalGeneIDRedirect,
        CollisionSourceSymbols,
        FinalDecision
      )
    ]
  )
}


# ============================================================
# 15. Final inclusion
# ============================================================

keep_classes <- c(
  "KEEP_UNIQUE_CURRENT_IDENTITY",
  "KEEP_STABLE_GENEID_SOURCE",
  "KEEP_CURRENT_PRIMARY_COLLISION_GROUP"
)

cand[
  ,
  IncludeCanonical_FINAL :=
    FinalDecision %chin%
      keep_classes
]

cand[
  ,
  FinalGeneID_FINAL :=
    fifelse(
      IncludeCanonical_FINAL,
      FinalResolvedGeneID,
      NA_character_
    )
]

cand[
  ,
  FinalSymbol_FINAL :=
    fifelse(
      IncludeCanonical_FINAL,
      FinalResolvedSymbol,
      NA_character_
    )
]


# ============================================================
# 16. Rebuild complete 13,769-row table
# ============================================================

candidate_final <- cand[
  ,
  .(
    FeatureOrder =
      as.integer(FeatureOrder),

    SourceSymbol,

    EvidenceClass,
    EvidenceDetail,

    OriginalCandidateGeneID,
    OriginalCandidateSymbol =
      CandidateSymbol,

    FinalResolvedGeneID,
    FinalResolvedSymbol,

    FinalGeneType,
    FinalDescription,

    DynamicGeneIDRedirect,
    HistoricalGeneIDRedirect,

    GeneIDSourceN,
    NonRedirectSourceN,
    CurrentPrimarySourceN,
    CollisionSourceSymbols,

    Decision =
      FinalDecision,

    DecisionReason =
      FinalDecisionReason,

    IncludeCanonical =
      IncludeCanonical_FINAL,

    FinalGeneID =
      FinalGeneID_FINAL,

    FinalSymbol =
      FinalSymbol_FINAL
  )
]


locked_final <- locked_exclusion[
  ,
  .(
    FeatureOrder =
      as.integer(FeatureOrder),

    SourceSymbol,

    EvidenceClass,
    EvidenceDetail,

    OriginalCandidateGeneID =
      CandidateGeneID,

    OriginalCandidateSymbol =
      CandidateSymbol,

    FinalResolvedGeneID =
      NA_character_,

    FinalResolvedSymbol =
      NA_character_,

    FinalGeneType =
      NA_character_,

    FinalDescription =
      NA_character_,

    DynamicGeneIDRedirect =
      NA,

    HistoricalGeneIDRedirect =
      NA,

    GeneIDSourceN =
      NA_integer_,

    NonRedirectSourceN =
      NA_integer_,

    CurrentPrimarySourceN =
      NA_integer_,

    CollisionSourceSymbols =
      NA_character_,

    Decision,
    DecisionReason,

    IncludeCanonical =
      FALSE,

    FinalGeneID =
      NA_character_,

    FinalSymbol =
      NA_character_
  )
]


final <- rbindlist(
  list(
    candidate_final,
    locked_final
  ),
  use.names = TRUE,
  fill = TRUE
)

setorder(
  final,
  FeatureOrder
)

stopifnot(
  nrow(final) == 13769,
  identical(
    final$FeatureOrder,
    1:13769
  ),
  !anyDuplicated(
    final$SourceSymbol
  )
)


# ============================================================
# 17. Final uniqueness checks
# ============================================================

included <- final[
  IncludeCanonical == TRUE
]

excluded <- final[
  IncludeCanonical == FALSE
]

final_review <- final[
  grepl(
    "^REVIEW_",
    Decision
  )
]

dup_gid <- included[
  duplicated(FinalGeneID) |
  duplicated(
    FinalGeneID,
    fromLast = TRUE
  )
]

dup_symbol <- included[
  duplicated(FinalSymbol) |
  duplicated(
    FinalSymbol,
    fromLast = TRUE
  )
]


cat("\n============================================================\n")
cat("FINAL DECISION SUMMARY\n")
cat("============================================================\n")

decision_summary <- final[
  ,
  .N,
  by = .(
    IncludeCanonical,
    Decision
  )
][
  order(
    -IncludeCanonical,
    -N
  )
]

print(
  decision_summary
)


cat(
  "\nIncluded:",
  nrow(included),
  "\n"
)

cat(
  "Excluded:",
  nrow(excluded),
  "\n"
)

cat(
  "Review:",
  nrow(final_review),
  "\n"
)

cat(
  "Unique FinalGeneID:",
  uniqueN(
    included$FinalGeneID
  ),
  "\n"
)

cat(
  "Unique FinalSymbol:",
  uniqueN(
    included$FinalSymbol
  ),
  "\n"
)


# ============================================================
# 18. Check FAM83F/SACK1F and IPW/SNHG14 explicitly
# ============================================================

cat("\n============================================================\n")
cat("KEY FINAL RECONCILIATIONS\n")
cat("============================================================\n")

print(
  final[
    SourceSymbol %chin% c(
      "FAM83F",
      "LOC100130899",
      "IPW",
      "SNHG14"
    ),
    .(
      SourceSymbol,
      EvidenceClass,
      OriginalCandidateGeneID,
      FinalResolvedGeneID,
      FinalResolvedSymbol,
      HistoricalGeneIDRedirect,
      Decision,
      IncludeCanonical,
      FinalGeneID,
      FinalSymbol
    )
  ]
)


# ============================================================
# 19. Save FINAL authoritative decision resources
# ============================================================

fwrite(
  final,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_decision_table_FINAL.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  included,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_included_features_FINAL.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  excluded,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_excluded_features_FINAL.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  collision,
  file.path(
    harm_dir,
    "GSE93624_FINAL_current_GeneID_collision_audit.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  collision_expr,
  file.path(
    harm_dir,
    "GSE93624_FINAL_collision_expression_diagnostics.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  decision_summary,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_decision_summary_FINAL.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 20. Final hard guards
# ============================================================

cat("\n============================================================\n")
cat("FINAL INTEGRITY CHECKS\n")
cat("============================================================\n")

cat("[PASS] 13,769 raw features accounted exactly once\n")


if (
  nrow(final_review) == 0 &&
  nrow(dup_gid) == 0 &&
  nrow(dup_symbol) == 0 &&
  nrow(included) ==
    uniqueN(included$FinalGeneID) &&
  nrow(included) ==
    uniqueN(included$FinalSymbol)
) {

  cat("[PASS] No review rows remain\n")
  cat("[PASS] FinalGeneID is unique\n")
  cat("[PASS] FinalSymbol is unique\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE93624_CANONICAL_GENE_HARMONIZATION_FINAL\n")

} else {

  if (nrow(final_review) > 0) {
    cat(
      "[REVIEW] Remaining review rows:",
      nrow(final_review),
      "\n"
    )
  }

  if (nrow(dup_gid) > 0) {
    cat(
      "[FAIL] Duplicate FinalGeneID rows:",
      nrow(dup_gid),
      "\n"
    )
  }

  if (nrow(dup_symbol) > 0) {
    cat(
      "[FAIL] Duplicate FinalSymbol rows:",
      nrow(dup_symbol),
      "\n"
    )
  }

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_REQUIRED_GSE93624_CANONICAL_GENE_HARMONIZATION\n")
}

cat("============================================================\n")
