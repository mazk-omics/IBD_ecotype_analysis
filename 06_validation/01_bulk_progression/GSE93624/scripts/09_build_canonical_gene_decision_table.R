suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 2F
#
# Final global gene-identity reconciliation +
# canonical gene decision table
#
# Inputs integrate all evidence from:
#   - authoritative raw symbol matrix
#   - org.Hs current-symbol / alias mapping
#   - punctuation rescue
#   - local NCBI gene_info
#   - expression collision audits
#   - external NCBI Datasets forensic
#
# PRINCIPLES:
#   - one canonical row = one current GeneID
#   - GeneID is the stable identity anchor
#   - current NCBI primary symbol is the canonical symbol
#   - no expression aggregation
#   - no arbitrary choice among ambiguous identities
#   - historical source row is excluded if current primary
#     source row already exists
#
# This script builds the DECISION TABLE only.
# It does not modify the expression matrix.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

harm_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

query_dir <- file.path(
  harm_dir,
  "ncbi_external_queries"
)

expr_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE93624_authoritative_expression_245.rds"
)

audit_file <- file.path(
  harm_dir,
  "GSE93624_gene_symbol_harmonization_audit.tsv"
)

unresolved_file <- file.path(
  harm_dir,
  "GSE93624_unresolved_symbol_structural_audit.tsv"
)

alias_validation_file <- file.path(
  harm_dir,
  "GSE93624_safe_alias_NCBI_validation.tsv"
)

punct_validation_file <- file.path(
  harm_dir,
  "GSE93624_safe_punctuation_NCBI_validation.tsv"
)

loc_file <- file.path(
  harm_dir,
  "GSE93624_LOC_style_NCBI_current_audit.tsv"
)

external_unique_file <- file.path(
  query_dir,
  "GSE93624_NCBI_external_unique_identity_candidates_CORRECTED.tsv"
)

external_unresolved_file <- file.path(
  query_dir,
  "GSE93624_NCBI_external_unresolved_queries_CORRECTED.tsv"
)


cat("============================================================\n")
cat("GSE93624 canonical gene decision table\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load authoritative expression
# ============================================================

obj <- readRDS(
  expr_file
)

expr <- obj$expression

raw_symbols <- rownames(expr)

stopifnot(
  nrow(expr) == 13769,
  ncol(expr) == 245,
  uniqueN(raw_symbols) == 13769
)

cat(
  "[PASS] Authoritative matrix:",
  nrow(expr),
  "x",
  ncol(expr),
  "\n"
)


# ============================================================
# 2. Load forensic resources
# ============================================================

a <- fread(
  audit_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

u <- fread(
  unresolved_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

alias_val <- fread(
  alias_validation_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

punct_val <- fread(
  punct_validation_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

loc <- fread(
  loc_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

ext_unique <- fread(
  external_unique_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

ext_unresolved <- fread(
  external_unresolved_file,
  colClasses = "character",
  na.strings = c("", "NA")
)


stopifnot(
  nrow(a) == 13769,
  nrow(u) == 776,
  nrow(alias_val) == 748,
  nrow(punct_val) == 197,
  nrow(loc) == 133,
  nrow(ext_unique) == 71,
  nrow(ext_unresolved) == 38
)

stopifnot(
  identical(
    a$GeneSymbol_raw,
    raw_symbols
  )
)

cat("[PASS] All Stage 2 forensic resources loaded\n\n")


# ============================================================
# 3. Candidate identity: 12,208 direct current symbols
# ============================================================

cand_direct <- a[
  PreliminaryIdentityStatus ==
    "CURRENT_IDENTITY_OK",
  .(
    SourceSymbol =
      GeneSymbol_raw,

    CandidateGeneID =
      Direct_EntrezIDs,

    CandidateSymbol =
      GeneSymbol_raw,

    EvidenceClass =
      "DIRECT_CURRENT_ORGHS",

    EvidenceDetail =
      "Exact current org.Hs SYMBOL with unique EntrezID"
  )
]

stopifnot(
  nrow(cand_direct) == 12208
)


# ============================================================
# 4. Candidate identity: 748 safe aliases
#
# IMPORTANT:
# Use NCBI current symbol from GeneID validation.
# This automatically handles the 10 symbol-version updates.
# ============================================================

cand_alias <- alias_val[
  ,
  .(
    SourceSymbol,

    CandidateGeneID,

    CandidateSymbol =
      NCBICurrentSymbol,

    EvidenceClass =
      "SAFE_ALIAS_GENEID_VALIDATED",

    EvidenceDetail =
      fifelse(
        ValidationStatus ==
          "PASS_CURRENT_NCBI_ID_SYMBOL",
        "org.Hs alias and current NCBI GeneID/symbol concordant",
        "GeneID validated; current NCBI primary symbol supersedes older candidate symbol"
      )
  )
]

stopifnot(
  nrow(cand_alias) == 748,
  all(
    !is.na(
      cand_alias$CandidateGeneID
    )
  ),
  all(
    !is.na(
      cand_alias$CandidateSymbol
    )
  )
)


# ============================================================
# 5. Candidate identity: 197 punctuation rescues
# ============================================================

cand_punct <- punct_val[
  ,
  .(
    SourceSymbol,

    CandidateGeneID,

    CandidateSymbol =
      NCBICurrentSymbol,

    EvidenceClass =
      "PUNCTUATION_RESCUE_GENEID_VALIDATED",

    EvidenceDetail =
      fifelse(
        ValidationStatus ==
          "PASS_CURRENT_NCBI_ID_SYMBOL",
        "Unique punctuation-only rescue confirmed by current NCBI GeneID/symbol",
        "GeneID validated; current NCBI primary symbol supersedes punctuation candidate symbol"
      )
  )
]

stopifnot(
  nrow(cand_punct) == 197,
  all(
    !is.na(
      cand_punct$CandidateGeneID
    )
  ),
  all(
    !is.na(
      cand_punct$CandidateSymbol
    )
  )
)


# ============================================================
# 6. Candidate identity: 122 current LOC GeneIDs
# ============================================================

cand_loc_current <- loc[
  NCBICurrentPresent == "TRUE",
  .(
    SourceSymbol =
      GeneSymbol_raw,

    CandidateGeneID =
      LOCNumericGeneID,

    CandidateSymbol =
      NCBICurrentSymbol,

    EvidenceClass =
      "LOC_NUMERIC_CURRENT_NCBI",

    EvidenceDetail =
      "LOC numeric suffix is a current NCBI GeneID"
  )
]

stopifnot(
  nrow(cand_loc_current) == 122
)


# ============================================================
# 7. Additional identities from external NCBI Datasets
#
# 71 unique external queries include:
#   12 previously-known alias/punctuation disagreement checks
#      -> ALREADY represented in cand_alias / cand_punct
#
# New source features only:
#   11 LOC redirects
#    2 alias-ambiguous -> unique synonym
#    1 TEC
#   45 unique historical simple-symbol queries
# --------------------------------------------
#   59 genuinely additional source features
# ============================================================

overlap_reasons <- c(
  "SAFE_ALIAS_CURRENT_SYMBOL_DISAGREEMENT",
  "PUNCTUATION_CURRENT_SYMBOL_DISAGREEMENT"
)

ext_new <- ext_unique[
  !Reason %chin%
    overlap_reasons
]

stopifnot(
  nrow(ext_new) == 59
)


ext_new[
  ,
  ExternalEvidenceClass :=
    fcase(

      Reason ==
        "LOC_GENEID_ABSENT_FROM_CURRENT_GENE_INFO" &
      DatasetsClass ==
        "REDIRECTED_TO_CURRENT_GENEID",
      "LOC_GENEID_REDIRECT",

      Reason ==
        "ALIAS_AMBIGUOUS_NOT_CURRENT_NCBI_PRIMARY" &
      DatasetsClass ==
        "UNIQUE_CURRENT_SYNONYM",
      "NCBI_UNIQUE_SYNONYM_RESCUE",

      Reason ==
        "DIRECT_AMBIGUOUS_NOT_UNIQUELY_RESOLVED_BY_NCBI_PRIMARY_SYMBOL" &
      DatasetsClass ==
        "UNIQUE_CURRENT_PRIMARY_SYMBOL",
      "NCBI_PRIMARY_SYMBOL_DISAMBIGUATION",

      Reason ==
        "SIMPLE_SYMBOL_NOT_CURRENT_NCBI_PRIMARY" &
      DatasetsClass ==
        "UNIQUE_RETURNED_GENE_UNCLEAR_MATCH_MODE",
      "NCBI_UNIQUE_SYMBOL_QUERY_RESCUE",

      default =
        "EXTERNAL_UNIQUE_OTHER"
    )
]


cand_external <- ext_new[
  ,
  .(
    SourceSymbol,

    CandidateGeneID =
      ResolvedGeneID,

    CandidateSymbol =
      ResolvedSymbol,

    EvidenceClass =
      ExternalEvidenceClass,

    EvidenceDetail =
      paste(
        Reason,
        DatasetsClass,
        sep = " | "
      )
  )
]


stopifnot(
  nrow(cand_external) == 59,
  all(
    !is.na(
      cand_external$CandidateGeneID
    )
  ),
  all(
    !is.na(
      cand_external$CandidateSymbol
    )
  )
)


# ============================================================
# 8. Combine all candidate identities
# ============================================================

candidates <- rbindlist(
  list(
    cand_direct,
    cand_alias,
    cand_punct,
    cand_loc_current,
    cand_external
  ),
  use.names = TRUE,
  fill = TRUE
)


stopifnot(
  nrow(candidates) == 13334,
  !anyDuplicated(
    candidates$SourceSymbol
  ),
  all(
    candidates$SourceSymbol %chin%
      raw_symbols
  ),
  all(
    !is.na(
      candidates$CandidateGeneID
    )
  ),
  all(
    !is.na(
      candidates$CandidateSymbol
    )
  )
)


cat("============================================================\n")
cat("CANDIDATE IDENTITY ACCOUNTING\n")
cat("============================================================\n")

print(
  candidates[
    ,
    .N,
    by = EvidenceClass
  ][
    order(-N)
  ]
)

cat(
  "\nTotal candidate source features:",
  nrow(candidates),
  "\n"
)


# ============================================================
# 9. Construct the 435 known exclusions
# ============================================================

# ------------------------------------------------------------
# 9A. Composite underscore features: 389
# ------------------------------------------------------------

ex_composite <- u[
  StructuralClass ==
    "COMPOSITE_UNDERSCORE",
  .(
    SourceSymbol =
      GeneSymbol_raw,

    CandidateGeneID =
      NA_character_,

    CandidateSymbol =
      NA_character_,

    EvidenceClass =
      "COMPOSITE_FEATURE",

    Decision =
      "EXCLUDE_COMPOSITE_FEATURE",

    DecisionReason =
      paste0(
        "Historical composite/multi-gene feature; ",
        "not split or assigned to a single canonical gene"
      )
  )
]

stopifnot(
  nrow(ex_composite) == 389
)


# ------------------------------------------------------------
# 9B. Alias target already exists as independent raw row: 6
#
# Expression collision forensic already showed all six pairs
# are DISTINCT_PROFILE.
# ------------------------------------------------------------

ex_target_present <- a[
  PreliminaryIdentityStatus ==
    "ALIAS_TARGET_ALREADY_PRESENT",
  .(
    SourceSymbol =
      GeneSymbol_raw,

    CandidateGeneID =
      CandidateEntrezID,

    CandidateSymbol =
      CandidateCurrentSymbol,

    EvidenceClass =
      "ALIAS_TARGET_ALREADY_PRESENT",

    Decision =
      "EXCLUDE_HISTORICAL_ALIAS_TARGET_PRESENT",

    DecisionReason =
      paste0(
        "Current target already exists as another raw feature; ",
        "expression profiles are distinct; no collapse permitted"
      )
  )
]

stopifnot(
  nrow(ex_target_present) == 6
)


# ------------------------------------------------------------
# 9C. Many-to-one historical aliases without safe choice: 2
# ------------------------------------------------------------

ex_many_old <- a[
  PreliminaryIdentityStatus ==
    "ALIAS_MANY_TO_ONE_COLLISION",
  .(
    SourceSymbol =
      GeneSymbol_raw,

    CandidateGeneID =
      CandidateEntrezID,

    CandidateSymbol =
      CandidateCurrentSymbol,

    EvidenceClass =
      "HISTORICAL_MANY_TO_ONE_COLLISION",

    Decision =
      "EXCLUDE_HISTORICAL_MANY_TO_ONE",

    DecisionReason =
      paste0(
        "Multiple distinct historical expression rows map to ",
        "one modern identity; no current-primary source row ",
        "provides a defensible row-level collapse"
      )
  )
]

stopifnot(
  nrow(ex_many_old) == 2
)


# ------------------------------------------------------------
# 9D. External NCBI no live symbol match: 12
# ------------------------------------------------------------

ex_no_live <- ext_unresolved[
  DatasetsClass ==
    "NO_LIVE_SYMBOL_MATCH",
  .(
    SourceSymbol,

    CandidateGeneID =
      NA_character_,

    CandidateSymbol =
      NA_character_,

    EvidenceClass =
      "NO_LIVE_NCBI_SYMBOL_MATCH",

    Decision =
      "EXCLUDE_NO_LIVE_CURRENT_IDENTITY",

    DecisionReason =
      "No unique live current NCBI Gene identity recovered"
  )
]

stopifnot(
  nrow(ex_no_live) == 12
)


# ------------------------------------------------------------
# 9E. External NCBI ambiguous multiple current genes: 26
# ------------------------------------------------------------

ex_ambiguous <- ext_unresolved[
  DatasetsClass ==
    "AMBIGUOUS_MULTIPLE_CURRENT_GENES",
  .(
    SourceSymbol,

    CandidateGeneID =
      ReturnedGeneIDs,

    CandidateSymbol =
      ReturnedSymbols,

    EvidenceClass =
      "AMBIGUOUS_MULTIPLE_CURRENT_GENES",

    Decision =
      "EXCLUDE_AMBIGUOUS_CURRENT_IDENTITY",

    DecisionReason =
      paste0(
        "Historical symbol maps to multiple current GeneIDs; ",
        "no arbitrary current identity selected"
      )
  )
]

stopifnot(
  nrow(ex_ambiguous) == 26
)


# ============================================================
# 10. Known exclusion accounting
# ============================================================

known_exclusions <- rbindlist(
  list(
    ex_composite,
    ex_target_present,
    ex_many_old,
    ex_no_live,
    ex_ambiguous
  ),
  use.names = TRUE,
  fill = TRUE
)


stopifnot(
  nrow(known_exclusions) == 435,
  !anyDuplicated(
    known_exclusions$SourceSymbol
  ),
  !any(
    known_exclusions$SourceSymbol %chin%
      candidates$SourceSymbol
  )
)


# ------------------------------------------------------------
# Full raw-feature accounting BEFORE final candidate collision
# audit.
# ------------------------------------------------------------

accounted_symbols <- c(
  candidates$SourceSymbol,
  known_exclusions$SourceSymbol
)

stopifnot(
  length(accounted_symbols) == 13769,
  uniqueN(accounted_symbols) == 13769,
  setequal(
    accounted_symbols,
    raw_symbols
  )
)


cat("\n============================================================\n")
cat("PRE-COLLISION FULL FEATURE ACCOUNTING\n")
cat("============================================================\n")

cat(
  "Candidate identities:",
  nrow(candidates),
  "\n"
)

cat(
  "Known exclusions:",
  nrow(known_exclusions),
  "\n"
)

cat(
  "Total:",
  nrow(candidates) +
    nrow(known_exclusions),
  "\n"
)

cat(
  "[PASS] All 13,769 raw features are accounted for\n"
)


# ============================================================
# 11. GLOBAL FinalGeneID collision audit
# ============================================================

gid_groups <- candidates[
  ,
  .(
    SourceFeatureN =
      .N,

    SourceSymbols =
      paste(
        SourceSymbol,
        collapse = "|"
      ),

    PrimarySourceN =
      sum(
        SourceSymbol ==
          CandidateSymbol
      ),

    PrimarySourceSymbols =
      paste(
        SourceSymbol[
          SourceSymbol ==
            CandidateSymbol
        ],
        collapse = "|"
      )
  ),
  by = .(
    CandidateGeneID,
    CandidateSymbol
  )
]


gid_idx <- match(
  candidates$CandidateGeneID,
  gid_groups$CandidateGeneID
)


candidates[
  ,
  CandidateGeneIDSourceN :=
    gid_groups$SourceFeatureN[
      gid_idx
    ]
]

candidates[
  ,
  CurrentPrimarySourceN :=
    gid_groups$PrimarySourceN[
      gid_idx
    ]
]

candidates[
  ,
  CollisionSourceSymbols :=
    gid_groups$SourceSymbols[
      gid_idx
    ]
]


# ============================================================
# 12. Preliminary candidate decisions
# ============================================================

candidates[
  ,
  Decision := fcase(

    CandidateGeneIDSourceN == 1,
    "KEEP_UNIQUE_CURRENT_IDENTITY",

    CandidateGeneIDSourceN > 1 &
      CurrentPrimarySourceN >= 1 &
      SourceSymbol ==
        CandidateSymbol,
    "KEEP_CURRENT_PRIMARY_COLLISION_GROUP",

    CandidateGeneIDSourceN > 1 &
      CurrentPrimarySourceN >= 1 &
      SourceSymbol !=
        CandidateSymbol,
    "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",

    CandidateGeneIDSourceN > 1 &
      CurrentPrimarySourceN == 0,
    "REVIEW_MANY_TO_ONE_NO_PRIMARY",

    default =
      "REVIEW_UNCLASSIFIED"
  )
]


candidates[
  ,
  DecisionReason := fcase(

    Decision ==
      "KEEP_UNIQUE_CURRENT_IDENTITY",
    "Unique current GeneID candidate",

    Decision ==
      "KEEP_CURRENT_PRIMARY_COLLISION_GROUP",
    paste0(
      "Multiple source rows map to same current GeneID; ",
      "retain the raw row whose symbol equals current primary symbol"
    ),

    Decision ==
      "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",
    paste0(
      "Historical/non-primary source maps to a current GeneID ",
      "whose current-primary raw row is already present; no collapse"
    ),

    Decision ==
      "REVIEW_MANY_TO_ONE_NO_PRIMARY",
    paste0(
      "Multiple source rows map to same current GeneID but no ",
      "current-primary raw row is present; arbitrary selection prohibited"
    ),

    default =
      "Unclassified candidate collision"
  )
]


# ============================================================
# 13. GeneID collision report
# ============================================================

gid_collision_rows <- candidates[
  CandidateGeneIDSourceN > 1
]


cat("\n============================================================\n")
cat("GLOBAL GeneID COLLISION AUDIT\n")
cat("============================================================\n")

cat(
  "Candidate rows involved in GeneID collisions:",
  nrow(gid_collision_rows),
  "\n"
)

cat(
  "Unique collided GeneIDs:",
  uniqueN(
    gid_collision_rows$CandidateGeneID
  ),
  "\n"
)

if (nrow(gid_collision_rows) > 0) {

  print(
    gid_collision_rows[
      ,
      .(
        SourceSymbol,
        CandidateGeneID,
        CandidateSymbol,
        EvidenceClass,
        CandidateGeneIDSourceN,
        CurrentPrimarySourceN,
        CollisionSourceSymbols,
        Decision
      )
    ][
      order(
        CandidateGeneID,
        SourceSymbol
      )
    ]
  )
}


# ============================================================
# 14. Final-symbol collision audit among provisionally kept rows
#
# A current symbol must not refer to >1 distinct GeneID in
# the final canonical expression matrix.
# ============================================================

keep_decisions <- c(
  "KEEP_UNIQUE_CURRENT_IDENTITY",
  "KEEP_CURRENT_PRIMARY_COLLISION_GROUP"
)


provisional_keep <- candidates[
  Decision %chin%
    keep_decisions
]


symbol_groups <- provisional_keep[
  ,
  .(
    DistinctGeneIDN =
      uniqueN(
        CandidateGeneID
      ),

    GeneIDs =
      paste(
        sort(
          unique(
            CandidateGeneID
          )
        ),
        collapse = "|"
      ),

    SourceSymbols =
      paste(
        SourceSymbol,
        collapse = "|"
      )
  ),
  by = CandidateSymbol
]


symbol_problem <- symbol_groups[
  DistinctGeneIDN > 1
]


if (nrow(symbol_problem) > 0) {

  problem_symbols <-
    symbol_problem$CandidateSymbol

  candidates[
    CandidateSymbol %chin%
      problem_symbols &
    Decision %chin%
      keep_decisions,
    `:=`(
      Decision =
        "REVIEW_FINAL_SYMBOL_COLLISION",

      DecisionReason =
        paste0(
          "Same proposed current symbol corresponds to multiple ",
          "distinct GeneIDs; requires explicit review"
        )
    )
  ]
}


cat("\n============================================================\n")
cat("FINAL SYMBOL COLLISION AUDIT\n")
cat("============================================================\n")

cat(
  "Current symbols mapping to >1 retained GeneID:",
  nrow(symbol_problem),
  "\n"
)

if (nrow(symbol_problem) > 0) {

  print(
    symbol_problem
  )
}


# ============================================================
# 15. Build full decision table in authoritative feature order
# ============================================================

candidate_decisions <- candidates[
  ,
  .(
    SourceSymbol,

    CandidateGeneID,
    CandidateSymbol,

    EvidenceClass,
    EvidenceDetail,

    Decision,
    DecisionReason,

    CandidateGeneIDSourceN,
    CurrentPrimarySourceN,
    CollisionSourceSymbols
  )
]


exclusion_decisions <- known_exclusions[
  ,
  .(
    SourceSymbol,

    CandidateGeneID,
    CandidateSymbol,

    EvidenceClass,

    EvidenceDetail =
      NA_character_,

    Decision,
    DecisionReason,

    CandidateGeneIDSourceN =
      NA_integer_,

    CurrentPrimarySourceN =
      NA_integer_,

    CollisionSourceSymbols =
      NA_character_
  )
]


decision <- rbindlist(
  list(
    candidate_decisions,
    exclusion_decisions
  ),
  use.names = TRUE,
  fill = TRUE
)


stopifnot(
  nrow(decision) == 13769,
  !anyDuplicated(
    decision$SourceSymbol
  ),
  setequal(
    decision$SourceSymbol,
    raw_symbols
  )
)


decision <- decision[
  match(
    raw_symbols,
    SourceSymbol
  )
]


stopifnot(
  identical(
    decision$SourceSymbol,
    raw_symbols
  )
)


decision[
  ,
  FeatureOrder :=
    seq_len(.N)
]


# ============================================================
# 16. Canonical inclusion flag
# ============================================================

decision[
  ,
  IncludeCanonical :=
    Decision %chin%
      keep_decisions
]


# ------------------------------------------------------------
# Final identity populated ONLY for included rows.
# ------------------------------------------------------------

decision[
  ,
  FinalGeneID :=
    fifelse(
      IncludeCanonical,
      CandidateGeneID,
      NA_character_
    )
]

decision[
  ,
  FinalSymbol :=
    fifelse(
      IncludeCanonical,
      CandidateSymbol,
      NA_character_
    )
]


# ============================================================
# 17. Final hard uniqueness checks on INCLUDED features
# ============================================================

included <- decision[
  IncludeCanonical == TRUE
]

review <- decision[
  grepl(
    "^REVIEW_",
    Decision
  )
]


included_geneid_unique <-
  !anyDuplicated(
    included$FinalGeneID
  )

included_symbol_unique <-
  !anyDuplicated(
    included$FinalSymbol
  )


cat("\n============================================================\n")
cat("FINAL DECISION SUMMARY\n")
cat("============================================================\n")

decision_summary <- decision[
  ,
  .N,
  by = .(
    IncludeCanonical,
    Decision
  )
][
  order(
    -IncludeCanonical,
    -N,
    Decision
  )
]

print(
  decision_summary
)


cat(
  "\nIncluded candidate rows:",
  nrow(included),
  "\n"
)

cat(
  "Excluded/review rows:",
  nrow(decision) -
    nrow(included),
  "\n"
)

cat(
  "Review-required rows:",
  nrow(review),
  "\n"
)

cat(
  "Unique included FinalGeneID:",
  uniqueN(
    included$FinalGeneID
  ),
  "\n"
)

cat(
  "Unique included FinalSymbol:",
  uniqueN(
    included$FinalSymbol
  ),
  "\n"
)


# ============================================================
# 18. Evidence summary for included features
# ============================================================

evidence_summary <- included[
  ,
  .N,
  by = EvidenceClass
][
  order(-N)
]


cat("\n============================================================\n")
cat("INCLUDED EVIDENCE SUMMARY\n")
cat("============================================================\n")

print(
  evidence_summary
)


# ============================================================
# 19. Save
# ============================================================

setcolorder(
  decision,
  c(
    "FeatureOrder",
    "SourceSymbol",
    "EvidenceClass",
    "EvidenceDetail",
    "CandidateGeneID",
    "CandidateSymbol",
    "CandidateGeneIDSourceN",
    "CurrentPrimarySourceN",
    "CollisionSourceSymbols",
    "Decision",
    "DecisionReason",
    "IncludeCanonical",
    "FinalGeneID",
    "FinalSymbol"
  )
)


fwrite(
  decision,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_decision_table.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  included,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_included_features.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  review,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_review_required.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  gid_collision_rows,
  file.path(
    harm_dir,
    "GSE93624_final_GeneID_collision_audit.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  symbol_problem,
  file.path(
    harm_dir,
    "GSE93624_final_symbol_collision_audit.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  decision_summary,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_decision_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  evidence_summary,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_evidence_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 20. Final status
# ============================================================

cat("\n============================================================\n")
cat("INTEGRITY CHECKS\n")
cat("============================================================\n")

cat("[PASS] 13,769 raw features accounted for exactly once\n")

if (
  nrow(review) == 0 &&
  included_geneid_unique &&
  included_symbol_unique
) {

  cat("[PASS] No unresolved candidate collision remains\n")
  cat("[PASS] Included FinalGeneID values are unique\n")
  cat("[PASS] Included FinalSymbol values are unique\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE93624_CANONICAL_GENE_DECISION_TABLE_BUILT\n")

} else {

  if (nrow(review) > 0) {

    cat(
      "[REVIEW] Candidate rows requiring review:",
      nrow(review),
      "\n"
    )
  }

  if (!included_geneid_unique) {

    cat(
      "[REVIEW] Duplicate FinalGeneID remains among included rows\n"
    )
  }

  if (!included_symbol_unique) {

    cat(
      "[REVIEW] Duplicate FinalSymbol remains among included rows\n"
    )
  }

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_REQUIRED_GSE93624_CANONICAL_GENE_DECISION_TABLE\n")
}

cat("============================================================\n")
