suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 2F-final
#
# Repair/finalize canonical gene decision table.
#
# WHY:
# Previous Stage 09 grouped candidate collisions by:
#   CandidateGeneID + CandidateSymbol
#
# This can miss:
#   same GeneID
#   + different historical/provisional symbols
#
# FINAL RULE:
#   GeneID is the identity anchor.
#   Current NCBI primary symbol is annotation/display identity.
#
# Therefore:
#   1. reconcile EVERY candidate GeneID against current NCBI
#   2. replace provisional symbol with current NCBI symbol
#   3. audit collisions by GeneID ONLY
#   4. if exactly one current-primary raw source exists:
#        keep it, exclude all non-primary historical rows
#   5. if no primary raw source exists:
#        REVIEW; do not arbitrarily select a row
#
# No expression values are modified.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

harm_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

old_decision_file <- file.path(
  harm_dir,
  "GSE93624_canonical_gene_decision_table.tsv"
)

gene_info_file <- paste0(
  "/home/mazekai/IBD_EcoTyper/06_validation/99_common/",
  "gene_annotation/NCBI_Gene_2026-09-21/",
  "Homo_sapiens.gene_info.gz"
)


cat("============================================================\n")
cat("GSE93624 FINAL GeneID-only canonical reconciliation\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load Stage 09 decision table
# ============================================================

d <- fread(
  old_decision_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  nrow(d) == 13769,
  !anyDuplicated(
    d$SourceSymbol
  )
)

cat(
  "[PASS] Loaded Stage 09 decision table:",
  nrow(d),
  "raw features\n"
)


# ============================================================
# 2. Identify the OLD residual duplicate GeneID
# ============================================================

old_included <- d[
  IncludeCanonical == "TRUE"
]

old_dup_geneids <- old_included[
  !is.na(FinalGeneID) &
  (
    duplicated(FinalGeneID) |
    duplicated(
      FinalGeneID,
      fromLast = TRUE
    )
  ),
  unique(FinalGeneID)
]


cat("\n============================================================\n")
cat("OLD STAGE 09 RESIDUAL GeneID DUPLICATE\n")
cat("============================================================\n")

cat(
  "Duplicated included GeneIDs:",
  length(old_dup_geneids),
  "\n"
)

if (length(old_dup_geneids) > 0) {

  old_dup_rows <- old_included[
    FinalGeneID %chin%
      old_dup_geneids,
    .(
      SourceSymbol,
      EvidenceClass,
      CandidateGeneID,
      CandidateSymbol,
      Decision,
      FinalGeneID,
      FinalSymbol
    )
  ][
    order(
      FinalGeneID,
      SourceSymbol
    )
  ]

  print(
    old_dup_rows
  )

} else {

  old_dup_rows <- data.table()
}


stopifnot(
  length(old_dup_geneids) == 1
)


# ============================================================
# 3. Define the 13,334 candidate-source rows
#
# These are the rows that participated in Stage 09 candidate
# collision reconciliation.
#
# Important:
# EXCLUDE_NONPRIMARY_COLLISION_SOURCE is still a valid
# candidate identity; it was excluded only because another
# source row represented the same GeneID.
# ============================================================

candidate_decision_classes <- c(
  "KEEP_UNIQUE_CURRENT_IDENTITY",
  "KEEP_CURRENT_PRIMARY_COLLISION_GROUP",
  "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",
  "REVIEW_MANY_TO_ONE_NO_PRIMARY",
  "REVIEW_UNCLASSIFIED",
  "REVIEW_FINAL_SYMBOL_COLLISION"
)


candidate <- copy(
  d[
    Decision %chin%
      candidate_decision_classes
  ]
)

stopifnot(
  nrow(candidate) == 13334,
  all(
    !is.na(
      candidate$CandidateGeneID
    )
  ),
  all(
    candidate$CandidateGeneID != ""
  ),
  !any(
    grepl(
      "\\|",
      candidate$CandidateGeneID
    )
  )
)

cat(
  "\n[PASS] Candidate identity rows:",
  nrow(candidate),
  "\n"
)


# ============================================================
# 4. Known non-candidate exclusions
# ============================================================

known_exclusion <- copy(
  d[
    !Decision %chin%
      candidate_decision_classes
  ]
)

stopifnot(
  nrow(known_exclusion) == 435
)

cat(
  "[PASS] Pre-existing non-candidate exclusions:",
  nrow(known_exclusion),
  "\n"
)


# ============================================================
# 5. Load current NCBI human Gene reference
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
  !anyDuplicated(
    gi$NCBI_GeneID
  )
)

cat(
  "[PASS] Current NCBI human Gene records:",
  nrow(gi),
  "\n"
)


# ============================================================
# 6. Reconcile ALL candidate GeneIDs against current NCBI
#
# This is the key correction.
# CandidateSymbol is NOT used as part of identity grouping.
# ============================================================

idx <- match(
  candidate$CandidateGeneID,
  gi$NCBI_GeneID
)

candidate[
  ,
  GeneIDPresentCurrentNCBI :=
    !is.na(idx)
]

candidate[
  ,
  NCBI_CurrentSymbol :=
    fifelse(
      !is.na(idx),
      gi$NCBI_CurrentSymbol[idx],
      NA_character_
    )
]

candidate[
  ,
  NCBI_GeneType :=
    fifelse(
      !is.na(idx),
      gi$NCBI_GeneType[idx],
      NA_character_
    )
]

candidate[
  ,
  NCBI_Description :=
    fifelse(
      !is.na(idx),
      gi$NCBI_Description[idx],
      NA_character_
    )
]


# ------------------------------------------------------------
# At this point all final candidate GeneIDs should be current.
#
# Historical LOC IDs were already redirected before Stage 09.
# ------------------------------------------------------------

missing_current <- candidate[
  GeneIDPresentCurrentNCBI == FALSE
]

cat("\n============================================================\n")
cat("CURRENT NCBI GeneID COVERAGE\n")
cat("============================================================\n")

cat(
  "Candidate rows absent from current gene_info:",
  nrow(missing_current),
  "\n"
)

if (nrow(missing_current) > 0) {

  print(
    missing_current[
      ,
      .(
        SourceSymbol,
        CandidateGeneID,
        CandidateSymbol,
        EvidenceClass
      )
    ]
  )

  stop(
    "Not all candidate GeneIDs are present in current NCBI gene_info."
  )
}


# ============================================================
# 7. Audit provisional-symbol updates
# ============================================================

candidate[
  ,
  SymbolChangedByFinalNCBI :=
    !is.na(CandidateSymbol) &
    CandidateSymbol !=
      NCBI_CurrentSymbol
]


symbol_updates <- candidate[
  SymbolChangedByFinalNCBI == TRUE
]


cat("\n============================================================\n")
cat("FINAL CURRENT-SYMBOL RECONCILIATION\n")
cat("============================================================\n")

cat(
  "Rows whose provisional CandidateSymbol differs from",
  " current NCBI primary symbol:",
  nrow(symbol_updates),
  "\n"
)

if (nrow(symbol_updates) > 0) {

  print(
    symbol_updates[
      ,
      .(
        SourceSymbol,
        CandidateGeneID,
        CandidateSymbol_old =
          CandidateSymbol,
        NCBI_CurrentSymbol,
        EvidenceClass
      )
    ]
  )
}


# ============================================================
# 8. GLOBAL collision audit BY GeneID ONLY
#
# This is the corrected grouping key.
# ============================================================

gid_group <- candidate[
  ,
  .(
    SourceFeatureN =
      .N,

    SourceSymbols =
      paste(
        SourceSymbol,
        collapse = "|"
      ),

    CurrentPrimarySourceN =
      sum(
        SourceSymbol ==
          NCBI_CurrentSymbol
      ),

    CurrentPrimarySourceSymbols =
      paste(
        SourceSymbol[
          SourceSymbol ==
            NCBI_CurrentSymbol
        ],
        collapse = "|"
      )
  ),
  by = .(
    CandidateGeneID
  )
]


idx_group <- match(
  candidate$CandidateGeneID,
  gid_group$CandidateGeneID
)

candidate[
  ,
  GeneIDSourceN :=
    gid_group$SourceFeatureN[
      idx_group
    ]
]

candidate[
  ,
  CurrentPrimarySourceN :=
    gid_group$CurrentPrimarySourceN[
      idx_group
    ]
]

candidate[
  ,
  CollisionSourceSymbols :=
    gid_group$SourceSymbols[
      idx_group
    ]
]


# ============================================================
# 9. Final candidate decision
# ============================================================

candidate[
  ,
  Decision_FINAL := fcase(

    GeneIDSourceN == 1,
    "KEEP_UNIQUE_CURRENT_IDENTITY",

    GeneIDSourceN > 1 &
      CurrentPrimarySourceN == 1 &
      SourceSymbol ==
        NCBI_CurrentSymbol,
    "KEEP_CURRENT_PRIMARY_COLLISION_GROUP",

    GeneIDSourceN > 1 &
      CurrentPrimarySourceN == 1 &
      SourceSymbol !=
        NCBI_CurrentSymbol,
    "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",

    GeneIDSourceN > 1 &
      CurrentPrimarySourceN == 0,
    "REVIEW_MANY_TO_ONE_NO_PRIMARY",

    GeneIDSourceN > 1 &
      CurrentPrimarySourceN > 1,
    "REVIEW_MULTIPLE_CURRENT_PRIMARY_SOURCES",

    default =
      "REVIEW_UNCLASSIFIED"
  )
]


candidate[
  ,
  DecisionReason_FINAL := fcase(

    Decision_FINAL ==
      "KEEP_UNIQUE_CURRENT_IDENTITY",
    "Only one source expression row maps to this current GeneID",

    Decision_FINAL ==
      "KEEP_CURRENT_PRIMARY_COLLISION_GROUP",
    paste0(
      "Multiple source rows map to this current GeneID; ",
      "retain the deposited row whose source symbol equals ",
      "the current NCBI primary symbol"
    ),

    Decision_FINAL ==
      "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",
    paste0(
      "Another deposited source row represents the same ",
      "current GeneID using the current NCBI primary symbol; ",
      "historical/non-primary row excluded without aggregation"
    ),

    Decision_FINAL ==
      "REVIEW_MANY_TO_ONE_NO_PRIMARY",
    paste0(
      "Multiple source expression rows map to the same current ",
      "GeneID but none uses the current primary symbol; ",
      "arbitrary selection prohibited"
    ),

    Decision_FINAL ==
      "REVIEW_MULTIPLE_CURRENT_PRIMARY_SOURCES",
    paste0(
      "More than one source row appears to use the current ",
      "primary symbol for the same GeneID"
    ),

    default =
      "Unclassified GeneID collision"
  )
]


# ============================================================
# 10. Corrected collision report
# ============================================================

collision <- candidate[
  GeneIDSourceN > 1
][
  order(
    CandidateGeneID,
    SourceSymbol
  )
]


cat("\n============================================================\n")
cat("CORRECTED GLOBAL GeneID-ONLY COLLISION AUDIT\n")
cat("============================================================\n")

cat(
  "Candidate rows involved:",
  nrow(collision),
  "\n"
)

cat(
  "Unique collided GeneIDs:",
  uniqueN(
    collision$CandidateGeneID
  ),
  "\n"
)

if (nrow(collision) > 0) {

  print(
    collision[
      ,
      .(
        SourceSymbol,
        CandidateGeneID,
        OldCandidateSymbol =
          CandidateSymbol,
        NCBI_CurrentSymbol,
        EvidenceClass,
        GeneIDSourceN,
        CurrentPrimarySourceN,
        CollisionSourceSymbols,
        Decision_FINAL
      )
    ]
  )
}


# ============================================================
# 11. Review rows
# ============================================================

review_classes <- c(
  "REVIEW_MANY_TO_ONE_NO_PRIMARY",
  "REVIEW_MULTIPLE_CURRENT_PRIMARY_SOURCES",
  "REVIEW_UNCLASSIFIED"
)

candidate_review <- candidate[
  Decision_FINAL %chin%
    review_classes
]


cat("\n============================================================\n")
cat("CANDIDATE REVIEW ROWS\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(candidate_review),
  "\n"
)

if (nrow(candidate_review) > 0) {

  print(
    candidate_review[
      ,
      .(
        SourceSymbol,
        CandidateGeneID,
        NCBI_CurrentSymbol,
        EvidenceClass,
        GeneIDSourceN,
        CurrentPrimarySourceN,
        CollisionSourceSymbols,
        Decision_FINAL
      )
    ]
  )
}


# ============================================================
# 12. Build final candidate rows
# ============================================================

keep_classes <- c(
  "KEEP_UNIQUE_CURRENT_IDENTITY",
  "KEEP_CURRENT_PRIMARY_COLLISION_GROUP"
)


candidate[
  ,
  IncludeCanonical_FINAL :=
    Decision_FINAL %chin%
      keep_classes
]


candidate[
  ,
  FinalGeneID_FINAL :=
    fifelse(
      IncludeCanonical_FINAL,
      CandidateGeneID,
      NA_character_
    )
]


candidate[
  ,
  FinalSymbol_FINAL :=
    fifelse(
      IncludeCanonical_FINAL,
      NCBI_CurrentSymbol,
      NA_character_
    )
]


# ============================================================
# 13. Rebuild full 13,769-row final decision table
# ============================================================

candidate_final <- candidate[
  ,
  .(
    FeatureOrder,
    SourceSymbol,
    EvidenceClass,
    EvidenceDetail,

    CandidateGeneID,

    CandidateSymbol_original =
      CandidateSymbol,

    NCBI_CurrentSymbol,

    NCBI_GeneType,
    NCBI_Description,

    GeneIDSourceN,
    CurrentPrimarySourceN,
    CollisionSourceSymbols,

    SymbolChangedByFinalNCBI,

    Decision =
      Decision_FINAL,

    DecisionReason =
      DecisionReason_FINAL,

    IncludeCanonical =
      IncludeCanonical_FINAL,

    FinalGeneID =
      FinalGeneID_FINAL,

    FinalSymbol =
      FinalSymbol_FINAL
  )
]


# ------------------------------------------------------------
# Preserve already-locked 435 non-candidate exclusions.
# ------------------------------------------------------------

known_exclusion_final <- known_exclusion[
  ,
  .(
    FeatureOrder,
    SourceSymbol,
    EvidenceClass,
    EvidenceDetail,

    CandidateGeneID,

    CandidateSymbol_original =
      CandidateSymbol,

    NCBI_CurrentSymbol =
      NA_character_,

    NCBI_GeneType =
      NA_character_,

    NCBI_Description =
      NA_character_,

    GeneIDSourceN =
      NA_integer_,

    CurrentPrimarySourceN =
      NA_integer_,

    CollisionSourceSymbols =
      NA_character_,

    SymbolChangedByFinalNCBI =
      NA,

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
    known_exclusion_final
  ),
  use.names = TRUE,
  fill = TRUE
)


final[
  ,
  FeatureOrder :=
    as.integer(
      FeatureOrder
    )
]

setorder(
  final,
  FeatureOrder
)


stopifnot(
  nrow(final) == 13769,
  !anyDuplicated(
    final$SourceSymbol
  ),
  identical(
    final$FeatureOrder,
    1:13769
  )
)


# ============================================================
# 14. FINAL uniqueness validation
# ============================================================

included <- final[
  IncludeCanonical == TRUE
]

excluded <- final[
  IncludeCanonical == FALSE
]

review <- final[
  grepl(
    "^REVIEW_",
    Decision
  )
]


dup_geneid_final <- included[
  duplicated(FinalGeneID) |
  duplicated(
    FinalGeneID,
    fromLast = TRUE
  )
]

dup_symbol_final <- included[
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
    -N,
    Decision
  )
]

print(
  decision_summary
)


cat(
  "\nIncluded rows:",
  nrow(included),
  "\n"
)

cat(
  "Excluded rows:",
  nrow(excluded),
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
# 15. Evidence summary
# ============================================================

evidence_summary <- included[
  ,
  .N,
  by = EvidenceClass
][
  order(-N)
]


cat("\n============================================================\n")
cat("FINAL INCLUDED EVIDENCE SUMMARY\n")
cat("============================================================\n")

print(
  evidence_summary
)


# ============================================================
# 16. Save authoritative FINAL outputs
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
    "GSE93624_final_GeneID_ONLY_collision_audit.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  symbol_updates,
  file.path(
    harm_dir,
    "GSE93624_final_NCBI_symbol_updates.tsv"
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

fwrite(
  evidence_summary,
  file.path(
    harm_dir,
    "GSE93624_canonical_gene_evidence_summary_FINAL.tsv"
  ),
  sep = "\t"
)

fwrite(
  old_dup_rows,
  file.path(
    harm_dir,
    "GSE93624_stage09_residual_duplicate_GeneID.tsv"
  ),
  sep = "\t",
  na = ""
)


# ============================================================
# 17. Final hard guards
# ============================================================

cat("\n============================================================\n")
cat("FINAL INTEGRITY CHECKS\n")
cat("============================================================\n")

cat("[PASS] 13,769 raw features accounted for exactly once\n")


if (nrow(review) > 0) {

  cat(
    "[REVIEW] Rows requiring manual review:",
    nrow(review),
    "\n"
  )
}

if (nrow(dup_geneid_final) > 0) {

  cat(
    "[FAIL] Duplicate FinalGeneID remains:",
    uniqueN(
      dup_geneid_final$FinalGeneID
    ),
    "\n"
  )

  print(
    dup_geneid_final
  )
}

if (nrow(dup_symbol_final) > 0) {

  cat(
    "[FAIL] Duplicate FinalSymbol remains:",
    uniqueN(
      dup_symbol_final$FinalSymbol
    ),
    "\n"
  )

  print(
    dup_symbol_final
  )
}


if (
  nrow(review) == 0 &&
  nrow(dup_geneid_final) == 0 &&
  nrow(dup_symbol_final) == 0 &&
  nrow(included) ==
    uniqueN(included$FinalGeneID) &&
  nrow(included) ==
    uniqueN(included$FinalSymbol)
) {

  cat("[PASS] No candidate review rows remain\n")
  cat("[PASS] Included FinalGeneID values are unique\n")
  cat("[PASS] Included FinalSymbol values are unique\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE93624_FINAL_CANONICAL_GENE_DECISION_TABLE\n")

} else {

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_REQUIRED_GSE93624_FINAL_CANONICAL_GENE_DECISION_TABLE\n")
}

cat("============================================================\n")
