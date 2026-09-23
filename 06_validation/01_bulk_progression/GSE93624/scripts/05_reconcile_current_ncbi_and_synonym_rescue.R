suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 2D
#
# Comprehensive current-NCBI identity reconciliation
#
# Goals:
#   1. Reconcile all already-known candidate GeneIDs against
#      current NCBI primary symbols.
#   2. Resolve org.Hs / NCBI symbol disagreements by GeneID.
#   3. Search current NCBI gene_info Synonyms for:
#         - 57 remaining simple symbols
#         - 28 alias-ambiguous symbols
#   4. Perform GLOBAL GeneID / symbol collision audit.
#   5. Generate only the truly residual external query set.
#
# IMPORTANT:
#   - No expression values are changed.
#   - No expression rows are renamed.
#   - No rows are merged or excluded yet.
#   - This is still an identity forensic stage.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

expr_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE93624_authoritative_expression_245.rds"
)

audit_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE93624_gene_symbol_harmonization_audit.tsv"
)

unresolved_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE93624_unresolved_symbol_structural_audit.tsv"
)

gene_info_file <- paste0(
  "/home/mazekai/IBD_EcoTyper/06_validation/99_common/",
  "gene_annotation/NCBI_Gene_2026-09-21/",
  "Homo_sapiens.gene_info.gz"
)

out_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

cat("============================================================\n")
cat("GSE93624 current-NCBI reconciliation + synonym rescue\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load resources
# ============================================================

obj <- readRDS(expr_file)
expr <- obj$expression

a <- fread(
  audit_file,
  na.strings = c("", "NA")
)

u <- fread(
  unresolved_file,
  na.strings = c("", "NA")
)

stopifnot(
  nrow(expr) == 13769,
  nrow(a) == 13769,
  nrow(u) == 776,
  identical(
    rownames(expr),
    a$GeneSymbol_raw
  )
)

raw_symbols <- rownames(expr)

cat("[PASS] Expression and forensic annotation resources aligned\n")


# ============================================================
# 2. Load current NCBI human gene_info
# ============================================================

gi <- fread(
  gene_info_file,
  header = TRUE,
  sep = "\t",
  select = c(
    "GeneID",
    "Symbol",
    "Synonyms",
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
    "Synonyms",
    "description",
    "type_of_gene"
  ),
  c(
    "NCBI_GeneID",
    "NCBI_Symbol",
    "NCBI_Synonyms",
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
  "Current NCBI human Gene records:",
  nrow(gi),
  "\n\n"
)


# ============================================================
# 3. Build existing candidate-identity table
# ============================================================

# ------------------------------------------------------------
# 3A. Direct current symbols
# ------------------------------------------------------------

cand_direct <- a[
  PreliminaryIdentityStatus ==
    "CURRENT_IDENTITY_OK",
  .(
    SourceSymbol =
      GeneSymbol_raw,

    CandidateGeneID =
      Direct_EntrezIDs,

    PreNCBISymbol =
      GeneSymbol_raw,

    EvidenceClass =
      "ORGHS_DIRECT_CURRENT"
  )
]

stopifnot(
  nrow(cand_direct) == 12208
)


# ------------------------------------------------------------
# 3B. Safe org.Hs alias candidates
# ------------------------------------------------------------

cand_alias <- a[
  PreliminaryIdentityStatus ==
    "SAFE_ALIAS_CANDIDATE",
  .(
    SourceSymbol =
      GeneSymbol_raw,

    CandidateGeneID =
      CandidateEntrezID,

    PreNCBISymbol =
      CandidateCurrentSymbol,

    EvidenceClass =
      "ORGHS_SAFE_ALIAS"
  )
]

stopifnot(
  nrow(cand_alias) == 748
)


# ------------------------------------------------------------
# 3C. Safe punctuation candidates
# ------------------------------------------------------------

cand_punct <- u[
  PunctuationRescueStatus ==
    "SAFE_PUNCTUATION_CANDIDATE",
  .(
    SourceSymbol =
      GeneSymbol_raw,

    CandidateGeneID =
      Punct_EntrezIDs,

    PreNCBISymbol =
      Punct_CurrentSymbols,

    EvidenceClass =
      "PUNCTUATION_RESCUE"
  )
]

stopifnot(
  nrow(cand_punct) == 197
)


# ------------------------------------------------------------
# 3D. Current LOC numeric GeneIDs
# ------------------------------------------------------------

loc_all <- u[
  StructuralClass ==
    "LOC_STYLE",
  .(
    SourceSymbol =
      GeneSymbol_raw
  )
]

stopifnot(
  nrow(loc_all) == 133
)

loc_all[
  ,
  CandidateGeneID :=
    sub(
      "^LOC",
      "",
      SourceSymbol
    )
]

loc_idx <- match(
  loc_all$CandidateGeneID,
  gi$NCBI_GeneID
)

loc_current <- loc_all[
  !is.na(loc_idx)
]

loc_current_idx <- match(
  loc_current$CandidateGeneID,
  gi$NCBI_GeneID
)

loc_current[
  ,
  PreNCBISymbol :=
    gi$NCBI_Symbol[
      loc_current_idx
    ]
]

loc_current[
  ,
  EvidenceClass :=
    "LOC_NUMERIC_CURRENT_GENEID"
]

stopifnot(
  nrow(loc_current) == 122
)


# ============================================================
# 4. NCBI synonym rescue for remaining 85 symbol problems
# ============================================================

simple_remaining <- u[
  StructuralClass ==
    "SIMPLE_SYMBOL_LIKE" &
  PunctuationRescueStatus ==
    "NO_PUNCTUATION_RESCUE",
  GeneSymbol_raw
]

alias_ambiguous <- a[
  PreliminaryIdentityStatus ==
    "ALIAS_AMBIGUOUS",
  GeneSymbol_raw
]

stopifnot(
  length(simple_remaining) == 57,
  length(alias_ambiguous) == 28
)

symbol_queries <- unique(
  c(
    simple_remaining,
    alias_ambiguous
  )
)

stopifnot(
  length(symbol_queries) == 85
)

cat(
  "Symbols entering local NCBI synonym audit:",
  length(symbol_queries),
  "\n"
)


# ------------------------------------------------------------
# Expand NCBI synonym field
# ------------------------------------------------------------

gi_with_syn <- gi[
  !is.na(NCBI_Synonyms) &
  NCBI_Synonyms != ""
]

cat(
  "NCBI records with synonyms:",
  nrow(gi_with_syn),
  "\n"
)

syn_long <- gi_with_syn[
  ,
  .(
    NCBI_Synonym =
      unlist(
        strsplit(
          NCBI_Synonyms,
          "|",
          fixed = TRUE
        )
      )
  ),
  by = .(
    NCBI_GeneID,
    NCBI_Symbol,
    NCBI_GeneType
  )
]

syn_long <- unique(
  syn_long[
    !is.na(NCBI_Synonym) &
    NCBI_Synonym != "" &
    NCBI_Synonym %chin%
      symbol_queries
  ]
)

cat(
  "Relevant exact NCBI synonym matches:",
  nrow(syn_long),
  "\n\n"
)


# ------------------------------------------------------------
# Summarize one row per queried symbol
# ------------------------------------------------------------

syn_summary <- syn_long[
  ,
  .(
    SynonymGeneID_N =
      uniqueN(NCBI_GeneID),

    SynonymCurrentSymbol_N =
      uniqueN(NCBI_Symbol),

    SynonymGeneIDs =
      paste(
        sort(
          unique(NCBI_GeneID)
        ),
        collapse = "|"
      ),

    SynonymCurrentSymbols =
      paste(
        sort(
          unique(NCBI_Symbol)
        ),
        collapse = "|"
      )
  ),
  by = .(
    QuerySymbol =
      NCBI_Synonym
  )
]


symbol_audit <- data.table(
  QuerySymbol =
    symbol_queries,

  OriginalProblemClass =
    fifelse(
      symbol_queries %chin%
        simple_remaining,
      "SIMPLE_REMAINING",
      "ALIAS_AMBIGUOUS"
    )
)

idx_syn <- match(
  symbol_audit$QuerySymbol,
  syn_summary$QuerySymbol
)

symbol_audit[
  ,
  SynonymGeneID_N :=
    fifelse(
      !is.na(idx_syn),
      syn_summary$SynonymGeneID_N[
        idx_syn
      ],
      0L
    )
]

symbol_audit[
  ,
  SynonymCurrentSymbol_N :=
    fifelse(
      !is.na(idx_syn),
      syn_summary$SynonymCurrentSymbol_N[
        idx_syn
      ],
      0L
    )
]

symbol_audit[
  ,
  SynonymGeneIDs :=
    fifelse(
      !is.na(idx_syn),
      syn_summary$SynonymGeneIDs[
        idx_syn
      ],
      NA_character_
    )
]

symbol_audit[
  ,
  SynonymCurrentSymbols :=
    fifelse(
      !is.na(idx_syn),
      syn_summary$SynonymCurrentSymbols[
        idx_syn
      ],
      NA_character_
    )
]


symbol_audit[
  ,
  LocalSynonymClass := fcase(

    SynonymGeneID_N == 1 &
    SynonymCurrentSymbol_N == 1,
    "UNIQUE_CURRENT_NCBI_SYNONYM",

    SynonymGeneID_N > 1 |
    SynonymCurrentSymbol_N > 1,
    "AMBIGUOUS_CURRENT_NCBI_SYNONYM",

    default =
      "NO_CURRENT_NCBI_SYNONYM"
  )
]


cat("============================================================\n")
cat("LOCAL NCBI SYNONYM AUDIT\n")
cat("============================================================\n")

print(
  symbol_audit[
    ,
    .N,
    by = .(
      OriginalProblemClass,
      LocalSynonymClass
    )
  ][
    order(
      OriginalProblemClass,
      -N
    )
  ]
)


# ============================================================
# 5. Build unique synonym candidate rows
# ============================================================

unique_syn <- symbol_audit[
  LocalSynonymClass ==
    "UNIQUE_CURRENT_NCBI_SYNONYM"
]

if (nrow(unique_syn) > 0) {

  cand_syn <- data.table(
    SourceSymbol =
      unique_syn$QuerySymbol,

    CandidateGeneID =
      unique_syn$SynonymGeneIDs,

    PreNCBISymbol =
      unique_syn$SynonymCurrentSymbols,

    EvidenceClass =
      fifelse(
        unique_syn$OriginalProblemClass ==
          "SIMPLE_REMAINING",
        "NCBI_EXACT_SYNONYM_SIMPLE",
        "NCBI_EXACT_SYNONYM_ALIAS_AMBIGUOUS"
      )
  )

} else {

  cand_syn <- data.table(
    SourceSymbol = character(),
    CandidateGeneID = character(),
    PreNCBISymbol = character(),
    EvidenceClass = character()
  )
}


# ============================================================
# 6. Combine ALL identity candidates
# ============================================================

candidates <- rbindlist(
  list(
    cand_direct,
    cand_alias,
    cand_punct,
    loc_current[
      ,
      .(
        SourceSymbol,
        CandidateGeneID,
        PreNCBISymbol,
        EvidenceClass
      )
    ],
    cand_syn
  ),
  use.names = TRUE
)

stopifnot(
  !anyDuplicated(
    candidates$SourceSymbol
  )
)

cat("\n============================================================\n")
cat("IDENTITY CANDIDATE ACCOUNTING\n")
cat("============================================================\n")

print(
  candidates[
    ,
    .N,
    by = EvidenceClass
  ][order(-N)]
)

cat(
  "\nTotal candidate source features:",
  nrow(candidates),
  "\n"
)


# ============================================================
# 7. Reconcile every candidate GeneID against current NCBI
# ============================================================

idx_gi <- match(
  candidates$CandidateGeneID,
  gi$NCBI_GeneID
)

candidates[
  ,
  CandidateGeneIDPresentCurrentNCBI :=
    !is.na(idx_gi)
]

candidates[
  ,
  NCBI_CurrentSymbol :=
    fifelse(
      !is.na(idx_gi),
      gi$NCBI_Symbol[
        idx_gi
      ],
      NA_character_
    )
]

candidates[
  ,
  NCBI_GeneType :=
    fifelse(
      !is.na(idx_gi),
      gi$NCBI_GeneType[
        idx_gi
      ],
      NA_character_
    )
]

candidates[
  ,
  NCBI_Description :=
    fifelse(
      !is.na(idx_gi),
      gi$NCBI_Description[
        idx_gi
      ],
      NA_character_
    )
]


# ------------------------------------------------------------
# Did the provisional symbol disagree with latest NCBI?
# ------------------------------------------------------------

candidates[
  ,
  CandidateSymbolDisagreesNCBI :=
    CandidateGeneIDPresentCurrentNCBI &
    !is.na(PreNCBISymbol) &
    PreNCBISymbol !=
      NCBI_CurrentSymbol
]


# ============================================================
# 8. Global GeneID multiplicity
# ============================================================

gid_mult <- candidates[
  CandidateGeneIDPresentCurrentNCBI,
  .(
    SourceFeatureN =
      .N,

    SourceSymbols =
      paste(
        SourceSymbol,
        collapse = "|"
      ),

    CurrentPrimarySourcePresent =
      any(
        SourceSymbol ==
          NCBI_CurrentSymbol
      )
  ),
  by = .(
    NCBI_GeneID =
      CandidateGeneID,

    NCBI_CurrentSymbol
  )
]

idx_mult <- match(
  candidates$CandidateGeneID,
  gid_mult$NCBI_GeneID
)

candidates[
  ,
  CurrentGeneIDSourceN :=
    fifelse(
      !is.na(idx_mult),
      gid_mult$SourceFeatureN[
        idx_mult
      ],
      NA_integer_
    )
]

candidates[
  ,
  CurrentPrimarySourcePresent :=
    fifelse(
      !is.na(idx_mult),
      gid_mult$CurrentPrimarySourcePresent[
        idx_mult
      ],
      FALSE
    )
]


# ============================================================
# 9. Does current NCBI primary symbol already exist as a raw
#    feature in the deposited matrix?
# ============================================================

candidates[
  ,
  CurrentSymbolAlreadyRaw :=
    CandidateGeneIDPresentCurrentNCBI &
    !is.na(NCBI_CurrentSymbol) &
    NCBI_CurrentSymbol %chin%
      raw_symbols &
    NCBI_CurrentSymbol !=
      SourceSymbol
]


# ============================================================
# 10. Global current-identity classification
# ============================================================

candidates[
  ,
  ReconciliationStatus := fcase(

    !CandidateGeneIDPresentCurrentNCBI,
    "CANDIDATE_GENEID_ABSENT_CURRENT_NCBI",

    CurrentGeneIDSourceN > 1 &
    SourceSymbol ==
      NCBI_CurrentSymbol,
    "KEEP_CURRENT_PRIMARY_IN_COLLISION_GROUP",

    CurrentGeneIDSourceN > 1 &
    CurrentPrimarySourcePresent &
    SourceSymbol !=
      NCBI_CurrentSymbol,
    "EXCLUDE_NONPRIMARY_COLLISION_SOURCE",

    CurrentGeneIDSourceN > 1 &
    !CurrentPrimarySourcePresent,
    "MANY_TO_ONE_NO_CURRENT_PRIMARY_REVIEW",

    CurrentGeneIDSourceN == 1 &
    CurrentSymbolAlreadyRaw,
    "CURRENT_TARGET_RAW_PRESENT_REVIEW",

    CurrentGeneIDSourceN == 1,
    "SAFE_CURRENT_IDENTITY",

    default =
      "REVIEW"
  )
]


cat("\n============================================================\n")
cat("CURRENT NCBI RECONCILIATION STATUS\n")
cat("============================================================\n")

reconciliation_summary <- candidates[
  ,
  .N,
  by = ReconciliationStatus
][order(-N)]

print(
  reconciliation_summary
)


# ============================================================
# 11. Symbol-disagreement audit
#
# This includes the previously observed 10 safe-alias +
# 2 punctuation disagreements, and catches any additional
# org.Hs -> current NCBI primary-symbol updates.
# ============================================================

symbol_disagreements <- candidates[
  CandidateSymbolDisagreesNCBI == TRUE
]

cat("\n============================================================\n")
cat("PROVISIONAL SYMBOL vs CURRENT NCBI SYMBOL DISAGREEMENTS\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(symbol_disagreements),
  "\n"
)

if (nrow(symbol_disagreements) > 0) {

  print(
    symbol_disagreements[
      ,
      .(
        SourceSymbol,
        CandidateGeneID,
        PreNCBISymbol,
        NCBI_CurrentSymbol,
        EvidenceClass,
        CurrentSymbolAlreadyRaw,
        CurrentGeneIDSourceN,
        ReconciliationStatus
      )
    ]
  )
}


# ============================================================
# 12. Candidate collision / review rows
# ============================================================

candidate_problems <- candidates[
  ReconciliationStatus !=
    "SAFE_CURRENT_IDENTITY"
]

cat("\n============================================================\n")
cat("CANDIDATE COLLISION / REVIEW ROWS\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(candidate_problems),
  "\n"
)

if (nrow(candidate_problems) > 0) {

  print(
    candidate_problems[
      ,
      .(
        SourceSymbol,
        CandidateGeneID,
        NCBI_CurrentSymbol,
        EvidenceClass,
        CurrentGeneIDSourceN,
        CurrentPrimarySourcePresent,
        CurrentSymbolAlreadyRaw,
        ReconciliationStatus
      )
    ]
  )
}


# ============================================================
# 13. Safe identity table
# ============================================================

safe_identity <- candidates[
  ReconciliationStatus ==
    "SAFE_CURRENT_IDENTITY",
  .(
    SourceSymbol,

    FinalGeneID =
      CandidateGeneID,

    FinalSymbol =
      NCBI_CurrentSymbol,

    NCBI_GeneType,

    EvidenceClass,

    PreNCBISymbol,

    SymbolUpdatedByNCBI =
      CandidateSymbolDisagreesNCBI
  )
]

stopifnot(
  !anyDuplicated(
    safe_identity$SourceSymbol
  ),
  !anyDuplicated(
    safe_identity$FinalGeneID
  ),
  !anyDuplicated(
    safe_identity$FinalSymbol
  )
)

cat("\n============================================================\n")
cat("SAFE CURRENT IDENTITY TABLE\n")
cat("============================================================\n")

cat(
  "Safe source features:",
  nrow(safe_identity),
  "\n"
)

cat(
  "Unique FinalGeneID:",
  uniqueN(
    safe_identity$FinalGeneID
  ),
  "\n"
)

cat(
  "Unique FinalSymbol:",
  uniqueN(
    safe_identity$FinalSymbol
  ),
  "\n"
)

cat(
  "Rows receiving a newer NCBI primary symbol:",
  sum(
    safe_identity$SymbolUpdatedByNCBI
  ),
  "\n"
)


# ============================================================
# 14. Residual external-query set
# ============================================================

# ------------------------------------------------------------
# 14A. LOC GeneIDs absent from current gene_info
# ------------------------------------------------------------

loc_external <- loc_all[
  !CandidateGeneID %chin%
    gi$NCBI_GeneID,
  .(
    QueryType =
      "GENE_ID",

    Query =
      CandidateGeneID,

    GeneSymbol_raw =
      SourceSymbol,

    Reason =
      "LOC_GENEID_ABSENT_FROM_CURRENT_GENE_INFO"
  )
]

stopifnot(
  nrow(loc_external) == 11
)


# ------------------------------------------------------------
# 14B. Remaining symbol queries not uniquely resolved by
#      current NCBI synonym table
# ------------------------------------------------------------

symbol_external <- symbol_audit[
  LocalSynonymClass !=
    "UNIQUE_CURRENT_NCBI_SYNONYM",
  .(
    QueryType =
      "SYMBOL",

    Query =
      QuerySymbol,

    GeneSymbol_raw =
      QuerySymbol,

    Reason =
      fifelse(
        LocalSynonymClass ==
          "AMBIGUOUS_CURRENT_NCBI_SYNONYM",
        "AMBIGUOUS_CURRENT_NCBI_SYNONYM",
        "NO_CURRENT_NCBI_SYNONYM"
      )
  )
]


# ------------------------------------------------------------
# 14C. TEC
# ------------------------------------------------------------

tec_external <- data.table(
  QueryType = "SYMBOL",
  Query = "TEC",
  GeneSymbol_raw = "TEC",
  Reason =
    "DIRECT_SYMBOL_MAPS_TO_TWO_CURRENT_NCBI_GENEIDS"
)


# ------------------------------------------------------------
# 14D. Candidate GeneIDs absent from current NCBI, if any
# ------------------------------------------------------------

candidate_id_external <- candidates[
  CandidateGeneIDPresentCurrentNCBI == FALSE,
  .(
    QueryType =
      "GENE_ID",

    Query =
      CandidateGeneID,

    GeneSymbol_raw =
      SourceSymbol,

    Reason =
      "CANDIDATE_GENEID_ABSENT_FROM_CURRENT_GENE_INFO"
  )
]


external_pending <- unique(
  rbindlist(
    list(
      loc_external,
      symbol_external,
      tec_external,
      candidate_id_external
    ),
    use.names = TRUE,
    fill = TRUE
  )
)

cat("\n============================================================\n")
cat("RESIDUAL EXTERNAL NCBI DATASETS QUERIES\n")
cat("============================================================\n")

cat(
  "Queries still required:",
  nrow(external_pending),
  "\n"
)

print(
  external_pending[
    ,
    .N,
    by = .(
      QueryType,
      Reason
    )
  ][order(QueryType, Reason)]
)


# ============================================================
# 15. Accounting of features not in safe table
# ============================================================

known_collision_symbols <- a[
  PreliminaryIdentityStatus %chin% c(
    "ALIAS_TARGET_ALREADY_PRESENT",
    "ALIAS_MANY_TO_ONE_COLLISION"
  ),
  GeneSymbol_raw
]

composite_symbols <- u[
  StructuralClass ==
    "COMPOSITE_UNDERSCORE",
  GeneSymbol_raw
]

tec_symbol <- "TEC"

accounting <- data.table(
  Category = c(
    "Total_raw_features",
    "Safe_current_identity",
    "Known_alias_collision_source_rows",
    "Composite_underscore_features",
    "Residual_external_query_features",
    "TEC_ambiguous"
  ),

  N = c(
    13769L,
    nrow(safe_identity),
    length(
      unique(
        known_collision_symbols
      )
    ),
    length(
      unique(
        composite_symbols
      )
    ),
    uniqueN(
      external_pending[
        GeneSymbol_raw != "TEC",
        GeneSymbol_raw
      ]
    ),
    1L
  )
)

cat("\n============================================================\n")
cat("FEATURE ACCOUNTING SNAPSHOT\n")
cat("============================================================\n")

print(accounting)


# ============================================================
# 16. Save
# ============================================================

fwrite(
  symbol_audit,
  file.path(
    out_dir,
    "GSE93624_local_NCBI_synonym_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  candidates,
  file.path(
    out_dir,
    "GSE93624_all_current_identity_candidates.tsv"
  ),
  sep = "\t"
)

fwrite(
  reconciliation_summary,
  file.path(
    out_dir,
    "GSE93624_current_identity_reconciliation_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  symbol_disagreements,
  file.path(
    out_dir,
    "GSE93624_current_symbol_disagreement_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  candidate_problems,
  file.path(
    out_dir,
    "GSE93624_current_identity_collision_review.tsv"
  ),
  sep = "\t"
)

fwrite(
  safe_identity,
  file.path(
    out_dir,
    "GSE93624_safe_current_identity_candidates.tsv"
  ),
  sep = "\t"
)

fwrite(
  external_pending,
  file.path(
    out_dir,
    "GSE93624_residual_external_NCBI_queries.tsv"
  ),
  sep = "\t"
)

fwrite(
  accounting,
  file.path(
    out_dir,
    "GSE93624_stage2D_feature_accounting.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 17. Final
# ============================================================

cat("\nIMPORTANT:\n")
cat(
  "- GeneID is treated as the stable identity anchor once a unique current GeneID is established.\n",
  "- Current NCBI primary symbols supersede older org.Hs candidate symbols only at the annotation layer.\n",
  "- No expression row has been renamed yet.\n",
  "- No collision row has been merged.\n",
  "- Composite features remain untouched and unresolved.\n",
  "- Only the residual query set should proceed to external NCBI Datasets forensic.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_CURRENT_NCBI_RECONCILIATION_AND_SYNONYM_AUDIT\n")
cat("============================================================\n")
