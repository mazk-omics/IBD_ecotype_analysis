suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ============================================================
# GSE93624
# Stage 2A: gene-symbol identity forensic audit
#
# IMPORTANT:
#   - Audit only.
#   - No expression rows renamed.
#   - No rows collapsed.
#   - No rows excluded.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

expr_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE93624_authoritative_expression_245.rds"
)

out_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE93624 gene-symbol harmonization audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load authoritative matrix
# ============================================================

obj <- readRDS(expr_file)

expr <- obj$expression

raw_symbols <- rownames(expr)

stopifnot(
  nrow(expr) == 13769,
  ncol(expr) == 245,
  length(raw_symbols) == 13769,
  uniqueN(raw_symbols) == 13769,
  !anyNA(raw_symbols)
)

cat("Raw features:", length(raw_symbols), "\n")
cat("Unique raw symbols:", uniqueN(raw_symbols), "\n\n")


# ============================================================
# 2. Current SYMBOL key universe
# ============================================================

current_symbol_keys <- keys(
  org.Hs.eg.db,
  keytype = "SYMBOL"
)

direct_flag <- raw_symbols %chin% current_symbol_keys

cat(
  "Exact current SYMBOL keys:",
  sum(direct_flag),
  "\n"
)

cat(
  "Not exact current SYMBOL keys:",
  sum(!direct_flag),
  "\n\n"
)


# ============================================================
# 3. Direct current-symbol mapping
# ============================================================

direct_symbols <- raw_symbols[
  direct_flag
]

direct_map <- as.data.table(
  suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = direct_symbols,
      keytype = "SYMBOL",
      columns = c(
        "SYMBOL",
        "ENTREZID"
      )
    )
  )
)

direct_map <- unique(
  direct_map[
    !is.na(SYMBOL) &
    !is.na(ENTREZID)
  ]
)

direct_summary <- direct_map[
  ,
  .(
    NEntrez =
      uniqueN(ENTREZID),

    EntrezIDs =
      paste(
        sort(
          unique(ENTREZID)
        ),
        collapse = "|"
      )
  ),
  by = SYMBOL
]


# ============================================================
# 4. Historical ALIAS mapping for non-current symbols
# ============================================================

alias_symbols <- raw_symbols[
  !direct_flag
]

if (length(alias_symbols) > 0) {

  alias_map <- as.data.table(
    suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = alias_symbols,
        keytype = "ALIAS",
        columns = c(
          "ALIAS",
          "ENTREZID",
          "SYMBOL"
        )
      )
    )
  )

  alias_map <- unique(
    alias_map[
      !is.na(ALIAS)
    ]
  )

} else {

  alias_map <- data.table(
    ALIAS = character(),
    ENTREZID = character(),
    SYMBOL = character()
  )
}


alias_summary <- alias_map[
  ,
  .(
    NEntrez =
      uniqueN(
        ENTREZID[
          !is.na(ENTREZID)
        ]
      ),

    NCurrentSymbol =
      uniqueN(
        SYMBOL[
          !is.na(SYMBOL)
        ]
      ),

    EntrezIDs =
      paste(
        sort(
          unique(
            ENTREZID[
              !is.na(ENTREZID)
            ]
          )
        ),
        collapse = "|"
      ),

    CurrentSymbols =
      paste(
        sort(
          unique(
            SYMBOL[
              !is.na(SYMBOL)
            ]
          )
        ),
        collapse = "|"
      )
  ),
  by = ALIAS
]


# ============================================================
# 5. Build one-row-per-original-feature audit table
# ============================================================

audit <- data.table(
  FeatureOrder =
    seq_along(raw_symbols),

  GeneSymbol_raw =
    raw_symbols
)

audit[
  ,
  ExactCurrentSymbol :=
    GeneSymbol_raw %chin%
    current_symbol_keys
]


# ------------------------------------------------------------
# Direct mapping info
# ------------------------------------------------------------

idx_direct <- match(
  audit$GeneSymbol_raw,
  direct_summary$SYMBOL
)

audit[
  ,
  Direct_NEntrez :=
    fifelse(
      !is.na(idx_direct),
      direct_summary$NEntrez[idx_direct],
      NA_integer_
    )
]

audit[
  ,
  Direct_EntrezIDs :=
    fifelse(
      !is.na(idx_direct),
      direct_summary$EntrezIDs[idx_direct],
      NA_character_
    )
]


# ------------------------------------------------------------
# Alias mapping info
# ------------------------------------------------------------

idx_alias <- match(
  audit$GeneSymbol_raw,
  alias_summary$ALIAS
)

audit[
  ,
  Alias_NEntrez :=
    fifelse(
      !is.na(idx_alias),
      alias_summary$NEntrez[idx_alias],
      0L
    )
]

audit[
  ,
  Alias_NCurrentSymbol :=
    fifelse(
      !is.na(idx_alias),
      alias_summary$NCurrentSymbol[idx_alias],
      0L
    )
]

audit[
  ,
  Alias_EntrezIDs :=
    fifelse(
      !is.na(idx_alias),
      alias_summary$EntrezIDs[idx_alias],
      NA_character_
    )
]

audit[
  ,
  Alias_CurrentSymbols :=
    fifelse(
      !is.na(idx_alias),
      alias_summary$CurrentSymbols[idx_alias],
      NA_character_
    )
]


# ============================================================
# 6. Preliminary identity classification
#
# Still AUDIT ONLY.
# ============================================================

audit[
  ,
  MappingClass := fcase(

    ExactCurrentSymbol &
      Direct_NEntrez == 1,
    "DIRECT_CURRENT_UNIQUE",

    ExactCurrentSymbol &
      Direct_NEntrez > 1,
    "DIRECT_CURRENT_AMBIGUOUS",

    !ExactCurrentSymbol &
      Alias_NEntrez == 1 &
      Alias_NCurrentSymbol == 1,
    "ALIAS_SINGLE_CURRENT",

    !ExactCurrentSymbol &
      (
        Alias_NEntrez > 1 |
        Alias_NCurrentSymbol > 1
      ),
    "ALIAS_AMBIGUOUS",

    default =
      "UNRESOLVED"
  )
]


# ============================================================
# 7. Candidate current identities
# ============================================================

audit[
  ,
  CandidateEntrezID := fcase(

    MappingClass ==
      "DIRECT_CURRENT_UNIQUE",
    Direct_EntrezIDs,

    MappingClass ==
      "ALIAS_SINGLE_CURRENT",
    Alias_EntrezIDs,

    default =
      NA_character_
  )
]

audit[
  ,
  CandidateCurrentSymbol := fcase(

    MappingClass ==
      "DIRECT_CURRENT_UNIQUE",
    GeneSymbol_raw,

    MappingClass ==
      "ALIAS_SINGLE_CURRENT",
    Alias_CurrentSymbols,

    default =
      NA_character_
  )
]


# ============================================================
# 8. Collision audits
# ============================================================

# Does an alias target already exist as another original feature?
audit[
  ,
  CandidateTargetAlreadyRaw :=
    !is.na(CandidateCurrentSymbol) &
    CandidateCurrentSymbol %chin%
      raw_symbols &
    CandidateCurrentSymbol !=
      GeneSymbol_raw
]


# How many original features map to same candidate GeneID?
target_geneid_mult <- audit[
  !is.na(CandidateEntrezID),
  .(
    SourceFeatureN = .N,
    SourceSymbols =
      paste(
        GeneSymbol_raw,
        collapse = "|"
      )
  ),
  by = CandidateEntrezID
]

idx_mult_id <- match(
  audit$CandidateEntrezID,
  target_geneid_mult$CandidateEntrezID
)

audit[
  ,
  CandidateGeneIDSourceN :=
    fifelse(
      !is.na(idx_mult_id),
      target_geneid_mult$SourceFeatureN[
        idx_mult_id
      ],
      NA_integer_
    )
]


# How many original features map to same current symbol?
target_symbol_mult <- audit[
  !is.na(CandidateCurrentSymbol),
  .(
    SourceFeatureN = .N,
    SourceSymbols =
      paste(
        GeneSymbol_raw,
        collapse = "|"
      )
  ),
  by = CandidateCurrentSymbol
]

idx_mult_symbol <- match(
  audit$CandidateCurrentSymbol,
  target_symbol_mult$CandidateCurrentSymbol
)

audit[
  ,
  CandidateSymbolSourceN :=
    fifelse(
      !is.na(idx_mult_symbol),
      target_symbol_mult$SourceFeatureN[
        idx_mult_symbol
      ],
      NA_integer_
    )
]


# ============================================================
# 9. Preliminary safety class
#
# Still no canonicalization performed.
# ============================================================

audit[
  ,
  PreliminaryIdentityStatus := fcase(

    MappingClass ==
      "DIRECT_CURRENT_UNIQUE",
    "CURRENT_IDENTITY_OK",

    MappingClass ==
      "ALIAS_SINGLE_CURRENT" &
      !CandidateTargetAlreadyRaw &
      CandidateGeneIDSourceN == 1 &
      CandidateSymbolSourceN == 1,
    "SAFE_ALIAS_CANDIDATE",

    MappingClass ==
      "ALIAS_SINGLE_CURRENT" &
      CandidateTargetAlreadyRaw,
    "ALIAS_TARGET_ALREADY_PRESENT",

    MappingClass ==
      "ALIAS_SINGLE_CURRENT" &
      (
        CandidateGeneIDSourceN > 1 |
        CandidateSymbolSourceN > 1
      ),
    "ALIAS_MANY_TO_ONE_COLLISION",

    MappingClass ==
      "DIRECT_CURRENT_AMBIGUOUS",
    "DIRECT_CURRENT_AMBIGUOUS",

    MappingClass ==
      "ALIAS_AMBIGUOUS",
    "ALIAS_AMBIGUOUS",

    default =
      "UNRESOLVED"
  )
]


# ============================================================
# 10. Summary
# ============================================================

cat("============================================================\n")
cat("MAPPING CLASS SUMMARY\n")
cat("============================================================\n")

mapping_summary <- audit[
  ,
  .N,
  by = MappingClass
][order(-N)]

print(
  mapping_summary
)


cat("\n============================================================\n")
cat("PRELIMINARY IDENTITY STATUS\n")
cat("============================================================\n")

status_summary <- audit[
  ,
  .N,
  by = PreliminaryIdentityStatus
][order(-N)]

print(
  status_summary
)


cat("\n============================================================\n")
cat("ALIAS TARGET ALREADY PRESENT\n")
cat("============================================================\n")

target_present <- audit[
  PreliminaryIdentityStatus ==
    "ALIAS_TARGET_ALREADY_PRESENT"
]

cat(
  "Rows:",
  nrow(target_present),
  "\n"
)

if (nrow(target_present) > 0) {

  print(
    target_present[
      ,
      .(
        GeneSymbol_raw,
        CandidateEntrezID,
        CandidateCurrentSymbol
      )
    ]
  )
}


cat("\n============================================================\n")
cat("MANY-TO-ONE COLLISIONS\n")
cat("============================================================\n")

collision <- audit[
  !is.na(CandidateGeneIDSourceN) &
    CandidateGeneIDSourceN > 1
]

cat(
  "Rows involved:",
  nrow(collision),
  "\n"
)

if (nrow(collision) > 0) {

  print(
    collision[
      ,
      .(
        GeneSymbol_raw,
        CandidateEntrezID,
        CandidateCurrentSymbol,
        CandidateGeneIDSourceN,
        CandidateSymbolSourceN
      )
    ]
  )
}


cat("\n============================================================\n")
cat("UNRESOLVED / AMBIGUOUS\n")
cat("============================================================\n")

problem <- audit[
  PreliminaryIdentityStatus %chin% c(
    "UNRESOLVED",
    "DIRECT_CURRENT_AMBIGUOUS",
    "ALIAS_AMBIGUOUS"
  )
]

cat(
  "Rows:",
  nrow(problem),
  "\n"
)

if (nrow(problem) > 0) {

  print(
    head(
      problem[
        ,
        .(
          GeneSymbol_raw,
          MappingClass,
          Alias_NEntrez,
          Alias_NCurrentSymbol,
          Alias_EntrezIDs,
          Alias_CurrentSymbols
        )
      ],
      100
    )
  )
}


# ============================================================
# 11. Version provenance
# ============================================================

pkg_version <- as.character(
  packageVersion(
    "org.Hs.eg.db"
  )
)

annotationdbi_version <- as.character(
  packageVersion(
    "AnnotationDbi"
  )
)

cat("\n============================================================\n")
cat("ANNOTATION DATABASE\n")
cat("============================================================\n")

cat(
  "org.Hs.eg.db:",
  pkg_version,
  "\n"
)

cat(
  "AnnotationDbi:",
  annotationdbi_version,
  "\n"
)


# ============================================================
# 12. Save
# ============================================================

fwrite(
  audit,
  file.path(
    out_dir,
    "GSE93624_gene_symbol_harmonization_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  mapping_summary,
  file.path(
    out_dir,
    "GSE93624_gene_symbol_mapping_class_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  status_summary,
  file.path(
    out_dir,
    "GSE93624_gene_symbol_identity_status_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  target_present,
  file.path(
    out_dir,
    "GSE93624_alias_target_already_present.tsv"
  ),
  sep = "\t"
)

fwrite(
  collision,
  file.path(
    out_dir,
    "GSE93624_gene_symbol_many_to_one_collisions.tsv"
  ),
  sep = "\t"
)

fwrite(
  problem,
  file.path(
    out_dir,
    "GSE93624_gene_symbol_unresolved_or_ambiguous.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 13. Final
# ============================================================

stopifnot(
  nrow(audit) == 13769,
  identical(
    audit$GeneSymbol_raw,
    rownames(expr)
  )
)

cat("\nIMPORTANT:\n")
cat(
  "- This stage performs annotation audit only.\n",
  "- No expression values were changed.\n",
  "- No feature was renamed, collapsed, excluded, or added.\n",
  "- Alias mappings are preliminary candidates, not automatic canonicalization rules.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_GENE_SYMBOL_HARMONIZATION_AUDIT\n")
cat("============================================================\n")
