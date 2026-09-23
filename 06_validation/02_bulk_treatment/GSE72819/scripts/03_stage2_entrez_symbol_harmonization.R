suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ============================================================
# GSE72819 Stage 2
# Lightweight Entrez GeneID -> current SYMBOL harmonization
#
# Input:
#   26,468 unique source Entrez GeneIDs
#
# Policy:
#   - stable Entrez GeneID is authoritative source identifier
#   - use installed org.Hs.eg.db only
#   - unique Entrez -> unique SYMBOL: keep candidate
#   - no mapping: exclude
#   - one Entrez -> multiple SYMBOLs: exclude
#   - multiple Entrez -> same SYMBOL: exclude all collision rows
#   - NO historical redirect forensic
#   - NO alias archaeology
#   - NO expression aggregation
#   - NO normalization
#
# Targeted rescue is deferred until final EcoTyper signatures
# exist and only if an excluded gene becomes relevant.
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

source_manifest_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE72819_source_gene_manifest.tsv"
)

rpkm_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE72819_authoritative_RPKM_26468x73.rds"
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


cat("============================================================\n")
cat("GSE72819 Stage 2 Entrez -> SYMBOL harmonization\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load frozen Stage-1 resources
# ============================================================

source_manifest <- fread(
  source_manifest_file,
  colClasses = "character"
)

rpkm <- readRDS(
  rpkm_file
)


stopifnot(
  nrow(source_manifest) == 26468,
  nrow(rpkm) == 26468,
  ncol(rpkm) == 73
)


required_cols <- c(
  "SourceName",
  "EntrezID",
  "ValidEntrez",
  "Width"
)

stopifnot(
  all(
    required_cols %in%
      names(source_manifest)
  )
)


source_manifest[
  ,
  SourceName :=
    as.character(SourceName)
]

source_manifest[
  ,
  EntrezID :=
    trimws(
      as.character(EntrezID)
    )
]


stopifnot(
  !anyNA(source_manifest$EntrezID),

  all(
    grepl(
      "^[0-9]+$",
      source_manifest$EntrezID
    )
  ),

  uniqueN(
    source_manifest$EntrezID
  ) == 26468,

  identical(
    source_manifest$SourceName,
    rownames(rpkm)
  ),

  identical(
    paste0(
      "GeneID:",
      source_manifest$EntrezID
    ),
    source_manifest$SourceName
  )
)


cat("[PASS] Stage-1 source manifest: 26,468 unique Entrez IDs\n")
cat("[PASS] Source manifest exactly aligned to authoritative RPKM\n")


# ============================================================
# 2. Record annotation resource version
# ============================================================

annotation_info <- data.table(
  Package = c(
    "AnnotationDbi",
    "org.Hs.eg.db"
  ),

  Version = c(
    as.character(
      packageVersion(
        "AnnotationDbi"
      )
    ),

    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    )
  )
)


cat("\n============================================================\n")
cat("ANNOTATION RESOURCE\n")
cat("============================================================\n")

print(
  annotation_info
)


db_metadata <- as.data.table(
  AnnotationDbi::metadata(
    org.Hs.eg.db
  )
)

fwrite(
  annotation_info,
  file.path(
    out_dir,
    "GSE72819_annotation_package_versions.tsv"
  ),
  sep = "\t"
)

fwrite(
  db_metadata,
  file.path(
    out_dir,
    "GSE72819_org_Hs_eg_db_metadata.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 3. Entrez -> SYMBOL query
#
# Do NOT use mapIds(..., multiVals="first"), because that would
# silently hide ambiguous mappings.
# ============================================================

cat("\nQuerying org.Hs.eg.db...\n")


map_raw <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,

    keys =
      source_manifest$EntrezID,

    keytype =
      "ENTREZID",

    columns = c(
      "SYMBOL",
      "GENENAME"
    )
  )
)


map_raw <- as.data.table(
  map_raw
)


map_raw[
  ,
  ENTREZID :=
    as.character(ENTREZID)
]

map_raw[
  ,
  SYMBOL :=
    as.character(SYMBOL)
]

map_raw[
  ,
  GENENAME :=
    as.character(GENENAME)
]


# Remove completely empty SYMBOL mappings from candidate table.
map_nonempty <- map_raw[
  !is.na(SYMBOL) &
  trimws(SYMBOL) != ""
]


# Remove exact duplicate annotation rows.
map_nonempty <- unique(
  map_nonempty,
  by = c(
    "ENTREZID",
    "SYMBOL",
    "GENENAME"
  )
)


# ============================================================
# 4. Summarize mapping multiplicity per Entrez ID
# ============================================================

map_summary <- map_nonempty[
  ,
  .(
    SymbolN =
      uniqueN(SYMBOL),

    SymbolSet =
      paste(
        sort(
          unique(SYMBOL)
        ),
        collapse = " | "
      ),

    CandidateSymbol =
      if (
        uniqueN(SYMBOL) == 1
      ) {
        unique(SYMBOL)[1]
      } else {
        NA_character_
      },

    GeneName =
      if (
        uniqueN(SYMBOL) == 1
      ) {

        gn <- unique(
          GENENAME[
            !is.na(GENENAME) &
            GENENAME != ""
          ]
        )

        if (length(gn) >= 1) {
          gn[1]
        } else {
          NA_character_
        }

      } else {

        NA_character_
      }
  ),

  by = ENTREZID
]


# ============================================================
# 5. Construct one-row-per-source-gene decision table
# ============================================================

idx <- match(
  source_manifest$EntrezID,
  map_summary$ENTREZID
)


decision <- data.table(
  SourceRow =
    seq_len(
      nrow(source_manifest)
    ),

  SourceName =
    source_manifest$SourceName,

  EntrezID =
    source_manifest$EntrezID,

  Width =
    suppressWarnings(
      as.numeric(
        source_manifest$Width
      )
    ),

  SymbolN =
    map_summary$SymbolN[
      idx
    ],

  SymbolSet =
    map_summary$SymbolSet[
      idx
    ],

  CandidateSymbol =
    map_summary$CandidateSymbol[
      idx
    ],

  GeneName =
    map_summary$GeneName[
      idx
    ]
)


decision[
  is.na(SymbolN),
  SymbolN := 0L
]


# ============================================================
# 6. Initial Entrez mapping decision
# ============================================================

decision[
  SymbolN == 0,
  InitialDecision :=
    "EXCLUDE_UNMAPPED_ENTREZ"
]

decision[
  SymbolN == 1,
  InitialDecision :=
    "CANDIDATE_UNIQUE_ENTREZ_TO_SYMBOL"
]

decision[
  SymbolN > 1,
  InitialDecision :=
    "EXCLUDE_AMBIGUOUS_ENTREZ_TO_SYMBOL"
]


stopifnot(
  !anyNA(
    decision$InitialDecision
  )
)


cat("\n============================================================\n")
cat("ENTREZ -> SYMBOL INITIAL MAPPING\n")
cat("============================================================\n")


initial_summary <- decision[
  ,
  .N,
  by = InitialDecision
][order(-N)]

print(
  initial_summary
)


# ============================================================
# 7. SYMBOL collision audit
#
# Multiple independent Entrez IDs mapping to the same symbol
# are NOT aggregated.
#
# Conservative rule:
#   exclude all members of such a symbol collision.
# ============================================================

candidate_symbol_counts <- decision[
  InitialDecision ==
    "CANDIDATE_UNIQUE_ENTREZ_TO_SYMBOL",

  .(
    SymbolMultiplicity =
      .N
  ),

  by = CandidateSymbol
]


sym_idx <- match(
  decision$CandidateSymbol,
  candidate_symbol_counts$CandidateSymbol
)


decision[
  ,
  SymbolMultiplicity :=
    candidate_symbol_counts$SymbolMultiplicity[
      sym_idx
    ]
]


decision[
  InitialDecision !=
    "CANDIDATE_UNIQUE_ENTREZ_TO_SYMBOL",

  SymbolMultiplicity :=
    NA_integer_
]


collision_symbols <- candidate_symbol_counts[
  SymbolMultiplicity > 1
]


collision_rows <- decision[
  InitialDecision ==
    "CANDIDATE_UNIQUE_ENTREZ_TO_SYMBOL" &
  SymbolMultiplicity > 1
]


cat("\n============================================================\n")
cat("SYMBOL COLLISION AUDIT\n")
cat("============================================================\n")

cat(
  "Candidate unique Entrez->SYMBOL rows:",
  sum(
    decision$InitialDecision ==
      "CANDIDATE_UNIQUE_ENTREZ_TO_SYMBOL"
  ),
  "\n"
)

cat(
  "Collided symbols:",
  nrow(collision_symbols),
  "\n"
)

cat(
  "Source rows involved in symbol collisions:",
  nrow(collision_rows),
  "\n"
)


if (nrow(collision_rows) > 0) {

  cat(
    "\nExamples of symbol collisions:\n"
  )

  print(
    head(
      collision_rows[
        order(
          CandidateSymbol,
          EntrezID
        ),
        .(
          EntrezID,
          CandidateSymbol,
          GeneName,
          SymbolMultiplicity
        )
      ],
      30
    )
  )
}


# ============================================================
# 8. Final decision
# ============================================================

decision[
  InitialDecision ==
    "EXCLUDE_UNMAPPED_ENTREZ",

  FinalDecision :=
    "EXCLUDE_UNMAPPED_ENTREZ"
]


decision[
  InitialDecision ==
    "EXCLUDE_AMBIGUOUS_ENTREZ_TO_SYMBOL",

  FinalDecision :=
    "EXCLUDE_AMBIGUOUS_ENTREZ_TO_SYMBOL"
]


decision[
  InitialDecision ==
    "CANDIDATE_UNIQUE_ENTREZ_TO_SYMBOL" &
  SymbolMultiplicity == 1,

  FinalDecision :=
    "KEEP_UNIQUE_ENTREZ_SYMBOL"
]


decision[
  InitialDecision ==
    "CANDIDATE_UNIQUE_ENTREZ_TO_SYMBOL" &
  SymbolMultiplicity > 1,

  FinalDecision :=
    "EXCLUDE_SYMBOL_COLLISION"
]


stopifnot(
  !anyNA(
    decision$FinalDecision
  )
)


decision[
  FinalDecision ==
    "KEEP_UNIQUE_ENTREZ_SYMBOL",

  `:=`(
    Included = TRUE,
    FinalGeneID = EntrezID,
    FinalSymbol = CandidateSymbol
  )
]


decision[
  FinalDecision !=
    "KEEP_UNIQUE_ENTREZ_SYMBOL",

  `:=`(
    Included = FALSE,
    FinalGeneID = NA_character_,
    FinalSymbol = NA_character_
  )
]


# ============================================================
# 9. Final mapping integrity
# ============================================================

included <- decision[
  Included == TRUE
]


stopifnot(
  nrow(included) > 0,

  uniqueN(
    included$FinalGeneID
  ) ==
    nrow(included),

  uniqueN(
    included$FinalSymbol
  ) ==
    nrow(included),

  !anyNA(
    included$FinalGeneID
  ),

  !anyNA(
    included$FinalSymbol
  ),

  all(
    grepl(
      "^[0-9]+$",
      included$FinalGeneID
    )
  )
)


# ============================================================
# 10. Final summaries
# ============================================================

final_summary <- decision[
  ,
  .N,
  by = FinalDecision
][order(-N)]


cat("\n============================================================\n")
cat("FINAL GENE HARMONIZATION DECISIONS\n")
cat("============================================================\n")

print(
  final_summary
)


cat(
  "\nSource genes:",
  nrow(decision),
  "\n"
)

cat(
  "Included canonical genes:",
  nrow(included),
  "\n"
)

cat(
  "Excluded genes:",
  sum(!decision$Included),
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
# 11. Unmapped examples
# ============================================================

unmapped <- decision[
  FinalDecision ==
    "EXCLUDE_UNMAPPED_ENTREZ"
]


if (nrow(unmapped) > 0) {

  cat("\n============================================================\n")
  cat("UNMAPPED ENTREZ EXAMPLES\n")
  cat("============================================================\n")

  print(
    head(
      unmapped[
        ,
        .(
          SourceName,
          EntrezID
        )
      ],
      30
    )
  )
}


# ============================================================
# 12. Save authoritative harmonization products
# ============================================================

fwrite(
  decision,
  file.path(
    out_dir,
    "GSE72819_canonical_gene_decision_table.tsv"
  ),
  sep = "\t",
  na = ""
)


fwrite(
  initial_summary,
  file.path(
    out_dir,
    "GSE72819_initial_mapping_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  final_summary,
  file.path(
    out_dir,
    "GSE72819_final_gene_harmonization_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  collision_rows,
  file.path(
    out_dir,
    "GSE72819_symbol_collision_rows.tsv"
  ),
  sep = "\t",
  na = ""
)


fwrite(
  unmapped,
  file.path(
    out_dir,
    "GSE72819_unmapped_entrez.tsv"
  ),
  sep = "\t",
  na = ""
)


# ============================================================
# 13. Mapping-method audit
# ============================================================

method_audit <- data.table(
  Field = c(
    "SourceIdentifier",
    "SourceGeneN",
    "AnnotationResource",
    "AnnotationResourceVersion",
    "UnmappedPolicy",
    "AmbiguousEntrezPolicy",
    "SymbolCollisionPolicy",
    "HistoricalRedirectForensic",
    "AliasRescue",
    "ExpressionAggregation",
    "TargetedFutureRescue"
  ),

  Value = c(
    "Entrez GeneID parsed from GeneID:<ID>",
    "26468",
    "org.Hs.eg.db",
    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    ),
    "Exclude",
    "Exclude",
    "Exclude all collision members",
    "Not performed",
    "Not performed",
    "Not performed",
    "Only if excluded genes later occur in final EcoTyper recovery signatures"
  )
)


fwrite(
  method_audit,
  file.path(
    out_dir,
    "GSE72819_gene_harmonization_method_audit.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 14. Final status
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")


critical_pass <- (
  nrow(decision) == 26468 &&

  nrow(included) > 0 &&

  uniqueN(
    included$FinalGeneID
  ) ==
    nrow(included) &&

  uniqueN(
    included$FinalSymbol
  ) ==
    nrow(included) &&

  !anyNA(
    decision$FinalDecision
  )
)


if (critical_pass) {

  cat("[PASS] 26,468 source Entrez IDs evaluated\n")
  cat("[PASS] Mapping performed using org.Hs.eg.db\n")
  cat("[PASS] Unmapped Entrez IDs excluded conservatively\n")
  cat("[PASS] Ambiguous Entrez mappings excluded conservatively\n")
  cat("[PASS] Symbol collisions excluded without aggregation\n")
  cat("[PASS] Included FinalGeneID values are unique\n")
  cat("[PASS] Included FinalSymbol values are unique\n")
  cat("[PASS] No historical/alias forensic performed\n")
  cat("[PASS] No expression values modified\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE2_GENE_HARMONIZATION\n")

} else {

  cat("[REVIEW] Gene harmonization requires review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE2_GENE_HARMONIZATION\n")
}

cat("============================================================\n")
