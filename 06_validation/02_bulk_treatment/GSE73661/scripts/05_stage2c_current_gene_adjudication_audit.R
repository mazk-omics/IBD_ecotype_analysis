suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 2C OPTIMIZED
# Category-aware current-gene adjudication audit
#
# Inputs:
#   Stage 2A:
#     ChipDb transcript-cluster -> Entrez mapping
#
#   Stage 2B:
#     GPL6244 gene_assignment parsing
#
#   NCBI:
#     Homo_sapiens.gene_info
#     gene_history
#
# Main optimization:
#   - resolve UNIQUE Entrez IDs only once
#   - vectorized gene_history traversal
#   - cached mapping of Entrez sets
#   - hash/match-based official-symbol lookup
#
# NO expression modification
# NO transcript-cluster collapsing
# NO final canonical matrix
# ============================================================


# ============================================================
# 0. Paths
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

harm_dir <- file.path(
  root,
  "prepared/02_gene_harmonization"
)

stage2a_file <- file.path(
  harm_dir,
  "GSE73661_stage2a_probe_gene_mapping_audit.tsv"
)

stage2b_parsed_file <- file.path(
  harm_dir,
  "GSE73661_stage2b_GPL6244_gene_assignment_parsed_33252.tsv"
)

stage2b_cmp_file <- file.path(
  harm_dir,
  "GSE73661_stage2b_ChipDb_vs_GPL6244_comparison.tsv"
)

ncbi_dir <- "/home/mazekai/IBD_EcoTyper/06_validation/99_common/gene_annotation/NCBI_Gene_2026-09-21"

gene_info_file <- file.path(
  ncbi_dir,
  "Homo_sapiens.gene_info.gz"
)

gene_history_file <- file.path(
  ncbi_dir,
  "gene_history.gz"
)

out_dir <- harm_dir


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
    x %in% c(
      "",
      "-",
      "---",
      "NA",
      "N/A",
      "na",
      "n/a",
      "NULL",
      "null"
    )

  x[bad] <- NA_character_

  x
}


collapse_unique <- function(x) {

  x <- clean_chr(x)

  x <- sort(
    unique(
      x[
        !is.na(x)
      ]
    )
  )

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = ";"
  )
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


extract_all_ids <- function(x) {

  x <- clean_chr(x)

  x <- x[
    !is.na(x)
  ]

  if (length(x) == 0L) {
    return(character())
  }

  z <- unlist(
    strsplit(
      x,
      ";",
      fixed = TRUE
    ),
    use.names = FALSE
  )

  z <- clean_chr(z)

  sort(
    unique(
      z[
        !is.na(z)
      ]
    )
  )
}


elapsed_now <- function(t0) {

  as.numeric(
    proc.time()[["elapsed"]] - t0
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 2C OPTIMIZED: current-gene adjudication\n")
cat("============================================================\n")

t0 <- proc.time()[["elapsed"]]


# ============================================================
# 1. Input existence
# ============================================================

for (f in c(
  stage2a_file,
  stage2b_parsed_file,
  stage2b_cmp_file,
  gene_info_file,
  gene_history_file
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
# 2. Load Stage 2A / 2B
# ============================================================

cat("\nLoading frozen Stage 2A / 2B audits...\n")

stage2a <- fread(
  stage2a_file,
  showProgress = FALSE
)

platform <- fread(
  stage2b_parsed_file,
  showProgress = FALSE
)

cmp <- fread(
  stage2b_cmp_file,
  showProgress = FALSE
)


check(
  nrow(stage2a) == 33252L,
  "Stage 2A contains 33,252 rows"
)

check(
  nrow(platform) == 33252L,
  "Stage 2B platform table contains 33,252 rows"
)

check(
  nrow(cmp) == 33252L,
  "Stage 2B comparison contains 33,252 rows"
)


check(
  identical(
    as.character(stage2a$SourceProbeID),
    as.character(platform$SourceProbeID)
  ),
  "Stage2A / platform source order identical"
)


check(
  identical(
    as.character(stage2a$SourceProbeID),
    as.character(cmp$SourceProbeID)
  ),
  "Stage2A / comparison source order identical"
)


cat(
  sprintf(
    "[INFO] Stage 2A/2B loaded in %.2f sec\n",
    elapsed_now(t0)
  )
)


# ============================================================
# 3. Load CURRENT human NCBI Gene
#
# Homo_sapiens.gene_info is already human-specific.
# Only first three fields are needed:
#   tax_id / GeneID / Symbol
# ============================================================

cat("\nLoading current NCBI Homo sapiens gene_info...\n")

gene_info <- fread(
  cmd = paste(
    "zcat",
    shQuote(gene_info_file)
  ),
  select = 1:3,
  showProgress = FALSE
)


tax_col_info <- names(gene_info)[1]


check(
  all(
    c(
      "GeneID",
      "Symbol"
    ) %in%
      names(gene_info)
  ),
  "gene_info contains GeneID and Symbol"
)


gene_info[
  ,
  GeneID :=
    clean_chr(GeneID)
]

gene_info[
  ,
  Symbol :=
    clean_chr(Symbol)
]


gene_info <- gene_info[
  !is.na(GeneID)
]


check(
  uniqueN(
    gene_info$GeneID
  ) ==
    nrow(gene_info),
  "Current human GeneIDs are unique"
)


current_gene_ids <- gene_info$GeneID


geneid_to_symbol <- setNames(
  gene_info$Symbol,
  gene_info$GeneID
)


cat(
  "[INFO] Current human GeneIDs: ",
  length(current_gene_ids),
  "\n",
  sep = ""
)


# ============================================================
# 4. Build exact OFFICIAL-symbol lookup
#
# Only exact current NCBI official symbols are accepted here.
# Synonyms/aliases are NOT used at Stage 2C.
# ============================================================

symbol_map <- gene_info[
  !is.na(Symbol),
  .(
    N_CurrentGeneIDs =
      uniqueN(GeneID),

    CurrentGeneIDs =
      collapse_unique(GeneID)
  ),
  by = Symbol
]


unique_symbol_map <- symbol_map[
  N_CurrentGeneIDs == 1L
]


symbol_to_geneid <- setNames(
  unique_symbol_map$CurrentGeneIDs,
  unique_symbol_map$Symbol
)


cat(
  "[INFO] Unique current official symbols: ",
  length(symbol_to_geneid),
  "\n",
  sep = ""
)


# ============================================================
# 5. Load human-relevant NCBI gene_history
#
# gene_history contains multiple species.
# Read only first three fields and retain tax_id 9606.
# ============================================================

cat("\nLoading NCBI gene_history...\n")

history <- fread(
  cmd = paste(
    "zcat",
    shQuote(gene_history_file)
  ),
  select = 1:3,
  showProgress = FALSE
)


tax_col_history <- names(history)[1]


check(
  all(
    c(
      "GeneID",
      "Discontinued_GeneID"
    ) %in%
      names(history)
  ),
  "gene_history contains GeneID and Discontinued_GeneID"
)


history <- history[
  suppressWarnings(
    as.integer(
      get(tax_col_history)
    )
  ) == 9606L
]


history[
  ,
  GeneID :=
    clean_chr(GeneID)
]

history[
  ,
  Discontinued_GeneID :=
    clean_chr(
      Discontinued_GeneID
    )
]


history <- history[
  !is.na(
    Discontinued_GeneID
  )
]


cat(
  "[INFO] Human gene_history rows: ",
  nrow(history),
  "\n",
  sep = ""
)


# ============================================================
# 6. Collapse gene_history to one record per discontinued ID
#
# Normally one discontinued ID has zero or one replacement.
# Any multi-replacement case is retained as ambiguous.
# ============================================================

history_map <- history[
  ,
  .(
    ReplacementIDs =
      list(
        sort(
          unique(
            GeneID[
              !is.na(GeneID)
            ]
          )
        )
      )
  ),
  by = Discontinued_GeneID
]


history_map[
  ,
  NReplacement :=
    lengths(
      ReplacementIDs
    )
]


history_map[
  ,
  ReplacementGeneID :=
    vapply(
      ReplacementIDs,
      function(z) {

        if (length(z) == 1L) {
          z
        } else {
          NA_character_
        }
      },
      character(1)
    )
]


history_map[
  ,
  ReplacementSet :=
    vapply(
      ReplacementIDs,
      collapse_unique,
      character(1)
    )
]


check(
  uniqueN(
    history_map$Discontinued_GeneID
  ) ==
    nrow(history_map),
  "gene_history lookup has unique discontinued GeneIDs"
)


cat(
  sprintf(
    "[INFO] NCBI lookup tables ready in %.2f sec\n",
    elapsed_now(t0)
  )
)


# ============================================================
# 7. Collect ALL Entrez IDs appearing anywhere in Stage2A/2B
# ============================================================

all_entrez_ids <- sort(
  unique(
    c(
      extract_all_ids(
        stage2a$ChipEntrez
      ),

      extract_all_ids(
        platform$PlatformEntrez
      )
    )
  )
)


cat(
  "\n[INFO] Unique Entrez IDs requiring current-NCBI resolution: ",
  length(all_entrez_ids),
  "\n",
  sep = ""
)


# ============================================================
# 8. FAST vectorized current-GeneID resolver
#
# Each unique Entrez ID is resolved ONCE.
# ============================================================

resolve_ids_fast <- function(
  ids,
  current_gene_ids,
  history_map,
  geneid_to_symbol,
  max_depth = 25L
) {

  ids <- clean_chr(ids)

  ids <- sort(
    unique(
      ids[
        !is.na(ids)
      ]
    )
  )


  out <- data.table(
    InputGeneID =
      ids,

    Cursor =
      ids,

    Depth =
      0L,

    Redirected =
      FALSE,

    ResolutionStatus =
      NA_character_,

    ResolvedGeneID =
      NA_character_
  )


  for (iteration in seq_len(max_depth + 1L)) {

    pending <- which(
      is.na(
        out$ResolutionStatus
      )
    )

    if (length(pending) == 0L) {
      break
    }


    cursor <- out$Cursor[
      pending
    ]


    # --------------------------------------------------------
    # Current GeneIDs
    # --------------------------------------------------------

    is_current <- cursor %chin%
      current_gene_ids


    if (any(is_current)) {

      idx <- pending[
        is_current
      ]


      out$ResolvedGeneID[
        idx
      ] <- out$Cursor[
        idx
      ]


      out$ResolutionStatus[
        idx
      ] <- ifelse(
        out$Redirected[
          idx
        ],
        "HISTORY_TO_CURRENT",
        "DIRECT_CURRENT"
      )
    }


    # --------------------------------------------------------
    # Non-current IDs
    # --------------------------------------------------------

    rem <- pending[
      !is_current
    ]


    if (length(rem) == 0L) {
      next
    }


    rem_cursor <- out$Cursor[
      rem
    ]


    hm <- match(
      rem_cursor,
      history_map$Discontinued_GeneID
    )


    no_history <- is.na(hm)


    if (any(no_history)) {

      out$ResolutionStatus[
        rem[
          no_history
        ]
      ] <- "NOT_CURRENT_NO_HISTORY"
    }


    rem2 <- rem[
      !no_history
    ]


    hm2 <- hm[
      !no_history
    ]


    if (length(rem2) == 0L) {
      next
    }


    nrep <- history_map$NReplacement[
      hm2
    ]


    # retired without replacement
    no_replacement <- nrep == 0L


    if (any(no_replacement)) {

      out$ResolutionStatus[
        rem2[
          no_replacement
        ]
      ] <- "RETIRED_NO_REPLACEMENT"
    }


    # multiple replacements
    multi_replacement <- nrep > 1L


    if (any(multi_replacement)) {

      out$ResolutionStatus[
        rem2[
          multi_replacement
        ]
      ] <- "HISTORY_MULTI_REPLACEMENT"
    }


    # unique replacement
    unique_replacement <- nrep == 1L


    if (any(unique_replacement)) {

      idx <- rem2[
        unique_replacement
      ]


      next_id <- history_map$ReplacementGeneID[
        hm2[
          unique_replacement
        ]
      ]


      cycle <- next_id ==
        out$Cursor[
          idx
        ]


      if (any(cycle)) {

        out$ResolutionStatus[
          idx[
            cycle
          ]
        ] <- "HISTORY_CYCLE"
      }


      move_idx <- idx[
        !cycle
      ]


      if (length(move_idx) > 0L) {

        out$Cursor[
          move_idx
        ] <- next_id[
          !cycle
        ]


        out$Depth[
          move_idx
        ] <- out$Depth[
          move_idx
        ] + 1L


        out$Redirected[
          move_idx
        ] <- TRUE
      }
    }
  }


  still_pending <- which(
    is.na(
      out$ResolutionStatus
    )
  )


  if (length(still_pending) > 0L) {

    out$ResolutionStatus[
      still_pending
    ] <- "HISTORY_MAX_DEPTH"
  }


  out[
    ,
    ResolvedSymbol :=
      unname(
        geneid_to_symbol[
          ResolvedGeneID
        ]
      )
  ]


  out[]
}


cat("\nResolving unique Entrez IDs...\n")

resolver <- resolve_ids_fast(
  ids = all_entrez_ids,
  current_gene_ids = current_gene_ids,
  history_map = history_map,
  geneid_to_symbol = geneid_to_symbol,
  max_depth = 25L
)


check(
  nrow(resolver) ==
    length(all_entrez_ids),
  "Every unique Entrez ID received a resolution record"
)


resolver_summary <- resolver[
  ,
  .N,
  by = ResolutionStatus
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("NCBI ENTREZ RESOLUTION SUMMARY\n")
cat("============================================================\n")

print(
  resolver_summary
)


resolved_gene_lookup <- setNames(
  resolver$ResolvedGeneID,
  resolver$InputGeneID
)


resolution_status_lookup <- setNames(
  resolver$ResolutionStatus,
  resolver$InputGeneID
)


cat(
  sprintf(
    "[INFO] Unique Entrez resolution completed in %.2f sec\n",
    elapsed_now(t0)
  )
)


# ============================================================
# 9. Cache resolver for semicolon-delimited Entrez sets
#
# Only UNIQUE raw strings are processed.
# ============================================================

resolve_set_column <- function(
  x,
  resolved_gene_lookup
) {

  x <- clean_chr(x)

  unique_strings <- sort(
    unique(
      x[
        !is.na(x)
      ]
    )
  )


  resolved_cache <- setNames(
    rep(
      NA_character_,
      length(unique_strings)
    ),
    unique_strings
  )


  n_resolved_cache <- setNames(
    integer(
      length(unique_strings)
    ),
    unique_strings
  )


  n_input_cache <- setNames(
    integer(
      length(unique_strings)
    ),
    unique_strings
  )


  n_unresolved_cache <- setNames(
    integer(
      length(unique_strings)
    ),
    unique_strings
  )


  for (s in unique_strings) {

    input_ids <- split_set(s)

    mapped <- unname(
      resolved_gene_lookup[
        input_ids
      ]
    )


    mapped <- clean_chr(
      mapped
    )


    mapped <- sort(
      unique(
        mapped[
          !is.na(mapped)
        ]
      )
    )


    resolved_cache[
      s
    ] <- collapse_unique(
      mapped
    )


    n_resolved_cache[
      s
    ] <- length(
      mapped
    )


    n_input_cache[
      s
    ] <- length(
      input_ids
    )


    n_unresolved_cache[
      s
    ] <- sum(
      is.na(
        resolved_gene_lookup[
          input_ids
        ]
      )
    )
  }


  list(
    resolved =
      unname(
        resolved_cache[
          x
        ]
      ),

    n_resolved =
      as.integer(
        unname(
          n_resolved_cache[
            x
          ]
        )
      ),

    n_input =
      as.integer(
        unname(
          n_input_cache[
            x
          ]
        )
      ),

    n_unresolved =
      as.integer(
        unname(
          n_unresolved_cache[
            x
          ]
        )
      )
  )
}


cat("\nResolving ChipDb Entrez sets...\n")

chip_resolved <- resolve_set_column(
  stage2a$ChipEntrez,
  resolved_gene_lookup
)


cat("Resolving GPL6244 Entrez sets...\n")

platform_resolved <- resolve_set_column(
  platform$PlatformEntrez,
  resolved_gene_lookup
)


# ============================================================
# 10. Exact CURRENT official-symbol lookup
# ============================================================

platform_symbols <- clean_chr(
  platform$PlatformSymbols
)


platform_symbol_current <- rep(
  NA_character_,
  length(platform_symbols)
)


unique_symbol_idx <- which(
  platform$N_PlatformSymbols == 1L &
    !is.na(
      platform_symbols
    )
)


if (length(unique_symbol_idx) > 0L) {

  platform_symbol_current[
    unique_symbol_idx
  ] <- unname(
    symbol_to_geneid[
      platform_symbols[
        unique_symbol_idx
      ]
    ]
  )
}


n_platform_symbol_current <- as.integer(
  !is.na(
    platform_symbol_current
  )
)


# ============================================================
# 11. Combined audit table
# ============================================================

audit <- data.table(
  SourceRow =
    seq_len(33252L),

  SourceProbeID =
    as.character(
      stage2a$SourceProbeID
    ),

  PlatformCategory =
    as.character(
      platform$Category
    ),

  MalformedAssignments =
    platform$MalformedAssignments,

  Stage2AStatus =
    as.character(
      stage2a$MappingStatus
    ),

  ChipEntrez =
    clean_chr(
      stage2a$ChipEntrez
    ),

  Stage2ACurrentSymbol =
    clean_chr(
      stage2a$CurrentSymbol
    ),

  ChipCurrentGeneIDs =
    chip_resolved$resolved,

  N_ChipCurrent =
    chip_resolved$n_resolved,

  N_ChipInputEntrez =
    chip_resolved$n_input,

  N_ChipUnresolvedEntrez =
    chip_resolved$n_unresolved,

  PlatformEntrez =
    clean_chr(
      platform$PlatformEntrez
    ),

  PlatformSymbols =
    platform_symbols,

  N_PlatformSymbols =
    platform$N_PlatformSymbols,

  PlatformCurrentGeneIDs =
    platform_resolved$resolved,

  N_PlatformCurrent =
    platform_resolved$n_resolved,

  N_PlatformInputEntrez =
    platform_resolved$n_input,

  N_PlatformUnresolvedEntrez =
    platform_resolved$n_unresolved,

  PlatformSymbolCurrentGeneID =
    platform_symbol_current,

  N_PlatformSymbolCurrent =
    n_platform_symbol_current,

  Stage2BEntrezConcordance =
    cmp$EntrezConcordance,

  Stage2BSymbolConcordance =
    cmp$SymbolConcordance
)


check(
  nrow(audit) == 33252L,
  "Combined adjudication table contains 33,252 rows"
)


check(
  identical(
    audit$SourceProbeID,
    as.character(
      stage2a$SourceProbeID
    )
  ),
  "Source transcript-cluster order preserved"
)


# ============================================================
# 12. Category structure
# ============================================================

cat("\n============================================================\n")
cat("CATEGORY x STAGE2A STATUS\n")
cat("============================================================\n")


category_status <- audit[
  ,
  .N,
  by = .(
    PlatformCategory,
    Stage2AStatus
  )
][
  order(
    PlatformCategory,
    -N
  )
]


print(
  category_status
)


check(
  audit[
    PlatformCategory == "main",
    .N
  ] == 28869L,
  "28,869 main-category transcript clusters confirmed"
)


check(
  audit[
    PlatformCategory != "main",
    .N
  ] == 4383L,
  "4,383 non-main platform features confirmed"
)


# ============================================================
# 13. Current GeneID set concordance
#
# Resolved sets are already sorted/canonicalized.
# ============================================================

audit[
  ,
  CurrentEntrezSetSame :=
    fcase(

      is.na(
        ChipCurrentGeneIDs
      ) &
        is.na(
          PlatformCurrentGeneIDs
        ),
      TRUE,

      !is.na(
        ChipCurrentGeneIDs
      ) &
        !is.na(
          PlatformCurrentGeneIDs
        ) &
        ChipCurrentGeneIDs ==
          PlatformCurrentGeneIDs,
      TRUE,

      default =
        FALSE
    )
]


# ============================================================
# 14. Conservative provisional adjudication
#
# NON-MAIN FEATURES ARE GATED FIRST.
#
# Malformed gene_assignment is isolated for manual audit.
# No final keep/exclude decision is frozen here.
# ============================================================

audit[
  ,
  AdjudicationClass :=
    fcase(

      # ------------------------------------------------------
      # Platform normalization / control / unmapped features
      # ------------------------------------------------------
      PlatformCategory != "main",
      "EXCLUDE_NON_MAIN_PLATFORM_CATEGORY",

      # ------------------------------------------------------
      # Parser anomaly
      # ------------------------------------------------------
      MalformedAssignments > 0L,
      "REVIEW_MALFORMED_PLATFORM_ASSIGNMENT",

      # ------------------------------------------------------
      # Strongest evidence:
      # both independent Entrez routes resolve to same current gene
      # ------------------------------------------------------
      N_ChipCurrent == 1L &
        N_PlatformCurrent == 1L &
        CurrentEntrezSetSame,
      "STRONG_CONCORDANT_CURRENT_GENE",

      # ------------------------------------------------------
      # ChipDb current ID + exact current official platform symbol
      # supporting same current gene
      # ------------------------------------------------------
      N_ChipCurrent == 1L &
        N_PlatformCurrent == 0L &
        N_PlatformSymbolCurrent == 1L &
        ChipCurrentGeneIDs ==
          PlatformSymbolCurrentGeneID,
      "CHIPDB_CURRENT_SYMBOL_SUPPORTED",

      # ------------------------------------------------------
      # ChipDb current ID only
      # ------------------------------------------------------
      N_ChipCurrent == 1L &
        N_PlatformCurrent == 0L &
        N_PlatformSymbolCurrent == 0L,
      "CHIPDB_CURRENT_ONLY",

      # ------------------------------------------------------
      # ChipDb current vs exact current platform symbol conflict
      # ------------------------------------------------------
      N_ChipCurrent == 1L &
        N_PlatformCurrent == 0L &
        N_PlatformSymbolCurrent == 1L &
        ChipCurrentGeneIDs !=
          PlatformSymbolCurrentGeneID,
      "REVIEW_CHIPDB_SYMBOL_CONFLICT",

      # ------------------------------------------------------
      # Both sides resolve uniquely but to different genes
      # ------------------------------------------------------
      N_ChipCurrent == 1L &
        N_PlatformCurrent == 1L &
        !CurrentEntrezSetSame,
      "REVIEW_UNIQUE_CURRENT_CONFLICT",

      # ------------------------------------------------------
      # ChipDb absent/unresolved; GPL6244 Entrez resolves uniquely
      # ------------------------------------------------------
      N_ChipCurrent == 0L &
        N_PlatformCurrent == 1L,
      "CANDIDATE_PLATFORM_ENTREZ_RESCUE",

      # ------------------------------------------------------
      # ChipDb ambiguous multi-gene mapping,
      # but GPL6244 resolves uniquely
      # ------------------------------------------------------
      N_ChipCurrent > 1L &
        N_PlatformCurrent == 1L,
      "CANDIDATE_PLATFORM_DISAMBIGUATION",

      # ------------------------------------------------------
      # Neither Entrez route resolves,
      # but GPL6244 has ONE exact CURRENT official symbol
      # ------------------------------------------------------
      N_ChipCurrent == 0L &
        N_PlatformCurrent == 0L &
        N_PlatformSymbolCurrent == 1L,
      "CANDIDATE_EXACT_CURRENT_SYMBOL_RESCUE",

      # ------------------------------------------------------
      # Both still multi-gene but agree exactly
      # ------------------------------------------------------
      N_ChipCurrent > 1L &
        N_PlatformCurrent > 1L &
        CurrentEntrezSetSame,
      "AMBIGUOUS_MULTI_CURRENT_EXACT_SET",

      # ------------------------------------------------------
      # At least one side remains multi-gene
      # ------------------------------------------------------
      N_ChipCurrent > 1L |
        N_PlatformCurrent > 1L,
      "AMBIGUOUS_MULTI_CURRENT",

      # ------------------------------------------------------
      # Nothing usable
      # ------------------------------------------------------
      default =
        "NO_CURRENT_GENE_IDENTITY"
    )
]


# ============================================================
# 15. Candidate current GeneID
#
# This is still PROVISIONAL.
# REVIEW / ambiguous classes remain NA.
# ============================================================

audit[
  ,
  CandidateCurrentGeneID :=
    fcase(

      AdjudicationClass ==
        "STRONG_CONCORDANT_CURRENT_GENE",
      ChipCurrentGeneIDs,

      AdjudicationClass ==
        "CHIPDB_CURRENT_SYMBOL_SUPPORTED",
      ChipCurrentGeneIDs,

      AdjudicationClass ==
        "CHIPDB_CURRENT_ONLY",
      ChipCurrentGeneIDs,

      AdjudicationClass ==
        "CANDIDATE_PLATFORM_ENTREZ_RESCUE",
      PlatformCurrentGeneIDs,

      AdjudicationClass ==
        "CANDIDATE_PLATFORM_DISAMBIGUATION",
      PlatformCurrentGeneIDs,

      AdjudicationClass ==
        "CANDIDATE_EXACT_CURRENT_SYMBOL_RESCUE",
      PlatformSymbolCurrentGeneID,

      default =
        NA_character_
    )
]


audit[
  ,
  CandidateCurrentSymbol :=
    unname(
      geneid_to_symbol[
        CandidateCurrentGeneID
      ]
    )
]


# ============================================================
# 16. Overall provisional adjudication
# ============================================================

adjudication_summary <- audit[
  ,
  .N,
  by = AdjudicationClass
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("PROVISIONAL CURRENT-GENE ADJUDICATION\n")
cat("============================================================\n")

print(
  adjudication_summary
)


# ============================================================
# 17. MAIN category only
# ============================================================

main_summary <- audit[
  PlatformCategory == "main",
  .N,
  by = AdjudicationClass
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("MAIN-CATEGORY ADJUDICATION ONLY\n")
cat("============================================================\n")

print(
  main_summary
)


# ============================================================
# 18. Malformed gene_assignment rows
# ============================================================

malformed <- audit[
  MalformedAssignments > 0L
]


cat("\n============================================================\n")
cat("MALFORMED gene_assignment ROWS\n")
cat("============================================================\n")

print(
  malformed[
    ,
    .(
      SourceProbeID,
      PlatformCategory,
      MalformedAssignments,
      Stage2AStatus,
      ChipEntrez,
      ChipCurrentGeneIDs,
      PlatformEntrez,
      PlatformSymbols,
      PlatformCurrentGeneIDs,
      AdjudicationClass
    )
  ]
)


check(
  nrow(malformed) == 5L,
  "Exactly five previously identified malformed assignments recovered"
)


# ============================================================
# 19. Unique current-gene conflicts
# ============================================================

conflicts <- audit[
  AdjudicationClass %in%
    c(
      "REVIEW_UNIQUE_CURRENT_CONFLICT",
      "REVIEW_CHIPDB_SYMBOL_CONFLICT"
    )
]


cat("\n============================================================\n")
cat("UNIQUE CURRENT-GENE CONFLICTS\n")
cat("============================================================\n")

cat(
  "Conflict rows: ",
  nrow(conflicts),
  "\n",
  sep = ""
)


if (nrow(conflicts) > 0L) {

  print(
    head(
      conflicts[
        ,
        .(
          SourceProbeID,
          PlatformCategory,
          Stage2AStatus,
          ChipEntrez,
          ChipCurrentGeneIDs,
          Stage2ACurrentSymbol,
          PlatformEntrez,
          PlatformCurrentGeneIDs,
          PlatformSymbols,
          PlatformSymbolCurrentGeneID,
          AdjudicationClass
        )
      ],
      30L
    )
  )
}


# ============================================================
# 20. Rescue / disambiguation summary
# ============================================================

rescue_classes <- c(
  "CANDIDATE_PLATFORM_ENTREZ_RESCUE",
  "CANDIDATE_PLATFORM_DISAMBIGUATION",
  "CANDIDATE_EXACT_CURRENT_SYMBOL_RESCUE"
)


rescue_summary <- audit[
  AdjudicationClass %in%
    rescue_classes,
  .N,
  by = .(
    PlatformCategory,
    AdjudicationClass
  )
][
  order(
    PlatformCategory,
    AdjudicationClass
  )
]


cat("\n============================================================\n")
cat("CANDIDATE RESCUES BY PLATFORM CATEGORY\n")
cat("============================================================\n")

print(
  rescue_summary
)


# ============================================================
# 21. Candidate current-gene collision preview
#
# Still NO expression collapse.
# ============================================================

candidate_gene_counts <- audit[
  PlatformCategory == "main" &
    !is.na(
      CandidateCurrentGeneID
    ),
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
            AdjudicationClass
          )
        ),
        collapse = ";"
      )
  ),
  by = .(
    CandidateCurrentGeneID,
    CandidateCurrentSymbol
  )
][
  order(
    -N_SourceTranscriptClusters,
    CandidateCurrentSymbol
  )
]


collision_summary <- data.table(
  Metric = c(
    "Candidate_current_genes",
    "Genes_with_one_source_transcript_cluster",
    "Genes_with_multiple_source_transcript_clusters",
    "Maximum_source_transcript_clusters_per_gene"
  ),

  N = c(
    nrow(
      candidate_gene_counts
    ),

    sum(
      candidate_gene_counts$N_SourceTranscriptClusters ==
        1L
    ),

    sum(
      candidate_gene_counts$N_SourceTranscriptClusters >
        1L
    ),

    if (
      nrow(candidate_gene_counts) > 0L
    ) {
      max(
        candidate_gene_counts$N_SourceTranscriptClusters
      )
    } else {
      0L
    }
  )
)


cat("\n============================================================\n")
cat("CURRENT-GENE COLLISION PREVIEW\n")
cat("============================================================\n")

print(
  collision_summary
)


cat("\nTranscript-cluster count distribution per candidate gene:\n")

print(
  candidate_gene_counts[
    ,
    .N,
    by = N_SourceTranscriptClusters
  ][
    order(
      N_SourceTranscriptClusters
    )
  ]
)


# ============================================================
# 22. Useful current-NCBI resolution metrics
# ============================================================

resolution_metrics <- data.table(
  Metric = c(
    "Unique_input_Entrez_IDs",
    "Direct_current",
    "History_to_current",
    "Retired_no_replacement",
    "Not_current_no_history",
    "History_multi_replacement",
    "History_cycle",
    "History_max_depth"
  ),

  N = c(
    nrow(resolver),

    resolver[
      ResolutionStatus ==
        "DIRECT_CURRENT",
      .N
    ],

    resolver[
      ResolutionStatus ==
        "HISTORY_TO_CURRENT",
      .N
    ],

    resolver[
      ResolutionStatus ==
        "RETIRED_NO_REPLACEMENT",
      .N
    ],

    resolver[
      ResolutionStatus ==
        "NOT_CURRENT_NO_HISTORY",
      .N
    ],

    resolver[
      ResolutionStatus ==
        "HISTORY_MULTI_REPLACEMENT",
      .N
    ],

    resolver[
      ResolutionStatus ==
        "HISTORY_CYCLE",
      .N
    ],

    resolver[
      ResolutionStatus ==
        "HISTORY_MAX_DEPTH",
      .N
    ]
  )
)


# ============================================================
# 23. Save artifacts
# ============================================================

fwrite(
  resolver,
  file.path(
    out_dir,
    "GSE73661_stage2c_NCBI_entrez_resolution.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  resolver_summary,
  file.path(
    out_dir,
    "GSE73661_stage2c_NCBI_entrez_resolution_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  resolution_metrics,
  file.path(
    out_dir,
    "GSE73661_stage2c_NCBI_resolution_metrics.tsv"
  ),
  sep = "\t"
)


fwrite(
  category_status,
  file.path(
    out_dir,
    "GSE73661_stage2c_category_by_stage2a_status.tsv"
  ),
  sep = "\t"
)


fwrite(
  audit,
  file.path(
    out_dir,
    "GSE73661_stage2c_current_gene_adjudication_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  adjudication_summary,
  file.path(
    out_dir,
    "GSE73661_stage2c_adjudication_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  main_summary,
  file.path(
    out_dir,
    "GSE73661_stage2c_main_category_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  malformed,
  file.path(
    out_dir,
    "GSE73661_stage2c_malformed_gene_assignment_rows.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  conflicts,
  file.path(
    out_dir,
    "GSE73661_stage2c_unique_current_gene_conflicts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  rescue_summary,
  file.path(
    out_dir,
    "GSE73661_stage2c_candidate_rescue_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  candidate_gene_counts,
  file.path(
    out_dir,
    "GSE73661_stage2c_candidate_current_gene_probe_counts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  collision_summary,
  file.path(
    out_dir,
    "GSE73661_stage2c_current_gene_collision_preview.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 24. Stage summary
# ============================================================

stage_summary <- data.table(
  Metric = c(
    "Source_features",
    "Main_features",
    "Non_main_features",
    "Malformed_platform_rows",
    "Candidate_current_gene_rows",
    "Candidate_unique_current_genes",
    "Candidate_genes_with_multiple_source_features",
    "Review_conflict_rows",
    "Candidate_rescue_rows",
    "Elapsed_seconds"
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

    nrow(malformed),

    audit[
      PlatformCategory == "main" &
        !is.na(
          CandidateCurrentGeneID
        ),
      .N
    ],

    nrow(
      candidate_gene_counts
    ),

    candidate_gene_counts[
      N_SourceTranscriptClusters >
        1L,
      .N
    ],

    nrow(conflicts),

    audit[
      AdjudicationClass %in%
        rescue_classes,
      .N
    ],

    round(
      elapsed_now(t0),
      2
    )
  )
)


cat("\n============================================================\n")
cat("STAGE 2C SUMMARY\n")
cat("============================================================\n")

print(
  stage_summary
)


fwrite(
  stage_summary,
  file.path(
    out_dir,
    "GSE73661_stage2c_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 25. Final structural gate
# ============================================================

critical_pass <- all(
  c(
    nrow(audit) ==
      33252L,

    identical(
      audit$SourceProbeID,
      as.character(
        stage2a$SourceProbeID
      )
    ),

    audit[
      PlatformCategory == "main",
      .N
    ] ==
      28869L,

    audit[
      PlatformCategory != "main",
      .N
    ] ==
      4383L,

    nrow(malformed) ==
      5L,

    nrow(resolver) ==
      length(
        all_entrez_ids
      )
  )
)


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] All 33,252 source features adjudicated\n")
  cat("[PASS] 28,869 main-category transcript clusters retained for gene-identity review\n")
  cat("[PASS] 4,383 non-main platform features separated\n")
  cat("[PASS] Unique Entrez universe resolved once against current NCBI Gene\n")
  cat("[PASS] NCBI gene_history redirects audited\n")
  cat("[PASS] Exact current official-symbol lookup completed\n")
  cat("[PASS] Candidate rescue/disambiguation classes identified\n")
  cat("[PASS] Five malformed GPL6244 assignments isolated\n")
  cat("[PASS] Current-gene collision structure previewed\n")
  cat("[PASS] No expression values modified\n")
  cat("[PASS] No transcript-cluster collapse performed\n")
  cat("[INFO] Final keep/exclude decision table remains unfrozen\n")

  cat(
    sprintf(
      "[INFO] Total elapsed time: %.2f seconds\n",
      elapsed_now(t0)
    )
  )

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE73661_STAGE2C_CURRENT_GENE_ADJUDICATION_AUDIT\n")

} else {

  cat("[REVIEW] Stage 2C structural checks require review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE73661_STAGE2C_CURRENT_GENE_ADJUDICATION_AUDIT\n")
}

cat("============================================================\n")
