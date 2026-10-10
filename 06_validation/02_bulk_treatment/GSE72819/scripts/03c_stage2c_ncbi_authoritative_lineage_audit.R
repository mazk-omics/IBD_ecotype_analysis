suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE72819 Stage 2C
# Authoritative NCBI GeneID lineage + collision audit
#
# Scope:
#   ALL 26,468 source Entrez GeneIDs
#
# Authority:
#   NCBI Homo_sapiens.gene_info
#   NCBI gene_history
#
# No expression values are touched.
# Existing Stage 2A table is NOT overwritten.
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

anno_root <- paste0(
  "/home/mazekai/IBD_EcoTyper/",
  "06_validation/99_common/gene_annotation/",
  "NCBI_Gene_2026-09-21"
)

stage2a_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE72819_canonical_gene_decision_table.tsv"
)

gene_info_file <- file.path(
  anno_root,
  "Homo_sapiens.gene_info.gz"
)

gene_history_file <- file.path(
  anno_root,
  "gene_history.gz"
)

out_dir <- file.path(
  root,
  "prepared/02_gene_harmonization"
)

stopifnot(
  file.exists(stage2a_file),
  file.exists(gene_info_file),
  file.exists(gene_history_file)
)

cat("============================================================\n")
cat("GSE72819 Stage 2C: authoritative NCBI lineage audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load original Stage 2A source table
# ============================================================

x <- fread(stage2a_file)

stopifnot(
  nrow(x) == 26468,
  all(c(
    "SourceRow",
    "SourceName",
    "EntrezID",
    "FinalSymbol",
    "Included"
  ) %in% names(x))
)

x[, SourceRow := as.integer(SourceRow)]
x[, SourceGeneID := as.character(EntrezID)]

setorder(x, SourceRow)

stopifnot(
  identical(x$SourceRow, seq_len(26468)),
  uniqueN(x$SourceGeneID) == 26468
)

cat("[PASS] Loaded 26,468 unique source Entrez GeneIDs\n")


# ============================================================
# 2. Load current human NCBI Gene
# ============================================================

gene_info <- fread(
  gene_info_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  na.strings = c("", "-")
)

setnames(
  gene_info,
  names(gene_info),
  sub("^#", "", names(gene_info))
)

stopifnot(
  all(c(
    "tax_id",
    "GeneID",
    "Symbol"
  ) %in% names(gene_info))
)

gene_info <- gene_info[tax_id == 9606]

gene_info[, GeneID := as.character(GeneID)]
gene_info[, Symbol := as.character(Symbol)]

stopifnot(
  uniqueN(gene_info$GeneID) == nrow(gene_info)
)

current_symbol <- setNames(
  gene_info$Symbol,
  gene_info$GeneID
)

cat(
  "[PASS] Current human NCBI GeneIDs:",
  nrow(gene_info),
  "\n"
)


# ============================================================
# 3. Load human NCBI gene_history
# ============================================================

history_cmd <- sprintf(
  "zcat %s | awk -F '\\t' 'NR==1 || $1==9606'",
  shQuote(gene_history_file)
)

hist <- fread(
  cmd = history_cmd,
  sep = "\t",
  header = TRUE,
  quote = "",
  na.strings = c("")
)

setnames(
  hist,
  names(hist),
  sub("^#", "", names(hist))
)

stopifnot(
  all(c(
    "GeneID",
    "Discontinued_GeneID",
    "Discontinued_Symbol",
    "Discontinue_Date"
  ) %in% names(hist))
)

hist[, GeneID := as.character(GeneID)]
hist[, Discontinued_GeneID := as.character(Discontinued_GeneID)]

hist[
  GeneID == "-" | GeneID == "",
  GeneID := NA_character_
]

cat(
  "[PASS] Human gene_history rows:",
  nrow(hist),
  "\n"
)


# ============================================================
# 4. Collapse history by discontinued GeneID
# ============================================================

hist_summary <- hist[
  ,
  .(
    HistoryRows = .N,

    ReplacementIDs = list(
      sort(unique(
        GeneID[
          !is.na(GeneID) &
          GeneID != ""
        ]
      ))
    ),

    DiscontinuedSymbols = paste(
      sort(unique(
        Discontinued_Symbol[
          !is.na(Discontinued_Symbol) &
          Discontinued_Symbol != ""
        ]
      )),
      collapse = ";"
    ),

    DiscontinueDates = paste(
      sort(unique(
        Discontinue_Date[
          !is.na(Discontinue_Date) &
          Discontinue_Date != ""
        ]
      )),
      collapse = ";"
    )
  ),
  by = Discontinued_GeneID
]

hist_index <- setNames(
  seq_len(nrow(hist_summary)),
  hist_summary$Discontinued_GeneID
)


# ============================================================
# 5. Recursive resolver
# ============================================================

resolve_geneid <- function(source_id, max_depth = 25L) {

  current <- source_id
  visited <- character()
  chain <- source_id
  depth <- 0L

  repeat {

    # Current GeneID wins immediately.
    symbol_now <- unname(current_symbol[current])

      if (
        length(symbol_now) == 1L &&
        !is.na(symbol_now) &&
        nzchar(symbol_now)
      ) {

      status <- if (depth == 0L) {
        "DIRECT_CURRENT"
      } else {
        "HISTORY_TO_CURRENT"
      }

      return(list(
        ResolveStatus = status,
        FinalGeneID = current,
        FinalSymbol = symbol_now,
        ChainDepth = depth,
        GeneIDChain = paste(chain, collapse = "->")
      ))
    }


    # Protect against cycles.
    if (current %in% visited) {

      return(list(
        ResolveStatus = "REVIEW_HISTORY_CYCLE",
        FinalGeneID = NA_character_,
        FinalSymbol = NA_character_,
        ChainDepth = depth,
        GeneIDChain = paste(chain, collapse = "->")
      ))
    }

    visited <- c(visited, current)


    # Look for current ID as a discontinued GeneID.
    idx <- match(
     current,
      hist_summary$Discontinued_GeneID
    )

    if (is.na(idx)) {
      return(list(
        ResolveStatus = "REVIEW_NOT_FOUND_CURRENT_OR_HISTORY",
        FinalGeneID = NA_character_,
        FinalSymbol = NA_character_,
        ChainDepth = depth,
        GeneIDChain = paste(chain, collapse = "->")
      ))
    }


    replacements <- hist_summary$ReplacementIDs[[idx]]

    # Explicitly discontinued without replacement.
    if (length(replacements) == 0L) {

      return(list(
        ResolveStatus = "RETIRED_NO_REPLACEMENT",
        FinalGeneID = NA_character_,
        FinalSymbol = NA_character_,
        ChainDepth = depth,
        GeneIDChain = paste(chain, collapse = "->")
      ))
    }


    # More than one replacement is not safe to resolve automatically.
    if (length(replacements) > 1L) {

      return(list(
        ResolveStatus = "REVIEW_MULTIPLE_REPLACEMENTS",
        FinalGeneID = NA_character_,
        FinalSymbol = NA_character_,
        ChainDepth = depth,
        GeneIDChain = paste(
          c(
            chain,
            paste(replacements, collapse = "|")
          ),
          collapse = "->"
        )
      ))
    }


    current <- replacements[[1]]
    chain <- c(chain, current)
    depth <- depth + 1L

    if (depth > max_depth) {

      return(list(
        ResolveStatus = "REVIEW_MAX_DEPTH_EXCEEDED",
        FinalGeneID = NA_character_,
        FinalSymbol = NA_character_,
        ChainDepth = depth,
        GeneIDChain = paste(chain, collapse = "->")
      ))
    }
  }
}


cat("\nResolving all 26,468 source GeneIDs...\n")

resolved_list <- lapply(
  x$SourceGeneID,
  resolve_geneid
)

resolved <- rbindlist(
  resolved_list,
  fill = TRUE
)

audit <- cbind(
  x[, .(
    SourceRow,
    SourceName,
    SourceGeneID,
    Stage2A_FinalSymbol = FinalSymbol,
    Stage2A_Included = Included
  )],
  resolved
)

stopifnot(
  nrow(audit) == 26468
)


# ============================================================
# 6. Resolution summary
# ============================================================

resolution_summary <- audit[
  ,
  .N,
  by = ResolveStatus
][order(-N)]

cat("\n============================================================\n")
cat("NCBI RESOLUTION SUMMARY\n")
cat("============================================================\n")

print(resolution_summary)


cat("\n============================================================\n")
cat("CHAIN DEPTH SUMMARY\n")
cat("============================================================\n")

print(
  audit[
    ,
    .N,
    by = ChainDepth
  ][order(ChainDepth)]
)


# ============================================================
# 7. Compare old Stage 2A against NCBI authority
# ============================================================

if (is.logical(audit$Stage2A_Included)) {

  audit[, Stage2A_IncludedLogical := Stage2A_Included]

} else {

  audit[
    ,
    Stage2A_IncludedLogical :=
      toupper(trimws(
        as.character(Stage2A_Included)
      )) == "TRUE"
  ]
}


audit[
  ,
  Stage2A_SymbolMatchesNCBI :=
    Stage2A_IncludedLogical &
    !is.na(FinalSymbol) &
    !is.na(Stage2A_FinalSymbol) &
    Stage2A_FinalSymbol == FinalSymbol
]


stage2a_compare <- data.table(
  Metric = c(
    "Stage2A_included",
    "Stage2A_excluded",
    "Stage2A_included_and_NCBI_direct_current",
    "Stage2A_included_not_NCBI_direct_current",
    "Stage2A_included_symbol_match_NCBI",
    "Stage2A_included_symbol_diff_NCBI"
  ),

  N = c(
    sum(audit$Stage2A_IncludedLogical),

    sum(!audit$Stage2A_IncludedLogical),

    sum(
      audit$Stage2A_IncludedLogical &
      audit$ResolveStatus == "DIRECT_CURRENT"
    ),

    sum(
      audit$Stage2A_IncludedLogical &
      audit$ResolveStatus != "DIRECT_CURRENT"
    ),

    sum(
      audit$Stage2A_SymbolMatchesNCBI,
      na.rm = TRUE
    ),

    sum(
      audit$Stage2A_IncludedLogical &
      audit$ResolveStatus == "DIRECT_CURRENT" &
      !audit$Stage2A_SymbolMatchesNCBI
    )
  )
)

cat("\n============================================================\n")
cat("STAGE 2A vs NCBI AUTHORITY\n")
cat("============================================================\n")

print(stage2a_compare)


# ============================================================
# 8. Global FinalGeneID collision structure
# ============================================================

target_stats <- audit[
  !is.na(FinalGeneID),
  .(
    SourceRowsPerTarget = .N,

    DirectCurrentSourcePresent =
      any(
        ResolveStatus == "DIRECT_CURRENT" &
        SourceGeneID == FinalGeneID
      ),

    SourceGeneIDs = paste(
      SourceGeneID,
      collapse = ";"
    )
  ),
  by = FinalGeneID
]

audit <- merge(
  audit,
  target_stats,
  by = "FinalGeneID",
  all.x = TRUE,
  sort = FALSE
)

setorder(audit, SourceRow)


# ============================================================
# 9. Preliminary GeneID-level decision
# ============================================================

audit[
  ,
  GeneIDDecision := fcase(

    ResolveStatus == "DIRECT_CURRENT",
    "KEEP_DIRECT_CURRENT",

    ResolveStatus == "HISTORY_TO_CURRENT" &
      DirectCurrentSourcePresent == TRUE,
    "EXCLUDE_HISTORY_TARGET_ALREADY_PRESENT",

    ResolveStatus == "HISTORY_TO_CURRENT" &
      DirectCurrentSourcePresent == FALSE &
      SourceRowsPerTarget > 1L,
    "EXCLUDE_HISTORY_TARGET_COLLISION",

    ResolveStatus == "HISTORY_TO_CURRENT" &
      DirectCurrentSourcePresent == FALSE &
      SourceRowsPerTarget == 1L,
    "KEEP_HISTORY_RESCUE_UNIQUE",

    ResolveStatus == "RETIRED_NO_REPLACEMENT",
    "EXCLUDE_RETIRED_NO_REPLACEMENT",

    grepl("^REVIEW_", ResolveStatus),
    ResolveStatus,

    default = "REVIEW_OTHER"
  )
]


audit[
  ,
  KeepBeforeSymbolAudit :=
    GeneIDDecision %in%
      c(
        "KEEP_DIRECT_CURRENT",
        "KEEP_HISTORY_RESCUE_UNIQUE"
      )
]


# ============================================================
# 10. Global FinalSymbol collision audit
# ============================================================

symbol_stats <- audit[
  KeepBeforeSymbolAudit == TRUE &
    !is.na(FinalSymbol),
  .(
    KeptRowsPerSymbol = .N,
    KeptGeneIDsPerSymbol = uniqueN(FinalGeneID),
    GeneIDs = paste(
      sort(unique(FinalGeneID)),
      collapse = ";"
    )
  ),
  by = FinalSymbol
]

audit <- merge(
  audit,
  symbol_stats,
  by = "FinalSymbol",
  all.x = TRUE,
  sort = FALSE
)

setorder(audit, SourceRow)


audit[
  ,
  FinalDecision := GeneIDDecision
]

audit[
  KeepBeforeSymbolAudit == TRUE &
    !is.na(KeptRowsPerSymbol) &
    (
      KeptRowsPerSymbol > 1L |
      KeptGeneIDsPerSymbol > 1L
    ),
  FinalDecision := "REVIEW_FINAL_SYMBOL_COLLISION"
]


audit[
  ,
  Included :=
    FinalDecision %in%
      c(
        "KEEP_DIRECT_CURRENT",
        "KEEP_HISTORY_RESCUE_UNIQUE"
      )
]


# ============================================================
# 11. Final summaries
# ============================================================

decision_summary <- audit[
  ,
  .N,
  by = FinalDecision
][order(-N)]

cat("\n============================================================\n")
cat("FINAL DECISION SUMMARY\n")
cat("============================================================\n")

print(decision_summary)


cat("\n============================================================\n")
cat("KEY COLLISION COUNTS\n")
cat("============================================================\n")

cat(
  "History rows whose current target already exists:",
  sum(
    audit$FinalDecision ==
      "EXCLUDE_HISTORY_TARGET_ALREADY_PRESENT"
  ),
  "\n"
)

cat(
  "History rows excluded by target collision:",
  sum(
    audit$FinalDecision ==
      "EXCLUDE_HISTORY_TARGET_COLLISION"
  ),
  "\n"
)

cat(
  "Unique history rescues retained:",
  sum(
    audit$FinalDecision ==
      "KEEP_HISTORY_RESCUE_UNIQUE"
  ),
  "\n"
)

cat(
  "Final symbol collision rows:",
  sum(
    audit$FinalDecision ==
      "REVIEW_FINAL_SYMBOL_COLLISION"
  ),
  "\n"
)

cat(
  "Total REVIEW rows:",
  sum(grepl(
    "^REVIEW_",
    audit$FinalDecision
  )),
  "\n"
)


# ============================================================
# 12. Final uniqueness checks
# ============================================================

kept <- audit[
  Included == TRUE
]

cat("\n============================================================\n")
cat("FINAL KEEP SET\n")
cat("============================================================\n")

cat("Kept rows:", nrow(kept), "\n")
cat(
  "Unique FinalGeneID:",
  uniqueN(kept$FinalGeneID),
  "\n"
)
cat(
  "Unique FinalSymbol:",
  uniqueN(kept$FinalSymbol),
  "\n"
)

stopifnot(
  uniqueN(kept$FinalGeneID) == nrow(kept),
  uniqueN(kept$FinalSymbol) == nrow(kept)
)


# ============================================================
# 13. Write outputs
# ============================================================

fwrite(
  audit,
  file.path(
    out_dir,
    "GSE72819_stage2c_ncbi_authoritative_audit.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  resolution_summary,
  file.path(
    out_dir,
    "GSE72819_stage2c_resolution_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  decision_summary,
  file.path(
    out_dir,
    "GSE72819_stage2c_decision_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  stage2a_compare,
  file.path(
    out_dir,
    "GSE72819_stage2c_stage2a_vs_ncbi.tsv"
  ),
  sep = "\t"
)

fwrite(
  audit[
    ResolveStatus == "HISTORY_TO_CURRENT"
  ],
  file.path(
    out_dir,
    "GSE72819_stage2c_history_resolved_rows.tsv"
  ),
  sep = "\t",
  na = ""
)

fwrite(
  audit[
    grepl("^REVIEW_", FinalDecision)
  ],
  file.path(
    out_dir,
    "GSE72819_stage2c_review_rows.tsv"
  ),
  sep = "\t",
  na = ""
)


# ============================================================
# 14. Final audit status
# ============================================================

review_n <- sum(
  grepl(
    "^REVIEW_",
    audit$FinalDecision
  )
)

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

cat("[PASS] All 26,468 source GeneIDs audited against NCBI\n")
cat("[PASS] Recursive GeneID history tracing completed\n")
cat("[PASS] Global current-GeneID collisions audited\n")
cat("[PASS] Global final-symbol collisions audited\n")
cat("[PASS] No expression values modified\n")

if (review_n == 0L) {

  cat("[PASS] No unresolved REVIEW rows remain\n")
  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE2C_NCBI_AUTHORITATIVE_AUDIT\n")

} else {

  cat(
    "[REVIEW] Rows requiring targeted review:",
    review_n,
    "\n"
  )

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE2C_NCBI_AUTHORITATIVE_AUDIT\n")
}

cat("============================================================\n")
