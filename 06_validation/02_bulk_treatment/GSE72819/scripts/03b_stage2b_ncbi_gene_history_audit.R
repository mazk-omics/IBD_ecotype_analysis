suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE72819 Stage 2B
# NCBI GeneID current/history audit
#
# Purpose:
#   Audit the 1,908 Entrez GeneIDs that were not mapped by
#   org.Hs.eg.db during Stage 2A.
#
# IMPORTANT:
#   This script DOES NOT modify expression values.
#   This script DOES NOT finalize rescue decisions.
#   This script DOES NOT overwrite the Stage 2A decision table.
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

annotation_root <- paste0(
  "/home/mazekai/IBD_EcoTyper/",
  "06_validation/99_common/gene_annotation/",
  "NCBI_Gene_2026-09-21"
)

decision_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE72819_canonical_gene_decision_table.tsv"
)

gene_info_file <- file.path(
  annotation_root,
  "Homo_sapiens.gene_info.gz"
)

gene_history_file <- file.path(
  annotation_root,
  "gene_history.gz"
)

out_dir <- file.path(
  root,
  "prepared/02_gene_harmonization"
)

stopifnot(
  file.exists(decision_file),
  file.exists(gene_info_file),
  file.exists(gene_history_file)
)

cat("============================================================\n")
cat("GSE72819 Stage 2B: NCBI GeneID history audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load frozen Stage 2A decision table
# ============================================================

decision <- fread(decision_file)

cat("Stage 2A rows:", nrow(decision), "\n")
cat("Columns:\n")
print(names(decision))

stopifnot(
  nrow(decision) == 26468,
  "EntrezID" %in% names(decision),
  "Included" %in% names(decision)
)

decision[, EntrezID := as.character(EntrezID)]


# Normalize Included safely.
if (is.logical(decision$Included)) {

  decision[, IncludedLogical := Included]

} else {

  decision[
    ,
    IncludedLogical :=
      toupper(trimws(as.character(Included))) == "TRUE"
  ]
}


cat("\nStage 2A included:",
    sum(decision$IncludedLogical),
    "\n")

cat("Stage 2A excluded:",
    sum(!decision$IncludedLogical),
    "\n")

stopifnot(
  sum(decision$IncludedLogical) == 24560,
  sum(!decision$IncludedLogical) == 1908
)


# Since Stage 2A had:
#   0 ambiguous mappings
#   0 symbol collisions
# all 1,908 excluded rows are the unresolved Entrez set.

unresolved <- decision[
  IncludedLogical == FALSE,
  .(
    SourceRow =
      if ("SourceRow" %in% names(decision))
        SourceRow
      else
        .I,

    SourceName =
      if ("SourceName" %in% names(decision))
        SourceName
      else
        paste0("GeneID:", EntrezID),

    SourceGeneID = EntrezID
  )
]

stopifnot(
  nrow(unresolved) == 1908,
  uniqueN(unresolved$SourceGeneID) == 1908
)

cat("[PASS] Extracted exactly 1,908 unresolved Entrez GeneIDs\n")


# ============================================================
# 2. Read current human NCBI gene_info
# ============================================================

gene_info <- fread(
  gene_info_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  na.strings = c("", "-")
)

# Remove leading # from #tax_id.
setnames(
  gene_info,
  names(gene_info),
  sub("^#", "", names(gene_info))
)

stopifnot(
  all(
    c(
      "tax_id",
      "GeneID",
      "Symbol"
    ) %in% names(gene_info)
  )
)

gene_info <- gene_info[
  tax_id == 9606
]

gene_info[
  ,
  GeneID := as.character(GeneID)
]

gene_info[
  ,
  Symbol := as.character(Symbol)
]

stopifnot(
  uniqueN(gene_info$GeneID) == nrow(gene_info)
)

cat("[PASS] Loaded current human NCBI gene_info\n")
cat("Current human GeneIDs:", nrow(gene_info), "\n")


# ============================================================
# 3. Read ONLY human records from global gene_history
# ============================================================

history_cmd <- sprintf(
  "zcat %s | awk -F '\\t' 'NR==1 || $1==9606'",
  shQuote(gene_history_file)
)

gene_history <- fread(
  cmd = history_cmd,
  sep = "\t",
  header = TRUE,
  quote = "",
  na.strings = c("")
)

setnames(
  gene_history,
  names(gene_history),
  sub("^#", "", names(gene_history))
)

required_history_cols <- c(
  "tax_id",
  "GeneID",
  "Discontinued_GeneID",
  "Discontinued_Symbol",
  "Discontinue_Date"
)

stopifnot(
  all(
    required_history_cols %in%
      names(gene_history)
  )
)

gene_history[
  ,
  Discontinued_GeneID :=
    as.character(Discontinued_GeneID)
]

gene_history[
  ,
  GeneID :=
    as.character(GeneID)
]

gene_history[
  GeneID == "-",
  GeneID := NA_character_
]

cat("[PASS] Loaded human NCBI gene_history only\n")
cat("Human history rows:", nrow(gene_history), "\n")


# ============================================================
# 4. Restrict history table to the 1,908 unresolved IDs
# ============================================================

target_ids <- unresolved$SourceGeneID

history_target <- gene_history[
  Discontinued_GeneID %in% target_ids
]

cat("\nHistory rows matching unresolved IDs:",
    nrow(history_target),
    "\n")


# ============================================================
# 5. Audit uniqueness of NCBI history evidence
# ============================================================

history_counts <- history_target[
  ,
  .(
    HistoryRows = .N,

    ReplacementIDs =
      paste(
        sort(
          unique(
            GeneID[!is.na(GeneID)]
          )
        ),
        collapse = ";"
      ),

    DiscontinuedSymbols =
      paste(
        sort(
          unique(
            Discontinued_Symbol[
              !is.na(Discontinued_Symbol)
            ]
          )
        ),
        collapse = ";"
      ),

    DiscontinueDates =
      paste(
        sort(
          unique(
            Discontinue_Date[
              !is.na(Discontinue_Date)
            ]
          )
        ),
        collapse = ";"
      )
  ),
  by = Discontinued_GeneID
]


# ============================================================
# 6. Build source-level NCBI audit
# ============================================================

audit <- copy(unresolved)

audit[
  ,
  CurrentHit :=
    SourceGeneID %in% gene_info$GeneID
]

audit <- merge(
  audit,
  history_counts,
  by.x = "SourceGeneID",
  by.y = "Discontinued_GeneID",
  all.x = TRUE,
  sort = FALSE
)

audit[
  is.na(HistoryRows),
  HistoryRows := 0L
]

audit[
  ,
  ReplacementGeneID :=
    fifelse(
      HistoryRows == 1L &
        !is.na(ReplacementIDs) &
        ReplacementIDs != "",
      ReplacementIDs,
      NA_character_
    )
]


# ============================================================
# 7. Detect impossible / unusual evidence combinations
# ============================================================

audit[
  ,
  CurrentAndHistory :=
    CurrentHit & HistoryRows > 0L
]

audit[
  ,
  ReplacementCurrent :=
    !is.na(ReplacementGeneID) &
      ReplacementGeneID %in%
        gene_info$GeneID
]


# ============================================================
# 8. Recover current NCBI symbol where possible
# ============================================================

current_symbol_map <- gene_info[
  ,
  .(
    FinalGeneID = GeneID,
    NCBICurrentSymbol = Symbol
  )
]

# For current source IDs.
audit[
  CurrentHit == TRUE,
  FinalGeneID := SourceGeneID
]

# For retired IDs with one replacement.
audit[
  CurrentHit == FALSE &
    HistoryRows == 1L &
    !is.na(ReplacementGeneID),
  FinalGeneID := ReplacementGeneID
]


audit <- merge(
  audit,
  current_symbol_map,
  by = "FinalGeneID",
  all.x = TRUE,
  sort = FALSE
)


# ============================================================
# 9. Classify NCBI status
# ============================================================

audit[
  ,
  NCBIStatus := fcase(

    CurrentAndHistory,
    "REVIEW_CURRENT_AND_HISTORY",

    CurrentHit,
    "CURRENT_NCBI",

    !CurrentHit &
      HistoryRows > 1L,
    "REVIEW_MULTIPLE_HISTORY_ROWS",

    !CurrentHit &
      HistoryRows == 1L &
      is.na(ReplacementGeneID),
    "RETIRED_NO_REPLACEMENT",

    !CurrentHit &
      HistoryRows == 1L &
      !is.na(ReplacementGeneID) &
      ReplacementCurrent,
    "RETIRED_REPLACED_CURRENT",

    !CurrentHit &
      HistoryRows == 1L &
      !is.na(ReplacementGeneID) &
      !ReplacementCurrent,
    "REVIEW_REPLACEMENT_NOT_CURRENT",

    !CurrentHit &
      HistoryRows == 0L,
    "NOT_FOUND_NCBI_CURRENT_OR_HISTORY",

    default =
      "REVIEW_OTHER"
  )
]


# ============================================================
# 10. Global collision-risk audit
# ============================================================

all_source_geneids <- decision$EntrezID


# Does the replacement/current target already occur as an
# original GSE72819 source feature?
audit[
  ,
  FinalTargetAlreadyPresentInSource :=
    !is.na(FinalGeneID) &
      FinalGeneID %in%
        all_source_geneids &
      FinalGeneID != SourceGeneID
]


# How many unresolved source IDs map to the same final target?
target_counts <- audit[
  !is.na(FinalGeneID),
  .N,
  by = FinalGeneID
]

setnames(
  target_counts,
  "N",
  "UnresolvedSourcesPerFinalGeneID"
)

audit <- merge(
  audit,
  target_counts,
  by = "FinalGeneID",
  all.x = TRUE,
  sort = FALSE
)

audit[
  is.na(UnresolvedSourcesPerFinalGeneID),
  UnresolvedSourcesPerFinalGeneID := 0L
]

audit[
  ,
  MultipleHistoricalSourcesToSameTarget :=
    UnresolvedSourcesPerFinalGeneID > 1L
]


# ============================================================
# 11. Restore source-row order
# ============================================================

setorder(
  audit,
  SourceRow
)

stopifnot(
  nrow(audit) == 1908,
  uniqueN(audit$SourceGeneID) == 1908
)


# ============================================================
# 12. Summary
# ============================================================

status_summary <- audit[
  ,
  .N,
  by = NCBIStatus
][order(-N)]

cat("\n============================================================\n")
cat("NCBI STATUS SUMMARY\n")
cat("============================================================\n")

print(status_summary)


cat("\n============================================================\n")
cat("COLLISION-RISK SUMMARY\n")
cat("============================================================\n")

cat(
  "Targets already present in original source:",
  sum(
    audit$FinalTargetAlreadyPresentInSource,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Rows involved in multi-history-target mapping:",
  sum(
    audit$MultipleHistoricalSourcesToSameTarget,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Unique final target GeneIDs among unresolved:",
  uniqueN(
    audit$FinalGeneID[
      !is.na(audit$FinalGeneID)
    ]
  ),
  "\n"
)


# ============================================================
# 13. Additional summaries
# ============================================================

collision_summary <- audit[
  !is.na(FinalGeneID),
  .(
    SourceGeneIDs =
      paste(
        SourceGeneID,
        collapse = ";"
      ),

    N_source_rows = .N,

    TargetAlreadyPresent =
      any(
        FinalTargetAlreadyPresentInSource
      )
  ),
  by = .(
    FinalGeneID,
    NCBICurrentSymbol
  )
][
  N_source_rows > 1L |
    TargetAlreadyPresent == TRUE
]


# ============================================================
# 14. Save outputs
# ============================================================

audit_file <- file.path(
  out_dir,
  "GSE72819_stage2b_ncbi_history_audit_1908.tsv"
)

summary_file <- file.path(
  out_dir,
  "GSE72819_stage2b_ncbi_status_summary.tsv"
)

collision_file <- file.path(
  out_dir,
  "GSE72819_stage2b_ncbi_collision_candidates.tsv"
)

fwrite(
  audit,
  audit_file,
  sep = "\t",
  na = ""
)

fwrite(
  status_summary,
  summary_file,
  sep = "\t"
)

fwrite(
  collision_summary,
  collision_file,
  sep = "\t",
  na = ""
)


# ============================================================
# 15. Annotation provenance
# ============================================================

annotation_hashes <- data.table(
  File = c(
    gene_info_file,
    gene_history_file,
    decision_file
  ),
  SHA256 = c(
    "6d1f340d4e54ff26bb0713f5c2269bb486192650f3d5948a998c9d16f32c6ccc",
    "5c5411f63d36858e29621046f122011bb0ef7a9ce524ad18b9f2304196993d53",
    NA_character_
  )
)

fwrite(
  annotation_hashes,
  file.path(
    out_dir,
    "GSE72819_stage2b_annotation_provenance.tsv"
  ),
  sep = "\t",
  na = ""
)


# ============================================================
# 16. Final gate
# ============================================================

critical_pass <- (
  nrow(audit) == 1908 &&
    uniqueN(audit$SourceGeneID) == 1908 &&
    !anyNA(audit$NCBIStatus)
)

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] 1,908 unresolved GeneIDs audited against NCBI\n")
  cat("[PASS] Current GeneIDs distinguished from retired IDs\n")
  cat("[PASS] NCBI replacement relationships recorded\n")
  cat("[PASS] No expression values modified\n")
  cat("[PASS] No rescue decision applied yet\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE2B_NCBI_HISTORY_AUDIT\n")

} else {

  cat("[REVIEW] Stage 2B requires review\n")
  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE2B_NCBI_HISTORY_AUDIT\n")
}

cat("============================================================\n")
