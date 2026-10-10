suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE72819 Stage 2D
# Finalize canonical gene decision table
#
# Inputs:
#   Stage 2A preliminary mapping table
#   Stage 2C authoritative NCBI audit
#
# This step:
#   - DOES NOT query any database
#   - DOES NOT alter expression
#   - DOES NOT make new rescue decisions
#   - ONLY freezes Stage 2C decisions into FINAL form
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

dir_harm <- file.path(
  root,
  "prepared/02_gene_harmonization"
)

stage2a_file <- file.path(
  dir_harm,
  "GSE72819_canonical_gene_decision_table.tsv"
)

stage2c_file <- file.path(
  dir_harm,
  "GSE72819_stage2c_ncbi_authoritative_audit.tsv"
)

final_file <- file.path(
  dir_harm,
  "GSE72819_canonical_gene_decision_table_FINAL.tsv"
)

change_file <- file.path(
  dir_harm,
  "GSE72819_stage2A_to_FINAL_changes.tsv"
)

summary_file <- file.path(
  dir_harm,
  "GSE72819_FINAL_gene_harmonization_summary.tsv"
)

stopifnot(
  file.exists(stage2a_file),
  file.exists(stage2c_file)
)

cat("============================================================\n")
cat("GSE72819 Stage 2D: finalize canonical gene decision table\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load Stage 2A provenance + Stage 2C authoritative result
# ============================================================

a <- fread(stage2a_file)
c <- fread(stage2c_file)

stopifnot(
  nrow(a) == 26468,
  nrow(c) == 26468
)

required_a <- c(
  "SourceRow",
  "SourceName",
  "EntrezID",
  "Included",
  "FinalSymbol"
)

required_c <- c(
  "SourceRow",
  "SourceName",
  "SourceGeneID",
  "ResolveStatus",
  "FinalGeneID",
  "FinalSymbol",
  "ChainDepth",
  "GeneIDChain",
  "FinalDecision",
  "Included"
)

stopifnot(
  all(required_a %in% names(a)),
  all(required_c %in% names(c))
)

a[, SourceRow := as.integer(SourceRow)]
c[, SourceRow := as.integer(SourceRow)]

setorder(a, SourceRow)
setorder(c, SourceRow)

stopifnot(
  identical(a$SourceRow, seq_len(26468)),
  identical(c$SourceRow, seq_len(26468)),
  identical(a$SourceName, c$SourceName),
  identical(as.character(a$EntrezID), as.character(c$SourceGeneID))
)

cat("[PASS] Stage 2A and Stage 2C align exactly by 26,468 source rows\n")


# ============================================================
# 2. Normalize logical fields
# ============================================================

to_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  toupper(trimws(as.character(x))) == "TRUE"
}

a[, Stage2A_IncludedLogical := to_logical(Included)]
c[, FinalIncludedLogical := to_logical(Included)]

stopifnot(
  sum(a$Stage2A_IncludedLogical) == 24560,
  sum(c$FinalIncludedLogical) == 24656
)


# ============================================================
# 3. Construct FINAL decision table
# ============================================================

final <- data.table(
  SourceRow = c$SourceRow,
  SourceName = c$SourceName,
  SourceGeneID = as.character(c$SourceGeneID),

  Width = if ("Width" %in% names(a)) a$Width else NA,

  Stage2A_CandidateSymbol =
    if ("CandidateSymbol" %in% names(a))
      as.character(a$CandidateSymbol)
    else
      NA_character_,

  Stage2A_FinalSymbol =
    as.character(a$FinalSymbol),

  Stage2A_Included =
    a$Stage2A_IncludedLogical,

  ResolveStatus =
    as.character(c$ResolveStatus),

  ChainDepth =
    as.integer(c$ChainDepth),

  GeneIDChain =
    as.character(c$GeneIDChain),

  FinalGeneID =
    as.character(c$FinalGeneID),

  FinalSymbol =
    as.character(c$FinalSymbol),

  FinalDecision =
    as.character(c$FinalDecision),

  Included =
    c$FinalIncludedLogical
)

stopifnot(
  nrow(final) == 26468,
  uniqueN(final$SourceRow) == 26468,
  uniqueN(final$SourceGeneID) == 26468
)


# ============================================================
# 4. Revalidate exact Stage 2C resolution counts
# ============================================================

resolution_summary <- final[
  ,
  .N,
  by = ResolveStatus
][order(-N)]

decision_summary <- final[
  ,
  .N,
  by = FinalDecision
][order(-N)]

cat("============================================================\n")
cat("FINAL RESOLUTION SUMMARY\n")
cat("============================================================\n")

print(resolution_summary)

cat("\n============================================================\n")
cat("FINAL DECISION SUMMARY\n")
cat("============================================================\n")

print(decision_summary)


# Frozen Stage 2C expectations
stopifnot(
  final[ResolveStatus == "DIRECT_CURRENT", .N] == 24555,
  final[ResolveStatus == "HISTORY_TO_CURRENT", .N] == 385,
  final[ResolveStatus == "RETIRED_NO_REPLACEMENT", .N] == 1528,

  final[FinalDecision == "KEEP_DIRECT_CURRENT", .N] == 24555,
  final[FinalDecision == "KEEP_HISTORY_RESCUE_UNIQUE", .N] == 101,
  final[FinalDecision == "EXCLUDE_HISTORY_TARGET_ALREADY_PRESENT", .N] == 278,
  final[FinalDecision == "EXCLUDE_HISTORY_TARGET_COLLISION", .N] == 6,
  final[FinalDecision == "EXCLUDE_RETIRED_NO_REPLACEMENT", .N] == 1528
)


# ============================================================
# 5. Final KEEP-set integrity
# ============================================================

keep <- final[
  Included == TRUE
]

exclude <- final[
  Included == FALSE
]

stopifnot(
  nrow(keep) == 24656,
  nrow(exclude) == 1812,

  !anyNA(keep$FinalGeneID),
  !anyNA(keep$FinalSymbol),

  uniqueN(keep$FinalGeneID) == 24656,
  uniqueN(keep$FinalSymbol) == 24656,

  !any(grepl("^REVIEW_", final$FinalDecision)),
  !any(grepl("^REVIEW_", final$ResolveStatus))
)

cat("\n============================================================\n")
cat("FINAL CANONICAL KEEP SET\n")
cat("============================================================\n")

cat("Source rows:             ", nrow(final), "\n")
cat("Included rows:           ", nrow(keep), "\n")
cat("Excluded rows:           ", nrow(exclude), "\n")
cat("Unique FinalGeneID:      ", uniqueN(keep$FinalGeneID), "\n")
cat("Unique FinalSymbol:      ", uniqueN(keep$FinalSymbol), "\n")
cat("Unresolved REVIEW rows:  ", sum(grepl("^REVIEW_", final$FinalDecision)), "\n")


# ============================================================
# 6. Audit Stage 2A → FINAL changes
# ============================================================

final[
  ,
  InclusionChange := fcase(
    Stage2A_Included & Included,
    "KEEP_TO_KEEP",

    Stage2A_Included & !Included,
    "KEEP_TO_EXCLUDE",

    !Stage2A_Included & Included,
    "EXCLUDE_TO_KEEP",

    !Stage2A_Included & !Included,
    "EXCLUDE_TO_EXCLUDE",

    default = "REVIEW"
  )
]

final[
  ,
  SymbolChangedFromStage2A :=
    Included == TRUE &
    Stage2A_Included == TRUE &
    !is.na(Stage2A_FinalSymbol) &
    !is.na(FinalSymbol) &
    Stage2A_FinalSymbol != FinalSymbol
]


change_summary <- final[
  ,
  .N,
  by = InclusionChange
][order(-N)]

cat("\n============================================================\n")
cat("STAGE 2A → FINAL INCLUSION CHANGES\n")
cat("============================================================\n")

print(change_summary)

cat(
  "\nIncluded rows with Stage2A vs NCBI symbol change:",
  sum(final$SymbolChangedFromStage2A),
  "\n"
)


# Keep only rows whose inclusion or symbol changed.
changes <- final[
  InclusionChange != "KEEP_TO_KEEP" |
    SymbolChangedFromStage2A == TRUE
]


# ============================================================
# 7. Summary table
# ============================================================

summary <- data.table(
  Metric = c(
    "Source_features",
    "Direct_current",
    "History_to_current",
    "Retired_no_replacement",
    "History_target_already_present",
    "History_target_collision",
    "History_unique_rescue",
    "Final_included",
    "Final_excluded",
    "Unique_FinalGeneID",
    "Unique_FinalSymbol",
    "Review_rows",
    "Stage2A_included",
    "Stage2A_excluded",
    "Stage2A_KEEP_to_FINAL_EXCLUDE",
    "Stage2A_EXCLUDE_to_FINAL_KEEP",
    "Included_symbol_changes_vs_Stage2A"
  ),

  N = c(
    nrow(final),
    final[ResolveStatus == "DIRECT_CURRENT", .N],
    final[ResolveStatus == "HISTORY_TO_CURRENT", .N],
    final[ResolveStatus == "RETIRED_NO_REPLACEMENT", .N],
    final[FinalDecision == "EXCLUDE_HISTORY_TARGET_ALREADY_PRESENT", .N],
    final[FinalDecision == "EXCLUDE_HISTORY_TARGET_COLLISION", .N],
    final[FinalDecision == "KEEP_HISTORY_RESCUE_UNIQUE", .N],
    nrow(keep),
    nrow(exclude),
    uniqueN(keep$FinalGeneID),
    uniqueN(keep$FinalSymbol),
    sum(grepl("^REVIEW_", final$FinalDecision)),
    sum(final$Stage2A_Included),
    sum(!final$Stage2A_Included),
    final[InclusionChange == "KEEP_TO_EXCLUDE", .N],
    final[InclusionChange == "EXCLUDE_TO_KEEP", .N],
    sum(final$SymbolChangedFromStage2A)
  )
)

cat("\n============================================================\n")
cat("FINAL HARMONIZATION SUMMARY\n")
cat("============================================================\n")

print(summary)


# ============================================================
# 8. Save FINAL outputs
# ============================================================

fwrite(
  final,
  final_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  changes,
  change_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  summary,
  summary_file,
  sep = "\t",
  quote = FALSE
)

fwrite(
  decision_summary,
  file.path(
    dir_harm,
    "GSE72819_FINAL_decision_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  resolution_summary,
  file.path(
    dir_harm,
    "GSE72819_FINAL_resolution_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 9. Round-trip re-read
# ============================================================

reload <- fread(final_file)

stopifnot(
  nrow(reload) == 26468,
  sum(to_logical(reload$Included)) == 24656,
  uniqueN(
    reload[to_logical(Included), FinalGeneID]
  ) == 24656,
  uniqueN(
    reload[to_logical(Included), FinalSymbol]
  ) == 24656
)

cat("\n[PASS] FINAL decision table round-trip validated\n")


# ============================================================
# 10. Hash FINAL decision table
# ============================================================

md5 <- tools::md5sum(final_file)

writeLines(
  paste(
    unname(md5),
    basename(final_file)
  ),
  file.path(
    dir_harm,
    "GSE72819_canonical_gene_decision_table_FINAL.md5"
  )
)


# ============================================================
# 11. Final gate
# ============================================================

critical_pass <- (
  nrow(final) == 26468 &&
  nrow(keep) == 24656 &&
  nrow(exclude) == 1812 &&
  uniqueN(keep$FinalGeneID) == 24656 &&
  uniqueN(keep$FinalSymbol) == 24656 &&
  sum(grepl("^REVIEW_", final$FinalDecision)) == 0
)

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] 26,468 source GeneIDs accounted for\n")
  cat("[PASS] 24,656 canonical genes retained\n")
  cat("[PASS] 1,812 source rows excluded\n")
  cat("[PASS] FinalGeneID values are unique\n")
  cat("[PASS] FinalSymbol values are unique\n")
  cat("[PASS] No unresolved REVIEW rows remain\n")
  cat("[PASS] Stage 2C decisions frozen without modification\n")
  cat("[PASS] No expression values accessed or modified\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE2D_FINAL_GENE_DECISION_TABLE\n")

} else {

  cat("[REVIEW] Final gene decision table requires review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE2D_FINAL_GENE_DECISION_TABLE\n")
}

cat("============================================================\n")
