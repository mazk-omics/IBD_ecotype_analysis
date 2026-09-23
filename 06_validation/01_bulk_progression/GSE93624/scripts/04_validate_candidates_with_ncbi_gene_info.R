suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 2C
#
# Cross-validate candidate gene identities against the current
# NCBI Homo_sapiens.gene_info reference already downloaded
# for the validation project.
#
# This stage:
#   - validates safe org.Hs alias candidates
#   - validates punctuation-only candidates
#   - checks LOC-style symbols by numeric GeneID
#   - checks unresolved simple symbols against current NCBI symbols
#   - investigates TEC current-symbol ambiguity
#
# This stage DOES NOT:
#   - modify expression
#   - rename expression rows
#   - collapse expression rows
#   - build the final canonical matrix
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

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
cat("GSE93624 NCBI gene_info identity validation\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load Stage 2A / 2B resources
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

stopifnot(
  nrow(a) == 13769,
  nrow(u) == 776
)

cat("[PASS] Stage 2A/2B audit resources loaded\n")


# ============================================================
# 2. Load current NCBI Gene reference
# ============================================================

if (!file.exists(gene_info_file)) {
  stop(
    "NCBI gene_info file not found: ",
    gene_info_file
  )
}

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
  old = c(
    "GeneID",
    "Symbol",
    "Synonyms",
    "description",
    "type_of_gene"
  ),
  new = c(
    "NCBI_GeneID",
    "NCBI_Symbol",
    "NCBI_Synonyms",
    "NCBI_Description",
    "NCBI_GeneType"
  )
)

stopifnot(
  !anyDuplicated(gi$NCBI_GeneID)
)

cat(
  "Current NCBI human Gene records:",
  nrow(gi),
  "\n\n"
)


# ============================================================
# 3. Helper: validate an EntrezID + symbol candidate
# ============================================================

validate_id_symbol <- function(
  source_symbol,
  candidate_id,
  candidate_symbol
) {

  hit <- gi[
    NCBI_GeneID == candidate_id
  ]

  if (nrow(hit) == 0) {

    return(
      data.table(
        SourceSymbol = source_symbol,
        CandidateGeneID = candidate_id,
        CandidateSymbol = candidate_symbol,
        NCBIPresent = FALSE,
        NCBICurrentSymbol = NA_character_,
        NCBIGeneType = NA_character_,
        IDSymbolConcordant = FALSE,
        ValidationStatus =
          "CANDIDATE_GENEID_ABSENT_FROM_CURRENT_NCBI"
      )
    )
  }

  stopifnot(
    nrow(hit) == 1
  )

  concordant <-
    identical(
      hit$NCBI_Symbol,
      candidate_symbol
    )

  data.table(
    SourceSymbol = source_symbol,
    CandidateGeneID = candidate_id,
    CandidateSymbol = candidate_symbol,
    NCBIPresent = TRUE,
    NCBICurrentSymbol =
      hit$NCBI_Symbol,
    NCBIGeneType =
      hit$NCBI_GeneType,
    IDSymbolConcordant =
      concordant,
    ValidationStatus =
      if (concordant) {
        "PASS_CURRENT_NCBI_ID_SYMBOL"
      } else {
        "NCBI_CURRENT_SYMBOL_DISAGREES"
      }
  )
}


# ============================================================
# 4. Validate 748 SAFE_ALIAS_CANDIDATE
# ============================================================

safe_alias <- a[
  PreliminaryIdentityStatus ==
    "SAFE_ALIAS_CANDIDATE"
]

cat("Safe alias candidates:", nrow(safe_alias), "\n")

stopifnot(
  nrow(safe_alias) == 748
)

safe_alias_validation <- rbindlist(
  lapply(
    seq_len(nrow(safe_alias)),
    function(i) {

      validate_id_symbol(
        safe_alias$GeneSymbol_raw[i],
        safe_alias$CandidateEntrezID[i],
        safe_alias$CandidateCurrentSymbol[i]
      )
    }
  )
)

cat("\n============================================================\n")
cat("SAFE ALIAS NCBI VALIDATION\n")
cat("============================================================\n")

print(
  safe_alias_validation[
    ,
    .N,
    by = ValidationStatus
  ][order(-N)]
)


# ============================================================
# 5. Validate 197 punctuation candidates
# ============================================================

safe_punct <- u[
  PunctuationRescueStatus ==
    "SAFE_PUNCTUATION_CANDIDATE"
]

cat(
  "\nSafe punctuation candidates:",
  nrow(safe_punct),
  "\n"
)

stopifnot(
  nrow(safe_punct) == 197
)

safe_punct_validation <- rbindlist(
  lapply(
    seq_len(nrow(safe_punct)),
    function(i) {

      validate_id_symbol(
        safe_punct$GeneSymbol_raw[i],
        safe_punct$Punct_EntrezIDs[i],
        safe_punct$Punct_CurrentSymbols[i]
      )
    }
  )
)

cat("\n============================================================\n")
cat("PUNCTUATION CANDIDATE NCBI VALIDATION\n")
cat("============================================================\n")

print(
  safe_punct_validation[
    ,
    .N,
    by = ValidationStatus
  ][order(-N)]
)


# ============================================================
# 6. LOC-style unresolved features
#
# LOC123456 -> candidate GeneID 123456
#
# This is only an identity forensic check.
# The numeric suffix is not automatically accepted unless the
# corresponding current NCBI Gene record exists.
# ============================================================

loc <- u[
  StructuralClass ==
    "LOC_STYLE"
]

stopifnot(
  nrow(loc) == 133
)

loc[
  ,
  LOCNumericGeneID :=
    sub(
      "^LOC",
      "",
      GeneSymbol_raw
    )
]

stopifnot(
  grepl(
    "^[0-9]+$",
    loc$LOCNumericGeneID
  )
)

loc_idx <- match(
  loc$LOCNumericGeneID,
  gi$NCBI_GeneID
)

loc[
  ,
  NCBICurrentPresent :=
    !is.na(loc_idx)
]

loc[
  ,
  NCBICurrentSymbol :=
    fifelse(
      !is.na(loc_idx),
      gi$NCBI_Symbol[loc_idx],
      NA_character_
    )
]

loc[
  ,
  NCBIGeneType :=
    fifelse(
      !is.na(loc_idx),
      gi$NCBI_GeneType[loc_idx],
      NA_character_
    )
]

loc[
  ,
  NCBIDescription :=
    fifelse(
      !is.na(loc_idx),
      gi$NCBI_Description[loc_idx],
      NA_character_
    )
]

loc[
  ,
  LOCStatus :=
    fifelse(
      NCBICurrentPresent,
      "CURRENT_NCBI_GENEID",
      "ABSENT_FROM_CURRENT_NCBI_GENE_INFO"
    )
]

cat("\n============================================================\n")
cat("LOC-STYLE CURRENT NCBI AUDIT\n")
cat("============================================================\n")

print(
  loc[
    ,
    .N,
    by = LOCStatus
  ][order(-N)]
)

cat("\nCurrent LOC mappings (first 50):\n")

print(
  head(
    loc[
      NCBICurrentPresent == TRUE,
      .(
        GeneSymbol_raw,
        LOCNumericGeneID,
        NCBICurrentSymbol,
        NCBIGeneType
      )
    ],
    50
  )
)


# ============================================================
# 7. Remaining 57 simple unresolved symbols
# ============================================================

simple_remaining <- u[
  StructuralClass ==
    "SIMPLE_SYMBOL_LIKE" &
  PunctuationRescueStatus ==
    "NO_PUNCTUATION_RESCUE"
]

stopifnot(
  nrow(simple_remaining) == 57
)

simple_idx <- match(
  simple_remaining$GeneSymbol_raw,
  gi$NCBI_Symbol
)

simple_remaining[
  ,
  NCBIExactCurrentSymbol :=
    !is.na(simple_idx)
]

simple_remaining[
  ,
  NCBICurrentGeneID :=
    fifelse(
      !is.na(simple_idx),
      gi$NCBI_GeneID[simple_idx],
      NA_character_
    )
]

simple_remaining[
  ,
  NCBICurrentSymbol :=
    fifelse(
      !is.na(simple_idx),
      gi$NCBI_Symbol[simple_idx],
      NA_character_
    )
]

simple_remaining[
  ,
  NCBIGeneType :=
    fifelse(
      !is.na(simple_idx),
      gi$NCBI_GeneType[simple_idx],
      NA_character_
    )
]

cat("\n============================================================\n")
cat("REMAINING SIMPLE SYMBOLS: EXACT CURRENT NCBI CHECK\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(simple_remaining),
  "\n"
)

cat(
  "Exact current NCBI symbols:",
  sum(
    simple_remaining$NCBIExactCurrentSymbol
  ),
  "\n"
)

cat(
  "Still unresolved:",
  sum(
    !simple_remaining$NCBIExactCurrentSymbol
  ),
  "\n"
)

if (
  any(
    simple_remaining$NCBIExactCurrentSymbol
  )
) {

  print(
    simple_remaining[
      NCBIExactCurrentSymbol == TRUE,
      .(
        GeneSymbol_raw,
        NCBICurrentGeneID,
        NCBICurrentSymbol,
        NCBIGeneType
      )
    ]
  )
}


# ============================================================
# 8. Resolve the DIRECT_CURRENT_AMBIGUOUS symbol TEC against
#    current NCBI primary symbols
# ============================================================

direct_ambiguous <- a[
  PreliminaryIdentityStatus ==
    "DIRECT_CURRENT_AMBIGUOUS"
]

stopifnot(
  nrow(direct_ambiguous) == 1
)

tec_symbol <-
  direct_ambiguous$GeneSymbol_raw[1]

tec_hits <- gi[
  NCBI_Symbol ==
    tec_symbol
]

cat("\n============================================================\n")
cat("DIRECT AMBIGUOUS: CURRENT NCBI PRIMARY-SYMBOL CHECK\n")
cat("============================================================\n")

cat(
  "Raw symbol:",
  tec_symbol,
  "\n"
)

cat(
  "Current NCBI primary-symbol records:",
  nrow(tec_hits),
  "\n"
)

print(
  tec_hits[
    ,
    .(
      NCBI_GeneID,
      NCBI_Symbol,
      NCBI_GeneType,
      NCBI_Description
    )
  ]
)


# ============================================================
# 9. Alias-ambiguous features:
#    check whether raw symbol itself is now a current NCBI
#    primary symbol.
# ============================================================

alias_ambiguous <- a[
  PreliminaryIdentityStatus ==
    "ALIAS_AMBIGUOUS"
]

stopifnot(
  nrow(alias_ambiguous) == 28
)

amb_idx <- match(
  alias_ambiguous$GeneSymbol_raw,
  gi$NCBI_Symbol
)

alias_ambiguous[
  ,
  NCBIExactCurrentSymbol :=
    !is.na(amb_idx)
]

alias_ambiguous[
  ,
  NCBICurrentGeneID :=
    fifelse(
      !is.na(amb_idx),
      gi$NCBI_GeneID[amb_idx],
      NA_character_
    )
]

alias_ambiguous[
  ,
  NCBICurrentSymbol :=
    fifelse(
      !is.na(amb_idx),
      gi$NCBI_Symbol[amb_idx],
      NA_character_
    )
]

alias_ambiguous[
  ,
  NCBIGeneType :=
    fifelse(
      !is.na(amb_idx),
      gi$NCBI_GeneType[amb_idx],
      NA_character_
    )
]

cat("\n============================================================\n")
cat("ALIAS-AMBIGUOUS: CURRENT NCBI PRIMARY-SYMBOL CHECK\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(alias_ambiguous),
  "\n"
)

cat(
  "Resolved by exact NCBI primary symbol:",
  sum(
    alias_ambiguous$NCBIExactCurrentSymbol
  ),
  "\n"
)

cat(
  "Still ambiguous:",
  sum(
    !alias_ambiguous$NCBIExactCurrentSymbol
  ),
  "\n"
)

if (
  any(
    alias_ambiguous$NCBIExactCurrentSymbol
  )
) {

  print(
    alias_ambiguous[
      NCBIExactCurrentSymbol == TRUE,
      .(
        GeneSymbol_raw,
        NCBICurrentGeneID,
        NCBICurrentSymbol,
        NCBIGeneType
      )
    ]
  )
}


# ============================================================
# 10. Composite feature policy inventory
# ============================================================

composite <- u[
  StructuralClass ==
    "COMPOSITE_UNDERSCORE"
]

stopifnot(
  nrow(composite) == 389
)

cat("\n============================================================\n")
cat("COMPOSITE FEATURES\n")
cat("============================================================\n")

cat(
  "Composite underscore features:",
  nrow(composite),
  "\n"
)

cat(
  "Policy: do not split / rename / aggregate into ordinary",
  " single-gene canonical identities.\n"
)


# ============================================================
# 11. Current resolved-identity accounting
# ============================================================

alias_pass <- safe_alias_validation[
  ValidationStatus ==
    "PASS_CURRENT_NCBI_ID_SYMBOL",
  .N
]

punct_pass <- safe_punct_validation[
  ValidationStatus ==
    "PASS_CURRENT_NCBI_ID_SYMBOL",
  .N
]

loc_current <- loc[
  NCBICurrentPresent == TRUE,
  .N
]

simple_current <- simple_remaining[
  NCBIExactCurrentSymbol == TRUE,
  .N
]

amb_current <- alias_ambiguous[
  NCBIExactCurrentSymbol == TRUE,
  .N
]

tec_current <- nrow(
  tec_hits
)

summary <- data.table(
  Category = c(
    "DIRECT_CURRENT_UNIQUE_orgHs",
    "SAFE_ALIAS_candidate_total",
    "SAFE_ALIAS_NCBI_validated",
    "SAFE_PUNCTUATION_candidate_total",
    "SAFE_PUNCTUATION_NCBI_validated",
    "LOC_STYLE_total",
    "LOC_STYLE_current_NCBI",
    "SIMPLE_remaining_total",
    "SIMPLE_exact_current_NCBI",
    "ALIAS_AMBIGUOUS_total",
    "ALIAS_AMBIGUOUS_exact_current_NCBI",
    "DIRECT_AMBIGUOUS_NCBI_primary_hits",
    "COMPOSITE_total"
  ),

  N = c(
    12208L,
    nrow(safe_alias),
    alias_pass,
    nrow(safe_punct),
    punct_pass,
    nrow(loc),
    loc_current,
    nrow(simple_remaining),
    simple_current,
    nrow(alias_ambiguous),
    amb_current,
    tec_current,
    nrow(composite)
  )
)

cat("\n============================================================\n")
cat("NCBI CROSS-VALIDATION SUMMARY\n")
cat("============================================================\n")

print(summary)


# ============================================================
# 12. Build list that still requires external NCBI Datasets
#     forensic after current gene_info matching
# ============================================================

pending_loc <- loc[
  NCBICurrentPresent == FALSE,
  .(
    QueryType =
      "GENE_ID",
    Query =
      LOCNumericGeneID,
    GeneSymbol_raw,
    Reason =
      "LOC_GENEID_ABSENT_FROM_CURRENT_GENE_INFO"
  )
]

pending_simple <- simple_remaining[
  NCBIExactCurrentSymbol == FALSE,
  .(
    QueryType =
      "SYMBOL",
    Query =
      GeneSymbol_raw,
    GeneSymbol_raw,
    Reason =
      "SIMPLE_SYMBOL_NOT_CURRENT_NCBI_PRIMARY"
  )
]

pending_alias_ambiguous <- alias_ambiguous[
  NCBIExactCurrentSymbol == FALSE,
  .(
    QueryType =
      "SYMBOL",
    Query =
      GeneSymbol_raw,
    GeneSymbol_raw,
    Reason =
      "ALIAS_AMBIGUOUS_NOT_CURRENT_NCBI_PRIMARY"
  )
]

pending_tec <- if (
  nrow(tec_hits) == 1
) {

  data.table(
    QueryType = character(),
    Query = character(),
    GeneSymbol_raw = character(),
    Reason = character()
  )

} else {

  data.table(
    QueryType = "SYMBOL",
    Query = tec_symbol,
    GeneSymbol_raw = tec_symbol,
    Reason =
      "DIRECT_AMBIGUOUS_NOT_UNIQUELY_RESOLVED_BY_NCBI_PRIMARY_SYMBOL"
  )
}

pending <- rbindlist(
  list(
    pending_loc,
    pending_simple,
    pending_alias_ambiguous,
    pending_tec
  ),
  use.names = TRUE,
  fill = TRUE
)

cat("\n============================================================\n")
cat("PENDING EXTERNAL NCBI DATASETS FORENSIC\n")
cat("============================================================\n")

cat(
  "Queries still required:",
  nrow(pending),
  "\n"
)

print(
  pending[
    ,
    .N,
    by = .(
      QueryType,
      Reason
    )
  ]
)


# ============================================================
# 13. Save
# ============================================================

fwrite(
  safe_alias_validation,
  file.path(
    out_dir,
    "GSE93624_safe_alias_NCBI_validation.tsv"
  ),
  sep = "\t"
)

fwrite(
  safe_punct_validation,
  file.path(
    out_dir,
    "GSE93624_safe_punctuation_NCBI_validation.tsv"
  ),
  sep = "\t"
)

fwrite(
  loc,
  file.path(
    out_dir,
    "GSE93624_LOC_style_NCBI_current_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  simple_remaining,
  file.path(
    out_dir,
    "GSE93624_remaining_simple_NCBI_current_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  alias_ambiguous,
  file.path(
    out_dir,
    "GSE93624_alias_ambiguous_NCBI_current_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  tec_hits,
  file.path(
    out_dir,
    "GSE93624_TEC_NCBI_current_records.tsv"
  ),
  sep = "\t"
)

fwrite(
  composite,
  file.path(
    out_dir,
    "GSE93624_composite_features.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary,
  file.path(
    out_dir,
    "GSE93624_NCBI_cross_validation_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  pending,
  file.path(
    out_dir,
    "GSE93624_pending_external_NCBI_queries.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 14. Final
# ============================================================

cat("\nIMPORTANT:\n")
cat(
  "- No expression values were changed.\n",
  "- Candidate mappings were cross-validated against current NCBI Gene records.\n",
  "- Composite features remain unresolved single-row features and are not split.\n",
  "- Collision features remain unmerged.\n",
  "- Only the residual pending set will require external NCBI Datasets queries.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_NCBI_GENE_INFO_CROSS_VALIDATION\n")
cat("============================================================\n")
