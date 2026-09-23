suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Build final external NCBI Datasets query manifest
#
# Query set:
#   A. 97 residual unresolved queries from Stage 2C
#   B. 10 safe-alias current-symbol disagreements
#   C. 2 punctuation current-symbol disagreements
#
# One QueryID = one independent NCBI request.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

harm_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

out_dir <- file.path(
  harm_dir,
  "ncbi_external_queries"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE93624 build external NCBI query manifest\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load Stage 2C residual 97 queries
# ============================================================

pending_file <- file.path(
  harm_dir,
  "GSE93624_pending_external_NCBI_queries.tsv"
)

alias_validation_file <- file.path(
  harm_dir,
  "GSE93624_safe_alias_NCBI_validation.tsv"
)

punct_validation_file <- file.path(
  harm_dir,
  "GSE93624_safe_punctuation_NCBI_validation.tsv"
)


pending <- fread(
  pending_file,
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

stopifnot(
  nrow(pending) == 97,
  nrow(alias_val) == 748,
  nrow(punct_val) == 197
)


# ============================================================
# 2. Normalize residual queries
# ============================================================

q_pending <- pending[
  ,
  .(
    QueryType,
    Query,
    SourceSymbol =
      GeneSymbol_raw,
    Reason,
    ExpectedCandidateGeneID =
      NA_character_,
    PreviousCandidateSymbol =
      NA_character_
  )
]


# ============================================================
# 3. Add 10 safe-alias GeneID disagreements
# ============================================================

alias_disagree <- alias_val[
  ValidationStatus ==
    "NCBI_CURRENT_SYMBOL_DISAGREES"
]

stopifnot(
  nrow(alias_disagree) == 10
)

q_alias <- alias_disagree[
  ,
  .(
    QueryType =
      "GENE_ID",

    Query =
      CandidateGeneID,

    SourceSymbol =
      SourceSymbol,

    Reason =
      "SAFE_ALIAS_CURRENT_SYMBOL_DISAGREEMENT",

    ExpectedCandidateGeneID =
      CandidateGeneID,

    PreviousCandidateSymbol =
      CandidateSymbol
  )
]


# ============================================================
# 4. Add 2 punctuation GeneID disagreements
# ============================================================

punct_disagree <- punct_val[
  ValidationStatus ==
    "NCBI_CURRENT_SYMBOL_DISAGREES"
]

stopifnot(
  nrow(punct_disagree) == 2
)

q_punct <- punct_disagree[
  ,
  .(
    QueryType =
      "GENE_ID",

    Query =
      CandidateGeneID,

    SourceSymbol =
      SourceSymbol,

    Reason =
      "PUNCTUATION_CURRENT_SYMBOL_DISAGREEMENT",

    ExpectedCandidateGeneID =
      CandidateGeneID,

    PreviousCandidateSymbol =
      CandidateSymbol
  )
]


# ============================================================
# 5. Combine
# ============================================================

manifest <- rbindlist(
  list(
    q_pending,
    q_alias,
    q_punct
  ),
  use.names = TRUE,
  fill = TRUE
)

stopifnot(
  nrow(manifest) == 109
)


# ------------------------------------------------------------
# No exact duplicate QueryType + Query + SourceSymbol + Reason
# ------------------------------------------------------------

stopifnot(
  !anyDuplicated(
    manifest[
      ,
      .(
        QueryType,
        Query,
        SourceSymbol,
        Reason
      )
    ]
  )
)


# ============================================================
# 6. Validate query types
# ============================================================

stopifnot(
  all(
    manifest$QueryType %chin%
      c(
        "GENE_ID",
        "SYMBOL"
      )
  )
)

stopifnot(
  all(
    !is.na(manifest$Query) &
    manifest$Query != ""
  )
)

# GeneID queries must be numeric.
stopifnot(
  all(
    grepl(
      "^[0-9]+$",
      manifest[
        QueryType == "GENE_ID",
        Query
      ]
    )
  )
)


# ============================================================
# 7. Stable QueryID
# ============================================================

manifest[
  ,
  QueryID :=
    sprintf(
      "Q%04d",
      seq_len(.N)
    )
]

setcolorder(
  manifest,
  c(
    "QueryID",
    "QueryType",
    "Query",
    "SourceSymbol",
    "Reason",
    "ExpectedCandidateGeneID",
    "PreviousCandidateSymbol"
  )
)


# ============================================================
# 8. Summary
# ============================================================

cat("Total external queries:", nrow(manifest), "\n\n")

cat("By query type:\n")

print(
  manifest[
    ,
    .N,
    by = QueryType
  ]
)

cat("\nBy reason:\n")

print(
  manifest[
    ,
    .N,
    by = .(
      QueryType,
      Reason
    )
  ][order(QueryType, Reason)]
)


# ============================================================
# 9. Save
# ============================================================

fwrite(
  manifest,
  file.path(
    out_dir,
    "GSE93624_NCBI_external_query_manifest.tsv"
  ),
  sep = "\t",
  na = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_NCBI_EXTERNAL_QUERY_MANIFEST_BUILT\n")
cat("============================================================\n")
