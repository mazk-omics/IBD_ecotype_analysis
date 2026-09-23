suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 2E: parse external NCBI Datasets forensic queries
#
# Inputs:
#   - 109-query manifest
#   - one independent Qxxxx.jsonl per query
#
# Parsing backend:
#   1. jsonlite if available
#   2. otherwise NCBI dataformat CLI
#
# IMPORTANT:
#   - NO Python required
#   - NO new NCBI queries performed
#   - query/output order is never assumed
#   - no expression values are modified
#   - no canonicalization is performed yet
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

query_dir <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "ncbi_external_queries"
)

manifest_file <- file.path(
  query_dir,
  "GSE93624_NCBI_external_query_manifest.tsv"
)

status_file <- file.path(
  query_dir,
  "GSE93624_NCBI_external_query_status.tsv"
)

json_dir <- file.path(
  query_dir,
  "raw_jsonl"
)

out_detail <- file.path(
  query_dir,
  "GSE93624_NCBI_external_results_detail.tsv"
)

out_class <- file.path(
  query_dir,
  "GSE93624_NCBI_external_query_classification.tsv"
)

out_summary <- file.path(
  query_dir,
  "GSE93624_NCBI_external_results_summary.tsv"
)

out_reason <- file.path(
  query_dir,
  "GSE93624_NCBI_external_results_by_reason.tsv"
)

out_disagreement <- file.path(
  query_dir,
  "GSE93624_NCBI_symbol_disagreement_external_audit.tsv"
)

out_unresolved <- file.path(
  query_dir,
  "GSE93624_NCBI_external_unresolved_queries.tsv"
)

out_unique <- file.path(
  query_dir,
  "GSE93624_NCBI_external_unique_identity_candidates.tsv"
)


cat("============================================================\n")
cat("GSE93624 NCBI external forensic parser\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load manifest
# ============================================================

manifest <- fread(
  manifest_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  nrow(manifest) == 109,
  !anyDuplicated(manifest$QueryID),
  all(
    manifest$QueryType %chin%
      c("GENE_ID", "SYMBOL")
  )
)

cat("[PASS] External query manifest contains 109 unique QueryIDs\n")


# ============================================================
# 2. Verify query completion status if status file exists
# ============================================================

if (file.exists(status_file)) {

  status <- fread(
    status_file,
    colClasses = "character",
    na.strings = c("", "NA")
  )

  pass_ids <- unique(
    status[
      Status == "PASS",
      QueryID
    ]
  )

  if (length(pass_ids) != 109) {

    stop(
      "Expected 109 successfully completed QueryIDs, found ",
      length(pass_ids)
    )
  }

  stopifnot(
    setequal(
      pass_ids,
      manifest$QueryID
    )
  )

  cat("[PASS] All 109 QueryIDs have successful NCBI query status\n")
}


# ============================================================
# 3. Verify JSONL files
# ============================================================

expected_json <- file.path(
  json_dir,
  paste0(
    manifest$QueryID,
    ".jsonl"
  )
)

missing_json <- expected_json[
  !file.exists(expected_json)
]

if (length(missing_json) > 0) {

  stop(
    "Missing JSONL files:\n",
    paste(
      missing_json,
      collapse = "\n"
    )
  )
}

cat("[PASS] All 109 per-query JSONL files are present\n\n")


# ============================================================
# 4. Select parser backend
# ============================================================

has_jsonlite <- requireNamespace(
  "jsonlite",
  quietly = TRUE
)

dataformat_path <- Sys.which(
  "dataformat"
)

has_dataformat <- nzchar(
  dataformat_path
)

if (has_jsonlite) {

  parser_backend <- "jsonlite"

} else if (has_dataformat) {

  parser_backend <- "NCBI_dataformat"

} else {

  stop(
    paste0(
      "Neither R package 'jsonlite' nor NCBI 'dataformat' ",
      "is available. At least one is required."
    )
  )
}

cat(
  "Parser backend:",
  parser_backend,
  "\n\n"
)


# ============================================================
# 5. Utility functions
# ============================================================

scalar_chr <- function(x) {

  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }

  as.character(x[[1]])
}


collapse_synonyms <- function(x) {

  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }

  x <- as.character(
    unlist(
      x,
      use.names = FALSE
    )
  )

  x <- x[
    !is.na(x) &
    nzchar(x)
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    unique(x),
    collapse = "|"
  )
}


split_synonyms <- function(x) {

  if (
    is.na(x) ||
    !nzchar(x)
  ) {
    return(character())
  }

  # JSON parser stores synonyms with |
  # dataformat normally emits comma-separated synonyms.
  out <- unlist(
    strsplit(
      x,
      "[|,]",
      perl = TRUE
    )
  )

  out <- trimws(out)

  unique(
    out[
      nzchar(out)
    ]
  )
}


empty_return_table <- function() {

  data.table(
    ReturnedGeneID = character(),
    ReturnedSymbol = character(),
    TaxID = character(),
    TaxName = character(),
    GeneType = character(),
    Description = character(),
    Synonyms = character(),
    ReplacedGeneID = character()
  )
}


# ============================================================
# 6A. jsonlite parser
# ============================================================

parse_with_jsonlite <- function(f) {

  lines <- readLines(
    f,
    warn = FALSE,
    encoding = "UTF-8"
  )

  lines <- lines[
    nzchar(
      trimws(lines)
    )
  ]

  if (length(lines) == 0) {
    return(
      empty_return_table()
    )
  }

  rows <- lapply(
    seq_along(lines),
    function(i) {

      z <- tryCatch(
        jsonlite::fromJSON(
          lines[i],
          simplifyVector = FALSE
        ),
        error = function(e) {

          stop(
            basename(f),
            ": invalid JSON on line ",
            i,
            ": ",
            conditionMessage(e)
          )
        }
      )

      data.table(
        ReturnedGeneID =
          scalar_chr(
            z$geneId
          ),

        ReturnedSymbol =
          scalar_chr(
            z$symbol
          ),

        TaxID =
          scalar_chr(
            z$taxId
          ),

        TaxName =
          scalar_chr(
            z$taxname
          ),

        GeneType =
          scalar_chr(
            z$type
          ),

        Description =
          scalar_chr(
            z$description
          ),

        Synonyms =
          collapse_synonyms(
            z$synonyms
          ),

        ReplacedGeneID =
          scalar_chr(
            z$replacedGeneId
          )
      )
    }
  )

  rbindlist(
    rows,
    use.names = TRUE,
    fill = TRUE
  )
}


# ============================================================
# 6B. NCBI dataformat fallback parser
# ============================================================

parse_with_dataformat <- function(f) {

  tmp <- tempfile(
    fileext = ".tsv"
  )

  err <- tempfile(
    fileext = ".log"
  )

  fields <- paste(
    c(
      "gene-id",
      "symbol",
      "tax-id",
      "tax-name",
      "gene-type",
      "description",
      "synonyms",
      "replaced-gene-id"
    ),
    collapse = ","
  )

  exit_code <- system2(
    command = dataformat_path,
    args = c(
      "tsv",
      "gene",
      "--inputfile",
      f,
      "--fields",
      fields,
      "--force"
    ),
    stdout = tmp,
    stderr = err
  )

  if (!identical(exit_code, 0L)) {

    err_txt <- paste(
      readLines(
        err,
        warn = FALSE
      ),
      collapse = "\n"
    )

    unlink(
      c(tmp, err)
    )

    stop(
      basename(f),
      ": dataformat failed:\n",
      err_txt
    )
  }

  if (
    !file.exists(tmp) ||
    file.info(tmp)$size == 0
  ) {

    unlink(
      c(tmp, err)
    )

    return(
      empty_return_table()
    )
  }

  x <- fread(
    tmp,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    colClasses = "character",
    na.strings = c("", "NA"),
    showProgress = FALSE
  )

  unlink(
    c(tmp, err)
  )

  if (nrow(x) == 0) {
    return(
      empty_return_table()
    )
  }

  if (ncol(x) != 8) {

    stop(
      basename(f),
      ": dataformat returned ",
      ncol(x),
      " columns instead of 8."
    )
  }

  # Column order is determined by --fields,
  # so do not depend on human-readable column headers.
  setnames(
    x,
    c(
      "ReturnedGeneID",
      "ReturnedSymbol",
      "TaxID",
      "TaxName",
      "GeneType",
      "Description",
      "Synonyms",
      "ReplacedGeneID"
    )
  )

  x[]
}


parse_one_file <- if (
  parser_backend == "jsonlite"
) {
  parse_with_jsonlite
} else {
  parse_with_dataformat
}


# ============================================================
# 7. Parse all 109 independent queries
# ============================================================

detail_list <- vector(
  "list",
  nrow(manifest)
)

class_list <- vector(
  "list",
  nrow(manifest)
)


for (i in seq_len(nrow(manifest))) {

  m <- manifest[i]

  f <- file.path(
    json_dir,
    paste0(
      m$QueryID,
      ".jsonl"
    )
  )

  x <- parse_one_file(f)


  # ----------------------------------------------------------
  # Retain human records only
  # ----------------------------------------------------------

  if (nrow(x) > 0) {

    x <- x[
      is.na(TaxID) |
      TaxID == "" |
      TaxID == "9606"
    ]
  }


  # ----------------------------------------------------------
  # Match-mode evidence
  # ----------------------------------------------------------

  if (nrow(x) > 0) {

    x[
      ,
      ExactSymbolMatch :=
        !is.na(ReturnedSymbol) &
        toupper(ReturnedSymbol) ==
          toupper(m$Query)
    ]

    x[
      ,
      QueryInSynonyms :=
        vapply(
          Synonyms,
          function(s) {

            m$Query %chin%
              split_synonyms(s)

          },
          logical(1)
        )
    ]

    detail_i <- data.table(
      QueryID =
        m$QueryID,

      QueryType =
        m$QueryType,

      Query =
        m$Query,

      SourceSymbol =
        m$SourceSymbol,

      Reason =
        m$Reason,

      ExpectedCandidateGeneID =
        m$ExpectedCandidateGeneID,

      PreviousCandidateSymbol =
        m$PreviousCandidateSymbol,

      ReturnedGeneID =
        x$ReturnedGeneID,

      ReturnedSymbol =
        x$ReturnedSymbol,

      TaxID =
        x$TaxID,

      TaxName =
        x$TaxName,

      GeneType =
        x$GeneType,

      Description =
        x$Description,

      Synonyms =
        x$Synonyms,

      ReplacedGeneID =
        x$ReplacedGeneID,

      ExactSymbolMatch =
        x$ExactSymbolMatch,

      QueryInSynonyms =
        x$QueryInSynonyms
    )

  } else {

    detail_i <- data.table(
      QueryID = character(),
      QueryType = character(),
      Query = character(),
      SourceSymbol = character(),
      Reason = character(),
      ExpectedCandidateGeneID = character(),
      PreviousCandidateSymbol = character(),
      ReturnedGeneID = character(),
      ReturnedSymbol = character(),
      TaxID = character(),
      TaxName = character(),
      GeneType = character(),
      Description = character(),
      Synonyms = character(),
      ReplacedGeneID = character(),
      ExactSymbolMatch = logical(),
      QueryInSynonyms = logical()
    )
  }

  detail_list[[i]] <- detail_i


  # ==========================================================
  # Query-level unique identities
  # ==========================================================

  returned_geneids <- sort(
    unique(
      detail_i$ReturnedGeneID[
        !is.na(detail_i$ReturnedGeneID) &
        nzchar(detail_i$ReturnedGeneID)
      ]
    )
  )

  returned_symbols <- sort(
    unique(
      detail_i$ReturnedSymbol[
        !is.na(detail_i$ReturnedSymbol) &
        nzchar(detail_i$ReturnedSymbol)
      ]
    )
  )

  n_geneids <- length(
    returned_geneids
  )

  n_symbols <- length(
    returned_symbols
  )


  # ==========================================================
  # Classification: GeneID query
  # ==========================================================

  if (m$QueryType == "GENE_ID") {

    if (nrow(detail_i) == 0) {

      datasets_class <-
        "NO_LIVE_METADATA_FROM_DATASETS"

    } else if (n_geneids == 1) {

      if (
        identical(
          returned_geneids[1],
          m$Query
        )
      ) {

        datasets_class <-
          "CURRENT_LIVE_GENEID"

      } else {

        datasets_class <-
          "REDIRECTED_TO_CURRENT_GENEID"
      }

    } else {

      datasets_class <-
        "MULTIPLE_GENE_RECORDS_FOR_GENEID_QUERY"
    }


  # ==========================================================
  # Classification: symbol query
  # ==========================================================

  } else {

    if (nrow(detail_i) == 0) {

      datasets_class <-
        "NO_LIVE_SYMBOL_MATCH"

    } else if (n_geneids == 1) {

      exact_match <- any(
        detail_i$ExactSymbolMatch,
        na.rm = TRUE
      )

      synonym_match <- any(
        detail_i$QueryInSynonyms,
        na.rm = TRUE
      )

      if (exact_match) {

        datasets_class <-
          "UNIQUE_CURRENT_PRIMARY_SYMBOL"

      } else if (synonym_match) {

        datasets_class <-
          "UNIQUE_CURRENT_SYNONYM"

      } else {

        datasets_class <-
          "UNIQUE_RETURNED_GENE_UNCLEAR_MATCH_MODE"
      }

    } else {

      datasets_class <-
        "AMBIGUOUS_MULTIPLE_CURRENT_GENES"
    }
  }


  # ==========================================================
  # Expected GeneID concordance
  #
  # Used primarily for the 12 current-symbol disagreements.
  # ==========================================================

  expected_concordant <- NA

  if (
    !is.na(m$ExpectedCandidateGeneID) &&
    nzchar(m$ExpectedCandidateGeneID)
  ) {

    expected_concordant <-
      n_geneids == 1 &&
      identical(
        returned_geneids[1],
        m$ExpectedCandidateGeneID
      )
  }


  class_list[[i]] <- data.table(
    QueryID =
      m$QueryID,

    QueryType =
      m$QueryType,

    Query =
      m$Query,

    SourceSymbol =
      m$SourceSymbol,

    Reason =
      m$Reason,

    ExpectedCandidateGeneID =
      m$ExpectedCandidateGeneID,

    PreviousCandidateSymbol =
      m$PreviousCandidateSymbol,

    NReturnedRecords =
      nrow(detail_i),

    NReturnedGeneIDs =
      n_geneids,

    ReturnedGeneIDs =
      if (n_geneids > 0) {
        paste(
          returned_geneids,
          collapse = "|"
        )
      } else {
        NA_character_
      },

    NReturnedSymbols =
      n_symbols,

    ReturnedSymbols =
      if (n_symbols > 0) {
        paste(
          returned_symbols,
          collapse = "|"
        )
      } else {
        NA_character_
      },

    ExpectedGeneIDConcordant =
      expected_concordant,

    DatasetsClass =
      datasets_class
  )
}


# ============================================================
# 8. Combine
# ============================================================

detail <- rbindlist(
  detail_list,
  use.names = TRUE,
  fill = TRUE
)

classification <- rbindlist(
  class_list,
  use.names = TRUE,
  fill = TRUE
)

stopifnot(
  nrow(classification) == 109,
  !anyDuplicated(
    classification$QueryID
  ),
  identical(
    classification$QueryID,
    manifest$QueryID
  )
)

cat(
  "[PASS] Parsed all 109 QueryIDs independently\n"
)

cat(
  "Returned human NCBI gene records:",
  nrow(detail),
  "\n\n"
)


# ============================================================
# 9. Overall classification
# ============================================================

summary <- classification[
  ,
  .N,
  by = DatasetsClass
][
  order(
    -N,
    DatasetsClass
  )
]

cat("============================================================\n")
cat("QUERY CLASSIFICATION\n")
cat("============================================================\n")

print(summary)


# ============================================================
# 10. Classification by original forensic reason
# ============================================================

reason_summary <- classification[
  ,
  .N,
  by = .(
    QueryType,
    Reason,
    DatasetsClass
  )
][
  order(
    QueryType,
    Reason,
    -N,
    DatasetsClass
  )
]

cat("\n============================================================\n")
cat("BY ORIGINAL REASON\n")
cat("============================================================\n")

print(
  reason_summary
)


# ============================================================
# 11. 12 known current-symbol disagreement queries
# ============================================================

disagreement <- classification[
  Reason %chin% c(
    "SAFE_ALIAS_CURRENT_SYMBOL_DISAGREEMENT",
    "PUNCTUATION_CURRENT_SYMBOL_DISAGREEMENT"
  )
]

stopifnot(
  nrow(disagreement) == 12
)

cat("\n============================================================\n")
cat("12 CURRENT-SYMBOL DISAGREEMENT QUERIES\n")
cat("============================================================\n")

print(
  disagreement[
    ,
    .(
      QueryID,
      SourceSymbol,
      ExpectedCandidateGeneID,
      PreviousCandidateSymbol,
      ReturnedGeneIDs,
      ReturnedSymbols,
      ExpectedGeneIDConcordant,
      DatasetsClass
    )
  ]
)


# ============================================================
# 12. TEC forensic
# ============================================================

tec <- classification[
  SourceSymbol == "TEC"
]

stopifnot(
  nrow(tec) == 1
)

cat("\n============================================================\n")
cat("TEC FORENSIC\n")
cat("============================================================\n")

print(
  tec
)

if (nrow(detail[SourceSymbol == "TEC"]) > 0) {

  cat("\nTEC returned NCBI records:\n")

  print(
    detail[
      SourceSymbol == "TEC",
      .(
        ReturnedGeneID,
        ReturnedSymbol,
        GeneType,
        Description,
        Synonyms,
        ReplacedGeneID
      )
    ]
  )
}


# ============================================================
# 13. Unique identity candidates from external queries
#
# Conservative:
# only one returned current GeneID is considered uniquely
# resolved at this stage.
# ============================================================

unique_identity <- classification[
  NReturnedGeneIDs == 1 &
  DatasetsClass %chin% c(
    "CURRENT_LIVE_GENEID",
    "REDIRECTED_TO_CURRENT_GENEID",
    "UNIQUE_CURRENT_PRIMARY_SYMBOL",
    "UNIQUE_CURRENT_SYNONYM",
    "UNIQUE_RETURNED_GENE_UNCLEAR_MATCH_MODE"
  )
]

unique_identity[
  ,
  ResolvedGeneID :=
    ReturnedGeneIDs
]

unique_identity[
  ,
  ResolvedSymbol :=
    ReturnedSymbols
]


cat("\n============================================================\n")
cat("UNIQUE EXTERNAL IDENTITY CANDIDATES\n")
cat("============================================================\n")

cat(
  "Queries with exactly one returned GeneID:",
  nrow(unique_identity),
  "\n"
)

print(
  unique_identity[
    ,
    .N,
    by = .(
      QueryType,
      Reason,
      DatasetsClass
    )
  ][
    order(
      Reason,
      DatasetsClass
    )
  ]
)


# ============================================================
# 14. Residual unresolved/ambiguous external queries
# ============================================================

unresolved_classes <- c(
  "NO_LIVE_METADATA_FROM_DATASETS",
  "NO_LIVE_SYMBOL_MATCH",
  "MULTIPLE_GENE_RECORDS_FOR_GENEID_QUERY",
  "AMBIGUOUS_MULTIPLE_CURRENT_GENES"
)

external_unresolved <- classification[
  DatasetsClass %chin%
    unresolved_classes
]


cat("\n============================================================\n")
cat("RESIDUAL UNRESOLVED / AMBIGUOUS QUERIES\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(external_unresolved),
  "\n"
)

if (nrow(external_unresolved) > 0) {

  print(
    external_unresolved[
      ,
      .(
        QueryID,
        QueryType,
        Query,
        SourceSymbol,
        Reason,
        NReturnedGeneIDs,
        ReturnedGeneIDs,
        ReturnedSymbols,
        DatasetsClass
      )
    ]
  )
}


# ============================================================
# 15. Save outputs
# ============================================================

fwrite(
  detail,
  out_detail,
  sep = "\t",
  na = ""
)

fwrite(
  classification,
  out_class,
  sep = "\t",
  na = ""
)

fwrite(
  summary,
  out_summary,
  sep = "\t"
)

fwrite(
  reason_summary,
  out_reason,
  sep = "\t"
)

fwrite(
  disagreement,
  out_disagreement,
  sep = "\t",
  na = ""
)

fwrite(
  unique_identity,
  out_unique,
  sep = "\t",
  na = ""
)

fwrite(
  external_unresolved,
  out_unresolved,
  sep = "\t",
  na = ""
)


# ============================================================
# 16. Final integrity checks
# ============================================================

stopifnot(
  sum(summary$N) == 109,
  nrow(classification) == 109
)

cat("\n============================================================\n")
cat("INTEGRITY CHECKS\n")
cat("============================================================\n")

cat("[PASS] 109 manifest queries classified exactly once\n")
cat("[PASS] QueryID provenance retained\n")
cat("[PASS] No batch output-order assumption used\n")
cat("[PASS] No expression values modified\n")


cat("\nIMPORTANT:\n")
cat(
  "- This script interprets the already-downloaded NCBI query results only.\n",
  "- No additional NCBI request was made.\n",
  "- A unique NCBI result is still an identity candidate until the final global collision audit.\n",
  "- Ambiguous and no-hit queries remain unresolved.\n",
  "- Final canonical inclusion/exclusion decisions are NOT made in this script.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_NCBI_EXTERNAL_RESULTS_PARSED_IN_R\n")
cat("============================================================\n")
