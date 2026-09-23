suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE93624
# Stage 2E corrected parser
#
# Parse the ALREADY-DOWNLOADED 109 NCBI Datasets JSONL files
# using NCBI's own dataformat utility.
#
# IMPORTANT:
#   - NO new NCBI queries
#   - NO Python
#   - NO direct JSON-field assumptions in R
#   - NO expression modification
#   - NO canonical decision yet
#
# This script supersedes:
#   08_parse_ncbi_external_forensics.R
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

tsv_dir <- file.path(
  query_dir,
  "dataformat_tsv"
)

dir.create(
  tsv_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_detail <- file.path(
  query_dir,
  "GSE93624_NCBI_external_results_detail_CORRECTED.tsv"
)

out_class <- file.path(
  query_dir,
  "GSE93624_NCBI_external_query_classification_CORRECTED.tsv"
)

out_summary <- file.path(
  query_dir,
  "GSE93624_NCBI_external_results_summary_CORRECTED.tsv"
)

out_reason <- file.path(
  query_dir,
  "GSE93624_NCBI_external_results_by_reason_CORRECTED.tsv"
)

out_unique <- file.path(
  query_dir,
  "GSE93624_NCBI_external_unique_identity_candidates_CORRECTED.tsv"
)

out_unresolved <- file.path(
  query_dir,
  "GSE93624_NCBI_external_unresolved_queries_CORRECTED.tsv"
)

out_disagreement <- file.path(
  query_dir,
  "GSE93624_NCBI_symbol_disagreement_external_audit_CORRECTED.tsv"
)


cat("============================================================\n")
cat("GSE93624 corrected NCBI external forensic parser\n")
cat("using official NCBI dataformat\n")
cat("============================================================\n\n")


# ============================================================
# 1. Tool check
# ============================================================

dataformat <- Sys.which(
  "dataformat"
)

if (!nzchar(dataformat)) {
  stop(
    "NCBI dataformat is not available in PATH."
  )
}

cat(
  "dataformat:",
  dataformat,
  "\n"
)


version_out <- tryCatch(
  system2(
    dataformat,
    "version",
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) {
    "version unavailable"
  }
)

cat(
  "dataformat version:",
  paste(version_out, collapse = " "),
  "\n\n"
)


# ============================================================
# 2. Manifest
# ============================================================

manifest <- fread(
  manifest_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  nrow(manifest) == 109,
  !anyDuplicated(
    manifest$QueryID
  ),
  all(
    manifest$QueryType %chin%
      c(
        "GENE_ID",
        "SYMBOL"
      )
  )
)

cat(
  "[PASS] Manifest contains 109 unique queries\n"
)


# ============================================================
# 3. Verify previous network-query success
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

  stopifnot(
    length(pass_ids) == 109,
    setequal(
      pass_ids,
      manifest$QueryID
    )
  )

  cat(
    "[PASS] Previous NCBI network queries were successful 109/109\n"
  )
}


# ============================================================
# 4. Verify JSONL resources
# ============================================================

json_files <- file.path(
  json_dir,
  paste0(
    manifest$QueryID,
    ".jsonl"
  )
)

stopifnot(
  all(
    file.exists(json_files)
  )
)

cat(
  "[PASS] All 109 original JSONL responses are present\n\n"
)


# ============================================================
# 5. Parse one NCBI gene JSONL using dataformat
# ============================================================

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


empty_gene_table <- function() {

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


parse_gene_jsonl <- function(
  input_file,
  query_id
) {

  output_file <- file.path(
    tsv_dir,
    paste0(
      query_id,
      ".tsv"
    )
  )

  error_file <- file.path(
    tsv_dir,
    paste0(
      query_id,
      ".stderr.log"
    )
  )


  # ----------------------------------------------------------
  # Let official NCBI dataformat parse the JSONL.
  # ----------------------------------------------------------

  exit_code <- system2(
    command = dataformat,
    args = c(
      "tsv",
      "gene",
      "--inputfile",
      input_file,
      "--fields",
      fields,
      "--force"
    ),
    stdout = output_file,
    stderr = error_file
  )


  if (!identical(exit_code, 0L)) {

    err <- if (
      file.exists(error_file)
    ) {
      paste(
        readLines(
          error_file,
          warn = FALSE
        ),
        collapse = "\n"
      )
    } else {
      ""
    }

    stop(
      query_id,
      ": dataformat failed with exit code ",
      exit_code,
      "\n",
      err
    )
  }


  if (
    !file.exists(output_file) ||
    file.info(output_file)$size == 0
  ) {

    return(
      empty_gene_table()
    )
  }


  lines <- readLines(
    output_file,
    warn = FALSE,
    encoding = "UTF-8"
  )

  if (length(lines) <= 1) {

    # Header only = no returned gene record
    return(
      empty_gene_table()
    )
  }


  x <- fread(
    output_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    colClasses = "character",
    na.strings = c("", "NA"),
    showProgress = FALSE
  )


  if (nrow(x) == 0) {

    return(
      empty_gene_table()
    )
  }


  if (ncol(x) != 8) {

    stop(
      query_id,
      ": expected 8 dataformat columns but received ",
      ncol(x),
      "\nHeaders: ",
      paste(
        names(x),
        collapse = " | "
      )
    )
  }


  # ----------------------------------------------------------
  # Do not depend on human-readable dataformat headers.
  #
  # Column order follows --fields exactly.
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # HARD CHECK:
  #
  # Any actual returned NCBI Gene record MUST have GeneID.
  #
  # This directly protects us against the failure in old 08.
  # ----------------------------------------------------------

  bad <- x[
    !is.na(ReturnedSymbol) &
    ReturnedSymbol != "" &
    (
      is.na(ReturnedGeneID) |
      ReturnedGeneID == ""
    )
  ]

  if (nrow(bad) > 0) {

    stop(
      query_id,
      ": dataformat returned gene symbols without GeneID. ",
      "This should not occur for a valid GeneDescriptor."
    )
  }


  x[]
}


# ============================================================
# 6. Synonym helper
#
# dataformat renders multi-valued synonyms comma-separated.
# ============================================================

split_synonyms <- function(x) {

  if (
    is.na(x) ||
    !nzchar(x)
  ) {
    return(
      character()
    )
  }

  z <- unlist(
    strsplit(
      x,
      ",",
      fixed = TRUE
    )
  )

  z <- trimws(z)

  unique(
    z[
      nzchar(z)
    ]
  )
}


# ============================================================
# 7. Parse + classify all 109 QueryIDs
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

  json_file <- file.path(
    json_dir,
    paste0(
      m$QueryID,
      ".jsonl"
    )
  )


  x <- parse_gene_jsonl(
    json_file,
    m$QueryID
  )


  # ----------------------------------------------------------
  # GeneID queries are globally unique.
  # Symbol queries were explicitly performed with --taxon human.
  #
  # Keep human records; explicitly retain unknown taxID only
  # for audit rather than silently deleting them.
  # ----------------------------------------------------------

  if (nrow(x) > 0) {

    nonhuman <- x[
      !is.na(TaxID) &
      TaxID != "" &
      TaxID != "9606"
    ]

    if (nrow(nonhuman) > 0) {

      warning(
        m$QueryID,
        ": non-human returned records detected."
      )
    }

    x <- x[
      is.na(TaxID) |
      TaxID == "" |
      TaxID == "9606"
    ]
  }


  # ----------------------------------------------------------
  # Match evidence
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
  # Unique returned identities
  # ==========================================================

  ids <- sort(
    unique(
      detail_i$ReturnedGeneID[
        !is.na(detail_i$ReturnedGeneID) &
        detail_i$ReturnedGeneID != ""
      ]
    )
  )

  syms <- sort(
    unique(
      detail_i$ReturnedSymbol[
        !is.na(detail_i$ReturnedSymbol) &
        detail_i$ReturnedSymbol != ""
      ]
    )
  )

  n_ids <- length(ids)
  n_syms <- length(syms)


  # ==========================================================
  # Classification: GeneID query
  # ==========================================================

  if (m$QueryType == "GENE_ID") {

    if (nrow(detail_i) == 0) {

      cls <-
        "NO_LIVE_METADATA_FROM_DATASETS"

    } else if (n_ids == 1) {

      if (
        identical(
          ids[1],
          m$Query
        )
      ) {

        cls <-
          "CURRENT_LIVE_GENEID"

      } else {

        cls <-
          "REDIRECTED_TO_CURRENT_GENEID"
      }

    } else {

      cls <-
        "MULTIPLE_GENE_RECORDS_FOR_GENEID_QUERY"
    }


  # ==========================================================
  # Classification: symbol query
  # ==========================================================

  } else {

    if (nrow(detail_i) == 0) {

      cls <-
        "NO_LIVE_SYMBOL_MATCH"

    } else if (n_ids == 1) {

      exact_match <- any(
        detail_i$ExactSymbolMatch,
        na.rm = TRUE
      )

      synonym_match <- any(
        detail_i$QueryInSynonyms,
        na.rm = TRUE
      )


      if (exact_match) {

        cls <-
          "UNIQUE_CURRENT_PRIMARY_SYMBOL"

      } else if (synonym_match) {

        cls <-
          "UNIQUE_CURRENT_SYNONYM"

      } else {

        cls <-
          "UNIQUE_RETURNED_GENE_UNCLEAR_MATCH_MODE"
      }

    } else {

      cls <-
        "AMBIGUOUS_MULTIPLE_CURRENT_GENES"
    }
  }


  # ==========================================================
  # Expected-candidate concordance
  # ==========================================================

  expected_concordant <- NA

  if (
    !is.na(
      m$ExpectedCandidateGeneID
    ) &&
    m$ExpectedCandidateGeneID != ""
  ) {

    expected_concordant <-
      n_ids == 1 &&
      identical(
        ids[1],
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
      n_ids,

    ReturnedGeneIDs =
      if (n_ids > 0) {
        paste(
          ids,
          collapse = "|"
        )
      } else {
        NA_character_
      },

    NReturnedSymbols =
      n_syms,

    ReturnedSymbols =
      if (n_syms > 0) {
        paste(
          syms,
          collapse = "|"
        )
      } else {
        NA_character_
      },

    ExpectedGeneIDConcordant =
      expected_concordant,

    DatasetsClass =
      cls
  )
}


# ============================================================
# 8. Assemble
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


# ============================================================
# 9. Critical corrected-parser integrity check
#
# Any returned Gene record must now carry GeneID.
# ============================================================

if (nrow(detail) > 0) {

  stopifnot(
    all(
      !is.na(
        detail$ReturnedGeneID
      )
    ),
    all(
      detail$ReturnedGeneID != ""
    )
  )
}

cat(
  "[PASS] Every returned NCBI Gene record has a GeneID\n"
)

cat(
  "[PASS] Parsed all 109 QueryIDs independently\n"
)

cat(
  "Returned human gene records:",
  nrow(detail),
  "\n\n"
)


# ============================================================
# 10. Overall classification
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
cat("QUERY CLASSIFICATION - CORRECTED\n")
cat("============================================================\n")

print(summary)


# ============================================================
# 11. By original reason
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
cat("BY ORIGINAL REASON - CORRECTED\n")
cat("============================================================\n")

print(
  reason_summary
)


# ============================================================
# 12. 12 known symbol disagreements
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
cat("12 CURRENT-SYMBOL DISAGREEMENTS - CORRECTED\n")
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
# 13. TEC
# ============================================================

tec <- classification[
  SourceSymbol == "TEC"
]

stopifnot(
  nrow(tec) == 1
)


cat("\n============================================================\n")
cat("TEC FORENSIC - CORRECTED\n")
cat("============================================================\n")

print(tec)


if (nrow(detail[SourceSymbol == "TEC"]) > 0) {

  cat("\nTEC returned records:\n")

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
# 14. Unique identity candidates
# ============================================================

unique_classes <- c(
  "CURRENT_LIVE_GENEID",
  "REDIRECTED_TO_CURRENT_GENEID",
  "UNIQUE_CURRENT_PRIMARY_SYMBOL",
  "UNIQUE_CURRENT_SYNONYM",
  "UNIQUE_RETURNED_GENE_UNCLEAR_MATCH_MODE"
)


unique_identity <- classification[
  NReturnedGeneIDs == 1 &
  DatasetsClass %chin%
    unique_classes
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
  "Rows:",
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
      QueryType,
      Reason,
      DatasetsClass
    )
  ]
)


# ============================================================
# 15. Residual unresolved
# ============================================================

unresolved_classes <- c(
  "NO_LIVE_METADATA_FROM_DATASETS",
  "NO_LIVE_SYMBOL_MATCH",
  "MULTIPLE_GENE_RECORDS_FOR_GENEID_QUERY",
  "AMBIGUOUS_MULTIPLE_CURRENT_GENES"
)


unresolved <- classification[
  DatasetsClass %chin%
    unresolved_classes
]


cat("\n============================================================\n")
cat("RESIDUAL UNRESOLVED / AMBIGUOUS - CORRECTED\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(unresolved),
  "\n"
)

if (nrow(unresolved) > 0) {

  print(
    unresolved[
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
# 16. Save corrected outputs
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
  unique_identity,
  out_unique,
  sep = "\t",
  na = ""
)

fwrite(
  unresolved,
  out_unresolved,
  sep = "\t",
  na = ""
)

fwrite(
  disagreement,
  out_disagreement,
  sep = "\t",
  na = ""
)


# ============================================================
# 17. Final integrity
# ============================================================

stopifnot(
  sum(
    summary$N
  ) == 109
)


cat("\n============================================================\n")
cat("FINAL INTEGRITY CHECKS\n")
cat("============================================================\n")

cat("[PASS] 109 queries classified exactly once\n")
cat("[PASS] Every returned gene record contains NCBI GeneID\n")
cat("[PASS] Official NCBI dataformat used for JSONL parsing\n")
cat("[PASS] No network query repeated\n")
cat("[PASS] No expression value modified\n")


cat("\nIMPORTANT:\n")
cat(
  "- Results from the previous 08 parser must NOT be used.\n",
  "- Files ending in _CORRECTED.tsv are authoritative for this external forensic stage.\n",
  "- Unique identities remain subject to one final global collision audit before canonical inclusion.\n",
  "- No feature has yet been renamed, merged, or excluded by this script.\n",
  sep = ""
)


cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_NCBI_EXTERNAL_RESULTS_CORRECTLY_PARSED\n")
cat("============================================================\n")
