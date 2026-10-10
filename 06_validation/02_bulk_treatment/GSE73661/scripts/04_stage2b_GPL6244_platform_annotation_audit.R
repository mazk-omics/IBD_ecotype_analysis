suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 2B v2
# GPL6244 official gene_assignment audit
#
# GPL6244 does NOT contain dedicated SYMBOL / ENTREZ columns.
# Gene information is encoded in multipart "gene_assignment":
#
# assignment1 /// assignment2 ...
#
# Each assignment is expected to contain fields separated by:
# " // "
#
# Relevant fields:
#   field 1 = transcript/public accession
#   field 2 = gene symbol
#   field 3 = gene description
#   field 4 = cytoband / other assignment metadata
#   field 5 = Entrez Gene ID, when available
#
# Goals:
#   1. Extract GPL6244 directly from GEO family SOFT
#   2. Align all 33,252 expression transcript clusters
#   3. Parse official gene_assignment
#   4. Compare GPL6244 vs Stage-2A ChipDb mappings
#   5. Audit unresolved / ambiguous probes
#
# NO expression modification
# NO probe collapse
# NO final inclusion/exclusion decisions
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

soft_file <- file.path(
  root,
  "00_provenance/geo_downloads",
  "GSE73661_family.soft.gz"
)

expr_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE73661_authoritative_log2_RMA_33252x178.rds"
)

stage2a_file <- file.path(
  root,
  "prepared/02_gene_harmonization",
  "GSE73661_stage2a_probe_gene_mapping_audit.tsv"
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

  cat("[PASS] ", message, "\n", sep = "")
}


clean_chr <- function(x) {

  x <- trimws(as.character(x))

  x[
    is.na(x) |
      x == "" |
      x %in% c(
        "---",
        "NA",
        "N/A",
        "na",
        "n/a",
        "NULL",
        "null"
      )
  ] <- NA_character_

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


split_semicolon_set <- function(x) {

  x <- clean_chr(x)

  if (length(x) == 0L || is.na(x[1])) {
    return(character())
  }

  z <- unlist(
    strsplit(
      x[1],
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


set_equal_semicolon <- function(x, y) {

  a <- split_semicolon_set(x)
  b <- split_semicolon_set(y)

  identical(
    a,
    b
  )
}


# ------------------------------------------------------------
# Robust parser for one GPL6244 gene_assignment cell
# ------------------------------------------------------------

parse_gene_assignment <- function(x) {

  x <- clean_chr(x)

  if (length(x) == 0L || is.na(x[1])) {

    return(
      list(
        NAssignments = 0L,
        Accessions = NA_character_,
        Symbols = NA_character_,
        Entrez = NA_character_,
        N_Symbols = 0L,
        N_Entrez = 0L,
        MalformedAssignments = 0L
      )
    )
  }

  # Multiple transcript/gene assignments are separated by ///
  assignments <- strsplit(
    x[1],
    "\\s*///\\s*",
    perl = TRUE
  )[[1]]

  assignments <- trimws(assignments)

  assignments <- assignments[
    nzchar(assignments)
  ]

  accessions <- character()
  symbols <- character()
  entrez <- character()

  malformed <- 0L

  for (a in assignments) {

    # Fields within an assignment are separated by //
    fields <- strsplit(
      a,
      "\\s*//\\s*",
      perl = TRUE
    )[[1]]

    fields <- trimws(fields)

    if (length(fields) < 2L) {
      malformed <- malformed + 1L
      next
    }

    # field 1: transcript/accession
    acc <- clean_chr(fields[1])

    if (!is.na(acc)) {
      accessions <- c(
        accessions,
        acc
      )
    }


    # field 2: gene symbol
    sym <- clean_chr(fields[2])

    if (!is.na(sym)) {
      symbols <- c(
        symbols,
        sym
      )
    }


    # field 5: Entrez Gene ID when supplied.
    #
    # We additionally require a numeric token because "---"
    # or other annotation text must not become an Entrez ID.
    if (length(fields) >= 5L) {

      e <- clean_chr(fields[5])

      if (
        !is.na(e) &&
        grepl(
          "^[0-9]+$",
          e
        )
      ) {

        entrez <- c(
          entrez,
          e
        )
      }
    }
  }


  accessions <- sort(
    unique(accessions)
  )

  symbols <- sort(
    unique(symbols)
  )

  entrez <- sort(
    unique(entrez)
  )


  list(
    NAssignments =
      length(assignments),

    Accessions =
      if (length(accessions) == 0L) {
        NA_character_
      } else {
        paste(
          accessions,
          collapse = ";"
        )
      },

    Symbols =
      if (length(symbols) == 0L) {
        NA_character_
      } else {
        paste(
          symbols,
          collapse = ";"
        )
      },

    Entrez =
      if (length(entrez) == 0L) {
        NA_character_
      } else {
        paste(
          entrez,
          collapse = ";"
        )
      },

    N_Symbols =
      length(symbols),

    N_Entrez =
      length(entrez),

    MalformedAssignments =
      malformed
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 2B v2: GPL6244 gene_assignment audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Inputs
# ============================================================

check(
  file.exists(soft_file),
  "GEO family SOFT exists"
)

check(
  file.exists(expr_file),
  "Frozen Stage-1 authoritative expression exists"
)

check(
  file.exists(stage2a_file),
  "Frozen Stage-2A mapping audit exists"
)


expr <- readRDS(
  expr_file
)

stage2a <- fread(
  stage2a_file
)


check(
  identical(
    dim(expr),
    c(
      33252L,
      178L
    )
  ),
  "Authoritative expression = 33,252 x 178"
)


source_ids <- as.character(
  rownames(expr)
)


check(
  nrow(stage2a) == 33252L,
  "Stage-2A audit contains 33,252 rows"
)

check(
  identical(
    as.character(
      stage2a$SourceProbeID
    ),
    source_ids
  ),
  "Stage-2A probe order matches expression"
)


# ============================================================
# 2. Extract GPL6244 platform table
# ============================================================

awk_file <- tempfile(
  fileext = ".awk"
)

platform_tmp <- tempfile(
  fileext = ".tsv"
)


writeLines(
  c(
    'tolower($0) ~ /^!platform_table_begin/ { in_table=1; next }',
    'tolower($0) ~ /^!platform_table_end/   { in_table=0; next }',
    'in_table { print }'
  ),
  awk_file
)


status <- system(
  paste(
    "zcat",
    shQuote(soft_file),
    "| awk -f",
    shQuote(awk_file),
    ">",
    shQuote(platform_tmp)
  )
)


check(
  status == 0,
  "GPL6244 platform table extracted"
)


platform <- fread(
  platform_tmp,
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = FALSE,
  na.strings = c(
    "",
    "---",
    "NA",
    "N/A"
  )
)


unlink(
  c(
    awk_file,
    platform_tmp
  )
)


cat("\n============================================================\n")
cat("GPL6244 PLATFORM STRUCTURE\n")
cat("============================================================\n")

cat(
  "Rows: ",
  nrow(platform),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(platform),
  "\n",
  sep = ""
)

print(
  data.table(
    Index =
      seq_along(
        names(platform)
      ),

    Column =
      names(platform)
  )
)


check(
  nrow(platform) == 33297L,
  "GPL6244 contains 33,297 platform rows"
)

check(
  "ID" %in%
    names(platform),
  "GPL6244 contains ID column"
)

check(
  "gene_assignment" %in%
    names(platform),
  "GPL6244 contains gene_assignment column"
)

check(
  "category" %in%
    names(platform),
  "GPL6244 contains category column"
)


platform[
  ,
  ID :=
    as.character(ID)
]


check(
  uniqueN(platform$ID) ==
    nrow(platform),
  "GPL6244 platform IDs are unique"
)


# ============================================================
# 3. Expression ↔ platform universe
# ============================================================

source_present <- source_ids %in%
  platform$ID


check(
  sum(source_present) == 33252L,
  "All 33,252 expression IDs occur in GPL6244"
)


extra_platform_ids <- setdiff(
  platform$ID,
  source_ids
)


cat(
  "[INFO] GPL6244 features absent from expression matrix: ",
  length(extra_platform_ids),
  "\n",
  sep = ""
)


check(
  length(extra_platform_ids) == 45L,
  "Exactly 45 GPL6244 features are unused by this expression matrix"
)


platform_source <- platform[
  match(
    source_ids,
    ID
  )
]


check(
  identical(
    platform_source$ID,
    source_ids
  ),
  "GPL6244 rows aligned exactly to expression order"
)


# ============================================================
# 4. Platform category audit
# ============================================================

category_summary <- platform_source[
  ,
  .N,
  by = category
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("GPL6244 CATEGORY DISTRIBUTION — EXPRESSION FEATURES\n")
cat("============================================================\n")

print(
  category_summary
)


# ============================================================
# 5. Parse all gene_assignment fields
# ============================================================

cat("\nParsing 33,252 gene_assignment fields...\n")


parsed_list <- lapply(
  platform_source$gene_assignment,
  parse_gene_assignment
)


parsed <- rbindlist(
  lapply(
    parsed_list,
    as.data.table
  )
)


check(
  nrow(parsed) == 33252L,
  "All 33,252 gene_assignment fields parsed"
)


platform_parsed <- data.table(
  SourceRow =
    seq_along(
      source_ids
    ),

  SourceProbeID =
    source_ids,

  Category =
    platform_source$category,

  GeneAssignmentRaw =
    platform_source$gene_assignment,

  N_Assignments =
    parsed$NAssignments,

  PlatformAccessions =
    parsed$Accessions,

  PlatformSymbols =
    parsed$Symbols,

  N_PlatformSymbols =
    parsed$N_Symbols,

  PlatformEntrez =
    parsed$Entrez,

  N_PlatformEntrez =
    parsed$N_Entrez,

  MalformedAssignments =
    parsed$MalformedAssignments
)


check(
  identical(
    platform_parsed$SourceProbeID,
    source_ids
  ),
  "Parsed platform annotation preserves source order"
)


# ============================================================
# 6. Parser QC
# ============================================================

parser_summary <- data.table(
  Metric = c(
    "Source_transcript_clusters",
    "No_gene_assignment",
    "Single_assignment",
    "Multiple_assignments",
    "No_symbol",
    "Unique_symbol",
    "Multiple_symbols",
    "No_Entrez",
    "Unique_Entrez",
    "Multiple_Entrez",
    "Rows_with_malformed_assignment"
  ),

  N = c(
    nrow(
      platform_parsed
    ),

    sum(
      platform_parsed$N_Assignments ==
        0L
    ),

    sum(
      platform_parsed$N_Assignments ==
        1L
    ),

    sum(
      platform_parsed$N_Assignments >
        1L
    ),

    sum(
      platform_parsed$N_PlatformSymbols ==
        0L
    ),

    sum(
      platform_parsed$N_PlatformSymbols ==
        1L
    ),

    sum(
      platform_parsed$N_PlatformSymbols >
        1L
    ),

    sum(
      platform_parsed$N_PlatformEntrez ==
        0L
    ),

    sum(
      platform_parsed$N_PlatformEntrez ==
        1L
    ),

    sum(
      platform_parsed$N_PlatformEntrez >
        1L
    ),

    sum(
      platform_parsed$MalformedAssignments >
        0L
    )
  )
)


cat("\n============================================================\n")
cat("GPL6244 gene_assignment PARSER SUMMARY\n")
cat("============================================================\n")

print(
  parser_summary
)


# ============================================================
# 7. Inspect parser examples
# ============================================================

cat("\n============================================================\n")
cat("PARSER EXAMPLES\n")
cat("============================================================\n")


example_unique <- platform_parsed[
  N_PlatformEntrez == 1L &
    N_PlatformSymbols == 1L
][1:min(.N, 3L)]


example_multi <- platform_parsed[
  N_PlatformSymbols > 1L |
    N_PlatformEntrez > 1L
][1:min(.N, 3L)]


cat("\nUnique examples:\n")

print(
  example_unique[
    ,
    .(
      SourceProbeID,
      PlatformSymbols,
      PlatformEntrez,
      N_Assignments
    )
  ]
)


cat("\nMultiple-assignment examples:\n")

print(
  example_multi[
    ,
    .(
      SourceProbeID,
      PlatformSymbols,
      PlatformEntrez,
      N_Assignments
    )
  ]
)


# ============================================================
# 8. Combine with Stage 2A
# ============================================================

comparison <- data.table(
  SourceRow =
    seq_along(
      source_ids
    ),

  SourceProbeID =
    source_ids,

  Stage2AStatus =
    stage2a$MappingStatus,

  ChipEntrez =
    clean_chr(
      stage2a$ChipEntrez
    ),

  N_ChipEntrez =
    stage2a$N_ChipEntrez,

  ChipSymbol =
    clean_chr(
      stage2a$ChipSymbol
    ),

  N_ChipSymbol =
    stage2a$N_ChipSymbol,

  PlatformCategory =
    platform_parsed$Category,

  N_PlatformAssignments =
    platform_parsed$N_Assignments,

  PlatformEntrez =
    platform_parsed$PlatformEntrez,

  N_PlatformEntrez =
    platform_parsed$N_PlatformEntrez,

  PlatformSymbols =
    platform_parsed$PlatformSymbols,

  N_PlatformSymbols =
    platform_parsed$N_PlatformSymbols
)


comparison[
  ,
  EntrezSetExact :=
    mapply(
      set_equal_semicolon,
      ChipEntrez,
      PlatformEntrez
    )
]


comparison[
  ,
  SymbolSetExact :=
    mapply(
      set_equal_semicolon,
      ChipSymbol,
      PlatformSymbols
    )
]


comparison[
  ,
  EntrezConcordance :=
    fcase(

      N_ChipEntrez == 0L &
        N_PlatformEntrez == 0L,
      "BOTH_NONE",

      N_ChipEntrez == 0L &
        N_PlatformEntrez > 0L,
      "PLATFORM_ONLY",

      N_ChipEntrez > 0L &
        N_PlatformEntrez == 0L,
      "CHIPDB_ONLY",

      EntrezSetExact,
      "EXACT_SET",

      default =
        "DIFFERENT_SET"
    )
]


comparison[
  ,
  SymbolConcordance :=
    fcase(

      N_ChipSymbol == 0L &
        N_PlatformSymbols == 0L,
      "BOTH_NONE",

      N_ChipSymbol == 0L &
        N_PlatformSymbols > 0L,
      "PLATFORM_ONLY",

      N_ChipSymbol > 0L &
        N_PlatformSymbols == 0L,
      "CHIPDB_ONLY",

      SymbolSetExact,
      "EXACT_SET",

      default =
        "DIFFERENT_SET"
    )
]


# ============================================================
# 9. Overall concordance
# ============================================================

entrez_concordance <- comparison[
  ,
  .N,
  by = EntrezConcordance
][
  order(
    -N
  )
]


symbol_concordance <- comparison[
  ,
  .N,
  by = SymbolConcordance
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("CHIPDB vs GPL6244 — ENTREZ CONCORDANCE\n")
cat("============================================================\n")

print(
  entrez_concordance
)


cat("\n============================================================\n")
cat("CHIPDB vs GPL6244 — SYMBOL CONCORDANCE\n")
cat("============================================================\n")

print(
  symbol_concordance
)


# ============================================================
# 10. Stage-2A confident mappings
# ============================================================

confident_audit <- comparison[
  Stage2AStatus ==
    "UNIQUE_ENTREZ_CURRENT_SYMBOL",
  .N,
  by = .(
    EntrezConcordance,
    SymbolConcordance
  )
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("STAGE-2A CONFIDENT MAPPINGS — GPL6244 CROSS-CHECK\n")
cat("============================================================\n")

print(
  confident_audit
)


# ============================================================
# 11. Stage-2A NO_ENTREZ
# ============================================================

no_entrez_audit <- comparison[
  Stage2AStatus ==
    "NO_ENTREZ_MAPPING",
  .N,
  by = .(
    N_PlatformAssignments,
    N_PlatformSymbols,
    N_PlatformEntrez,
    EntrezConcordance,
    SymbolConcordance
  )
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("STAGE-2A NO_ENTREZ — GPL6244 RESOLUTION AUDIT\n")
cat("============================================================\n")

print(
  no_entrez_audit
)


# ============================================================
# 12. Stage-2A MULTI_ENTREZ
# ============================================================

multi_entrez_audit <- comparison[
  Stage2AStatus ==
    "MULTI_ENTREZ_AMBIGUOUS",
  .N,
  by = .(
    N_ChipEntrez,
    N_PlatformEntrez,
    EntrezConcordance,
    SymbolConcordance
  )
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("STAGE-2A MULTI_ENTREZ — GPL6244 CROSS-CHECK\n")
cat("============================================================\n")

print(
  multi_entrez_audit
)


# ============================================================
# 13. Important potentially rescuable subsets
# ============================================================

rescue_summary <- data.table(
  Metric = c(
    "Stage2A_NO_ENTREZ_total",
    "NO_ENTREZ_platform_unique_symbol",
    "NO_ENTREZ_platform_unique_Entrez",
    "NO_ENTREZ_platform_multi_symbol",
    "Stage2A_MULTI_ENTREZ_total",
    "MULTI_ENTREZ_platform_unique_symbol",
    "MULTI_ENTREZ_platform_unique_Entrez",
    "MULTI_ENTREZ_exact_Entrez_set"
  ),

  N = c(
    comparison[
      Stage2AStatus ==
        "NO_ENTREZ_MAPPING",
      .N
    ],

    comparison[
      Stage2AStatus ==
        "NO_ENTREZ_MAPPING" &
        N_PlatformSymbols ==
          1L,
      .N
    ],

    comparison[
      Stage2AStatus ==
        "NO_ENTREZ_MAPPING" &
        N_PlatformEntrez ==
          1L,
      .N
    ],

    comparison[
      Stage2AStatus ==
        "NO_ENTREZ_MAPPING" &
        N_PlatformSymbols >
          1L,
      .N
    ],

    comparison[
      Stage2AStatus ==
        "MULTI_ENTREZ_AMBIGUOUS",
      .N
    ],

    comparison[
      Stage2AStatus ==
        "MULTI_ENTREZ_AMBIGUOUS" &
        N_PlatformSymbols ==
          1L,
      .N
    ],

    comparison[
      Stage2AStatus ==
        "MULTI_ENTREZ_AMBIGUOUS" &
        N_PlatformEntrez ==
          1L,
      .N
    ],

    comparison[
      Stage2AStatus ==
        "MULTI_ENTREZ_AMBIGUOUS" &
        EntrezConcordance ==
          "EXACT_SET",
      .N
    ]
  )
)


cat("\n============================================================\n")
cat("POTENTIALLY RESCUABLE / DISAMBIGUATING SUBSETS\n")
cat("============================================================\n")

print(
  rescue_summary
)


# ============================================================
# 14. Save outputs
# ============================================================

fwrite(
  category_summary,
  file.path(
    out_dir,
    "GSE73661_stage2b_GPL6244_category_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  platform_parsed,
  file.path(
    out_dir,
    "GSE73661_stage2b_GPL6244_gene_assignment_parsed_33252.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  parser_summary,
  file.path(
    out_dir,
    "GSE73661_stage2b_gene_assignment_parser_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  comparison,
  file.path(
    out_dir,
    "GSE73661_stage2b_ChipDb_vs_GPL6244_comparison.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  entrez_concordance,
  file.path(
    out_dir,
    "GSE73661_stage2b_entrez_concordance_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  symbol_concordance,
  file.path(
    out_dir,
    "GSE73661_stage2b_symbol_concordance_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  confident_audit,
  file.path(
    out_dir,
    "GSE73661_stage2b_confident_mapping_crosscheck.tsv"
  ),
  sep = "\t"
)


fwrite(
  no_entrez_audit,
  file.path(
    out_dir,
    "GSE73661_stage2b_no_entrez_resolution_audit.tsv"
  ),
  sep = "\t"
)


fwrite(
  multi_entrez_audit,
  file.path(
    out_dir,
    "GSE73661_stage2b_multi_entrez_crosscheck.tsv"
  ),
  sep = "\t"
)


fwrite(
  rescue_summary,
  file.path(
    out_dir,
    "GSE73661_stage2b_rescue_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 15. Final summary
# ============================================================

summary <- data.table(
  Metric = c(
    "GPL6244_total_platform_rows",
    "Expression_source_rows",
    "GPL6244_unused_rows",
    "Source_IDs_found_in_GPL6244",
    "Platform_no_symbol",
    "Platform_unique_symbol",
    "Platform_multi_symbol",
    "Platform_no_Entrez",
    "Platform_unique_Entrez",
    "Platform_multi_Entrez",
    "Rows_with_malformed_assignments",
    "Stage2A_confident_unique",
    "Stage2A_no_Entrez",
    "Stage2A_multi_Entrez",
    "Entrez_exact_set",
    "Entrez_different_set",
    "Entrez_platform_only",
    "Entrez_chipdb_only",
    "Entrez_both_none"
  ),

  N = c(
    nrow(platform),
    length(source_ids),
    length(extra_platform_ids),
    sum(source_present),

    sum(
      platform_parsed$N_PlatformSymbols ==
        0L
    ),

    sum(
      platform_parsed$N_PlatformSymbols ==
        1L
    ),

    sum(
      platform_parsed$N_PlatformSymbols >
        1L
    ),

    sum(
      platform_parsed$N_PlatformEntrez ==
        0L
    ),

    sum(
      platform_parsed$N_PlatformEntrez ==
        1L
    ),

    sum(
      platform_parsed$N_PlatformEntrez >
        1L
    ),

    sum(
      platform_parsed$MalformedAssignments >
        0L
    ),

    comparison[
      Stage2AStatus ==
        "UNIQUE_ENTREZ_CURRENT_SYMBOL",
      .N
    ],

    comparison[
      Stage2AStatus ==
        "NO_ENTREZ_MAPPING",
      .N
    ],

    comparison[
      Stage2AStatus ==
        "MULTI_ENTREZ_AMBIGUOUS",
      .N
    ],

    comparison[
      EntrezConcordance ==
        "EXACT_SET",
      .N
    ],

    comparison[
      EntrezConcordance ==
        "DIFFERENT_SET",
      .N
    ],

    comparison[
      EntrezConcordance ==
        "PLATFORM_ONLY",
      .N
    ],

    comparison[
      EntrezConcordance ==
        "CHIPDB_ONLY",
      .N
    ],

    comparison[
      EntrezConcordance ==
        "BOTH_NONE",
      .N
    ]
  )
)


cat("\n============================================================\n")
cat("STAGE 2B SUMMARY\n")
cat("============================================================\n")

print(
  summary
)


fwrite(
  summary,
  file.path(
    out_dir,
    "GSE73661_stage2b_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 16. Final gate
# ============================================================

critical_pass <- all(
  c(
    nrow(platform) ==
      33297L,

    length(source_ids) ==
      33252L,

    sum(source_present) ==
      33252L,

    length(extra_platform_ids) ==
      45L,

    nrow(platform_parsed) ==
      33252L,

    identical(
      platform_parsed$SourceProbeID,
      source_ids
    ),

    nrow(comparison) ==
      33252L,

    identical(
      comparison$SourceProbeID,
      source_ids
    ),

    identical(
      dim(expr),
      c(
        33252L,
        178L
      )
    )
  )
)


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] GPL6244 official platform table parsed\n")
  cat("[PASS] 33,297 GPL6244 platform rows verified\n")
  cat("[PASS] All 33,252 expression IDs aligned to GPL6244\n")
  cat("[PASS] 45 unused GPL6244 features identified\n")
  cat("[PASS] Multipart gene_assignment parsed\n")
  cat("[PASS] Platform gene symbols audited\n")
  cat("[PASS] Platform Entrez IDs audited\n")
  cat("[PASS] ChipDb vs GPL6244 mapping cross-check completed\n")
  cat("[PASS] Stage-2A NO_ENTREZ cases independently audited\n")
  cat("[PASS] Stage-2A MULTI_ENTREZ cases independently audited\n")
  cat("[PASS] No expression values modified\n")
  cat("[PASS] No probe collapse performed\n")
  cat("[INFO] Final gene-mapping decisions remain unfrozen\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE73661_STAGE2B_GPL6244_GENE_ASSIGNMENT_AUDIT\n")

} else {

  cat("[REVIEW] Stage 2B requires review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE73661_STAGE2B_GPL6244_GENE_ASSIGNMENT_AUDIT\n")
}

cat("============================================================\n")
