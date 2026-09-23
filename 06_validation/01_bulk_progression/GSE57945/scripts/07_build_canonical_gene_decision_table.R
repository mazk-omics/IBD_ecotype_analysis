suppressPackageStartupMessages({
  library(data.table)
})

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

harm_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_gene_id_harmonization.tsv.gz"
)

rescue_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_ncbi_geneid_rescue_classified.tsv"
)

out_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

cat("============================================================\n")
cat("GSE57945 canonical gene decision table\n")
cat("============================================================\n\n")

harm <- fread(
  harm_file,
  colClasses = "character"
)

rescue <- fread(
  rescue_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(nrow(harm) == 36372)

# ------------------------------------------------------------
# Base decision table
# ------------------------------------------------------------

d <- copy(harm)

d[
  ,
  `:=`(
    FinalGeneID = NA_character_,
    FinalSymbol = NA_character_,
    CanonicalAction = NA_character_,
    CanonicalReason = NA_character_
  )
]


# ------------------------------------------------------------
# 1. Direct current mappings
# ------------------------------------------------------------

direct_idx <- which(
  d$MappingStatus != "ENTREZ_NOT_CURRENT" &
  grepl("^[0-9]+$", d$GeneID_raw) &
  !is.na(d$CanonicalSymbol) &
  d$CanonicalSymbol != ""
)

d[
  direct_idx,
  `:=`(
    FinalGeneID = GeneID_raw,
    FinalSymbol = CanonicalSymbol,
    CanonicalAction = "KEEP_DIRECT_CURRENT",
    CanonicalReason =
      "Direct current GeneID with unique canonical symbol"
  )
]


# ------------------------------------------------------------
# 2. Numeric IDs absent from org.Hs.eg.db
# ------------------------------------------------------------

r <- rescue

# redirect multiplicity
r[
  !is.na(ReturnedGeneID),
  RedirectMultiplicityAll :=
    .N,
  by = ReturnedGeneID
]

# current targets already represented by clean direct rows
direct_geneids <- d[
  CanonicalAction == "KEEP_DIRECT_CURRENT",
  GeneID_raw
]

direct_symbols <- d[
  CanonicalAction == "KEEP_DIRECT_CURRENT",
  FinalSymbol
]


# ------------------------------------------------------------
# Candidate safe redirects
# ------------------------------------------------------------

r[
  ,
  SafeRedirectCandidate :=
    RescueClass == "REDIRECTED_TO_CURRENT" &
    !is.na(ReturnedGeneID) &
    !ReturnedGeneID %chin% d$GeneID_raw &
    RedirectMultiplicityAll == 1 &
    !is.na(Symbol) &
    Symbol != ""
]


# symbol collision with clean direct core
r[
  ,
  SymbolCollisionWithDirect :=
    SafeRedirectCandidate &
    Symbol %chin% direct_symbols
]

r[
  ,
  SafeRedirect :=
    SafeRedirectCandidate &
    !SymbolCollisionWithDirect
]


# Apply numeric decisions
for (i in seq_len(nrow(r))) {

  old <- r$QueryGeneID[i]

  idx <- which(d$GeneID_raw == old)

  if (length(idx) != 1)
    stop("Cannot uniquely locate GeneID: ", old)

  if (isTRUE(r$SafeRedirect[i])) {

    d[
      idx,
      `:=`(
        FinalGeneID = r$ReturnedGeneID[i],
        FinalSymbol = r$Symbol[i],
        CanonicalAction = "KEEP_SAFE_REDIRECT",
        CanonicalReason =
          "Singleton old-to-current redirect; target absent; no symbol collision"
      )
    ]

  } else if (
    r$RescueClass[i] == "REDIRECTED_TO_CURRENT" &&
    r$ReturnedGeneID[i] %chin% d$GeneID_raw
  ) {

    d[
      idx,
      `:=`(
        CanonicalAction = "EXCLUDE_REDIRECT_TARGET_ALREADY_PRESENT",
        CanonicalReason =
          "Current target already represented by another original matrix row"
      )
    ]

  } else if (
    r$RescueClass[i] == "REDIRECTED_TO_CURRENT" &&
    !is.na(r$RedirectMultiplicityAll[i]) &&
    r$RedirectMultiplicityAll[i] > 1
  ) {

    d[
      idx,
      `:=`(
        CanonicalAction = "EXCLUDE_REDIRECT_MANY_TO_ONE",
        CanonicalReason =
          "Multiple historical GeneIDs redirect to the same current GeneID"
      )
    ]

  } else if (
    r$RescueClass[i] == "NO_LIVE_METADATA_FROM_DATASETS"
  ) {

    d[
      idx,
      `:=`(
        CanonicalAction = "EXCLUDE_NO_LIVE_METADATA",
        CanonicalReason =
          "No live current gene metadata returned by NCBI Datasets"
      )
    ]

  } else {

    d[
      idx,
      `:=`(
        CanonicalAction = "EXCLUDE_UNRESOLVED_NUMERIC",
        CanonicalReason =
          "Numeric GeneID not safely resolvable to a unique live current gene"
      )
    ]
  }
}


# ------------------------------------------------------------
# 3. Locus-suffixed IDs
# ------------------------------------------------------------

locus_idx <- which(
  grepl(
    "^[0-9]+\\.chr([0-9]+|X|Y)(\\.[0-9]+)?$",
    d$GeneID_raw
  )
)

d[
  locus_idx,
  `:=`(
    CanonicalAction = "EXCLUDE_LOCUS_SUFFIXED",
    CanonicalReason =
      "Locus-specific feature; multiple rows per base GeneID and no justified gene-level collapse"
  )
]


# ------------------------------------------------------------
# Sanity
# ------------------------------------------------------------

if (any(is.na(d$CanonicalAction))) {

  cat("\nUnclassified rows:\n")

  print(
    d[
      is.na(CanonicalAction),
      .(
        GeneID_raw,
        GeneSymbol_raw,
        MappingStatus
      )
    ]
  )

  stop("Some rows remain unclassified.")
}


# ------------------------------------------------------------
# Included canonical set
# ------------------------------------------------------------

included <- d[
  CanonicalAction %chin% c(
    "KEEP_DIRECT_CURRENT",
    "KEEP_SAFE_REDIRECT"
  )
]

cat("Decision counts:\n")

print(
  d[
    ,
    .N,
    by = CanonicalAction
  ][order(-N)]
)

cat("\nIncluded rows:", nrow(included), "\n")
cat(
  "Unique FinalGeneID:",
  uniqueN(included$FinalGeneID),
  "\n"
)

cat(
  "Unique FinalSymbol:",
  uniqueN(included$FinalSymbol),
  "\n"
)


# ------------------------------------------------------------
# Critical uniqueness checks
# ------------------------------------------------------------

stopifnot(
  nrow(included) ==
    uniqueN(included$FinalGeneID)
)

stopifnot(
  nrow(included) ==
    uniqueN(included$FinalSymbol)
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fwrite(
  d,
  file.path(
    out_dir,
    "GSE57945_canonical_gene_decision_table.tsv"
  ),
  sep = "\t"
)

fwrite(
  included,
  file.path(
    out_dir,
    "GSE57945_canonical_gene_included_features.tsv"
  ),
  sep = "\t"
)


summary <- d[
  ,
  .N,
  by = CanonicalAction
][order(-N)]

fwrite(
  summary,
  file.path(
    out_dir,
    "GSE57945_canonical_gene_decision_summary.tsv"
  ),
  sep = "\t"
)


cat("\nFINAL STATUS:\n")
cat("PASS_CANONICAL_GENE_DECISION_TABLE_BUILT\n")
cat("============================================================\n")
