suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

mapping_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_EcoTyper_NMF_gene_identity_mapping.tsv"
)

decision_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_canonical_gene_decision_table.tsv"
)

out_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

cat("============================================================\n")
cat("EcoTyper unresolved symbol alias audit\n")
cat("============================================================\n\n")

m <- fread(
  mapping_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

d <- fread(
  decision_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

included <- d[
  CanonicalAction %chin% c(
    "KEEP_DIRECT_CURRENT",
    "KEEP_SAFE_REDIRECT"
  )
]

missing <- unique(
  m[
    MatchClass == "ABSENT_FROM_GSE57945",
    EcoTyperGene
  ]
)

cat("Input unresolved genes:", length(missing), "\n\n")

# ------------------------------------------------------------
# Alias lookup
# ------------------------------------------------------------

alias_map <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = missing,
    keytype = "ALIAS",
    columns = c(
      "ENTREZID",
      "SYMBOL"
    )
  )
)

alias_map <- as.data.table(alias_map)

setnames(
  alias_map,
  "ALIAS",
  "EcoTyperGene"
)

# remove failed mappings
alias_map <- alias_map[
  !is.na(ENTREZID)
]

# Is mapped Entrez ID actually represented in canonical GSE57945?
alias_map[
  ,
  InCanonical :=
    ENTREZID %chin% included$FinalGeneID
]

# attach canonical target
canon_lookup <- included[
  ,
  .(
    ENTREZID = FinalGeneID,
    FinalSymbol,
    SourceGeneID = GeneID_raw,
    SourceRawSymbol = GeneSymbol_raw
  )
]

alias_map <- merge(
  alias_map,
  canon_lookup,
  by = "ENTREZID",
  all.x = TRUE
)

# ------------------------------------------------------------
# Classification per EcoTyper gene
# ------------------------------------------------------------

result_list <- vector(
  "list",
  length(missing)
)

for (i in seq_along(missing)) {

  g <- missing[i]

  x <- alias_map[
    EcoTyperGene == g &
    InCanonical == TRUE
  ]

  targets <- unique(
    x[
      !is.na(ENTREZID),
      ENTREZID
    ]
  )

  if (length(targets) == 1) {

    z <- x[
      ENTREZID == targets[1]
    ][1]

    result_list[[i]] <- data.table(
      EcoTyperGene = g,
      AliasClass = "SAFE_ALIAS_TO_CANONICAL",
      TargetGeneID = z$ENTREZID,
      TargetCanonicalSymbol = z$FinalSymbol,
      OrgHsSymbol = z$SYMBOL,
      SourceRawSymbol = z$SourceRawSymbol
    )

  } else if (length(targets) > 1) {

    result_list[[i]] <- data.table(
      EcoTyperGene = g,
      AliasClass = "AMBIGUOUS_ALIAS",
      TargetGeneID = paste(
        sort(targets),
        collapse = "|"
      ),
      TargetCanonicalSymbol = paste(
        sort(unique(x$FinalSymbol)),
        collapse = "|"
      ),
      OrgHsSymbol = paste(
        sort(unique(x$SYMBOL)),
        collapse = "|"
      ),
      SourceRawSymbol = paste(
        sort(unique(x$SourceRawSymbol)),
        collapse = "|"
      )
    )

  } else {

    any_alias <- alias_map[
      EcoTyperGene == g
    ]

    if (nrow(any_alias) > 0) {

      result_list[[i]] <- data.table(
        EcoTyperGene = g,
        AliasClass = "ALIAS_NOT_IN_CANONICAL",
        TargetGeneID = paste(
          sort(unique(any_alias$ENTREZID)),
          collapse = "|"
        ),
        TargetCanonicalSymbol = NA_character_,
        OrgHsSymbol = paste(
          sort(unique(any_alias$SYMBOL)),
          collapse = "|"
        ),
        SourceRawSymbol = NA_character_
      )

    } else {

      result_list[[i]] <- data.table(
        EcoTyperGene = g,
        AliasClass = "NO_ALIAS_MAPPING",
        TargetGeneID = NA_character_,
        TargetCanonicalSymbol = NA_character_,
        OrgHsSymbol = NA_character_,
        SourceRawSymbol = NA_character_
      )
    }
  }
}

result <- rbindlist(
  result_list,
  fill = TRUE
)

cat("============================================================\n")
cat("ALIAS CLASS SUMMARY\n")
cat("============================================================\n")

print(
  result[
    ,
    .N,
    by = AliasClass
  ][order(-N)]
)

cat("\n============================================================\n")
cat("DETAIL\n")
cat("============================================================\n")

print(result, nrows = 100)

safe <- result[
  AliasClass == "SAFE_ALIAS_TO_CANONICAL"
]

cat("\nSafe alias rescues:", nrow(safe), "\n")

cat(
  "Remaining unresolved:",
  nrow(result) - nrow(safe),
  "\n"
)

# check target uniqueness
dup <- safe[
  ,
  .N,
  by = TargetGeneID
][N > 1]

cat(
  "Alias targets receiving >1 EcoTyper gene:",
  nrow(dup),
  "\n"
)

fwrite(
  result,
  file.path(
    out_dir,
    "GSE57945_EcoTyper_missing_alias_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  safe,
  file.path(
    out_dir,
    "GSE57945_EcoTyper_safe_alias_rescues.tsv"
  ),
  sep = "\t"
)

fwrite(
  dup,
  file.path(
    out_dir,
    "GSE57945_EcoTyper_alias_duplicate_targets.tsv"
  ),
  sep = "\t"
)

cat("\nFINAL STATUS:\n")
cat("PASS_ECOTYPER_MISSING_ALIAS_AUDIT\n")
cat("============================================================\n")
