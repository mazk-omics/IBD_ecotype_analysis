suppressPackageStartupMessages({
  library(data.table)
})

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

harm_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_gene_id_harmonization.tsv.gz"
)

query_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_unresolved_numeric_single_query.tsv"
)

out_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

cat("============================================================\n")
cat("GSE57945 NCBI GeneID rescue classification\n")
cat("============================================================\n\n")

harm <- fread(
  harm_file,
  colClasses = "character"
)

q <- fread(
  query_file,
  colClasses = "character",
  na.strings = c("NA", "")
)

stopifnot(
  uniqueN(q$QueryGeneID) == nrow(q)
)

# ------------------------------------------------------------
# Classification
# ------------------------------------------------------------

q[
  ,
  RescueClass := fifelse(
    ReturnedGeneID == "QUERY_FAILED",
    "QUERY_FAILED",

    fifelse(
      !is.na(ReturnedGeneID) &
      ReturnedGeneID != QueryGeneID &
      TaxID == "9606" &
      !is.na(Symbol),
      "REDIRECTED_TO_CURRENT",

      fifelse(
        ReturnedGeneID == QueryGeneID &
        TaxID == "9606" &
        !is.na(Symbol),
        "CURRENT_IN_NCBI",

        fifelse(
          ReturnedGeneID == QueryGeneID &
          is.na(Symbol) &
          is.na(TaxID),
          "NO_LIVE_METADATA_FROM_DATASETS",
          "OTHER"
        )
      )
    )
  )
]

cat("Rescue classes:\n")
print(
  q[, .N, by = RescueClass][order(-N)]
)


# ------------------------------------------------------------
# Does target GeneID already exist in original matrix?
# ------------------------------------------------------------

original_ids <- harm$GeneID_raw

q[
  ,
  TargetAlreadyInMatrix :=
    RescueClass == "REDIRECTED_TO_CURRENT" &
    ReturnedGeneID %chin% original_ids
]


# ------------------------------------------------------------
# Attach target annotation already present in harmonization
# ------------------------------------------------------------

idx_target <- match(
  q$ReturnedGeneID,
  harm$GeneID_raw
)

q[
  ,
  ExistingTargetRawSymbol :=
    harm$GeneSymbol_raw[idx_target]
]

q[
  ,
  ExistingTargetCanonicalSymbol :=
    harm$CanonicalSymbol[idx_target]
]

q[
  ,
  ExistingTargetMappingStatus :=
    harm$MappingStatus[idx_target]
]


# ------------------------------------------------------------
# Number of old IDs mapping into same current target
# ------------------------------------------------------------

q[
  RescueClass == "REDIRECTED_TO_CURRENT",
  RedirectMultiplicity :=
    .N,
  by = ReturnedGeneID
]


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fwrite(
  q,
  file.path(
    out_dir,
    "GSE57945_ncbi_geneid_rescue_classified.tsv"
  ),
  sep = "\t"
)

redirect <- q[
  RescueClass == "REDIRECTED_TO_CURRENT"
][
  order(
    ReturnedGeneID,
    QueryGeneID
  )
]

fwrite(
  redirect,
  file.path(
    out_dir,
    "GSE57945_redirected_geneids.tsv"
  ),
  sep = "\t"
)

collision <- redirect[
  TargetAlreadyInMatrix == TRUE |
  RedirectMultiplicity > 1
]

fwrite(
  collision,
  file.path(
    out_dir,
    "GSE57945_geneid_rescue_collisions.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("REDIRECTION SUMMARY\n")
cat("============================================================\n")

cat(
  "Redirected old GeneIDs          :",
  nrow(redirect),
  "\n"
)

cat(
  "Unique current targets          :",
  uniqueN(redirect$ReturnedGeneID),
  "\n"
)

cat(
  "Target already exists in matrix :",
  sum(
    redirect$TargetAlreadyInMatrix,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Targets receiving >1 old ID     :",
  uniqueN(
    redirect[
      RedirectMultiplicity > 1,
      ReturnedGeneID
    ]
  ),
  "\n"
)

cat("\nExamples of collisions:\n")

print(
  head(
    collision[
      ,
      .(
        QueryGeneID,
        ReturnedGeneID,
        Symbol,
        TargetAlreadyInMatrix,
        ExistingTargetRawSymbol,
        ExistingTargetCanonicalSymbol,
        RedirectMultiplicity
      )
    ],
    30
  )
)

cat("\nIMPORTANT:\n")
cat(
  "- No expression matrix was modified.\n",
  "- No redirected GeneID was automatically collapsed.\n",
  "- Collisions must be evaluated before generating canonical expression matrix.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_NCBI_GENEID_RESCUE_CLASSIFIED\n")
