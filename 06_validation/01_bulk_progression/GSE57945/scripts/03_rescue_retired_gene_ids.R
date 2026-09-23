suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE57945
# Stage 2B: retired / non-current GeneID forensic rescue
#
# Uses:
#   NCBI current Homo_sapiens.gene_info
#   NCBI gene_history
#
# IMPORTANT:
#   - expression matrix is NOT modified
#   - no genes are dropped
#   - raw symbols are NOT used as authoritative rescue
# ============================================================

gse <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

ref <- "/home/mazekai/IBD_EcoTyper/06_validation/99_common/gene_annotation/NCBI_Gene_2026-09-21"

harm_file <- file.path(
  gse,
  "prepared",
  "02_gene_harmonization",
  "GSE57945_gene_id_harmonization.tsv.gz"
)

gene_info_file <- file.path(
  ref,
  "Homo_sapiens.gene_info.gz"
)

gene_history_file <- file.path(
  ref,
  "gene_history.gz"
)

out_dir <- file.path(
  gse,
  "prepared",
  "02_gene_harmonization"
)

cat("============================================================\n")
cat("GSE57945 retired GeneID forensic rescue\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load previous harmonization
# ============================================================

harm <- fread(
  harm_file,
  colClasses = "character"
)

cat("[1] Previous harmonization\n")
cat("Rows:", nrow(harm), "\n")

unresolved <- harm[
  MappingStatus == "ENTREZ_NOT_CURRENT"
]

cat("ENTREZ_NOT_CURRENT:", nrow(unresolved), "\n\n")


# ============================================================
# 2. Numeric / non-numeric split
# ============================================================

unresolved[
  ,
  NumericGeneID :=
    grepl("^[0-9]+$", GeneID_raw)
]

numeric_missing <- unresolved[
  NumericGeneID == TRUE
]

nonnumeric_missing <- unresolved[
  NumericGeneID == FALSE
]

cat("[2] Unresolved ID classes\n")
cat(
  "Numeric       :",
  nrow(numeric_missing),
  "\n"
)

cat(
  "Non-numeric   :",
  nrow(nonnumeric_missing),
  "\n\n"
)


# ============================================================
# 3. Load current human gene_info
# ============================================================

cat("[3] Loading current human gene_info...\n")

gi <- fread(
  gene_info_file,
  colClasses = "character"
)

setnames(
  gi,
  sub("^#", "", names(gi))
)

if (!all(
  c(
    "tax_id",
    "GeneID",
    "Symbol",
    "Synonyms",
    "description"
  ) %in% names(gi)
)) {
  stop("Unexpected gene_info schema.")
}

gi <- gi[
  tax_id == "9606"
]

cat(
  "Current human GeneIDs:",
  uniqueN(gi$GeneID),
  "\n\n"
)


# Preferred current symbol:
# use nomenclature authority symbol when available,
# otherwise NCBI Symbol.

official_symbol_col <-
  intersect(
    c(
      "Symbol_from_nomenclature_authority",
      "Symbol_from_nomenclature"
    ),
    names(gi)
  )

if (length(official_symbol_col) >= 1) {

  osc <- official_symbol_col[1]

  gi[
    ,
    CurrentSymbol :=
      fifelse(
        !is.na(get(osc)) &
        get(osc) != "" &
        get(osc) != "-",
        get(osc),
        Symbol
      )
  ]

} else {

  gi[
    ,
    CurrentSymbol := Symbol
  ]
}


# ============================================================
# 4. Load NCBI gene_history
# ============================================================

cat("[4] Loading gene_history...\n")

gh <- fread(
  gene_history_file,
  colClasses = "character"
)

setnames(
  gh,
  sub("^#", "", names(gh))
)

required_history <- c(
  "tax_id",
  "GeneID",
  "Discontinued_GeneID",
  "Discontinued_Symbol",
  "Discontinue_Date"
)

if (!all(
  required_history %in% names(gh)
)) {
  stop(
    "Unexpected gene_history schema.\nObserved columns:\n",
    paste(names(gh), collapse = ", ")
  )
}

gh_human <- gh[
  tax_id == "9606"
]

cat(
  "Human gene-history records:",
  nrow(gh_human),
  "\n\n"
)


# ============================================================
# 5. Build discontinued -> replacement map
# ============================================================

history_map <- gh_human[
  ,
  .(
    Discontinued_GeneID,
    Replacement_GeneID = GeneID,
    Discontinued_Symbol,
    Discontinue_Date
  )
]

# GeneID == "-" means discontinued without replacement.

history_map[
  ,
  HasReplacement :=
    !is.na(Replacement_GeneID) &
    Replacement_GeneID != "" &
    Replacement_GeneID != "-"
]


# ============================================================
# 6. Resolve replacement chains
# ============================================================

current_ids <- unique(
  gi$GeneID
)

replacement_lookup <- setNames(
  history_map$Replacement_GeneID,
  history_map$Discontinued_GeneID
)


resolve_geneid <- function(id, max_steps = 20L) {

  if (is.na(id) || id == "") {
    return(
      list(
        final_id = NA_character_,
        status = "INVALID_ID",
        steps = 0L,
        chain = NA_character_
      )
    )
  }

  original <- id
  chain <- id
  steps <- 0L
  seen <- character()

  repeat {

    if (id %chin% current_ids) {

      return(
        list(
          final_id = id,
          status =
            if (steps == 0L)
              "CURRENT"
            else
              "HISTORY_REPLACED_CURRENT",
          steps = steps,
          chain = paste(
            chain,
            collapse = ">"
          )
        )
      )
    }

    if (id %chin% seen) {

      return(
        list(
          final_id = NA_character_,
          status = "HISTORY_CYCLE",
          steps = steps,
          chain = paste(
            chain,
            collapse = ">"
          )
        )
      )
    }

    seen <- c(seen, id)

    next_id <- replacement_lookup[[id]]

    if (
      is.null(next_id) ||
      is.na(next_id) ||
      next_id == ""
    ) {

      return(
        list(
          final_id = NA_character_,
          status = "NOT_FOUND_IN_HISTORY",
          steps = steps,
          chain = paste(
            chain,
            collapse = ">"
          )
        )
      )
    }

    if (next_id == "-") {

      return(
        list(
          final_id = NA_character_,
          status = "DISCONTINUED_NO_REPLACEMENT",
          steps = steps + 1L,
          chain = paste(
            c(chain, "-"),
            collapse = ">"
          )
        )
      )
    }

    id <- next_id
    steps <- steps + 1L
    chain <- c(
      chain,
      id
    )

    if (steps >= max_steps) {

      return(
        list(
          final_id = NA_character_,
          status = "MAX_HISTORY_DEPTH",
          steps = steps,
          chain = paste(
            chain,
            collapse = ">"
          )
        )
      )
    }
  }
}


cat("[5] Resolving numeric discontinued GeneIDs...\n")

resolved_list <- lapply(
  numeric_missing$GeneID_raw,
  resolve_geneid
)

resolved_numeric <- data.table(
  GeneID_raw =
    numeric_missing$GeneID_raw,

  Rescue_GeneID =
    vapply(
      resolved_list,
      `[[`,
      character(1),
      "final_id"
    ),

  RescueStatus =
    vapply(
      resolved_list,
      `[[`,
      character(1),
      "status"
    ),

  ReplacementSteps =
    vapply(
      resolved_list,
      `[[`,
      integer(1),
      "steps"
    ),

  ReplacementChain =
    vapply(
      resolved_list,
      `[[`,
      character(1),
      "chain"
    )
)


# ============================================================
# 7. Attach current NCBI annotation to rescued IDs
# ============================================================

idx <- match(
  resolved_numeric$Rescue_GeneID,
  gi$GeneID
)

resolved_numeric[
  ,
  RescueSymbol :=
    gi$CurrentSymbol[idx]
]

resolved_numeric[
  ,
  RescueNCBISymbol :=
    gi$Symbol[idx]
]

resolved_numeric[
  ,
  RescueDescription :=
    gi$description[idx]
]


# Attach original raw symbol
idx_raw <- match(
  resolved_numeric$GeneID_raw,
  harm$GeneID_raw
)

resolved_numeric[
  ,
  GeneSymbol_raw :=
    harm$GeneSymbol_raw[idx_raw]
]


# ============================================================
# 8. Attach original gene_history record
# ============================================================

first_hist <- history_map[
  match(
    resolved_numeric$GeneID_raw,
    history_map$Discontinued_GeneID
  )
]

resolved_numeric[
  ,
  HistoryFirstReplacement :=
    first_hist$Replacement_GeneID
]

resolved_numeric[
  ,
  HistoryDiscontinuedSymbol :=
    first_hist$Discontinued_Symbol
]

resolved_numeric[
  ,
  HistoryDiscontinueDate :=
    first_hist$Discontinue_Date
]


# ============================================================
# 9. Non-numeric IDs
# ============================================================

nonnumeric_out <- nonnumeric_missing[
  ,
  .(
    GeneID_raw,
    GeneSymbol_raw
  )
]

nonnumeric_out[
  ,
  RescueStatus :=
    "NONNUMERIC_GENE_ID"
]


# ============================================================
# 10. Summary
# ============================================================

status_summary <- resolved_numeric[
  ,
  .N,
  by = RescueStatus
][
  order(-N)
]

cat("\n")
cat("============================================================\n")
cat("NUMERIC GeneID RESCUE STATUS\n")
cat("============================================================\n")

print(status_summary)


n_rescued <- resolved_numeric[
  RescueStatus ==
    "HISTORY_REPLACED_CURRENT",
  .N
]

n_no_replacement <- resolved_numeric[
  RescueStatus ==
    "DISCONTINUED_NO_REPLACEMENT",
  .N
]

n_not_found <- resolved_numeric[
  RescueStatus ==
    "NOT_FOUND_IN_HISTORY",
  .N
]


summary_dt <- data.table(
  Metric = c(
    "Total_not_current",
    "Numeric_not_current",
    "NonNumeric_not_current",
    "Numeric_rescued_to_current",
    "Numeric_discontinued_no_replacement",
    "Numeric_not_found_in_history"
  ),

  Value = c(
    nrow(unresolved),
    nrow(numeric_missing),
    nrow(nonnumeric_missing),
    n_rescued,
    n_no_replacement,
    n_not_found
  )
)


# ============================================================
# 11. Save
# ============================================================

fwrite(
  resolved_numeric,
  file.path(
    out_dir,
    "GSE57945_retired_numeric_entrez_rescue.tsv"
  ),
  sep = "\t"
)

fwrite(
  nonnumeric_out,
  file.path(
    out_dir,
    "GSE57945_nonnumeric_gene_ids.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary_dt,
  file.path(
    out_dir,
    "GSE57945_retired_gene_id_rescue_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 12. Examples
# ============================================================

cat("\nRescued examples:\n")

print(
  head(
    resolved_numeric[
      RescueStatus ==
        "HISTORY_REPLACED_CURRENT",
      .(
        GeneID_raw,
        GeneSymbol_raw,
        Rescue_GeneID,
        RescueSymbol,
        ReplacementChain,
        HistoryDiscontinueDate
      )
    ],
    20
  )
)

cat("\nNon-numeric examples:\n")

print(
  head(
    nonnumeric_out,
    30
  )
)


cat("\n")
cat("============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n")

print(summary_dt)

cat("\nIMPORTANT:\n")
cat(
  "- No expression values were changed.\n",
  "- Gene-history replacement is treated separately from symbol-based rescue.\n",
  "- Non-numeric and no-replacement entries remain unresolved for now.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_RETIRED_GENE_ID_FORENSIC_AUDIT_COMPLETED\n")
cat("============================================================\n")
