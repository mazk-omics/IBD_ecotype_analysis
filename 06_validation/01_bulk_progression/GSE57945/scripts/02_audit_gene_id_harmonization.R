suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ============================================================
# GSE57945
# Stage 2A: Gene ID identity + canonical symbol harmonization
#
# IMPORTANT:
#   - expression matrix is NOT changed
#   - genes are NOT filtered
#   - duplicate symbols are NOT collapsed
#   - raw annotation is preserved
#
# Purpose:
#   Determine whether GeneID_raw represents Entrez Gene IDs,
#   map valid current Entrez IDs to current SYMBOL annotation,
#   and audit all mapping failures / symbol changes / duplicates.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

input_file <- file.path(
  base,
  "prepared",
  "01_authoritative_matrix",
  "GSE57945_authoritative_expression_322.rds"
)

out_dir <- file.path(
  base,
  "prepared",
  "02_gene_harmonization"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE57945 Gene ID harmonization audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load authoritative object
# ============================================================

obj <- readRDS(input_file)

expr <- obj$expression
ann  <- as.data.table(
  copy(obj$gene_annotation)
)

stopifnot(
  nrow(expr) == 36372,
  ncol(expr) == 322,
  nrow(ann) == 36372
)

ann[, GeneID_raw := as.character(GeneID_raw)]
ann[, GeneSymbol_raw := as.character(GeneSymbol_raw)]

if (!identical(
  rownames(expr),
  ann$GeneID_raw
)) {
  stop(
    "Gene annotation and expression matrix are not aligned."
  )
}

cat("[PASS] authoritative matrix loaded\n")
cat("[PASS] 36,372 genes x 322 samples\n\n")


# ============================================================
# 2. Basic Gene ID identity audit
# ============================================================

ann[, GeneID_numeric_only :=
      grepl("^[0-9]+$", GeneID_raw)]

ann[, GeneID_has_whitespace :=
      grepl("^[[:space:]]|[[:space:]]$", GeneID_raw)]

ann[, GeneID_empty :=
      is.na(GeneID_raw) | GeneID_raw == ""]

cat("[1] Raw Gene ID structure\n")
cat(
  "Numeric-only IDs       :",
  sum(ann$GeneID_numeric_only),
  "/",
  nrow(ann),
  "\n"
)

cat(
  "Whitespace IDs         :",
  sum(ann$GeneID_has_whitespace),
  "\n"
)

cat(
  "Empty IDs              :",
  sum(ann$GeneID_empty),
  "\n"
)

cat(
  "Unique IDs             :",
  uniqueN(ann$GeneID_raw),
  "\n\n"
)


# ============================================================
# 3. Audit org.Hs.eg.db database
# ============================================================

cat("[2] Annotation database\n")

cat(
  "org.Hs.eg.db version   :",
  as.character(packageVersion("org.Hs.eg.db")),
  "\n"
)

cat(
  "AnnotationDbi version  :",
  as.character(packageVersion("AnnotationDbi")),
  "\n"
)

available_keytypes <- keytypes(org.Hs.eg.db)

if (!"ENTREZID" %in% available_keytypes) {
  stop("ENTREZID is not available in org.Hs.eg.db.")
}

db_entrez <- keys(
  org.Hs.eg.db,
  keytype = "ENTREZID"
)

cat(
  "Current ENTREZID keys  :",
  length(db_entrez),
  "\n\n"
)


# ============================================================
# 4. Determine current Entrez ID presence
# ============================================================

ann[
  ,
  CurrentEntrezPresent :=
    GeneID_raw %chin% db_entrez
]

cat("[3] Entrez ID database overlap\n")

cat(
  "Present in org.Hs.eg.db:",
  sum(ann$CurrentEntrezPresent),
  "\n"
)

cat(
  "Not present            :",
  sum(!ann$CurrentEntrezPresent),
  "\n"
)

cat(
  "Coverage               :",
  sprintf(
    "%.2f%%",
    100 * mean(ann$CurrentEntrezPresent)
  ),
  "\n\n"
)


# ============================================================
# 5. Map valid ENTREZIDs to SYMBOL
#    Keep ALL mappings; do not arbitrarily take first.
# ============================================================

valid_ids <- ann[
  CurrentEntrezPresent == TRUE,
  GeneID_raw
]

cat("[4] Mapping ENTREZID -> SYMBOL...\n")

symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = valid_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "CharacterList"
)

symbol_list <- as.list(symbol_map)

symbol_dt <- rbindlist(
  lapply(
    names(symbol_list),
    function(id) {

      vals <- as.character(
        symbol_list[[id]]
      )

      vals <- vals[
        !is.na(vals) &
        vals != ""
      ]

      vals <- sort(
        unique(vals)
      )

      data.table(
        GeneID_raw = id,
        CanonicalSymbol_n =
          length(vals),

        CanonicalSymbol_all =
          if (length(vals) == 0)
            NA_character_
          else
            paste(
              vals,
              collapse = "|"
            ),

        CanonicalSymbol =
          if (length(vals) == 1)
            vals
          else
            NA_character_
      )
    }
  )
)


# ============================================================
# 6. Optional GENENAME mapping
# ============================================================

cat("[5] Mapping ENTREZID -> GENENAME...\n")

genename_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = valid_ids,
  column = "GENENAME",
  keytype = "ENTREZID",
  multiVals = "CharacterList"
)

genename_list <- as.list(genename_map)

genename_dt <- rbindlist(
  lapply(
    names(genename_list),
    function(id) {

      vals <- as.character(
        genename_list[[id]]
      )

      vals <- vals[
        !is.na(vals) &
        vals != ""
      ]

      vals <- sort(
        unique(vals)
      )

      data.table(
        GeneID_raw = id,

        GeneName_n =
          length(vals),

        GeneName_current =
          if (length(vals) == 1)
            vals
          else if (length(vals) > 1)
            paste(
              vals,
              collapse = "|"
            )
          else
            NA_character_
      )
    }
  )
)


# ============================================================
# 7. Add mappings back in original gene order
# ============================================================

idx_symbol <- match(
  ann$GeneID_raw,
  symbol_dt$GeneID_raw
)

ann[
  ,
  CanonicalSymbol_n :=
    symbol_dt$CanonicalSymbol_n[
      idx_symbol
    ]
]

ann[
  ,
  CanonicalSymbol_all :=
    symbol_dt$CanonicalSymbol_all[
      idx_symbol
    ]
]

ann[
  ,
  CanonicalSymbol :=
    symbol_dt$CanonicalSymbol[
      idx_symbol
    ]
]

idx_name <- match(
  ann$GeneID_raw,
  genename_dt$GeneID_raw
)

ann[
  ,
  GeneName_n :=
    genename_dt$GeneName_n[
      idx_name
    ]
]

ann[
  ,
  GeneName_current :=
    genename_dt$GeneName_current[
      idx_name
    ]
]


# ============================================================
# 8. Mapping status
# ============================================================

ann[
  ,
  MappingStatus :=
    fifelse(
      !CurrentEntrezPresent,
      "ENTREZ_NOT_CURRENT",

      fifelse(
        is.na(CanonicalSymbol_n) |
        CanonicalSymbol_n == 0,
        "CURRENT_ENTREZ_NO_SYMBOL",

        fifelse(
          CanonicalSymbol_n == 1,
          "UNIQUE_SYMBOL",
          "AMBIGUOUS_SYMBOL"
        )
      )
    )
]


# ============================================================
# 9. Raw vs current symbol comparison
# ============================================================

ann[
  ,
  RawSymbolExactMatch :=
    MappingStatus == "UNIQUE_SYMBOL" &
    !is.na(GeneSymbol_raw) &
    GeneSymbol_raw == CanonicalSymbol
]

ann[
  ,
  RawSymbolChanged :=
    MappingStatus == "UNIQUE_SYMBOL" &
    !is.na(GeneSymbol_raw) &
    GeneSymbol_raw != CanonicalSymbol
]


# ============================================================
# 10. Detect obvious historical spreadsheet/date-like symbols
# ============================================================

date_patterns <- c(
  "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}$",
  "^[0-9]{1,2}-[A-Za-z]{3}$",
  "^[A-Za-z]{3}-[0-9]{1,2}$",
  "^[0-9]{1,2}/[A-Za-z]{3}$",
  "^[A-Za-z]{3}/[0-9]{1,2}$"
)

date_regex <- paste(
  date_patterns,
  collapse = "|"
)

ann[
  ,
  RawSymbolLooksDateLike :=
    grepl(
      date_regex,
      GeneSymbol_raw
    )
]


# ============================================================
# 11. Canonical symbol duplicate audit
# ============================================================

ann[
  ,
  CanonicalSymbolMultiplicity :=
    NA_integer_
]

ann[
  !is.na(CanonicalSymbol),
  CanonicalSymbolMultiplicity := .N,
  by = CanonicalSymbol
]

ann[
  ,
  CanonicalSymbolDuplicated :=
    !is.na(CanonicalSymbolMultiplicity) &
    CanonicalSymbolMultiplicity > 1
]


# ============================================================
# 12. Summary statistics
# ============================================================

n_total <- nrow(ann)

n_numeric <-
  sum(ann$GeneID_numeric_only)

n_current_entrez <-
  sum(ann$CurrentEntrezPresent)

n_not_current <-
  sum(!ann$CurrentEntrezPresent)

n_unique_symbol <-
  sum(
    ann$MappingStatus ==
      "UNIQUE_SYMBOL"
  )

n_ambiguous_symbol <-
  sum(
    ann$MappingStatus ==
      "AMBIGUOUS_SYMBOL"
  )

n_current_no_symbol <-
  sum(
    ann$MappingStatus ==
      "CURRENT_ENTREZ_NO_SYMBOL"
  )

n_exact <-
  sum(
    ann$RawSymbolExactMatch,
    na.rm = TRUE
  )

n_changed <-
  sum(
    ann$RawSymbolChanged,
    na.rm = TRUE
  )

n_date_like <-
  sum(
    ann$RawSymbolLooksDateLike,
    na.rm = TRUE
  )

n_date_like_rescued <-
  sum(
    ann$RawSymbolLooksDateLike &
    ann$MappingStatus == "UNIQUE_SYMBOL",
    na.rm = TRUE
  )

n_unique_canonical_symbols <-
  uniqueN(
    ann[
      !is.na(CanonicalSymbol),
      CanonicalSymbol
    ]
  )

n_duplicate_canonical_names <-
  uniqueN(
    ann[
      CanonicalSymbolDuplicated == TRUE,
      CanonicalSymbol
    ]
  )

n_rows_duplicate_canonical <-
  sum(
    ann$CanonicalSymbolDuplicated,
    na.rm = TRUE
  )


summary_dt <- data.table(
  Metric = c(
    "Total_gene_rows",
    "Numeric_only_GeneID",
    "Current_Entrez_present",
    "Current_Entrez_absent",
    "Unique_symbol_mapping",
    "Ambiguous_symbol_mapping",
    "Current_Entrez_without_symbol",
    "Raw_symbol_exact_current_symbol",
    "Raw_symbol_changed",
    "Raw_symbol_date_like",
    "Date_like_symbol_with_current_mapping",
    "Unique_current_symbols",
    "Duplicated_current_symbol_names",
    "Rows_in_duplicated_current_symbol_groups"
  ),

  Value = c(
    n_total,
    n_numeric,
    n_current_entrez,
    n_not_current,
    n_unique_symbol,
    n_ambiguous_symbol,
    n_current_no_symbol,
    n_exact,
    n_changed,
    n_date_like,
    n_date_like_rescued,
    n_unique_canonical_symbols,
    n_duplicate_canonical_names,
    n_rows_duplicate_canonical
  )
)


# ============================================================
# 13. Save full harmonization table
# ============================================================

full_file <- file.path(
  out_dir,
  "GSE57945_gene_id_harmonization.tsv.gz"
)

fwrite(
  ann,
  full_file,
  sep = "\t",
  compress = "gzip"
)

fwrite(
  summary_dt,
  file.path(
    out_dir,
    "GSE57945_gene_id_harmonization_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 14. Save forensic subsets
# ============================================================

fwrite(
  ann[
    CurrentEntrezPresent == FALSE
  ],
  file.path(
    out_dir,
    "GSE57945_entrez_not_current.tsv"
  ),
  sep = "\t"
)

fwrite(
  ann[
    MappingStatus ==
      "CURRENT_ENTREZ_NO_SYMBOL"
  ],
  file.path(
    out_dir,
    "GSE57945_current_entrez_without_symbol.tsv"
  ),
  sep = "\t"
)

fwrite(
  ann[
    MappingStatus ==
      "AMBIGUOUS_SYMBOL"
  ],
  file.path(
    out_dir,
    "GSE57945_ambiguous_symbol_mapping.tsv"
  ),
  sep = "\t"
)

fwrite(
  ann[
    RawSymbolChanged == TRUE,
    .(
      GeneID_raw,
      GeneSymbol_raw,
      CanonicalSymbol,
      GeneName_current,
      RawSymbolLooksDateLike
    )
  ],
  file.path(
    out_dir,
    "GSE57945_raw_vs_current_symbol_changes.tsv"
  ),
  sep = "\t"
)

fwrite(
  ann[
    RawSymbolLooksDateLike == TRUE,
    .(
      GeneID_raw,
      GeneSymbol_raw,
      MappingStatus,
      CanonicalSymbol,
      GeneName_current
    )
  ],
  file.path(
    out_dir,
    "GSE57945_date_like_raw_symbols.tsv"
  ),
  sep = "\t"
)

fwrite(
  ann[
    CanonicalSymbolDuplicated == TRUE
  ][
    order(
      CanonicalSymbol,
      GeneID_raw
    )
  ],
  file.path(
    out_dir,
    "GSE57945_duplicated_canonical_symbols.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 15. Known-ID sanity checks
# ============================================================

sanity_ids <- c(
  "366",   # AQP9
  "3576",  # historic IL8 / current CXCL8
  "4314"   # MMP3
)

sanity <- ann[
  GeneID_raw %chin% sanity_ids,
  .(
    GeneID_raw,
    GeneSymbol_raw,
    CanonicalSymbol,
    GeneName_current,
    MappingStatus
  )
]

fwrite(
  sanity,
  file.path(
    out_dir,
    "GSE57945_known_gene_mapping_sanity.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 16. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GENE HARMONIZATION SUMMARY\n")
cat("============================================================\n")

print(summary_dt)

cat("\nKnown-gene sanity check:\n")
print(sanity)

cat("\nMapping status distribution:\n")
print(
  ann[
    ,
    .N,
    by = MappingStatus
  ][
    order(-N)
  ]
)

cat("\n")
cat("IMPORTANT:\n")
cat(
  "- Expression matrix was NOT modified.\n",
  "- No genes were removed.\n",
  "- No duplicate symbols were collapsed.\n",
  "- CanonicalSymbol is assigned only for unique SYMBOL mappings.\n",
  sep = ""
)

cat("\n")
cat("FINAL STATUS:\n")
cat("PASS_GENE_ID_HARMONIZATION_AUDIT_COMPLETED\n")
cat("============================================================\n")
