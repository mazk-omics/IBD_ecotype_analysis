suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

obj_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE57945_authoritative_expression_322.rds"
)

harm_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_gene_id_harmonization.tsv.gz"
)

out_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

cat("============================================================\n")
cat("GSE57945 locus-suffixed GeneID audit\n")
cat("============================================================\n\n")

obj <- readRDS(obj_file)
expr <- obj$expression

harm <- fread(
  harm_file,
  colClasses = "character"
)

# ============================================================
# 1. Identify non-numeric IDs
# ============================================================

nn <- harm[
  !grepl("^[0-9]+$", GeneID_raw)
]

cat("Non-numeric rows:", nrow(nn), "\n")

# Expected form:
# 12345.chr1.1
# 12345.chrX
# 12345.chrY
# etc.

nn[
  ,
  BaseGeneID :=
    sub("\\.chr.*$", "", GeneID_raw)
]

nn[
  ,
  LocusSuffix :=
    sub("^[0-9]+\\.", "", GeneID_raw)
]

nn[
  ,
  ParsedSuccessfully :=
    grepl(
      "^[0-9]+\\.chr([0-9]+|X|Y)(\\.[0-9]+)?$",
      GeneID_raw
    )
]

cat(
  "Successfully parsed:",
  sum(nn$ParsedSuccessfully),
  "/",
  nrow(nn),
  "\n\n"
)

if (!all(nn$ParsedSuccessfully)) {
  cat("Unparsed IDs:\n")
  print(nn[ParsedSuccessfully == FALSE])
}


# ============================================================
# 2. Base GeneID mapping
# ============================================================

db_ids <- keys(
  org.Hs.eg.db,
  keytype = "ENTREZID"
)

nn[
  ,
  BasePresentInOrgHs :=
    BaseGeneID %chin% db_ids
]

valid_base <- unique(
  nn[
    BasePresentInOrgHs == TRUE,
    BaseGeneID
  ]
)

sym <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = valid_base,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

nn[
  ,
  BaseCanonicalSymbol :=
    unname(sym[BaseGeneID])
]

nn[
  ,
  BaseNumericRowPresent :=
    BaseGeneID %in% rownames(expr)
]


# ============================================================
# 3. Group structure
# ============================================================

group_summary <- nn[
  ,
  .(
    NumLocusRows = .N,
    RawSymbols =
      paste(
        sort(unique(GeneSymbol_raw)),
        collapse = "|"
      ),
    LocusIDs =
      paste(
        GeneID_raw,
        collapse = "|"
      ),
    BasePresentInOrgHs =
      any(BasePresentInOrgHs),
    BaseNumericRowPresent =
      any(BaseNumericRowPresent),
    CanonicalSymbol =
      paste(
        unique(
          BaseCanonicalSymbol[
            !is.na(BaseCanonicalSymbol)
          ]
        ),
        collapse = "|"
      )
  ),
  by = BaseGeneID
]

setorder(
  group_summary,
  -NumLocusRows,
  BaseGeneID
)


# ============================================================
# 4. Pairwise expression within locus groups
# ============================================================

safe_cor <- function(x, y) {

  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x,
      y,
      method = "pearson"
    )
  )
}

pair_list <- list()
k <- 0L

for (base_id in unique(nn$BaseGeneID)) {

  ids <- nn[
    BaseGeneID == base_id,
    GeneID_raw
  ]

  if (length(ids) < 2)
    next

  cmb <- combn(
    ids,
    2,
    simplify = FALSE
  )

  for (pair in cmb) {

    x <- expr[pair[1], ]
    y <- expr[pair[2], ]

    k <- k + 1L

    pair_list[[k]] <- data.table(
      BaseGeneID = base_id,
      Row1 = pair[1],
      Row2 = pair[2],

      ExactIdentical =
        identical(
          as.numeric(x),
          as.numeric(y)
        ),

      Pearson =
        safe_cor(x, y),

      MaxAbsDifference =
        max(abs(x - y)),

      MedianAbsDifference =
        median(abs(x - y)),

      ZeroPatternAgreement =
        mean(
          (x == 0) ==
          (y == 0)
        )
    )
  }
}

pairwise <- if (length(pair_list)) {
  rbindlist(pair_list)
} else {
  data.table()
}


# ============================================================
# 5. Compare locus rows with unsuffixed numeric row
# ============================================================

base_compare_list <- list()
k <- 0L

for (i in seq_len(nrow(nn))) {

  if (!nn$BaseNumericRowPresent[i])
    next

  locus_id <- nn$GeneID_raw[i]
  base_id  <- nn$BaseGeneID[i]

  x <- expr[locus_id, ]
  y <- expr[base_id, ]

  k <- k + 1L

  base_compare_list[[k]] <- data.table(
    LocusGeneID = locus_id,
    BaseGeneID = base_id,

    ExactIdentical =
      identical(
        as.numeric(x),
        as.numeric(y)
      ),

    Pearson =
      safe_cor(x, y),

    MaxAbsDifference =
      max(abs(x - y)),

    MedianAbsDifference =
      median(abs(x - y))
  )
}

base_compare <- if (length(base_compare_list)) {
  rbindlist(base_compare_list)
} else {
  data.table()
}


# ============================================================
# 6. Save
# ============================================================

fwrite(
  nn,
  file.path(
    out_dir,
    "GSE57945_locus_suffixed_geneids_annotated.tsv"
  ),
  sep = "\t"
)

fwrite(
  group_summary,
  file.path(
    out_dir,
    "GSE57945_locus_suffixed_group_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  pairwise,
  file.path(
    out_dir,
    "GSE57945_locus_suffixed_pairwise_expression.tsv"
  ),
  sep = "\t"
)

fwrite(
  base_compare,
  file.path(
    out_dir,
    "GSE57945_locus_vs_unsuffixed_expression.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 7. Report
# ============================================================

cat("============================================================\n")
cat("LOCUS ID STRUCTURE\n")
cat("============================================================\n")

cat(
  "Locus rows                  :",
  nrow(nn),
  "\n"
)

cat(
  "Unique base GeneIDs         :",
  uniqueN(nn$BaseGeneID),
  "\n"
)

cat(
  "Base IDs in org.Hs.eg.db    :",
  uniqueN(
    nn[
      BasePresentInOrgHs == TRUE,
      BaseGeneID
    ]
  ),
  "\n"
)

cat(
  "Base numeric row also exists:",
  uniqueN(
    nn[
      BaseNumericRowPresent == TRUE,
      BaseGeneID
    ]
  ),
  "\n"
)

cat("\nGroup-size distribution:\n")

print(
  group_summary[
    ,
    .N,
    by = NumLocusRows
  ][
    order(NumLocusRows)
  ]
)


cat("\n============================================================\n")
cat("LOCUS-LOCUS EXPRESSION RELATIONSHIP\n")
cat("============================================================\n")

if (nrow(pairwise) > 0) {

  cat(
    "Pairwise comparisons:",
    nrow(pairwise),
    "\n"
  )

  cat(
    "Exact duplicates    :",
    sum(pairwise$ExactIdentical),
    "\n"
  )

  valid <- pairwise[
    !is.na(Pearson)
  ]

  if (nrow(valid) > 0) {

    cat(
      "Pearson min        :",
      min(valid$Pearson),
      "\n"
    )

    cat(
      "Pearson median     :",
      median(valid$Pearson),
      "\n"
    )

    cat(
      "Pearson max        :",
      max(valid$Pearson),
      "\n"
    )
  }
}


cat("\n============================================================\n")
cat("LOCUS vs UNSUFFIXED BASE ROW\n")
cat("============================================================\n")

cat(
  "Comparisons:",
  nrow(base_compare),
  "\n"
)

if (nrow(base_compare) > 0) {

  cat(
    "Exact duplicates:",
    sum(base_compare$ExactIdentical),
    "\n"
  )

  valid2 <- base_compare[
    !is.na(Pearson)
  ]

  if (nrow(valid2) > 0) {

    cat(
      "Pearson median:",
      median(valid2$Pearson),
      "\n"
    )
  }
}


cat("\nExamples:\n")

print(
  head(
    group_summary,
    30
  )
)

cat("\nIMPORTANT:\n")
cat(
  "- No .chr suffix was removed from the expression matrix.\n",
  "- No locus rows were summed or collapsed.\n",
  "- This is a forensic audit only.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_LOCUS_SUFFIXED_GENEID_AUDIT_COMPLETED\n")
cat("============================================================\n")
