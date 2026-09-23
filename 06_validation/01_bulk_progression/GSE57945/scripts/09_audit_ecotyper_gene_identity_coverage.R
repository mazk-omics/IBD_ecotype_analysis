suppressPackageStartupMessages({
  library(data.table)
})

gse_base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

eco_base <- paste0(
  "/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/",
  "MSCCR_CORE12_CTSPEC_FORMAL/",
  "Cell_type_specific_genes/Cell_States/discovery"
)

decision_file <- file.path(
  gse_base,
  "prepared/02_gene_harmonization",
  "GSE57945_canonical_gene_decision_table.tsv"
)

out_dir <- file.path(
  gse_base,
  "prepared/02_gene_harmonization"
)

cell_types <- c(
  "B",
  "CD4_T",
  "CD8_T",
  "Endothelial",
  "Enteroendocrine",
  "Fibroblast",
  "Glial",
  "Goblet",
  "ILC",
  "Mast",
  "Monocyte_Macrophage",
  "Pericyte"
)

cat("============================================================\n")
cat("EcoTyper NMF gene identity-aware coverage audit\n")
cat("============================================================\n\n")

decision <- fread(
  decision_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

included <- decision[
  CanonicalAction %chin% c(
    "KEEP_DIRECT_CURRENT",
    "KEEP_SAFE_REDIRECT"
  )
]

stopifnot(
  uniqueN(included$FinalGeneID) == nrow(included),
  uniqueN(included$FinalSymbol) == nrow(included)
)

canonical_symbols <- included$FinalSymbol


# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

read_genes <- function(path) {

  x <- fread(
    path,
    select = 1,
    colClasses = "character",
    showProgress = FALSE,
    check.names = FALSE
  )

  unique(
    trimws(
      x[[1]]
    )
  )
}


# ------------------------------------------------------------
# Union of EcoTyper current NMF-input genes
# ------------------------------------------------------------

gene_by_ct <- list()

for (ct in cell_types) {

  f <- file.path(
    eco_base,
    ct,
    "expression_top_genes_scaled.txt"
  )

  gene_by_ct[[ct]] <- read_genes(f)
}

eco_genes <- sort(
  unique(
    unlist(
      gene_by_ct,
      use.names = FALSE
    )
  )
)

cat("EcoTyper union genes:", length(eco_genes), "\n\n")


# ------------------------------------------------------------
# Map each EcoTyper symbol to GSE57945 gene identity
# ------------------------------------------------------------

map_list <- vector(
  "list",
  length(eco_genes)
)

for (i in seq_along(eco_genes)) {

  g <- eco_genes[i]

  # ----------------------------------------------------------
  # A. Exact current canonical symbol
  # ----------------------------------------------------------

  hit_final <- included[
    FinalSymbol == g
  ]

  if (nrow(hit_final) == 1) {

    map_list[[i]] <- data.table(
      EcoTyperGene = g,
      MatchClass = "EXACT_CANONICAL_MATCH",
      FinalGeneID = hit_final$FinalGeneID,
      FinalSymbol = hit_final$FinalSymbol,
      SourceGeneID = hit_final$GeneID_raw,
      SourceRawSymbol = hit_final$GeneSymbol_raw
    )

    next
  }


  # ----------------------------------------------------------
  # B. Old/raw symbol belonging to an included canonical row
  # ----------------------------------------------------------

  hit_raw_included <- included[
    GeneSymbol_raw == g
  ]

  # only allow a unique 1:1 identity
  if (
    nrow(hit_raw_included) == 1 &&
    !is.na(hit_raw_included$FinalGeneID) &&
    !is.na(hit_raw_included$FinalSymbol)
  ) {

    map_list[[i]] <- data.table(
      EcoTyperGene = g,
      MatchClass = "SAFE_SYMBOL_RENAME",
      FinalGeneID = hit_raw_included$FinalGeneID,
      FinalSymbol = hit_raw_included$FinalSymbol,
      SourceGeneID = hit_raw_included$GeneID_raw,
      SourceRawSymbol = hit_raw_included$GeneSymbol_raw
    )

    next
  }


  # ----------------------------------------------------------
  # C. Locus-suffixed only
  # ----------------------------------------------------------

  locus_hit <- decision[
    CanonicalAction == "EXCLUDE_LOCUS_SUFFIXED" &
    GeneSymbol_raw == g
  ]

  if (nrow(locus_hit) > 0) {

    map_list[[i]] <- data.table(
      EcoTyperGene = g,
      MatchClass = "LOCUS_ONLY",
      FinalGeneID = NA_character_,
      FinalSymbol = NA_character_,
      SourceGeneID = paste(
        locus_hit$GeneID_raw,
        collapse = "|"
      ),
      SourceRawSymbol = g
    )

    next
  }


  # ----------------------------------------------------------
  # D. Other excluded rows carrying that symbol
  # ----------------------------------------------------------

  excluded_hit <- decision[
    !CanonicalAction %chin% c(
      "KEEP_DIRECT_CURRENT",
      "KEEP_SAFE_REDIRECT"
    ) &
    (
      GeneSymbol_raw == g |
      CanonicalSymbol == g
    )
  ]

  if (nrow(excluded_hit) > 0) {

    map_list[[i]] <- data.table(
      EcoTyperGene = g,
      MatchClass = "EXCLUDED_AMBIGUOUS_OR_OBSOLETE",
      FinalGeneID = NA_character_,
      FinalSymbol = NA_character_,
      SourceGeneID = paste(
        unique(excluded_hit$GeneID_raw),
        collapse = "|"
      ),
      SourceRawSymbol = paste(
        unique(excluded_hit$GeneSymbol_raw),
        collapse = "|"
      )
    )

    next
  }


  # ----------------------------------------------------------
  # E. Truly absent
  # ----------------------------------------------------------

  map_list[[i]] <- data.table(
    EcoTyperGene = g,
    MatchClass = "ABSENT_FROM_GSE57945",
    FinalGeneID = NA_character_,
    FinalSymbol = NA_character_,
    SourceGeneID = NA_character_,
    SourceRawSymbol = NA_character_
  )
}

mapping <- rbindlist(
  map_list,
  fill = TRUE
)


# ------------------------------------------------------------
# Global summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("GLOBAL IDENTITY MATCH CLASSES\n")
cat("============================================================\n")

global_summary <- mapping[
  ,
  .N,
  by = MatchClass
][order(-N)]

print(global_summary)


usable_classes <- c(
  "EXACT_CANONICAL_MATCH",
  "SAFE_SYMBOL_RENAME"
)

usable <- mapping[
  MatchClass %chin% usable_classes
]

cat("\nIdentity-resolved usable genes:",
    nrow(usable), "\n")

cat(
  "Identity-aware coverage:",
  sprintf(
    "%.4f",
    nrow(usable) / nrow(mapping)
  ),
  "\n"
)


# ------------------------------------------------------------
# Critical uniqueness
# ------------------------------------------------------------

dup_target <- usable[
  ,
  .N,
  by = FinalGeneID
][
  N > 1
]

cat(
  "Canonical targets receiving >1 EcoTyper symbol:",
  nrow(dup_target),
  "\n"
)


# ------------------------------------------------------------
# Per cell type
# ------------------------------------------------------------

ct_summary_list <- list()

for (ct in cell_types) {

  genes <- gene_by_ct[[ct]]

  m <- mapping[
    EcoTyperGene %chin% genes
  ]

  ct_summary_list[[ct]] <- data.table(
    CellType = ct,
    TotalGenes = length(genes),

    ExactCanonical =
      sum(
        m$MatchClass ==
          "EXACT_CANONICAL_MATCH"
      ),

    SafeSymbolRename =
      sum(
        m$MatchClass ==
          "SAFE_SYMBOL_RENAME"
      ),

    LocusOnly =
      sum(
        m$MatchClass ==
          "LOCUS_ONLY"
      ),

    ExcludedOther =
      sum(
        m$MatchClass ==
          "EXCLUDED_AMBIGUOUS_OR_OBSOLETE"
      ),

    TrulyAbsent =
      sum(
        m$MatchClass ==
          "ABSENT_FROM_GSE57945"
      ),

    IdentityResolved =
      sum(
        m$MatchClass %chin%
          usable_classes
      ),

    IdentityCoverage =
      sum(
        m$MatchClass %chin%
          usable_classes
      ) / length(genes)
  )
}

ct_summary <- rbindlist(ct_summary_list)


cat("\n============================================================\n")
cat("PER-CELL-TYPE IDENTITY-AWARE COVERAGE\n")
cat("============================================================\n")

print(ct_summary)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fwrite(
  mapping,
  file.path(
    out_dir,
    "GSE57945_EcoTyper_NMF_gene_identity_mapping.tsv"
  ),
  sep = "\t"
)

fwrite(
  global_summary,
  file.path(
    out_dir,
    "GSE57945_EcoTyper_NMF_gene_identity_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  ct_summary,
  file.path(
    out_dir,
    "GSE57945_EcoTyper_NMF_identity_coverage_by_celltype.tsv"
  ),
  sep = "\t"
)

fwrite(
  dup_target,
  file.path(
    out_dir,
    "GSE57945_EcoTyper_NMF_duplicate_canonical_targets.tsv"
  ),
  sep = "\t"
)


cat("\n============================================================\n")
cat("UNRESOLVED GENES\n")
cat("============================================================\n")

print(
  mapping[
    !MatchClass %chin% usable_classes,
    .(
      EcoTyperGene,
      MatchClass,
      SourceGeneID,
      SourceRawSymbol
    )
  ],
  nrows = 300
)


cat("\nIMPORTANT:\n")
cat(
  "- SAFE_SYMBOL_RENAME means the EcoTyper symbol uniquely identifies\n",
  "  an included GSE57945 gene row whose current canonical symbol differs.\n",
  "- LOCUS_ONLY features remain unresolved at gene level.\n",
  "- No expression matrix was modified.\n",
  "- Final recovery-signature audit is still required after EcoTyper discovery.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_ECOTYPER_GENE_IDENTITY_COVERAGE_AUDIT\n")
cat("============================================================\n")
