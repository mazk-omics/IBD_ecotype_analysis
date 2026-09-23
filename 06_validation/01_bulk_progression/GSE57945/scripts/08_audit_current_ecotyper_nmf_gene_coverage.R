suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE57945
# Stage 2F:
# Audit coverage of CURRENT EcoTyper pre-NMF feature sets
#
# This is NOT a final state-signature / recovery audit.
# EcoTyper discovery is still running.
# ============================================================

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
cat("GSE57945 vs CURRENT EcoTyper pre-NMF genes\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# Helper: read ONLY first column (gene names)
# ------------------------------------------------------------

read_gene_column <- function(path) {

  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }

  x <- fread(
    path,
    select = 1,
    colClasses = "character",
    showProgress = FALSE,
    check.names = FALSE
  )

  genes <- trimws(x[[1]])

  genes <- genes[
    !is.na(genes) &
    genes != ""
  ]

  genes
}


# ------------------------------------------------------------
# Load GSE57945 decision table
# ------------------------------------------------------------

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

canonical_symbols <- included$FinalSymbol

stopifnot(
  length(canonical_symbols) ==
    uniqueN(canonical_symbols)
)

cat(
  "GSE57945 canonical symbols:",
  length(canonical_symbols),
  "\n\n"
)


# ------------------------------------------------------------
# Audit each EcoTyper cell type
# ------------------------------------------------------------

summary_list <- list()
missing_list <- list()
all_top_genes <- list()

for (ct in cell_types) {

  cat("------------------------------------------------------------\n")
  cat("Cell type:", ct, "\n")

  ct_dir <- file.path(
    eco_base,
    ct
  )

  f_log2 <- file.path(
    ct_dir,
    "expression_top_genes_log2.txt"
  )

  f_scaled <- file.path(
    ct_dir,
    "expression_top_genes_scaled.txt"
  )

  genes_log2 <- read_gene_column(f_log2)
  genes_scaled <- read_gene_column(f_scaled)

  # Structural QC
  n_log2 <- length(genes_log2)
  n_scaled <- length(genes_scaled)

  u_log2 <- uniqueN(genes_log2)
  u_scaled <- uniqueN(genes_scaled)

  duplicate_log2 <- n_log2 - u_log2
  duplicate_scaled <- n_scaled - u_scaled

  exact_same_order <- identical(
    genes_log2,
    genes_scaled
  )

  same_gene_set <- setequal(
    genes_log2,
    genes_scaled
  )

  if (!same_gene_set) {
    stop(
      "Gene set mismatch between log2 and scaled files for: ",
      ct
    )
  }

  if (duplicate_log2 > 0 || duplicate_scaled > 0) {
    warning(
      "Duplicated genes detected for: ",
      ct
    )
  }

  # Use unique gene universe for coverage audit
  top_genes <- unique(genes_scaled)

  all_top_genes[[ct]] <- top_genes

  present <- top_genes %chin% canonical_symbols
  missing <- top_genes[!present]

  summary_list[[ct]] <- data.table(
    CellType = ct,

    Log2Rows = n_log2,
    ScaledRows = n_scaled,

    Log2UniqueGenes = u_log2,
    ScaledUniqueGenes = u_scaled,

    DuplicateGenesLog2 = duplicate_log2,
    DuplicateGenesScaled = duplicate_scaled,

    Log2ScaledSameSet = same_gene_set,
    Log2ScaledSameOrder = exact_same_order,

    NMFInputUniqueGenes = length(top_genes),

    PresentInGSE57945 = sum(present),
    MissingFromGSE57945 = length(missing),

    CoverageFraction =
      sum(present) / length(top_genes)
  )

  if (length(missing) > 0) {

    missing_list[[ct]] <- data.table(
      CellType = ct,
      Gene = missing
    )
  }

  cat(
    "NMF input genes:",
    length(top_genes),
    "\n"
  )

  cat(
    "Present:",
    sum(present),
    "\n"
  )

  cat(
    "Missing:",
    length(missing),
    "\n"
  )

  cat(
    "Coverage:",
    sprintf(
      "%.4f",
      sum(present) / length(top_genes)
    ),
    "\n"
  )
}


summary_dt <- rbindlist(
  summary_list,
  use.names = TRUE,
  fill = TRUE
)

missing_dt <- if (length(missing_list) > 0) {
  rbindlist(
    missing_list,
    use.names = TRUE,
    fill = TRUE
  )
} else {
  data.table(
    CellType = character(),
    Gene = character()
  )
}


# ------------------------------------------------------------
# Union of CURRENT pre-NMF genes across 12 cell types
# ------------------------------------------------------------

union_top <- sort(
  unique(
    unlist(
      all_top_genes,
      use.names = FALSE
    )
  )
)

union_present <-
  union_top %chin% canonical_symbols

union_missing <-
  union_top[!union_present]


cat("\n============================================================\n")
cat("UNION OF CURRENT NMF INPUT GENES\n")
cat("============================================================\n")

cat(
  "Unique EcoTyper NMF-input genes:",
  length(union_top),
  "\n"
)

cat(
  "Present in GSE57945 canonical:",
  sum(union_present),
  "\n"
)

cat(
  "Missing:",
  length(union_missing),
  "\n"
)

cat(
  "Coverage:",
  sprintf(
    "%.4f",
    sum(union_present) / length(union_top)
  ),
  "\n"
)


# ------------------------------------------------------------
# Trace missing genes back to original GSE57945 annotation
# ------------------------------------------------------------

if (length(union_missing) > 0) {

  # A missing EcoTyper gene may correspond to:
  # 1. raw symbol in excluded feature
  # 2. canonical symbol in excluded feature
  # 3. not present anywhere in GSE57945

  trace_list <- vector(
    "list",
    length(union_missing)
  )

  for (i in seq_along(union_missing)) {

    g <- union_missing[i]

    hit <- decision[
      GeneSymbol_raw == g |
      CanonicalSymbol == g |
      FinalSymbol == g
    ]

    if (nrow(hit) == 0) {

      trace_list[[i]] <- data.table(
        Gene = g,
        FoundInOriginalGSE57945 = FALSE,
        GeneID_raw = NA_character_,
        GeneSymbol_raw = NA_character_,
        CanonicalSymbol = NA_character_,
        CanonicalAction =
          "ABSENT_FROM_ORIGINAL_GSE57945",
        CanonicalReason =
          "Gene symbol not found in original GSE57945 annotation"
      )

    } else {

      trace_list[[i]] <- hit[
        ,
        .(
          Gene = g,
          FoundInOriginalGSE57945 = TRUE,
          GeneID_raw,
          GeneSymbol_raw,
          CanonicalSymbol,
          CanonicalAction,
          CanonicalReason
        )
      ]
    }
  }

  missing_trace <- rbindlist(
    trace_list,
    use.names = TRUE,
    fill = TRUE
  )

} else {

  missing_trace <- data.table(
    Gene = character(),
    FoundInOriginalGSE57945 = logical(),
    GeneID_raw = character(),
    GeneSymbol_raw = character(),
    CanonicalSymbol = character(),
    CanonicalAction = character(),
    CanonicalReason = character()
  )
}


# ------------------------------------------------------------
# Missing-gene decision-class summary
# ------------------------------------------------------------

if (nrow(missing_trace) > 0) {

  missing_reason_summary <- missing_trace[
    ,
    .(
      Rows = .N,
      UniqueGenes = uniqueN(Gene)
    ),
    by = CanonicalAction
  ][
    order(-UniqueGenes)
  ]

} else {

  missing_reason_summary <- data.table(
    CanonicalAction = character(),
    Rows = integer(),
    UniqueGenes = integer()
  )
}


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fwrite(
  summary_dt,
  file.path(
    out_dir,
    "GSE57945_current_EcoTyper_NMF_gene_coverage_by_celltype.tsv"
  ),
  sep = "\t"
)

fwrite(
  data.table(
    Gene = union_top,
    PresentInGSE57945Canonical =
      union_top %chin% canonical_symbols
  ),
  file.path(
    out_dir,
    "GSE57945_current_EcoTyper_NMF_gene_union.tsv"
  ),
  sep = "\t"
)

fwrite(
  missing_dt,
  file.path(
    out_dir,
    "GSE57945_current_EcoTyper_NMF_missing_by_celltype.tsv"
  ),
  sep = "\t"
)

fwrite(
  missing_trace,
  file.path(
    out_dir,
    "GSE57945_current_EcoTyper_NMF_missing_gene_trace.tsv"
  ),
  sep = "\t"
)

fwrite(
  missing_reason_summary,
  file.path(
    out_dir,
    "GSE57945_current_EcoTyper_NMF_missing_reason_summary.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PER-CELL-TYPE COVERAGE\n")
cat("============================================================\n")

print(
  summary_dt[
    ,
    .(
      CellType,
      NMFInputUniqueGenes,
      PresentInGSE57945,
      MissingFromGSE57945,
      CoverageFraction,
      Log2ScaledSameSet,
      Log2ScaledSameOrder
    )
  ]
)

cat("\n============================================================\n")
cat("MISSING-GENE REASONS\n")
cat("============================================================\n")

print(missing_reason_summary)


cat("\nIMPORTANT:\n")
cat(
  "- This audits the CURRENT pre-NMF EcoTyper feature sets.\n",
  "- These are NOT final state-marker or recovery signatures.\n",
  "- Final recovery coverage must be re-audited after discovery completes.\n",
  "- No GSE57945 expression matrix was modified.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_CURRENT_ECOTYPER_NMF_GENE_COVERAGE_AUDIT\n")
cat("============================================================\n")
