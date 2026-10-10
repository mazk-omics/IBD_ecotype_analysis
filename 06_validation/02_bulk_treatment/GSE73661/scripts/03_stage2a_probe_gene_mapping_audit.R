suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(hugene10sttranscriptcluster.db)
})

# ============================================================
# GSE73661 Stage 2A
# Affymetrix transcript-cluster -> gene mapping audit
#
# Source:
#   33,252 GPL6244 transcript-cluster IDs
#
# Strategy:
#   1. ChipDb establishes probe/transcript-cluster -> Entrez
#   2. org.Hs.eg.db establishes current Entrez -> symbol
#   3. classify ambiguity/collisions
#
# NO expression collapse
# NO averaging
# NO probe selection
# NO expression modification
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

expr_file <- file.path(
  root,
  "prepared/01_authoritative_matrix",
  "GSE73661_authoritative_log2_RMA_33252x178.rds"
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
      x == ""
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


count_unique_nonmissing <- function(x) {

  x <- clean_chr(x)

  length(
    unique(
      x[
        !is.na(x)
      ]
    )
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 2A: probe-to-gene mapping audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load frozen authoritative matrix
# ============================================================

expr <- readRDS(
  expr_file
)

check(
  is.matrix(expr),
  "Authoritative expression is a matrix"
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

probe_ids <- as.character(
  rownames(expr)
)

check(
  length(probe_ids) == 33252,
  "33,252 source transcript-cluster IDs recovered"
)

check(
  uniqueN(probe_ids) == 33252,
  "All source transcript-cluster IDs are unique"
)


# ============================================================
# 2. Package provenance
# ============================================================

package_versions <- data.table(
  Package = c(
    "R",
    "AnnotationDbi",
    "org.Hs.eg.db",
    "hugene10sttranscriptcluster.db"
  ),

  Version = c(
    paste(
      R.version$major,
      R.version$minor,
      sep = "."
    ),

    as.character(
      packageVersion(
        "AnnotationDbi"
      )
    ),

    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    ),

    as.character(
      packageVersion(
        "hugene10sttranscriptcluster.db"
      )
    )
  )
)


cat("============================================================\n")
cat("ANNOTATION PACKAGE VERSIONS\n")
cat("============================================================\n")

print(
  package_versions
)


fwrite(
  package_versions,
  file.path(
    out_dir,
    "GSE73661_annotation_package_versions.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 3. Verify ChipDb understands source IDs
# ============================================================

chip_keytypes <- AnnotationDbi::keytypes(
  hugene10sttranscriptcluster.db
)

chip_columns <- AnnotationDbi::columns(
  hugene10sttranscriptcluster.db
)


check(
  "PROBEID" %in% chip_keytypes,
  "ChipDb supports PROBEID keytype"
)

check(
  "ENTREZID" %in% chip_columns,
  "ChipDb provides ENTREZID"
)

check(
  "SYMBOL" %in% chip_columns,
  "ChipDb provides SYMBOL"
)


chip_keys <- AnnotationDbi::keys(
  hugene10sttranscriptcluster.db,
  keytype = "PROBEID"
)


check(
  length(chip_keys) > 0,
  "ChipDb probe universe loaded"
)


source_in_chip <- probe_ids %in%
  chip_keys


cat(
  "[INFO] Source IDs present in ChipDb: ",
  sum(source_in_chip),
  " / ",
  length(probe_ids),
  "\n",
  sep = ""
)


# ============================================================
# 4. Probe -> historical ChipDb gene annotation
# ============================================================

chip_map <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  keytype = "PROBEID",
  columns = c(
    "ENTREZID",
    "SYMBOL"
  )
)

chip_map <- as.data.table(
  chip_map
)

setnames(
  chip_map,
  "PROBEID",
  "SourceProbeID"
)


chip_map[
  ,
  SourceProbeID :=
    as.character(
      SourceProbeID
    )
]

chip_map[
  ,
  ENTREZID :=
    clean_chr(
      ENTREZID
    )
]

chip_map[
  ,
  SYMBOL :=
    clean_chr(
      SYMBOL
    )
]


# ============================================================
# 5. Summarize mappings per source probe
# ============================================================

probe_map <- chip_map[
  ,
  .(
    ChipEntrez =
      collapse_unique(
        ENTREZID
      ),

    N_ChipEntrez =
      count_unique_nonmissing(
        ENTREZID
      ),

    ChipSymbol =
      collapse_unique(
        SYMBOL
      ),

    N_ChipSymbol =
      count_unique_nonmissing(
        SYMBOL
      )
  ),
  by = SourceProbeID
]


source_manifest <- data.table(
  SourceRow =
    seq_along(
      probe_ids
    ),

  SourceProbeID =
    probe_ids
)


probe_map <- merge(
  source_manifest,
  probe_map,
  by = "SourceProbeID",
  all.x = TRUE,
  sort = FALSE
)


setorder(
  probe_map,
  SourceRow
)


check(
  nrow(probe_map) == 33252,
  "Mapping audit contains all 33,252 source probes"
)

check(
  identical(
    probe_map$SourceProbeID,
    probe_ids
  ),
  "Mapping audit preserves source probe order"
)


probe_map[
  is.na(N_ChipEntrez),
  N_ChipEntrez := 0L
]

probe_map[
  is.na(N_ChipSymbol),
  N_ChipSymbol := 0L
]


# ============================================================
# 6. Resolve unique Entrez IDs against current org.Hs.eg.db
# ============================================================

unique_entrez <- unique(
  probe_map[
    N_ChipEntrez == 1L &
      !is.na(ChipEntrez),
    ChipEntrez
  ]
)


current_map_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique_entrez,
  keytype = "ENTREZID",
  columns = c(
    "SYMBOL",
    "GENENAME"
  )
)


current_map_raw <- as.data.table(
  current_map_raw
)


current_map <- current_map_raw[
  ,
  .(
    CurrentSymbol =
      collapse_unique(
        SYMBOL
      ),

    N_CurrentSymbol =
      count_unique_nonmissing(
        SYMBOL
      ),

    CurrentGeneName =
      collapse_unique(
        GENENAME
      )
  ),
  by = ENTREZID
]


probe_map[
  current_map,
  on = c(
    "ChipEntrez" =
      "ENTREZID"
  ),
  `:=`(
    CurrentSymbol =
      i.CurrentSymbol,

    N_CurrentSymbol =
      i.N_CurrentSymbol,

    CurrentGeneName =
      i.CurrentGeneName
  )
]


probe_map[
  is.na(N_CurrentSymbol),
  N_CurrentSymbol := 0L
]


# ============================================================
# 7. Mapping status
# ============================================================

probe_map[
  ,
  MappingStatus :=
    fcase(

      N_ChipEntrez == 0L,
      "NO_ENTREZ_MAPPING",

      N_ChipEntrez > 1L,
      "MULTI_ENTREZ_AMBIGUOUS",

      N_ChipEntrez == 1L &
        N_CurrentSymbol == 0L,
      "UNIQUE_ENTREZ_NO_CURRENT_SYMBOL",

      N_ChipEntrez == 1L &
        N_CurrentSymbol == 1L,
      "UNIQUE_ENTREZ_CURRENT_SYMBOL",

      default =
        "REVIEW"
    )
]


status_summary <- probe_map[
  ,
  .N,
  by = MappingStatus
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("PROBE MAPPING STATUS\n")
cat("============================================================\n")

print(
  status_summary
)


check(
  sum(
    probe_map$MappingStatus ==
      "REVIEW"
  ) == 0,
  "No source probes fall into unclassified REVIEW status"
)


# ============================================================
# 8. Compare historical chip symbol with current symbol
# ============================================================

probe_map[
  ,
  SymbolConcordance :=
    fcase(

      MappingStatus !=
        "UNIQUE_ENTREZ_CURRENT_SYMBOL",
      NA_character_,

      is.na(ChipSymbol),
      "CHIP_SYMBOL_MISSING",

      ChipSymbol ==
        CurrentSymbol,
      "SAME",

      default =
        "CHANGED"
    )
]


symbol_summary <- probe_map[
  !is.na(SymbolConcordance),
  .N,
  by = SymbolConcordance
][
  order(
    -N
  )
]


cat("\n============================================================\n")
cat("CHIP SYMBOL vs CURRENT SYMBOL\n")
cat("============================================================\n")

print(
  symbol_summary
)


# ============================================================
# 9. Gene-level collision audit
#
# Multiple source probe IDs mapping to the same current gene.
# DO NOT collapse them at Stage 2A.
# ============================================================

current_gene_probe_counts <- probe_map[
  MappingStatus ==
    "UNIQUE_ENTREZ_CURRENT_SYMBOL",
  .(
    N_SourceProbes =
      .N,

    SourceProbeIDs =
      paste(
        SourceProbeID,
        collapse = ";"
      )
  ),
  by = .(
    ChipEntrez,
    CurrentSymbol
  )
][
  order(
    -N_SourceProbes,
    CurrentSymbol
  )
]


gene_collision_summary <- data.table(
  Metric = c(
    "Unique_current_genes",
    "Genes_with_1_source_probe",
    "Genes_with_multiple_source_probes",
    "Maximum_source_probes_per_gene"
  ),

  N = c(
    nrow(
      current_gene_probe_counts
    ),

    sum(
      current_gene_probe_counts$N_SourceProbes ==
        1L
    ),

    sum(
      current_gene_probe_counts$N_SourceProbes >
        1L
    ),

    max(
      current_gene_probe_counts$N_SourceProbes
    )
  )
)


cat("\n============================================================\n")
cat("MULTIPLE-PROBE PER GENE AUDIT\n")
cat("============================================================\n")

print(
  gene_collision_summary
)


cat("\nProbe-count distribution per mapped gene:\n")

print(
  current_gene_probe_counts[
    ,
    .N,
    by = N_SourceProbes
  ][
    order(
      N_SourceProbes
    )
  ]
)


# ============================================================
# 10. Source expression remains untouched
# ============================================================

check(
  identical(
    rownames(expr),
    probe_ids
  ),
  "Authoritative expression probe IDs remained unchanged"
)

check(
  identical(
    dim(expr),
    c(
      33252L,
      178L
    )
  ),
  "Authoritative expression dimensions remained unchanged"
)


# ============================================================
# 11. Save audits
# ============================================================

fwrite(
  probe_map,
  file.path(
    out_dir,
    "GSE73661_stage2a_probe_gene_mapping_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  status_summary,
  file.path(
    out_dir,
    "GSE73661_stage2a_mapping_status_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  symbol_summary,
  file.path(
    out_dir,
    "GSE73661_stage2a_symbol_concordance_summary.tsv"
  ),
  sep = "\t",
  na = "NA"
)


fwrite(
  current_gene_probe_counts,
  file.path(
    out_dir,
    "GSE73661_stage2a_current_gene_probe_counts.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  gene_collision_summary,
  file.path(
    out_dir,
    "GSE73661_stage2a_gene_collision_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 12. Final Stage-2A summary
# ============================================================

summary <- data.table(
  Metric = c(
    "Source_probes",
    "Source_probes_in_ChipDb",
    "No_Entrez_mapping",
    "Multi_Entrez_ambiguous",
    "Unique_Entrez_no_current_symbol",
    "Unique_Entrez_current_symbol",
    "Unique_current_genes",
    "Genes_with_multiple_source_probes",
    "Maximum_source_probes_per_gene",
    "Review_rows"
  ),

  N = c(
    length(probe_ids),

    sum(source_in_chip),

    probe_map[
      MappingStatus ==
        "NO_ENTREZ_MAPPING",
      .N
    ],

    probe_map[
      MappingStatus ==
        "MULTI_ENTREZ_AMBIGUOUS",
      .N
    ],

    probe_map[
      MappingStatus ==
        "UNIQUE_ENTREZ_NO_CURRENT_SYMBOL",
      .N
    ],

    probe_map[
      MappingStatus ==
        "UNIQUE_ENTREZ_CURRENT_SYMBOL",
      .N
    ],

    nrow(
      current_gene_probe_counts
    ),

    sum(
      current_gene_probe_counts$N_SourceProbes >
        1L
    ),

    max(
      current_gene_probe_counts$N_SourceProbes
    ),

    probe_map[
      MappingStatus ==
        "REVIEW",
      .N
    ]
  )
)


cat("\n============================================================\n")
cat("STAGE 2A SUMMARY\n")
cat("============================================================\n")

print(
  summary
)


fwrite(
  summary,
  file.path(
    out_dir,
    "GSE73661_stage2a_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 13. Final gate
# ============================================================

critical_pass <- all(
  c(
    nrow(probe_map) == 33252,
    identical(
      probe_map$SourceProbeID,
      probe_ids
    ),
    probe_map[
      MappingStatus ==
        "REVIEW",
      .N
    ] == 0,
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

  cat("[PASS] All 33,252 source probes audited\n")
  cat("[PASS] Probe order preserved\n")
  cat("[PASS] Probe-to-Entrez mapping structure classified\n")
  cat("[PASS] Current gene symbols audited\n")
  cat("[PASS] Multiple-probe-per-gene collisions quantified\n")
  cat("[PASS] No expression values modified\n")
  cat("[PASS] No probe collapsing performed\n")
  cat("[PASS] No representative probe selected yet\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE73661_STAGE2A_PROBE_GENE_MAPPING_AUDIT\n")

} else {

  cat("[REVIEW] Stage 2A requires review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE73661_STAGE2A_PROBE_GENE_MAPPING_AUDIT\n")
}

cat("============================================================\n")
