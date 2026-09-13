library(data.table)

# ============================================================
# GSE282122 final_analysis -> reference identity mapping
#
# Current ontology:
#   18 KEEP identities
#   2 CANDIDATE identities: Glial, Tuft
#
# REVIEW:
#   Cycling stroma
#   Cycling MNP
#
# EXCLUDE:
#   gd SOX4pos T
#   gd T
#   MAIT
#   Non ileal M like
# ============================================================


# ------------------------------------------------------------
# 1. Input / output
# ------------------------------------------------------------

input_file <- "GSE282122_celltype_mapping_template.tsv"

output_file <- "GSE282122_celltype_mapping.tsv"


# ------------------------------------------------------------
# 2. Frozen KEEP ontology
# ------------------------------------------------------------

keep_identities <- c(
  "B",
  "Plasma",
  "CD4_T",
  "CD8_T",
  "NK",
  "ILC",
  "Monocyte_Macrophage",
  "DC",
  "Mast",
  "Fibroblast",
  "Pericyte",
  "Endothelial",
  "Ileal_absorptive",
  "Colonic_absorptive",
  "Goblet",
  "Stem_TA_progenitor",
  "Paneth",
  "Enteroendocrine"
)

candidate_identities <- c(
  "Glial",
  "Tuft"
)


# ------------------------------------------------------------
# 3. Helper function
# ------------------------------------------------------------

make_map <- function(labels,
                     reference_identity,
                     status = "KEEP",
                     action = "merge",
                     note = "") {

  data.table(
    final_analysis = labels,
    reference_identity = reference_identity,
    status = status,
    action = action,
    note = note
  )
}


# ------------------------------------------------------------
# 4. KEEP mappings
# ------------------------------------------------------------

mapping_list <- list(

  # -------------------------
  # B cells
  # -------------------------

  make_map(
    c(
      "Memory B",
      "Follicular B",
      "Intermediate B",
      "Cycling B",
      "FCRL4pos memory B",
      "IFN resp memory B",
      "CCL22pos memory B",
      "GClike B",
      "IgGhi IgAhi memory B",
      "CCL3pos CCL4pos follicular B",
      "IGLC6hi memory B"
    ),
    "B"
  ),


  # -------------------------
  # Plasma cells
  # -------------------------

  make_map(
    c(
      "IgApos plasma",
      "IgGpos CXCR4hi plasma",
      "IgGpos CXCR4lo plasma",
      "Plasmablast",
      "IgApos IFN resp plasma"
    ),
    "Plasma"
  ),


  # -------------------------
  # CD4 T cells
  # -------------------------

  make_map(
    c(
      "CD4 FOSpos T",
      "CD4 FOShi T",
      "CD4 naive T",
      "Th22",
      "Tph Tfh CXCL13pos",
      "Tph Tfh",
      "CD4 KLF2hi T",
      "CD4 IKZF2hi TNFRSF18lo Treg",
      "CD4 IKZF2lo LAG3pos Treg",
      "CD4 KLF2int T",
      "Th17",
      "CD4 IKZF2hi TNFRSF18hi Treg",
      "Th1 17 GZMAhi",
      "Th1 17 GZMApos",
      "CD4 TWIST1 Treg",
      "Th1",
      "CD4 HSPhi Treg",
      "CD4 TNFSF13Bhi T",
      "CD4 HSPhi CD70pos Treg",
      "Th1 17 22"
    ),
    "CD4_T"
  ),


  # -------------------------
  # CD8 T cells
  # -------------------------

  make_map(
    c(
      "CD8 IL7Rhi T",
      "CD8 GZMKhi T",
      "CD8 GZMKint T",
      "CD8 EGR1hi T",
      "CD8 naive",
      "CD8 FGFBP2pos T",
      "CD8 IL17Apos IL26pos IL23Rpos T",
      "CD8 CTLA4hi TIGIThi T",
      "CD8 TNFhi IFNGhi IL2pos T",
      "CD8 cycling T"
    ),
    "CD8_T"
  ),


  # -------------------------
  # NK
  # -------------------------

  make_map(
    "NK",
    "NK",
    action = "keep"
  ),


  # -------------------------
  # ILC
  # -------------------------

  make_map(
    "ILC",
    "ILC",
    action = "keep"
  ),


  # -------------------------
  # Monocyte / macrophage
  # -------------------------

  make_map(
    c(
      "S100A8 A9hi mono",
      "C1Qhi IL1Bhi macro",
      "C1Qhi IL1Blo macro",
      "S100A8 A9hi TNFhi IL6pos mono"
    ),
    "Monocyte_Macrophage"
  ),


  # -------------------------
  # DC
  # -------------------------

  make_map(
    c(
      "CD1Chi DC",
      "XCR1pos DC",
      "LAMP3pos DC",
      "pDC",
      "LAMP3pos IL1Bpos DC"
    ),
    "DC"
  ),


  # -------------------------
  # Mast
  # -------------------------

  make_map(
    "Mast",
    "Mast",
    action = "keep"
  ),


  # -------------------------
  # Fibroblast
  # -------------------------

  make_map(
    c(
      "ABCA8pos WNT2Bpos FOSlo fibroblast",
      "SOX6pos POSTNpos fibroblast",
      "ABCA8pos WNT2Bpos FOShi fibroblast",
      "C3hi RSPO3pos fibroblast",
      "SOX6pos POSTNpos NRG1hi NPYpos fibroblast",
      "Myofibroblast",
      "C3hi CCL19pos fibroblast",
      "THY1pos FAPpos PDPNpos fibroblast"
    ),
    "Fibroblast"
  ),


  # -------------------------
  # Pericyte
  # -------------------------

  make_map(
    c(
      "CD74hi HLADRB1hi arterial pericyte",
      "NOTCH3hi TNClo pericyte",
      "NOTCH3hi MYH11pos pericyte",
      "NOTCH3hi TNChi LOXL2pos pericyte",
      "NOTCH3hi TNCint CCL19pos pericyte",
      "CD74hi HLADRB1hi venous pericyte"
    ),
    "Pericyte"
  ),


  # -------------------------
  # Endothelial
  # -------------------------

  make_map(
    c(
      "Arterial endothelium",
      "Venous endothelium",
      "Lymphatic endothelium"
    ),
    "Endothelial"
  ),


  # -------------------------
  # Ileal absorptive
  # -------------------------

  make_map(
    c(
      "Ileal enterocyte",
      "Ileal enterocyte PLCG2hi",
      "Ileal RPShi enterocyte",
      "Ileal enterocyte DUOX2pos LCN2pos",
      "Ileal BEST4 OTOP2",
      "Ileal enterocyte DUOX2pos MUC6pos PGCpos"
    ),
    "Ileal_absorptive"
  ),


  # -------------------------
  # Colonic absorptive
  # -------------------------

  make_map(
    c(
      "Non ileal undiff enterocyte",
      "Non ileal CT colonocyte",
      "Non ileal enterocyte PLCG2hi",
      "Non ileal enterocyte DUOXA2pos",
      "Non ileal BEST4 OTOP2",
      "Non ileal enterocyte DUOXA2pos CXCL11pos"
    ),
    "Colonic_absorptive"
  ),


  # -------------------------
  # Goblet
  # -------------------------

  make_map(
    c(
      "Non ileal goblet",
      "Ileal goblet"
    ),
    "Goblet"
  ),


  # -------------------------
  # Stem / TA progenitor
  # -------------------------

  make_map(
    c(
      "Non ileal TA",
      "Non ileal LGR5pos stem",
      "Ileal TA",
      "Ileal LGR5pos stem"
    ),
    "Stem_TA_progenitor"
  ),


  # -------------------------
  # Paneth
  # -------------------------

  make_map(
    c(
      "Non ileal paneth",
      "Ileal paneth"
    ),
    "Paneth"
  ),


  # -------------------------
  # Enteroendocrine
  # -------------------------

  make_map(
    c(
      "Non ileal EEC DDCpos CHGApos",
      "Non ileal EEC GCGpos PYYpos",
      "Non ileal EEC NEUROG3pos",
      "Ileal EEC CHGAhi CHGBhi",
      "Ileal EEC CHGAlo CHGBlo"
    ),
    "Enteroendocrine"
  ),


  # ==========================================================
  # 5. Candidate identities
  # ==========================================================

  make_map(
    "Glia",
    "Glial",
    status = "CANDIDATE",
    action = "hold",
    note = "Candidate identity; cross-dataset stability test pending"
  ),

  make_map(
    c(
      "Non ileal tuft",
      "Ileal tuft"
    ),
    "Tuft",
    status = "CANDIDATE",
    action = "hold",
    note = "Candidate identity; cross-dataset stability test pending"
  ),


  # ==========================================================
  # 6. REVIEW
  # ==========================================================

  make_map(
    "Cycling stroma",
    NA_character_,
    status = "REVIEW",
    action = "review",
    note = "Cycling stromal identity; parent lineage not sufficiently resolved"
  ),

  make_map(
    "Cycling MNP",
    NA_character_,
    status = "REVIEW",
    action = "review",
    note = "Cycling MNP; parent monocyte/macrophage/DC identity requires review"
  ),


  # ==========================================================
  # 7. EXCLUDE
  # ==========================================================

  make_map(
    c(
      "gd SOX4pos T",
      "gd T",
      "MAIT",
      "Non ileal M like"
    ),
    NA_character_,
    status = "EXCLUDE",
    action = "exclude",
    note = "Not represented in current frozen KEEP ontology"
  )
)


# ------------------------------------------------------------
# 8. Assemble mapping dictionary
# ------------------------------------------------------------

mapping <- rbindlist(mapping_list, use.names = TRUE)


cat("========================================\n")
cat("MAPPING DICTIONARY SELF-CHECK\n")
cat("========================================\n\n")


# ------------------------------------------------------------
# 9. Mapping dictionary must contain exactly 109 source labels
# ------------------------------------------------------------

cat("Mapping dictionary rows:", nrow(mapping), "\n")

if (nrow(mapping) != 109L) {
  stop(
    "ERROR: mapping dictionary contains ",
    nrow(mapping),
    " rows; expected exactly 109."
  )
}


# ------------------------------------------------------------
# 10. No duplicate source labels
# ------------------------------------------------------------

dup <- mapping[
  duplicated(final_analysis) |
    duplicated(final_analysis, fromLast = TRUE),
  final_analysis
]

if (length(dup) > 0L) {

  cat("\nDuplicated mapping labels:\n")
  print(unique(dup))

  stop("ERROR: duplicated source labels in mapping dictionary.")
}

cat("PASS: 109 unique mapping labels\n")


# ------------------------------------------------------------
# 11. Validate KEEP / CANDIDATE identities
# ------------------------------------------------------------

bad_keep <- mapping[
  status == "KEEP" &
    !reference_identity %in% keep_identities
]

if (nrow(bad_keep) > 0L) {

  print(bad_keep)

  stop(
    "ERROR: invalid identity found among KEEP mappings."
  )
}


bad_candidate <- mapping[
  status == "CANDIDATE" &
    !reference_identity %in% candidate_identities
]

if (nrow(bad_candidate) > 0L) {

  print(bad_candidate)

  stop(
    "ERROR: invalid candidate identity."
  )
}

cat("PASS: ontology whitelist validation\n")


# ------------------------------------------------------------
# 12. Read H5AD-derived template
# ------------------------------------------------------------

template <- fread(
  input_file,
  na.strings = c("", "NA")
)

required_columns <- c(
  "code",
  "final_analysis",
  "n_cells"
)

missing_columns <- setdiff(
  required_columns,
  colnames(template)
)

if (length(missing_columns) > 0L) {

  stop(
    "ERROR: template missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 13. Validate template size
# ------------------------------------------------------------

if (nrow(template) != 109L) {

  stop(
    "ERROR: input template contains ",
    nrow(template),
    " rows; expected 109."
  )
}

cat("PASS: input template contains 109 rows\n")


# ------------------------------------------------------------
# 14. Compare dictionary vs actual H5AD labels
# ------------------------------------------------------------

missing_from_mapping <- setdiff(
  template$final_analysis,
  mapping$final_analysis
)

extra_in_mapping <- setdiff(
  mapping$final_analysis,
  template$final_analysis
)


if (length(missing_from_mapping) > 0L) {

  cat("\nLabels present in H5AD/template but missing from mapping:\n")

  print(missing_from_mapping)

  stop("ERROR: mapping dictionary is incomplete.")
}


if (length(extra_in_mapping) > 0L) {

  cat("\nLabels present in mapping but absent from H5AD/template:\n")

  print(extra_in_mapping)

  stop("ERROR: mapping dictionary contains invalid labels.")
}


cat("PASS: mapping labels match H5AD labels exactly (109/109)\n")


# ------------------------------------------------------------
# 15. Join mapping onto template
# ------------------------------------------------------------

# Remove empty mapping columns that may already exist in template
drop_cols <- intersect(
  c(
    "reference_identity",
    "status",
    "action",
    "note"
  ),
  colnames(template)
)

if (length(drop_cols) > 0L) {
  template[, (drop_cols) := NULL]
}


result <- merge(
  template,
  mapping,
  by = "final_analysis",
  all.x = TRUE,
  sort = FALSE
)


# Restore original template order using code
setorder(result, code)


# ------------------------------------------------------------
# 16. Final completeness checks
# ------------------------------------------------------------

if (anyNA(result$status)) {

  print(
    result[
      is.na(status),
      .(code, final_analysis)
    ]
  )

  stop("ERROR: unmapped rows remain after merge.")
}


if (nrow(result) != 109L) {

  stop(
    "ERROR: output mapping does not contain exactly 109 rows."
  )
}


# ------------------------------------------------------------
# 17. Cell-count preservation check
# ------------------------------------------------------------

expected_cells <- 987743L

observed_cells <- sum(result$n_cells)

cat(
  "Cells represented by mapping:",
  observed_cells,
  "\n"
)

if (observed_cells != expected_cells) {

  stop(
    "ERROR: cell count mismatch. Expected ",
    expected_cells,
    ", observed ",
    observed_cells,
    "."
  )
}

cat("PASS: all 987743 cells accounted for\n")


# ------------------------------------------------------------
# 18. Anatomical sanity checks
# ------------------------------------------------------------

bad_ileal <- result[
  status == "KEEP" &
    grepl("^Ileal ", final_analysis) &
    grepl(
      "enterocyte|BEST4",
      final_analysis,
      ignore.case = TRUE
    ) &
    reference_identity != "Ileal_absorptive"
]

if (nrow(bad_ileal) > 0L) {

  print(bad_ileal)

  stop(
    "ERROR: possible ileal absorptive misclassification."
  )
}


bad_colonic <- result[
  status == "KEEP" &
    grepl("^Non ileal ", final_analysis) &
    grepl(
      "enterocyte|colonocyte|BEST4",
      final_analysis,
      ignore.case = TRUE
    ) &
    reference_identity != "Colonic_absorptive"
]

if (nrow(bad_colonic) > 0L) {

  print(bad_colonic)

  stop(
    "ERROR: possible non-ileal/colonic absorptive misclassification."
  )
}

cat("PASS: ileal/non-ileal absorptive sanity check\n")


# ------------------------------------------------------------
# 19. Status/action consistency
# ------------------------------------------------------------

bad_candidate_action <- result[
  status == "CANDIDATE" &
    action != "hold"
]

bad_exclude_action <- result[
  status == "EXCLUDE" &
    action != "exclude"
]

bad_review_action <- result[
  status == "REVIEW" &
    action != "review"
]

if (
  nrow(bad_candidate_action) > 0L ||
  nrow(bad_exclude_action) > 0L ||
  nrow(bad_review_action) > 0L
) {

  stop(
    "ERROR: inconsistent status/action combination."
  )
}

cat("PASS: status/action consistency\n")


# ------------------------------------------------------------
# 20. Write final mapping
# ------------------------------------------------------------

fwrite(
  result,
  output_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 21. Audit summaries
# ------------------------------------------------------------

cat("\n========================================\n")
cat("STATUS SUMMARY\n")
cat("========================================\n")

print(
  result[
    ,
    .(
      n_source_labels = .N,
      n_cells = sum(n_cells)
    ),
    by = status
  ][order(status)]
)


cat("\n========================================\n")
cat("KEEP IDENTITY SUMMARY\n")
cat("========================================\n")

keep_summary <- result[
  status == "KEEP",
  .(
    n_source_labels = .N,
    n_cells = sum(n_cells)
  ),
  by = reference_identity
][order(-n_cells)]

print(keep_summary)


cat("\n========================================\n")
cat("CANDIDATE / REVIEW / EXCLUDE\n")
cat("========================================\n")

print(
  result[
    status != "KEEP",
    .(
      code,
      final_analysis,
      n_cells,
      reference_identity,
      status,
      action,
      note
    )
  ]
)


cat("\n========================================\n")
cat("SUCCESS\n")
cat("========================================\n")

cat("Final mapping written to:\n")
cat(output_file, "\n")
cat("\nAll 109 source labels accounted for.\n")