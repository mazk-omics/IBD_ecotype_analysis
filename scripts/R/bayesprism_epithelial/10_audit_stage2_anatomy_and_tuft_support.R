#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1
)

# ==============================================================================
# IBD EcoTyper
# Stage 2 epithelial BayesPrism
#
# QC-05
#   1. Validate official GSE193677 metadata
#   2. Audit anatomical specificity of epithelial components
#   3. Audit Tuft theta across anatomical regions
#   4. Audit disease / inflammation structure
#   5. Compare Tuft theta against Tuft-specific transcript support
#
# IMPORTANT
#   - Diagnostic QC only
#   - No BayesPrism fitting
#   - No modification of existing posterior results
# ==============================================================================


# ==============================================================================
# 1. Paths
# ==============================================================================

PROJECT_ROOT <- "/home/mazekai/IBD_EcoTyper"

RUN_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "20_bayesprism_epithelial_full_run"
)

INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "17_bayesprism_epithelial_input"
)

QC_DIR <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "21_bayesprism_epithelial_qc"
)

SOFT_FILE <- file.path(
  PROJECT_ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "15_bayesprism_broad14_sanity_audit",
  "GSE193677_family.soft.gz"
)

RESULT_RDS <- file.path(
  RUN_DIR,
  "bayesprism_epithelial7_result.rds"
)

MIXTURE_RDS <- file.path(
  INPUT_DIR,
  "bayesprism_stage2_epithelial_pseudomixture.rds"
)

MASS_QC_TSV <- file.path(
  QC_DIR,
  "stage2_fractional_mass_audit_by_sample.tsv"
)

TUFT_GENE_QC_TSV <- file.path(
  QC_DIR,
  "stage2_tuft_discriminative_support_by_gene.tsv"
)


dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Outputs
# ==============================================================================

OUT_METADATA <- file.path(
  QC_DIR,
  "GSE193677_metadata_used_for_stage2_qc.tsv"
)

OUT_SAMPLE <- file.path(
  QC_DIR,
  "stage2_epithelial_anatomy_sample_level.tsv"
)

OUT_REGION <- file.path(
  QC_DIR,
  "stage2_epithelial_theta_by_region.tsv"
)

OUT_DISEASE <- file.path(
  QC_DIR,
  "stage2_epithelial_theta_by_disease.tsv"
)

OUT_INFLAMMATION <- file.path(
  QC_DIR,
  "stage2_epithelial_theta_by_inflammation.tsv"
)

OUT_REGION_DISEASE <- file.path(
  QC_DIR,
  "stage2_epithelial_theta_by_region_disease.tsv"
)

OUT_REGION_INFLAMMATION <- file.path(
  QC_DIR,
  "stage2_epithelial_theta_by_region_inflammation.tsv"
)

OUT_MARKER_COR <- file.path(
  QC_DIR,
  "stage2_tuft_marker_support_correlations.tsv"
)

OUT_LOG <- file.path(
  QC_DIR,
  "stage2_anatomy_tuft_support_audit.log"
)

OUT_SESSION <- file.path(
  QC_DIR,
  "stage2_anatomy_tuft_support_sessionInfo.txt"
)


# ==============================================================================
# 3. Logging
# ==============================================================================

log_con <- file(
  OUT_LOG,
  open = "wt"
)

sink(
  log_con,
  type = "output",
  split = TRUE
)

sink(
  log_con,
  type = "message"
)

on.exit({

  while (sink.number(type = "message") > 0) {
    sink(type = "message")
  }

  while (sink.number(type = "output") > 0) {
    sink(type = "output")
  }

  close(log_con)

}, add = TRUE)


# ==============================================================================
# 4. Helpers
# ==============================================================================

assert_true <- function(condition, message) {

  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}


write_tsv <- function(x, path) {

  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


normalize_key <- function(x) {

  x <- tolower(
    trimws(x)
  )

  x <- gsub(
    "[^a-z0-9]+",
    "_",
    x
  )

  x <- gsub(
    "^_|_$",
    "",
    x
  )

  x
}


safe_cpm <- function(numerator, denominator) {

  out <- rep(
    NA_real_,
    length(denominator)
  )

  valid <- is.finite(denominator) &
    denominator > 0

  out[valid] <- numerator[valid] /
    denominator[valid] * 1e6

  out
}


cat("============================================================\n")
cat("Stage 2 anatomy + Tuft support QC\n")
cat("============================================================\n\n")

cat(
  "START:",
  format(Sys.time()),
  "\n\n"
)


# ==============================================================================
# 5. Load BayesPrism result
# ==============================================================================

suppressPackageStartupMessages(
  library(BayesPrism)
)


assert_true(
  file.exists(RESULT_RDS),
  paste(
    "Missing Stage 2 BayesPrism result:",
    RESULT_RDS
  )
)


bp <- readRDS(
  RESULT_RDS
)


assert_true(
  is(bp, "BayesPrism"),
  "Stage 2 result is not a BayesPrism object."
)


theta_initial <- bp@posterior.initial.cellType@theta
theta_final <- bp@posterior.theta_f@theta


assert_true(
  identical(
    rownames(theta_initial),
    rownames(theta_final)
  ),
  "Initial/final theta sample IDs differ."
)


assert_true(
  identical(
    colnames(theta_initial),
    colnames(theta_final)
  ),
  "Initial/final theta component names differ."
)


components <- colnames(
  theta_final
)


bp_ids <- rownames(
  theta_final
)


cat(
  "[PASS] Theta:",
  nrow(theta_final),
  "samples x",
  ncol(theta_final),
  "components\n\n"
)


# ==============================================================================
# 6. Load Stage 1 epithelial pseudo-mixture
# ==============================================================================

assert_true(
  file.exists(MIXTURE_RDS),
  paste(
    "Missing Stage 2 pseudo-mixture:",
    MIXTURE_RDS
  )
)


mixture <- readRDS(
  MIXTURE_RDS
)


assert_true(
  identical(
    rownames(mixture),
    bp_ids
  ),
  "Pseudo-mixture sample IDs/order differ from BayesPrism result."
)


cat(
  "[PASS] Pseudo-mixture:",
  nrow(mixture),
  "x",
  ncol(mixture),
  "\n\n"
)


# ==============================================================================
# 7. Validate and parse the actual GSE193677 family SOFT
# ==============================================================================

assert_true(
  file.exists(SOFT_FILE),
  paste(
    "Missing GEO family SOFT:",
    SOFT_FILE
  )
)


cat(
  "Using GEO SOFT:\n",
  SOFT_FILE,
  "\n\n"
)


cat(
  "Parsing GEO sample metadata...\n"
)


con <- gzfile(
  SOFT_FILE,
  open = "rt"
)

soft_lines <- readLines(
  con,
  warn = FALSE
)

close(con)


sample_starts <- grep(
  "^\\^SAMPLE = ",
  soft_lines
)


assert_true(
  length(sample_starts) == 2490L,
  paste(
    "Expected 2490 GEO SAMPLE blocks but found",
    length(sample_starts)
  )
)


sample_ends <- c(
  sample_starts[-1] - 1L,
  length(soft_lines)
)


records <- vector(
  "list",
  length(sample_starts)
)


for (i in seq_along(sample_starts)) {

  block <- soft_lines[
    sample_starts[i]:sample_ends[i]
  ]


  gsm <- sub(
    "^\\^SAMPLE =\\s*",
    "",
    block[1]
  )


  title_lines <- grep(
    "^!Sample_title = ",
    block,
    value = TRUE
  )


  title <- if (length(title_lines) > 0) {

    sub(
      "^!Sample_title =\\s*",
      "",
      title_lines[1]
    )

  } else {

    NA_character_
  }


  characteristic_lines <- grep(
    "^!Sample_characteristics_ch1 = ",
    block,
    value = TRUE
  )


  characteristic_values <- sub(
    "^!Sample_characteristics_ch1 =\\s*",
    "",
    characteristic_lines
  )


  rec <- list(
    gsm = gsm,
    title = title
  )


  for (value in characteristic_values) {

    if (!grepl(":", value, fixed = TRUE)) {
      next
    }


    key <- normalize_key(
      sub(
        ":.*$",
        "",
        value
      )
    )


    field_value <- trimws(
      sub(
        "^[^:]+:",
        "",
        value
      )
    )


    if (is.null(getElement(rec, key))) {
      rec <- append(
        rec,
        setNames(
          list(field_value),
          key
        )
      )
    }
  }


  records[i] <- list(rec)
}


all_keys <- unique(
  unlist(
    lapply(
      records,
      names
    )
  )
)


meta_columns <- lapply(
  all_keys,
  function(key) {

    vapply(
      records,
      function(rec) {

        value <- getElement(
          rec,
          key
        )

        if (
          is.null(value) ||
          length(value) == 0
        ) {
          return(NA_character_)
        }

        as.character(
          value[1]
        )

      },
      character(1)
    )
  }
)


names(meta_columns) <- all_keys


meta <- as.data.frame(
  meta_columns,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "[PASS] Parsed GEO samples:",
  nrow(meta),
  "\n"
)

cat(
  "Metadata fields:\n",
  paste(
    colnames(meta),
    collapse = ", "
  ),
  "\n\n"
)


# ==============================================================================
# 8. Validate actual metadata schema
# ==============================================================================

required_fields <- c(
  "gsm",
  "title",
  "regionre",
  "diseasetypere",
  "ibd_disease",
  "typere"
)


missing_fields <- setdiff(
  required_fields,
  colnames(meta)
)


assert_true(
  length(missing_fields) == 0,
  paste(
    "Missing required GEO fields:",
    paste(
      missing_fields,
      collapse = ", "
    )
  )
)


cat(
  "[PASS] Actual GEO metadata schema validated\n\n"
)


# ==============================================================================
# 9. Derive exact BayesPrism sample ID from GEO title
#
# Actual title examples:
#
# MSCCR_reGRID_10_Biopsy_19, CD participants,LeftColon NonI tissue
# MSCCR_reGRID_1323_Biopsy_2057, CD participants,Ileum I tissue
#
# Therefore the string before the first comma is the exact sample ID.
# ==============================================================================

geo_sample_id <- trimws(
  sub(
    ",.*$",
    "",
    meta$title
  )
)


assert_true(
  !anyNA(geo_sample_id),
  "NA GEO sample IDs derived from title."
)


assert_true(
  !anyDuplicated(geo_sample_id),
  "Duplicated GEO sample IDs derived from title."
)


cat(
  "First 10 GEO sample IDs derived from title:\n"
)

print(
  head(
    geo_sample_id,
    10
  )
)


cat("\n")


# ==============================================================================
# 10. Exact match to BayesPrism samples
# ==============================================================================

idx <- match(
  bp_ids,
  geo_sample_id
)


n_matched <- sum(
  !is.na(idx)
)


cat(
  "Exact GEO-title/BayesPrism matches:",
  n_matched,
  "/",
  length(bp_ids),
  "\n\n"
)


assert_true(
  n_matched == length(bp_ids),
  paste(
    length(bp_ids) - n_matched,
    "BayesPrism samples failed exact GEO-title matching."
  )
)


meta <- meta[
  idx,
  ,
  drop = FALSE
]


geo_sample_id <- geo_sample_id[
  idx
]


assert_true(
  identical(
    geo_sample_id,
    bp_ids
  ),
  "GEO metadata ordering differs from BayesPrism sample ordering."
)


meta$sample_id <- bp_ids


cat(
  "[PASS] GEO metadata matched exactly:",
  nrow(meta),
  "/",
  length(bp_ids),
  "\n\n"
)


# ==============================================================================
# 11. Clean biological metadata
# ==============================================================================

region_raw <- trimws(
  meta$regionre
)


region_lower <- tolower(
  region_raw
)


region_broad <- rep(
  "Other",
  length(region_raw)
)


region_broad[
  grepl(
    "ileum",
    region_lower
  )
] <- "Ileum"


colon_pattern <- paste(
  c(
    "colon",
    "rectum",
    "rectal",
    "sigmoid",
    "cecum",
    "caecum",
    "transverse",
    "ascending",
    "descending"
  ),
  collapse = "|"
)


region_broad[
  grepl(
    colon_pattern,
    region_lower
  )
] <- "Colon_rectum"


disease <- trimws(
  meta$ibd_disease
)


disease_type_re <- trimws(
  meta$diseasetypere
)


type_re <- trimws(
  meta$typere
)


inflammation <- ifelse(
  type_re == "I",
  "Inflamed",
  ifelse(
    type_re == "NonI",
    "Non-inflamed",
    type_re
  )
)


assert_true(
  all(
    disease %in%
      c(
        "CD",
        "UC",
        "Control"
      )
  ),
  paste(
    "Unexpected ibd_disease values:",
    paste(
      setdiff(
        unique(disease),
        c(
          "CD",
          "UC",
          "Control"
        )
      ),
      collapse = ", "
    )
  )
)


assert_true(
  all(
    type_re %in%
      c(
        "I",
        "NonI"
      )
  ),
  paste(
    "Unexpected typere values:",
    paste(
      setdiff(
        unique(type_re),
        c(
          "I",
          "NonI"
        )
      ),
      collapse = ", "
    )
  )
)


has_patient_id <- grepl(
  "reGRID_[0-9]+",
  bp_ids
)


assert_true(
  all(has_patient_id),
  "Could not extract reGRID patient ID from all samples."
)


patient_id <- sub(
  ".*(reGRID_[0-9]+).*",
  "\\1",
  bp_ids
)


# ==============================================================================
# 12. Metadata audit
# ==============================================================================

cat("============================================================\n")
cat("GSE193677 METADATA AUDIT\n")
cat("============================================================\n\n")


cat("First 10 matched records:\n\n")

print(
  data.frame(
    sample_id = head(bp_ids, 10),
    title = head(meta$title, 10),
    region = head(region_raw, 10),
    disease = head(disease, 10),
    typere = head(type_re, 10),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


cat("\nRaw region counts:\n")

print(
  sort(
    table(
      region_raw,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


cat("\nBroad region counts:\n")

print(
  sort(
    table(
      region_broad,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


cat("\nDisease counts (ibd_disease):\n")

print(
  sort(
    table(
      disease,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


cat("\nInflammation counts (typere):\n")

print(
  sort(
    table(
      inflammation,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


cat("\nComposite diseasetypere counts:\n")

print(
  sort(
    table(
      disease_type_re,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


cat("\nRegion x disease:\n")

print(
  addmargins(
    table(
      region_broad,
      disease
    )
  )
)


cat("\nRegion x inflammation:\n")

print(
  addmargins(
    table(
      region_broad,
      inflammation
    )
  )
)


cat("\n")


# ==============================================================================
# 13. Save exact metadata used
# ==============================================================================

metadata_used <- data.frame(
  sample_id = bp_ids,
  gsm = meta$gsm,
  title = meta$title,
  patient_id = patient_id,
  region = region_raw,
  region_broad = region_broad,
  disease = disease,
  disease_type_re = disease_type_re,
  typere = type_re,
  inflammation = inflammation,
  stringsAsFactors = FALSE
)


write_tsv(
  metadata_used,
  OUT_METADATA
)


# ==============================================================================
# 14. Fractional pseudo-mixture reliability
# ==============================================================================

assert_true(
  file.exists(MASS_QC_TSV),
  paste(
    "Missing fractional-mass QC:",
    MASS_QC_TSV
  )
)


mass_qc <- read.delim(
  MASS_QC_TSV,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


mass_ids <- getElement(
  mass_qc,
  "sample_id"
)


mass_idx <- match(
  bp_ids,
  mass_ids
)


assert_true(
  !anyNA(mass_idx),
  "Could not match all samples to fractional-mass QC."
)


relative_loss <- getElement(
  mass_qc,
  "relative_mass_loss"
)[mass_idx]


allocated_mass <- getElement(
  mass_qc,
  "posterior_allocated_mass"
)[mass_idx]


original_mass <- getElement(
  mass_qc,
  "original_mass"
)[mass_idx]


reliable <- (
  allocated_mass > 0 &
    relative_loss <= 0.01
)


cat(
  "Numerically reliable Stage 2 samples:",
  sum(reliable),
  "/",
  length(reliable),
  "\n"
)

cat(
  "Excluded by diagnostic <=1% mass-loss rule:",
  sum(!reliable),
  "\n\n"
)


# ==============================================================================
# 15. Tuft marker support in Stage 1 epithelial pseudo-mixture
# ==============================================================================

epi_mass <- rowSums(
  mixture
)


assert_true(
  all(
    is.finite(epi_mass) &
      epi_mass > 0
  ),
  "Some Stage 1 epithelial pseudo-mixtures have zero/non-finite total mass."
)


canonical_tuft <- c(
  "POU2F3",
  "PLCB2",
  "GFI1B",
  "AVIL",
  "PTGS1",
  "SH2D6"
)


canonical_present <- intersect(
  canonical_tuft,
  colnames(mixture)
)


assert_true(
  length(canonical_present) >= 4,
  "Too few canonical Tuft markers present in Stage 2 gene universe."
)


canonical_mass <- rowSums(
  mixture[
    ,
    canonical_present,
    drop = FALSE
  ]
)


canonical_CPM <- safe_cpm(
  canonical_mass,
  epi_mass
)


# ==============================================================================
# 16. Data-driven Tuft-specific gene support
# ==============================================================================

assert_true(
  file.exists(TUFT_GENE_QC_TSV),
  paste(
    "Missing Tuft gene QC table:",
    TUFT_GENE_QC_TSV
  )
)


tuft_gene_qc <- read.delim(
  TUFT_GENE_QC_TSV,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


gene_qc_gene <- getElement(
  tuft_gene_qc,
  "gene"
)

reference_share <- getElement(
  tuft_gene_qc,
  "tuft_reference_share"
)

reference_fold <- getElement(
  tuft_gene_qc,
  "tuft_vs_max_other_fold"
)

reference_top <- getElement(
  tuft_gene_qc,
  "reference_top_component"
)


strong80 <- gene_qc_gene[
  is.finite(reference_share) &
    reference_share >= 0.80
]


strong5x <- gene_qc_gene[
  is.finite(reference_fold) &
    reference_top == "Tuft" &
    reference_fold >= 5
]


strong80 <- intersect(
  strong80,
  colnames(mixture)
)


strong5x <- intersect(
  strong5x,
  colnames(mixture)
)


assert_true(
  length(strong80) > 0,
  "No Tuft reference-share >=80% genes found."
)


assert_true(
  length(strong5x) > 0,
  "No >=5x Tuft-specific genes found."
)


strong80_mass <- rowSums(
  mixture[
    ,
    strong80,
    drop = FALSE
  ]
)


strong5x_mass <- rowSums(
  mixture[
    ,
    strong5x,
    drop = FALSE
  ]
)


strong80_CPM <- safe_cpm(
  strong80_mass,
  epi_mass
)


strong5x_CPM <- safe_cpm(
  strong5x_mass,
  epi_mass
)


cat(
  "Canonical Tuft genes:",
  paste(
    canonical_present,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Reference-share >=80% Tuft genes:",
  length(strong80),
  "\n"
)

cat(
  ">=5x Tuft-specific genes:",
  length(strong5x),
  "\n\n"
)


# ==============================================================================
# 17. Construct sample-level audit table
# ==============================================================================

sample_table <- data.frame(
  sample_id = bp_ids,
  gsm = meta$gsm,
  patient_id = patient_id,
  region = region_raw,
  region_broad = region_broad,
  disease = disease,
  disease_type_re = disease_type_re,
  inflammation = inflammation,
  reliable_fractional_mass = reliable,
  stage1_epithelial_mass = original_mass,
  stage2_integerized_mass = allocated_mass,
  stage2_relative_mass_loss = relative_loss,
  canonical_Tuft_CPM = canonical_CPM,
  strong80_Tuft_CPM = strong80_CPM,
  strong5x_Tuft_CPM = strong5x_CPM,
  stringsAsFactors = FALSE
)


for (component in components) {

  initial_col <- paste0(
    "initial_",
    component
  )

  final_col <- paste0(
    "final_",
    component
  )

  sample_table[
    ,
    initial_col
  ] <- theta_initial[
    ,
    component
  ]

  sample_table[
    ,
    final_col
  ] <- theta_final[
    ,
    component
  ]
}


write_tsv(
  sample_table,
  OUT_SAMPLE
)


# ==============================================================================
# 18. Generic grouped theta summary
# ==============================================================================

summarise_theta <- function(df, group_cols) {

  complete_group <- complete.cases(
    df[
      ,
      group_cols,
      drop = FALSE
    ]
  )


  df2 <- df[
    complete_group,
    ,
    drop = FALSE
  ]


  interaction_args <- c(
    as.list(
      df2[
        ,
        group_cols,
        drop = FALSE
      ]
    ),
    list(
      drop = TRUE,
      sep = " | ",
      lex.order = TRUE
    )
  )


  group_key <- do.call(
    interaction,
    interaction_args
  )


  groups <- split(
    seq_len(nrow(df2)),
    group_key
  )


  output_rows <- list()

  counter <- 1L


  for (group_indices in groups) {

    sub_df <- df2[
      group_indices,
      ,
      drop = FALSE
    ]


    for (stage in c(
      "initial",
      "final"
    )) {

      for (component in components) {

        value_col <- paste0(
          stage,
          "_",
          component
        )


        x <- sub_df[
          ,
          value_col,
          drop = TRUE
        ]


        row_out <- data.frame(
          stage = stage,
          component = component,
          n_samples = nrow(sub_df),
          n_patients = length(
            unique(
              sub_df$patient_id
            )
          ),
          mean = mean(
            x,
            na.rm = TRUE
          ),
          q25 = unname(
            quantile(
              x,
              0.25,
              na.rm = TRUE
            )
          ),
          median = median(
            x,
            na.rm = TRUE
          ),
          q75 = unname(
            quantile(
              x,
              0.75,
              na.rm = TRUE
            )
          ),
          stringsAsFactors = FALSE
        )


        for (group_col in group_cols) {

          row_out[
            ,
            group_col
          ] <- sub_df[
            1,
            group_col,
            drop = TRUE
          ]
        }


        output_rows[counter] <- list(
          row_out
        )

        counter <- counter + 1L
      }
    }
  }


  result <- do.call(
    rbind,
    output_rows
  )


  result[
    ,
    c(
      group_cols,
      "stage",
      "component",
      "n_samples",
      "n_patients",
      "mean",
      "q25",
      "median",
      "q75"
    ),
    drop = FALSE
  ]
}


# ==============================================================================
# 19. Restrict biological interpretation to reliable samples
# ==============================================================================

reliable_df <- sample_table[
  sample_table$reliable_fractional_mass,
  ,
  drop = FALSE
]


assert_true(
  nrow(reliable_df) == sum(reliable),
  "Reliable-sample filtering mismatch."
)


region_summary <- summarise_theta(
  reliable_df,
  "region_broad"
)


disease_summary <- summarise_theta(
  reliable_df,
  "disease"
)


inflammation_summary <- summarise_theta(
  reliable_df,
  "inflammation"
)


region_disease_summary <- summarise_theta(
  reliable_df,
  c(
    "region_broad",
    "disease"
  )
)


region_inflammation_summary <- summarise_theta(
  reliable_df,
  c(
    "region_broad",
    "inflammation"
  )
)


write_tsv(
  region_summary,
  OUT_REGION
)

write_tsv(
  disease_summary,
  OUT_DISEASE
)

write_tsv(
  inflammation_summary,
  OUT_INFLAMMATION
)

write_tsv(
  region_disease_summary,
  OUT_REGION_DISEASE
)

write_tsv(
  region_inflammation_summary,
  OUT_REGION_INFLAMMATION
)


# ==============================================================================
# 20. Tuft theta vs marker-support correlation
# ==============================================================================

rel_idx <- which(
  reliable
)


make_cor_row <- function(
  stage_name,
  marker_name,
  theta_values,
  marker_values,
  region_name = "ALL"
) {

  valid <- is.finite(theta_values) &
    is.finite(marker_values)


  data.frame(
    region = region_name,
    theta = stage_name,
    marker_score = marker_name,
    n = sum(valid),
    spearman = cor(
      theta_values[valid],
      marker_values[valid],
      method = "spearman"
    ),
    stringsAsFactors = FALSE
  )
}


cor_rows <- list()

counter <- 1L


marker_scores <- list(
  canonical_Tuft_CPM = canonical_CPM,
  strong80_Tuft_CPM = strong80_CPM,
  strong5x_Tuft_CPM = strong5x_CPM
)


for (marker_name in names(marker_scores)) {

  marker_values <- getElement(
    marker_scores,
    marker_name
  )


  cor_rows[counter] <- list(
    make_cor_row(
      "initial_Tuft",
      marker_name,
      theta_initial[
        rel_idx,
        "Tuft"
      ],
      marker_values[
        rel_idx
      ]
    )
  )

  counter <- counter + 1L


  cor_rows[counter] <- list(
    make_cor_row(
      "final_Tuft",
      marker_name,
      theta_final[
        rel_idx,
        "Tuft"
      ],
      marker_values[
        rel_idx
      ]
    )
  )

  counter <- counter + 1L
}


for (region_name in sort(
  unique(
    reliable_df$region_broad
  )
)) {

  region_idx <- which(
    reliable &
      region_broad == region_name
  )


  if (length(region_idx) < 10) {
    next
  }


  for (marker_name in names(marker_scores)) {

    marker_values <- getElement(
      marker_scores,
      marker_name
    )


    cor_rows[counter] <- list(
      make_cor_row(
        "final_Tuft",
        marker_name,
        theta_final[
          region_idx,
          "Tuft"
        ],
        marker_values[
          region_idx
        ],
        region_name
      )
    )

    counter <- counter + 1L
  }
}


cor_table <- do.call(
  rbind,
  cor_rows
)


write_tsv(
  cor_table,
  OUT_MARKER_COR
)


# ==============================================================================
# 21. Console: final theta by anatomical region
# ==============================================================================

cat("============================================================\n")
cat("FINAL THETA MEDIANS BY BROAD REGION\n")
cat("Reliable samples only\n")
cat("============================================================\n\n")


final_region <- region_summary[
  region_summary$stage == "final",
  ,
  drop = FALSE
]


for (region_name in sort(
  unique(
    final_region$region_broad
  )
)) {

  cat(
    "\nREGION:",
    region_name,
    "\n"
  )


  region_rows <- final_region[
    final_region$region_broad == region_name,
    ,
    drop = FALSE
  ]


  print(
    region_rows[
      ,
      c(
        "component",
        "n_samples",
        "n_patients",
        "median",
        "q25",
        "q75"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


# ==============================================================================
# 22. Console: Tuft by region
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("TUFT FINAL THETA BY REGION\n")
cat("============================================================\n\n")


for (region_name in sort(
  unique(
    reliable_df$region_broad
  )
)) {

  x <- reliable_df[
    reliable_df$region_broad == region_name,
    "final_Tuft",
    drop = TRUE
  ]


  cat(
    region_name,
    ": n=",
    length(x),
    " median=",
    signif(
      median(x),
      6
    ),
    " IQR=[",
    signif(
      quantile(
        x,
        0.25
      ),
      6
    ),
    ", ",
    signif(
      quantile(
        x,
        0.75
      ),
      6
    ),
    "] >10%=",
    signif(
      mean(
        x > 0.10
      ),
      4
    ),
    " >20%=",
    signif(
      mean(
        x > 0.20
      ),
      4
    ),
    "\n",
    sep = ""
  )
}


# ==============================================================================
# 23. Console: Tuft by disease
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("TUFT FINAL THETA BY DISEASE\n")
cat("============================================================\n\n")


for (group_name in c(
  "CD",
  "UC",
  "Control"
)) {

  x <- reliable_df[
    reliable_df$disease == group_name,
    "final_Tuft",
    drop = TRUE
  ]


  if (length(x) == 0) {
    next
  }


  cat(
    group_name,
    ": n=",
    length(x),
    " median=",
    signif(
      median(x),
      6
    ),
    " IQR=[",
    signif(
      quantile(
        x,
        0.25
      ),
      6
    ),
    ", ",
    signif(
      quantile(
        x,
        0.75
      ),
      6
    ),
    "]\n",
    sep = ""
  )
}


# ==============================================================================
# 24. Console: Tuft / Stem-TA by inflammation
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("TUFT / STEM-TA BY INFLAMMATION\n")
cat("============================================================\n\n")


for (group_name in c(
  "Inflamed",
  "Non-inflamed"
)) {

  sub_df <- reliable_df[
    reliable_df$inflammation == group_name,
    ,
    drop = FALSE
  ]


  if (nrow(sub_df) == 0) {
    next
  }


  cat(
    group_name,
    ": n=",
    nrow(sub_df),
    " Tuft median=",
    signif(
      median(
        sub_df$final_Tuft
      ),
      6
    ),
    " StemTA median=",
    signif(
      median(
        sub_df$final_Stem_TA_progenitor
      ),
      6
    ),
    "\n",
    sep = ""
  )
}


# ==============================================================================
# 25. Console: marker support correlations
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("TUFT THETA vs TUFT-SPECIFIC MARKER SUPPORT\n")
cat("============================================================\n\n")


print(
  cor_table,
  row.names = FALSE
)


# ==============================================================================
# 26. Directional anatomical checks
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("ANATOMICAL DIRECTIONAL CHECKS\n")
cat("============================================================\n\n")


available_regions <- unique(
  reliable_df$region_broad
)


if (
  all(
    c(
      "Ileum",
      "Colon_rectum"
    ) %in%
      available_regions
  )
) {

  ileum <- reliable_df[
    reliable_df$region_broad == "Ileum",
    ,
    drop = FALSE
  ]


  colon <- reliable_df[
    reliable_df$region_broad == "Colon_rectum",
    ,
    drop = FALSE
  ]


  cat(
    "Ileal_absorptive median:\n"
  )

  cat(
    "  Ileum        = ",
    median(
      ileum$final_Ileal_absorptive
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Colon/rectum = ",
    median(
      colon$final_Ileal_absorptive
    ),
    "\n\n",
    sep = ""
  )


  cat(
    "Paneth median:\n"
  )

  cat(
    "  Ileum        = ",
    median(
      ileum$final_Paneth
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Colon/rectum = ",
    median(
      colon$final_Paneth
    ),
    "\n\n",
    sep = ""
  )


  cat(
    "Colonic_absorptive median:\n"
  )

  cat(
    "  Ileum        = ",
    median(
      ileum$final_Colonic_absorptive
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Colon/rectum = ",
    median(
      colon$final_Colonic_absorptive
    ),
    "\n\n",
    sep = ""
  )


  cat(
    "Goblet median:\n"
  )

  cat(
    "  Ileum        = ",
    median(
      ileum$final_Goblet
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Colon/rectum = ",
    median(
      colon$final_Goblet
    ),
    "\n\n",
    sep = ""
  )


  cat(
    "Tuft median:\n"
  )

  cat(
    "  Ileum        = ",
    median(
      ileum$final_Tuft
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Colon/rectum = ",
    median(
      colon$final_Tuft
    ),
    "\n\n",
    sep = ""
  )
}


# ==============================================================================
# 27. Final audit summary
# ==============================================================================

cat("============================================================\n")
cat("QC-05 SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Total samples                 :",
  nrow(sample_table),
  "\n"
)

cat(
  "Reliable samples              :",
  nrow(reliable_df),
  "\n"
)

cat(
  "Unique patients               :",
  length(
    unique(
      sample_table$patient_id
    )
  ),
  "\n"
)

cat(
  "Broad regions                 :",
  paste(
    sort(
      unique(
        sample_table$region_broad
      )
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Diseases                      :",
  paste(
    sort(
      unique(
        sample_table$disease
      )
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Inflammation states           :",
  paste(
    sort(
      unique(
        sample_table$inflammation
      )
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Canonical Tuft marker genes   :",
  length(canonical_present),
  "\n"
)

cat(
  "Tuft share>=80 genes          :",
  length(strong80),
  "\n"
)

cat(
  "Tuft >=5x genes               :",
  length(strong5x),
  "\n"
)


capture.output(
  sessionInfo(),
  file = OUT_SESSION
)


cat("\nOutput directory:\n")
cat(QC_DIR, "\n")


cat(
  "\nEND:",
  format(Sys.time()),
  "\n"
)

