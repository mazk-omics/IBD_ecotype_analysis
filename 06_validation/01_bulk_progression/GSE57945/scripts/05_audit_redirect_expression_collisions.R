suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE57945
# Stage 2C: expression forensic audit for redirected GeneIDs
#
# Important:
#   - NO expression values are modified
#   - NO rows are collapsed
#   - this is diagnostic only
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE57945"

obj_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE57945_authoritative_expression_322.rds"
)

rescue_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE57945_ncbi_geneid_rescue_classified.tsv"
)

out_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

cat("============================================================\n")
cat("GSE57945 redirect expression collision audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load data
# ============================================================

obj <- readRDS(obj_file)
expr <- obj$expression

rescue <- fread(
  rescue_file,
  colClasses = "character"
)

redirect <- rescue[
  RescueClass == "REDIRECTED_TO_CURRENT"
]

cat("Redirected old GeneIDs :", nrow(redirect), "\n")
cat(
  "Unique current targets :",
  uniqueN(redirect$ReturnedGeneID),
  "\n\n"
)


# ============================================================
# 2. Helper
# ============================================================

safe_cor <- function(x, y, method = "pearson") {

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
      method = method,
      use = "pairwise.complete.obs"
    )
  )
}


# ============================================================
# 3. Compare each redirected row to current target
# ============================================================

audit_list <- vector(
  "list",
  nrow(redirect)
)

for (i in seq_len(nrow(redirect))) {

  old_id <- redirect$QueryGeneID[i]
  new_id <- redirect$ReturnedGeneID[i]

  old_present <- old_id %in% rownames(expr)
  new_present <- new_id %in% rownames(expr)

  if (!old_present) {
    stop(
      "Old GeneID missing from authoritative matrix: ",
      old_id
    )
  }

  x <- expr[old_id, ]

  if (new_present) {

    y <- expr[new_id, ]

    abs_diff <- abs(x - y)

    denom <- pmax(
      abs(x),
      abs(y),
      1e-12
    )

    relative_diff <- abs_diff / denom

    pearson <- safe_cor(
      x,
      y,
      "pearson"
    )

    spearman <- safe_cor(
      x,
      y,
      "spearman"
    )

    exact_identical <- identical(
      as.numeric(x),
      as.numeric(y)
    )

    max_abs_diff <- max(abs_diff)
    median_abs_diff <- median(abs_diff)

    mean_abs_diff <- mean(abs_diff)

    median_relative_diff <-
      median(relative_diff)

    zero_pattern_agreement <-
      mean(
        (x == 0) ==
        (y == 0)
      )

    # Diagnostic classification only
    relationship <- fifelse(
      exact_identical,
      "EXACT_DUPLICATE",

      fifelse(
        !is.na(pearson) &&
        pearson >= 0.9999 &&
        median_relative_diff <= 1e-4,
        "NEAR_DUPLICATE",

        fifelse(
          !is.na(pearson) &&
          pearson >= 0.95,
          "HIGHLY_CORRELATED_DISTINCT",
          "DISTINCT_PROFILE"
        )
      )
    )

  } else {

    pearson <- NA_real_
    spearman <- NA_real_

    exact_identical <- NA
    max_abs_diff <- NA_real_
    median_abs_diff <- NA_real_
    mean_abs_diff <- NA_real_
    median_relative_diff <- NA_real_
    zero_pattern_agreement <- NA_real_

    relationship <-
      "TARGET_NOT_IN_ORIGINAL_MATRIX"
  }

  audit_list[[i]] <- data.table(
    OldGeneID = old_id,
    TargetGeneID = new_id,
    TargetSymbol = redirect$Symbol[i],

    TargetPresentInMatrix = new_present,

    Pearson = pearson,
    Spearman = spearman,

    ExactIdentical = exact_identical,

    MaxAbsDifference = max_abs_diff,
    MeanAbsDifference = mean_abs_diff,
    MedianAbsDifference = median_abs_diff,

    MedianRelativeDifference =
      median_relative_diff,

    ZeroPatternAgreement =
      zero_pattern_agreement,

    RelationshipClass =
      relationship
  )
}

audit <- rbindlist(
  audit_list,
  fill = TRUE
)


# ============================================================
# 4. Add redirect multiplicity
# ============================================================

audit[
  ,
  RedirectMultiplicity :=
    .N,
  by = TargetGeneID
]


# ============================================================
# 5. Group-level collision summary
# ============================================================

group_summary <- audit[
  ,
  .(
    NumOldGeneIDs = .N,

    TargetPresentInMatrix =
      any(TargetPresentInMatrix),

    NumExactDuplicates =
      sum(
        RelationshipClass ==
          "EXACT_DUPLICATE"
      ),

    NumNearDuplicates =
      sum(
        RelationshipClass ==
          "NEAR_DUPLICATE"
      ),

    NumHighlyCorrelatedDistinct =
      sum(
        RelationshipClass ==
          "HIGHLY_CORRELATED_DISTINCT"
      ),

    NumDistinctProfiles =
      sum(
        RelationshipClass ==
          "DISTINCT_PROFILE"
      ),

    NumTargetAbsent =
      sum(
        RelationshipClass ==
          "TARGET_NOT_IN_ORIGINAL_MATRIX"
      ),

    MinPearson =
      if (all(is.na(Pearson)))
        NA_real_
      else
        min(Pearson, na.rm = TRUE),

    MedianPearson =
      if (all(is.na(Pearson)))
        NA_real_
      else
        median(Pearson, na.rm = TRUE)

  ),
  by = .(
    TargetGeneID,
    TargetSymbol
  )
]


# ============================================================
# 6. Save
# ============================================================

fwrite(
  audit,
  file.path(
    out_dir,
    "GSE57945_redirect_expression_pairwise_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  group_summary,
  file.path(
    out_dir,
    "GSE57945_redirect_expression_group_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 7. Console summary
# ============================================================

cat("============================================================\n")
cat("PAIRWISE RELATIONSHIP CLASSES\n")
cat("============================================================\n")

print(
  audit[
    ,
    .N,
    by = RelationshipClass
  ][
    order(-N)
  ]
)

cat("\n")
cat("============================================================\n")
cat("TARGET PRESENCE\n")
cat("============================================================\n")

print(
  audit[
    ,
    .N,
    by = TargetPresentInMatrix
  ]
)

cat("\n")
cat("============================================================\n")
cat("MULTI-OLD-ID TARGETS\n")
cat("============================================================\n")

multi <- group_summary[
  NumOldGeneIDs > 1
]

cat(
  "Targets receiving >1 old GeneID:",
  nrow(multi),
  "\n"
)

print(
  head(
    multi[
      order(
        -NumOldGeneIDs,
        TargetGeneID
      )
    ],
    30
  )
)


cat("\n")
cat("============================================================\n")
cat("CORRELATION SUMMARY FOR TARGETS ALREADY PRESENT\n")
cat("============================================================\n")

present <- audit[
  TargetPresentInMatrix == TRUE &
  !is.na(Pearson)
]

if (nrow(present) > 0) {

  cat(
    "N comparisons       :",
    nrow(present),
    "\n"
  )

  cat(
    "Pearson min         :",
    min(present$Pearson),
    "\n"
  )

  cat(
    "Pearson Q1          :",
    quantile(
      present$Pearson,
      0.25
    ),
    "\n"
  )

  cat(
    "Pearson median      :",
    median(present$Pearson),
    "\n"
  )

  cat(
    "Pearson Q3          :",
    quantile(
      present$Pearson,
      0.75
    ),
    "\n"
  )

  cat(
    "Pearson max         :",
    max(present$Pearson),
    "\n"
  )
}


cat("\nIMPORTANT:\n")
cat(
  "- This audit does NOT decide collapse strategy.\n",
  "- Exact or correlated profiles are descriptive only.\n",
  "- No expression matrix was modified.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_REDIRECT_EXPRESSION_COLLISION_AUDIT_COMPLETED\n")
cat("============================================================\n")
