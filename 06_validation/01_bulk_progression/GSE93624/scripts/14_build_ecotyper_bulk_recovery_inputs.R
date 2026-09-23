suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# EcoTyper bulk recovery input preparation
#
# Datasets:
#   GSE57945
#   GSE93624
#
# PRINCIPLES
#   - use frozen canonical matrices
#   - NO new expression transformation
#   - NO z-scoring
#   - NO rank transformation
#   - NO cross-cohort harmonization
#   - NO batch correction
#   - first column = GENES
#   - remaining columns = sample IDs
#   - annotation first column = ID
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression"

base57945 <- file.path(root, "GSE57945")
base93624 <- file.path(root, "GSE93624")

out57945 <- file.path(
  base57945,
  "prepared/06_ecotyper_recovery_input"
)

out93624 <- file.path(
  base93624,
  "prepared/06_ecotyper_recovery_input"
)

dir.create(out57945, recursive = TRUE, showWarnings = FALSE)
dir.create(out93624, recursive = TRUE, showWarnings = FALSE)


cat("============================================================\n")
cat("EcoTyper bulk recovery input preparation\n")
cat("============================================================\n\n")


# ============================================================
# Helper: robust expression loader
# ============================================================

get_expression <- function(obj) {

  if (is.matrix(obj)) {
    return(obj)
  }

  if (
    is.list(obj) &&
    !is.null(obj$expression)
  ) {
    return(
      as.matrix(obj$expression)
    )
  }

  stop(
    "Object is neither an expression matrix nor a list ",
    "containing $expression"
  )
}


# ============================================================
# Helper: expression matrix validation
# ============================================================

validate_expression <- function(x, dataset) {

  stopifnot(
    is.matrix(x),
    nrow(x) > 0,
    ncol(x) > 0,
    !is.null(rownames(x)),
    !is.null(colnames(x)),
    !anyDuplicated(rownames(x)),
    !anyDuplicated(colnames(x)),
    !anyNA(x),
    all(is.finite(x))
  )

  # EcoTyper uses read.delim; sample IDs should survive
  # make.names unchanged whenever possible.
  changed <- make.names(
    colnames(x)
  ) != colnames(x)

  if (any(changed)) {

    stop(
      dataset,
      ": some sample IDs would be modified by make.names(): ",
      paste(
        colnames(x)[changed],
        collapse = ", "
      )
    )
  }

  cat(
    "[PASS]",
    dataset,
    "expression validation:",
    nrow(x),
    "x",
    ncol(x),
    "\n"
  )
}


# ============================================================
# 1. GSE57945
# ============================================================

cat("\n============================================================\n")
cat("GSE57945\n")
cat("============================================================\n")

f57945 <- file.path(
  base57945,
  "prepared/03_canonical_expression",
  "GSE57945_canonical_expression_34368x322.rds"
)

obj57945 <- readRDS(f57945)
x57945 <- get_expression(obj57945)

stopifnot(
  nrow(x57945) == 34368,
  ncol(x57945) == 322
)

validate_expression(
  x57945,
  "GSE57945"
)


# ------------------------------------------------------------
# Export expression exactly as frozen.
# ------------------------------------------------------------

export57945 <- data.table(
  GENES = rownames(x57945)
)

export57945 <- cbind(
  export57945,
  as.data.table(x57945)
)

expr57945_file <- file.path(
  out57945,
  "GSE57945_EcoTyper_recovery_expression.txt"
)

fwrite(
  export57945,
  expr57945_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# Minimal annotation.
#
# Clinical metadata is NOT required for state/ecotype recovery.
# We deliberately keep this input minimal and attach clinical
# information to recovery outputs downstream.
# ------------------------------------------------------------

ann57945 <- data.table(
  ID = colnames(x57945)
)

ann57945_file <- file.path(
  out57945,
  "GSE57945_EcoTyper_recovery_annotation.txt"
)

fwrite(
  ann57945,
  ann57945_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 2. GSE93624
# ============================================================

cat("\n============================================================\n")
cat("GSE93624\n")
cat("============================================================\n")

f93624 <- file.path(
  base93624,
  "prepared/03_canonical_expression",
  "GSE93624_canonical_expression_13323x245.rds"
)

meta93624_file <- file.path(
  base93624,
  "prepared/04_sample_metadata",
  "GSE93624_canonical_sample_metadata.tsv"
)

obj93624 <- readRDS(f93624)
x93624 <- get_expression(obj93624)

stopifnot(
  nrow(x93624) == 13323,
  ncol(x93624) == 245
)

validate_expression(
  x93624,
  "GSE93624"
)


meta93624 <- fread(
  meta93624_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  nrow(meta93624) == 245,
  identical(
    meta93624$ExpressionSampleID,
    colnames(x93624)
  )
)


# ------------------------------------------------------------
# Export expression EXACTLY AS FROZEN.
#
# Important:
# Negative values are retained.
# We do not exponentiate / shift / renormalize.
# ------------------------------------------------------------

export93624 <- data.table(
  GENES = rownames(x93624)
)

export93624 <- cbind(
  export93624,
  as.data.table(x93624)
)

expr93624_file <- file.path(
  out93624,
  "GSE93624_EcoTyper_recovery_expression.txt"
)

fwrite(
  export93624,
  expr93624_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# Recovery annotation.
#
# ID must exactly equal expression sample column names.
# Progression is intentionally included because this cohort
# will later be used for clinical progression analyses.
# ------------------------------------------------------------

ann93624 <- meta93624[
  ,
  .(
    ID =
      ExpressionSampleID,

    GSM,
    Diagnosis,
    Progression3yr,
    Sex,
    Tissue,
    ParisAge,
    Ancestry
  )
]


stopifnot(
  identical(
    ann93624$ID,
    colnames(x93624)
  )
)


ann93624_file <- file.path(
  out93624,
  "GSE93624_EcoTyper_recovery_annotation.txt"
)

fwrite(
  ann93624,
  ann93624_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 3. Reload exported expression files
#
# Check that text export has not altered dimensions,
# row identities, sample identities, or numeric values.
# ============================================================

reload_expr <- function(path) {

  y <- fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )

  stopifnot(
    colnames(y)[1] == "GENES"
  )

  genes <- y[[1]]

  m <- as.matrix(
    y[, -1, drop = FALSE]
  )

  storage.mode(m) <- "double"

  rownames(m) <- genes

  m
}


r57945 <- reload_expr(
  expr57945_file
)

r93624 <- reload_expr(
  expr93624_file
)


stopifnot(
  identical(
    dim(r57945),
    dim(x57945)
  ),
  identical(
    rownames(r57945),
    rownames(x57945)
  ),
  identical(
    colnames(r57945),
    colnames(x57945)
  ),
  isTRUE(
    all.equal(
      unname(r57945),
      unname(x57945),
      tolerance = 0,
      check.attributes = FALSE
    )
  )
)


stopifnot(
  identical(
    dim(r93624),
    dim(x93624)
  ),
  identical(
    rownames(r93624),
    rownames(x93624)
  ),
  identical(
    colnames(r93624),
    colnames(x93624)
  ),
  isTRUE(
    all.equal(
      unname(r93624),
      unname(x93624),
      tolerance = 0,
      check.attributes = FALSE
    )
  )
)


cat(
  "\n[PASS] Text exports reproduce frozen canonical matrices exactly\n"
)


# ============================================================
# 4. EcoTyper automatic preprocessing audit
#
# state_recovery_bulk.R does:
#
#   if(all(data >= 0) && max(data) > 50)
#       data = log2(data + 1)
#
# Then it gene-wise scales across samples.
#
# We only document what EcoTyper itself will do.
# ============================================================

auto_log57945 <- (
  all(x57945 >= 0) &&
  max(x57945) > 50
)

auto_log93624 <- (
  all(x93624 >= 0) &&
  max(x93624) > 50
)


preprocess_audit <- data.table(
  Dataset = c(
    "GSE57945",
    "GSE93624"
  ),

  Genes = c(
    nrow(x57945),
    nrow(x93624)
  ),

  Samples = c(
    ncol(x57945),
    ncol(x93624)
  ),

  Min = c(
    min(x57945),
    min(x93624)
  ),

  Max = c(
    max(x57945),
    max(x93624)
  ),

  NegativeValues = c(
    sum(x57945 < 0),
    sum(x93624 < 0)
  ),

  InputScale = c(
    "RPKM",
    paste0(
      "GEO processed continuous normalized expression; ",
      "exact deposited transformation unresolved"
    )
  ),

  ExpressionModifiedDuringInputBuild = c(
    FALSE,
    FALSE
  ),

  EcoTyperWillAutoLog2Plus1 = c(
    auto_log57945,
    auto_log93624
  ),

  EcoTyperWillGeneWiseScale = c(
    TRUE,
    TRUE
  )
)


cat("\n============================================================\n")
cat("ECOTYPER PREPROCESSING AUDIT\n")
cat("============================================================\n")

print(
  preprocess_audit
)


# ============================================================
# 5. Annotation audits
# ============================================================

ann_audit <- data.table(
  Dataset = c(
    "GSE57945",
    "GSE93624"
  ),

  ExpressionSamples = c(
    ncol(x57945),
    ncol(x93624)
  ),

  AnnotationRows = c(
    nrow(ann57945),
    nrow(ann93624)
  ),

  UniqueIDs = c(
    uniqueN(ann57945$ID),
    uniqueN(ann93624$ID)
  ),

  ExactOrderMatch = c(
    identical(
      ann57945$ID,
      colnames(x57945)
    ),
    identical(
      ann93624$ID,
      colnames(x93624)
    )
  )
)


cat("\n============================================================\n")
cat("ANNOTATION AUDIT\n")
cat("============================================================\n")

print(
  ann_audit
)


stopifnot(
  all(
    ann_audit$ExpressionSamples ==
      ann_audit$AnnotationRows
  ),
  all(
    ann_audit$ExpressionSamples ==
      ann_audit$UniqueIDs
  ),
  all(
    ann_audit$ExactOrderMatch
  )
)


# ============================================================
# 6. Build final input manifest
# ============================================================

manifest <- data.table(
  Dataset = c(
    "GSE57945",
    "GSE93624"
  ),

  Role = c(
    "Recovery robustness / reproducibility",
    "Clinical progression recovery / association"
  ),

  ExpressionFile = c(
    expr57945_file,
    expr93624_file
  ),

  AnnotationFile = c(
    ann57945_file,
    ann93624_file
  ),

  GeneRows = c(
    34368,
    13323
  ),

  Samples = c(
    322,
    245
  ),

  InputExpressionChanged = c(
    FALSE,
    FALSE
  ),

  ReadyForEcoTyperInput = c(
    TRUE,
    TRUE
  )
)


fwrite(
  preprocess_audit,
  file.path(
    out93624,
    "EcoTyper_recovery_input_preprocessing_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  ann_audit,
  file.path(
    out93624,
    "EcoTyper_recovery_input_annotation_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  manifest,
  file.path(
    out93624,
    "EcoTyper_bulk_recovery_input_manifest.tsv"
  ),
  sep = "\t"
)


# Copy manifest into GSE57945 side as well.
fwrite(
  manifest,
  file.path(
    out57945,
    "EcoTyper_bulk_recovery_input_manifest.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 7. MD5 provenance
# ============================================================

files <- c(
  expr57945_file,
  ann57945_file,
  expr93624_file,
  ann93624_file
)

hash <- tools::md5sum(files)

hash_table <- data.table(
  File = names(hash),
  MD5 = unname(hash)
)

fwrite(
  hash_table,
  file.path(
    out93624,
    "EcoTyper_bulk_recovery_input_md5.tsv"
  ),
  sep = "\t"
)

fwrite(
  hash_table,
  file.path(
    out57945,
    "EcoTyper_bulk_recovery_input_md5.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 8. Final
# ============================================================

cat("\n============================================================\n")
cat("FINAL INPUT CHECKS\n")
cat("============================================================\n")

cat("[PASS] GSE57945 recovery expression exported unchanged\n")
cat("[PASS] GSE57945 annotation IDs exactly match expression\n")

cat("[PASS] GSE93624 recovery expression exported unchanged\n")
cat("[PASS] GSE93624 annotation IDs exactly match expression\n")

cat("[PASS] Gene symbols unique in both datasets\n")
cat("[PASS] Sample IDs unique in both datasets\n")
cat("[PASS] Sample IDs are make.names-safe\n")
cat("[PASS] No cross-cohort normalization performed\n")
cat("[PASS] No expression transformation performed during input build\n")

cat("\nFINAL STATUS:\n")
cat("PASS_BULK_ECOTYPER_RECOVERY_INPUTS_BUILT\n")
cat("============================================================\n")
