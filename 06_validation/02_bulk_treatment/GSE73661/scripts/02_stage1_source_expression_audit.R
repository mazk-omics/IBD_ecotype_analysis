suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# GSE73661 Stage 1 v2
# Authoritative processed-expression audit
#
# Frozen facts already established:
#   GEO platform: GPL6244
#   Series Matrix: 33,252 probe sets x 178 samples
#   GEO processing: RMA via aroma.affymetrix
#   Output scale: log2 expression value per probe set
#
# This stage:
#   - validates Series Matrix
#   - validates all SOFT sample tables
#   - performs value-level SOFT vs Series Matrix checks
#   - freezes author-deposited processed expression
#
# NO CEL reprocessing
# NO normalization
# NO transformation
# NO probe-to-gene mapping
# ============================================================

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE73661"

matrix_file <- file.path(
  root,
  "00_provenance/geo_downloads",
  "GSE73661_series_matrix.txt.gz"
)

soft_file <- file.path(
  root,
  "00_provenance/geo_downloads",
  "GSE73661_family.soft.gz"
)

raw_file <- file.path(
  root,
  "00_provenance/geo_downloads",
  "GSE73661_RAW.tar"
)

cohort_file <- file.path(
  root,
  "prepared/00_cohort_reconstruction",
  "GSE73661_cohort_reconstruction_178_samples.tsv"
)

out_dir <- file.path(
  root,
  "prepared/01_authoritative_matrix"
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


same_chr <- function(x, y) {

  x <- as.character(x)
  y <- as.character(y)

  if (length(x) != length(y)) {
    return(FALSE)
  }

  all(
    !is.na(x) &
      !is.na(y) &
      x == y
  )
}


shell_count <- function(pattern, file) {

  cmd <- paste(
    "zgrep -Eic",
    shQuote(pattern),
    shQuote(file)
  )

  out <- system(
    cmd,
    intern = TRUE
  )

  suppressWarnings(
    as.integer(out[1])
  )
}


cat("============================================================\n")
cat("GSE73661 Stage 1 v2: authoritative processed-expression audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Input existence
# ============================================================

check(
  file.exists(matrix_file),
  "GEO Series Matrix exists"
)

check(
  file.exists(soft_file),
  "GEO SOFT exists"
)

check(
  file.exists(raw_file),
  "RAW CEL archive exists"
)

check(
  file.exists(cohort_file),
  "Stage-0B cohort reconstruction exists"
)


# ============================================================
# 2. Extract Series Matrix expression block
# ============================================================

matrix_awk <- tempfile(fileext = ".awk")
matrix_tmp <- tempfile(fileext = ".tsv")

writeLines(
  c(
    'tolower($0) ~ /^!series_matrix_table_begin/ { in_table=1; next }',
    'tolower($0) ~ /^!series_matrix_table_end/   { in_table=0; next }',
    'in_table { print }'
  ),
  matrix_awk
)

status <- system(
  paste(
    "zcat",
    shQuote(matrix_file),
    "| awk -f",
    shQuote(matrix_awk),
    ">",
    shQuote(matrix_tmp)
  )
)

check(
  status == 0,
  "Series Matrix expression block extracted"
)

check(
  file.exists(matrix_tmp) &&
    file.info(matrix_tmp)$size > 0,
  "Extracted Series Matrix table is non-empty"
)

series <- fread(
  matrix_tmp,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  na.strings = c(
    "null",
    "NULL",
    "NA"
  )
)

unlink(
  c(
    matrix_awk,
    matrix_tmp
  )
)


check(
  ncol(series) == 179,
  "Series Matrix contains ID_REF + 178 samples"
)

check(
  names(series)[1] %in%
    c(
      "ID_REF",
      "ID"
    ),
  "Series Matrix first column is probe-set identifier"
)


probe_ids <- as.character(
  series[[1]]
)

sample_ids <- as.character(
  names(series)[-1]
)


check(
  nrow(series) == 33252,
  "Series Matrix contains 33,252 probe sets"
)

check(
  length(sample_ids) == 178,
  "Series Matrix contains 178 samples"
)

check(
  uniqueN(probe_ids) == 33252,
  "All probe-set IDs are unique"
)

check(
  uniqueN(sample_ids) == 178,
  "All sample IDs are unique"
)


# ============================================================
# 3. Numeric expression matrix
# ============================================================

value_dt <- series[
  ,
  -1,
  with = FALSE
]

for (nm in names(value_dt)) {

  if (!is.numeric(value_dt[[nm]])) {

    old <- value_dt[[nm]]

    new <- suppressWarnings(
      as.numeric(old)
    )

    invalid <- is.na(new) &
      !is.na(old)

    check(
      !any(invalid),
      paste0(
        nm,
        ": values convert cleanly to numeric"
      )
    )

    set(
      value_dt,
      j = nm,
      value = new
    )
  }
}


expr <- as.matrix(
  value_dt
)

storage.mode(expr) <- "double"

rownames(expr) <- probe_ids
colnames(expr) <- sample_ids


check(
  identical(
    dim(expr),
    c(
      33252L,
      178L
    )
  ),
  "Authoritative expression dimensions = 33,252 x 178"
)

check(
  !anyNA(expr),
  "Expression matrix contains no NA"
)

check(
  all(is.finite(expr)),
  "Expression matrix contains only finite values"
)


# ============================================================
# 4. Cohort alignment
# ============================================================

cohort <- fread(
  cohort_file
)

check(
  nrow(cohort) == 178,
  "Cohort reconstruction contains 178 samples"
)

check(
  uniqueN(cohort$GSM) == 178,
  "Cohort reconstruction contains 178 unique GSM IDs"
)

check(
  setequal(
    sample_ids,
    cohort$GSM
  ),
  "Series Matrix and Stage-0B contain identical GSM universe"
)

cohort <- cohort[
  match(
    sample_ids,
    GSM
  )
]

check(
  same_chr(
    cohort$GSM,
    sample_ids
  ),
  "Cohort reconstruction aligned to Series Matrix order"
)


# ============================================================
# 5. GEO processing metadata
# ============================================================

processing_lines <- system(
  paste(
    "zgrep",
    shQuote("^!Sample_data_processing ="),
    shQuote(soft_file),
    "| sort -u"
  ),
  intern = TRUE
)

check(
  length(processing_lines) > 0,
  "GEO processing metadata recovered"
)

rma_declared <- any(
  grepl(
    "RMA|robust multichip average|robust multi-array average",
    processing_lines,
    ignore.case = TRUE
  )
)

log2_declared <- any(
  grepl(
    "log2",
    processing_lines,
    ignore.case = TRUE
  )
)


cat("\n============================================================\n")
cat("GEO PROCESSING METADATA\n")
cat("============================================================\n")

cat(
  paste(
    processing_lines,
    collapse = "\n"
  ),
  "\n"
)

cat(
  "\nRMA declared: ",
  rma_declared,
  "\n",
  sep = ""
)

cat(
  "Log2 explicitly declared: ",
  log2_declared,
  "\n",
  sep = ""
)


check(
  rma_declared,
  "GEO explicitly declares RMA processing"
)

check(
  log2_declared,
  "GEO explicitly declares log2 expression values"
)


writeLines(
  processing_lines,
  file.path(
    out_dir,
    "GSE73661_GEO_data_processing_metadata.txt"
  )
)


# ============================================================
# 6. Count SOFT sample-table markers BEFORE parsing
# ============================================================

begin_n <- shell_count(
  "^!sample_table_begin[[:space:]]*$",
  soft_file
)

end_n <- shell_count(
  "^!sample_table_end[[:space:]]*$",
  soft_file
)


cat("\n============================================================\n")
cat("SOFT SAMPLE-TABLE MARKER AUDIT\n")
cat("============================================================\n")

cat(
  "sample_table_begin markers: ",
  begin_n,
  "\n",
  sep = ""
)

cat(
  "sample_table_end markers:   ",
  end_n,
  "\n",
  sep = ""
)


check(
  !is.na(begin_n),
  "sample_table_begin marker count recovered"
)

check(
  !is.na(end_n),
  "sample_table_end marker count recovered"
)

check(
  begin_n == 178,
  "Exactly 178 sample_table_begin markers found"
)

check(
  end_n == 178,
  "Exactly 178 sample_table_end markers found"
)


# ============================================================
# 7. Full SOFT table audit
#
# Header is split into two explicit fields to avoid embedding
# a TAB inside a TSV column.
# ============================================================

audit_awk <- tempfile(fileext = ".awk")
audit_tmp <- tempfile(fileext = ".tsv")

writeLines(
  c(
    'BEGIN {',
    '  OFS="\\t"',
    '  print "GSM","TableRows","HeaderCol1","HeaderCol2"',
    '}',
    '',
    '/^\\^SAMPLE = / {',
    '  gsm=$0',
    '  sub(/^\\^SAMPLE = /,"",gsm)',
    '  next',
    '}',
    '',
    'tolower($0) ~ /^!sample_table_begin/ {',
    '  in_table=1',
    '  header_pending=1',
    '  rows=0',
    '  h1=""',
    '  h2=""',
    '  next',
    '}',
    '',
    'in_table && header_pending {',
    '  split($0,a,"\\t")',
    '  h1=a[1]',
    '  h2=a[2]',
    '  gsub(/^"|"$/,"",h1)',
    '  gsub(/^"|"$/,"",h2)',
    '  header_pending=0',
    '  next',
    '}',
    '',
    'tolower($0) ~ /^!sample_table_end/ {',
    '  if (in_table) {',
    '    print gsm,rows,h1,h2',
    '  }',
    '  in_table=0',
    '  header_pending=0',
    '  next',
    '}',
    '',
    'in_table {',
    '  rows++',
    '}'
  ),
  audit_awk
)


status <- system(
  paste(
    "zcat",
    shQuote(soft_file),
    "| awk -f",
    shQuote(audit_awk),
    ">",
    shQuote(audit_tmp)
  )
)

check(
  status == 0,
  "SOFT processed tables scanned"
)


soft_table_audit <- fread(
  audit_tmp,
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = FALSE
)

unlink(
  c(
    audit_awk,
    audit_tmp
  )
)


cat("\n============================================================\n")
cat("SOFT SAMPLE-TABLE SCHEMA SUMMARY\n")
cat("============================================================\n")

print(
  soft_table_audit[
    ,
    .N,
    by = .(
      TableRows,
      HeaderCol1,
      HeaderCol2
    )
  ][
    order(
      TableRows,
      HeaderCol1,
      HeaderCol2
    )
  ]
)


check(
  nrow(soft_table_audit) == 178,
  "Exactly 178 SOFT processed-data tables parsed"
)

check(
  uniqueN(soft_table_audit$GSM) == 178,
  "All 178 SOFT table GSM IDs are unique"
)

check(
  setequal(
    soft_table_audit$GSM,
    sample_ids
  ),
  "SOFT tables and Series Matrix contain identical GSM universe"
)

check(
  all(
    soft_table_audit$TableRows ==
      33252
  ),
  "Every SOFT sample table contains 33,252 rows"
)

check(
  all(
    soft_table_audit$HeaderCol1 ==
      "ID_REF"
  ),
  "Every SOFT sample table has ID_REF as first column"
)

check(
  all(
    soft_table_audit$HeaderCol2 ==
      "VALUE"
  ),
  "Every SOFT sample table has VALUE as second column"
)


fwrite(
  soft_table_audit,
  file.path(
    out_dir,
    "GSE73661_SOFT_sample_table_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 8. Representative sample selection
# ============================================================

pick_first <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      nzchar(x)
  ]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  x[1]
}


selected <- unique(
  c(
    pick_first(
      cohort[
        Branch == "CONTROL",
        GSM
      ]
    ),

    pick_first(
      cohort[
        Branch == "IFX" &
          Timepoint == "W0",
        GSM
      ]
    ),

    pick_first(
      cohort[
        Branch == "IFX" &
          Timepoint == "W4_6",
        GSM
      ]
    ),

    pick_first(
      cohort[
        Branch == "VDZ_TRIAL" &
          Timepoint == "W0" &
          ExposureClass == "VDZ_exposed",
        GSM
      ]
    ),

    pick_first(
      cohort[
        Branch == "VDZ_TRIAL" &
          Timepoint == "W0" &
          ExposureClass == "Placebo_only",
        GSM
      ]
    ),

    pick_first(
      cohort[
        Branch == "VDZ_TRIAL" &
          Timepoint == "W6" &
          ExposureClass == "VDZ_exposed",
        GSM
      ]
    ),

    pick_first(
      cohort[
        Branch == "VDZ_TRIAL" &
          Timepoint == "W12" &
          ExposureClass == "VDZ_exposed",
        GSM
      ]
    ),

    pick_first(
      cohort[
        Branch == "VDZ_TRIAL" &
          Timepoint == "W52" &
          ExposureClass == "VDZ_exposed",
        GSM
      ]
    )
  )
)


selected <- selected[
  !is.na(selected)
]


if (length(selected) < 8L) {

  filler <- setdiff(
    sample_ids,
    selected
  )

  selected <- c(
    selected,
    head(
      filler,
      8L - length(selected)
    )
  )
}


check(
  length(selected) == 8,
  "Eight representative samples selected"
)


cat("\nRepresentative GSMs:\n")
print(selected)


# ============================================================
# 9. Extract representative SOFT tables
# ============================================================

selected_file <- tempfile(fileext = ".txt")
extract_awk <- tempfile(fileext = ".awk")
extract_tmp <- tempfile(fileext = ".tsv")

writeLines(
  selected,
  selected_file
)


writeLines(
  c(
    'BEGIN {',
    '  OFS="\\t"',
    '  while ((getline z < selected_file) > 0) {',
    '    gsub(/\\r/,"",z)',
    '    if (z != "") selected[z]=1',
    '  }',
    '  print "GSM","Row","ID_REF","VALUE"',
    '}',
    '',
    '/^\\^SAMPLE = / {',
    '  gsm=$0',
    '  sub(/^\\^SAMPLE = /,"",gsm)',
    '  next',
    '}',
    '',
    'tolower($0) ~ /^!sample_table_begin/ {',
    '  in_table=1',
    '  header_pending=1',
    '  row_number=0',
    '  active=(gsm in selected)',
    '  next',
    '}',
    '',
    'in_table && header_pending {',
    '  header_pending=0',
    '  next',
    '}',
    '',
    'tolower($0) ~ /^!sample_table_end/ {',
    '  in_table=0',
    '  header_pending=0',
    '  active=0',
    '  next',
    '}',
    '',
    'in_table && active {',
    '  row_number++',
    '  split($0,a,"\\t")',
    '  id=a[1]',
    '  val=a[2]',
    '  gsub(/^"|"$/,"",id)',
    '  gsub(/^"|"$/,"",val)',
    '  print gsm,row_number,id,val',
    '}'
  ),
  extract_awk
)


status <- system(
  paste(
    "zcat",
    shQuote(soft_file),
    "| awk",
  paste0("-v selected_file=", shQuote(selected_file)),
    "-f",
    shQuote(extract_awk),
    ">",
    shQuote(extract_tmp)
  )
)

check(
  status == 0,
  "Representative SOFT tables extracted"
)


soft_selected <- fread(
  extract_tmp,
  sep = "\t",
  header = TRUE,
  quote = "",
  na.strings = c(
    "NA",
    "null",
    "NULL"
  )
)


unlink(
  c(
    selected_file,
    extract_awk,
    extract_tmp
  )
)


check(
  nrow(soft_selected) ==
    8L * 33252L,
  "Representative SOFT extraction has expected row count"
)

check(
  uniqueN(soft_selected$GSM) == 8,
  "All eight representative GSM tables recovered"
)


soft_selected[
  ,
  VALUE_NUM :=
    suppressWarnings(
      as.numeric(VALUE)
    )
]


check(
  !anyNA(
    soft_selected$VALUE_NUM
  ),
  "Representative SOFT values are numeric"
)


# ============================================================
# 10. SOFT vs Series Matrix value audit
# ============================================================

comparison_list <- vector(
  "list",
  length(selected)
)


for (i in seq_along(selected)) {

  gsm <- selected[i]

  s <- soft_selected[
    GSM == gsm
  ]

  check(
    nrow(s) == 33252,
    paste0(
      gsm,
      ": 33,252 SOFT rows"
    )
  )

  probe_match <- same_chr(
    s$ID_REF,
    probe_ids
  )

  check(
    probe_match,
    paste0(
      gsm,
      ": probe order matches Series Matrix"
    )
  )

  matrix_values <- expr[
    ,
    gsm
  ]

  diff <- abs(
    s$VALUE_NUM -
      matrix_values
  )

  comparison_list[[i]] <- data.table(
    GSM =
      gsm,

    SOFTRows =
      nrow(s),

    ProbeOrderExact =
      probe_match,

    MaxAbsoluteDifference =
      max(diff),

    MeanAbsoluteDifference =
      mean(diff),

    ExactNumericEquality =
      all(
        s$VALUE_NUM ==
          matrix_values
      ),

    Within1e12 =
      all(
        diff <=
          1e-12
      )
  )
}


soft_vs_series <- rbindlist(
  comparison_list
)


cat("\n============================================================\n")
cat("SOFT vs SERIES MATRIX VALUE AUDIT\n")
cat("============================================================\n")

print(
  soft_vs_series
)


check(
  all(
    soft_vs_series$ProbeOrderExact
  ),
  "All representative probe orders match"
)

check(
  all(
    soft_vs_series$Within1e12
  ),
  "Representative SOFT values reproduce Series Matrix within 1e-12"
)


fwrite(
  soft_vs_series,
  file.path(
    out_dir,
    "GSE73661_SOFT_vs_series_matrix_value_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 11. Authoritative expression distribution
# ============================================================

expression_summary <- data.table(
  Metric = c(
    "ProbeSets",
    "Samples",
    "Values",
    "Min",
    "Q01",
    "Q25",
    "Median",
    "Mean",
    "Q75",
    "Q99",
    "Max",
    "Zero_N",
    "Negative_N",
    "NA_N",
    "Nonfinite_N",
    "UniqueProbeSets",
    "UniqueSamples",
    "RMA_declared",
    "Log2_declared"
  ),

  Value = c(
    nrow(expr),
    ncol(expr),
    length(expr),

    min(expr),

    quantile(
      expr,
      0.01,
      names = FALSE
    ),

    quantile(
      expr,
      0.25,
      names = FALSE
    ),

    median(expr),

    mean(expr),

    quantile(
      expr,
      0.75,
      names = FALSE
    ),

    quantile(
      expr,
      0.99,
      names = FALSE
    ),

    max(expr),

    sum(expr == 0),
    sum(expr < 0),
    sum(is.na(expr)),
    sum(!is.finite(expr)),
    uniqueN(rownames(expr)),
    uniqueN(colnames(expr)),
    rma_declared,
    log2_declared
  )
)


cat("\n============================================================\n")
cat("AUTHORITATIVE EXPRESSION SUMMARY\n")
cat("============================================================\n")

print(
  expression_summary
)


# ============================================================
# 12. Per-sample descriptive QC
# ============================================================

sample_qc <- data.table(
  GSM =
    sample_ids,

  Min =
    apply(
      expr,
      2,
      min
    ),

  Q25 =
    apply(
      expr,
      2,
      quantile,
      probs = 0.25,
      names = FALSE
    ),

  Median =
    apply(
      expr,
      2,
      median
    ),

  Mean =
    colMeans(expr),

  Q75 =
    apply(
      expr,
      2,
      quantile,
      probs = 0.75,
      names = FALSE
    ),

  Max =
    apply(
      expr,
      2,
      max
    ),

  SD =
    apply(
      expr,
      2,
      sd
    )
)


sample_qc <- cbind(
  cohort[
    ,
    .(
      GSM,
      Branch,
      PatientKey,
      PatientNumber,
      Timepoint,
      ExposureClass,
      ResponseAtVisit,
      TreatmentSchedule
    )
  ],

  sample_qc[
    ,
    -"GSM"
  ]
)


check(
  same_chr(
    sample_qc$GSM,
    sample_ids
  ),
  "Per-sample QC matches authoritative matrix order"
)


# ============================================================
# 13. Save authoritative matrix
# ============================================================

authoritative_file <- file.path(
  out_dir,
  "GSE73661_authoritative_log2_RMA_33252x178.rds"
)


saveRDS(
  expr,
  authoritative_file,
  compress = FALSE
)


check(
  file.exists(authoritative_file),
  "Authoritative log2-RMA RDS written"
)


reload <- readRDS(
  authoritative_file
)


check(
  identical(
    dim(reload),
    c(
      33252L,
      178L
    )
  ),
  "RDS dimensions round-trip exactly"
)

check(
  same_chr(
    rownames(reload),
    rownames(expr)
  ),
  "RDS probe order round-trip exactly"
)

check(
  same_chr(
    colnames(reload),
    colnames(expr)
  ),
  "RDS sample order round-trip exactly"
)

check(
  identical(
    as.vector(reload),
    as.vector(expr)
  ),
  "RDS numeric payload round-trip exactly"
)


# ============================================================
# 14. Source sample manifest
# ============================================================

sample_manifest <- copy(
  cohort
)

sample_manifest[
  ,
  MatrixColumn :=
    seq_len(.N)
]


setcolorder(
  sample_manifest,
  c(
    "GSM",
    "MatrixColumn",
    setdiff(
      names(sample_manifest),
      c(
        "GSM",
        "MatrixColumn"
      )
    )
  )
)


check(
  same_chr(
    sample_manifest$GSM,
    sample_ids
  ),
  "Sample manifest matches authoritative matrix order"
)


# ============================================================
# 15. Save artifacts
# ============================================================

fwrite(
  sample_manifest,
  file.path(
    out_dir,
    "GSE73661_source_sample_manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  expression_summary,
  file.path(
    out_dir,
    "GSE73661_authoritative_expression_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

fwrite(
  sample_qc,
  file.path(
    out_dir,
    "GSE73661_source_sample_expression_qc.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


selected_manifest <- cohort[
  GSM %in%
    selected
]

selected_manifest[
  ,
  SelectionOrder :=
    match(
      GSM,
      selected
    )
]

setorder(
  selected_manifest,
  SelectionOrder
)


fwrite(
  selected_manifest,
  file.path(
    out_dir,
    "GSE73661_SOFT_value_audit_selected_samples.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


# ============================================================
# 16. Stage summary
# ============================================================

stage_summary <- data.table(
  Metric = c(
    "Probe_sets",
    "Samples",
    "SOFT_begin_markers",
    "SOFT_end_markers",
    "SOFT_tables_parsed",
    "SOFT_tables_with_33252_rows",
    "SOFT_tables_ID_REF_VALUE_schema",
    "Representative_samples_checked",
    "Representative_samples_matching_series",
    "RMA_declared",
    "Log2_declared",
    "NA_values",
    "Nonfinite_values",
    "Authoritative_scale",
    "Authoritative_source"
  ),

  Value = c(
    "33252",
    "178",
    as.character(begin_n),
    as.character(end_n),
    as.character(
      nrow(
        soft_table_audit
      )
    ),
    as.character(
      sum(
        soft_table_audit$TableRows ==
          33252
      )
    ),
    as.character(
      sum(
        soft_table_audit$HeaderCol1 ==
          "ID_REF" &
        soft_table_audit$HeaderCol2 ==
          "VALUE"
      )
    ),
    as.character(
      nrow(
        soft_vs_series
      )
    ),
    as.character(
      sum(
        soft_vs_series$Within1e12
      )
    ),
    as.character(
      rma_declared
    ),
    as.character(
      log2_declared
    ),
    as.character(
      sum(
        is.na(expr)
      )
    ),
    as.character(
      sum(
        !is.finite(expr)
      )
    ),
    "author-provided log2 RMA",
    "GEO Series Matrix / deposited Sample tables"
  )
)


cat("\n============================================================\n")
cat("STAGE 1 SUMMARY\n")
cat("============================================================\n")

print(
  stage_summary
)


fwrite(
  stage_summary,
  file.path(
    out_dir,
    "GSE73661_stage1_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 17. Hashes
# ============================================================

hashes <- tools::md5sum(
  c(
    matrix_file,
    soft_file,
    raw_file,
    cohort_file,
    authoritative_file
  )
)


fwrite(
  data.table(
    File =
      names(hashes),

    MD5 =
      unname(hashes)
  ),
  file.path(
    out_dir,
    "GSE73661_stage1_artifact_md5.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 18. Final gate
# ============================================================

critical_pass <- all(
  c(
    nrow(expr) == 33252,

    ncol(expr) == 178,

    uniqueN(
      rownames(expr)
    ) == 33252,

    uniqueN(
      colnames(expr)
    ) == 178,

    begin_n == 178,

    end_n == 178,

    nrow(
      soft_table_audit
    ) == 178,

    all(
      soft_table_audit$TableRows ==
        33252
    ),

    all(
      soft_table_audit$HeaderCol1 ==
        "ID_REF"
    ),

    all(
      soft_table_audit$HeaderCol2 ==
        "VALUE"
    ),

    all(
      soft_vs_series$ProbeOrderExact
    ),

    all(
      soft_vs_series$Within1e12
    ),

    rma_declared,

    log2_declared,

    !anyNA(expr),

    all(
      is.finite(expr)
    ),

    identical(
      as.vector(reload),
      as.vector(expr)
    )
  )
)


cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (critical_pass) {

  cat("[PASS] 33,252 probe sets retained\n")
  cat("[PASS] 178 samples retained\n")
  cat("[PASS] All 178 SOFT Sample tables verified\n")
  cat("[PASS] Every SOFT Sample table contains 33,252 rows\n")
  cat("[PASS] Every SOFT Sample table uses ID_REF / VALUE schema\n")
  cat("[PASS] Representative SOFT values reproduce Series Matrix\n")
  cat("[PASS] GEO explicitly documents RMA processing\n")
  cat("[PASS] GEO explicitly documents log2 expression values\n")
  cat("[PASS] Authoritative numeric payload preserved exactly\n")
  cat("[PASS] No CEL reprocessing performed\n")
  cat("[PASS] No normalization performed\n")
  cat("[PASS] No inverse-log transformation performed\n")
  cat("[PASS] No probe-to-gene mapping performed\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE73661_STAGE1_AUTHORITATIVE_LOG2_RMA_EXPRESSION\n")

} else {

  cat("[REVIEW] Stage 1 structural checks failed\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE73661_STAGE1_AUTHORITATIVE_EXPRESSION\n")
}

cat("============================================================\n")
