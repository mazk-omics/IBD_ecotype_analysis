suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

root <- "/home/mazekai/IBD_EcoTyper/06_validation/02_bulk_treatment/GSE72819"

soft_file <- file.path(
  root,
  "00_provenance/geo_downloads/GSE72819_family.soft.gz"
)

matrix_file <- file.path(
  root,
  "00_provenance/geo_downloads/GSE72819_series_matrix.txt.gz"
)

meta_file <- file.path(
  root,
  "01_source_metadata/original/metadata_geo_GSE72819.xlsx"
)

raw_dir <- file.path(
  root,
  "02_source_expression/raw_tar_contents"
)

out_dir <- file.path(
  root,
  "00_provenance"
)

cat("============================================================\n")
cat("GSE72819 Stage 0 provenance audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Local curated metadata
# ============================================================

meta <- as.data.table(
  read_excel(
    meta_file,
    sheet = "GEO数据整理"
  )
)

gsm_meta <- sort(
  unique(
    as.character(meta$GSM编号)
  )
)

stopifnot(
  length(gsm_meta) == 73
)

cat("[PASS] Local metadata GSM:", length(gsm_meta), "\n")


# ============================================================
# 2. SOFT GSM
# ============================================================

soft_lines <- readLines(
  gzfile(soft_file),
  warn = FALSE
)

soft_sample_lines <- grep(
  "^\\^SAMPLE = GSM[0-9]+",
  soft_lines,
  value = TRUE
)

gsm_soft <- sub(
  "^\\^SAMPLE = ",
  "",
  soft_sample_lines
)

gsm_soft <- sort(
  unique(gsm_soft)
)

cat("[PASS] SOFT GSM:", length(gsm_soft), "\n")


# ============================================================
# 3. Series Matrix GSM
# ============================================================

matrix_lines <- readLines(
  gzfile(matrix_file),
  warn = FALSE
)

geo_line <- grep(
  "^!Sample_geo_accession",
  matrix_lines,
  value = TRUE
)

if (length(geo_line) != 1) {
  stop(
    "Expected exactly one !Sample_geo_accession line in series matrix"
  )
}

gsm_matrix <- regmatches(
  geo_line,
  gregexpr(
    "GSM[0-9]+",
    geo_line
  )
)[[1]]

gsm_matrix <- sort(
  unique(gsm_matrix)
)

cat("[PASS] Series Matrix GSM:", length(gsm_matrix), "\n")


# ============================================================
# 4. RAW processed files
# ============================================================

files <- list.files(
  raw_dir,
  pattern = "\\.counts_gene_exonic\\.tab\\.gz$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop(
    "No extracted processed files found in: ",
    raw_dir
  )
}

fn <- basename(files)

gsm_raw <- sub(
  "^(GSM[0-9]+)_.*$",
  "\\1",
  fn
)

stopifnot(
  all(
    grepl(
      "^GSM[0-9]+$",
      gsm_raw
    )
  )
)

cat("[PASS] RAW processed files:", length(files), "\n")
cat("[PASS] Unique RAW GSM:", uniqueN(gsm_raw), "\n")


# ============================================================
# 5. Technical fields encoded in filenames
# ============================================================

extract_group <- function(x, pattern) {

  m <- regexec(
    pattern,
    x
  )

  z <- regmatches(
    x,
    m
  )

  vapply(
    z,
    function(y) {
      if (length(y) >= 2) y[2] else NA_character_
    },
    character(1)
  )
}


raw_manifest <- data.table(
  GSM = gsm_raw,
  File = fn,

  Run = extract_group(
    fn,
    "_(R[0-9]+)_LIB"
  ),

  Library = extract_group(
    fn,
    "_(LIB[0-9]+)_SAM"
  ),

  SAM = extract_group(
    fn,
    "_(SAM[0-9]+)_L[0-9]+"
  ),

  Lane = extract_group(
    fn,
    "_(L[0-9]+)_NXG"
  ),

  NXG = extract_group(
    fn,
    "_(NXG[0-9]+)\\."
  )
)


stopifnot(
  nrow(raw_manifest) == 73,
  uniqueN(raw_manifest$GSM) == 73
)


cat("\n============================================================\n")
cat("TECHNICAL FILE-NAME STRUCTURE\n")
cat("============================================================\n")

cat("\nRun:\n")
print(
  raw_manifest[
    ,
    .N,
    by = Run
  ][order(Run)]
)

cat("\nLane:\n")
print(
  raw_manifest[
    ,
    .N,
    by = Lane
  ][order(Lane)]
)

cat("\nRun x Lane:\n")
print(
  dcast(
    raw_manifest,
    Run ~ Lane,
    value.var = "GSM",
    fun.aggregate = length
  )
)


# ============================================================
# 6. Four-way GSM concordance
# ============================================================

all_gsm <- sort(
  unique(
    c(
      gsm_meta,
      gsm_soft,
      gsm_matrix,
      gsm_raw
    )
  )
)

crosswalk <- data.table(
  GSM = all_gsm
)

crosswalk[
  ,
  In_LocalMetadata :=
    GSM %chin% gsm_meta
]

crosswalk[
  ,
  In_SOFT :=
    GSM %chin% gsm_soft
]

crosswalk[
  ,
  In_SeriesMatrix :=
    GSM %chin% gsm_matrix
]

crosswalk[
  ,
  In_RAW :=
    GSM %chin% gsm_raw
]

crosswalk[
  ,
  AllFour :=
    In_LocalMetadata &
    In_SOFT &
    In_SeriesMatrix &
    In_RAW
]


cat("\n============================================================\n")
cat("FOUR-WAY GSM CONCORDANCE\n")
cat("============================================================\n")

summary <- data.table(
  Source = c(
    "Local metadata",
    "SOFT",
    "Series Matrix",
    "RAW processed files",
    "Union",
    "Present in all four"
  ),

  N = c(
    length(gsm_meta),
    length(gsm_soft),
    length(gsm_matrix),
    uniqueN(gsm_raw),
    length(all_gsm),
    sum(crosswalk$AllFour)
  )
)

print(summary)


discordant <- crosswalk[
  AllFour == FALSE
]

cat(
  "\nDiscordant GSM:",
  nrow(discordant),
  "\n"
)

if (nrow(discordant) > 0) {
  print(discordant)
}


# ============================================================
# 7. Processed file schema audit
#
# Read only header + first few rows here.
# Full expression-value audit comes in Stage 1.
# ============================================================

schema_list <- lapply(
  files,
  function(f) {

    con <- gzfile(
      f,
      open = "rt"
    )

    on.exit(
      close(con),
      add = TRUE
    )

    x <- read.delim(
      con,
      nrows = 5,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    data.table(
      GSM =
        sub(
          "^(GSM[0-9]+)_.*$",
          "\\1",
          basename(f)
        ),

      File =
        basename(f),

      NColumns =
        ncol(x),

      ColumnNames =
        paste(
          names(x),
          collapse = " | "
        )
    )
  }
)

schema <- rbindlist(
  schema_list
)


cat("\n============================================================\n")
cat("PROCESSED FILE SCHEMA\n")
cat("============================================================\n")

schema_summary <- schema[
  ,
  .N,
  by = .(
    NColumns,
    ColumnNames
  )
][order(-N)]

print(schema_summary)


# ============================================================
# 8. Series Matrix: expression present or metadata only?
# ============================================================

table_begin <- grep(
  "^!series_matrix_table_begin",
  matrix_lines
)

table_end <- grep(
  "^!series_matrix_table_end",
  matrix_lines
)

matrix_data_rows <- 0L

if (
  length(table_begin) == 1 &&
  length(table_end) == 1 &&
  table_end > table_begin
) {
  matrix_data_rows <-
    table_end -
    table_begin -
    1L
}


cat("\n============================================================\n")
cat("SERIES MATRIX CONTENT\n")
cat("============================================================\n")

cat(
  "Rows between table_begin and table_end:",
  matrix_data_rows,
  "\n"
)


# ============================================================
# 9. Save outputs
# ============================================================

fwrite(
  raw_manifest,
  file.path(
    out_dir,
    "GSE72819_RAW_processed_file_manifest.tsv"
  ),
  sep = "\t"
)

fwrite(
  crosswalk,
  file.path(
    out_dir,
    "GSE72819_four_way_GSM_crosswalk.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary,
  file.path(
    out_dir,
    "GSE72819_four_way_GSM_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  schema,
  file.path(
    out_dir,
    "GSE72819_processed_file_schema.tsv"
  ),
  sep = "\t"
)

fwrite(
  schema_summary,
  file.path(
    out_dir,
    "GSE72819_processed_file_schema_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 10. Final status
# ============================================================

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n")

if (
  length(gsm_meta) == 73 &&
  length(gsm_soft) == 73 &&
  length(gsm_matrix) == 73 &&
  uniqueN(gsm_raw) == 73 &&
  length(all_gsm) == 73 &&
  sum(crosswalk$AllFour) == 73 &&
  nrow(schema_summary) == 1
) {

  cat("[PASS] 73 GSM in local metadata\n")
  cat("[PASS] 73 GSM in SOFT\n")
  cat("[PASS] 73 GSM in Series Matrix\n")
  cat("[PASS] 73 GSM in processed RAW files\n")
  cat("[PASS] All four sources contain exactly the same GSM set\n")
  cat("[PASS] All 73 processed files share one schema\n")

  cat("\nFINAL STATUS:\n")
  cat("PASS_GSE72819_STAGE0_PROVENANCE_AUDIT\n")

} else {

  cat("[REVIEW] One or more provenance checks require review\n")

  cat("\nFINAL STATUS:\n")
  cat("REVIEW_GSE72819_STAGE0_PROVENANCE_AUDIT\n")
}

cat("============================================================\n")
