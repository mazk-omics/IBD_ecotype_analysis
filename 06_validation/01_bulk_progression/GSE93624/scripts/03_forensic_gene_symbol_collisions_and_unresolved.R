suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ============================================================
# GSE93624
# Stage 2B
#
# 1. Audit expression profiles for alias collisions
# 2. Classify unresolved symbol structures
# 3. Identify conservative punctuation-only current-symbol
#    rescue candidates
#
# IMPORTANT:
#   - No expression row is renamed.
#   - No expression row is merged.
#   - No expression row is excluded.
# ============================================================

base <- "/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

expr_file <- file.path(
  base,
  "prepared/01_authoritative_matrix",
  "GSE93624_authoritative_expression_245.rds"
)

audit_file <- file.path(
  base,
  "prepared/02_gene_harmonization",
  "GSE93624_gene_symbol_harmonization_audit.tsv"
)

out_dir <- file.path(
  base,
  "prepared/02_gene_harmonization"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("GSE93624 collision + unresolved-symbol forensic audit\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load
# ============================================================

obj <- readRDS(expr_file)
expr <- obj$expression

a <- fread(
  audit_file,
  colClasses = "character",
  na.strings = c("", "NA")
)

stopifnot(
  nrow(expr) == 13769,
  ncol(expr) == 245,
  nrow(a) == 13769,
  identical(
    rownames(expr),
    a$GeneSymbol_raw
  )
)

raw_symbols <- rownames(expr)

cat("[PASS] Expression and annotation audit aligned\n\n")


# ============================================================
# 2. Helper for pairwise expression comparison
# ============================================================

compare_profiles <- function(symbol1, symbol2) {

  x <- as.numeric(
    expr[symbol1, ]
  )

  y <- as.numeric(
    expr[symbol2, ]
  )

  exact <- identical(
    x,
    y
  )

  max_abs_diff <- max(
    abs(x - y)
  )

  mean_abs_diff <- mean(
    abs(x - y)
  )

  rmse <- sqrt(
    mean(
      (x - y)^2
    )
  )

  pearson <- suppressWarnings(
    cor(
      x,
      y,
      method = "pearson"
    )
  )

  spearman <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )

  if (exact) {

    relation <- "EXACT_DUPLICATE"

  } else if (!is.finite(pearson)) {

    relation <- "CORRELATION_UNDEFINED"

  } else if (pearson >= 0.99) {

    relation <- "VERY_HIGH_CORRELATED_DISTINCT"

  } else if (pearson >= 0.95) {

    relation <- "HIGH_CORRELATED_DISTINCT"

  } else {

    relation <- "DISTINCT_PROFILE"
  }

  data.table(
    Symbol1 = symbol1,
    Symbol2 = symbol2,
    ExactDuplicate = exact,
    Pearson = pearson,
    Spearman = spearman,
    MeanAbsDiff = mean_abs_diff,
    RMSE = rmse,
    MaxAbsDiff = max_abs_diff,
    Relationship = relation
  )
}


# ============================================================
# 3. ALIAS_TARGET_ALREADY_PRESENT expression audit
# ============================================================

target_present <- a[
  PreliminaryIdentityStatus ==
    "ALIAS_TARGET_ALREADY_PRESENT"
]

collision_target_list <- list()

if (nrow(target_present) > 0) {

  for (i in seq_len(nrow(target_present))) {

    old_symbol <-
      target_present$GeneSymbol_raw[i]

    current_symbol <-
      target_present$CandidateCurrentSymbol[i]

    stopifnot(
      old_symbol %chin% raw_symbols,
      current_symbol %chin% raw_symbols
    )

    x <- compare_profiles(
      old_symbol,
      current_symbol
    )

    x[
      ,
      `:=`(
        CandidateEntrezID =
          target_present$CandidateEntrezID[i],

        CollisionType =
          "ALIAS_TARGET_ALREADY_PRESENT"
      )
    ]

    collision_target_list[[i]] <- x
  }
}

if (length(collision_target_list)) {

  target_collision_audit <-
    rbindlist(
      collision_target_list
    )

} else {

  target_collision_audit <- data.table()
}


# ============================================================
# 4. Many-to-one alias collision audit
#
# Specifically audit groups where multiple historical symbols
# map to the same candidate GeneID.
# ============================================================

candidate_collision_rows <- a[
  !is.na(CandidateEntrezID) &
    !is.na(CandidateGeneIDSourceN) &
    as.integer(CandidateGeneIDSourceN) > 1
]

collision_ids <- unique(
  candidate_collision_rows$CandidateEntrezID
)

many_to_one_list <- list()
counter <- 0L

for (gid in collision_ids) {

  symbols <- candidate_collision_rows[
    CandidateEntrezID == gid,
    GeneSymbol_raw
  ]

  symbols <- unique(symbols)

  if (length(symbols) < 2)
    next

  pairs <- combn(
    symbols,
    2,
    simplify = FALSE
  )

  for (pair in pairs) {

    counter <- counter + 1L

    x <- compare_profiles(
      pair[1],
      pair[2]
    )

    x[
      ,
      `:=`(
        CandidateEntrezID = gid,
        CollisionType =
          "MULTIPLE_SOURCE_SYMBOLS_TO_ONE_CURRENT_GENE"
      )
    ]

    many_to_one_list[[counter]] <- x
  }
}

if (length(many_to_one_list)) {

  many_to_one_audit <-
    rbindlist(
      many_to_one_list
    )

} else {

  many_to_one_audit <- data.table()
}


# ============================================================
# 5. Direct-current ambiguous feature
# ============================================================

direct_ambiguous <- a[
  MappingClass ==
    "DIRECT_CURRENT_AMBIGUOUS"
]

cat("============================================================\n")
cat("DIRECT CURRENT AMBIGUOUS\n")
cat("============================================================\n")

cat(
  "Rows:",
  nrow(direct_ambiguous),
  "\n"
)

if (nrow(direct_ambiguous) > 0) {

  print(
    direct_ambiguous[
      ,
      .(
        GeneSymbol_raw,
        Direct_NEntrez,
        Direct_EntrezIDs
      )
    ]
  )
}


# ============================================================
# 6. Unresolved structural classification
# ============================================================

u <- copy(
  a[
    MappingClass ==
      "UNRESOLVED"
  ]
)

stopifnot(
  nrow(u) == 776
)


u[
  ,
  StructuralClass := fcase(

    grepl(
      "_",
      GeneSymbol_raw,
      fixed = TRUE
    ),
    "COMPOSITE_UNDERSCORE",

    grepl(
      "^LOC[0-9]+$",
      GeneSymbol_raw
    ),
    "LOC_STYLE",

    grepl(
      "^[A-Za-z0-9.-]+$",
      GeneSymbol_raw
    ),
    "SIMPLE_SYMBOL_LIKE",

    default =
      "OTHER_COMPLEX"
  )
]


# ============================================================
# 7. Conservative punctuation-only rescue
#
# Principle:
#   Normalize ONLY by removing punctuation from current symbols.
#
# Example:
#   ACTA2-AS1 -> ACTA2AS1
#   NKX2-3    -> NKX23
#
# Candidate accepted only if:
#   - raw feature is SIMPLE_SYMBOL_LIKE
#   - normalized key matches exactly ONE current SYMBOL
#   - maps to exactly ONE EntrezID
#   - candidate differs only by punctuation/case
#
# Still AUDIT ONLY.
# ============================================================

normalize_token <- function(x) {

  x <- toupper(x)

  gsub(
    "[^A-Z0-9]",
    "",
    x
  )
}


current_symbols <- keys(
  org.Hs.eg.db,
  keytype = "SYMBOL"
)

current_map <- as.data.table(
  suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = current_symbols,
      keytype = "SYMBOL",
      columns = c(
        "SYMBOL",
        "ENTREZID"
      )
    )
  )
)

current_map <- unique(
  current_map[
    !is.na(SYMBOL) &
      !is.na(ENTREZID)
  ]
)

current_map[
  ,
  NormalizedToken :=
    normalize_token(SYMBOL)
]

current_norm <- current_map[
  ,
  .(
    NCurrentSymbols =
      uniqueN(SYMBOL),

    NEntrez =
      uniqueN(ENTREZID),

    CurrentSymbols =
      paste(
        sort(
          unique(SYMBOL)
        ),
        collapse = "|"
      ),

    EntrezIDs =
      paste(
        sort(
          unique(ENTREZID)
        ),
        collapse = "|"
      )
  ),
  by = NormalizedToken
]


u[
  ,
  NormalizedRawToken :=
    normalize_token(
      GeneSymbol_raw
    )
]

idx_norm <- match(
  u$NormalizedRawToken,
  current_norm$NormalizedToken
)

u[
  ,
  Punct_NCurrentSymbols :=
    fifelse(
      !is.na(idx_norm),
      current_norm$NCurrentSymbols[
        idx_norm
      ],
      0L
    )
]

u[
  ,
  Punct_NEntrez :=
    fifelse(
      !is.na(idx_norm),
      current_norm$NEntrez[
        idx_norm
      ],
      0L
    )
]

u[
  ,
  Punct_CurrentSymbols :=
    fifelse(
      !is.na(idx_norm),
      current_norm$CurrentSymbols[
        idx_norm
      ],
      NA_character_
    )
]

u[
  ,
  Punct_EntrezIDs :=
    fifelse(
      !is.na(idx_norm),
      current_norm$EntrezIDs[
        idx_norm
      ],
      NA_character_
    )
]


# only simple symbols are eligible
u[
  ,
  PunctuationOnlyCandidate :=
    StructuralClass ==
      "SIMPLE_SYMBOL_LIKE" &
    Punct_NCurrentSymbols == 1 &
    Punct_NEntrez == 1 &
    !is.na(Punct_CurrentSymbols) &
    GeneSymbol_raw !=
      Punct_CurrentSymbols
]


# ------------------------------------------------------------
# Does target already exist as an original raw feature?
# ------------------------------------------------------------

u[
  ,
  PunctTargetAlreadyRaw :=
    PunctuationOnlyCandidate &
    Punct_CurrentSymbols %chin%
      raw_symbols
]


# ------------------------------------------------------------
# Candidate multiplicity among unresolved raw features
# ------------------------------------------------------------

punct_mult <- u[
  PunctuationOnlyCandidate == TRUE,
  .(
    SourceFeatureN = .N,
    SourceSymbols =
      paste(
        GeneSymbol_raw,
        collapse = "|"
      )
  ),
  by = Punct_EntrezIDs
]

idx_pm <- match(
  u$Punct_EntrezIDs,
  punct_mult$Punct_EntrezIDs
)

u[
  ,
  PunctSourceFeatureN :=
    fifelse(
      PunctuationOnlyCandidate &
        !is.na(idx_pm),
      punct_mult$SourceFeatureN[
        idx_pm
      ],
      NA_integer_
    )
]


# ------------------------------------------------------------
# Preliminary punctuation rescue status
# ------------------------------------------------------------

u[
  ,
  PunctuationRescueStatus := fcase(

    !PunctuationOnlyCandidate,
    "NO_PUNCTUATION_RESCUE",

    PunctuationOnlyCandidate &
      PunctTargetAlreadyRaw,
    "PUNCT_TARGET_ALREADY_PRESENT",

    PunctuationOnlyCandidate &
      PunctSourceFeatureN > 1,
    "PUNCT_MANY_TO_ONE",

    PunctuationOnlyCandidate &
      !PunctTargetAlreadyRaw &
      PunctSourceFeatureN == 1,
    "SAFE_PUNCTUATION_CANDIDATE",

    default =
      "PUNCTUATION_REVIEW"
  )
]


# ============================================================
# 8. Summaries
# ============================================================

cat("\n============================================================\n")
cat("ALIAS TARGET COLLISION EXPRESSION AUDIT\n")
cat("============================================================\n")

cat(
  "Pairs:",
  nrow(target_collision_audit),
  "\n"
)

if (nrow(target_collision_audit) > 0) {

  print(
    target_collision_audit[
      ,
      .(
        Symbol1,
        Symbol2,
        CandidateEntrezID,
        Pearson,
        Spearman,
        MeanAbsDiff,
        RMSE,
        MaxAbsDiff,
        Relationship
      )
    ]
  )
}


cat("\n============================================================\n")
cat("MANY-TO-ONE EXPRESSION AUDIT\n")
cat("============================================================\n")

cat(
  "Pairwise comparisons:",
  nrow(many_to_one_audit),
  "\n"
)

if (nrow(many_to_one_audit) > 0) {

  print(
    many_to_one_audit[
      ,
      .(
        Symbol1,
        Symbol2,
        CandidateEntrezID,
        Pearson,
        Spearman,
        MeanAbsDiff,
        RMSE,
        MaxAbsDiff,
        Relationship
      )
    ]
  )
}


cat("\n============================================================\n")
cat("UNRESOLVED STRUCTURAL CLASSES\n")
cat("============================================================\n")

structural_summary <- u[
  ,
  .N,
  by = StructuralClass
][order(-N)]

print(
  structural_summary
)


cat("\n============================================================\n")
cat("PUNCTUATION-ONLY RESCUE SUMMARY\n")
cat("============================================================\n")

punct_summary <- u[
  ,
  .N,
  by = PunctuationRescueStatus
][order(-N)]

print(
  punct_summary
)


cat("\n============================================================\n")
cat("SAFE PUNCTUATION CANDIDATES\n")
cat("============================================================\n")

safe_punct <- u[
  PunctuationRescueStatus ==
    "SAFE_PUNCTUATION_CANDIDATE"
]

cat(
  "Rows:",
  nrow(safe_punct),
  "\n"
)

if (nrow(safe_punct) > 0) {

  print(
    safe_punct[
      ,
      .(
        GeneSymbol_raw,
        Punct_CurrentSymbols,
        Punct_EntrezIDs,
        StructuralClass
      )
    ]
  )
}


cat("\n============================================================\n")
cat("PUNCTUATION COLLISIONS\n")
cat("============================================================\n")

punct_problem <- u[
  PunctuationRescueStatus %chin% c(
    "PUNCT_TARGET_ALREADY_PRESENT",
    "PUNCT_MANY_TO_ONE",
    "PUNCTUATION_REVIEW"
  )
]

cat(
  "Rows:",
  nrow(punct_problem),
  "\n"
)

if (nrow(punct_problem) > 0) {

  print(
    punct_problem[
      ,
      .(
        GeneSymbol_raw,
        Punct_CurrentSymbols,
        Punct_EntrezIDs,
        PunctTargetAlreadyRaw,
        PunctSourceFeatureN,
        PunctuationRescueStatus
      )
    ]
  )
}


# ============================================================
# 9. Save
# ============================================================

fwrite(
  target_collision_audit,
  file.path(
    out_dir,
    "GSE93624_alias_target_present_expression_collision_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  many_to_one_audit,
  file.path(
    out_dir,
    "GSE93624_many_to_one_expression_collision_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  direct_ambiguous,
  file.path(
    out_dir,
    "GSE93624_direct_current_ambiguous.tsv"
  ),
  sep = "\t"
)

fwrite(
  u,
  file.path(
    out_dir,
    "GSE93624_unresolved_symbol_structural_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  structural_summary,
  file.path(
    out_dir,
    "GSE93624_unresolved_structural_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  punct_summary,
  file.path(
    out_dir,
    "GSE93624_punctuation_rescue_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  safe_punct,
  file.path(
    out_dir,
    "GSE93624_safe_punctuation_candidates.tsv"
  ),
  sep = "\t"
)

fwrite(
  punct_problem,
  file.path(
    out_dir,
    "GSE93624_punctuation_collision_candidates.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 10. Final
# ============================================================

stopifnot(
  nrow(u) == 776
)

cat("\nIMPORTANT:\n")
cat(
  "- Expression collisions were audited but not resolved by aggregation.\n",
  "- Punctuation-only mappings are candidate identity rescues only.\n",
  "- Composite underscore features are not split into component genes.\n",
  "- No feature was renamed, merged, excluded, or added.\n",
  sep = ""
)

cat("\nFINAL STATUS:\n")
cat("PASS_GSE93624_COLLISION_AND_UNRESOLVED_FORENSIC_AUDIT\n")
cat("============================================================\n")
