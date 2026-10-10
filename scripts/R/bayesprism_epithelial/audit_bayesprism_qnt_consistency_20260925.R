#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

suppressPackageStartupMessages({
  library(BayesPrism)
  library(data.table)
})

# ==============================================================
# BayesPrism -> EcoTyper q_nt consistency audit
# Version: 2026-09-25 v2
#
# PURPOSE
# -------
#
# Stage 1:
#
#   q_nt =
#     sum_g Z_ntg /
#     (R_n * theta0_nt)
#
# Stage 2 conditional:
#
#   q_nt^(S2|Epi) =
#     sum_g Z_ntg /
#     (R_S2,n * theta0_S2,nt)
#
# Stage 2 whole-biopsy hierarchical:
#
#   q_nt^(hier) =
#     sum_g Z_ntg /
#     (
#       R_S1,n *
#       theta0_S1,n,Epithelial *
#       theta0_S2,nt
#     )
#
# Two depth definitions:
#
#   R_raw = sum_g X_ng
#   R_eff = sum_g floor(X_ng)
#
# Existing CORE12 matrix:
#
#   G_current =
#     Z / sum_g(Z) * 1e6
#
# Eq2-style analogue:
#
#   G_Eq2 =
#     Z / (R * theta0) * 1e6
#
# Therefore:
#
#   G_Eq2 / G_current = q_nt
#
#
# IMPORTANT TECHNICAL-ZERO RULE
# -----------------------------
#
# Stage2 pseudo-mixture contains fractional posterior values.
# For the four predefined technical-zero samples:
#
#   R_raw > 0
#
# may coexist with:
#
#   R_eff = sum_g floor(X_g) = 0
#
# Therefore:
#
# - R_eff <= 0 is allowed ONLY for samples already present in
#   technical_exclusion_samples.tsv.
#
# - For those samples, q_cond_eff is mathematically undefined
#   and MUST remain NA.
#
# - No epsilon is added.
#
# - Any R_eff <= 0 in a NON-technical sample is a hard failure.
#
#
# READ-ONLY GUARANTEE
# -------------------
#
# This script does NOT modify:
#
# - BayesPrism outputs
# - G-like matrices
# - EcoTyper inputs
# - EcoTyper runtime
# - running NMF jobs
#
# It writes only into:
#
#   03_reference/combined_reference/qc/
#       bayesprism_qnt_audit/<timestamp>/
#
# ==============================================================


# ==============================================================
# Project paths
# ==============================================================

ROOT <- "/home/mazekai/IBD_EcoTyper"


# --------------------------------------------------------------
# Stage 1 Broad14
# --------------------------------------------------------------

S1_RESULT <- file.path(
  ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "14_bayesprism_broad_run",
  "bayesprism_broad14_result.rds"
)

S1_Z_DIR <- file.path(
  ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "14_bayesprism_broad_run",
  "celltype_expression"
)


# --------------------------------------------------------------
# Stage 2 epithelial7
# --------------------------------------------------------------

S2_RESULT <- file.path(
  ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "20_bayesprism_epithelial_full_run",
  "bayesprism_epithelial7_result.rds"
)

S2_Z_DIR <- file.path(
  ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "20_bayesprism_epithelial_full_run",
  "celltype_expression"
)

S2_PSEUDOMIXTURE <- file.path(
  ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "17_bayesprism_epithelial_input",
  "bayesprism_stage2_epithelial_pseudomixture.rds"
)


# --------------------------------------------------------------
# Existing formal CORE12 G-like input
# --------------------------------------------------------------

CORE12_GLIKE_DIR <- file.path(
  ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "27_ecotyper_core12_glike",
  "expression_matrices"
)


# --------------------------------------------------------------
# Predefined technical exclusions
# --------------------------------------------------------------

TECH_FILE <- file.path(
  ROOT,
  "03_reference",
  "combined_reference",
  "output",
  "26_panibd_panel_classification",
  "technical_exclusion_samples.tsv"
)


# --------------------------------------------------------------
# Audit output
# --------------------------------------------------------------

AUDIT_ROOT <- file.path(
  ROOT,
  "03_reference",
  "combined_reference",
  "qc",
  "bayesprism_qnt_audit"
)

RUN_ID <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

OUT <- file.path(
  AUDIT_ROOT,
  RUN_ID
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================
# Locked cell-type panels
# ==============================================================

CORE10 <- c(
  "B",
  "CD4_T",
  "CD8_T",
  "Endothelial",
  "Fibroblast",
  "Glial",
  "ILC",
  "Mast",
  "Monocyte_Macrophage",
  "Pericyte"
)

STAGE2_CORE2 <- c(
  "Goblet",
  "Enteroendocrine"
)

CORE12 <- c(
  CORE10,
  STAGE2_CORE2
)

HIER_PARENT <- "Epithelial"


# ==============================================================
# Generic helpers
# ==============================================================

stopf <- function(cond, msg) {

  if (!isTRUE(cond)) {

    stop(
      msg,
      call. = FALSE
    )
  }
}


write_tsv <- function(x, path) {

  fwrite(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}


safe_ratio <- function(num, den) {

  out <- rep(
    NA_real_,
    length(num)
  )

  ok <- (
    is.finite(num) &
      is.finite(den) &
      den > 0
  )

  out[ok] <- (
    num[ok] /
      den[ok]
  )

  out
}


get_z_path <- function(
    z_dir,
    cell_type) {

  file.path(
    z_dir,
    paste0(
      "bayesprism_expression_",
      cell_type,
      ".rds"
    )
  )
}


# ==============================================================
# Input path audit
# ==============================================================

audit_paths <- function() {

  required_files <- c(
    S1_RESULT,
    S2_RESULT,
    S2_PSEUDOMIXTURE,
    TECH_FILE
  )

  required_dirs <- c(
    S1_Z_DIR,
    S2_Z_DIR,
    CORE12_GLIKE_DIR
  )

  for (f in required_files) {

    stopf(
      file.exists(f),
      paste0(
        "[FAIL] Required file missing:\n",
        f
      )
    )
  }

  for (d in required_dirs) {

    stopf(
      dir.exists(d),
      paste0(
        "[FAIL] Required directory missing:\n",
        d
      )
    )
  }

  cat(
    "[PASS] All required input paths exist.\n"
  )
}


# ==============================================================
# Effective-depth validation
#
# Non-positive effective depth is allowed ONLY for samples
# explicitly listed in the technical-exclusion manifest.
# ==============================================================

audit_effective_depth <- function(
    label,
    sample_ids,
    R_raw,
    R_eff,
    technical_exclusion) {

  stopf(
    length(sample_ids) ==
      length(R_raw),
    paste0(
      label,
      ": sample/R_raw length mismatch."
    )
  )

  stopf(
    length(sample_ids) ==
      length(R_eff),
    paste0(
      label,
      ": sample/R_eff length mismatch."
    )
  )

  stopf(
    length(sample_ids) ==
      length(technical_exclusion),
    paste0(
      label,
      ": technical flag length mismatch."
    )
  )


  # ------------------------------------------------------------
  # Raw depth
  # ------------------------------------------------------------

  bad_raw <- (
    !is.finite(R_raw) |
      R_raw <= 0
  )

  if (any(bad_raw)) {

    stop(
      paste0(
        label,
        ": non-positive/non-finite R_raw detected in: ",
        paste(
          sample_ids[
            bad_raw
          ],
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Effective depth
  # ------------------------------------------------------------

  zero_eff <- (
    !is.finite(R_eff) |
      R_eff <= 0
  )

  zero_eff_main <- (
    zero_eff &
      !technical_exclusion
  )

  zero_eff_technical <- (
    zero_eff &
      technical_exclusion
  )


  # ------------------------------------------------------------
  # HARD FAIL:
  # effective depth <= 0 in a formal sample
  # ------------------------------------------------------------

  if (any(zero_eff_main)) {

    stop(
      paste0(
        label,
        ": non-positive R_eff detected in NON-technical samples: ",
        paste(
          sample_ids[
            zero_eff_main
          ],
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # ALLOWED:
  # predefined technical-zero samples
  # ------------------------------------------------------------

  if (any(zero_eff_technical)) {

    cat(
      "[INFO] ",
      label,
      ": ",
      sum(zero_eff_technical),
      " predefined technical-exclusion sample(s) ",
      "have R_eff <= 0.\n",
      sep = ""
    )

    cat(
      "       These samples will retain q_eff = NA and ",
      "will NOT enter primary summaries.\n"
    )

    cat(
      "       Samples: ",
      paste(
        sample_ids[
          zero_eff_technical
        ],
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }


  # ------------------------------------------------------------
  # Return flags for downstream forensic output
  # ------------------------------------------------------------

  data.frame(
    sample_id = sample_ids,

    R_raw = R_raw,

    R_eff = R_eff,

    R_eff_over_R_raw =
      safe_ratio(
        R_eff,
        R_raw
      ),

    technical_exclusion =
      technical_exclusion,

    effective_zero =
      zero_eff,

    technical_effective_zero =
      zero_eff_technical,

    stringsAsFactors = FALSE
  )
}


# ==============================================================
# Load one posterior-Z matrix
# ==============================================================

load_z <- function(
    z_dir,
    cell_type,
    expected_samples,
    expected_genes) {

  path <- get_z_path(
    z_dir,
    cell_type
  )

  stopf(
    file.exists(path),
    paste0(
      "Missing posterior-Z file:\n",
      path
    )
  )

  z <- readRDS(path)

  stopf(
    is.matrix(z) ||
      inherits(z, "Matrix"),
    paste0(
      cell_type,
      ": Z is not matrix/Matrix."
    )
  )

  stopf(
    nrow(z) ==
      length(expected_samples),
    paste0(
      cell_type,
      ": unexpected sample count."
    )
  )

  stopf(
    ncol(z) ==
      length(expected_genes),
    paste0(
      cell_type,
      ": unexpected gene count."
    )
  )

  stopf(
    identical(
      rownames(z),
      expected_samples
    ),
    paste0(
      cell_type,
      ": sample names/order mismatch."
    )
  )

  stopf(
    identical(
      colnames(z),
      expected_genes
    ),
    paste0(
      cell_type,
      ": gene names/order mismatch."
    )
  )

  stopf(
    all(is.finite(z)),
    paste0(
      cell_type,
      ": non-finite Z detected."
    )
  )

  stopf(
    min(z) >= -1e-10,
    paste0(
      cell_type,
      ": negative posterior mass detected."
    )
  )

  z
}


# ==============================================================
# Effective Dirichlet prior after state -> type aggregation
# ==============================================================

get_alpha_info <- function(bp) {

  map <- bp@prism@map

  state_names <- rownames(
    bp@prism@phi_cellState@phi
  )

  alpha <- bp@
    control_param$
    gibbs.control$
    alpha

  stopf(
    !is.null(alpha),
    "alpha missing from BayesPrism object."
  )

  if (length(alpha) == 1) {

    alpha_state <- rep(
      as.numeric(alpha),
      length(state_names)
    )

  } else {

    stopf(
      length(alpha) ==
        length(state_names),
      paste0(
        "alpha length = ",
        length(alpha),
        "; states = ",
        length(state_names)
      )
    )

    alpha_state <- as.numeric(alpha)
  }

  names(alpha_state) <- state_names

  stopf(
    setequal(
      unlist(map),
      state_names
    ),
    paste0(
      "State names in map do not match ",
      "phi_cellState."
    )
  )

  alpha_type <- sapply(
    names(map),
    function(ct) {

      sum(
        alpha_state[
          map[[ct]]
        ]
      )
    }
  )

  list(
    alpha_state =
      alpha_state,

    alpha_type =
      alpha_type,

    alpha_total =
      sum(alpha_state),

    n_states =
      length(state_names)
  )
}


# ==============================================================
# Validate/extract BayesPrism model metadata
# ==============================================================

extract_bp_metadata <- function(
    bp,
    label) {

  stopf(
    is(bp, "BayesPrism"),
    paste0(
      label,
      ": object is not BayesPrism."
    )
  )

  mixture <- bp@prism@mixture

  theta <- bp@
    posterior.initial.cellType@
    theta

  z_obj <- bp@
    posterior.initial.cellType@
    Z

  stopf(
    identical(
      rownames(mixture),
      rownames(theta)
    ),
    paste0(
      label,
      ": mixture/theta sample mismatch."
    )
  )

  stopf(
    identical(
      dimnames(z_obj)[[1]],
      rownames(mixture)
    ),
    paste0(
      label,
      ": Z/mixture sample mismatch."
    )
  )

  stopf(
    identical(
      dimnames(z_obj)[[2]],
      colnames(mixture)
    ),
    paste0(
      label,
      ": Z/mixture gene mismatch."
    )
  )

  stopf(
    identical(
      dimnames(z_obj)[[3]],
      colnames(theta)
    ),
    paste0(
      label,
      ": Z/theta cell-type mismatch."
    )
  )

  stopf(
    all(is.finite(mixture)),
    paste0(
      label,
      ": non-finite mixture values."
    )
  )

  stopf(
    min(mixture) >= 0,
    paste0(
      label,
      ": negative mixture values."
    )
  )

  list(
    mixture =
      mixture,

    theta =
      theta,

    map =
      bp@prism@map,

    alpha_info =
      get_alpha_info(bp),

    z_dim =
      dim(z_obj)
  )
}


# ==============================================================
# q summary
# ==============================================================

summarize_q <- function(
    df,
    q_col,
    model_name,
    denominator_name,
    technical = FALSE) {

  d0 <- if (technical) {

    df[
      df$technical_exclusion,
      ,
      drop = FALSE
    ]

  } else {

    df[
      !df$technical_exclusion,
      ,
      drop = FALSE
    ]
  }

  if (nrow(d0) == 0) {
    return(data.frame())
  }

  groups <- c(
    "ALL",
    unique(d0$cell_type)
  )

  strata <- c(
    "all",
    "Zmass_ge_100",
    "Zmass_ge_1000"
  )

  output <- list()

  ii <- 1L

  for (grp in groups) {

    d1 <- if (grp == "ALL") {

      d0

    } else {

      d0[
        d0$cell_type == grp,
        ,
        drop = FALSE
      ]
    }

    for (stratum in strata) {

      keep <- switch(
        stratum,

        all =
          rep(
            TRUE,
            nrow(d1)
          ),

        Zmass_ge_100 =
          d1$posterior_mass >= 100,

        Zmass_ge_1000 =
          d1$posterior_mass >= 1000
      )

      x <- d1[
        keep,
        q_col
      ]

      x <- x[
        is.finite(x)
      ]

      if (length(x) == 0) {
        next
      }

      abs_err <- abs(
        x - 1
      )

      qx <- quantile(
        x,
        probs = c(
          0.01,
          0.05,
          0.50,
          0.95,
          0.99
        ),
        names = FALSE,
        na.rm = TRUE
      )

      qe <- quantile(
        abs_err,
        probs = c(
          0.50,
          0.95,
          0.99
        ),
        names = FALSE,
        na.rm = TRUE
      )

      output[[ii]] <- data.frame(
        model =
          model_name,

        denominator =
          denominator_name,

        technical_exclusion_group =
          technical,

        cell_type =
          grp,

        mass_stratum =
          stratum,

        n =
          length(x),

        q_min =
          min(x),

        q_q01 =
          qx[1],

        q_q05 =
          qx[2],

        q_median =
          qx[3],

        q_q95 =
          qx[4],

        q_q99 =
          qx[5],

        q_max =
          max(x),

        abs_q_minus_1_median =
          qe[1],

        abs_q_minus_1_q95 =
          qe[2],

        abs_q_minus_1_q99 =
          qe[3],

        abs_q_minus_1_max =
          max(abs_err),

        pct_within_0_1pct =
          mean(
            abs_err <= 0.001
          ) * 100,

        pct_within_1pct =
          mean(
            abs_err <= 0.01
          ) * 100,

        pct_within_5pct =
          mean(
            abs_err <= 0.05
          ) * 100,

        stringsAsFactors = FALSE
      )

      ii <- ii + 1L
    }
  }

  if (length(output) == 0) {
    return(data.frame())
  }

  do.call(
    rbind,
    output
  )
}


# ==============================================================
# Theta reconstruction summary
# ==============================================================

summarize_theta_error <- function(
    df,
    error_col,
    model_name,
    denominator_name) {

  d0 <- df[
    !df$technical_exclusion,
    ,
    drop = FALSE
  ]

  groups <- c(
    "ALL",
    unique(d0$cell_type)
  )

  strata <- c(
    "all",
    "Zmass_ge_100",
    "Zmass_ge_1000"
  )

  output <- list()

  ii <- 1L

  for (grp in groups) {

    d1 <- if (grp == "ALL") {

      d0

    } else {

      d0[
        d0$cell_type == grp,
        ,
        drop = FALSE
      ]
    }

    for (stratum in strata) {

      keep <- switch(
        stratum,

        all =
          rep(
            TRUE,
            nrow(d1)
          ),

        Zmass_ge_100 =
          d1$posterior_mass >= 100,

        Zmass_ge_1000 =
          d1$posterior_mass >= 1000
      )

      x <- d1[
        keep,
        error_col
      ]

      x <- x[
        is.finite(x)
      ]

      if (length(x) == 0) {
        next
      }

      qx <- quantile(
        x,
        probs = c(
          0.50,
          0.95,
          0.99
        ),
        names = FALSE,
        na.rm = TRUE
      )

      output[[ii]] <- data.frame(
        model =
          model_name,

        denominator =
          denominator_name,

        cell_type =
          grp,

        mass_stratum =
          stratum,

        n =
          length(x),

        abs_error_median =
          qx[1],

        abs_error_q95 =
          qx[2],

        abs_error_q99 =
          qx[3],

        abs_error_max =
          max(x),

        stringsAsFactors = FALSE
      )

      ii <- ii + 1L
    }
  }

  do.call(
    rbind,
    output
  )
}


# ==============================================================
# Chunked exact matrix comparison
# ==============================================================

max_abs_diff_chunked <- function(
    a,
    b,
    chunk_size = 500L) {

  stopf(
    identical(
      dim(a),
      dim(b)
    ),
    "Matrices have different dimensions."
  )

  max_diff <- 0

  starts <- seq(
    1L,
    ncol(a),
    by = chunk_size
  )

  for (s in starts) {

    e <- min(
      ncol(a),
      s + chunk_size - 1L
    )

    d <- max(
      abs(
        a[
          ,
          s:e,
          drop = FALSE
        ] -
          b[
            ,
            s:e,
            drop = FALSE
          ]
      )
    )

    if (d > max_diff) {
      max_diff <- d
    }
  }

  max_diff
}


# ==============================================================
# Existing G-like formula audit
#
# Formal files should represent:
#
#   G_current =
#     Z / rowSums(Z) * 1e6
#
# We validate:
#
# - sample universe
# - gene universe
# - deterministic numerical spot-check
#
# ==============================================================

spotcheck_glike <- function(
    cell_type,
    z,
    expected_formal_samples,
    expected_genes) {

  glike_file <- file.path(
    CORE12_GLIKE_DIR,
    paste0(
      cell_type,
      ".txt"
    )
  )

  stopf(
    file.exists(glike_file),
    paste0(
      "Missing G-like file:\n",
      glike_file
    )
  )


  # ------------------------------------------------------------
  # Header / sample universe
  # ------------------------------------------------------------

  header <- names(
    fread(
      glike_file,
      nrows = 0,
      showProgress = FALSE
    )
  )

  stopf(
    length(header) >= 2,
    paste0(
      cell_type,
      ": malformed G-like header."
    )
  )

  stopf(
    header[1] == "Gene",
    paste0(
      cell_type,
      ": first G-like column is not Gene."
    )
  )

  observed_samples <-
    header[-1]

  stopf(
    setequal(
      observed_samples,
      expected_formal_samples
    ),
    paste0(
      cell_type,
      ": formal sample universe differs ",
      "from expected 2486 sample set."
    )
  )

  sample_order_identical <- identical(
    observed_samples,
    expected_formal_samples
  )


  # ------------------------------------------------------------
  # Deterministic spot-check set
  # ------------------------------------------------------------

  n_sample_check <- min(
    20L,
    length(expected_formal_samples)
  )

  n_gene_check <- min(
    50L,
    length(expected_genes)
  )

  sample_idx <- unique(
    round(
      seq(
        1,
        length(expected_formal_samples),
        length.out =
          n_sample_check
      )
    )
  )

  gene_idx <- unique(
    round(
      seq(
        1,
        length(expected_genes),
        length.out =
          n_gene_check
      )
    )
  )

  sample_keep <-
    expected_formal_samples[
      sample_idx
    ]

  gene_keep <-
    expected_genes[
      gene_idx
    ]


  # ------------------------------------------------------------
  # Load only selected sample columns
  # ------------------------------------------------------------

  dt <- fread(
    glike_file,
    select = c(
      "Gene",
      sample_keep
    ),
    showProgress = FALSE
  )

  observed_gene_universe <-
    dt$Gene

  stopf(
    setequal(
      observed_gene_universe,
      expected_genes
    ),
    paste0(
      cell_type,
      ": formal G-like gene universe differs ",
      "from expected BayesPrism gene universe."
    )
  )

  gene_order_identical <- identical(
    observed_gene_universe,
    expected_genes
  )

  dt_small <- dt[
    match(
      gene_keep,
      Gene
    )
  ]

  stopf(
    !anyNA(dt_small$Gene),
    paste0(
      cell_type,
      ": one or more spot-check genes missing."
    )
  )


  # ------------------------------------------------------------
  # Expected current G-like
  # ------------------------------------------------------------

  posterior_mass <- rowSums(z)

  stopf(
    all(
      posterior_mass[
        sample_keep
      ] > 0
    ),
    paste0(
      cell_type,
      ": selected formal sample has ",
      "non-positive posterior mass."
    )
  )

  expected <- z[
    sample_keep,
    gene_keep,
    drop = FALSE
  ]

  expected <- sweep(
    expected,
    MARGIN = 1,
    STATS =
      posterior_mass[
        sample_keep
      ],
    FUN = "/"
  )

  expected <- expected * 1e6


  # ------------------------------------------------------------
  # Observed formal G-like
  # ------------------------------------------------------------

  observed <- as.matrix(
    dt_small[
      ,
      ..sample_keep
    ]
  )

  rownames(observed) <-
    dt_small$Gene

  observed <- t(
    observed
  )

  observed <- observed[
    sample_keep,
    gene_keep,
    drop = FALSE
  ]


  # ------------------------------------------------------------
  # Difference
  # ------------------------------------------------------------

  abs_diff <- abs(
    observed -
      expected
  )

  data.frame(
    cell_type =
      cell_type,

    formal_n_samples =
      length(observed_samples),

    formal_n_genes =
      length(observed_gene_universe),

    sample_order_identical =
      sample_order_identical,

    gene_order_identical =
      gene_order_identical,

    n_samples_checked =
      length(sample_keep),

    n_genes_checked =
      length(gene_keep),

    n_values_checked =
      length(abs_diff),

    max_abs_diff =
      max(abs_diff),

    median_abs_diff =
      median(abs_diff),

    q99_abs_diff =
      quantile(
        abs_diff,
        0.99,
        names = FALSE
      ),

    stringsAsFactors = FALSE
  )
}


# ==============================================================
# Main
# ==============================================================

main <- function() {

  cat(
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  cat(
    "BayesPrism q_nt consistency audit\n"
  )

  cat(
    "Run ID: ",
    RUN_ID,
    "\n",
    sep = ""
  )

  cat(
    "Audit output: ",
    OUT,
    "\n",
    sep = ""
  )

  cat(
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n\n",
    sep = ""
  )


  # ============================================================
  # Paths
  # ==============================================================

  audit_paths()


  # ============================================================
  # Technical exclusion manifest
  # ==============================================================

  cat(
    "[INFO] Reading technical-exclusion manifest.\n"
  )

  tech <- fread(
    TECH_FILE,
    data.table = FALSE
  )

  stopf(
    "sample_id" %in%
      colnames(tech),
    "TECH_FILE has no sample_id column."
  )

  tech_ids <- unique(
    as.character(
      tech$sample_id
    )
  )

  stopf(
    length(tech_ids) > 0,
    "Technical-exclusion manifest is empty."
  )


  # ============================================================
  # STAGE 1
  # ==============================================================

  cat(
    "[INFO] Loading Stage 1 BayesPrism result.\n"
  )

  bp1 <- readRDS(
    S1_RESULT
  )

  s1 <- extract_bp_metadata(
    bp1,
    "Stage1"
  )

  samples1 <- rownames(
    s1$mixture
  )

  genes1 <- colnames(
    s1$mixture
  )

  theta1 <- s1$theta

  stopf(
    all(
      c(
        CORE10,
        HIER_PARENT
      ) %in%
        colnames(theta1)
    ),
    "Required Stage1 cell type missing."
  )

  stopf(
    all(
      tech_ids %in%
        samples1
    ),
    paste0(
      "Technical-exclusion manifest contains ",
      "samples absent from Stage1."
    )
  )

  tech_flag1 <- (
    samples1 %in%
      tech_ids
  )

  R1_raw <- rowSums(
    s1$mixture
  )

  R1_eff <- rowSums(
    floor(
      s1$mixture
    )
  )

  depth1 <- audit_effective_depth(
    label =
      "Stage1",

    sample_ids =
      samples1,

    R_raw =
      R1_raw,

    R_eff =
      R1_eff,

    technical_exclusion =
      tech_flag1
  )

  write_tsv(
    depth1,
    file.path(
      OUT,
      "stage1_depth_raw_vs_effective.tsv"
    )
  )

  theta1_epi <- theta1[
    ,
    HIER_PARENT
  ]

  alpha1_type <- s1$
    alpha_info$
    alpha_type

  alpha1_total <- s1$
    alpha_info$
    alpha_total


  # ------------------------------------------------------------
  # Release full Stage1 BayesPrism result
  # ------------------------------------------------------------

  rm(bp1)

  s1$mixture <- NULL

  gc()


  cat(
    "[INFO] Auditing Stage 1 q_nt, ",
    "R_raw/R_eff, and theta reconstruction.\n",
    sep = ""
  )


  # ============================================================
  # Stage1 q
  # ==============================================================

  stage1_types <- unique(
    c(
      CORE10,
      HIER_PARENT
    )
  )

  stage1_rows <- vector(
    "list",
    length(stage1_types)
  )

  glike_rows <- list()

  glike_idx <- 1L


  for (i in seq_along(stage1_types)) {

    ct <- stage1_types[i]

    z <- load_z(
      S1_Z_DIR,
      ct,
      samples1,
      genes1
    )

    S <- rowSums(z)

    th <- theta1[
      ,
      ct
    ]

    alpha_ct <- alpha1_type[
      ct
    ]

    q_raw <- safe_ratio(
      S,
      R1_raw * th
    )

    q_eff <- safe_ratio(
      S,
      R1_eff * th
    )


    # ----------------------------------------------------------
    # For Stage1, any R_eff zero should already have been
    # validated above.
    # ----------------------------------------------------------

    theta_hat_raw <- (
      S + alpha_ct
    ) / (
      R1_raw +
        alpha1_total
    )

    theta_hat_eff <- (
      S + alpha_ct
    ) / (
      R1_eff +
        alpha1_total
    )

    stage1_rows[[i]] <- data.frame(
      sample_id =
        samples1,

      cell_type =
        ct,

      panel_role =
        ifelse(
          ct %in%
            CORE10,
          "CORE10",
          "Hierarchy_parent"
        ),

      posterior_mass =
        S,

      theta0 =
        th,

      R_raw =
        R1_raw,

      R_eff =
        R1_eff,

      effective_zero =
        depth1$effective_zero,

      q_raw =
        q_raw,

      q_eff =
        q_eff,

      G_Eq2_over_Gcurrent_raw =
        q_raw,

      G_Eq2_over_Gcurrent_eff =
        q_eff,

      theta_hat_raw =
        theta_hat_raw,

      theta_hat_eff =
        theta_hat_eff,

      theta_abs_error_raw =
        abs(
          theta_hat_raw -
            th
        ),

      theta_abs_error_eff =
        abs(
          theta_hat_eff -
            th
        ),

      technical_exclusion =
        tech_flag1,

      stringsAsFactors = FALSE
    )


    # ----------------------------------------------------------
    # Formal G-like audit for CORE10 only
    # ----------------------------------------------------------

    if (ct %in% CORE10) {

      expected_formal_samples <-
        samples1[
          !tech_flag1
        ]

      glike_rows[[glike_idx]] <-
        spotcheck_glike(
          cell_type =
            ct,

          z =
            z,

          expected_formal_samples =
            expected_formal_samples,

          expected_genes =
            genes1
        )

      glike_idx <-
        glike_idx + 1L
    }


    rm(
      z,
      S,
      th,
      q_raw,
      q_eff,
      theta_hat_raw,
      theta_hat_eff
    )

    gc()
  }


  stage1_long <- do.call(
    rbind,
    stage1_rows
  )

  write_tsv(
    stage1_long,
    file.path(
      OUT,
      "stage1_qnt_long.tsv"
    )
  )

  rm(stage1_rows)

  gc()


  # ============================================================
  # STAGE 2
  # ==============================================================

  cat(
    "[INFO] Loading Stage 2 BayesPrism result.\n"
  )

  bp2 <- readRDS(
    S2_RESULT
  )

  s2 <- extract_bp_metadata(
    bp2,
    "Stage2"
  )

  samples2 <- rownames(
    s2$mixture
  )

  genes2 <- colnames(
    s2$mixture
  )

  theta2 <- s2$theta

  stopf(
    identical(
      samples2,
      samples1
    ),
    "Stage1/Stage2 sample order mismatch."
  )

  stopf(
    identical(
      genes2,
      genes1
    ),
    "Stage1/Stage2 gene order mismatch."
  )

  stopf(
    all(
      STAGE2_CORE2 %in%
        colnames(theta2)
    ),
    "Goblet or Enteroendocrine missing from Stage2."
  )

  tech_flag2 <- (
    samples2 %in%
      tech_ids
  )

  R2_raw <- rowSums(
    s2$mixture
  )

  R2_eff <- rowSums(
    floor(
      s2$mixture
    )
  )


  # ------------------------------------------------------------
  # THIS IS THE KEY V2 FIX.
  #
  # R_eff == 0 is allowed for predefined technical exclusions.
  # ------------------------------------------------------------

  depth2 <- audit_effective_depth(
    label =
      "Stage2_conditional",

    sample_ids =
      samples2,

    R_raw =
      R2_raw,

    R_eff =
      R2_eff,

    technical_exclusion =
      tech_flag2
  )

  write_tsv(
    depth2,
    file.path(
      OUT,
      "stage2_depth_raw_vs_effective.tsv"
    )
  )


  alpha2_type <- s2$
    alpha_info$
    alpha_type

  alpha2_total <- s2$
    alpha_info$
    alpha_total


  # ============================================================
  # Stage2 pseudo-mixture triple consistency
  # ==============================================================

  cat(
    "[INFO] Auditing Stage 2 pseudo-mixture consistency.\n"
  )


  # ------------------------------------------------------------
  # Stored pseudo-mixture
  # ------------------------------------------------------------

  pseudo_external <- readRDS(
    S2_PSEUDOMIXTURE
  )

  stopf(
    identical(
      dim(pseudo_external),
      dim(s2$mixture)
    ),
    "Stored pseudo-mixture dimension mismatch."
  )

  stopf(
    identical(
      rownames(pseudo_external),
      samples2
    ),
    "Stored pseudo-mixture sample mismatch."
  )

  stopf(
    identical(
      colnames(pseudo_external),
      genes2
    ),
    "Stored pseudo-mixture gene mismatch."
  )


  diff_external_internal <-
    max_abs_diff_chunked(
      pseudo_external,
      s2$mixture
    )


  # ------------------------------------------------------------
  # Stage1 Epithelial Z
  # ------------------------------------------------------------

  epi_z <- load_z(
    S1_Z_DIR,
    HIER_PARENT,
    samples1,
    genes1
  )

  diff_epi_s2 <-
    max_abs_diff_chunked(
      epi_z,
      s2$mixture
    )

  diff_epi_rowsum <-
    max(
      abs(
        rowSums(epi_z) -
          R2_raw
      )
    )


  pseudomix_check <- data.frame(
    check = c(
      "stored_stage2_pseudomixture_vs_stage2_prism_mixture",
      "stage1_epithelial_Z_vs_stage2_prism_mixture",
      "stage1_epithelial_Z_rowsum_vs_stage2_Rraw"
    ),

    max_abs_difference = c(
      diff_external_internal,
      diff_epi_s2,
      diff_epi_rowsum
    ),

    stringsAsFactors = FALSE
  )

  write_tsv(
    pseudomix_check,
    file.path(
      OUT,
      "stage2_pseudomixture_consistency.tsv"
    )
  )


  rm(
    pseudo_external,
    epi_z
  )

  gc()


  # ------------------------------------------------------------
  # Release full Stage2 result before loading individual Z files
  # ------------------------------------------------------------

  rm(bp2)

  s2$mixture <- NULL

  gc()


  cat(
    "[INFO] Auditing Stage 2 conditional q_nt, ",
    "R_raw/R_eff, and theta reconstruction.\n",
    sep = ""
  )


  # ============================================================
  # Stage2 q
  # ==============================================================

  stage2_types <- colnames(
    theta2
  )

  stage2_rows <- vector(
    "list",
    length(stage2_types)
  )


  for (i in seq_along(stage2_types)) {

    ct <- stage2_types[i]

    z <- load_z(
      S2_Z_DIR,
      ct,
      samples2,
      genes2
    )

    S <- rowSums(z)

    th2 <- theta2[
      ,
      ct
    ]

    alpha_ct <- alpha2_type[
      ct
    ]


    # ==========================================================
    # Stage2 conditional q
    # ==========================================================

    q_cond_raw <- safe_ratio(
      S,
      R2_raw * th2
    )


    # ----------------------------------------------------------
    # For technical-zero samples:
    #
    # R2_eff == 0
    #
    # safe_ratio() deliberately returns NA.
    #
    # No epsilon.
    # No artificial replacement.
    # ----------------------------------------------------------

    q_cond_eff <- safe_ratio(
      S,
      R2_eff * th2
    )


    # ==========================================================
    # Theta reconstruction
    #
    # Even when R_eff == 0:
    #
    #   theta_hat =
    #     alpha_type / alpha_total
    #
    # remains mathematically meaningful under the prior.
    # ==========================================================

    theta2_hat_raw <- (
      S + alpha_ct
    ) / (
      R2_raw +
        alpha2_total
    )

    theta2_hat_eff <- (
      S + alpha_ct
    ) / (
      R2_eff +
        alpha2_total
    )


    # ==========================================================
    # Whole-biopsy hierarchical q
    # ==========================================================

    theta_hier <- (
      theta1_epi *
        th2
    )

    q_hier_raw <- safe_ratio(
      S,
      R1_raw *
        theta_hier
    )

    q_hier_eff <- safe_ratio(
      S,
      R1_eff *
        theta_hier
    )


    stage2_rows[[i]] <- data.frame(
      sample_id =
        samples2,

      cell_type =
        ct,

      panel_role =
        ifelse(
          ct %in%
            STAGE2_CORE2,
          "CORE12_stage2",
          "Stage2_other"
        ),

      posterior_mass =
        S,

      theta2_conditional =
        th2,

      theta1_epithelial =
        theta1_epi,

      theta_hierarchical =
        theta_hier,

      R2_raw =
        R2_raw,

      R2_eff =
        R2_eff,

      R2_effective_zero =
        depth2$effective_zero,

      R2_technical_effective_zero =
        depth2$technical_effective_zero,

      q_cond_raw =
        q_cond_raw,

      q_cond_eff =
        q_cond_eff,

      G_Eq2_cond_over_Gcurrent_raw =
        q_cond_raw,

      G_Eq2_cond_over_Gcurrent_eff =
        q_cond_eff,

      theta2_hat_raw =
        theta2_hat_raw,

      theta2_hat_eff =
        theta2_hat_eff,

      theta2_abs_error_raw =
        abs(
          theta2_hat_raw -
            th2
        ),

      theta2_abs_error_eff =
        abs(
          theta2_hat_eff -
            th2
        ),

      R1_raw =
        R1_raw,

      R1_eff =
        R1_eff,

      q_hier_raw =
        q_hier_raw,

      q_hier_eff =
        q_hier_eff,

      G_Eq2_hier_over_Gcurrent_raw =
        q_hier_raw,

      G_Eq2_hier_over_Gcurrent_eff =
        q_hier_eff,

      technical_exclusion =
        tech_flag2,

      stringsAsFactors = FALSE
    )


    # ----------------------------------------------------------
    # Formal G-like audit for the two Stage2 CORE12 types
    # ----------------------------------------------------------

    if (ct %in% STAGE2_CORE2) {

      expected_formal_samples <-
        samples2[
          !tech_flag2
        ]

      glike_rows[[glike_idx]] <-
        spotcheck_glike(
          cell_type =
            ct,

          z =
            z,

          expected_formal_samples =
            expected_formal_samples,

          expected_genes =
            genes2
        )

      glike_idx <-
        glike_idx + 1L
    }


    rm(
      z,
      S,
      th2,
      theta_hier,
      q_cond_raw,
      q_cond_eff,
      q_hier_raw,
      q_hier_eff,
      theta2_hat_raw,
      theta2_hat_eff
    )

    gc()
  }


  stage2_long <- do.call(
    rbind,
    stage2_rows
  )

  write_tsv(
    stage2_long,
    file.path(
      OUT,
      "stage2_qnt_long.tsv"
    )
  )


  # ============================================================
  # Explicit technical-zero forensic table
  # ==============================================================

  stage2_technical_zero <- stage2_long[
    stage2_long$
      R2_technical_effective_zero,
    ,
    drop = FALSE
  ]

  if (nrow(stage2_technical_zero) > 0) {

    write_tsv(
      stage2_technical_zero,
      file.path(
        OUT,
        "stage2_technical_zero_forensic.tsv"
      )
    )
  }


  rm(stage2_rows)

  gc()


  # ============================================================
  # Existing G-like formula audit
  # ==============================================================

  glike_check <- do.call(
    rbind,
    glike_rows
  )

  write_tsv(
    glike_check,
    file.path(
      OUT,
      "core12_glike_formula_spotcheck.tsv"
    )
  )


  # ============================================================
  # Combined CORE12
  # ==============================================================

  core12_s1 <- stage1_long[
    stage1_long$cell_type %in%
      CORE10,
    ,
    drop = FALSE
  ]

  core12_s1_out <- data.frame(
    sample_id =
      core12_s1$sample_id,

    cell_type =
      core12_s1$cell_type,

    source =
      "Stage1",

    posterior_mass =
      core12_s1$posterior_mass,

    q_raw =
      core12_s1$q_raw,

    q_eff =
      core12_s1$q_eff,

    technical_exclusion =
      core12_s1$technical_exclusion,

    stringsAsFactors = FALSE
  )


  core12_s2 <- stage2_long[
    stage2_long$cell_type %in%
      STAGE2_CORE2,
    ,
    drop = FALSE
  ]

  core12_s2_out <- data.frame(
    sample_id =
      core12_s2$sample_id,

    cell_type =
      core12_s2$cell_type,

    source =
      "Stage2_hierarchical",

    posterior_mass =
      core12_s2$posterior_mass,

    q_raw =
      core12_s2$q_hier_raw,

    q_eff =
      core12_s2$q_hier_eff,

    technical_exclusion =
      core12_s2$technical_exclusion,

    stringsAsFactors = FALSE
  )


  core12_long <- rbind(
    core12_s1_out,
    core12_s2_out
  )

  write_tsv(
    core12_long,
    file.path(
      OUT,
      "core12_qnt_long.tsv"
    )
  )


  # ============================================================
  # Main q summaries
  # ==============================================================

  summary_list <- list(
    summarize_q(
      stage1_long[
        stage1_long$cell_type %in%
          CORE10,
        ,
        drop = FALSE
      ],
      "q_raw",
      "Stage1_CORE10",
      "R_raw"
    ),

    summarize_q(
      stage1_long[
        stage1_long$cell_type %in%
          CORE10,
        ,
        drop = FALSE
      ],
      "q_eff",
      "Stage1_CORE10",
      "R_eff"
    ),

    summarize_q(
      stage2_long,
      "q_cond_raw",
      "Stage2_conditional",
      "R2_raw"
    ),

    summarize_q(
      stage2_long,
      "q_cond_eff",
      "Stage2_conditional",
      "R2_eff"
    ),

    summarize_q(
      stage2_long,
      "q_hier_raw",
      "Stage2_hierarchical",
      "R1_raw"
    ),

    summarize_q(
      stage2_long,
      "q_hier_eff",
      "Stage2_hierarchical",
      "R1_eff"
    ),

    summarize_q(
      core12_long,
      "q_raw",
      "CORE12",
      "raw"
    ),

    summarize_q(
      core12_long,
      "q_eff",
      "CORE12",
      "effective"
    )
  )


  summary_list <- summary_list[
    vapply(
      summary_list,
      nrow,
      integer(1)
    ) > 0
  ]


  q_summary <- do.call(
    rbind,
    summary_list
  )

  write_tsv(
    q_summary,
    file.path(
      OUT,
      "q_summary_main.tsv"
    )
  )


  # ============================================================
  # Technical-exclusion summary
  # ==============================================================

  technical_summary_list <- list(
    summarize_q(
      stage1_long[
        stage1_long$cell_type %in%
          CORE10,
        ,
        drop = FALSE
      ],
      "q_raw",
      "Stage1_CORE10",
      "R_raw",
      technical = TRUE
    ),

    summarize_q(
      stage2_long,
      "q_cond_raw",
      "Stage2_conditional",
      "R2_raw",
      technical = TRUE
    ),

    summarize_q(
      stage2_long,
      "q_cond_eff",
      "Stage2_conditional",
      "R2_eff",
      technical = TRUE
    ),

    summarize_q(
      stage2_long,
      "q_hier_raw",
      "Stage2_hierarchical",
      "R1_raw",
      technical = TRUE
    ),

    summarize_q(
      stage2_long,
      "q_hier_eff",
      "Stage2_hierarchical",
      "R1_eff",
      technical = TRUE
    )
  )


  technical_summary_list <-
    technical_summary_list[
      vapply(
        technical_summary_list,
        nrow,
        integer(1)
      ) > 0
    ]


  if (length(technical_summary_list) > 0) {

    technical_summary <- do.call(
      rbind,
      technical_summary_list
    )

    write_tsv(
      technical_summary,
      file.path(
        OUT,
        "q_summary_technical_exclusions.tsv"
      )
    )
  }


  # ============================================================
  # Theta reconstruction
  # ==============================================================

  theta_error <- rbind(
    summarize_theta_error(
      stage1_long,
      "theta_abs_error_eff",
      "Stage1",
      "R_eff"
    ),

    summarize_theta_error(
      stage2_long,
      "theta2_abs_error_eff",
      "Stage2_conditional",
      "R2_eff"
    )
  )

  write_tsv(
    theta_error,
    file.path(
      OUT,
      "theta0_reconstruction_error.tsv"
    )
  )


  # ============================================================
  # Depth summaries
  # ==============================================================

  summarize_depth <- function(
      df,
      stage_name) {

    output <- list()

    for (
      group_name in c(
        "formal",
        "technical_exclusion"
      )
    ) {

      d <- if (
        group_name == "formal"
      ) {

        df[
          !df$technical_exclusion,
          ,
          drop = FALSE
        ]

      } else {

        df[
          df$technical_exclusion,
          ,
          drop = FALSE
        ]
      }

      if (nrow(d) == 0) {
        next
      }

      x <- d$
        R_eff_over_R_raw

      x <- x[
        is.finite(x)
      ]

      if (length(x) == 0) {
        next
      }

      qx <- quantile(
        x,
        probs = c(
          0,
          0.01,
          0.05,
          0.50,
          0.95,
          0.99,
          1
        ),
        names = FALSE
      )

      output[[group_name]] <- data.frame(
        stage =
          stage_name,

        sample_group =
          group_name,

        n =
          length(x),

        n_R_eff_zero =
          sum(
            d$effective_zero
          ),

        min =
          qx[1],

        q01 =
          qx[2],

        q05 =
          qx[3],

        median =
          qx[4],

        q95 =
          qx[5],

        q99 =
          qx[6],

        max =
          qx[7],

        stringsAsFactors = FALSE
      )
    }

    do.call(
      rbind,
      output
    )
  }


  depth_summary <- rbind(
    summarize_depth(
      depth1,
      "Stage1"
    ),

    summarize_depth(
      depth2,
      "Stage2"
    )
  )

  write_tsv(
    depth_summary,
    file.path(
      OUT,
      "depth_effective_vs_raw_summary.tsv"
    )
  )


  # ============================================================
  # Model metadata
  # ==============================================================

  model_metadata <- data.frame(
    stage = c(
      "Stage1",
      "Stage2"
    ),

    samples = c(
      length(samples1),
      length(samples2)
    ),

    genes = c(
      length(genes1),
      length(genes2)
    ),

    cell_types = c(
      ncol(theta1),
      ncol(theta2)
    ),

    cell_states = c(
      s1$alpha_info$n_states,
      s2$alpha_info$n_states
    ),

    alpha_total = c(
      alpha1_total,
      alpha2_total
    ),

    BayesPrism_version =
      as.character(
        packageVersion(
          "BayesPrism"
        )
      ),

    stringsAsFactors = FALSE
  )

  write_tsv(
    model_metadata,
    file.path(
      OUT,
      "model_metadata.tsv"
    )
  )


  # ============================================================
  # PRIMARY CORE12 RESULT
  #
  # IMPORTANT:
  # only formal samples
  # only posterior mass >= 1000
  # ==============================================================

  cat(
    "\n",
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  cat(
    "PRIMARY CORE12 RESULT\n"
  )

  cat(
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )


  primary <- core12_long[
    !core12_long$technical_exclusion &
      core12_long$posterior_mass >= 1000 &
      is.finite(
        core12_long$q_eff
      ),
    ,
    drop = FALSE
  ]


  if (nrow(primary) == 0) {

    cat(
      "[WARN] No eligible CORE12 profiles ",
      "with posterior mass >= 1000.\n"
    )

  } else {

    abs_err <- abs(
      primary$q_eff -
        1
    )

    qr <- quantile(
      primary$q_eff,
      c(
        0.01,
        0.99
      ),
      names = FALSE
    )

    cat(
      "Eligible profiles            : ",
      nrow(primary),
      "\n",
      sep = ""
    )

    cat(
      "q_eff median                 : ",
      signif(
        median(
          primary$q_eff
        ),
        10
      ),
      "\n",
      sep = ""
    )

    cat(
      "q_eff Q01                    : ",
      signif(
        qr[1],
        10
      ),
      "\n",
      sep = ""
    )

    cat(
      "q_eff Q99                    : ",
      signif(
        qr[2],
        10
      ),
      "\n",
      sep = ""
    )

    cat(
      "|q_eff - 1| median           : ",
      signif(
        median(
          abs_err
        ),
        8
      ),
      "\n",
      sep = ""
    )

    cat(
      "|q_eff - 1| Q95              : ",
      signif(
        quantile(
          abs_err,
          0.95,
          names = FALSE
        ),
        8
      ),
      "\n",
      sep = ""
    )

    cat(
      "|q_eff - 1| Q99              : ",
      signif(
        quantile(
          abs_err,
          0.99,
          names = FALSE
        ),
        8
      ),
      "\n",
      sep = ""
    )

    cat(
      "% within 0.1%               : ",
      round(
        mean(
          abs_err <= 0.001
        ) * 100,
        4
      ),
      "%\n",
      sep = ""
    )

    cat(
      "% within 1%                 : ",
      round(
        mean(
          abs_err <= 0.01
        ) * 100,
        4
      ),
      "%\n",
      sep = ""
    )

    cat(
      "% within 5%                 : ",
      round(
        mean(
          abs_err <= 0.05
        ) * 100,
        4
      ),
      "%\n",
      sep = ""
    )
  }


  # ============================================================
  # Per-cell-type CORE12
  # ==============================================================

  cat(
    "\n",
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  cat(
    "CORE12 PER-CELL-TYPE q_eff | Z mass >= 1000\n"
  )

  cat(
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )


  for (ct in CORE12) {

    d <- core12_long[
      core12_long$cell_type == ct &
        !core12_long$technical_exclusion &
        core12_long$posterior_mass >= 1000 &
        is.finite(
          core12_long$q_eff
        ),
      ,
      drop = FALSE
    ]

    if (nrow(d) == 0) {

      cat(
        sprintf(
          "%-24s n=0\n",
          ct
        )
      )

      next
    }

    abs_err <- abs(
      d$q_eff -
        1
    )

    cat(
      sprintf(
        paste0(
          "%-24s ",
          "n=%4d  ",
          "median=%.8f  ",
          "Q01=%.8f  ",
          "Q99=%.8f  ",
          "abs.err.Q95=%.6g\n"
        ),

        ct,

        nrow(d),

        median(
          d$q_eff
        ),

        quantile(
          d$q_eff,
          0.01,
          names = FALSE
        ),

        quantile(
          d$q_eff,
          0.99,
          names = FALSE
        ),

        quantile(
          abs_err,
          0.95,
          names = FALSE
        )
      )
    )
  }


  # ============================================================
  # Stage2 technical-zero result
  # ==============================================================

  cat(
    "\n",
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  cat(
    "STAGE2 TECHNICAL-ZERO AUDIT\n"
  )

  cat(
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )


  zero_sample_table <- unique(
    stage2_long[
      stage2_long$
        R2_technical_effective_zero,
      c(
        "sample_id",
        "R2_raw",
        "R2_eff",
        "technical_exclusion"
      ),
      drop = FALSE
    ]
  )

  if (nrow(zero_sample_table) == 0) {

    cat(
      "No technical sample with R_eff <= 0 detected.\n"
    )

  } else {

    print(
      zero_sample_table,
      row.names = FALSE
    )
  }


  # ============================================================
  # G-like check
  # ==============================================================

  cat(
    "\n",
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  cat(
    "EXISTING CORE12 G-LIKE FORMULA AUDIT\n"
  )

  cat(
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  print(
    glike_check,
    row.names = FALSE
  )


  # ============================================================
  # Stage2 pseudo-mixture check
  # ==============================================================

  cat(
    "\n",
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  cat(
    "STAGE2 PSEUDO-MIXTURE CONSISTENCY\n"
  )

  cat(
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  print(
    pseudomix_check,
    row.names = FALSE
  )


  # ============================================================
  # Session information
  # ==============================================================

  capture.output(
    sessionInfo(),
    file = file.path(
      OUT,
      "sessionInfo.txt"
    )
  )


  # ============================================================
  # Finish
  # ==============================================================

  cat(
    "\n",
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  cat(
    "AUDIT COMPLETE\n"
  )

  cat(
    paste0(
      rep("=", 120),
      collapse = ""
    ),
    "\n",
    sep = ""
  )

  cat(
    "Output directory:\n",
    OUT,
    "\n\n",
    sep = ""
  )

  cat(
    "Key files:\n"
  )

  cat(
    "  q_summary_main.tsv\n"
  )

  cat(
    "  core12_qnt_long.tsv\n"
  )

  cat(
    "  theta0_reconstruction_error.tsv\n"
  )

  cat(
    "  depth_effective_vs_raw_summary.tsv\n"
  )

  cat(
    "  stage2_pseudomixture_consistency.tsv\n"
  )

  cat(
    "  stage2_technical_zero_forensic.tsv\n"
  )

  cat(
    "  core12_glike_formula_spotcheck.tsv\n"
  )

  cat(
    "\nEND: ",
    format(Sys.time()),
    "\n",
    sep = ""
  )
}


main()