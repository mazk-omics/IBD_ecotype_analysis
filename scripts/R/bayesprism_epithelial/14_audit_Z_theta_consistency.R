#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages(
  library(BayesPrism)
)

ROOT <- "/home/mazekai/IBD_EcoTyper"

S1_RESULT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/14_bayesprism_broad_run",
  "bayesprism_broad14_result.rds"
)

S1_Z_DIR <- file.path(
  ROOT,
  "03_reference/combined_reference/output/14_bayesprism_broad_run",
  "celltype_expression"
)

S2_RESULT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/20_bayesprism_epithelial_full_run",
  "bayesprism_epithelial7_result.rds"
)

S2_Z_DIR <- file.path(
  ROOT,
  "03_reference/combined_reference/output/20_bayesprism_epithelial_full_run",
  "celltype_expression"
)

OUT <- file.path(
  ROOT,
  "03_reference/combined_reference/output/25_Z_theta_consistency_audit"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

bp1 <- readRDS(S1_RESULT)
bp2 <- readRDS(S2_RESULT)

theta1 <- bp1@posterior.initial.cellType@theta
theta2 <- bp2@posterior.initial.cellType@theta


audit_stage <- function(
  theta,
  z_dir,
  stage_name
) {

  types <- colnames(theta)
  sample_ids <- rownames(theta)

  masses <- matrix(
    0,
    nrow = nrow(theta),
    ncol = length(types),
    dimnames = list(
      sample_ids,
      types
    )
  )

  for (ct in types) {

    path <- file.path(
      z_dir,
      paste0(
        "bayesprism_expression_",
        ct,
        ".rds"
      )
    )

    if (!file.exists(path)) {
      stop(
        stage_name,
        ": missing Z file for ",
        ct
      )
    }

    z <- readRDS(path)

    stopifnot(
      setequal(
        rownames(z),
        sample_ids
      )
    )

    z <- z[
      sample_ids,
      ,
      drop = FALSE
    ]

    masses[, ct] <- rowSums(z)

    rm(z)
    gc()
  }

  total_mass <- rowSums(masses)

  valid <- total_mass > 0

  mass_fraction <- matrix(
    NA_real_,
    nrow = nrow(masses),
    ncol = ncol(masses),
    dimnames = dimnames(masses)
  )

  mass_fraction[valid, ] <-
    masses[valid, , drop = FALSE] /
    total_mass[valid]

  rows <- list()

  for (i in seq_along(types)) {

    ct <- types[i]

    x <- mass_fraction[valid, ct]
    y <- theta[valid, ct]

    rows[[i]] <- data.frame(
      stage = stage_name,
      cell_type = ct,
      n_valid = sum(valid),

      pearson = cor(
        x,
        y,
        method = "pearson"
      ),

      spearman = cor(
        x,
        y,
        method = "spearman"
      ),

      MAE = mean(
        abs(x - y)
      ),

      max_abs_diff = max(
        abs(x - y)
      ),

      median_Z_mass = median(
        masses[, ct]
      ),

      zero_mass_samples = sum(
        masses[, ct] <= 0
      ),

      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(
    rbind,
    rows
  )

  list(
    summary = summary,
    masses = masses,
    total_mass = total_mass,
    valid = valid
  )
}


cat("Auditing Stage1...\n")

a1 <- audit_stage(
  theta1,
  S1_Z_DIR,
  "Stage1"
)

cat("Auditing Stage2...\n")

a2 <- audit_stage(
  theta2,
  S2_Z_DIR,
  "Stage2"
)

summary_all <- rbind(
  a1$summary,
  a2$summary
)

write.table(
  summary_all,
  file.path(
    OUT,
    "Z_theta_consistency_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================\n")
cat("Z / THETA CONSISTENCY\n")
cat("============================================\n\n")

print(
  summary_all,
  row.names = FALSE
)

pass <- all(
  summary_all$MAE < 1e-8
) &&
  all(
    summary_all$max_abs_diff < 1e-6
  )

cat(
  "\nFINAL STATUS:",
  ifelse(
    pass,
    "PASS",
    "REVIEW"
  ),
  "\n"
)

cat(
  "Output:",
  OUT,
  "\n"
)

if (!pass) {
  quit(
    save = "no",
    status = 1
  )
}

