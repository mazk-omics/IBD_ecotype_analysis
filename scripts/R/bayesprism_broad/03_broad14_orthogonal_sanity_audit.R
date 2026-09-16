options(stringsAsFactors = FALSE)

project_root <- "/home/mazekai/IBD_EcoTyper"

bp_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/14_bayesprism_broad_run"
)

input_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/11_bayesprism_full_input"
)

out_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/15_bayesprism_broad14_sanity_audit"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

fraction_file <- file.path(
    bp_dir,
    "bayesprism_broad14_fraction_final.tsv"
)

mixture_file <- file.path(
    input_dir,
    "bayesprism_full_msccr_mixture.rds"
)

soft_file <- file.path(
    out_dir,
    "GSE193677_family.soft.gz"
)

cat("============================================================\n")
cat("BROAD14 ORTHOGONAL BIOLOGICAL SANITY AUDIT\n")
cat("============================================================\n")

# ============================================================
# 1. LOAD BAYESPRISM FRACTIONS + BULK
# ============================================================

cat("\n[1] LOAD BROAD14 FRACTIONS + MSCCR BULK\n")

stopifnot(
    file.exists(fraction_file),
    file.exists(mixture_file)
)

frac <- read.delim(
    fraction_file,
    check.names = FALSE
)

stopifnot(
    nrow(frac) == 2490L,
    "sample_id" %in% colnames(frac)
)

mix_obj <- readRDS(mixture_file)

stopifnot(
    is.list(mix_obj),
    "mixture" %in% names(mix_obj)
)

mixture <- mix_obj$mixture

stopifnot(
    is.matrix(mixture),
    identical(
        dim(mixture),
        c(2490L, 53326L)
    )
)

cat(
    "Fractions:",
    nrow(frac),
    "samples x",
    ncol(frac) - 1L,
    "cell types\n"
)

cat(
    "Bulk:",
    nrow(mixture),
    "samples x",
    ncol(mixture),
    "genes\n"
)

# Fraction file and mixture originate from same prism ordering.
mix_idx <- match(
    frac$sample_id,
    rownames(mixture)
)

if (anyNA(mix_idx)) {
    stop(
        "Fraction sample IDs do not map 1:1 to mixture rownames."
    )
}

mixture <- mixture[
    mix_idx,
    ,
    drop = FALSE
]

stopifnot(
    identical(
        rownames(mixture),
        frac$sample_id
    )
)

cat("Fraction/bulk sample alignment: PASS\n")

# ============================================================
# 2. DOWNLOAD OFFICIAL GEO SOFT METADATA
# ============================================================

cat("\n[2] GEO METADATA\n")

if (!file.exists(soft_file)) {

    geo_url <-
        paste0(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/",
            "GSE193nnn/GSE193677/soft/",
            "GSE193677_family.soft.gz"
        )

    cat("Downloading official GSE193677 SOFT metadata...\n")

    download.file(
        geo_url,
        soft_file,
        mode = "wb",
        method = "libcurl",
        quiet = FALSE
    )
}

stopifnot(
    file.exists(soft_file),
    file.info(soft_file)$size > 0
)

soft <- readLines(
    gzfile(soft_file),
    warn = FALSE
)

sample_starts <- grep(
    "^\\^SAMPLE = ",
    soft
)

stopifnot(
    length(sample_starts) == 2490L
)

sample_ends <- c(
    sample_starts[-1L] - 1L,
    length(soft)
)

extract_field <- function(block, prefix) {

    hit <- grep(
        paste0("^", prefix),
        block,
        value = TRUE
    )

    if (length(hit) == 0L) {
        return(NA_character_)
    }

    sub(
        paste0("^", prefix, "\\s*"),
        "",
        hit[1L]
    )
}

extract_characteristic <- function(block, key) {

    x <- grep(
        "^!Sample_characteristics_ch1 = ",
        block,
        value = TRUE
    )

    x <- sub(
        "^!Sample_characteristics_ch1 = ",
        "",
        x
    )

    target <- paste0(
        "^",
        key,
        ":\\s*"
    )

    hit <- grep(
        target,
        x,
        value = TRUE
    )

    if (length(hit) == 0L) {
        return(NA_character_)
    }

    sub(
        target,
        "",
        hit[1L]
    )
}

meta_list <- vector(
    "list",
    length(sample_starts)
)

for (i in seq_along(sample_starts)) {

    block <- soft[
        sample_starts[i]:
        sample_ends[i]
    ]

    gsm <- sub(
        "^\\^SAMPLE = ",
        "",
        block[1L]
    )

    title <- extract_field(
        block,
        "!Sample_title ="
    )

    sample_core <- sub(
        ",.*$",
        "",
        title
    )

    patient_id <- sub(
        "^MSCCR_reGRID_([^_]+)_Biopsy_.*$",
        "\\1",
        sample_core
    )

    meta_list[[i]] <- data.frame(
        gsm = gsm,
        title = title,
        sample_core = sample_core,
        patient_id = patient_id,
        typere = extract_characteristic(
            block,
            "typere"
        ),
        region = extract_characteristic(
            block,
            "regionre"
        ),
        ibd_disease = extract_characteristic(
            block,
            "ibd_disease"
        ),
        disease_type_region = extract_characteristic(
            block,
            "diseasetypere"
        ),
        endoremiss = extract_characteristic(
            block,
            "endoremiss"
        ),
        historemiss = extract_characteristic(
            block,
            "historemiss"
        ),
        stringsAsFactors = FALSE
    )
}

meta <- do.call(
    rbind,
    meta_list
)

stopifnot(
    nrow(meta) == 2490L,
    !anyDuplicated(meta$gsm),
    !anyDuplicated(meta$sample_core)
)

cat("GEO samples parsed:", nrow(meta), "\n")

cat("\nInflammation labels from GEO:\n")
print(
    table(
        meta$typere,
        useNA = "ifany"
    )
)

cat("\nDisease labels:\n")
print(
    table(
        meta$ibd_disease,
        useNA = "ifany"
    )
)

cat("\nRegions:\n")
print(
    sort(
        table(
            meta$region,
            useNA = "ifany"
        ),
        decreasing = TRUE
    )
)

# ============================================================
# 3. MATCH GEO METADATA TO BAYESPRISM SAMPLES
# ============================================================

cat("\n[3] MATCH GEO METADATA TO BAYESPRISM SAMPLES\n")

bp_ids <- frac$sample_id

candidate_ids <- list(
    gsm = meta$gsm,
    sample_core = meta$sample_core,
    title = meta$title
)

match_counts <- vapply(
    candidate_ids,
    function(x) {
        sum(bp_ids %in% x)
    },
    integer(1)
)

print(match_counts)

best_id <- names(
    which.max(match_counts)
)

cat(
    "Best metadata identifier:",
    best_id,
    "\n"
)

cat(
    "Matched samples:",
    max(match_counts),
    "/",
    length(bp_ids),
    "\n"
)

if (max(match_counts) != 2490L) {
    stop(
        paste0(
            "Could not obtain exact 2490/2490 GEO metadata match. ",
            "Do not continue QC."
        )
    )
}

meta_idx <- match(
    bp_ids,
    meta[[best_id]]
)

meta <- meta[
    meta_idx,
    ,
    drop = FALSE
]

stopifnot(
    !anyNA(meta_idx)
)

cat("GEO/BayesPrism sample alignment: PASS\n")

# ============================================================
# 4. DERIVE BROAD BIOLOGICAL COMPARTMENTS
# ============================================================

cat("\n[4] DERIVE BROAD COMPARTMENTS\n")

immune_types <- c(
    "B",
    "CD4_T",
    "CD8_T",
    "DC",
    "ILC",
    "Mast",
    "Monocyte_Macrophage",
    "NK",
    "Plasma"
)

stromal_types <- c(
    "Fibroblast",
    "Endothelial",
    "Pericyte",
    "Glial"
)

stopifnot(
    all(
        c(
            immune_types,
            stromal_types,
            "Epithelial"
        ) %in% colnames(frac)
    )
)

frac$immune_total <- rowSums(
    frac[
        ,
        immune_types,
        drop = FALSE
    ]
)

frac$stromal_total <- rowSums(
    frac[
        ,
        stromal_types,
        drop = FALSE
    ]
)

frac$T_total <-
    frac$CD4_T +
    frac$CD8_T

frac$myeloid_total <-
    frac$DC +
    frac$Monocyte_Macrophage

cat(
    "Whole-cohort immune mean:",
    mean(frac$immune_total),
    "\n"
)

cat(
    "Whole-cohort immune median:",
    median(frac$immune_total),
    "\n"
)

# ============================================================
# 5. BULK MARKER SCORES
#
# Raw MSCCR counts -> CPM -> log2(CPM+1)
# Used ONLY as orthogonal QC, not formal inference.
# ============================================================

cat("\n[5] BUILD BULK MARKER SCORES\n")

marker_sets <- list(

    Immune = c(
        "PTPRC"
    ),

    Epithelial = c(
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19"
    ),

    T_cell = c(
        "CD3D",
        "CD3E",
        "TRBC1"
    ),

    B_cell = c(
        "MS4A1",
        "CD79A"
    ),

    MonoMac = c(
        "LST1",
        "TYROBP",
        "FCER1G",
        "C1QC"
    ),

    Plasma = c(
        "MZB1",
        "JCHAIN",
        "SDC1"
    ),

    Neutrophil = c(
        "FCGR3B",
        "CSF3R",
        "CXCR1",
        "CXCR2",
        "FPR1",
        "S100A8",
        "S100A9"
    )
)

marker_presence <- do.call(
    rbind,
    lapply(
        names(marker_sets),
        function(set_name) {

            genes <- marker_sets[[set_name]]

            data.frame(
                marker_set = set_name,
                gene = genes,
                present = genes %in%
                    colnames(mixture),
                stringsAsFactors = FALSE
            )
        }
    )
)

print(marker_presence)

write.table(
    marker_presence,
    file.path(
        out_dir,
        "marker_gene_presence.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

libsize <- rowSums(
    mixture
)

stopifnot(
    all(
        is.finite(libsize)
    ),
    all(
        libsize > 0
    )
)

marker_scores <- data.frame(
    sample_id = rownames(mixture),
    stringsAsFactors = FALSE
)

for (set_name in names(marker_sets)) {

    genes <- intersect(
        marker_sets[[set_name]],
        colnames(mixture)
    )

    score_name <- paste0(
        set_name,
        "_score"
    )

    if (length(genes) == 0L) {

        warning(
            paste(
                "No genes available for marker set:",
                set_name
            )
        )

        marker_scores[[score_name]] <- NA_real_

        next
    }

    x <- mixture[
        ,
        genes,
        drop = FALSE
    ]

    cpm <- sweep(
        x,
        1L,
        libsize,
        "/"
    ) * 1e6

    log_cpm <- log2(
        cpm + 1
    )

    marker_scores[[score_name]] <- rowMeans(
        log_cpm
    )

    cat(
        set_name,
        ":",
        paste(genes, collapse = ", "),
        "\n"
    )
}

# ============================================================
# 6. MASTER SAMPLE QC TABLE
# ============================================================

cat("\n[6] BUILD MASTER QC TABLE\n")

qc <- cbind(
    meta,
    frac,
    marker_scores[
        ,
        setdiff(
            colnames(marker_scores),
            "sample_id"
        ),
        drop = FALSE
    ]
)

stopifnot(
    nrow(qc) == 2490L
)

write.table(
    qc,
    file.path(
        out_dir,
        "broad14_sample_level_qc.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. ORTHOGONAL MARKER-FRACTION CORRELATIONS
# ============================================================

cat("\n[7] MARKER-FRACTION CORRELATIONS\n")

cor_pairs <- list(

    c(
        "immune_total",
        "Immune_score"
    ),

    c(
        "Epithelial",
        "Epithelial_score"
    ),

    c(
        "T_total",
        "T_cell_score"
    ),

    c(
        "B",
        "B_cell_score"
    ),

    c(
        "Monocyte_Macrophage",
        "MonoMac_score"
    ),

    c(
        "Plasma",
        "Plasma_score"
    )
)

cor_tab <- do.call(
    rbind,
    lapply(
        cor_pairs,
        function(p) {

            x <- qc[[p[1L]]]
            y <- qc[[p[2L]]]

            ok <- is.finite(x) &
                is.finite(y)

            data.frame(
                fraction = p[1L],
                marker_score = p[2L],
                n = sum(ok),
                spearman = cor(
                    x[ok],
                    y[ok],
                    method = "spearman"
                ),
                pearson = cor(
                    x[ok],
                    y[ok],
                    method = "pearson"
                ),
                stringsAsFactors = FALSE
            )
        }
    )
)

print(
    cor_tab,
    row.names = FALSE
)

write.table(
    cor_tab,
    file.path(
        out_dir,
        "marker_fraction_correlations.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. WHERE DOES NEUTROPHIL-LIKE BULK SIGNAL GO?
# ============================================================

cat("\n[8] NEUTROPHIL SIGNAL vs BROAD14 FRACTIONS\n")

broad_types <- c(
    "B",
    "CD4_T",
    "CD8_T",
    "DC",
    "Endothelial",
    "Epithelial",
    "Fibroblast",
    "Glial",
    "ILC",
    "Mast",
    "Monocyte_Macrophage",
    "NK",
    "Pericyte",
    "Plasma"
)

neut_score <- qc$Neutrophil_score

neut_cor <- data.frame(
    cell_type = broad_types,
    spearman = vapply(
        broad_types,
        function(ct) {

            cor(
                qc[[ct]],
                neut_score,
                use = "complete.obs",
                method = "spearman"
            )
        },
        numeric(1)
    ),
    stringsAsFactors = FALSE
)

neut_cor <- neut_cor[
    order(
        neut_cor$spearman,
        decreasing = TRUE
    ),
]

print(
    neut_cor,
    row.names = FALSE
)

write.table(
    neut_cor,
    file.path(
        out_dir,
        "neutrophil_signal_vs_broad14.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 9. INFLAMED vs NON-INFLAMED
#
# Primary descriptive analysis restricted to IBD only.
# Controls are excluded to avoid disease-status confounding.
# ============================================================

cat("\n[9] IBD INFLAMED vs NON-INFLAMED\n")

ibd <- qc[
    qc$ibd_disease %in%
        c("CD", "UC") &
    qc$typere %in%
        c("I", "NonI"),
    ,
    drop = FALSE
]

cat(
    "IBD biopsies with I/NonI labels:",
    nrow(ibd),
    "\n"
)

cat("\nIBD inflammation counts:\n")

print(
    table(
        ibd$typere
    )
)

metrics <- c(
    "immune_total",
    "Epithelial",
    "stromal_total",
    "T_total",
    "B",
    "Monocyte_Macrophage",
    "DC",
    "Plasma",
    "Immune_score",
    "Epithelial_score",
    "T_cell_score",
    "MonoMac_score",
    "Neutrophil_score"
)

group_summary <- do.call(
    rbind,
    lapply(
        metrics,
        function(metric) {

            x_i <- ibd[
                ibd$typere == "I",
                metric
            ]

            x_n <- ibd[
                ibd$typere == "NonI",
                metric
            ]

            data.frame(
                metric = metric,
                inflamed_n =
                    sum(is.finite(x_i)),
                noninflamed_n =
                    sum(is.finite(x_n)),
                inflamed_mean =
                    mean(
                        x_i,
                        na.rm = TRUE
                    ),
                noninflamed_mean =
                    mean(
                        x_n,
                        na.rm = TRUE
                    ),
                inflamed_median =
                    median(
                        x_i,
                        na.rm = TRUE
                    ),
                noninflamed_median =
                    median(
                        x_n,
                        na.rm = TRUE
                    ),
                median_difference_I_minus_NonI =
                    median(
                        x_i,
                        na.rm = TRUE
                    ) -
                    median(
                        x_n,
                        na.rm = TRUE
                    ),
                stringsAsFactors = FALSE
            )
        }
    )
)

print(
    group_summary,
    row.names = FALSE
)

write.table(
    group_summary,
    file.path(
        out_dir,
        "IBD_inflamed_vs_noninflamed_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 10. WITHIN-PATIENT PAIRED SANITY CHECK
#
# Handles repeated biopsies more appropriately:
# mean multiple biopsies within patient/status, then compare
# I - NonI in patients who have both.
# ============================================================

cat("\n[10] WITHIN-PATIENT PAIRED CHECK\n")

paired_results <- list()

for (metric in metrics) {

    tmp <- ibd[
        is.finite(
            ibd[[metric]]
        ),
        c(
            "patient_id",
            "typere",
            metric
        ),
        drop = FALSE
    ]

    agg <- aggregate(
        tmp[[metric]],
        by = list(
            patient_id =
                tmp$patient_id,
            typere =
                tmp$typere
        ),
        FUN = mean
    )

    colnames(agg)[3L] <- "value"

    wide <- reshape(
        agg,
        idvar = "patient_id",
        timevar = "typere",
        direction = "wide"
    )

    if (
        !all(
            c(
                "value.I",
                "value.NonI"
            ) %in% colnames(wide)
        )
    ) {
        next
    }

    wide <- wide[
        is.finite(wide$value.I) &
        is.finite(wide$value.NonI),
        ,
        drop = FALSE
    ]

    delta <-
        wide$value.I -
        wide$value.NonI

    p <- if (length(delta) >= 3L) {

        suppressWarnings(
            wilcox.test(
                wide$value.I,
                wide$value.NonI,
                paired = TRUE,
                exact = FALSE
            )$p.value
        )

    } else {
        NA_real_
    }

    paired_results[[metric]] <- data.frame(
        metric = metric,
        paired_patients = length(delta),
        median_within_patient_difference =
            median(
                delta,
                na.rm = TRUE
            ),
        mean_within_patient_difference =
            mean(
                delta,
                na.rm = TRUE
            ),
        proportion_I_greater_than_NonI =
            mean(
                delta > 0,
                na.rm = TRUE
            ),
        paired_wilcox_p = p,
        stringsAsFactors = FALSE
    )
}

paired_tab <- do.call(
    rbind,
    paired_results
)

print(
    paired_tab,
    row.names = FALSE
)

write.table(
    paired_tab,
    file.path(
        out_dir,
        "IBD_within_patient_paired_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 11. REGION-STRATIFIED DESCRIPTIVE CHECK
# ============================================================

cat("\n[11] REGION-STRATIFIED IMMUNE / EPITHELIAL CHECK\n")

regions_keep <- names(
    which(
        table(
            ibd$region,
            ibd$typere
        )[, "I"] >= 10 &
        table(
            ibd$region,
            ibd$typere
        )[, "NonI"] >= 10
    )
)

region_tab <- do.call(
    rbind,
    lapply(
        regions_keep,
        function(region_name) {

            d <- ibd[
                ibd$region == region_name,
                ,
                drop = FALSE
            ]

            data.frame(
                region = region_name,
                n_I =
                    sum(d$typere == "I"),
                n_NonI =
                    sum(d$typere == "NonI"),

                immune_median_I =
                    median(
                        d$immune_total[
                            d$typere == "I"
                        ],
                        na.rm = TRUE
                    ),

                immune_median_NonI =
                    median(
                        d$immune_total[
                            d$typere == "NonI"
                        ],
                        na.rm = TRUE
                    ),

                epithelial_median_I =
                    median(
                        d$Epithelial[
                            d$typere == "I"
                        ],
                        na.rm = TRUE
                    ),

                epithelial_median_NonI =
                    median(
                        d$Epithelial[
                            d$typere == "NonI"
                        ],
                        na.rm = TRUE
                    ),

                PTPRC_median_I =
                    median(
                        d$Immune_score[
                            d$typere == "I"
                        ],
                        na.rm = TRUE
                    ),

                PTPRC_median_NonI =
                    median(
                        d$Immune_score[
                            d$typere == "NonI"
                        ],
                        na.rm = TRUE
                    ),

                neutrophil_score_median_I =
                    median(
                        d$Neutrophil_score[
                            d$typere == "I"
                        ],
                        na.rm = TRUE
                    ),

                neutrophil_score_median_NonI =
                    median(
                        d$Neutrophil_score[
                            d$typere == "NonI"
                        ],
                        na.rm = TRUE
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)

print(
    region_tab,
    row.names = FALSE
)

write.table(
    region_tab,
    file.path(
        out_dir,
        "IBD_region_stratified_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 12. KEY SANITY FLAGS
#
# These are descriptive QC signals, NOT formal pass/fail rules
# from BayesPrism.
# ============================================================

cat("\n[12] KEY SANITY SIGNALS\n")

get_cor <- function(name) {
    cor_tab$spearman[
        cor_tab$fraction == name
    ]
}

immune_r <- get_cor(
    "immune_total"
)

epi_r <- get_cor(
    "Epithelial"
)

immune_group <- group_summary[
    group_summary$metric ==
        "immune_total",
    ,
    drop = FALSE
]

epi_group <- group_summary[
    group_summary$metric ==
        "Epithelial",
    ,
    drop = FALSE
]

neut_group <- group_summary[
    group_summary$metric ==
        "Neutrophil_score",
    ,
    drop = FALSE
]

cat(
    "Spearman immune_total vs PTPRC:",
    immune_r,
    "\n"
)

cat(
    "Spearman Epithelial vs epithelial marker score:",
    epi_r,
    "\n"
)

cat(
    "Immune median I - NonI:",
    immune_group$
        median_difference_I_minus_NonI,
    "\n"
)

cat(
    "Epithelial median I - NonI:",
    epi_group$
        median_difference_I_minus_NonI,
    "\n"
)

cat(
    "Neutrophil marker median I - NonI:",
    neut_group$
        median_difference_I_minus_NonI,
    "\n"
)

cat("\n============================================================\n")
cat("ORTHOGONAL SANITY AUDIT COMPLETE\n")
cat("============================================================\n")

cat("\nPrimary outputs:\n")
cat(
    "  ",
    file.path(
        out_dir,
        "marker_fraction_correlations.tsv"
    ),
    "\n"
)

cat(
    "  ",
    file.path(
        out_dir,
        "neutrophil_signal_vs_broad14.tsv"
    ),
    "\n"
)

cat(
    "  ",
    file.path(
        out_dir,
        "IBD_inflamed_vs_noninflamed_summary.tsv"
    ),
    "\n"
)

cat(
    "  ",
    file.path(
        out_dir,
        "IBD_within_patient_paired_summary.tsv"
    ),
    "\n"
)

cat(
    "  ",
    file.path(
        out_dir,
        "IBD_region_stratified_summary.tsv"
    ),
    "\n"
)
