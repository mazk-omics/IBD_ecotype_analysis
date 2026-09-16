suppressPackageStartupMessages({
    library(BayesPrism)
})

options(stringsAsFactors = FALSE)

project_root <- "/home/mazekai/IBD_EcoTyper"

bp_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/14_bayesprism_broad_run"
)

qc_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/15_bayesprism_broad14_sanity_audit"
)

expression_dir <- file.path(
    bp_dir,
    "celltype_expression"
)

prism_rds <- file.path(
    bp_dir,
    "bayesprism_broad14_prism_object.rds"
)

sample_qc_file <- file.path(
    qc_dir,
    "broad14_sample_level_qc.tsv"
)

out_dir <- file.path(
    project_root,
    "03_reference/combined_reference/output/16_bayesprism_broad14_neutrophil_allocation_audit"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

cat("============================================================\n")
cat("BROAD14 NEUTROPHIL POSTERIOR ALLOCATION AUDIT\n")
cat("============================================================\n")

# ============================================================
# 1. DEFINITIONS
# ============================================================

cell_types <- c(
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

# Relatively specific neutrophil / granulocyte markers.
primary_neutrophil <- c(
    "FCGR3B",
    "CSF3R",
    "CXCR1",
    "CXCR2",
    "FPR1"
)

# Important inflammatory granulocyte-associated genes,
# but substantially less cell-type specific.
secondary_neutrophil <- c(
    "S100A8",
    "S100A9"
)

# Internal controls:
# these help confirm that Z allocation behaves biologically
# as expected for lineages represented in the reference.
control_genes <- c(
    "PTPRC",
    "LST1",
    "TYROBP",
    "C1QC",
    "CD3D",
    "CD3E",
    "MS4A1",
    "EPCAM"
)

genes_requested <- unique(
    c(
        primary_neutrophil,
        secondary_neutrophil,
        control_genes
    )
)

gene_class <- c(
    setNames(
        rep(
            "primary_neutrophil",
            length(primary_neutrophil)
        ),
        primary_neutrophil
    ),
    setNames(
        rep(
            "secondary_neutrophil",
            length(secondary_neutrophil)
        ),
        secondary_neutrophil
    ),
    setNames(
        rep(
            "control",
            length(control_genes)
        ),
        control_genes
    )
)

# ============================================================
# 2. LOAD PRISM + SAMPLE METADATA
# ============================================================

cat("\n[1] LOAD PRISM + SAMPLE QC\n")

stopifnot(
    file.exists(prism_rds),
    file.exists(sample_qc_file),
    dir.exists(expression_dir)
)

prism_obj <- readRDS(
    prism_rds
)

mixture <- prism_obj@mixture

stopifnot(
    nrow(mixture) == 2490L,
    ncol(mixture) == 13668L
)

qc <- read.delim(
    sample_qc_file,
    check.names = FALSE
)

stopifnot(
    nrow(qc) == 2490L,
    "sample_id" %in% colnames(qc),
    "typere" %in% colnames(qc),
    "ibd_disease" %in% colnames(qc)
)

qc_idx <- match(
    rownames(mixture),
    qc$sample_id
)

if (anyNA(qc_idx)) {
    stop(
        "Sample QC metadata does not align with prism mixture."
    )
}

qc <- qc[
    qc_idx,
    ,
    drop = FALSE
]

stopifnot(
    identical(
        qc$sample_id,
        rownames(mixture)
    )
)

cat(
    "Prism mixture:",
    nrow(mixture),
    "samples x",
    ncol(mixture),
    "genes\n"
)

cat("Sample metadata alignment: PASS\n")

# ============================================================
# 3. CHECK TARGET GENES
# ============================================================

cat("\n[2] CHECK TARGET GENES\n")

genes_present <- intersect(
    genes_requested,
    colnames(mixture)
)

genes_missing <- setdiff(
    genes_requested,
    genes_present
)

cat(
    "Requested genes:",
    length(genes_requested),
    "\n"
)

cat(
    "Present in prism gene space:",
    length(genes_present),
    "\n"
)

if (length(genes_missing) > 0L) {

    cat(
        "Missing genes:",
        paste(
            genes_missing,
            collapse = ", "
        ),
        "\n"
    )
}

if (
    length(
        intersect(
            primary_neutrophil,
            genes_present
        )
    ) == 0L
) {
    stop(
        "No primary neutrophil markers remain in prism gene space."
    )
}

cat(
    "Genes analyzed:",
    paste(
        genes_present,
        collapse = ", "
    ),
    "\n"
)

# ============================================================
# 4. LOAD ONLY SELECTED GENES FROM EACH Z MATRIX
#
# Do NOT hold all 14 full 2490 x 13668 matrices in memory.
# Each RDS is loaded, subset immediately, then released.
# ============================================================

cat("\n[3] LOAD SELECTED GENES FROM 14 POSTERIOR Z MATRICES\n")

n_samples <- nrow(mixture)
n_genes <- length(genes_present)
n_types <- length(cell_types)

z_array <- array(
    0,
    dim = c(
        n_samples,
        n_genes,
        n_types
    ),
    dimnames = list(
        rownames(mixture),
        genes_present,
        cell_types
    )
)

for (i in seq_along(cell_types)) {

    ct <- cell_types[i]

    safe_ct <- gsub(
        "[^A-Za-z0-9_.-]",
        "_",
        ct
    )

    z_file <- file.path(
        expression_dir,
        paste0(
            "bayesprism_expression_",
            safe_ct,
            ".rds"
        )
    )

    if (!file.exists(z_file)) {
        stop(
            "Missing Z matrix: ",
            z_file
        )
    }

    cat(
        "Loading:",
        ct,
        "\n"
    )

    z <- readRDS(
        z_file
    )

    if (!is.matrix(z)) {
        z <- as.matrix(z)
    }

    stopifnot(
        identical(
            dim(z),
            c(2490L, 13668L)
        ),
        identical(
            rownames(z),
            rownames(mixture)
        ),
        identical(
            colnames(z),
            colnames(mixture)
        )
    )

    z_array[
        ,
        ,
        i
    ] <- z[
        ,
        genes_present,
        drop = FALSE
    ]

    rm(z)

    gc(verbose = FALSE)
}

cat("Selected-gene Z extraction: PASS\n")

# ============================================================
# 5. CONSERVATION CHECK
#
# For each sample/gene:
#
# sum(cell-type Z) should recover the corresponding posterior
# allocated read mass for the mixture gene.
#
# This is checked descriptively rather than imposing an
# arbitrary biological threshold.
# ============================================================

cat("\n[4] Z READ-MASS CONSERVATION\n")

z_total_by_sample_gene <- apply(
    z_array,
    c(1L, 2L),
    sum
)

bulk_selected <- mixture[
    ,
    genes_present,
    drop = FALSE
]

conservation <- do.call(
    rbind,
    lapply(
        genes_present,
        function(g) {

            z_total <- z_total_by_sample_gene[
                ,
                g
            ]

            bulk_total <- bulk_selected[
                ,
                g
            ]

            total_bulk_mass <- sum(
                bulk_total
            )

            total_z_mass <- sum(
                z_total
            )

            data.frame(
                gene = g,
                marker_class =
                    unname(
                        gene_class[g]
                    ),
                bulk_total_mass =
                    total_bulk_mass,
                z_total_mass =
                    total_z_mass,
                total_Z_to_bulk_ratio =
                    if (
                        total_bulk_mass > 0
                    ) {
                        total_z_mass /
                            total_bulk_mass
                    } else {
                        NA_real_
                    },
                samplewise_pearson =
                    if (
                        sd(z_total) > 0 &&
                        sd(bulk_total) > 0
                    ) {
                        cor(
                            z_total,
                            bulk_total,
                            method = "pearson"
                        )
                    } else {
                        NA_real_
                    },
                mean_absolute_difference =
                    mean(
                        abs(
                            z_total -
                            bulk_total
                        )
                    ),
                stringsAsFactors = FALSE
            )
        }
    )
)

print(
    conservation,
    row.names = FALSE
)

write.table(
    conservation,
    file.path(
        out_dir,
        "selected_gene_Z_conservation.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. WHO GETS EACH GENE?
#
# A[t,g] =
#   sum_n Z[n,g,t] /
#   sum_t sum_n Z[n,g,t]
# ============================================================

cat("\n[5] OVERALL POSTERIOR READ ALLOCATION\n")

overall_mass <- matrix(
    0,
    nrow = n_genes,
    ncol = n_types,
    dimnames = list(
        genes_present,
        cell_types
    )
)

for (i in seq_along(cell_types)) {

    overall_mass[
        ,
        i
    ] <- colSums(
        z_array[
            ,
            ,
            i,
            drop = FALSE
        ]
    )
}

gene_total_mass <- rowSums(
    overall_mass
)

overall_fraction <- sweep(
    overall_mass,
    1L,
    gene_total_mass,
    "/"
)

overall_fraction[
    !is.finite(overall_fraction)
] <- NA_real_

allocation_long <- do.call(
    rbind,
    lapply(
        genes_present,
        function(g) {

            data.frame(
                gene = g,
                marker_class =
                    unname(
                        gene_class[g]
                    ),
                cell_type =
                    cell_types,
                posterior_Z_mass =
                    overall_mass[
                        g,
                        cell_types
                    ],
                allocation_fraction =
                    overall_fraction[
                        g,
                        cell_types
                    ],
                stringsAsFactors = FALSE
            )
        }
    )
)

allocation_long <- allocation_long[
    order(
        allocation_long$gene,
        -allocation_long$allocation_fraction
    ),
    ,
    drop = FALSE
]

write.table(
    allocation_long,
    file.path(
        out_dir,
        "selected_gene_allocation_long.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. COMPACT PER-GENE SUMMARY
# ============================================================

cat("\n[6] COMPACT GENE SUMMARY\n")

gene_summary_list <- lapply(
    genes_present,
    function(g) {

        f <- overall_fraction[
            g,
            ,
            drop = TRUE
        ]

        ord <- order(
            f,
            decreasing = TRUE,
            na.last = TRUE
        )

        data.frame(
            gene = g,
            marker_class =
                unname(
                    gene_class[g]
                ),

            top1_type =
                names(f)[ord[1L]],

            top1_fraction =
                unname(f[ord[1L]]),

            top2_type =
                names(f)[ord[2L]],

            top2_fraction =
                unname(f[ord[2L]]),

            top3_type =
                names(f)[ord[3L]],

            top3_fraction =
                unname(f[ord[3L]]),

            DC_fraction =
                unname(f["DC"]),

            Monocyte_Macrophage_fraction =
                unname(
                    f[
                        "Monocyte_Macrophage"
                    ]
                ),

            Endothelial_fraction =
                unname(
                    f["Endothelial"]
                ),

            CD4_T_fraction =
                unname(
                    f["CD4_T"]
                ),

            Fibroblast_fraction =
                unname(
                    f["Fibroblast"]
                ),

            Epithelial_fraction =
                unname(
                    f["Epithelial"]
                ),

            stringsAsFactors = FALSE
        )
    }
)

gene_summary <- do.call(
    rbind,
    gene_summary_list
)

print(
    gene_summary,
    row.names = FALSE
)

write.table(
    gene_summary,
    file.path(
        out_dir,
        "selected_gene_allocation_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. IBD INFLAMED vs NON-INFLAMED ALLOCATION
#
# Important:
# This asks whether the SAME gene is assigned differently
# across broad cell types in inflamed vs non-inflamed IBD.
# ============================================================

cat("\n[7] IBD INFLAMED vs NON-INFLAMED ALLOCATION\n")

idx_I <- which(
    qc$ibd_disease %in% c(
        "CD",
        "UC"
    ) &
    qc$typere == "I"
)

idx_NonI <- which(
    qc$ibd_disease %in% c(
        "CD",
        "UC"
    ) &
    qc$typere == "NonI"
)

cat(
    "Inflamed IBD:",
    length(idx_I),
    "\n"
)

cat(
    "Non-inflamed IBD:",
    length(idx_NonI),
    "\n"
)

group_allocation <- list()

for (group_name in c("I", "NonI")) {

    idx <- if (group_name == "I") {
        idx_I
    } else {
        idx_NonI
    }

    group_mass <- matrix(
        0,
        nrow = n_genes,
        ncol = n_types,
        dimnames = list(
            genes_present,
            cell_types
        )
    )

    for (i in seq_along(cell_types)) {

        group_mass[, i] <- colSums(
            z_array[
                idx,
                ,
                i,
                drop = FALSE
            ]
        )
    }

    group_total <- rowSums(
        group_mass
    )

    group_frac <- sweep(
        group_mass,
        1L,
        group_total,
        "/"
    )

    group_frac[
        !is.finite(group_frac)
    ] <- NA_real_

    group_allocation[[group_name]] <- group_frac
}

stopifnot(
    all(
        dim(
            group_allocation[["I"]]
        ) == c(
            n_genes,
            n_types
        )
    ),
    all(
        dim(
            group_allocation[["NonI"]]
        ) == c(
            n_genes,
            n_types
        )
    )
)

group_long <- do.call(
    rbind,
    lapply(
        genes_present,
        function(g) {

            do.call(
                rbind,
                lapply(
                    cell_types,
                    function(ct) {

                        frac_I <-
                            group_allocation[["I"]][
                                g,
                                ct
                            ]

                        frac_NonI <-
                            group_allocation[["NonI"]][
                                g,
                                ct
                            ]

                        data.frame(
                            gene = g,
                            marker_class =
                                unname(
                                    gene_class[g]
                                ),
                            cell_type = ct,
                            allocation_fraction_I =
                                frac_I,
                            allocation_fraction_NonI =
                                frac_NonI,
                            difference_I_minus_NonI =
                                frac_I -
                                frac_NonI,
                            stringsAsFactors = FALSE
                        )
                    }
                )
            )
        }
    )
)

stopifnot(
    nrow(group_long) ==
        length(genes_present) *
        length(cell_types)
)

write.table(
    group_long,
    file.path(
        out_dir,
        "IBD_inflammation_stratified_gene_allocation.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "Inflammation-stratified allocation table:",
    nrow(group_long),
    "rows\n"
)

cat("Inflammation-stratified allocation: PASS\n")

# ============================================================
# 9. PRIMARY NEUTROPHIL MARKER REPORT
# ============================================================

cat("\n[8] PRIMARY NEUTROPHIL MARKER REPORT\n")

primary_present <- intersect(
    primary_neutrophil,
    genes_present
)

primary_report <- gene_summary[
    gene_summary$gene %in%
        primary_present,
    ,
    drop = FALSE
]

primary_report <- primary_report[
    match(
        primary_present,
        primary_report$gene
    ),
    ,
    drop = FALSE
]

print(
    primary_report,
    row.names = FALSE
)

# Aggregate all primary neutrophil marker Z mass together.
primary_mass_by_type <- colSums(
    overall_mass[
        primary_present,
        ,
        drop = FALSE
    ]
)

primary_fraction_by_type <-
    primary_mass_by_type /
    sum(primary_mass_by_type)

primary_aggregate <- data.frame(
    cell_type =
        names(
            primary_fraction_by_type
        ),
    posterior_Z_mass =
        unname(
            primary_mass_by_type
        ),
    allocation_fraction =
        unname(
            primary_fraction_by_type
        ),
    stringsAsFactors = FALSE
)

primary_aggregate <- primary_aggregate[
    order(
        primary_aggregate$
            allocation_fraction,
        decreasing = TRUE
    ),
    ,
    drop = FALSE
]

cat(
    "\nAggregate allocation of primary neutrophil-marker Z mass:\n"
)

print(
    primary_aggregate,
    row.names = FALSE
)

write.table(
    primary_report,
    file.path(
        out_dir,
        "primary_neutrophil_gene_allocation_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    primary_aggregate,
    file.path(
        out_dir,
        "primary_neutrophil_marker_aggregate_allocation.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 10. INTERNAL CONTROL REPORT
# ============================================================

cat("\n[9] INTERNAL CONTROL GENES\n")

control_present <- intersect(
    control_genes,
    genes_present
)

control_report <- gene_summary[
    gene_summary$gene %in%
        control_present,
    ,
    drop = FALSE
]

control_report <- control_report[
    match(
        control_present,
        control_report$gene
    ),
    ,
    drop = FALSE
]

print(
    control_report,
    row.names = FALSE
)

# ============================================================
# 11. FINAL OUTPUT SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("NEUTROPHIL POSTERIOR ALLOCATION AUDIT COMPLETE\n")
cat("============================================================\n")

cat("\nPrimary interpretation files:\n")

cat(
    "  ",
    file.path(
        out_dir,
        "primary_neutrophil_gene_allocation_summary.tsv"
    ),
    "\n"
)

cat(
    "  ",
    file.path(
        out_dir,
        "primary_neutrophil_marker_aggregate_allocation.tsv"
    ),
    "\n"
)

cat(
    "  ",
    file.path(
        out_dir,
        "selected_gene_Z_conservation.tsv"
    ),
    "\n"
)

cat(
    "  ",
    file.path(
        out_dir,
        "IBD_inflammation_stratified_gene_allocation.tsv"
    ),
    "\n"
)

cat("\nDONE\n")
