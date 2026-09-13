#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
15_build_reference_availability_audit.py

Final-reference availability audit before balanced sampling.

Purpose
-------
Freeze the current 20-class reference ontology and describe the
ACTUAL available cell distribution across:

    dataset
    donor
    biological sample
    identity
    disease
    site
    inflammation

This script DOES NOT:
    - sample cells
    - determine donor caps
    - force equal dataset contribution
    - read expression matrices
    - modify previous audit outputs

This script DOES:
    - inherit canonical audited metadata
    - promote Glial/Tuft into the final availability universe
    - reproduce prior cross-dataset identity-count checkpoints
    - quantify donor and dataset concentration
    - generate tables needed to design balanced sampling
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd


# ============================================================
# 0. Paths
# ============================================================

ROOT = Path("/home/mazekai/IBD_EcoTyper")

OUT = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "02_reference_availability"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ------------------------------------------------------------
# Established harmonization outputs
# ------------------------------------------------------------

HARMONIZATION = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "00_cross_dataset_harmonization"
)

CORE_SUMMARY = (
    HARMONIZATION
    / "core_identity_cross_dataset_summary.tsv"
)

CANDIDATE_SUMMARY = (
    HARMONIZATION
    / "candidate_identity_cross_dataset_summary.tsv"
)


# ------------------------------------------------------------
# Canonical audited cell metadata
# ------------------------------------------------------------

GSE_META = (
    ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "04_identity_coverage_v1"
    / "mapped_candidate_cells_v1.parquet"
)

SCP1884_META = (
    ROOT
    / "03_reference"
    / "SCP1884"
    / "audit_v1"
    / "output"
    / "02_candidate_taxonomy"
    / "scp1884_candidate_cells_provisional_taxonomy.parquet"
)

SCP259_META = (
    ROOT
    / "03_reference"
    / "SCP259"
    / "audit_v1"
    / "output"
    / "01_candidate_taxonomy_coverage"
    / "scp259_candidate_cells_provisional_taxonomy.parquet"
)


# ============================================================
# 1. Final ontology
# ============================================================

CORE_IDENTITIES = [
    "B",
    "Plasma",
    "CD4_T",
    "CD8_T",
    "NK",
    "ILC",
    "Monocyte_Macrophage",
    "DC",
    "Mast",
    "Fibroblast",
    "Pericyte",
    "Endothelial",
    "Ileal_absorptive",
    "Colonic_absorptive",
    "Goblet",
    "Stem_TA_progenitor",
    "Paneth",
    "Enteroendocrine",
]

PROMOTED_IDENTITIES = [
    "Glial",
    "Tuft",
]

FINAL_IDENTITIES = (
    CORE_IDENTITIES
    + PROMOTED_IDENTITIES
)

TISSUE_RESTRICTED = {
    "Ileal_absorptive",
    "Paneth",
}

DATASETS = [
    "GSE282122",
    "SCP1884",
    "SCP259",
]


# ============================================================
# 2. Utilities
# ============================================================

def fail(msg):
    print(
        f"\nERROR: {msg}",
        file=sys.stderr,
    )
    sys.exit(1)


def require_file(path):
    if not Path(path).exists():
        fail(
            f"Required file not found:\n{path}"
        )


def string_col(x):
    return (
        x.astype("string")
        .str.strip()
    )


def donor_hhi(counts):
    """
    Herfindahl-Hirschman Index of donor contribution.
    """
    counts = np.asarray(
        counts,
        dtype=float,
    )

    total = counts.sum()

    if total <= 0:
        return np.nan

    p = counts / total

    return float(
        np.sum(p ** 2)
    )


def concentration_metrics(counts):
    """
    Basic concentration metrics from positive cell counts.
    """

    x = np.asarray(
        counts,
        dtype=float,
    )

    if len(x) == 0:
        return {
            "min": 0,
            "q25": 0,
            "median": 0,
            "q75": 0,
            "q90": 0,
            "max": 0,
            "max_fraction": np.nan,
            "HHI": np.nan,
        }

    return {
        "min":
            int(np.min(x)),

        "q25":
            float(np.quantile(x, 0.25)),

        "median":
            float(np.median(x)),

        "q75":
            float(np.quantile(x, 0.75)),

        "q90":
            float(np.quantile(x, 0.90)),

        "max":
            int(np.max(x)),

        "max_fraction":
            float(
                np.max(x)
                / np.sum(x)
            ),

        "HHI":
            donor_hhi(x),
    }


# ============================================================
# 3. Required files
# ============================================================

for path in [
    CORE_SUMMARY,
    CANDIDATE_SUMMARY,
    GSE_META,
    SCP1884_META,
    SCP259_META,
]:
    require_file(path)


print("=" * 78)
print("FINAL REFERENCE AVAILABILITY AUDIT")
print("=" * 78)


# ============================================================
# 4. Read previous harmonization checkpoints
# ============================================================

core_summary = pd.read_csv(
    CORE_SUMMARY,
    sep="\t",
)

candidate_summary = pd.read_csv(
    CANDIDATE_SUMMARY,
    sep="\t",
)


if set(
    core_summary["identity"]
) != set(
    CORE_IDENTITIES
):
    fail(
        "18-class core ontology does not "
        "match harmonization output."
    )


candidate_rows = (
    candidate_summary.loc[
        candidate_summary[
            "identity"
        ].isin(
            PROMOTED_IDENTITIES
        )
    ]
)


if set(
    candidate_rows["identity"]
) != set(
    PROMOTED_IDENTITIES
):
    fail(
        "Glial/Tuft missing from "
        "cross-dataset harmonization."
    )


print(
    "PASS: prior cross-dataset "
    "harmonization outputs loaded."
)


# ============================================================
# 5. Freeze ontology table
# ============================================================

ontology_rows = []

for i, identity in enumerate(
    FINAL_IDENTITIES,
    start=1,
):

    ontology_rows.append(
        {
            "ontology_order": i,
            "reference_identity":
                identity,

            "final_reference_action":
                "KEEP_FINAL",

            "ontology_origin":
                (
                    "PROMOTED_CANDIDATE"
                    if identity
                    in PROMOTED_IDENTITIES
                    else
                    "ESTABLISHED_CORE"
                ),

            "source_action_required":
                (
                    "HOLD_CROSS_DATASET"
                    if identity
                    in PROMOTED_IDENTITIES
                    else
                    "KEEP"
                ),

            "tissue_restricted":
                identity
                in TISSUE_RESTRICTED,
        }
    )


ontology = pd.DataFrame(
    ontology_rows
)


ontology.to_csv(
    OUT
    / "final_reference_ontology.tsv",
    sep="\t",
    index=False,
)


print(
    "PASS: 20-class ontology "
    "recorded for availability audit."
)


# ============================================================
# 6. Read and standardize three datasets
# ============================================================

def read_gse282122():

    x = pd.read_parquet(
        GSE_META,
        columns=[
            "cell_id",
            "sample_id",
            "Patient",
            "Disease",
            "Site",
            "Inflammation",
            "coarse_identity",
            "primary_reference_action",
        ],
    )

    return pd.DataFrame(
        {
            "dataset":
                "GSE282122",

            "cell_id":
                string_col(
                    x["cell_id"]
                ),

            "donor_id":
                string_col(
                    x["Patient"]
                ),

            "sample_id":
                string_col(
                    x["sample_id"]
                ),

            "disease":
                string_col(
                    x["Disease"]
                ),

            "site":
                string_col(
                    x["Site"]
                ),

            "inflammation":
                string_col(
                    x["Inflammation"]
                ),

            "identity":
                string_col(
                    x["coarse_identity"]
                ),

            "source_action":
                string_col(
                    x[
                        "primary_reference_action"
                    ]
                ),
        }
    )


def read_scp1884():

    x = pd.read_parquet(
        SCP1884_META,
        columns=[
            "cell_id",
            "v2_donor_id",
            "biological_sample_id",
            "disease",
            "site_canonical",
            "inflammation",
            "provisional_coarse_identity",
            "primary_reference_action",
        ],
    )

    return pd.DataFrame(
        {
            "dataset":
                "SCP1884",

            "cell_id":
                string_col(
                    x["cell_id"]
                ),

            "donor_id":
                string_col(
                    x["v2_donor_id"]
                ),

            "sample_id":
                string_col(
                    x[
                        "biological_sample_id"
                    ]
                ),

            "disease":
                string_col(
                    x["disease"]
                ),

            "site":
                string_col(
                    x[
                        "site_canonical"
                    ]
                ),

            "inflammation":
                string_col(
                    x["inflammation"]
                ),

            "identity":
                string_col(
                    x[
                        "provisional_coarse_identity"
                    ]
                ),

            "source_action":
                string_col(
                    x[
                        "primary_reference_action"
                    ]
                ),
        }
    )


def read_scp259():

    x = pd.read_parquet(
        SCP259_META,
        columns=[
            "cell_id",
            "donor_id",
            "sample_id",
            "disease",
            "site",
            "inflammation",
            "coarse_identity",
            "primary_reference_action",
        ],
    )

    return pd.DataFrame(
        {
            "dataset":
                "SCP259",

            "cell_id":
                string_col(
                    x["cell_id"]
                ),

            "donor_id":
                string_col(
                    x["donor_id"]
                ),

            "sample_id":
                string_col(
                    x["sample_id"]
                ),

            "disease":
                string_col(
                    x["disease"]
                ),

            "site":
                string_col(
                    x["site"]
                ),

            "inflammation":
                string_col(
                    x["inflammation"]
                ),

            "identity":
                string_col(
                    x["coarse_identity"]
                ),

            "source_action":
                string_col(
                    x[
                        "primary_reference_action"
                    ]
                ),
        }
    )


tables = []

for name, func in [
    (
        "GSE282122",
        read_gse282122,
    ),
    (
        "SCP1884",
        read_scp1884,
    ),
    (
        "SCP259",
        read_scp259,
    ),
]:

    x = func()

    print(
        f"{name}: "
        f"{len(x):,} audited candidate cells"
    )

    tables.append(x)


cells = pd.concat(
    tables,
    ignore_index=True,
)


# ============================================================
# 7. Structural hard checks
# ============================================================

for col in [
    "cell_id",
    "donor_id",
    "sample_id",
    "identity",
    "source_action",
]:

    missing = (
        cells[col].isna()
        |
        (
            cells[col]
            .astype("string")
            .str.len()
            == 0
        )
    )

    if missing.any():
        fail(
            f"{col}: "
            f"{int(missing.sum()):,} "
            "missing values."
        )


dup = cells.duplicated(
    subset=[
        "dataset",
        "cell_id",
    ],
    keep=False,
)

if dup.any():

    fail(
        "Duplicate dataset × cell_id "
        "detected."
    )


cells["global_cell_id"] = (
    cells["dataset"]
    .astype(str)
    + "::"
    + cells["cell_id"]
    .astype(str)
)


if cells[
    "global_cell_id"
].duplicated().any():

    fail(
        "global_cell_id is not unique."
    )


print(
    "PASS: cell/donor/sample structural "
    "checks."
)


# ============================================================
# 8. Construct final eligible universe
# ============================================================

core_mask = (
    cells["identity"].isin(
        CORE_IDENTITIES
    )
    &
    (
        cells["source_action"]
        == "KEEP"
    )
)


promoted_mask = (
    cells["identity"].isin(
        PROMOTED_IDENTITIES
    )
    &
    (
        cells["source_action"]
        == "HOLD_CROSS_DATASET"
    )
)


final_cells = (
    cells.loc[
        core_mask
        | promoted_mask
    ]
    .copy()
)


final_cells[
    "reference_identity"
] = final_cells["identity"]


final_cells[
    "final_reference_action"
] = "KEEP_FINAL"


if set(
    final_cells[
        "reference_identity"
    ].unique()
) != set(
    FINAL_IDENTITIES
):

    missing = (
        set(FINAL_IDENTITIES)
        -
        set(
            final_cells[
                "reference_identity"
            ].unique()
        )
    )

    fail(
        "Final universe missing identities: "
        f"{sorted(missing)}"
    )


print(
    f"Final 20-class eligible universe: "
    f"{len(final_cells):,} cells"
)


# ============================================================
# 9. HARD checkpoint against previous harmonization
# ============================================================

expected = {}


for _, row in (
    core_summary.iterrows()
):

    identity = row["identity"]

    for dataset in DATASETS:

        expected[
            (
                dataset,
                identity,
            )
        ] = int(
            row[
                f"{dataset}__n_cells"
            ]
        )


for _, row in (
    candidate_rows.iterrows()
):

    identity = row["identity"]

    for dataset in DATASETS:

        expected[
            (
                dataset,
                identity,
            )
        ] = int(
            row[
                f"{dataset}__n_cells"
            ]
        )


observed = (
    final_cells
    .groupby(
        [
            "dataset",
            "reference_identity",
        ],
        observed=True,
    )
    .size()
)


checkpoint_rows = []


for identity in FINAL_IDENTITIES:

    for dataset in DATASETS:

        obs = int(
            observed.get(
                (
                    dataset,
                    identity,
                ),
                0,
            )
        )

        exp = int(
            expected.get(
                (
                    dataset,
                    identity,
                ),
                0,
            )
        )


        checkpoint_rows.append(
            {
                "dataset":
                    dataset,

                "reference_identity":
                    identity,

                "expected_n_cells":
                    exp,

                "observed_n_cells":
                    obs,

                "difference":
                    obs - exp,

                "status":
                    (
                        "PASS"
                        if obs == exp
                        else "FAIL"
                    ),
            }
        )


checkpoint = pd.DataFrame(
    checkpoint_rows
)


checkpoint.to_csv(
    OUT
    / "identity_count_checkpoint.tsv",
    sep="\t",
    index=False,
)


if (
    checkpoint["status"]
    != "PASS"
).any():

    bad = checkpoint.loc[
        checkpoint["status"]
        != "PASS"
    ]

    print(
        bad.to_string(
            index=False
        )
    )

    fail(
        "Final eligible universe does not "
        "reproduce prior harmonization counts."
    )


print(
    "PASS: all dataset × identity cell "
    "counts reproduce prior audit."
)


# ============================================================
# 10. Positive donor × identity availability
# ============================================================

donor_identity = (
    final_cells
    .groupby(
        [
            "dataset",
            "donor_id",
            "reference_identity",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),

        n_samples=(
            "sample_id",
            "nunique",
        ),

        n_sites=(
            "site",
            "nunique",
        ),

        n_disease_labels=(
            "disease",
            "nunique",
        ),

        n_inflammation_labels=(
            "inflammation",
            "nunique",
        ),
    )
    .reset_index()
)


donor_identity.to_csv(
    OUT
    / "dataset_donor_identity_availability_positive.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 11. COMPLETE donor × identity grid, including zeros
# ============================================================

all_grid_parts = []


for dataset in DATASETS:

    donors = sorted(
        final_cells.loc[
            final_cells["dataset"]
            == dataset,
            "donor_id",
        ].unique()
    )


    grid = pd.MultiIndex.from_product(
        [
            [dataset],
            donors,
            FINAL_IDENTITIES,
        ],
        names=[
            "dataset",
            "donor_id",
            "reference_identity",
        ],
    ).to_frame(
        index=False
    )


    all_grid_parts.append(
        grid
    )


full_grid = pd.concat(
    all_grid_parts,
    ignore_index=True,
)


full_grid = full_grid.merge(
    donor_identity,
    on=[
        "dataset",
        "donor_id",
        "reference_identity",
    ],
    how="left",
    validate="one_to_one",
)


for col in [
    "n_cells",
    "n_samples",
    "n_sites",
    "n_disease_labels",
    "n_inflammation_labels",
]:

    full_grid[col] = (
        full_grid[col]
        .fillna(0)
        .astype(int)
    )


full_grid[
    "identity_present"
] = (
    full_grid["n_cells"]
    > 0
)


full_grid.to_csv(
    OUT
    / "dataset_donor_identity_availability_complete.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 12. Wide donor × identity matrix
# ============================================================

wide = (
    full_grid
    .pivot(
        index=[
            "dataset",
            "donor_id",
        ],
        columns="reference_identity",
        values="n_cells",
    )
    .reset_index()
)


# Preserve ontology ordering.

ordered_cols = [
    "dataset",
    "donor_id",
] + FINAL_IDENTITIES


wide = wide[
    ordered_cols
]


wide.to_csv(
    OUT
    / "dataset_donor_identity_matrix.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 13. Donor × SAMPLE × identity availability
# ============================================================

sample_identity = (
    final_cells
    .groupby(
        [
            "dataset",
            "donor_id",
            "sample_id",
            "reference_identity",
            "disease",
            "site",
            "inflammation",
        ],
        dropna=False,
        observed=True,
    )
    .size()
    .rename(
        "n_cells"
    )
    .reset_index()
)


sample_identity.to_csv(
    OUT
    / "dataset_donor_sample_identity_availability.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 14. Dataset × identity donor-balance summary
# ============================================================

dataset_identity_rows = []


for dataset in DATASETS:

    for identity in FINAL_IDENTITIES:

        x = donor_identity.loc[
            (
                donor_identity["dataset"]
                == dataset
            )
            &
            (
                donor_identity[
                    "reference_identity"
                ]
                == identity
            ),
            "n_cells",
        ]


        total = int(
            x.sum()
        )


        metrics = (
            concentration_metrics(
                x
            )
        )


        dataset_identity_rows.append(
            {
                "dataset":
                    dataset,

                "reference_identity":
                    identity,

                "tissue_restricted":
                    identity
                    in TISSUE_RESTRICTED,

                "n_cells":
                    total,

                "n_donors_with_cells":
                    int(len(x)),

                "min_cells_per_positive_donor":
                    metrics["min"],

                "q25_cells_per_positive_donor":
                    metrics["q25"],

                "median_cells_per_positive_donor":
                    metrics["median"],

                "q75_cells_per_positive_donor":
                    metrics["q75"],

                "q90_cells_per_positive_donor":
                    metrics["q90"],

                "max_cells_per_positive_donor":
                    metrics["max"],

                "max_donor_fraction":
                    metrics[
                        "max_fraction"
                    ],

                "donor_HHI":
                    metrics["HHI"],
            }
        )


dataset_identity_summary = pd.DataFrame(
    dataset_identity_rows
)


dataset_identity_summary.to_csv(
    OUT
    / "dataset_identity_donor_balance_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 15. Overall identity availability across datasets
# ============================================================

identity_rows = []


for identity in FINAL_IDENTITIES:

    sub = final_cells.loc[
        final_cells[
            "reference_identity"
        ]
        == identity
    ]


    donor_counts = (
        sub
        .groupby(
            [
                "dataset",
                "donor_id",
            ],
            observed=True,
        )
        .size()
    )


    dataset_counts = (
        sub
        .groupby(
            "dataset",
            observed=True,
        )
        .size()
    )


    total = len(
        sub
    )


    donor_metrics = (
        concentration_metrics(
            donor_counts
        )
    )


    if total > 0:

        dataset_fractions = (
            dataset_counts
            / total
        )

        max_dataset_fraction = float(
            dataset_fractions.max()
        )

        dominant_dataset = (
            dataset_fractions.idxmax()
        )

    else:

        max_dataset_fraction = np.nan
        dominant_dataset = "NA"


    identity_rows.append(
        {
            "reference_identity":
                identity,

            "tissue_restricted":
                identity
                in TISSUE_RESTRICTED,

            "n_cells_total":
                total,

            "n_datasets_with_cells":
                int(
                    dataset_counts.size
                ),

            "n_dataset_donor_units":
                int(
                    donor_counts.size
                ),

            "n_unique_samples":
                int(
                    sub[
                        [
                            "dataset",
                            "sample_id",
                        ]
                    ]
                    .drop_duplicates()
                    .shape[0]
                ),

            "median_cells_per_dataset_donor":
                donor_metrics[
                    "median"
                ],

            "q75_cells_per_dataset_donor":
                donor_metrics[
                    "q75"
                ],

            "q90_cells_per_dataset_donor":
                donor_metrics[
                    "q90"
                ],

            "max_cells_per_dataset_donor":
                donor_metrics[
                    "max"
                ],

            "max_dataset_donor_fraction":
                donor_metrics[
                    "max_fraction"
                ],

            "dataset_donor_HHI":
                donor_metrics[
                    "HHI"
                ],

            "dominant_dataset":
                dominant_dataset,

            "max_dataset_fraction":
                max_dataset_fraction,
        }
    )


identity_summary = pd.DataFrame(
    identity_rows
)


ontology_order = {
    x: i
    for i, x
    in enumerate(
        FINAL_IDENTITIES,
        start=1,
    )
}


identity_summary[
    "ontology_order"
] = (
    identity_summary[
        "reference_identity"
    ]
    .map(
        ontology_order
    )
)


identity_summary = (
    identity_summary
    .sort_values(
        "ontology_order"
    )
)


identity_summary.to_csv(
    OUT
    / "identity_availability_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 16. Dataset contribution by identity
# ============================================================

dataset_contribution = (
    final_cells
    .groupby(
        [
            "reference_identity",
            "dataset",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),

        n_donors=(
            "donor_id",
            "nunique",
        ),

        n_samples=(
            "sample_id",
            "nunique",
        ),
    )
    .reset_index()
)


identity_totals = (
    dataset_contribution
    .groupby(
        "reference_identity"
    )["n_cells"]
    .transform(
        "sum"
    )
)


dataset_contribution[
    "fraction_within_identity"
] = (
    dataset_contribution[
        "n_cells"
    ]
    /
    identity_totals
)


dataset_contribution.to_csv(
    OUT
    / "identity_dataset_contribution.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 17. Disease / site / inflammation coverage
# ============================================================

for field in [
    "disease",
    "site",
    "inflammation",
]:

    summary = (
        final_cells
        .groupby(
            [
                "reference_identity",
                "dataset",
                field,
            ],
            dropna=False,
            observed=True,
        )
        .agg(
            n_cells=(
                "cell_id",
                "size",
            ),

            n_donors=(
                "donor_id",
                "nunique",
            ),

            n_samples=(
                "sample_id",
                "nunique",
            ),
        )
        .reset_index()
    )


    summary.to_csv(
        OUT
        / (
            "identity_dataset_"
            f"{field}_coverage.tsv"
        ),
        sep="\t",
        index=False,
    )


# ============================================================
# 18. Dataset-level overall summary
# ============================================================

dataset_summary = (
    final_cells
    .groupby(
        "dataset",
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),

        n_donors=(
            "donor_id",
            "nunique",
        ),

        n_samples=(
            "sample_id",
            "nunique",
        ),

        n_identities=(
            "reference_identity",
            "nunique",
        ),
    )
    .reset_index()
)


dataset_summary[
    "fraction_final_universe"
] = (
    dataset_summary[
        "n_cells"
    ]
    /
    len(final_cells)
)


dataset_summary.to_csv(
    OUT
    / "dataset_availability_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 19. Human-readable report
# ============================================================

report = (
    OUT
    / "reference_availability_audit_summary.txt"
)


with report.open(
    "w",
    encoding="utf-8",
) as f:

    f.write(
        "IBD EcoTyper combined reference\n"
    )

    f.write(
        "Pre-sampling availability audit\n"
    )

    f.write(
        "=" * 78
        + "\n\n"
    )


    f.write(
        "Final ontology\n"
    )

    f.write(
        "-" * 78
        + "\n"
    )

    f.write(
        "20 identities\n"
    )

    f.write(
        "18 established core identities\n"
    )

    f.write(
        "Glial + Tuft promoted into the "
        "final reference universe\n\n"
    )


    f.write(
        "Important\n"
    )

    f.write(
        "-" * 78
        + "\n"
    )

    f.write(
        "No cell sampling was performed.\n"
    )

    f.write(
        "No donor cap was imposed.\n"
    )

    f.write(
        "No dataset balancing was imposed.\n"
    )

    f.write(
        "All numbers represent the actual "
        "availability inherited from the "
        "canonical audits.\n\n"
    )


    f.write(
        "Dataset summary\n"
    )

    f.write(
        "-" * 78
        + "\n"
    )

    f.write(
        dataset_summary.to_string(
            index=False
        )
    )

    f.write(
        "\n\n"
    )


    f.write(
        "Identity availability\n"
    )

    f.write(
        "-" * 78
        + "\n"
    )

    f.write(
        identity_summary.to_string(
            index=False
        )
    )

    f.write("\n")


# ============================================================
# 20. Console output
# ============================================================

print()
print("=" * 78)
print("DATASET AVAILABILITY")
print("=" * 78)

print(
    dataset_summary.to_string(
        index=False
    )
)


print()
print("=" * 78)
print("IDENTITY AVAILABILITY")
print("=" * 78)

print(
    identity_summary[
        [
            "reference_identity",
            "n_cells_total",
            "n_datasets_with_cells",
            "n_dataset_donor_units",
            "median_cells_per_dataset_donor",
            "q75_cells_per_dataset_donor",
            "q90_cells_per_dataset_donor",
            "max_cells_per_dataset_donor",
            "max_dataset_donor_fraction",
            "dominant_dataset",
            "max_dataset_fraction",
        ]
    ].to_string(
        index=False
    )
)


print()
print("=" * 78)
print("DATASET × IDENTITY DONOR BALANCE")
print("=" * 78)

print(
    dataset_identity_summary[
        [
            "dataset",
            "reference_identity",
            "n_cells",
            "n_donors_with_cells",
            "median_cells_per_positive_donor",
            "q75_cells_per_positive_donor",
            "q90_cells_per_positive_donor",
            "max_cells_per_positive_donor",
            "max_donor_fraction",
        ]
    ].to_string(
        index=False
    )
)


print()
print("=" * 78)
print("OUTPUT DIRECTORY")
print("=" * 78)

print(OUT)


print()
print("=" * 78)
print("SUCCESS")
print("=" * 78)

print(
    "Final 20-class reference availability "
    "audit completed."
)

print(
    "No sampling rule has been applied."
)

print(
    "Use these distributions to design "
    "the balanced sampling strategy."
)