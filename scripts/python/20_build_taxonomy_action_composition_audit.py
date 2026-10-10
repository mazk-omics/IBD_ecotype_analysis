#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# Configuration
# ============================================================

ROOT = Path(
    "/home/mazekai/IBD_EcoTyper"
)

OUT = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "00_taxonomy_v2_audit"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ------------------------------------------------------------
# Canonical audited cell-level metadata
# ------------------------------------------------------------

GSE_FILE = (
    ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "04_identity_coverage_v1"
    / "mapped_candidate_cells_v1.parquet"
)

SCP1884_FILE = (
    ROOT
    / "03_reference"
    / "SCP1884"
    / "audit_v1"
    / "output"
    / "02_candidate_taxonomy"
    / "scp1884_candidate_cells_provisional_taxonomy.parquet"
)

SCP259_FILE = (
    ROOT
    / "03_reference"
    / "SCP259"
    / "audit_v1"
    / "output"
    / "01_candidate_taxonomy_coverage"
    / "scp259_candidate_cells_provisional_taxonomy.parquet"
)


EXPECTED_N = {
    "GSE282122": 505_591,
    "SCP1884": 720_633,
    "SCP259": 365_492,
}


CURRENT_REFERENCE20 = {
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
    "Glial",
    "Tuft",
}


CURRENT_ECOTYPER_CORE12 = {
    "B",
    "CD4_T",
    "CD8_T",
    "Endothelial",
    "Enteroendocrine",
    "Fibroblast",
    "Glial",
    "Goblet",
    "ILC",
    "Mast",
    "Monocyte_Macrophage",
    "Pericyte",
}


# ============================================================
# Helpers
# ============================================================

def fail(msg):
    raise RuntimeError(msg)


def require_file(path):

    if not path.is_file():
        fail(
            f"Missing required input:\n{path}"
        )


def pick_col(
    df,
    candidates,
    required=True,
):

    for x in candidates:

        if x in df.columns:
            return x

    if required:
        fail(
            "None of the candidate columns exists: "
            + ", ".join(candidates)
        )

    return None


def clean_string(x):

    return (
        x
        .astype("string")
        .fillna("<NA>")
        .replace("", "<EMPTY>")
    )


def unique_join(x):

    values = sorted(
        set(
            clean_string(x)
        )
    )

    return " | ".join(values)


# ============================================================
# Dataset readers
# ============================================================

def load_gse282122():

    x = pd.read_parquet(
        GSE_FILE
    )

    if len(x) != EXPECTED_N["GSE282122"]:

        fail(
            "Unexpected GSE282122 row count: "
            f"{len(x):,}"
        )

    required = [
        "cell_id",
        "sample_id",
        "Patient",
        "coarse_identity",
        "primary_reference_action",
    ]

    missing = [
        c
        for c in required
        if c not in x.columns
    ]

    if missing:

        fail(
            "GSE282122 missing columns: "
            + ", ".join(missing)
        )

    annotation_col = pick_col(
        x,
        [
            "final_analysis",
            "minor",
        ],
    )

    out = pd.DataFrame({

        "dataset":
            "GSE282122",

        "cell_id":
            clean_string(
                x["cell_id"]
            ),

        "donor_id":
            clean_string(
                x["Patient"]
            ),

        "sample_id":
            clean_string(
                x["sample_id"]
            ),

        "coarse_identity":
            clean_string(
                x["coarse_identity"]
            ),

        "primary_reference_action":
            clean_string(
                x[
                    "primary_reference_action"
                ]
            ),

        "source_annotation":
            clean_string(
                x[annotation_col]
            ),
    })

    return out


def load_scp1884():

    x = pd.read_parquet(
        SCP1884_FILE
    )

    if len(x) != EXPECTED_N["SCP1884"]:

        fail(
            "Unexpected SCP1884 row count: "
            f"{len(x):,}"
        )

    donor_col = pick_col(
        x,
        [
            "v2_donor_id",
            "PubID",
        ],
    )

    required = [
        "cell_id",
        "biological_sample_id",
        "anno2",
        "provisional_coarse_identity",
        "primary_reference_action",
    ]

    missing = [
        c
        for c in required
        if c not in x.columns
    ]

    if missing:

        fail(
            "SCP1884 missing columns: "
            + ", ".join(missing)
        )

    out = pd.DataFrame({

        "dataset":
            "SCP1884",

        "cell_id":
            clean_string(
                x["cell_id"]
            ),

        "donor_id":
            clean_string(
                x[donor_col]
            ),

        "sample_id":
            clean_string(
                x[
                    "biological_sample_id"
                ]
            ),

        "coarse_identity":
            clean_string(
                x[
                    "provisional_coarse_identity"
                ]
            ),

        "primary_reference_action":
            clean_string(
                x[
                    "primary_reference_action"
                ]
            ),

        "source_annotation":
            clean_string(
                x["anno2"]
            ),
    })

    return out


def load_scp259():

    x = pd.read_parquet(
        SCP259_FILE
    )

    if len(x) != EXPECTED_N["SCP259"]:

        fail(
            "Unexpected SCP259 row count: "
            f"{len(x):,}"
        )

    cell_col = pick_col(
        x,
        [
            "cell_id",
            "NAME",
        ],
    )

    annotation_col = pick_col(
        x,
        [
            "author_annotation",
            "Cluster",
        ],
    )

    required = [
        "donor_id",
        "sample_id",
        "coarse_identity",
        "primary_reference_action",
    ]

    missing = [
        c
        for c in required
        if c not in x.columns
    ]

    if missing:

        fail(
            "SCP259 missing columns: "
            + ", ".join(missing)
        )

    out = pd.DataFrame({

        "dataset":
            "SCP259",

        "cell_id":
            clean_string(
                x[cell_col]
            ),

        "donor_id":
            clean_string(
                x["donor_id"]
            ),

        "sample_id":
            clean_string(
                x["sample_id"]
            ),

        "coarse_identity":
            clean_string(
                x["coarse_identity"]
            ),

        "primary_reference_action":
            clean_string(
                x[
                    "primary_reference_action"
                ]
            ),

        "source_annotation":
            clean_string(
                x[annotation_col]
            ),
    })

    return out


# ============================================================
# Input checks
# ============================================================

for f in [
    GSE_FILE,
    SCP1884_FILE,
    SCP259_FILE,
]:

    require_file(f)


print(
    "=" * 78
)

print(
    "IBD SIGNATURE TAXONOMY V2"
)

print(
    "ACTION COMPOSITION AUDIT"
)

print(
    "=" * 78
)


tables = [

    load_gse282122(),
    load_scp1884(),
    load_scp259(),

]


for tab in tables:

    dataset = (
        tab["dataset"]
        .iloc[0]
    )

    print(
        f"{dataset}: "
        f"{len(tab):,} cells"
    )


cells = pd.concat(
    tables,
    axis=0,
    ignore_index=True,
)


# ============================================================
# Global IDs
# ============================================================

cells["global_cell_id"] = (
    cells["dataset"]
    + "::"
    + cells["cell_id"]
)

cells["global_sample_id"] = (
    cells["dataset"]
    + "::"
    + cells["sample_id"]
)

cells["global_donor_id"] = (
    cells["dataset"]
    + "::"
    + cells["donor_id"]
)


if cells[
    "global_cell_id"
].duplicated().any():

    fail(
        "Duplicated global cell IDs."
    )


# ============================================================
# Dataset × identity totals
#
# Denominator for action fractions.
# ============================================================

identity_dataset_total = (

    cells.groupby(
        [
            "dataset",
            "coarse_identity",
        ],
        observed=True,
        dropna=False,
    )

    .agg(

        identity_dataset_n_cells=(
            "global_cell_id",
            "size",
        ),

        identity_dataset_n_samples=(
            "global_sample_id",
            "nunique",
        ),

        identity_dataset_n_donors=(
            "global_donor_id",
            "nunique",
        ),
    )

    .reset_index()
)


# ============================================================
# Dataset × identity × action
# ============================================================

rows = []


for (
    dataset,
    identity,
    action,
), sub in cells.groupby(
    [
        "dataset",
        "coarse_identity",
        "primary_reference_action",
    ],
    observed=True,
    dropna=False,
):

    rows.append({

        "dataset":
            dataset,

        "coarse_identity":
            identity,

        "primary_reference_action":
            action,

        "n_cells":
            len(sub),

        "n_samples":
            sub[
                "global_sample_id"
            ].nunique(),

        "n_donors":
            sub[
                "global_donor_id"
            ].nunique(),

        "source_annotations":
            unique_join(
                sub[
                    "source_annotation"
                ]
            ),
    })


action = pd.DataFrame(
    rows
)


# ============================================================
# Add denominators + fractions
# ============================================================

action = action.merge(

    identity_dataset_total,

    on=[
        "dataset",
        "coarse_identity",
    ],

    how="left",

    validate="many_to_one",
)


action[
    "fraction_within_identity_dataset"
] = (

    action["n_cells"]
    /
    action[
        "identity_dataset_n_cells"
    ]
)


action[
    "fraction_within_identity_dataset"
] = (

    action[
        "fraction_within_identity_dataset"
    ]
    .round(6)
)


action[
    "percent_within_identity_dataset"
] = (

    100
    *
    action[
        "fraction_within_identity_dataset"
    ]
).round(3)


# ============================================================
# Current design flags
# ============================================================

action[
    "current_reference20_member"
] = (

    action[
        "coarse_identity"
    ]
    .isin(
        CURRENT_REFERENCE20
    )
)


action[
    "current_ecotyper_core12_member"
] = (

    action[
        "coarse_identity"
    ]
    .isin(
        CURRENT_ECOTYPER_CORE12
    )
)


# ============================================================
# Sanity check:
# fractions must sum to 1 per dataset × identity
# ============================================================

fraction_check = (

    action.groupby(
        [
            "dataset",
            "coarse_identity",
        ],
        observed=True,
    )[
        "fraction_within_identity_dataset"
    ]

    .sum()

    .reset_index(
        name="fraction_sum"
    )
)


fraction_check[
    "abs_error_from_1"
] = (

    fraction_check[
        "fraction_sum"
    ]
    - 1
).abs()


bad_fraction = (

    fraction_check[
        "abs_error_from_1"
    ]
    > 1e-5
)


if bad_fraction.any():

    bad = fraction_check.loc[
        bad_fraction
    ]

    bad.to_csv(

        OUT
        / "ERROR_action_fraction_not_one.tsv",

        sep="\t",
        index=False,
    )

    fail(
        "Action fractions do not sum to 1 "
        "for some dataset × identity units."
    )


# ============================================================
# Order columns
# ============================================================

column_order = [

    "coarse_identity",
    "dataset",
    "primary_reference_action",

    "current_reference20_member",
    "current_ecotyper_core12_member",

    "n_cells",
    "fraction_within_identity_dataset",
    "percent_within_identity_dataset",

    "n_samples",
    "n_donors",

    "identity_dataset_n_cells",
    "identity_dataset_n_samples",
    "identity_dataset_n_donors",

    "source_annotations",
]


action = action[
    column_order
]


action = (

    action.sort_values(

        [
            "coarse_identity",
            "dataset",
            "fraction_within_identity_dataset",
        ],

        ascending=[
            True,
            True,
            False,
        ],
    )
)


# ============================================================
# Save primary audit table
# ============================================================

OUT_FILE = (

    OUT
    / "06_identity_action_composition_by_dataset.tsv"
)


action.to_csv(

    OUT_FILE,

    sep="\t",
    index=False,
)


# ============================================================
# Build compact global identity summary
#
# Useful for quickly identifying identities with substantial
# excluded/held populations.
# ============================================================

global_action = (

    cells.groupby(
        [
            "coarse_identity",
            "primary_reference_action",
        ],
        observed=True,
        dropna=False,
    )

    .agg(

        n_cells=(
            "global_cell_id",
            "size",
        ),

        n_samples=(
            "global_sample_id",
            "nunique",
        ),

        n_donors=(
            "global_donor_id",
            "nunique",
        ),

        n_datasets=(
            "dataset",
            "nunique",
        ),
    )

    .reset_index()
)


identity_total = (

    cells.groupby(
        "coarse_identity",
        observed=True,
    )

    .size()

    .rename(
        "identity_total_cells"
    )

    .reset_index()
)


global_action = global_action.merge(

    identity_total,

    on="coarse_identity",

    how="left",

    validate="many_to_one",
)


global_action[
    "fraction_within_identity_global"
] = (

    global_action["n_cells"]
    /
    global_action[
        "identity_total_cells"
    ]
)


global_action[
    "percent_within_identity_global"
] = (

    100
    *
    global_action[
        "fraction_within_identity_global"
    ]
).round(3)


global_action = (

    global_action.sort_values(

        [
            "coarse_identity",
            "fraction_within_identity_global",
        ],

        ascending=[
            True,
            False,
        ],
    )
)


global_action.to_csv(

    OUT
    / "07_identity_action_composition_global.tsv",

    sep="\t",
    index=False,
)


# ============================================================
# Focus table:
# identities containing any non-KEEP cells
# ============================================================

nonkeep_identities = set(

    action.loc[
        action[
            "primary_reference_action"
        ]
        != "KEEP",
        "coarse_identity",
    ]
)


focus = (

    action.loc[
        action[
            "coarse_identity"
        ]
        .isin(
            nonkeep_identities
        )
    ]

    .copy()
)


focus.to_csv(

    OUT
    / "08_nonkeep_identity_action_focus.tsv",

    sep="\t",
    index=False,
)


# ============================================================
# Console summary
# ============================================================

print()
print(
    "=" * 78
)

print(
    "ACTION COMPOSITION AUDIT COMPLETE"
)

print(
    "=" * 78
)

print(
    f"Combined candidate cells: "
    f"{len(cells):,}"
)

print(
    f"Dataset × identity × action rows: "
    f"{len(action):,}"
)

print(
    f"Identities with non-KEEP populations: "
    f"{len(nonkeep_identities):,}"
)

print()

print(
    "Primary output:"
)

print(
    OUT_FILE
)

print()

print(
    "Additional outputs:"
)

print(
    OUT
    / "07_identity_action_composition_global.tsv"
)

print(
    OUT
    / "08_nonkeep_identity_action_focus.tsv"
)

print()

print(
    "FINAL STATUS: PASS"
)
