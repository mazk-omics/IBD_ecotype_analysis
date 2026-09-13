#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
18_select_final_reference_cells.py

Final cell selection for the IBD EcoTyper single-cell reference.

FROZEN SAMPLING RULE
--------------------
For each final reference identity i:

    cap_i = min(round(Q75_i), 200)

where Q75_i is calculated from positive:

    dataset × donor × identity

cell availability.

For every dataset × donor × identity unit:

    n_selected = min(n_available, cap_i)

No additional dataset-level equalization is applied.

Within each donor × identity unit:
    - biological samples are represented as evenly as possible;
    - exact cell IDs are then sampled deterministically;
    - random seeds are stable and reproducible.

IMPORTANT
---------
This script:
    - DOES select real cell IDs;
    - DOES generate the final cell manifest;
    - DOES perform post-sampling audit;
    - DOES NOT read expression matrices;
    - DOES NOT harmonize genes;
    - DOES NOT modify upstream audit files.

The resulting selected_cell_manifest.parquet is the handoff
to the subsequent raw-count extraction stage.
"""

from pathlib import Path
import hashlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


# ============================================================
# 0. Configuration
# ============================================================

ROOT = Path(
    "/home/mazekai/IBD_EcoTyper"
)

GLOBAL_SEED = 20260913

MAX_DONOR_CAP = 200


# ------------------------------------------------------------
# Upstream availability audit
# ------------------------------------------------------------

AVAIL_DIR = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "02_reference_availability"
)

AVAIL_FILE = (
    AVAIL_DIR
    / "dataset_donor_identity_availability_positive.tsv"
)

ONTOLOGY_FILE = (
    AVAIL_DIR
    / "final_reference_ontology.tsv"
)


# ------------------------------------------------------------
# Frozen sampling simulation
# ------------------------------------------------------------

SIM_DIR = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "03_sampling_simulation"
)

SIM_CAP_FILE = (
    SIM_DIR
    / "strategy_identity_caps.tsv"
)

SIM_SELECTED_FILE = (
    SIM_DIR
    / "comparison_selected_cells_by_identity.tsv"
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


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUT = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "05_final_cell_selection"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Frozen ontology
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


def clean_string(x):
    return (
        x.astype("string")
        .str.strip()
    )


def stable_seed(*parts):
    """
    Deterministic seed independent of Python's
    session-specific hash implementation.
    """

    text = (
        str(GLOBAL_SEED)
        + "||"
        + "||".join(
            str(x)
            for x in parts
        )
    )

    digest = hashlib.sha256(
        text.encode("utf-8")
    ).hexdigest()

    return int(
        digest[:8],
        16,
    )


def hhi(values):
    x = np.asarray(
        values,
        dtype=float,
    )

    total = x.sum()

    if total <= 0:
        return np.nan

    p = x / total

    return float(
        np.sum(p ** 2)
    )


def effective_n(values):
    value = hhi(values)

    if (
        not np.isfinite(value)
        or value <= 0
    ):
        return np.nan

    return float(
        1.0 / value
    )


def parquet_columns(path):
    return set(
        pq.ParquetFile(
            path
        ).schema.names
    )


# ============================================================
# 3. Required inputs
# ============================================================

for path in [
    AVAIL_FILE,
    ONTOLOGY_FILE,
    SIM_CAP_FILE,
    SIM_SELECTED_FILE,
    GSE_META,
    SCP1884_META,
    SCP259_META,
]:
    require_file(path)


print("=" * 82)
print("FINAL IBD REFERENCE CELL SELECTION")
print("=" * 82)


# ============================================================
# 4. Load frozen ontology
# ============================================================

ontology = pd.read_csv(
    ONTOLOGY_FILE,
    sep="\t",
).sort_values(
    "ontology_order"
)


observed_ontology = (
    ontology[
        "reference_identity"
    ]
    .tolist()
)


if observed_ontology != FINAL_IDENTITIES:
    fail(
        "Final ontology does not exactly match "
        "the frozen 20-class identity order."
    )


print(
    "PASS: frozen 20-class ontology loaded."
)


# ============================================================
# 5. Load availability audit
# ============================================================

availability = pd.read_csv(
    AVAIL_FILE,
    sep="\t",
)


required_avail = {
    "dataset",
    "donor_id",
    "reference_identity",
    "n_cells",
}


missing = (
    required_avail
    - set(
        availability.columns
    )
)


if missing:
    fail(
        "Availability table missing columns: "
        f"{sorted(missing)}"
    )


if availability.duplicated(
    subset=[
        "dataset",
        "donor_id",
        "reference_identity",
    ]
).any():

    fail(
        "Duplicate dataset × donor × identity "
        "availability rows."
    )


if (
    availability[
        "n_cells"
    ] <= 0
).any():

    fail(
        "Positive availability table contains "
        "non-positive counts."
    )


expected_eligible_total = int(
    availability[
        "n_cells"
    ].sum()
)


print(
    f"Availability universe: "
    f"{expected_eligible_total:,} cells"
)


# ============================================================
# 6. Load frozen simulation caps
# ============================================================

sim_caps = pd.read_csv(
    SIM_CAP_FILE,
    sep="\t",
)


sim_caps = sim_caps.loc[
    sim_caps[
        "strategy"
    ]
    == "IDENTITY_Q75_MAX200"
].copy()


if set(
    sim_caps[
        "reference_identity"
    ]
) != set(
    FINAL_IDENTITIES
):

    fail(
        "Simulation cap table does not contain "
        "all 20 final identities."
    )


sim_cap_map = dict(
    zip(
        sim_caps[
            "reference_identity"
        ],
        sim_caps[
            "cap"
        ].astype(int),
    )
)


# ============================================================
# 7. Independently recalculate Q75_MAX200 caps
# ============================================================

recalculated_cap_map = {}


for identity in FINAL_IDENTITIES:

    x = availability.loc[
        availability[
            "reference_identity"
        ]
        == identity,
        "n_cells",
    ]


    if len(x) == 0:
        fail(
            f"No positive availability "
            f"for {identity}."
        )


    q75 = float(
        x.quantile(
            0.75
        )
    )


    rounded_q75 = max(
        1,
        int(
            np.rint(
                q75
            )
        ),
    )


    cap = min(
        rounded_q75,
        MAX_DONOR_CAP,
    )


    recalculated_cap_map[
        identity
    ] = cap


    if (
        cap
        != sim_cap_map[
            identity
        ]
    ):

        fail(
            f"{identity}: recalculated cap "
            f"{cap} != frozen simulation cap "
            f"{sim_cap_map[identity]}"
        )


print(
    "PASS: all identity-specific caps exactly "
    "reproduce Q75_MAX200 simulation."
)


# ============================================================
# 8. Load expected identity counts from simulation
# ============================================================

sim_selected = pd.read_csv(
    SIM_SELECTED_FILE,
    sep="\t",
)


if (
    "IDENTITY_Q75_MAX200"
    not in sim_selected.columns
):

    fail(
        "Simulation selected-cell table lacks "
        "IDENTITY_Q75_MAX200."
    )


expected_identity_selected = dict(
    zip(
        sim_selected[
            "reference_identity"
        ],
        sim_selected[
            "IDENTITY_Q75_MAX200"
        ].astype(int),
    )
)


if set(
    expected_identity_selected.keys()
) != set(
    FINAL_IDENTITIES
):

    fail(
        "Simulation identity count checkpoint "
        "does not match final ontology."
    )


expected_total_selected = int(
    sum(
        expected_identity_selected.values()
    )
)


print(
    f"Frozen expected selected cells: "
    f"{expected_total_selected:,}"
)


# ============================================================
# 9. Metadata readers
# ============================================================

def read_gse282122():

    available_cols = parquet_columns(
        GSE_META
    )

    required = [
        "cell_id",
        "sample_id",
        "Patient",
        "Disease",
        "Site",
        "Inflammation",
        "coarse_identity",
        "primary_reference_action",
    ]

    for col in required:
        if col not in available_cols:
            fail(
                f"GSE282122 metadata missing: {col}"
            )


    optional = []

    if "obs_name" in available_cols:
        optional.append(
            "obs_name"
        )


    x = pd.read_parquet(
        GSE_META,
        columns=(
            required
            + optional
        ),
    )


    source_cell_key = (
        clean_string(
            x["obs_name"]
        )
        if "obs_name" in x.columns
        else
        clean_string(
            x["cell_id"]
        )
    )


    return pd.DataFrame(
        {
            "dataset":
                "GSE282122",

            "cell_id":
                clean_string(
                    x["cell_id"]
                ),

            "source_cell_key":
                source_cell_key,

            "source_bundle":
                "H5AD",

            "donor_id":
                clean_string(
                    x["Patient"]
                ),

            "sample_id":
                clean_string(
                    x["sample_id"]
                ),

            "disease":
                clean_string(
                    x["Disease"]
                ),

            "site":
                clean_string(
                    x["Site"]
                ),

            "inflammation":
                clean_string(
                    x["Inflammation"]
                ),

            "identity":
                clean_string(
                    x["coarse_identity"]
                ),

            "source_action":
                clean_string(
                    x[
                        "primary_reference_action"
                    ]
                ),
        }
    )


def read_scp1884():

    available_cols = parquet_columns(
        SCP1884_META
    )

    required = [
        "cell_id",
        "v2_donor_id",
        "biological_sample_id",
        "disease",
        "site_canonical",
        "inflammation",
        "expression_bundle",
        "provisional_coarse_identity",
        "primary_reference_action",
    ]


    for col in required:
        if col not in available_cols:
            fail(
                f"SCP1884 metadata missing: {col}"
            )


    x = pd.read_parquet(
        SCP1884_META,
        columns=required,
    )


    return pd.DataFrame(
        {
            "dataset":
                "SCP1884",

            "cell_id":
                clean_string(
                    x["cell_id"]
                ),

            "source_cell_key":
                clean_string(
                    x["cell_id"]
                ),

            "source_bundle":
                clean_string(
                    x["expression_bundle"]
                ),

            "donor_id":
                clean_string(
                    x["v2_donor_id"]
                ),

            "sample_id":
                clean_string(
                    x[
                        "biological_sample_id"
                    ]
                ),

            "disease":
                clean_string(
                    x["disease"]
                ),

            "site":
                clean_string(
                    x[
                        "site_canonical"
                    ]
                ),

            "inflammation":
                clean_string(
                    x["inflammation"]
                ),

            "identity":
                clean_string(
                    x[
                        "provisional_coarse_identity"
                    ]
                ),

            "source_action":
                clean_string(
                    x[
                        "primary_reference_action"
                    ]
                ),
        }
    )


def read_scp259():

    available_cols = parquet_columns(
        SCP259_META
    )

    required = [
        "cell_id",
        "donor_id",
        "sample_id",
        "disease",
        "site",
        "inflammation",
        "expression_bundle",
        "coarse_identity",
        "primary_reference_action",
    ]


    for col in required:
        if col not in available_cols:
            fail(
                f"SCP259 metadata missing: {col}"
            )


    x = pd.read_parquet(
        SCP259_META,
        columns=required,
    )


    return pd.DataFrame(
        {
            "dataset":
                "SCP259",

            "cell_id":
                clean_string(
                    x["cell_id"]
                ),

            "source_cell_key":
                clean_string(
                    x["cell_id"]
                ),

            "source_bundle":
                clean_string(
                    x["expression_bundle"]
                ),

            "donor_id":
                clean_string(
                    x["donor_id"]
                ),

            "sample_id":
                clean_string(
                    x["sample_id"]
                ),

            "disease":
                clean_string(
                    x["disease"]
                ),

            "site":
                clean_string(
                    x["site"]
                ),

            "inflammation":
                clean_string(
                    x["inflammation"]
                ),

            "identity":
                clean_string(
                    x["coarse_identity"]
                ),

            "source_action":
                clean_string(
                    x[
                        "primary_reference_action"
                    ]
                ),
        }
    )


# ============================================================
# 10. Read canonical metadata
# ============================================================

print()
print("=" * 82)
print("READING CANONICAL AUDITED CELL METADATA")
print("=" * 82)


parts = []


for dataset, reader in [
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

    x = reader()

    print(
        f"{dataset}: "
        f"{len(x):,} audited candidate cells"
    )

    parts.append(x)


cells = pd.concat(
    parts,
    ignore_index=True,
)


# ============================================================
# 11. Structural metadata audit
# ============================================================

for col in [
    "cell_id",
    "source_cell_key",
    "source_bundle",
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
            f"{int(missing.sum()):,} missing values."
        )


if cells.duplicated(
    subset=[
        "dataset",
        "cell_id",
    ],
    keep=False,
).any():

    fail(
        "Duplicate dataset × cell_id detected."
    )


cells[
    "global_cell_id"
] = (
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
    "PASS: metadata structural audit."
)


# ============================================================
# 12. Reconstruct final eligible universe
# ============================================================

core_mask = (
    cells[
        "identity"
    ].isin(
        CORE_IDENTITIES
    )
    &
    (
        cells[
            "source_action"
        ]
        == "KEEP"
    )
)


promoted_mask = (
    cells[
        "identity"
    ].isin(
        PROMOTED_IDENTITIES
    )
    &
    (
        cells[
            "source_action"
        ]
        == "HOLD_CROSS_DATASET"
    )
)


eligible = (
    cells.loc[
        core_mask
        | promoted_mask
    ]
    .copy()
)


eligible[
    "reference_identity"
] = (
    eligible[
        "identity"
    ]
)


eligible[
    "promotion_status"
] = np.where(
    eligible[
        "reference_identity"
    ].isin(
        PROMOTED_IDENTITIES
    ),
    "PROMOTED_FROM_HOLD",
    "ESTABLISHED_CORE",
)


eligible[
    "final_reference_action"
] = "KEEP_FINAL"


if set(
    eligible[
        "reference_identity"
    ].unique()
) != set(
    FINAL_IDENTITIES
):

    fail(
        "Eligible cell universe does not "
        "contain all final identities."
    )


if len(
    eligible
) != expected_eligible_total:

    fail(
        "Eligible metadata cell total does not "
        "match availability audit: "
        f"{len(eligible):,} != "
        f"{expected_eligible_total:,}"
    )


print(
    f"PASS: final eligible universe "
    f"reproduced: {len(eligible):,} cells."
)


# ============================================================
# 13. Reproduce availability table exactly
# ============================================================

observed_availability = (
    eligible
    .groupby(
        [
            "dataset",
            "donor_id",
            "reference_identity",
        ],
        observed=True,
    )
    .size()
    .rename(
        "observed_n_cells"
    )
    .reset_index()
)


avail_check = (
    availability
    .merge(
        observed_availability,
        on=[
            "dataset",
            "donor_id",
            "reference_identity",
        ],
        how="outer",
        validate="one_to_one",
    )
)


if (
    avail_check[
        "n_cells"
    ].isna().any()
    or
    avail_check[
        "observed_n_cells"
    ].isna().any()
):

    fail(
        "Availability unit set differs from "
        "canonical eligible metadata."
    )


avail_check[
    "observed_n_cells"
] = (
    avail_check[
        "observed_n_cells"
    ].astype(int)
)


if not np.array_equal(
    avail_check[
        "n_cells"
    ].to_numpy(
        dtype=int
    ),
    avail_check[
        "observed_n_cells"
    ].to_numpy(
        dtype=int
    ),
):

    bad = avail_check.loc[
        avail_check[
            "n_cells"
        ]
        !=
        avail_check[
            "observed_n_cells"
        ]
    ]

    print(
        bad.head(50).to_string(
            index=False
        )
    )

    fail(
        "Dataset × donor × identity "
        "availability checkpoint failed."
    )


print(
    "PASS: all availability units exactly "
    "reproduce pre-sampling audit."
)


# ============================================================
# 14. Build exact frozen sampling plan
# ============================================================

plan = availability.copy()


plan[
    "identity_donor_cap"
] = (
    plan[
        "reference_identity"
    ]
    .map(
        recalculated_cap_map
    )
)


plan[
    "n_selected_planned"
] = np.minimum(
    plan[
        "n_cells"
    ].to_numpy(
        dtype=int
    ),
    plan[
        "identity_donor_cap"
    ].to_numpy(
        dtype=int
    ),
)


plan[
    "fraction_retained"
] = (
    plan[
        "n_selected_planned"
    ]
    /
    plan[
        "n_cells"
    ]
)


# ------------------------------------------------------------
# Identity checkpoint against simulation
# ------------------------------------------------------------

planned_identity_counts = (
    plan
    .groupby(
        "reference_identity",
        observed=True,
    )[
        "n_selected_planned"
    ]
    .sum()
)


for identity in FINAL_IDENTITIES:

    observed = int(
        planned_identity_counts.loc[
            identity
        ]
    )

    expected = int(
        expected_identity_selected[
            identity
        ]
    )


    if observed != expected:

        fail(
            f"{identity}: planned selection "
            f"{observed:,} != frozen simulation "
            f"{expected:,}"
        )


planned_total = int(
    plan[
        "n_selected_planned"
    ].sum()
)


if (
    planned_total
    != expected_total_selected
):

    fail(
        "Planned total does not reproduce "
        "simulation."
    )


print(
    f"PASS: exact sampling plan reproduces "
    f"simulation: {planned_total:,} cells."
)


plan.to_csv(
    OUT
    / "dataset_donor_identity_sampling_plan.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 15. Sample-aware quota allocator
# ============================================================

def allocate_sample_quotas(
    sample_counts,
    target,
    seed,
):
    """
    Allocate target cells across biological samples.

    Strategy
    --------
    Round-robin allocation across samples with available cells.

    Consequences
    ------------
    - avoids one biopsy dominating a donor-level quota;
    - preserves as many biological samples as possible;
    - never exceeds sample availability;
    - exact total equals target;
    - deterministic via a stable seed.

    This does NOT change the donor-level cap.
    """

    counts = {
        str(k): int(v)
        for k, v
        in sample_counts.items()
    }


    total_available = int(
        sum(
            counts.values()
        )
    )


    if target > total_available:
        fail(
            "Sample allocation target exceeds "
            "available cells."
        )


    if target == total_available:
        return counts.copy()


    rng = np.random.default_rng(
        seed
    )


    samples = sorted(
        counts.keys()
    )


    rng.shuffle(
        samples
    )


    allocated = {
        s: 0
        for s in samples
    }


    remaining = int(
        target
    )


    while remaining > 0:

        active = [
            s
            for s in samples
            if allocated[s]
            < counts[s]
        ]


        if not active:
            fail(
                "No remaining sample capacity "
                "during quota allocation."
            )


        for sample in active:

            if remaining == 0:
                break

            allocated[
                sample
            ] += 1

            remaining -= 1


    result = {
        s: n
        for s, n
        in allocated.items()
        if n > 0
    }


    if sum(
        result.values()
    ) != target:

        fail(
            "Sample quota allocation total "
            "does not match target."
        )


    return result


# ============================================================
# 16. Deterministic actual cell selection
# ============================================================

print()
print("=" * 82)
print("SELECTING FINAL REFERENCE CELLS")
print("=" * 82)


eligible = (
    eligible
    .sort_values(
        [
            "reference_identity",
            "dataset",
            "donor_id",
            "sample_id",
            "cell_id",
        ]
    )
    .reset_index(
        drop=True
    )
)


selected_parts = []
sample_quota_rows = []


grouped = eligible.groupby(
    [
        "dataset",
        "donor_id",
        "reference_identity",
    ],
    sort=True,
    observed=True,
)


plan_lookup = (
    plan
    .set_index(
        [
            "dataset",
            "donor_id",
            "reference_identity",
        ]
    )
)


for (
    dataset,
    donor,
    identity,
), unit in grouped:

    key = (
        dataset,
        donor,
        identity,
    )


    if key not in (
        plan_lookup.index
    ):
        fail(
            f"Sampling unit absent from plan: "
            f"{key}"
        )


    target = int(
        plan_lookup.loc[
            key,
            "n_selected_planned",
        ]
    )


    n_available = len(
        unit
    )


    expected_available = int(
        plan_lookup.loc[
            key,
            "n_cells",
        ]
    )


    if (
        n_available
        != expected_available
    ):

        fail(
            f"{key}: unit availability changed."
        )


    unit_seed = stable_seed(
        "UNIT",
        dataset,
        donor,
        identity,
    )


    sample_counts = (
        unit[
            "sample_id"
        ]
        .value_counts()
        .to_dict()
    )


    quotas = allocate_sample_quotas(
        sample_counts,
        target,
        unit_seed,
    )


    chosen_parts = []


    for sample_id in sorted(
        quotas.keys()
    ):

        quota = int(
            quotas[
                sample_id
            ]
        )


        sample_pool = (
            unit.loc[
                unit[
                    "sample_id"
                ]
                == sample_id
            ]
            .sort_values(
                "cell_id"
            )
            .reset_index(
                drop=True
            )
        )


        sample_available = len(
            sample_pool
        )


        if quota > sample_available:
            fail(
                f"{key}/{sample_id}: quota "
                "exceeds availability."
            )


        cell_seed = stable_seed(
            "CELL",
            dataset,
            donor,
            identity,
            sample_id,
        )


        rng = np.random.default_rng(
            cell_seed
        )


        if quota == sample_available:

            chosen = (
                sample_pool.copy()
            )

        else:

            indices = rng.choice(
                sample_available,
                size=quota,
                replace=False,
            )


            chosen = (
                sample_pool
                .iloc[
                    np.sort(
                        indices
                    )
                ]
                .copy()
            )


        chosen[
            "unit_n_available"
        ] = n_available

        chosen[
            "identity_donor_cap"
        ] = int(
            recalculated_cap_map[
                identity
            ]
        )

        chosen[
            "unit_n_selected"
        ] = target

        chosen[
            "unit_sampling_seed"
        ] = unit_seed

        chosen[
            "cell_sampling_seed"
        ] = cell_seed


        chosen_parts.append(
            chosen
        )


        sample_quota_rows.append(
            {
                "dataset":
                    dataset,

                "donor_id":
                    donor,

                "reference_identity":
                    identity,

                "sample_id":
                    sample_id,

                "sample_n_available":
                    sample_available,

                "sample_n_selected":
                    quota,

                "unit_n_available":
                    n_available,

                "unit_n_selected":
                    target,

                "identity_donor_cap":
                    int(
                        recalculated_cap_map[
                            identity
                        ]
                    ),

                "unit_sampling_seed":
                    unit_seed,

                "cell_sampling_seed":
                    cell_seed,
            }
        )


    chosen_unit = pd.concat(
        chosen_parts,
        ignore_index=True,
    )


    if len(
        chosen_unit
    ) != target:

        fail(
            f"{key}: selected "
            f"{len(chosen_unit)} != "
            f"target {target}"
        )


    selected_parts.append(
        chosen_unit
    )


selected = pd.concat(
    selected_parts,
    ignore_index=True,
)


sample_quota_table = pd.DataFrame(
    sample_quota_rows
)


print(
    f"Selected raw manifest: "
    f"{len(selected):,} cells"
)


# ============================================================
# 17. Hard post-selection audit
# ============================================================

print()
print("=" * 82)
print("POST-SAMPLING HARD AUDIT")
print("=" * 82)


# ------------------------------------------------------------
# A. Total count
# ------------------------------------------------------------

if len(
    selected
) != expected_total_selected:

    fail(
        "Actual selected total does not "
        "match frozen simulation: "
        f"{len(selected):,} != "
        f"{expected_total_selected:,}"
    )


print(
    f"PASS: total selected cells = "
    f"{len(selected):,}"
)


# ------------------------------------------------------------
# B. Cell uniqueness
# ------------------------------------------------------------

if selected[
    "global_cell_id"
].duplicated().any():

    fail(
        "Duplicate global_cell_id detected "
        "after sampling."
    )


print(
    "PASS: all selected cell IDs are unique."
)


# ------------------------------------------------------------
# C. Selected cells are subset of eligible universe
# ------------------------------------------------------------

eligible_ids = set(
    eligible[
        "global_cell_id"
    ]
)


selected_ids = set(
    selected[
        "global_cell_id"
    ]
)


if not selected_ids.issubset(
    eligible_ids
):

    fail(
        "Selected cells include IDs outside "
        "eligible universe."
    )


print(
    "PASS: selected cells are a strict subset "
    "of audited eligible universe."
)


# ------------------------------------------------------------
# D. All 20 identities present
# ------------------------------------------------------------

if set(
    selected[
        "reference_identity"
    ].unique()
) != set(
    FINAL_IDENTITIES
):

    fail(
        "Selected manifest does not contain "
        "all 20 identities."
    )


print(
    "PASS: all 20 final identities retained."
)


# ------------------------------------------------------------
# E. Exact identity count checkpoint
# ------------------------------------------------------------

actual_identity_counts = (
    selected
    .groupby(
        "reference_identity",
        observed=True,
    )
    .size()
)


for identity in FINAL_IDENTITIES:

    actual = int(
        actual_identity_counts.loc[
            identity
        ]
    )

    expected = int(
        expected_identity_selected[
            identity
        ]
    )


    if actual != expected:

        fail(
            f"{identity}: actual selected "
            f"{actual:,} != expected "
            f"{expected:,}"
        )


print(
    "PASS: every identity exactly reproduces "
    "the frozen simulation count."
)


# ------------------------------------------------------------
# F. Exact dataset × donor × identity plan
# ------------------------------------------------------------

actual_units = (
    selected
    .groupby(
        [
            "dataset",
            "donor_id",
            "reference_identity",
        ],
        observed=True,
    )
    .size()
    .rename(
        "n_selected_actual"
    )
    .reset_index()
)


unit_check = (
    plan
    .merge(
        actual_units,
        on=[
            "dataset",
            "donor_id",
            "reference_identity",
        ],
        how="left",
        validate="one_to_one",
    )
)


unit_check[
    "n_selected_actual"
] = (
    unit_check[
        "n_selected_actual"
    ]
    .fillna(0)
    .astype(int)
)


unit_check[
    "difference"
] = (
    unit_check[
        "n_selected_actual"
    ]
    -
    unit_check[
        "n_selected_planned"
    ]
)


unit_check[
    "status"
] = np.where(
    unit_check[
        "difference"
    ]
    == 0,
    "PASS",
    "FAIL",
)


if (
    unit_check[
        "status"
    ]
    != "PASS"
).any():

    print(
        unit_check.loc[
            unit_check[
                "status"
            ]
            != "PASS"
        ]
        .head(50)
        .to_string(
            index=False
        )
    )

    fail(
        "Actual unit-level sampling differs "
        "from frozen plan."
    )


print(
    "PASS: every dataset × donor × identity "
    "unit exactly matches sampling plan."
)


# ------------------------------------------------------------
# G. Donor cap enforcement
# ------------------------------------------------------------

actual_with_caps = (
    actual_units.copy()
)


actual_with_caps[
    "identity_donor_cap"
] = (
    actual_with_caps[
        "reference_identity"
    ]
    .map(
        recalculated_cap_map
    )
)


if (
    actual_with_caps[
        "n_selected_actual"
    ]
    >
    actual_with_caps[
        "identity_donor_cap"
    ]
).any():

    fail(
        "At least one sampling unit exceeds "
        "its frozen identity donor cap."
    )


print(
    "PASS: all donor-level caps enforced."
)


# ============================================================
# 18. Final manifest
# ============================================================

manifest_cols = [
    "global_cell_id",
    "dataset",
    "cell_id",
    "source_cell_key",
    "source_bundle",
    "donor_id",
    "sample_id",
    "disease",
    "site",
    "inflammation",
    "reference_identity",
    "identity",
    "source_action",
    "promotion_status",
    "final_reference_action",
    "unit_n_available",
    "identity_donor_cap",
    "unit_n_selected",
    "unit_sampling_seed",
    "cell_sampling_seed",
]


selected = (
    selected[
        manifest_cols
    ]
    .sort_values(
        [
            "reference_identity",
            "dataset",
            "donor_id",
            "sample_id",
            "cell_id",
        ]
    )
    .reset_index(
        drop=True
    )
)


selected.to_parquet(
    OUT
    / "selected_cell_manifest.parquet",
    index=False,
)


selected.to_csv(
    OUT
    / "selected_cell_manifest.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)


sample_quota_table.to_csv(
    OUT
    / "sample_level_selection_quotas.tsv",
    sep="\t",
    index=False,
)


unit_check.to_csv(
    OUT
    / "post_sampling_unit_checkpoint.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 19. Post-sampling identity summary
# ============================================================

identity_rows = []


for identity in FINAL_IDENTITIES:

    x = selected.loc[
        selected[
            "reference_identity"
        ]
        == identity
    ]


    donor_units = (
        x
        .groupby(
            [
                "dataset",
                "donor_id",
            ],
            observed=True,
        )
        .size()
    )


    datasets = (
        x
        .groupby(
            "dataset",
            observed=True,
        )
        .size()
    )


    total = len(x)


    identity_rows.append(
        {
            "reference_identity":
                identity,

            "n_selected_cells":
                total,

            "n_datasets":
                int(
                    datasets.size
                ),

            "n_dataset_donor_units":
                int(
                    donor_units.size
                ),

            "n_unique_samples":
                int(
                    x[
                        [
                            "dataset",
                            "sample_id",
                        ]
                    ]
                    .drop_duplicates()
                    .shape[0]
                ),

            "max_dataset_donor_fraction":
                float(
                    donor_units.max()
                    / total
                ),

            "dataset_donor_HHI":
                hhi(
                    donor_units.values
                ),

            "effective_dataset_donor_units":
                effective_n(
                    donor_units.values
                ),

            "dominant_dataset":
                datasets.idxmax(),

            "max_dataset_fraction":
                float(
                    datasets.max()
                    / total
                ),
        }
    )


identity_summary = pd.DataFrame(
    identity_rows
)


identity_summary[
    "ontology_order"
] = (
    identity_summary[
        "reference_identity"
    ].map(
        {
            identity: i
            for i, identity
            in enumerate(
                FINAL_IDENTITIES,
                start=1,
            )
        }
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
    / "post_sampling_identity_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 20. Dataset × identity summary
# ============================================================

dataset_identity = (
    selected
    .groupby(
        [
            "reference_identity",
            "dataset",
        ],
        observed=True,
    )
    .agg(
        n_selected_cells=(
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


totals = (
    dataset_identity
    .groupby(
        "reference_identity"
    )[
        "n_selected_cells"
    ]
    .transform(
        "sum"
    )
)


dataset_identity[
    "fraction_within_identity"
] = (
    dataset_identity[
        "n_selected_cells"
    ]
    /
    totals
)


dataset_identity.to_csv(
    OUT
    / "post_sampling_dataset_identity_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 21. Dataset overall summary
# ============================================================

dataset_summary = (
    selected
    .groupby(
        "dataset",
        observed=True,
    )
    .agg(
        n_selected_cells=(
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
    "fraction_final_reference"
] = (
    dataset_summary[
        "n_selected_cells"
    ]
    /
    len(selected)
)


dataset_summary.to_csv(
    OUT
    / "post_sampling_dataset_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 22. Disease / site / inflammation coverage
# ============================================================

for field in [
    "disease",
    "site",
    "inflammation",
]:

    temp = (
        selected
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
            n_selected_cells=(
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


    temp.to_csv(
        OUT
        / (
            "post_sampling_identity_dataset_"
            f"{field}_coverage.tsv"
        ),
        sep="\t",
        index=False,
    )


# ============================================================
# 23. Source-bundle audit for future extraction
# ============================================================

source_bundle_summary = (
    selected
    .groupby(
        [
            "dataset",
            "source_bundle",
            "reference_identity",
        ],
        observed=True,
    )
    .size()
    .rename(
        "n_selected_cells"
    )
    .reset_index()
)


source_bundle_summary.to_csv(
    OUT
    / "selected_cells_by_source_bundle.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 24. Human-readable summary
# ============================================================

summary_file = (
    OUT
    / "final_cell_selection_summary.txt"
)


with summary_file.open(
    "w",
    encoding="utf-8",
) as f:

    f.write(
        "IBD EcoTyper Final Reference Cell Selection\n"
    )

    f.write(
        "=" * 82
        + "\n\n"
    )


    f.write(
        "Frozen sampling rule\n"
    )

    f.write(
        "-" * 82
        + "\n"
    )

    f.write(
        "cap_i = min(round(Q75_i), 200)\n"
    )

    f.write(
        "Unit = dataset × donor × identity\n"
    )

    f.write(
        "No additional dataset-level "
        "equalization.\n"
    )

    f.write(
        "Within-unit biological samples "
        "balanced by deterministic "
        "round-robin allocation.\n"
    )

    f.write(
        f"Global seed = {GLOBAL_SEED}\n\n"
    )


    f.write(
        "Hard checkpoints\n"
    )

    f.write(
        "-" * 82
        + "\n"
    )

    f.write(
        f"Eligible cells: "
        f"{len(eligible):,}\n"
    )

    f.write(
        f"Selected cells: "
        f"{len(selected):,}\n"
    )

    f.write(
        "20/20 reference identities retained\n"
    )

    f.write(
        "All identity counts reproduce "
        "frozen simulation\n"
    )

    f.write(
        "All dataset×donor×identity units "
        "reproduce frozen plan\n"
    )

    f.write(
        "All selected global cell IDs unique\n"
    )

    f.write(
        "No donor unit exceeds frozen cap\n\n"
    )


    f.write(
        "Identity summary\n"
    )

    f.write(
        "-" * 82
        + "\n"
    )

    f.write(
        identity_summary.to_string(
            index=False
        )
    )

    f.write(
        "\n\n"
    )


    f.write(
        "Dataset summary\n"
    )

    f.write(
        "-" * 82
        + "\n"
    )

    f.write(
        dataset_summary.to_string(
            index=False
        )
    )

    f.write(
        "\n"
    )


# ============================================================
# 25. Console summary
# ============================================================

print()
print("=" * 82)
print("FINAL IDENTITY SUMMARY")
print("=" * 82)

print(
    identity_summary[
        [
            "reference_identity",
            "n_selected_cells",
            "n_datasets",
            "n_dataset_donor_units",
            "n_unique_samples",
            "max_dataset_donor_fraction",
            "effective_dataset_donor_units",
            "dominant_dataset",
            "max_dataset_fraction",
        ]
    ].to_string(
        index=False
    )
)


print()
print("=" * 82)
print("FINAL DATASET SUMMARY")
print("=" * 82)

print(
    dataset_summary.to_string(
        index=False
    )
)


print()
print("=" * 82)
print("SOURCE BUNDLE SUMMARY")
print("=" * 82)

print(
    source_bundle_summary.to_string(
        index=False
    )
)


print()
print("=" * 82)
print("OUTPUT DIRECTORY")
print("=" * 82)

print(OUT)


print()
print("=" * 82)
print("SUCCESS")
print("=" * 82)

print(
    f"Final reference cell manifest created: "
    f"{len(selected):,} cells."
)

print(
    "All frozen sampling-plan checkpoints passed."
)

print(
    "Expression matrices were NOT accessed."
)

print(
    "The next stage is raw-count extraction "
    "using selected_cell_manifest.parquet."
)