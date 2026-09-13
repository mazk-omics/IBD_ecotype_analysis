#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
16_simulate_reference_sampling_strategies.py

Simulate donor-level sampling strategies for the final 20-class
IBD single-cell reference.

IMPORTANT
---------
This script performs MATHEMATICAL SIMULATION ONLY.

It does NOT:
    - select any cell IDs
    - modify the final ontology
    - modify previous audit outputs
    - read expression matrices

Input
-----
dataset_donor_identity_availability_positive.tsv

Each row represents:
    dataset × donor × identity

with the actual number of available cells.

Strategies
----------
RAW_UNCAPPED
    no sampling; baseline

CAP_100
    max 100 cells / dataset × donor × identity

CAP_200
    max 200 cells / dataset × donor × identity

CAP_300
    max 300 cells / dataset × donor × identity

IDENTITY_Q75
    identity-specific cap = round(Q75 of positive donor-unit availability)

IDENTITY_Q75_MAX200
    identity-specific cap = min(round(Q75), 200)

The goal is NOT to automatically choose a winner.
The goal is to expose tradeoffs between:

    - retained cell number
    - donor dominance
    - dataset composition
    - low-abundance identity preservation

so the final sampling rule can be chosen explicitly.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd


# ============================================================
# 0. Paths
# ============================================================

ROOT = Path(
    "/home/mazekai/IBD_EcoTyper"
)

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

OUT = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "03_sampling_simulation"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Strategies
# ============================================================

STRATEGIES = [
    {
        "strategy": "RAW_UNCAPPED",
        "type": "raw",
    },
    {
        "strategy": "CAP_100",
        "type": "fixed",
        "cap": 100,
    },
    {
        "strategy": "CAP_200",
        "type": "fixed",
        "cap": 200,
    },
    {
        "strategy": "CAP_300",
        "type": "fixed",
        "cap": 300,
    },
    {
        "strategy": "IDENTITY_Q75",
        "type": "q75",
    },
    {
        "strategy": "IDENTITY_Q75_MAX200",
        "type": "q75_max",
        "max_cap": 200,
    },
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


def hhi(counts):
    """
    Herfindahl-Hirschman concentration index.
    """

    x = np.asarray(
        counts,
        dtype=float,
    )

    total = x.sum()

    if total <= 0:
        return np.nan

    p = x / total

    return float(
        np.sum(p ** 2)
    )


def effective_n(counts):
    """
    Effective number of equally contributing units.

        N_eff = 1 / HHI

    Example:
        10 perfectly equal donors -> N_eff = 10
    """

    value = hhi(counts)

    if (
        not np.isfinite(value)
        or value <= 0
    ):
        return np.nan

    return float(
        1.0 / value
    )


# ============================================================
# 3. Load inputs
# ============================================================

require_file(
    AVAIL_FILE
)

require_file(
    ONTOLOGY_FILE
)


availability = pd.read_csv(
    AVAIL_FILE,
    sep="\t",
)

ontology = pd.read_csv(
    ONTOLOGY_FILE,
    sep="\t",
)


required_cols = {
    "dataset",
    "donor_id",
    "reference_identity",
    "n_cells",
}

missing = (
    required_cols
    - set(
        availability.columns
    )
)

if missing:
    fail(
        "Availability table missing columns: "
        f"{sorted(missing)}"
    )


if (
    availability["n_cells"]
    <= 0
).any():
    fail(
        "Positive availability table contains "
        "zero or negative n_cells."
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
        "rows detected."
    )


if (
    "ontology_order"
    not in ontology.columns
    or
    "reference_identity"
    not in ontology.columns
):
    fail(
        "Ontology file lacks required columns."
    )


ontology = ontology.sort_values(
    "ontology_order"
)

IDENTITIES = (
    ontology[
        "reference_identity"
    ]
    .tolist()
)


if len(IDENTITIES) != 20:
    fail(
        f"Expected 20 identities, "
        f"found {len(IDENTITIES)}."
    )


if set(
    availability[
        "reference_identity"
    ]
) != set(
    IDENTITIES
):
    fail(
        "Identity set in availability table "
        "does not match final ontology."
    )


DATASETS = sorted(
    availability[
        "dataset"
    ].unique()
)


raw_total_cells = int(
    availability[
        "n_cells"
    ].sum()
)


print("=" * 80)
print("REFERENCE SAMPLING STRATEGY SIMULATION")
print("=" * 80)
print()
print(
    f"Raw final-reference universe: "
    f"{raw_total_cells:,} cells"
)
print(
    f"Datasets: {', '.join(DATASETS)}"
)
print(
    f"Final identities: {len(IDENTITIES)}"
)


# ============================================================
# 4. Identity-specific availability quantiles
# ============================================================

identity_stats = (
    availability
    .groupby(
        "reference_identity",
        observed=True,
    )["n_cells"]
    .agg(
        n_dataset_donor_units="size",
        total_available="sum",
        median_available="median",
        q75_available=lambda x:
            float(
                x.quantile(0.75)
            ),
        q90_available=lambda x:
            float(
                x.quantile(0.90)
            ),
        max_available="max",
    )
    .reset_index()
)


identity_stats[
    "q75_cap_rounded"
] = (
    identity_stats[
        "q75_available"
    ]
    .apply(
        lambda x:
            max(
                1,
                int(
                    np.rint(x)
                ),
            )
    )
)


q75_cap_map = dict(
    zip(
        identity_stats[
            "reference_identity"
        ],
        identity_stats[
            "q75_cap_rounded"
        ],
    )
)


# ============================================================
# 5. Create strategy-specific caps
# ============================================================

cap_rows = []


for strategy in STRATEGIES:

    name = strategy[
        "strategy"
    ]

    kind = strategy[
        "type"
    ]


    for identity in IDENTITIES:

        q75_cap = int(
            q75_cap_map[
                identity
            ]
        )


        if kind == "raw":

            cap = np.nan
            rule = "NO_CAP"


        elif kind == "fixed":

            cap = int(
                strategy["cap"]
            )

            rule = (
                f"FIXED_{cap}"
            )


        elif kind == "q75":

            cap = q75_cap

            rule = (
                "ROUND_IDENTITY_Q75"
            )


        elif kind == "q75_max":

            cap = min(
                q75_cap,
                int(
                    strategy[
                        "max_cap"
                    ]
                ),
            )

            rule = (
                "MIN_ROUND_IDENTITY_Q75_"
                f"{strategy['max_cap']}"
            )


        else:
            fail(
                f"Unknown strategy type: {kind}"
            )


        cap_rows.append(
            {
                "strategy":
                    name,

                "reference_identity":
                    identity,

                "raw_identity_q75":
                    identity_stats.loc[
                        identity_stats[
                            "reference_identity"
                        ]
                        == identity,
                        "q75_available",
                    ].iloc[0],

                "cap":
                    cap,

                "cap_rule":
                    rule,
            }
        )


cap_table = pd.DataFrame(
    cap_rows
)


cap_table.to_csv(
    OUT
    / "strategy_identity_caps.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 6. Run simulations
# ============================================================

simulated_unit_rows = []

identity_summary_rows = []
dataset_identity_rows = []
dataset_summary_rows = []


for strategy in STRATEGIES:

    name = strategy[
        "strategy"
    ]

    kind = strategy[
        "type"
    ]


    print()
    print(
        f"Simulating: {name}"
    )


    sim = availability.copy()


    if kind == "raw":

        sim["cap"] = np.nan

        sim[
            "n_selected_simulated"
        ] = sim[
            "n_cells"
        ].astype(int)


    else:

        strategy_caps = (
            cap_table.loc[
                cap_table[
                    "strategy"
                ]
                == name,
                [
                    "reference_identity",
                    "cap",
                ],
            ]
        )


        sim = sim.merge(
            strategy_caps,
            on="reference_identity",
            how="left",
            validate="many_to_one",
        )


        if sim[
            "cap"
        ].isna().any():
            fail(
                f"{name}: missing simulated cap."
            )


        sim[
            "n_selected_simulated"
        ] = np.minimum(
            sim[
                "n_cells"
            ].to_numpy(),
            sim[
                "cap"
            ].to_numpy(),
        ).astype(int)


    sim[
        "strategy"
    ] = name


    sim[
        "retained_fraction_unit"
    ] = (
        sim[
            "n_selected_simulated"
        ]
        /
        sim[
            "n_cells"
        ]
    )


    simulated_unit_rows.append(
        sim[
            [
                "strategy",
                "dataset",
                "donor_id",
                "reference_identity",
                "n_cells",
                "cap",
                "n_selected_simulated",
                "retained_fraction_unit",
            ]
        ].copy()
    )


    # ========================================================
    # 6A. Identity-level summary
    # ========================================================

    for identity in IDENTITIES:

        x = sim.loc[
            sim[
                "reference_identity"
            ]
            == identity
        ].copy()


        selected_total = int(
            x[
                "n_selected_simulated"
            ].sum()
        )

        available_total = int(
            x[
                "n_cells"
            ].sum()
        )


        donor_counts = (
            x[
                "n_selected_simulated"
            ]
            .to_numpy()
        )


        dataset_counts = (
            x
            .groupby(
                "dataset",
                observed=True,
            )[
                "n_selected_simulated"
            ]
            .sum()
        )


        dataset_fractions = (
            dataset_counts
            / selected_total
        )


        dominant_dataset = (
            dataset_fractions.idxmax()
        )


        max_dataset_fraction = float(
            dataset_fractions.max()
        )


        identity_summary_rows.append(
            {
                "strategy":
                    name,

                "reference_identity":
                    identity,

                "n_available_cells":
                    available_total,

                "n_selected_simulated":
                    selected_total,

                "fraction_retained":
                    (
                        selected_total
                        / available_total
                    ),

                "n_dataset_donor_units":
                    len(x),

                "max_dataset_donor_fraction":
                    float(
                        donor_counts.max()
                        / selected_total
                    ),

                "dataset_donor_HHI":
                    hhi(
                        donor_counts
                    ),

                "effective_dataset_donor_units":
                    effective_n(
                        donor_counts
                    ),

                "dominant_dataset":
                    dominant_dataset,

                "max_dataset_fraction":
                    max_dataset_fraction,

                "n_datasets_with_cells":
                    int(
                        (
                            dataset_counts
                            > 0
                        ).sum()
                    ),
            }
        )


    # ========================================================
    # 6B. Dataset × identity donor balance
    # ========================================================

    for dataset in DATASETS:

        for identity in IDENTITIES:

            x = sim.loc[
                (
                    sim["dataset"]
                    == dataset
                )
                &
                (
                    sim[
                        "reference_identity"
                    ]
                    == identity
                )
            ].copy()


            if len(x) == 0:

                dataset_identity_rows.append(
                    {
                        "strategy":
                            name,

                        "dataset":
                            dataset,

                        "reference_identity":
                            identity,

                        "n_available_cells":
                            0,

                        "n_selected_simulated":
                            0,

                        "fraction_retained":
                            np.nan,

                        "n_donors":
                            0,

                        "median_selected_per_donor":
                            0,

                        "q75_selected_per_donor":
                            0,

                        "max_selected_per_donor":
                            0,

                        "max_donor_fraction":
                            np.nan,

                        "donor_HHI":
                            np.nan,

                        "effective_donors":
                            np.nan,
                    }
                )

                continue


            selected = (
                x[
                    "n_selected_simulated"
                ]
                .to_numpy()
            )


            selected_total = int(
                selected.sum()
            )

            available_total = int(
                x[
                    "n_cells"
                ].sum()
            )


            dataset_identity_rows.append(
                {
                    "strategy":
                        name,

                    "dataset":
                        dataset,

                    "reference_identity":
                        identity,

                    "n_available_cells":
                        available_total,

                    "n_selected_simulated":
                        selected_total,

                    "fraction_retained":
                        (
                            selected_total
                            / available_total
                        ),

                    "n_donors":
                        len(x),

                    "median_selected_per_donor":
                        float(
                            np.median(
                                selected
                            )
                        ),

                    "q75_selected_per_donor":
                        float(
                            np.quantile(
                                selected,
                                0.75,
                            )
                        ),

                    "max_selected_per_donor":
                        int(
                            selected.max()
                        ),

                    "max_donor_fraction":
                        float(
                            selected.max()
                            / selected_total
                        ),

                    "donor_HHI":
                        hhi(
                            selected
                        ),

                    "effective_donors":
                        effective_n(
                            selected
                        ),
                }
            )


    # ========================================================
    # 6C. Overall dataset contribution
    # ========================================================

    dataset_totals = (
        sim
        .groupby(
            "dataset",
            observed=True,
        )
        .agg(
            n_available_cells=(
                "n_cells",
                "sum",
            ),

            n_selected_simulated=(
                "n_selected_simulated",
                "sum",
            ),
        )
        .reset_index()
    )


    total_selected = int(
        dataset_totals[
            "n_selected_simulated"
        ].sum()
    )


    for _, row in (
        dataset_totals.iterrows()
    ):

        dataset_summary_rows.append(
            {
                "strategy":
                    name,

                "dataset":
                    row[
                        "dataset"
                    ],

                "n_available_cells":
                    int(
                        row[
                            "n_available_cells"
                        ]
                    ),

                "n_selected_simulated":
                    int(
                        row[
                            "n_selected_simulated"
                        ]
                    ),

                "fraction_dataset_cells_retained":
                    (
                        row[
                            "n_selected_simulated"
                        ]
                        /
                        row[
                            "n_available_cells"
                        ]
                    ),

                "fraction_of_simulated_reference":
                    (
                        row[
                            "n_selected_simulated"
                        ]
                        /
                        total_selected
                    ),
            }
        )


# ============================================================
# 7. Combine simulation outputs
# ============================================================

simulated_units = pd.concat(
    simulated_unit_rows,
    ignore_index=True,
)


identity_summary = pd.DataFrame(
    identity_summary_rows
)


dataset_identity_summary = pd.DataFrame(
    dataset_identity_rows
)


dataset_summary = pd.DataFrame(
    dataset_summary_rows
)


# ============================================================
# 8. Add RAW baseline deltas
# ============================================================

raw_identity = (
    identity_summary.loc[
        identity_summary[
            "strategy"
        ]
        == "RAW_UNCAPPED",
        [
            "reference_identity",
            "max_dataset_donor_fraction",
            "dataset_donor_HHI",
            "effective_dataset_donor_units",
            "max_dataset_fraction",
        ],
    ]
    .rename(
        columns={
            "max_dataset_donor_fraction":
                "raw_max_dataset_donor_fraction",

            "dataset_donor_HHI":
                "raw_dataset_donor_HHI",

            "effective_dataset_donor_units":
                "raw_effective_dataset_donor_units",

            "max_dataset_fraction":
                "raw_max_dataset_fraction",
        }
    )
)


identity_summary = (
    identity_summary
    .merge(
        raw_identity,
        on="reference_identity",
        how="left",
        validate="many_to_one",
    )
)


identity_summary[
    "change_max_dataset_donor_fraction"
] = (
    identity_summary[
        "max_dataset_donor_fraction"
    ]
    -
    identity_summary[
        "raw_max_dataset_donor_fraction"
    ]
)


identity_summary[
    "change_dataset_donor_HHI"
] = (
    identity_summary[
        "dataset_donor_HHI"
    ]
    -
    identity_summary[
        "raw_dataset_donor_HHI"
    ]
)


identity_summary[
    "change_effective_dataset_donor_units"
] = (
    identity_summary[
        "effective_dataset_donor_units"
    ]
    -
    identity_summary[
        "raw_effective_dataset_donor_units"
    ]
)


identity_summary[
    "change_max_dataset_fraction"
] = (
    identity_summary[
        "max_dataset_fraction"
    ]
    -
    identity_summary[
        "raw_max_dataset_fraction"
    ]
)


# ============================================================
# 9. Strategy-wide overview
# ============================================================

overview_rows = []


for strategy in [
    x["strategy"]
    for x in STRATEGIES
]:

    identity_sub = (
        identity_summary.loc[
            identity_summary[
                "strategy"
            ]
            == strategy
        ]
    )


    dataset_identity_sub = (
        dataset_identity_summary.loc[
            (
                dataset_identity_summary[
                    "strategy"
                ]
                == strategy
            )
            &
            (
                dataset_identity_summary[
                    "n_selected_simulated"
                ]
                > 0
            )
        ]
    )


    total_selected = int(
        identity_sub[
            "n_selected_simulated"
        ].sum()
    )


    overview_rows.append(
        {
            "strategy":
                strategy,

            "total_selected_cells":
                total_selected,

            "fraction_of_raw_cells_retained":
                (
                    total_selected
                    / raw_total_cells
                ),

            "median_identity_fraction_retained":
                float(
                    identity_sub[
                        "fraction_retained"
                    ].median()
                ),

            "minimum_identity_fraction_retained":
                float(
                    identity_sub[
                        "fraction_retained"
                    ].min()
                ),

            "median_identity_max_dataset_donor_fraction":
                float(
                    identity_sub[
                        "max_dataset_donor_fraction"
                    ].median()
                ),

            "worst_identity_max_dataset_donor_fraction":
                float(
                    identity_sub[
                        "max_dataset_donor_fraction"
                    ].max()
                ),

            "median_effective_dataset_donor_units":
                float(
                    identity_sub[
                        "effective_dataset_donor_units"
                    ].median()
                ),

            "median_identity_max_dataset_fraction":
                float(
                    identity_sub[
                        "max_dataset_fraction"
                    ].median()
                ),

            "worst_identity_max_dataset_fraction":
                float(
                    identity_sub[
                        "max_dataset_fraction"
                    ].max()
                ),

            "median_within_dataset_max_donor_fraction":
                float(
                    dataset_identity_sub[
                        "max_donor_fraction"
                    ].median()
                ),

            "worst_within_dataset_max_donor_fraction":
                float(
                    dataset_identity_sub[
                        "max_donor_fraction"
                    ].max()
                ),

            "n_identities_selected_lt_1000":
                int(
                    (
                        identity_sub[
                            "n_selected_simulated"
                        ]
                        < 1000
                    ).sum()
                ),
        }
    )


strategy_overview = pd.DataFrame(
    overview_rows
)


# ============================================================
# 10. Wide comparison tables
# ============================================================

selected_cells_wide = (
    identity_summary
    .pivot(
        index="reference_identity",
        columns="strategy",
        values="n_selected_simulated",
    )
    .reindex(
        IDENTITIES
    )
    .reset_index()
)


retention_wide = (
    identity_summary
    .pivot(
        index="reference_identity",
        columns="strategy",
        values="fraction_retained",
    )
    .reindex(
        IDENTITIES
    )
    .reset_index()
)


donor_dominance_wide = (
    identity_summary
    .pivot(
        index="reference_identity",
        columns="strategy",
        values="max_dataset_donor_fraction",
    )
    .reindex(
        IDENTITIES
    )
    .reset_index()
)


dataset_dominance_wide = (
    identity_summary
    .pivot(
        index="reference_identity",
        columns="strategy",
        values="max_dataset_fraction",
    )
    .reindex(
        IDENTITIES
    )
    .reset_index()
)


# ============================================================
# 11. Worst donor-concentration cases
# ============================================================

worst_rows = []


for strategy in [
    x["strategy"]
    for x in STRATEGIES
]:

    temp = (
        dataset_identity_summary.loc[
            (
                dataset_identity_summary[
                    "strategy"
                ]
                == strategy
            )
            &
            (
                dataset_identity_summary[
                    "n_selected_simulated"
                ]
                > 0
            )
        ]
        .sort_values(
            "max_donor_fraction",
            ascending=False,
        )
        .head(15)
        .copy()
    )


    temp[
        "rank_within_strategy"
    ] = np.arange(
        1,
        len(temp) + 1,
    )


    worst_rows.append(
        temp
    )


worst_cases = pd.concat(
    worst_rows,
    ignore_index=True,
)


# ============================================================
# 12. Save outputs
# ============================================================

simulated_units.to_csv(
    OUT
    / "simulated_dataset_donor_identity_units.tsv",
    sep="\t",
    index=False,
)


identity_summary.to_csv(
    OUT
    / "strategy_identity_summary.tsv",
    sep="\t",
    index=False,
)


dataset_identity_summary.to_csv(
    OUT
    / "strategy_dataset_identity_summary.tsv",
    sep="\t",
    index=False,
)


dataset_summary.to_csv(
    OUT
    / "strategy_dataset_summary.tsv",
    sep="\t",
    index=False,
)


strategy_overview.to_csv(
    OUT
    / "strategy_overview.tsv",
    sep="\t",
    index=False,
)


selected_cells_wide.to_csv(
    OUT
    / "comparison_selected_cells_by_identity.tsv",
    sep="\t",
    index=False,
)


retention_wide.to_csv(
    OUT
    / "comparison_retention_by_identity.tsv",
    sep="\t",
    index=False,
)


donor_dominance_wide.to_csv(
    OUT
    / "comparison_donor_dominance_by_identity.tsv",
    sep="\t",
    index=False,
)


dataset_dominance_wide.to_csv(
    OUT
    / "comparison_dataset_dominance_by_identity.tsv",
    sep="\t",
    index=False,
)


worst_cases.to_csv(
    OUT
    / "worst_within_dataset_donor_dominance_cases.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 13. Console report
# ============================================================

print()
print("=" * 80)
print("STRATEGY OVERVIEW")
print("=" * 80)

print(
    strategy_overview.to_string(
        index=False,
    )
)


print()
print("=" * 80)
print("IDENTITY — SIMULATED SELECTED CELL COUNTS")
print("=" * 80)

print(
    selected_cells_wide.to_string(
        index=False,
    )
)


print()
print("=" * 80)
print("IDENTITY — MAX SINGLE DATASET×DONOR FRACTION")
print("=" * 80)

print(
    donor_dominance_wide.to_string(
        index=False,
    )
)


print()
print("=" * 80)
print("IDENTITY — MAX DATASET FRACTION")
print("=" * 80)

print(
    dataset_dominance_wide.to_string(
        index=False,
    )
)


print()
print("=" * 80)
print("WORST WITHIN-DATASET DONOR DOMINANCE")
print("=" * 80)

for strategy in [
    x["strategy"]
    for x in STRATEGIES
]:

    print()
    print(
        f"--- {strategy} ---"
    )

    temp = worst_cases.loc[
        worst_cases[
            "strategy"
        ]
        == strategy,
        [
            "rank_within_strategy",
            "dataset",
            "reference_identity",
            "n_selected_simulated",
            "n_donors",
            "max_donor_fraction",
            "effective_donors",
        ],
    ]

    print(
        temp.to_string(
            index=False,
        )
    )


print()
print("=" * 80)
print("OUTPUT DIRECTORY")
print("=" * 80)

print(OUT)


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)

print(
    "Sampling strategies simulated."
)

print(
    "No individual cell has been selected."
)

print(
    "Use the simulation results to choose "
    "the final donor-balancing rule."
)