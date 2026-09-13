#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
17_simulate_dataset_cap_sensitivity.py

Sensitivity analysis for dataset-level balancing after the primary
donor-balancing rule has been fixed.

PRIMARY RULE
------------
For every identity:

    donor_cap_i = min(round(Q75_i), 200)

where Q75_i is calculated from positive
dataset × donor × identity availability.

SECONDARY RULES TESTED
----------------------
BASE_Q75_MAX200
    No dataset-level correction.

DATASET_MAX60
    For non-tissue-restricted identities only:
    maximum single-dataset contribution <= 60%.

DATASET_MAX65
    Same, threshold = 65%.

DATASET_MAX70
    Same, threshold = 70%.

Tissue-restricted identities:
    Ileal_absorptive
    Paneth

are NEVER dataset-balanced.

IMPORTANT
---------
This script performs mathematical simulation only.

It does NOT:
    - select real cell IDs
    - modify ontology
    - modify prior audit outputs
    - read expression matrices
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
    / "04_dataset_cap_sensitivity"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Configuration
# ============================================================

MAX_DONOR_CAP = 200

TISSUE_RESTRICTED = {
    "Ileal_absorptive",
    "Paneth",
}


STRATEGIES = [
    {
        "strategy": "BASE_Q75_MAX200",
        "dataset_cap": None,
    },
    {
        "strategy": "DATASET_MAX60",
        "dataset_cap": 0.60,
    },
    {
        "strategy": "DATASET_MAX65",
        "dataset_cap": 0.65,
    },
    {
        "strategy": "DATASET_MAX70",
        "dataset_cap": 0.70,
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
            f"Required file missing:\n{path}"
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
        np.sum(
            p ** 2
        )
    )


def effective_n(values):

    value = hhi(
        values
    )

    if (
        not np.isfinite(value)
        or value <= 0
    ):
        return np.nan

    return float(
        1 / value
    )


# ============================================================
# 3. Integer allocation helper
# ============================================================

def allocate_integer_preserve_donors(
    counts,
    target,
    labels,
):
    """
    Reduce a vector of donor-level cell counts to an
    exact integer target.

    Principles
    ----------
    1. Never exceed original donor-level counts.
    2. If target >= number of positive donors,
       retain at least 1 cell from every donor.
    3. Remaining cells are allocated proportionally
       to remaining donor capacity.
    4. Hamilton / largest-remainder allocation makes
       the final sum exactly equal to target.

    This is simulation only; no actual cell IDs are selected.
    """

    counts = np.asarray(
        counts,
        dtype=int,
    )

    labels = np.asarray(
        labels,
        dtype=str,
    )

    if len(counts) != len(labels):
        fail(
            "Allocation labels/counts length mismatch."
        )


    total = int(
        counts.sum()
    )

    target = int(
        target
    )


    if target < 0:
        fail(
            "Negative allocation target."
        )


    if target >= total:
        return counts.copy()


    if target == 0:
        return np.zeros(
            len(counts),
            dtype=int,
        )


    positive = (
        counts > 0
    )

    n_positive = int(
        positive.sum()
    )


    allocation = np.zeros(
        len(counts),
        dtype=int,
    )


    # --------------------------------------------------------
    # Prefer preserving all positive donors.
    # --------------------------------------------------------

    if target >= n_positive:

        allocation[
            positive
        ] = 1

        remaining = (
            target
            - n_positive
        )

        capacity = (
            counts
            - allocation
        )


    else:

        # Very unlikely in this project, but defined.
        # If target is smaller than number of donors,
        # retain donors with highest availability.
        order = sorted(
            range(len(counts)),
            key=lambda i: (
                -counts[i],
                labels[i],
            ),
        )

        chosen = order[
            :target
        ]

        allocation[
            chosen
        ] = 1

        return allocation


    if remaining == 0:
        return allocation


    capacity_total = int(
        capacity.sum()
    )


    if remaining > capacity_total:
        fail(
            "Allocation target exceeds capacity."
        )


    if capacity_total == 0:
        return allocation


    raw_extra = (
        capacity
        / capacity_total
        * remaining
    )


    extra_floor = np.floor(
        raw_extra
    ).astype(int)


    # Never exceed donor capacity.
    extra_floor = np.minimum(
        extra_floor,
        capacity,
    )


    allocation += extra_floor


    left = (
        target
        - int(
            allocation.sum()
        )
    )


    if left > 0:

        fractions = (
            raw_extra
            - extra_floor
        )


        candidates = [
            i
            for i in range(
                len(counts)
            )
            if allocation[i]
            < counts[i]
        ]


        candidates = sorted(
            candidates,
            key=lambda i: (
                -fractions[i],
                labels[i],
            ),
        )


        for i in candidates:

            if left == 0:
                break

            if (
                allocation[i]
                < counts[i]
            ):
                allocation[i] += 1
                left -= 1


    # Rare fallback if numerical/tie behavior
    # leaves anything unresolved.
    if left > 0:

        candidates = sorted(
            [
                i
                for i in range(
                    len(counts)
                )
                if allocation[i]
                < counts[i]
            ],
            key=lambda i: (
                -(
                    counts[i]
                    - allocation[i]
                ),
                labels[i],
            ),
        )


        for i in candidates:

            while (
                left > 0
                and
                allocation[i]
                < counts[i]
            ):

                allocation[i] += 1
                left -= 1


            if left == 0:
                break


    if int(
        allocation.sum()
    ) != target:

        fail(
            "Integer allocation failed to "
            f"reach target: "
            f"{allocation.sum()} != {target}"
        )


    if np.any(
        allocation > counts
    ):
        fail(
            "Allocation exceeds donor availability."
        )


    return allocation


# ============================================================
# 4. Load inputs
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


required = {
    "dataset",
    "donor_id",
    "reference_identity",
    "n_cells",
}


missing = (
    required
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
        "rows detected."
    )


if (
    availability[
        "n_cells"
    ] <= 0
).any():

    fail(
        "Positive availability table contains "
        "non-positive cell counts."
    )


ontology = ontology.sort_values(
    "ontology_order"
)


IDENTITIES = (
    ontology[
        "reference_identity"
    ].tolist()
)


if len(IDENTITIES) != 20:

    fail(
        f"Expected 20 identities; "
        f"found {len(IDENTITIES)}."
    )


DATASETS = sorted(
    availability[
        "dataset"
    ].unique()
)


RAW_TOTAL = int(
    availability[
        "n_cells"
    ].sum()
)


# ============================================================
# 5. Reproduce primary Q75-MAX200 rule
# ============================================================

identity_q75 = (
    availability
    .groupby(
        "reference_identity",
        observed=True,
    )[
        "n_cells"
    ]
    .quantile(
        0.75
    )
)


identity_caps = {}


for identity in IDENTITIES:

    q75 = float(
        identity_q75.loc[
            identity
        ]
    )


    q75_rounded = max(
        1,
        int(
            np.rint(q75)
        ),
    )


    identity_caps[
        identity
    ] = min(
        q75_rounded,
        MAX_DONOR_CAP,
    )


primary = availability.copy()


primary[
    "identity_donor_cap"
] = (
    primary[
        "reference_identity"
    ]
    .map(
        identity_caps
    )
)


primary[
    "n_after_donor_cap"
] = np.minimum(
    primary[
        "n_cells"
    ].to_numpy(),
    primary[
        "identity_donor_cap"
    ].to_numpy(),
).astype(int)


primary_total = int(
    primary[
        "n_after_donor_cap"
    ].sum()
)


print("=" * 82)
print("DATASET-CAP SENSITIVITY SIMULATION")
print("=" * 82)

print(
    f"Raw reference universe: "
    f"{RAW_TOTAL:,}"
)

print(
    f"After primary Q75_MAX200 donor cap: "
    f"{primary_total:,}"
)


# ============================================================
# 6. Simulate dataset-level corrections
# ============================================================

all_units = []
correction_rows = []


for strategy_cfg in STRATEGIES:

    strategy = (
        strategy_cfg[
            "strategy"
        ]
    )

    threshold = (
        strategy_cfg[
            "dataset_cap"
        ]
    )


    sim = primary.copy()


    sim[
        "n_selected_simulated"
    ] = (
        sim[
            "n_after_donor_cap"
        ]
        .astype(int)
    )


    # --------------------------------------------------------
    # Process identity by identity.
    # --------------------------------------------------------

    for identity in IDENTITIES:

        identity_mask = (
            sim[
                "reference_identity"
            ]
            == identity
        )


        subset = sim.loc[
            identity_mask
        ]


        dataset_totals = (
            subset
            .groupby(
                "dataset",
                observed=True,
            )[
                "n_selected_simulated"
            ]
            .sum()
        )


        identity_total_before = int(
            dataset_totals.sum()
        )


        dominant_dataset = (
            dataset_totals.idxmax()
        )


        dominant_before = int(
            dataset_totals.loc[
                dominant_dataset
            ]
        )


        dominant_fraction_before = (
            dominant_before
            / identity_total_before
        )


        tissue_restricted = (
            identity
            in TISSUE_RESTRICTED
        )


        triggered = False
        target_dominant = (
            dominant_before
        )


        # ----------------------------------------------------
        # Secondary correction.
        # ----------------------------------------------------

        if (
            threshold is not None
            and
            not tissue_restricted
            and
            dominant_fraction_before
            > threshold
        ):

            other_total = (
                identity_total_before
                - dominant_before
            )


            if other_total <= 0:

                fail(
                    f"{strategy}/{identity}: "
                    "cannot dataset-balance an "
                    "identity present in only "
                    "one dataset."
                )


            # Solve:
            #
            # D' / (D' + O) <= threshold
            #
            # D' <= threshold * O / (1-threshold)

            target_dominant = int(
                np.floor(
                    threshold
                    * other_total
                    /
                    (
                        1
                        - threshold
                    )
                )
            )


            target_dominant = min(
                target_dominant,
                dominant_before,
            )


            target_dominant = max(
                target_dominant,
                1,
            )


            dominant_mask = (
                identity_mask
                &
                (
                    sim["dataset"]
                    == dominant_dataset
                )
            )


            donor_rows = (
                sim.loc[
                    dominant_mask
                ]
                .copy()
            )


            allocations = (
                allocate_integer_preserve_donors(
                    donor_rows[
                        "n_selected_simulated"
                    ].to_numpy(),
                    target_dominant,
                    donor_rows[
                        "donor_id"
                    ].astype(str).to_numpy(),
                )
            )


            sim.loc[
                donor_rows.index,
                "n_selected_simulated",
            ] = allocations


            triggered = True


        # ----------------------------------------------------
        # Recalculate post-correction composition.
        # ----------------------------------------------------

        post_subset = sim.loc[
            identity_mask
        ]


        post_dataset_totals = (
            post_subset
            .groupby(
                "dataset",
                observed=True,
            )[
                "n_selected_simulated"
            ]
            .sum()
        )


        identity_total_after = int(
            post_dataset_totals.sum()
        )


        post_dominant_dataset = (
            post_dataset_totals.idxmax()
        )


        post_dominant_n = int(
            post_dataset_totals.loc[
                post_dominant_dataset
            ]
        )


        post_fraction = (
            post_dominant_n
            / identity_total_after
        )


        if (
            threshold is not None
            and
            triggered
            and
            post_fraction
            > threshold + 1e-6
        ):

            fail(
                f"{strategy}/{identity}: "
                "dataset correction failed. "
                f"Post fraction={post_fraction:.6f}, "
                f"threshold={threshold:.6f}"
            )


        correction_rows.append(
            {
                "strategy":
                    strategy,

                "dataset_fraction_threshold":
                    threshold,

                "reference_identity":
                    identity,

                "tissue_restricted":
                    tissue_restricted,

                "correction_triggered":
                    triggered,

                "dominant_dataset_before":
                    dominant_dataset,

                "dominant_dataset_cells_before":
                    dominant_before,

                "max_dataset_fraction_before":
                    dominant_fraction_before,

                "dominant_dataset_target":
                    target_dominant,

                "identity_cells_before_dataset_correction":
                    identity_total_before,

                "identity_cells_after_dataset_correction":
                    identity_total_after,

                "cells_removed_by_dataset_correction":
                    (
                        identity_total_before
                        - identity_total_after
                    ),

                "dominant_dataset_after":
                    post_dominant_dataset,

                "max_dataset_fraction_after":
                    post_fraction,
            }
        )


    sim[
        "strategy"
    ] = strategy


    sim[
        "dataset_cap_threshold"
    ] = threshold


    all_units.append(
        sim
    )


simulated_units = pd.concat(
    all_units,
    ignore_index=True,
)


corrections = pd.DataFrame(
    correction_rows
)


# ============================================================
# 7. Identity summaries
# ============================================================

identity_rows = []


for strategy_cfg in STRATEGIES:

    strategy = (
        strategy_cfg[
            "strategy"
        ]
    )


    s = simulated_units.loc[
        simulated_units[
            "strategy"
        ]
        == strategy
    ]


    for identity in IDENTITIES:

        x = s.loc[
            s[
                "reference_identity"
            ]
            == identity
        ]


        selected = (
            x[
                "n_selected_simulated"
            ].to_numpy()
        )


        selected_total = int(
            selected.sum()
        )


        raw_total = int(
            x[
                "n_cells"
            ].sum()
        )


        primary_total_identity = int(
            x[
                "n_after_donor_cap"
            ].sum()
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


        identity_rows.append(
            {
                "strategy":
                    strategy,

                "reference_identity":
                    identity,

                "tissue_restricted":
                    identity
                    in TISSUE_RESTRICTED,

                "n_raw_cells":
                    raw_total,

                "n_after_primary_donor_cap":
                    primary_total_identity,

                "n_selected_simulated":
                    selected_total,

                "fraction_of_raw_retained":
                    (
                        selected_total
                        / raw_total
                    ),

                "fraction_of_primary_retained":
                    (
                        selected_total
                        / primary_total_identity
                    ),

                "n_positive_dataset_donor_units":
                    int(
                        (
                            selected > 0
                        ).sum()
                    ),

                "max_dataset_donor_fraction":
                    float(
                        selected.max()
                        / selected_total
                    ),

                "dataset_donor_HHI":
                    hhi(
                        selected
                    ),

                "effective_dataset_donor_units":
                    effective_n(
                        selected
                    ),

                "dominant_dataset":
                    dominant_dataset,

                "max_dataset_fraction":
                    float(
                        dataset_fractions.max()
                    ),
            }
        )


identity_summary = pd.DataFrame(
    identity_rows
)


# ============================================================
# 8. Dataset × identity donor-balance summary
# ============================================================

dataset_identity_rows = []


for strategy_cfg in STRATEGIES:

    strategy = (
        strategy_cfg[
            "strategy"
        ]
    )


    s = simulated_units.loc[
        simulated_units[
            "strategy"
        ]
        == strategy
    ]


    for identity in IDENTITIES:

        for dataset in DATASETS:

            x = s.loc[
                (
                    s[
                        "reference_identity"
                    ]
                    == identity
                )
                &
                (
                    s["dataset"]
                    == dataset
                )
            ]


            if len(x) == 0:

                dataset_identity_rows.append(
                    {
                        "strategy":
                            strategy,

                        "reference_identity":
                            identity,

                        "dataset":
                            dataset,

                        "n_selected_simulated":
                            0,

                        "n_positive_donors":
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


            vals = (
                x[
                    "n_selected_simulated"
                ].to_numpy()
            )


            vals_positive = (
                vals[
                    vals > 0
                ]
            )


            total = int(
                vals_positive.sum()
            )


            dataset_identity_rows.append(
                {
                    "strategy":
                        strategy,

                    "reference_identity":
                        identity,

                    "dataset":
                        dataset,

                    "n_selected_simulated":
                        total,

                    "n_positive_donors":
                        len(
                            vals_positive
                        ),

                    "max_donor_fraction":
                        (
                            float(
                                vals_positive.max()
                                / total
                            )
                            if total > 0
                            else np.nan
                        ),

                    "donor_HHI":
                        hhi(
                            vals_positive
                        ),

                    "effective_donors":
                        effective_n(
                            vals_positive
                        ),
                }
            )


dataset_identity_summary = (
    pd.DataFrame(
        dataset_identity_rows
    )
)


# ============================================================
# 9. Strategy overview
# ============================================================

overview_rows = []


for strategy_cfg in STRATEGIES:

    strategy = (
        strategy_cfg[
            "strategy"
        ]
    )

    threshold = (
        strategy_cfg[
            "dataset_cap"
        ]
    )


    x = identity_summary.loc[
        identity_summary[
            "strategy"
        ]
        == strategy
    ]


    corr = corrections.loc[
        corrections[
            "strategy"
        ]
        == strategy
    ]


    non_tissue = x.loc[
        ~x[
            "tissue_restricted"
        ]
    ]


    worst_non_tissue_idx = (
        non_tissue[
            "max_dataset_fraction"
        ].idxmax()
    )


    worst_non_tissue_identity = (
        non_tissue.loc[
            worst_non_tissue_idx,
            "reference_identity",
        ]
    )


    worst_non_tissue_fraction = float(
        non_tissue.loc[
            worst_non_tissue_idx,
            "max_dataset_fraction",
        ]
    )


    total_selected = int(
        x[
            "n_selected_simulated"
        ].sum()
    )


    overview_rows.append(
        {
            "strategy":
                strategy,

            "dataset_cap_threshold":
                threshold,

            "total_selected_cells":
                total_selected,

            "fraction_raw_reference_retained":
                (
                    total_selected
                    / RAW_TOTAL
                ),

            "fraction_primary_Q75MAX200_retained":
                (
                    total_selected
                    / primary_total
                ),

            "n_non_tissue_identities_corrected":
                int(
                    (
                        corr[
                            "correction_triggered"
                        ]
                    ).sum()
                ),

            "cells_removed_by_dataset_correction":
                int(
                    corr[
                        "cells_removed_by_dataset_correction"
                    ].sum()
                ),

            "median_max_dataset_fraction":
                float(
                    x[
                        "max_dataset_fraction"
                    ].median()
                ),

            "worst_non_tissue_identity":
                worst_non_tissue_identity,

            "worst_non_tissue_max_dataset_fraction":
                worst_non_tissue_fraction,

            "median_max_dataset_donor_fraction":
                float(
                    x[
                        "max_dataset_donor_fraction"
                    ].median()
                ),

            "worst_max_dataset_donor_fraction":
                float(
                    x[
                        "max_dataset_donor_fraction"
                    ].max()
                ),

            "median_effective_dataset_donor_units":
                float(
                    x[
                        "effective_dataset_donor_units"
                    ].median()
                ),

            "minimum_identity_selected_cells":
                int(
                    x[
                        "n_selected_simulated"
                    ].min()
                ),
        }
    )


strategy_overview = pd.DataFrame(
    overview_rows
)


# ============================================================
# 10. Wide comparison tables
# ============================================================

selected_wide = (
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


dataset_fraction_wide = (
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


donor_fraction_wide = (
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


# ============================================================
# 11. Save outputs
# ============================================================

simulated_units.to_csv(
    OUT
    / "simulated_dataset_donor_identity_units.tsv",
    sep="\t",
    index=False,
)


corrections.to_csv(
    OUT
    / "dataset_correction_details.tsv",
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


strategy_overview.to_csv(
    OUT
    / "strategy_overview.tsv",
    sep="\t",
    index=False,
)


selected_wide.to_csv(
    OUT
    / "comparison_selected_cells.tsv",
    sep="\t",
    index=False,
)


dataset_fraction_wide.to_csv(
    OUT
    / "comparison_max_dataset_fraction.tsv",
    sep="\t",
    index=False,
)


donor_fraction_wide.to_csv(
    OUT
    / "comparison_max_dataset_donor_fraction.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 12. Console report
# ============================================================

print()
print("=" * 82)
print("STRATEGY OVERVIEW")
print("=" * 82)

print(
    strategy_overview.to_string(
        index=False,
    )
)


print()
print("=" * 82)
print("TRIGGERED DATASET CORRECTIONS")
print("=" * 82)


triggered = corrections.loc[
    corrections[
        "correction_triggered"
    ]
].copy()


if len(triggered) == 0:

    print(
        "No corrections triggered."
    )

else:

    print(
        triggered[
            [
                "strategy",
                "reference_identity",
                "dominant_dataset_before",
                "max_dataset_fraction_before",
                "dominant_dataset_cells_before",
                "dominant_dataset_target",
                "cells_removed_by_dataset_correction",
                "max_dataset_fraction_after",
            ]
        ].to_string(
            index=False,
        )
    )


print()
print("=" * 82)
print("IDENTITY — SELECTED CELL COUNTS")
print("=" * 82)

print(
    selected_wide.to_string(
        index=False,
    )
)


print()
print("=" * 82)
print("IDENTITY — MAX DATASET FRACTION")
print("=" * 82)

print(
    dataset_fraction_wide.to_string(
        index=False,
    )
)


print()
print("=" * 82)
print("IDENTITY — MAX DATASET×DONOR FRACTION")
print("=" * 82)

print(
    donor_fraction_wide.to_string(
        index=False,
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
    "Dataset-cap sensitivity analysis completed."
)

print(
    "No real cell IDs were selected."
)

print(
    "Primary donor rule remained Q75_MAX200 "
    "for all strategies."
)