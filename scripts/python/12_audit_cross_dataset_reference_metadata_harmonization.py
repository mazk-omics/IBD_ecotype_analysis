#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_audit_cross_dataset_reference_harmonization.py

Cross-dataset harmonization audit for the IBD EcoTyper
single-cell reference.

Datasets:
    - GSE282122
    - SCP1884
    - SCP259

This script DOES NOT:
    - re-read raw expression matrices
    - repeat dataset-level audits
    - modify old audit outputs
    - promote Glial/Tuft into the final ontology automatically

This script DOES:
    - reuse established audit outputs
    - harmonize identity/action tables
    - validate the frozen 18-class core taxonomy
    - evaluate cross-dataset representation
    - summarize Glial/Tuft candidates
    - inventory all other HOLD identities
    - produce a formal hand-off table for final ontology decisions
"""

from pathlib import Path
import sys

import pandas as pd


# ============================================================
# 0. Configuration
# ============================================================

PROJECT_ROOT = Path("/home/mazekai/IBD_EcoTyper")

OUTPUT_DIR = (
    PROJECT_ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "00_cross_dataset_harmonization"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ------------------------------------------------------------
# Established dataset-level audit outputs
# ------------------------------------------------------------

INPUT_FILES = {
    "GSE282122": (
        PROJECT_ROOT
        / "03_reference"
        / "GSE282122"
        / "audit_v2"
        / "output"
        / "04_identity_coverage_v1"
        / "all_candidate_cells_by_identity_action.tsv"
    ),
    "SCP1884": (
        PROJECT_ROOT
        / "03_reference"
        / "SCP1884"
        / "audit_v1"
        / "output"
        / "02_candidate_taxonomy"
        / "provisional_identity_counts.tsv"
    ),
    "SCP259": (
        PROJECT_ROOT
        / "03_reference"
        / "SCP259"
        / "audit_v1"
        / "output"
        / "01_candidate_taxonomy_coverage"
        / "identity_coverage.tsv"
    ),
}


# ------------------------------------------------------------
# Frozen 18-class core ontology
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Candidate additional identities
# ------------------------------------------------------------

CANDIDATE_IDENTITIES = [
    "Glial",
    "Tuft",
]


# ------------------------------------------------------------
# Identities for which absence from a colon-only dataset
# must not be interpreted as taxonomy failure.
# ------------------------------------------------------------

TISSUE_RESTRICTED_IDENTITIES = {
    "Ileal_absorptive",
    "Paneth",
}


# ------------------------------------------------------------
# Allowed action vocabulary inherited from old audits
# ------------------------------------------------------------

ALLOWED_ACTIONS = {
    "KEEP",
    "HOLD_CROSS_DATASET",
    "EXCLUDE_PRIMARY",
}


# ------------------------------------------------------------
# Expected core representation from established audits
# ------------------------------------------------------------

EXPECTED_CORE_BY_DATASET = {
    "GSE282122": set(CORE_IDENTITIES),
    "SCP1884": set(CORE_IDENTITIES),

    # SCP259 is a colon-focused dataset.
    # Ileal absorptive and Paneth are not expected.
    "SCP259": (
        set(CORE_IDENTITIES)
        - {
            "Ileal_absorptive",
            "Paneth",
        }
    ),
}


# ------------------------------------------------------------
# Historical hard checkpoints for GSE282122
# ------------------------------------------------------------

GSE282122_EXPECTED = {
    "KEEP": 484_979,
    "HOLD_CROSS_DATASET": 16_238,
    "EXCLUDE_PRIMARY": 4_374,
    "TOTAL": 505_591,
}


# ============================================================
# 1. Utilities
# ============================================================

def fail(message: str) -> None:
    """Stop with a clear error message."""
    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )
    sys.exit(1)


def detect_column(
    df: pd.DataFrame,
    candidates,
    role: str,
    dataset: str,
):
    """
    Detect one structural column conservatively.

    Stops rather than guessing when zero or multiple
    candidate columns are present.
    """

    hits = [
        x
        for x in candidates
        if x in df.columns
    ]

    if len(hits) == 0:
        fail(
            f"{dataset}: could not identify "
            f"{role} column.\n"
            f"Tried: {candidates}\n"
            f"Observed: {list(df.columns)}"
        )

    if len(hits) > 1:
        fail(
            f"{dataset}: ambiguous {role} columns: "
            f"{hits}"
        )

    return hits[0]


def read_action_table(
    dataset: str,
    path: Path,
) -> pd.DataFrame:
    """
    Read and standardize an established action-level
    audit table.
    """

    if not path.exists():
        fail(
            f"{dataset}: input file not found:\n"
            f"{path}"
        )

    df = pd.read_csv(
        path,
        sep="\t",
        dtype_backend="numpy_nullable",
    )

    print(
        f"{dataset}: "
        f"{len(df):,} rows read"
    )

    identity_col = detect_column(
        df,
        [
            "coarse_identity",
            "provisional_coarse_identity",
            "identity",
            "reference_identity",
        ],
        "identity",
        dataset,
    )

    action_col = detect_column(
        df,
        [
            "primary_reference_action",
            "reference_action",
            "action",
        ],
        "action",
        dataset,
    )

    count_col = detect_column(
        df,
        [
            "n_cells",
            "cell_count",
            "count",
        ],
        "cell-count",
        dataset,
    )

    sample_col = detect_column(
        df,
        [
            "n_samples",
            "sample_count",
        ],
        "sample-count",
        dataset,
    )

    donor_col = detect_column(
        df,
        [
            "n_donors",
            "donor_count",
        ],
        "donor-count",
        dataset,
    )

    out = pd.DataFrame(
        {
            "dataset": dataset,
            "identity": (
                df[identity_col]
                .astype("string")
                .str.strip()
            ),
            "action": (
                df[action_col]
                .astype("string")
                .str.strip()
            ),
            "n_cells": pd.to_numeric(
                df[count_col],
                errors="raise",
            ),
            "n_samples": pd.to_numeric(
                df[sample_col],
                errors="raise",
            ),
            "n_donors": pd.to_numeric(
                df[donor_col],
                errors="raise",
            ),
        }
    )

    if out[
        [
            "identity",
            "action",
            "n_cells",
            "n_samples",
            "n_donors",
        ]
    ].isna().any().any():
        fail(
            f"{dataset}: missing structural "
            f"values after standardization."
        )

    unexpected_actions = (
        set(out["action"])
        - ALLOWED_ACTIONS
    )

    if unexpected_actions:
        fail(
            f"{dataset}: unexpected action values: "
            f"{sorted(unexpected_actions)}"
        )

    if (
        (out["n_cells"] < 0).any()
        or (out["n_samples"] < 0).any()
        or (out["n_donors"] < 0).any()
    ):
        fail(
            f"{dataset}: negative counts detected."
        )

    duplicated = out.duplicated(
        subset=[
            "dataset",
            "identity",
            "action",
        ],
        keep=False,
    )

    if duplicated.any():
        print(
            out.loc[
                duplicated
            ].sort_values(
                [
                    "identity",
                    "action",
                ]
            )
        )

        fail(
            f"{dataset}: duplicate "
            "identity × action rows detected."
        )

    return out


# ============================================================
# 2. Start
# ============================================================

print("=" * 72)
print("IBD REFERENCE CROSS-DATASET HARMONIZATION")
print("=" * 72)
print()


# ============================================================
# 3. Read established audit tables
# ============================================================

tables = []

for dataset, path in INPUT_FILES.items():

    print(
        f"Reading {dataset}..."
    )

    tab = read_action_table(
        dataset,
        path,
    )

    tables.append(tab)

inventory = pd.concat(
    tables,
    ignore_index=True,
)

print()
print(
    "PASS: all established audit tables loaded"
)


# ============================================================
# 4. Save normalized action inventory
# ============================================================

inventory = inventory.sort_values(
    [
        "dataset",
        "action",
        "identity",
    ]
).reset_index(drop=True)

inventory.to_csv(
    OUTPUT_DIR
    / "normalized_identity_action_inventory.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 5. Dataset-level action totals
# ============================================================

action_totals = (
    inventory
    .groupby(
        [
            "dataset",
            "action",
        ],
        as_index=False,
        dropna=False,
    )
    .agg(
        n_cells=(
            "n_cells",
            "sum",
        ),
    )
)

action_totals.to_csv(
    OUTPUT_DIR
    / "dataset_action_totals.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 72)
print("DATASET ACTION TOTALS")
print("=" * 72)

print(
    action_totals.to_string(
        index=False
    )
)


# ============================================================
# 6. Reproduce historical GSE282122 checkpoints
# ============================================================

gse = inventory.loc[
    inventory["dataset"]
    == "GSE282122"
].copy()

for action in [
    "KEEP",
    "HOLD_CROSS_DATASET",
    "EXCLUDE_PRIMARY",
]:

    observed = int(
        gse.loc[
            gse["action"] == action,
            "n_cells",
        ].sum()
    )

    expected = (
        GSE282122_EXPECTED[action]
    )

    if observed != expected:
        fail(
            "GSE282122 historical checkpoint "
            f"failed for {action}: "
            f"expected {expected:,}, "
            f"observed {observed:,}."
        )

gse_total = int(
    gse["n_cells"].sum()
)

if (
    gse_total
    != GSE282122_EXPECTED["TOTAL"]
):
    fail(
        "GSE282122 total candidate-cell "
        "checkpoint failed: "
        f"expected "
        f"{GSE282122_EXPECTED['TOTAL']:,}, "
        f"observed {gse_total:,}."
    )

print()
print(
    "PASS: GSE282122 historical "
    "505,591-cell audit reproduced"
)


# ============================================================
# 7. Validate KEEP ontology per dataset
# ============================================================

print()
print("=" * 72)
print("CORE KEEP ONTOLOGY VALIDATION")
print("=" * 72)

for dataset in INPUT_FILES:

    keep_set = set(
        inventory.loc[
            (
                inventory["dataset"]
                == dataset
            )
            & (
                inventory["action"]
                == "KEEP"
            ),
            "identity",
        ]
    )

    expected = (
        EXPECTED_CORE_BY_DATASET[
            dataset
        ]
    )

    unexpected = (
        keep_set
        - set(CORE_IDENTITIES)
    )

    if unexpected:
        fail(
            f"{dataset}: unexpected identities "
            f"in KEEP pool: "
            f"{sorted(unexpected)}"
        )

    if keep_set != expected:

        missing = sorted(
            expected
            - keep_set
        )

        extra = sorted(
            keep_set
            - expected
        )

        fail(
            f"{dataset}: KEEP ontology differs "
            "from established expectation.\n"
            f"Missing: {missing}\n"
            f"Extra: {extra}"
        )

    print(
        f"{dataset}: PASS "
        f"({len(keep_set)} core identities)"
    )


# ============================================================
# 8. Build complete dataset × core-identity grid
# ============================================================

datasets = list(
    INPUT_FILES.keys()
)

core_grid = pd.MultiIndex.from_product(
    [
        datasets,
        CORE_IDENTITIES,
    ],
    names=[
        "dataset",
        "identity",
    ],
).to_frame(index=False)

keep_inventory = (
    inventory.loc[
        inventory["action"]
        == "KEEP"
    ]
    [
        [
            "dataset",
            "identity",
            "action",
            "n_cells",
            "n_samples",
            "n_donors",
        ]
    ]
    .copy()
)

core_long = core_grid.merge(
    keep_inventory,
    on=[
        "dataset",
        "identity",
    ],
    how="left",
    validate="one_to_one",
)

core_long["present"] = (
    core_long["n_cells"]
    .fillna(0)
    > 0
)

for col in [
    "n_cells",
    "n_samples",
    "n_donors",
]:
    core_long[col] = (
        core_long[col]
        .fillna(0)
        .astype(int)
    )

core_long["action"] = (
    core_long["action"]
    .fillna("ABSENT")
)


# ============================================================
# 9. Core identity cross-dataset summary
# ============================================================

core_summary_rows = []

for identity in CORE_IDENTITIES:

    sub = core_long.loc[
        core_long["identity"]
        == identity
    ].copy()

    present_datasets = (
        sub.loc[
            sub["present"],
            "dataset",
        ]
        .tolist()
    )

    n_present = len(
        present_datasets
    )

    tissue_restricted = (
        identity
        in TISSUE_RESTRICTED_IDENTITIES
    )

    if n_present == 3:

        decision = (
            "KEEP_CORE_CROSS_DATASET"
        )

    elif (
        tissue_restricted
        and n_present >= 2
    ):

        decision = (
            "KEEP_CORE_"
            "TISSUE_RESTRICTED"
        )

    else:

        decision = (
            "REVIEW_CORE_COVERAGE"
        )

    row = {
        "identity": identity,
        "core_taxonomy_member": True,
        "tissue_restricted":
            tissue_restricted,
        "n_datasets_present":
            n_present,
        "datasets_present":
            " | ".join(
                present_datasets
            ),
        "total_cells_across_datasets":
            int(
                sub["n_cells"].sum()
            ),
        "sum_dataset_donors":
            int(
                sub["n_donors"].sum()
            ),
        "harmonization_decision":
            decision,
    }

    for dataset in datasets:

        x = sub.loc[
            sub["dataset"]
            == dataset
        ].iloc[0]

        row[
            f"{dataset}__present"
        ] = bool(
            x["present"]
        )

        row[
            f"{dataset}__n_cells"
        ] = int(
            x["n_cells"]
        )

        row[
            f"{dataset}__n_samples"
        ] = int(
            x["n_samples"]
        )

        row[
            f"{dataset}__n_donors"
        ] = int(
            x["n_donors"]
        )

    core_summary_rows.append(
        row
    )

core_summary = pd.DataFrame(
    core_summary_rows
)

core_summary.to_csv(
    OUTPUT_DIR
    / "core_identity_cross_dataset_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 10. Hard validation of core decisions
# ============================================================

bad_core = core_summary.loc[
    core_summary[
        "harmonization_decision"
    ]
    == "REVIEW_CORE_COVERAGE"
]

if len(bad_core) > 0:

    print()
    print(
        "Core identities requiring review:"
    )

    print(
        bad_core[
            [
                "identity",
                "n_datasets_present",
                "datasets_present",
            ]
        ].to_string(
            index=False
        )
    )

    fail(
        "At least one frozen core identity "
        "has unexpected cross-dataset coverage."
    )

print()
print(
    "PASS: frozen 18-class core ontology "
    "is compatible across the three audits"
)


# ============================================================
# 11. Candidate Glial / Tuft assessment
# ============================================================

candidate_rows = []

for identity in CANDIDATE_IDENTITIES:

    sub = inventory.loc[
        inventory["identity"]
        == identity
    ].copy()

    rows = []

    for dataset in datasets:

        x = sub.loc[
            sub["dataset"]
            == dataset
        ]

        if len(x) == 0:

            rows.append(
                {
                    "dataset":
                        dataset,
                    "action":
                        "ABSENT",
                    "n_cells":
                        0,
                    "n_samples":
                        0,
                    "n_donors":
                        0,
                }
            )

        elif len(x) == 1:

            rows.append(
                {
                    "dataset":
                        dataset,
                    "action":
                        x.iloc[0][
                            "action"
                        ],
                    "n_cells":
                        int(
                            x.iloc[0][
                                "n_cells"
                            ]
                        ),
                    "n_samples":
                        int(
                            x.iloc[0][
                                "n_samples"
                            ]
                        ),
                    "n_donors":
                        int(
                            x.iloc[0][
                                "n_donors"
                            ]
                        ),
                }
            )

        else:

            fail(
                f"{dataset}: multiple action "
                f"rows detected for candidate "
                f"{identity}."
            )

    temp = pd.DataFrame(
        rows
    )

    n_present = int(
        (temp["n_cells"] > 0).sum()
    )

    all_hold = bool(
        (
            temp["action"]
            == "HOLD_CROSS_DATASET"
        ).all()
    )

    if (
        n_present == 3
        and all_hold
    ):
        decision = (
            "READY_FOR_"
            "EXPRESSION_VALIDATION"
        )
    else:
        decision = (
            "REVIEW_CANDIDATE_"
            "COVERAGE"
        )

    row = {
        "identity": identity,
        "n_datasets_present":
            n_present,
        "total_cells_across_datasets":
            int(
                temp["n_cells"].sum()
            ),
        "sum_dataset_donors":
            int(
                temp["n_donors"].sum()
            ),
        "all_datasets_hold_status":
            all_hold,
        "harmonization_decision":
            decision,
    }

    for _, x in temp.iterrows():

        dataset = x["dataset"]

        row[
            f"{dataset}__action"
        ] = x["action"]

        row[
            f"{dataset}__n_cells"
        ] = int(
            x["n_cells"]
        )

        row[
            f"{dataset}__n_samples"
        ] = int(
            x["n_samples"]
        )

        row[
            f"{dataset}__n_donors"
        ] = int(
            x["n_donors"]
        )

    candidate_rows.append(
        row
    )

candidate_summary = pd.DataFrame(
    candidate_rows
)

candidate_summary.to_csv(
    OUTPUT_DIR
    / "candidate_identity_cross_dataset_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 12. Validate candidate representation
# ============================================================

bad_candidate = (
    candidate_summary.loc[
        candidate_summary[
            "harmonization_decision"
        ]
        != (
            "READY_FOR_"
            "EXPRESSION_VALIDATION"
        )
    ]
)

if len(bad_candidate) > 0:

    print()
    print(
        "WARNING: candidate identity "
        "coverage requires review:"
    )

    print(
        bad_candidate.to_string(
            index=False
        )
    )

else:

    print()
    print(
        "PASS: Glial and Tuft are "
        "represented as HOLD identities "
        "in all three datasets"
    )


# ============================================================
# 13. Inventory all HOLD identities
# ============================================================

hold_inventory = (
    inventory.loc[
        inventory["action"]
        == "HOLD_CROSS_DATASET"
    ]
    .copy()
)

hold_inventory[
    "candidate_for_final_test"
] = (
    hold_inventory["identity"]
    .isin(
        CANDIDATE_IDENTITIES
    )
)

hold_inventory = (
    hold_inventory
    .sort_values(
        [
            "identity",
            "dataset",
        ]
    )
    .reset_index(drop=True)
)

hold_inventory.to_csv(
    OUTPUT_DIR
    / "hold_identity_inventory.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 14. HOLD presence matrix
# ============================================================

all_hold_identities = sorted(
    hold_inventory[
        "identity"
    ].unique()
)

hold_matrix_rows = []

for identity in all_hold_identities:

    row = {
        "identity": identity,
        "candidate_for_final_test":
            identity
            in CANDIDATE_IDENTITIES,
    }

    n_present = 0
    total_cells = 0

    for dataset in datasets:

        x = hold_inventory.loc[
            (
                hold_inventory[
                    "identity"
                ]
                == identity
            )
            & (
                hold_inventory[
                    "dataset"
                ]
                == dataset
            )
        ]

        if len(x) == 1:

            n = int(
                x.iloc[0][
                    "n_cells"
                ]
            )

            row[
                f"{dataset}__n_cells"
            ] = n

            row[
                f"{dataset}__n_donors"
            ] = int(
                x.iloc[0][
                    "n_donors"
                ]
            )

            row[
                f"{dataset}__present"
            ] = True

            n_present += 1
            total_cells += n

        elif len(x) == 0:

            row[
                f"{dataset}__n_cells"
            ] = 0

            row[
                f"{dataset}__n_donors"
            ] = 0

            row[
                f"{dataset}__present"
            ] = False

        else:

            fail(
                "Unexpected duplicate HOLD "
                f"rows for {dataset} / "
                f"{identity}"
            )

    row[
        "n_datasets_present"
    ] = n_present

    row[
        "total_cells_across_datasets"
    ] = total_cells

    hold_matrix_rows.append(
        row
    )

hold_matrix = pd.DataFrame(
    hold_matrix_rows
)

hold_matrix.to_csv(
    OUTPUT_DIR
    / "hold_identity_cross_dataset_matrix.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 15. Build taxonomy review table
# ============================================================

review_rows = []

for _, row in core_summary.iterrows():

    review_rows.append(
        {
            "identity":
                row["identity"],
            "taxonomy_scope":
                "CORE",
            "current_status":
                "KEEP",
            "n_datasets_present":
                row[
                    "n_datasets_present"
                ],
            "total_cells_across_datasets":
                row[
                    "total_cells_across_datasets"
                ],
            "next_step":
                "PROCEED_TO_BALANCED_SAMPLING",
            "decision_note":
                row[
                    "harmonization_decision"
                ],
        }
    )

for _, row in candidate_summary.iterrows():

    review_rows.append(
        {
            "identity":
                row["identity"],
            "taxonomy_scope":
                "CANDIDATE",
            "current_status":
                "HOLD_CROSS_DATASET",
            "n_datasets_present":
                row[
                    "n_datasets_present"
                ],
            "total_cells_across_datasets":
                row[
                    "total_cells_across_datasets"
                ],
            "next_step":
                "EXPRESSION_LEVEL_VALIDATION",
            "decision_note":
                row[
                    "harmonization_decision"
                ],
        }
    )


other_holds = (
    hold_matrix.loc[
        ~hold_matrix[
            "identity"
        ].isin(
            CANDIDATE_IDENTITIES
        )
    ]
)

for _, row in other_holds.iterrows():

    review_rows.append(
        {
            "identity":
                row["identity"],
            "taxonomy_scope":
                "OTHER_HOLD",
            "current_status":
                "HOLD_CROSS_DATASET",
            "n_datasets_present":
                row[
                    "n_datasets_present"
                ],
            "total_cells_across_datasets":
                row[
                    "total_cells_across_datasets"
                ],
            "next_step":
                "DO_NOT_INCLUDE_IN_PRIMARY_REFERENCE",
            "decision_note":
                "Retain for provenance; "
                "not part of current "
                "18-class core or "
                "Glial/Tuft candidate test.",
        }
    )


taxonomy_review = pd.DataFrame(
    review_rows
).sort_values(
    [
        "taxonomy_scope",
        "identity",
    ]
)

taxonomy_review.to_csv(
    OUTPUT_DIR
    / "taxonomy_decision_review_table.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 16. Human-readable summary
# ============================================================

summary_file = (
    OUTPUT_DIR
    / "cross_dataset_harmonization_summary.txt"
)

with summary_file.open(
    "w",
    encoding="utf-8",
) as f:

    f.write(
        "IBD EcoTyper reference\n"
    )

    f.write(
        "Cross-dataset taxonomy "
        "harmonization audit\n"
    )

    f.write(
        "=" * 72
        + "\n\n"
    )

    f.write(
        "Datasets\n"
    )

    f.write(
        "-" * 72
        + "\n"
    )

    for dataset in datasets:
        f.write(
            f"{dataset}: "
            f"{INPUT_FILES[dataset]}\n"
        )

    f.write("\n")

    f.write(
        "Core ontology\n"
    )

    f.write(
        "-" * 72
        + "\n"
    )

    f.write(
        f"Frozen core identities: "
        f"{len(CORE_IDENTITIES)}\n"
    )

    f.write(
        "GSE282122 core KEEP: 18/18\n"
    )

    f.write(
        "SCP1884 core KEEP: 18/18\n"
    )

    f.write(
        "SCP259 expected core KEEP: "
        "16/18\n"
    )

    f.write(
        "SCP259 expected absent "
        "tissue-restricted identities: "
        "Ileal_absorptive, Paneth\n\n"
    )

    f.write(
        "Candidate identities\n"
    )

    f.write(
        "-" * 72
        + "\n"
    )

    for _, row in (
        candidate_summary.iterrows()
    ):

        f.write(
            f"{row['identity']}: "
            f"{row['n_datasets_present']}/3 "
            "datasets, "
            f"{int(row['total_cells_across_datasets']):,} "
            "cells, "
            f"decision="
            f"{row['harmonization_decision']}\n"
        )

    f.write("\n")

    f.write(
        "Interpretation\n"
    )

    f.write(
        "-" * 72
        + "\n"
    )

    f.write(
        "1. Existing dataset-level audits "
        "are reused as the canonical baseline.\n"
    )

    f.write(
        "2. The frozen 18-class core "
        "taxonomy is retained.\n"
    )

    f.write(
        "3. Absence of Ileal_absorptive "
        "and Paneth from SCP259 is treated "
        "as tissue/sampling structure, "
        "not taxonomy failure.\n"
    )

    f.write(
        "4. Glial and Tuft are not "
        "automatically promoted to KEEP.\n"
    )

    f.write(
        "5. Glial/Tuft proceed to "
        "expression-level cross-dataset "
        "validation before final ontology "
        "freeze.\n"
    )

    f.write(
        "6. Other HOLD identities remain "
        "outside the primary reference "
        "unless separately justified.\n"
    )


# ============================================================
# 17. Console report
# ============================================================

print()
print("=" * 72)
print("CORE IDENTITY SUMMARY")
print("=" * 72)

print(
    core_summary[
        [
            "identity",
            "n_datasets_present",
            "total_cells_across_datasets",
            "harmonization_decision",
        ]
    ].to_string(
        index=False
    )
)


print()
print("=" * 72)
print("CANDIDATE IDENTITY SUMMARY")
print("=" * 72)

print(
    candidate_summary[
        [
            "identity",
            "n_datasets_present",
            "total_cells_across_datasets",
            "sum_dataset_donors",
            "harmonization_decision",
        ]
    ].to_string(
        index=False
    )
)


print()
print("=" * 72)
print("OTHER HOLD IDENTITIES")
print("=" * 72)

if len(other_holds) == 0:

    print("None")

else:

    print(
        other_holds[
            [
                "identity",
                "n_datasets_present",
                "total_cells_across_datasets",
            ]
        ].to_string(
            index=False
        )
    )


print()
print("=" * 72)
print("OUTPUT DIRECTORY")
print("=" * 72)

print(OUTPUT_DIR)


print()
print("=" * 72)
print("SUCCESS")
print("=" * 72)

print(
    "Cross-dataset reference "
    "harmonization audit completed."
)

print(
    "The 18-class core can proceed "
    "toward balanced sampling."
)

print(
    "Glial and Tuft require only the "
    "final expression-level validation."
)