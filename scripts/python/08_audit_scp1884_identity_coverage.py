#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


PROJECT = Path(
    "/home/mazekai/IBD_EcoTyper"
)

ENV_ROOT = Path(
    "/home/mazekai/miniconda3/envs/ibd_refaudit"
)

SOURCE = (
    PROJECT
    / "03_reference/SCP1884/audit_v1/output/"
    "02_candidate_taxonomy/"
    "scp1884_candidate_cells_provisional_taxonomy.parquet"
)

OUT = (
    PROJECT
    / "03_reference/SCP1884/audit_v1/output/"
    "03_identity_coverage"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


EXPECTED_PRIMARY = [
    "B",
    "Plasma",
    "CD4_T",
    "CD8_T",
    "NK",
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
    "ILC",
]


def fail(msg):
    raise RuntimeError(msg)


# ============================================================
# 0. Environment
# ============================================================

python_exe = Path(
    sys.executable
).resolve()

if ENV_ROOT not in python_exe.parents:
    fail(
        f"Wrong environment: {python_exe}"
    )


# ============================================================
# 1. Load candidate taxonomy
# ============================================================

df = pd.read_parquet(
    SOURCE
)

required = [
    "cell_id",
    "PubID",
    "biological_sample_id",
    "disease",
    "inflammation",
    "Site",
    "provisional_coarse_identity",
    "primary_reference_action",
]

missing = [
    x
    for x in required
    if x not in df.columns
]

if missing:
    fail(
        f"Missing columns: {missing}"
    )

if len(df) != 720_633:
    fail(
        f"Unexpected cell count: {len(df):,}"
    )


# ============================================================
# 2. Primary-reference cells
# ============================================================

primary = (
    df.loc[
        df[
            "primary_reference_action"
        ].eq(
            "KEEP"
        )
    ]
    .copy()
)

if len(primary) != 667_483:
    fail(
        f"KEEP cells={len(primary):,}; "
        "expected=667,483."
    )

observed_primary = set(
    primary[
        "provisional_coarse_identity"
    ].unique()
)

expected_primary = set(
    EXPECTED_PRIMARY
)

missing_primary = (
    expected_primary
    - observed_primary
)

unexpected_primary = (
    observed_primary
    - expected_primary
)

if missing_primary:
    fail(
        "Missing primary identities: "
        + ", ".join(
            sorted(
                missing_primary
            )
        )
    )

if unexpected_primary:
    fail(
        "Unexpected KEEP identities: "
        + ", ".join(
            sorted(
                unexpected_primary
            )
        )
    )


# ============================================================
# 3. Identity-level donor balance
# ============================================================

rows = []

for identity, sub in primary.groupby(
    "provisional_coarse_identity",
    observed=True,
):

    donor_counts = (
        sub[
            "PubID"
        ]
        .value_counts()
    )

    total = int(
        donor_counts.sum()
    )

    donor_frac = (
        donor_counts
        / total
    )

    rows.append(
        {
            "identity":
                identity,

            "n_cells":
                len(sub),

            "n_samples":
                sub[
                    "biological_sample_id"
                ].nunique(),

            "n_donors":
                sub[
                    "PubID"
                ].nunique(),

            "median_cells_per_donor":
                float(
                    donor_counts.median()
                ),

            "min_cells_per_donor":
                int(
                    donor_counts.min()
                ),

            "max_cells_per_donor":
                int(
                    donor_counts.max()
                ),

            "max_donor_fraction":
                float(
                    donor_frac.max()
                ),

            "donor_HHI":
                float(
                    (
                        donor_frac ** 2
                    ).sum()
                ),

            "n_CO_cells":
                int(
                    sub[
                        "Site"
                    ]
                    .eq(
                        "CO"
                    )
                    .sum()
                ),

            "n_TI_cells":
                int(
                    sub[
                        "Site"
                    ]
                    .eq(
                        "TI"
                    )
                    .sum()
                ),

            "n_SB_cells":
                int(
                    sub[
                        "Site"
                    ]
                    .eq(
                        "SB"
                    )
                    .sum()
                ),

            "n_CD_cells":
                int(
                    sub[
                        "disease"
                    ]
                    .eq(
                        "CD"
                    )
                    .sum()
                ),

            "n_Healthy_cells":
                int(
                    sub[
                        "disease"
                    ]
                    .eq(
                        "Healthy"
                    )
                    .sum()
                ),
        }
    )


coverage = pd.DataFrame(
    rows
)

coverage[
    "identity"
] = pd.Categorical(
    coverage[
        "identity"
    ],
    categories=EXPECTED_PRIMARY,
    ordered=True,
)

coverage = (
    coverage
    .sort_values(
        "identity"
    )
    .reset_index(
        drop=True
    )
)

coverage.to_csv(
    OUT
    / "primary_identity_coverage.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 4. Identity x site
# ============================================================

identity_site = (
    primary.groupby(
        [
            "provisional_coarse_identity",
            "Site",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),
        n_samples=(
            "biological_sample_id",
            "nunique",
        ),
        n_donors=(
            "PubID",
            "nunique",
        ),
    )
    .reset_index()
)

identity_site.to_csv(
    OUT
    / "primary_identity_by_site.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 5. Identity x disease / inflammation
# ============================================================

identity_condition = (
    primary.groupby(
        [
            "provisional_coarse_identity",
            "disease",
            "inflammation",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),
        n_samples=(
            "biological_sample_id",
            "nunique",
        ),
        n_donors=(
            "PubID",
            "nunique",
        ),
    )
    .reset_index()
)

identity_condition.to_csv(
    OUT
    / "primary_identity_by_condition.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 6. HOLD / EXCLUDE audit
# ============================================================

non_primary = (
    df.loc[
        ~df[
            "primary_reference_action"
        ].eq(
            "KEEP"
        )
    ]
    .groupby(
        [
            "provisional_coarse_identity",
            "primary_reference_action",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),
        n_samples=(
            "biological_sample_id",
            "nunique",
        ),
        n_donors=(
            "PubID",
            "nunique",
        ),
    )
    .reset_index()
    .sort_values(
        "n_cells",
        ascending=False,
    )
)

non_primary.to_csv(
    OUT
    / "hold_exclude_coverage.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 7. Simple audit flags
#
# These are descriptive QC flags, not automatic biological
# exclusion criteria.
# ============================================================

coverage[
    "coverage_flag"
] = "PASS"

coverage.loc[
    coverage[
        "n_donors"
    ] < 20,
    "coverage_flag",
] = "REVIEW_DONOR_COVERAGE"

coverage.loc[
    coverage[
        "max_donor_fraction"
    ] > 0.20,
    "coverage_flag",
] = "REVIEW_DONOR_DOMINANCE"


# Paneth remains a deliberate biological review target,
# independent of donor coverage.

coverage.loc[
    coverage[
        "identity"
    ].astype(str).eq(
        "Paneth"
    ),
    "coverage_flag",
] = "REVIEW_CROSS_DATASET"


coverage.to_csv(
    OUT
    / "primary_identity_coverage_with_flags.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 8. Console report
# ============================================================

print()
print("=" * 78)
print("SCP1884 PHASE 3 IDENTITY COVERAGE")
print("=" * 78)

print()
print(
    f"Candidate cells       : {len(df):,}"
)
print(
    f"Primary KEEP cells    : {len(primary):,}"
)
print(
    f"Primary identities    : {len(observed_primary)} / 18"
)

print()
print(
    coverage[
        [
            "identity",
            "n_cells",
            "n_samples",
            "n_donors",
            "median_cells_per_donor",
            "max_donor_fraction",
            "donor_HHI",
            "coverage_flag",
        ]
    ].to_string(
        index=False
    )
)

print()
print("HOLD / EXCLUDE:")
print(
    non_primary.to_string(
        index=False
    )
)

print()
print(
    "Output directory:"
)
print(
    OUT
)
