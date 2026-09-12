#!/usr/bin/env python3

from pathlib import Path
import sys
import importlib.metadata as md

import numpy as np
import pandas as pd


# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT = Path("/home/mazekai/IBD_EcoTyper")

EXPECTED_ENV_ROOT = Path(
    "/home/mazekai/miniconda3/envs/ibd_refaudit"
)

CANDIDATE_FILE = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "01_candidate_pool"
    / "candidate_cells.parquet"
)

MAPPING_FILE = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "03_taxonomy_mapping_v1"
    / "taxonomy_mapping_v1.tsv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "04_identity_coverage_v1"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "identity_coverage_v1_summary.txt"
)

PATH_COLS = [
    "bucket",
    "sub_bucket",
    "major",
    "minor",
    "final_analysis",
]

CORE_IDENTITIES = [
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

EXPECTED_CANDIDATE_CELLS = 505_591
EXPECTED_PRIMARY_CELLS = 484_979


# ============================================================
# Helpers
# ============================================================

def package_version(name):
    try:
        return md.version(name)
    except md.PackageNotFoundError:
        return "NOT INSTALLED"


def write_summary(lines):
    SUMMARY_FILE.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def fail(message, summary):
    summary.append(f"[FAIL] {message}")
    write_summary(summary)
    raise RuntimeError(message)


def calc_balance_table(df, group_col):
    """
    Per-identity donor/sample balance metrics.
    """

    rows = []

    for identity, sub in df.groupby(
        group_col,
        observed=True,
    ):

        donor_counts = (
            sub.groupby(
                "Patient",
                observed=True,
            )
            .size()
        )

        donor_fraction = (
            donor_counts
            / donor_counts.sum()
        )

        sample_counts = (
            sub.groupby(
                "sample_id",
                observed=True,
            )
            .size()
        )

        rows.append(
            {
                group_col: identity,

                "n_cells":
                    len(sub),

                "n_samples":
                    sub["sample_id"].nunique(),

                "n_donors":
                    sub["Patient"].nunique(),

                "median_cells_per_donor":
                    float(donor_counts.median()),

                "min_cells_per_donor":
                    int(donor_counts.min()),

                "max_cells_per_donor":
                    int(donor_counts.max()),

                "q25_cells_per_donor":
                    float(
                        donor_counts.quantile(0.25)
                    ),

                "q75_cells_per_donor":
                    float(
                        donor_counts.quantile(0.75)
                    ),

                "max_donor_fraction":
                    float(
                        donor_fraction.max()
                    ),

                "donor_HHI":
                    float(
                        np.square(
                            donor_fraction
                        ).sum()
                    ),

                "median_cells_per_sample":
                    float(sample_counts.median()),

                "min_cells_per_sample":
                    int(sample_counts.min()),

                "max_cells_per_sample":
                    int(sample_counts.max()),
            }
        )

    out = pd.DataFrame(rows)

    out["max_donor_fraction"] = (
        out["max_donor_fraction"]
        .round(4)
    )

    out["donor_HHI"] = (
        out["donor_HHI"]
        .round(4)
    )

    return out


def composition_table(
    df,
    identity_col,
    category_col,
    prefix,
):
    """
    Wide table:
    one row per identity,
    counts for each value of category_col.
    """

    tab = (
        df.groupby(
            [
                identity_col,
                category_col,
            ],
            observed=True,
            dropna=False,
        )
        .size()
        .unstack(
            fill_value=0
        )
        .reset_index()
    )

    tab.columns.name = None

    rename = {}

    for col in tab.columns:

        if col == identity_col:
            continue

        rename[col] = (
            f"{prefix}_{col}"
            .replace(" ", "_")
        )

    return tab.rename(
        columns=rename
    )


def coverage_long(
    df,
    identity_col,
    category_col,
):
    """
    Long-format coverage table including cells,
    samples and donors.
    """

    return (
        df.groupby(
            [
                identity_col,
                category_col,
            ],
            observed=True,
            dropna=False,
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
            n_donors=(
                "Patient",
                "nunique",
            ),
        )
        .reset_index()
    )


# ============================================================
# 0. Environment preflight
# ============================================================

summary = []

summary.append(
    "GSE282122 Reference Audit v2"
)
summary.append(
    "Phase 3A: Identity coverage and balance audit"
)
summary.append("=" * 78)
summary.append("")

python_exe = Path(
    sys.executable
).resolve()

summary.append("Environment")
summary.append("-" * 78)

summary.append(
    f"Python executable: {python_exe}"
)

if EXPECTED_ENV_ROOT not in python_exe.parents:
    fail(
        "Script is not running inside ibd_refaudit.",
        summary,
    )

summary.append(
    "[PASS] Correct Conda environment detected."
)
summary.append(
    f"Python: {sys.version.split()[0]}"
)
summary.append(
    f"pandas: {package_version('pandas')}"
)
summary.append(
    f"numpy: {package_version('numpy')}"
)
summary.append("")


# ============================================================
# 1. Load inputs
# ============================================================

for f in [
    CANDIDATE_FILE,
    MAPPING_FILE,
]:
    if not f.is_file():
        fail(
            f"Missing input: {f}",
            summary,
        )

candidate = pd.read_parquet(
    CANDIDATE_FILE
)

mapping = pd.read_csv(
    MAPPING_FILE,
    sep="\t",
)

summary.append("Inputs")
summary.append("-" * 78)

summary.append(
    f"Candidate cells: {len(candidate):,}"
)
summary.append(
    f"Mapping paths: {len(mapping):,}"
)
summary.append("")


# ============================================================
# 2. Input checkpoints
# ============================================================

if len(candidate) != EXPECTED_CANDIDATE_CELLS:
    fail(
        "Candidate pool does not contain "
        f"{EXPECTED_CANDIDATE_CELLS:,} cells.",
        summary,
    )

if len(mapping) != 110:
    fail(
        "Taxonomy mapping does not contain 110 paths.",
        summary,
    )

if candidate["cell_id"].duplicated().any():
    fail(
        "Candidate cell_id is not unique.",
        summary,
    )

required_mapping_cols = (
    PATH_COLS
    + [
        "path_id",
        "coarse_identity",
        "identity_status",
        "core_taxonomy_member",
        "primary_reference_action",
        "state_flag",
        "mapping_rationale",
        "review_notes",
    ]
)

missing_mapping = [
    c
    for c in required_mapping_cols
    if c not in mapping.columns
]

if missing_mapping:
    fail(
        "Mapping columns missing: "
        + ", ".join(missing_mapping),
        summary,
    )


# ============================================================
# 3. Ensure annotation paths are unique in mapping
# ============================================================

if mapping.duplicated(
    subset=PATH_COLS
).any():

    dup = mapping.loc[
        mapping.duplicated(
            subset=PATH_COLS,
            keep=False,
        ),
        PATH_COLS + ["path_id"],
    ]

    dup.to_csv(
        OUTPUT_DIR
        / "ERROR_duplicate_mapping_paths.tsv",
        sep="\t",
        index=False,
    )

    fail(
        "Mapping contains duplicate annotation paths.",
        summary,
    )

summary.append("Mapping validation")
summary.append("-" * 78)
summary.append(
    "[PASS] 110 unique annotation paths."
)
summary.append("")


# ============================================================
# 4. Apply mapping to every candidate cell
# ============================================================

mapping_cols = (
    PATH_COLS
    + [
        "path_id",
        "coarse_identity",
        "identity_status",
        "core_taxonomy_member",
        "primary_reference_action",
        "state_flag",
        "mapping_rationale",
        "review_notes",
    ]
)

mapped = candidate.merge(
    mapping[
        mapping_cols
    ],
    on=PATH_COLS,
    how="left",
    validate="many_to_one",
)


# ============================================================
# 5. Mapping integrity
# ============================================================

unmapped = mapped[
    "path_id"
].isna()

if unmapped.any():

    (
        mapped.loc[
            unmapped,
            PATH_COLS,
        ]
        .drop_duplicates()
        .to_csv(
            OUTPUT_DIR
            / "ERROR_unmapped_annotation_paths.tsv",
            sep="\t",
            index=False,
        )
    )

    fail(
        f"{int(unmapped.sum()):,} cells "
        "failed taxonomy mapping.",
        summary,
    )

if len(mapped) != len(candidate):
    fail(
        "Row count changed during taxonomy mapping.",
        summary,
    )

if not mapped["cell_id"].is_unique:
    fail(
        "cell_id uniqueness lost after mapping.",
        summary,
    )

summary.append("Cell-level mapping")
summary.append("-" * 78)

summary.append(
    f"[PASS] {len(mapped):,}/{len(candidate):,} "
    "cells mapped."
)

summary.append(
    "[PASS] No candidate cell duplicated or lost."
)
summary.append("")


# ============================================================
# 6. Save all mapped candidate metadata
# ============================================================

mapped.to_parquet(
    OUTPUT_DIR
    / "mapped_candidate_cells_v1.parquet",
    index=False,
    engine="pyarrow",
)


# ============================================================
# 7. Build primary KEEP reference metadata
# ============================================================

primary = (
    mapped.loc[
        mapped[
            "primary_reference_action"
        ] == "KEEP"
    ]
    .copy()
)

if len(primary) != EXPECTED_PRIMARY_CELLS:
    fail(
        "Primary KEEP pool mismatch: "
        f"expected {EXPECTED_PRIMARY_CELLS:,}, "
        f"observed {len(primary):,}.",
        summary,
    )

observed_core = set(
    primary[
        "coarse_identity"
    ]
    .dropna()
    .astype(str)
    .unique()
)

expected_core = set(
    CORE_IDENTITIES
)

missing_core = (
    expected_core
    - observed_core
)

extra_primary = (
    observed_core
    - expected_core
)

if missing_core:
    fail(
        "Core identities missing from KEEP pool: "
        + ", ".join(
            sorted(missing_core)
        ),
        summary,
    )

if extra_primary:
    fail(
        "Unexpected identities present in KEEP pool: "
        + ", ".join(
            sorted(extra_primary)
        ),
        summary,
    )

primary.to_parquet(
    OUTPUT_DIR
    / "primary_reference_cells_v1.parquet",
    index=False,
    engine="pyarrow",
)

summary.append("Primary reference candidate")
summary.append("-" * 78)

summary.append(
    f"KEEP cells: {len(primary):,}"
)

summary.append(
    f"Core identities: {len(observed_core)}"
)

summary.append(
    "[PASS] All 18 core identities are represented."
)
summary.append("")


# ============================================================
# 8. Identity balance metrics
# ============================================================

identity_balance = calc_balance_table(
    primary,
    "coarse_identity",
)

identity_balance.to_csv(
    OUTPUT_DIR
    / "identity_balance_primary_keep.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 9. Disease composition
# ============================================================

disease_long = coverage_long(
    primary,
    "coarse_identity",
    "Disease",
)

disease_long.to_csv(
    OUTPUT_DIR
    / "identity_coverage_by_disease.tsv",
    sep="\t",
    index=False,
)

disease_wide = composition_table(
    primary,
    "coarse_identity",
    "Disease",
    "n_cells",
)


# ============================================================
# 10. Site composition
# ============================================================

site_long = coverage_long(
    primary,
    "coarse_identity",
    "Site",
)

site_long.to_csv(
    OUTPUT_DIR
    / "identity_coverage_by_site.tsv",
    sep="\t",
    index=False,
)

site_wide = composition_table(
    primary,
    "coarse_identity",
    "Site",
    "n_cells",
)


# ============================================================
# 11. Inflammation composition
# ============================================================

inflammation_long = coverage_long(
    primary,
    "coarse_identity",
    "Inflammation",
)

inflammation_long.to_csv(
    OUTPUT_DIR
    / "identity_coverage_by_inflammation.tsv",
    sep="\t",
    index=False,
)

inflammation_wide = composition_table(
    primary,
    "coarse_identity",
    "Inflammation",
    "n_cells",
)


# ============================================================
# 12. Disease x site composition
# ============================================================

disease_site = (
    primary.groupby(
        [
            "coarse_identity",
            "Disease",
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
            "sample_id",
            "nunique",
        ),
        n_donors=(
            "Patient",
            "nunique",
        ),
    )
    .reset_index()
)

disease_site.to_csv(
    OUTPUT_DIR
    / "identity_coverage_by_disease_site.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 13. Donor-level cell counts
# ============================================================

donor_counts = (
    primary.groupby(
        [
            "coarse_identity",
            "Patient",
            "Disease",
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
    )
    .reset_index()
)

donor_total = (
    donor_counts.groupby(
        "coarse_identity",
        observed=True,
    )["n_cells"]
    .transform("sum")
)

donor_counts[
    "donor_fraction"
] = (
    donor_counts["n_cells"]
    / donor_total
)

donor_counts[
    "donor_fraction"
] = (
    donor_counts[
        "donor_fraction"
    ]
    .round(6)
)

donor_counts.to_csv(
    OUTPUT_DIR
    / "identity_counts_by_donor.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 14. Sample-level cell counts
# ============================================================

sample_counts = (
    primary.groupby(
        [
            "coarse_identity",
            "sample_id",
            "Patient",
            "Disease",
            "Site",
            "Inflammation",
        ],
        observed=True,
    )
    .size()
    .rename(
        "n_cells"
    )
    .reset_index()
)

sample_counts.to_csv(
    OUTPUT_DIR
    / "identity_counts_by_sample.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 15. Dataset sampling opportunity table
#
# Important for tissue-aware interpretation.
# ============================================================

sampling_opportunity = (
    primary[
        [
            "sample_id",
            "Patient",
            "Disease",
            "Site",
            "Inflammation",
        ]
    ]
    .drop_duplicates()
    .groupby(
        [
            "Disease",
            "Site",
            "Inflammation",
        ],
        observed=True,
    )
    .agg(
        n_samples=(
            "sample_id",
            "nunique",
        ),
        n_donors=(
            "Patient",
            "nunique",
        ),
    )
    .reset_index()
)

sampling_opportunity.to_csv(
    OUTPUT_DIR
    / "sampling_opportunity_by_disease_site_inflammation.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 16. Combined identity overview
# ============================================================

overview = (
    identity_balance
    .merge(
        disease_wide,
        on="coarse_identity",
        how="left",
        validate="one_to_one",
    )
    .merge(
        site_wide,
        on="coarse_identity",
        how="left",
        validate="one_to_one",
    )
    .merge(
        inflammation_wide,
        on="coarse_identity",
        how="left",
        validate="one_to_one",
    )
)

overview = overview.sort_values(
    "n_cells",
    ascending=False,
)

overview.to_csv(
    OUTPUT_DIR
    / "identity_coverage_overview_primary_keep.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 17. Action-level audit on ALL candidate cells
# ============================================================

action_audit = (
    mapped.groupby(
        [
            "primary_reference_action",
            "coarse_identity",
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
        n_donors=(
            "Patient",
            "nunique",
        ),
    )
    .reset_index()
)

action_audit.to_csv(
    OUTPUT_DIR
    / "all_candidate_cells_by_identity_action.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 18. Console summary
# ============================================================

summary.append(
    "Primary identity overview"
)
summary.append("-" * 78)

for row in overview.itertuples(
    index=False
):

    summary.append(
        f"{row.coarse_identity}: "
        f"{row.n_cells:,} cells, "
        f"{row.n_samples} samples, "
        f"{row.n_donors} donors, "
        f"max donor fraction="
        f"{row.max_donor_fraction:.4f}, "
        f"HHI={row.donor_HHI:.4f}"
    )

summary.append("")

summary.append(
    "Interpretation note"
)
summary.append("-" * 78)

summary.append(
    "No Green/Yellow/Red quality tier is assigned "
    "in Phase 3A."
)

summary.append(
    "Tissue-restricted identities must be evaluated "
    "against appropriate sampling opportunities rather "
    "than all 41 donors or all 116 biopsies."
)

summary.append(
    "Final quality tiers will be assigned only after "
    "reviewing disease, site, inflammation and donor "
    "coverage together."
)

summary.append("")

summary.append(
    "Final Phase-3A status"
)
summary.append("-" * 78)

summary.append(
    "[PASS] 505,591 candidate cells mapped."
)
summary.append(
    "[PASS] 484,979 primary KEEP cells retained."
)
summary.append(
    "[PASS] All 18 core identities represented."
)
summary.append(
    "[PASS] Donor balance metrics calculated."
)
summary.append(
    "[PASS] Disease/site/inflammation coverage calculated."
)
summary.append(
    "[PASS] Tissue-aware sampling opportunity table generated."
)
summary.append(
    "[PASS] No expression matrix was accessed."
)

write_summary(
    summary
)


# ============================================================
# Console output
# ============================================================

print()
print("=" * 78)
print(
    "GSE282122 Phase 3A identity coverage audit completed"
)
print("=" * 78)

print()
print(
    f"Mapped candidate cells : {len(mapped):,}"
)
print(
    f"Primary KEEP cells     : {len(primary):,}"
)
print(
    f"Core identities        : {len(observed_core)}"
)

print()
print(
    overview[
        [
            "coarse_identity",
            "n_cells",
            "n_samples",
            "n_donors",
            "median_cells_per_donor",
            "max_donor_fraction",
            "donor_HHI",
        ]
    ].to_string(
        index=False
    )
)

print()
print("Outputs:")
print(
    f"  {OUTPUT_DIR}"
)

print()
print("Summary:")
print(
    f"  {SUMMARY_FILE}"
)