#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd
import numpy as np


# ============================================================
# Configuration
# ============================================================

ROOT = Path("/home/mazekai/IBD_EcoTyper")

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
# Canonical upstream audited metadata
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


# ============================================================
# Expected upstream checkpoints
# ============================================================

EXPECTED_N = {
    "GSE282122": 505_591,
    "SCP1884": 720_633,
    "SCP259": 365_492,
}

EXPECTED_TOTAL = sum(
    EXPECTED_N.values()
)


# ============================================================
# Current reference ontology
#
# IMPORTANT:
# These are NOT treated as the v2 answer.
# They are recorded only so the audit can tell us what the
# previous pipeline kept.
# ============================================================

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


# Current downstream EcoTyper CORE12.
# Again: informational only.
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


ANATOMY_ENCODED_IDENTITIES = {
    "Ileal_absorptive",
    "Colonic_absorptive",
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


def pick_col(df, candidates, required=True):
    """
    Return first existing column name.
    """

    for x in candidates:
        if x in df.columns:
            return x

    if required:
        fail(
            "None of these columns exists: "
            + ", ".join(candidates)
        )

    return None


def as_string(series):
    """
    Convert to explicit audit-safe strings.
    """

    return (
        series
        .astype("string")
        .fillna("<NA>")
        .replace("", "<EMPTY>")
    )


def optional_series(df, candidates, default=""):
    col = pick_col(
        df,
        candidates,
        required=False,
    )

    if col is None:
        return pd.Series(
            [default] * len(df),
            index=df.index,
            dtype="string",
        )

    return as_string(
        df[col]
    )


def join_annotation_path(df, columns):
    """
    Combine hierarchical source annotations.
    """

    existing = [
        x for x in columns
        if x in df.columns
    ]

    if len(existing) == 0:
        return pd.Series(
            [""] * len(df),
            index=df.index,
            dtype="string",
        )

    parts = [
        as_string(df[x])
        for x in existing
    ]

    result = parts[0]

    for x in parts[1:]:
        result = (
            result
            + " > "
            + x
        )

    return result


def count_string(series):
    """
    Compact composition:
        CD=100; UC=80; Healthy=20
    """

    s = (
        series
        .astype("string")
        .fillna("<NA>")
    )

    tab = (
        s.value_counts(
            dropna=False
        )
    )

    return "; ".join(
        f"{idx}={int(n)}"
        for idx, n in tab.items()
    )


def unique_string(series):
    vals = sorted(
        set(
            series
            .astype("string")
            .fillna("<NA>")
        )
    )

    return " | ".join(vals)


def donor_balance(sub):
    """
    Basic concentration metrics.
    """

    counts = (
        sub["global_donor_id"]
        .value_counts()
    )

    if len(counts) == 0:
        return np.nan, np.nan

    frac = (
        counts
        / counts.sum()
    )

    return (
        float(frac.max()),
        float(
            np.square(frac).sum()
        ),
    )


# ============================================================
# GSE282122
# ============================================================

def load_gse282122():

    x = pd.read_parquet(
        GSE_FILE
    )

    if len(x) != EXPECTED_N["GSE282122"]:
        fail(
            "GSE282122 row count mismatch: "
            f"{len(x):,}"
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

    missing = [
        c for c in required
        if c not in x.columns
    ]

    if missing:
        fail(
            "GSE282122 missing columns: "
            + ", ".join(missing)
        )

    out = pd.DataFrame({
        "dataset":
            "GSE282122",

        "cell_id":
            as_string(
                x["cell_id"]
            ),

        "donor_id":
            as_string(
                x["Patient"]
            ),

        "sample_id":
            as_string(
                x["sample_id"]
            ),

        "disease":
            as_string(
                x["Disease"]
            ),

        "site":
            as_string(
                x["Site"]
            ),

        "inflammation":
            as_string(
                x["Inflammation"]
            ),

        # Original author annotation hierarchy
        "author_major_annotation":
            optional_series(
                x,
                ["major", "bucket"],
            ),

        "author_fine_annotation":
            optional_series(
                x,
                ["final_analysis", "minor"],
            ),

        "author_annotation_path":
            join_annotation_path(
                x,
                [
                    "bucket",
                    "sub_bucket",
                    "major",
                    "minor",
                    "final_analysis",
                ],
            ),

        "portal_or_secondary_annotation":
            "",

        # Current harmonized taxonomy
        "coarse_identity":
            as_string(
                x["coarse_identity"]
            ),

        "primary_reference_action":
            as_string(
                x[
                    "primary_reference_action"
                ]
            ),

        "state_flag":
            optional_series(
                x,
                ["state_flag"],
            ),

        "mapping_rationale":
            optional_series(
                x,
                ["mapping_rationale"],
            ),

        "review_notes":
            optional_series(
                x,
                ["review_notes"],
            ),
    })

    return out


# ============================================================
# SCP1884
# ============================================================

def load_scp1884():

    x = pd.read_parquet(
        SCP1884_FILE
    )

    if len(x) != EXPECTED_N["SCP1884"]:
        fail(
            "SCP1884 row count mismatch: "
            f"{len(x):,}"
        )

    donor_col = pick_col(
        x,
        [
            "v2_donor_id",
            "PubID",
        ],
    )

    site_col = pick_col(
        x,
        [
            "site_canonical",
            "Site",
        ],
    )

    required = [
        "cell_id",
        "biological_sample_id",
        "disease",
        "inflammation",
        "provisional_coarse_identity",
        "primary_reference_action",
        "anno2",
    ]

    missing = [
        c for c in required
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
            as_string(
                x["cell_id"]
            ),

        "donor_id":
            as_string(
                x[donor_col]
            ),

        "sample_id":
            as_string(
                x[
                    "biological_sample_id"
                ]
            ),

        "disease":
            as_string(
                x["disease"]
            ),

        "site":
            as_string(
                x[site_col]
            ),

        "inflammation":
            as_string(
                x["inflammation"]
            ),

        "author_major_annotation":
            optional_series(
                x,
                ["anno_overall"],
            ),

        "author_fine_annotation":
            as_string(
                x["anno2"]
            ),

        "author_annotation_path":
            join_annotation_path(
                x,
                [
                    "anno_overall",
                    "anno2",
                ],
            ),

        "portal_or_secondary_annotation":
            optional_series(
                x,
                [
                    "v2_celltype",
                    "Celltype",
                ],
            ),

        "coarse_identity":
            as_string(
                x[
                    "provisional_coarse_identity"
                ]
            ),

        "primary_reference_action":
            as_string(
                x[
                    "primary_reference_action"
                ]
            ),

        "state_flag":
            optional_series(
                x,
                ["state_flag"],
            ),

        "mapping_rationale":
            optional_series(
                x,
                ["mapping_rationale"],
            ),

        "review_notes":
            optional_series(
                x,
                ["review_notes"],
            ),
    })

    return out


# ============================================================
# SCP259
# ============================================================

def load_scp259():

    x = pd.read_parquet(
        SCP259_FILE
    )

    if len(x) != EXPECTED_N["SCP259"]:
        fail(
            "SCP259 row count mismatch: "
            f"{len(x):,}"
        )

    cell_col = pick_col(
        x,
        [
            "cell_id",
            "NAME",
        ],
    )

    author_col = pick_col(
        x,
        [
            "author_annotation",
            "Cluster",
        ],
    )

    required = [
        "donor_id",
        "sample_id",
        "disease",
        "site",
        "inflammation",
        "coarse_identity",
        "primary_reference_action",
    ]

    missing = [
        c for c in required
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
            as_string(
                x[cell_col]
            ),

        "donor_id":
            as_string(
                x["donor_id"]
            ),

        "sample_id":
            as_string(
                x["sample_id"]
            ),

        "disease":
            as_string(
                x["disease"]
            ),

        "site":
            as_string(
                x["site"]
            ),

        "inflammation":
            as_string(
                x["inflammation"]
            ),

        # SCP259 uses author Cluster as its fine annotation.
        # tissue_fraction is useful as the higher-level context.
        "author_major_annotation":
            optional_series(
                x,
                [
                    "tissue_fraction",
                    "expression_bundle",
                ],
            ),

        "author_fine_annotation":
            as_string(
                x[author_col]
            ),

        "author_annotation_path":
            (
                optional_series(
                    x,
                    [
                        "tissue_fraction",
                        "expression_bundle",
                    ],
                )
                + " > "
                + as_string(
                    x[author_col]
                )
            ),

        "portal_or_secondary_annotation":
            "",

        "coarse_identity":
            as_string(
                x["coarse_identity"]
            ),

        "primary_reference_action":
            as_string(
                x[
                    "primary_reference_action"
                ]
            ),

        "state_flag":
            optional_series(
                x,
                ["state_flag"],
            ),

        "mapping_rationale":
            optional_series(
                x,
                ["mapping_rationale"],
            ),

        "review_notes":
            optional_series(
                x,
                ["review_notes"],
            ),
    })

    return out


# ============================================================
# Load + standardize
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
    "IBD SIGNATURE TAXONOMY V2 AUDIT"
)
print(
    "=" * 78
)


tables = [
    load_gse282122(),
    load_scp1884(),
    load_scp259(),
]


for x in tables:
    name = x["dataset"].iloc[0]

    print(
        f"{name}: "
        f"{len(x):,} candidate cells"
    )


cells = pd.concat(
    tables,
    axis=0,
    ignore_index=True,
)


if len(cells) != EXPECTED_TOTAL:
    fail(
        "Combined candidate-cell count mismatch: "
        f"{len(cells):,} != "
        f"{EXPECTED_TOTAL:,}"
    )


# ============================================================
# Global identifiers
#
# Important because donor/sample IDs may theoretically overlap
# across independent studies.
# ============================================================

cells["global_cell_id"] = (
    cells["dataset"]
    + "::"
    + cells["cell_id"]
)

cells["global_donor_id"] = (
    cells["dataset"]
    + "::"
    + cells["donor_id"]
)

cells["global_sample_id"] = (
    cells["dataset"]
    + "::"
    + cells["sample_id"]
)


if cells[
    "global_cell_id"
].duplicated().any():

    fail(
        "Duplicated global cell IDs "
        "after dataset combination."
    )


# ============================================================
# Informational flags
# ============================================================

cells[
    "current_reference20_member"
] = (
    cells[
        "coarse_identity"
    ].isin(
        CURRENT_REFERENCE20
    )
)

cells[
    "current_ecotyper_core12_member"
] = (
    cells[
        "coarse_identity"
    ].isin(
        CURRENT_ECOTYPER_CORE12
    )
)

cells[
    "anatomy_encoded_identity"
] = (
    cells[
        "coarse_identity"
    ].isin(
        ANATOMY_ENCODED_IDENTITIES
    )
)


# ============================================================
# 1. SOURCE ANNOTATION AUDIT
#
# One row per:
# dataset × original author annotation × current mapping
# ============================================================

SOURCE_GROUP = [
    "dataset",
    "author_major_annotation",
    "author_fine_annotation",
    "author_annotation_path",
    "portal_or_secondary_annotation",
    "coarse_identity",
    "primary_reference_action",
    "state_flag",
    "current_reference20_member",
    "current_ecotyper_core12_member",
    "anatomy_encoded_identity",
]


source_rows = []


for keys, sub in cells.groupby(
    SOURCE_GROUP,
    observed=True,
    dropna=False,
):

    row = dict(
        zip(
            SOURCE_GROUP,
            keys,
        )
    )

    max_frac, hhi = (
        donor_balance(sub)
    )

    row.update({
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

        "max_donor_fraction":
            max_frac,

        "donor_HHI":
            hhi,

        "disease_levels":
            unique_string(
                sub["disease"]
            ),

        "disease_cell_counts":
            count_string(
                sub["disease"]
            ),

        "site_levels":
            unique_string(
                sub["site"]
            ),

        "site_cell_counts":
            count_string(
                sub["site"]
            ),

        "inflammation_levels":
            unique_string(
                sub[
                    "inflammation"
                ]
            ),

        "inflammation_cell_counts":
            count_string(
                sub[
                    "inflammation"
                ]
            ),

        "mapping_rationale":
            unique_string(
                sub[
                    "mapping_rationale"
                ]
            ),

        "review_notes":
            unique_string(
                sub[
                    "review_notes"
                ]
            ),
    })

    source_rows.append(row)


source_audit = pd.DataFrame(
    source_rows
)


source_audit[
    "max_donor_fraction"
] = (
    source_audit[
        "max_donor_fraction"
    ]
    .round(4)
)

source_audit[
    "donor_HHI"
] = (
    source_audit[
        "donor_HHI"
    ]
    .round(4)
)


source_audit = (
    source_audit
    .sort_values(
        [
            "coarse_identity",
            "dataset",
            "primary_reference_action",
            "n_cells",
        ],
        ascending=[
            True,
            True,
            True,
            False,
        ],
    )
)


source_audit.to_csv(
    OUT
    / "01_source_annotation_audit.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 2. DATASET × IDENTITY AUDIT
# ============================================================

identity_rows = []


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

    max_frac, hhi = (
        donor_balance(sub)
    )

    identity_rows.append({
        "dataset":
            dataset,

        "coarse_identity":
            identity,

        "primary_reference_action":
            action,

        "current_reference20_member":
            identity
            in CURRENT_REFERENCE20,

        "current_ecotyper_core12_member":
            identity
            in CURRENT_ECOTYPER_CORE12,

        "anatomy_encoded_identity":
            identity
            in ANATOMY_ENCODED_IDENTITIES,

        "n_source_annotations":
            sub[
                "author_annotation_path"
            ].nunique(),

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

        "max_donor_fraction":
            max_frac,

        "donor_HHI":
            hhi,

        "disease_levels":
            unique_string(
                sub["disease"]
            ),

        "site_levels":
            unique_string(
                sub["site"]
            ),

        "inflammation_levels":
            unique_string(
                sub["inflammation"]
            ),

        "source_annotations":
            unique_string(
                sub[
                    "author_fine_annotation"
                ]
            ),
    })


dataset_identity = pd.DataFrame(
    identity_rows
)


dataset_identity[
    "max_donor_fraction"
] = (
    dataset_identity[
        "max_donor_fraction"
    ]
    .round(4)
)

dataset_identity[
    "donor_HHI"
] = (
    dataset_identity[
        "donor_HHI"
    ]
    .round(4)
)


dataset_identity = (
    dataset_identity
    .sort_values(
        [
            "coarse_identity",
            "dataset",
            "primary_reference_action",
        ]
    )
)


dataset_identity.to_csv(
    OUT
    / "02_dataset_identity_audit.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 3. CROSS-DATASET IDENTITY AUDIT
#
# One row per biological coarse identity.
# This is the main table for taxonomy v2 review.
# ============================================================

cross_rows = []


for identity, sub in cells.groupby(
    "coarse_identity",
    observed=True,
    dropna=False,
):

    datasets_present = sorted(
        sub["dataset"].unique()
    )

    row = {
        "coarse_identity":
            identity,

        "current_reference20_member":
            identity
            in CURRENT_REFERENCE20,

        "current_ecotyper_core12_member":
            identity
            in CURRENT_ECOTYPER_CORE12,

        "anatomy_encoded_identity":
            identity
            in ANATOMY_ENCODED_IDENTITIES,

        "n_datasets_present":
            len(
                datasets_present
            ),

        "datasets_present":
            " | ".join(
                datasets_present
            ),

        "n_cells_total":
            len(sub),

        "n_samples_total":
            sub[
                "global_sample_id"
            ].nunique(),

        "n_donors_total":
            sub[
                "global_donor_id"
            ].nunique(),

        "n_source_annotations_total":
            (
                sub[
                    [
                        "dataset",
                        "author_annotation_path",
                    ]
                ]
                .drop_duplicates()
                .shape[0]
            ),

        "disease_levels":
            unique_string(
                sub["disease"]
            ),

        "site_levels":
            unique_string(
                sub["site"]
            ),

        "inflammation_levels":
            unique_string(
                sub[
                    "inflammation"
                ]
            ),
    }


    # Per-dataset coverage
    for dataset in [
        "GSE282122",
        "SCP1884",
        "SCP259",
    ]:

        temp = sub.loc[
            sub["dataset"]
            == dataset
        ]

        prefix = dataset

        row[
            f"{prefix}__present"
        ] = len(temp) > 0

        row[
            f"{prefix}__n_cells"
        ] = len(temp)

        row[
            f"{prefix}__n_samples"
        ] = (
            temp[
                "global_sample_id"
            ].nunique()
        )

        row[
            f"{prefix}__n_donors"
        ] = (
            temp[
                "global_donor_id"
            ].nunique()
        )

        if len(temp) > 0:

            row[
                f"{prefix}__actions"
            ] = unique_string(
                temp[
                    "primary_reference_action"
                ]
            )

            row[
                f"{prefix}__source_annotations"
            ] = unique_string(
                temp[
                    "author_fine_annotation"
                ]
            )

        else:

            row[
                f"{prefix}__actions"
            ] = "ABSENT"

            row[
                f"{prefix}__source_annotations"
            ] = ""


    cross_rows.append(row)


cross = pd.DataFrame(
    cross_rows
)


cross = (
    cross
    .sort_values(
        [
            "current_reference20_member",
            "n_datasets_present",
            "n_cells_total",
        ],
        ascending=[
            False,
            False,
            False,
        ],
    )
)


cross.to_csv(
    OUT
    / "03_cross_dataset_identity_audit.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 4. REVIEW FOCUS TABLE
#
# No biological KEEP/DROP decision is made here.
# These are only objective reasons to inspect an identity.
# ============================================================

review = cross.copy()


review_reasons = []


for _, row in review.iterrows():

    reasons = []

    if not row[
        "current_reference20_member"
    ]:
        reasons.append(
            "NOT_IN_CURRENT_REFERENCE20"
        )

    if (
        row[
            "n_datasets_present"
        ]
        < 2
    ):
        reasons.append(
            "LIMITED_CROSS_DATASET_COVERAGE"
        )

    if row[
        "anatomy_encoded_identity"
    ]:
        reasons.append(
            "ANATOMY_ENCODED_IDENTITY"
        )

    # Any upstream non-KEEP mapping
    action_cols = [
        "GSE282122__actions",
        "SCP1884__actions",
        "SCP259__actions",
    ]

    actions = " | ".join(
        str(
            row[x]
        )
        for x in action_cols
    )

    if (
        "EXCLUDE" in actions
        or
        "HOLD" in actions
        or
        "REVIEW" in actions
    ):
        reasons.append(
            "UPSTREAM_NON_KEEP_MAPPING"
        )

    review_reasons.append(
        " | ".join(reasons)
    )


review[
    "audit_review_reason"
] = review_reasons


review = (
    review.loc[
        review[
            "audit_review_reason"
        ] != ""
    ]
    .copy()
)


review.to_csv(
    OUT
    / "04_review_focus.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 5. Action summary
# ============================================================

action_summary = (
    cells.groupby(
        [
            "dataset",
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
        n_identities=(
            "coarse_identity",
            "nunique",
        ),
    )
    .reset_index()
)


action_summary.to_csv(
    OUT
    / "05_mapping_action_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 6. Summary
# ============================================================

summary_lines = [
    "IBD signature taxonomy v2 audit",
    "================================",
    "",
    f"Total candidate cells: {len(cells):,}",
    "",
]

for dataset in [
    "GSE282122",
    "SCP1884",
    "SCP259",
]:

    temp = cells.loc[
        cells["dataset"]
        == dataset
    ]

    summary_lines.extend([
        dataset,
        "-" * len(dataset),

        f"Cells: {len(temp):,}",

        "Donors: "
        f"{temp['global_donor_id'].nunique():,}",

        "Samples: "
        f"{temp['global_sample_id'].nunique():,}",

        "Coarse identities: "
        f"{temp['coarse_identity'].nunique():,}",

        "Author fine annotations: "
        f"{temp['author_fine_annotation'].nunique():,}",

        "",
    ])


summary_lines.extend([
    "Combined:",
    f"Coarse identities: "
    f"{cells['coarse_identity'].nunique():,}",

    f"Current reference20 members represented: "
    f"{len(set(cells['coarse_identity']) & CURRENT_REFERENCE20)}",

    f"Current EcoTyper CORE12 members represented: "
    f"{len(set(cells['coarse_identity']) & CURRENT_ECOTYPER_CORE12)}",

    "",
    "Outputs:",
    "01_source_annotation_audit.tsv",
    "02_dataset_identity_audit.tsv",
    "03_cross_dataset_identity_audit.tsv",
    "04_review_focus.tsv",
    "05_mapping_action_summary.tsv",
    "",
    "FINAL STATUS: PASS",
])


(
    OUT
    / "taxonomy_v2_audit_summary.txt"
).write_text(
    "\n".join(
        summary_lines
    )
    + "\n",
    encoding="utf-8",
)


print()
print(
    "=" * 78
)
print(
    "AUDIT COMPLETE"
)
print(
    "=" * 78
)

print(
    f"Total candidate cells: "
    f"{len(cells):,}"
)

print(
    f"Coarse identities: "
    f"{cells['coarse_identity'].nunique()}"
)

print(
    f"Source-annotation audit rows: "
    f"{len(source_audit):,}"
)

print(
    f"Cross-dataset identities: "
    f"{len(cross):,}"
)

print()
print(
    "Output directory:"
)
print(
    OUT
)

print()
print(
    "FINAL STATUS: PASS"
)
