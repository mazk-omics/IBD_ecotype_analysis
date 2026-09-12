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

INPUT_PARQUET = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "01_candidate_pool"
    / "candidate_cells.parquet"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "02_taxonomy_mapping"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "taxonomy_mapping_audit_summary.txt"
)


# ============================================================
# Author annotation hierarchy
# ============================================================

PATH_COLS = [
    "bucket",
    "sub_bucket",
    "major",
    "minor",
    "final_analysis",
]


# ============================================================
# Primary candidate taxonomy
#
# These are allowed downstream mapping targets.
# No automatic mapping is performed in Phase 2A.
# ============================================================

ALLOWED_IDENTITIES = [
    ("B", "B lymphocytes"),
    ("Plasma", "Plasma cells / plasmablasts"),
    ("CD4_T", "Conventional CD4 T cells"),
    ("CD8_T", "Conventional CD8 T cells"),
    ("NK", "Natural killer cells"),
    (
        "Monocyte_Macrophage",
        "Monocytes and macrophages",
    ),
    ("DC", "Dendritic cells"),
    ("Mast", "Mast cells"),
    (
        "Fibroblast",
        "Fibroblast lineage",
    ),
    ("Pericyte", "Pericytes"),
    (
        "Endothelial",
        "Blood and lymphatic endothelial cells",
    ),
    (
        "Ileal_absorptive",
        "Ileal absorptive epithelial cells",
    ),
    (
        "Colonic_absorptive",
        "Non-ileal / colonic absorptive epithelial cells",
    ),
    ("Goblet", "Goblet cells"),
    (
        "Stem_TA_progenitor",
        "Stem and transit-amplifying progenitor cells",
    ),
    ("Paneth", "Paneth cells"),
    (
        "Enteroendocrine",
        "Enteroendocrine cells",
    ),
    (
        "ILC",
        "Innate lymphoid cells excluding NK",
    ),
    (
        "EXCLUDE",
        "Explicitly excluded from primary reference",
    ),
    (
        "REVIEW",
        "Requires manual biological review",
    ),
    (
        "UNMAPPED",
        "Not yet assigned",
    ),
]


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
    summary.append(
        f"[FAIL] {message}"
    )
    write_summary(summary)
    raise RuntimeError(message)


# ============================================================
# 0. Environment preflight
# ============================================================

summary = []

summary.append(
    "GSE282122 Reference Audit v2"
)
summary.append(
    "Phase 2A: Annotation-path audit and taxonomy mapping template"
)
summary.append("=" * 76)
summary.append("")

python_exe = Path(
    sys.executable
).resolve()

summary.append("Environment")
summary.append("-" * 76)

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

summary.append(
    f"pyarrow: {package_version('pyarrow')}"
)

summary.append("")


# ============================================================
# 1. Load candidate pool
# ============================================================

summary.append("Input")
summary.append("-" * 76)

summary.append(
    f"Candidate metadata: {INPUT_PARQUET}"
)

if not INPUT_PARQUET.is_file():
    fail(
        f"Candidate pool not found: {INPUT_PARQUET}",
        summary,
    )

df = pd.read_parquet(
    INPUT_PARQUET
)

summary.append(
    f"Rows: {len(df):,}"
)

summary.append(
    f"Columns: {len(df.columns):,}"
)

summary.append("")


# ============================================================
# 2. Input validation
# ============================================================

required = [
    "cell_id",
    "sample_id",
    "Patient",
    "Disease",
    "Site",
    "Inflammation",
] + PATH_COLS

missing = [
    x for x in required
    if x not in df.columns
]

if missing:
    fail(
        "Required columns missing: "
        + ", ".join(missing),
        summary,
    )

if df["cell_id"].isna().any():
    fail(
        "cell_id contains missing values.",
        summary,
    )

if not df["cell_id"].is_unique:
    fail(
        "cell_id is not unique.",
        summary,
    )

for col in PATH_COLS:

    if df[col].isna().any():
        fail(
            f"{col} contains missing annotations.",
            summary,
        )

summary.append("Input validation")
summary.append("-" * 76)

summary.append(
    "[PASS] Required metadata fields present."
)

summary.append(
    "[PASS] cell_id is complete and unique."
)

summary.append(
    "[PASS] Author annotation hierarchy is complete."
)

summary.append("")


# ============================================================
# 3. Basic candidate-pool checkpoint
# ============================================================

n_cells = len(df)
n_samples = df["sample_id"].nunique()
n_donors = df["Patient"].nunique()

summary.append("Candidate pool")
summary.append("-" * 76)

summary.append(
    f"Cells: {n_cells:,}"
)

summary.append(
    f"Samples: {n_samples:,}"
)

summary.append(
    f"Donors: {n_donors:,}"
)

if n_cells != 505_591:
    fail(
        "Candidate pool no longer contains 505,591 cells.",
        summary,
    )

if n_samples != 116:
    fail(
        "Candidate pool no longer contains 116 samples.",
        summary,
    )

if n_donors != 41:
    fail(
        "Candidate pool no longer contains 41 donors.",
        summary,
    )

summary.append(
    "[PASS] Phase-1 candidate pool checkpoints retained."
)

summary.append("")


# ============================================================
# 4. Annotation path statistics
# ============================================================

path_stats = (
    df
    .groupby(
        PATH_COLS,
        observed=True,
        dropna=False,
    )
    .agg(
        n_cells=("cell_id", "size"),
        n_samples=("sample_id", "nunique"),
        n_donors=("Patient", "nunique"),
    )
    .reset_index()
)


# ============================================================
# 5. Disease composition
# ============================================================

disease_counts = (
    df
    .groupby(
        PATH_COLS + ["Disease"],
        observed=True,
        dropna=False,
    )
    .size()
    .unstack(
        fill_value=0
    )
    .reset_index()
)

disease_counts.columns.name = None

for disease in [
    "CD",
    "UC",
    "Healthy",
]:

    if disease not in disease_counts.columns:
        disease_counts[disease] = 0

disease_counts = disease_counts.rename(
    columns={
        "CD": "n_CD",
        "UC": "n_UC",
        "Healthy": "n_Healthy",
    }
)

path_stats = path_stats.merge(
    disease_counts[
        PATH_COLS
        + [
            "n_CD",
            "n_UC",
            "n_Healthy",
        ]
    ],
    on=PATH_COLS,
    how="left",
    validate="one_to_one",
)


# ============================================================
# 6. Inflammation composition
# ============================================================

inflammation_counts = (
    df
    .groupby(
        PATH_COLS + ["Inflammation"],
        observed=True,
        dropna=False,
    )
    .size()
    .unstack(
        fill_value=0
    )
    .reset_index()
)

inflammation_counts.columns.name = None

for state in [
    "Inflamed",
    "Non_Inflamed",
    "Healthy",
]:

    if state not in inflammation_counts.columns:
        inflammation_counts[state] = 0

inflammation_counts = (
    inflammation_counts.rename(
        columns={
            "Inflamed":
                "n_Inflamed",
            "Non_Inflamed":
                "n_NonInflamed",
            "Healthy":
                "n_HealthyInflammation",
        }
    )
)

path_stats = path_stats.merge(
    inflammation_counts[
        PATH_COLS
        + [
            "n_Inflamed",
            "n_NonInflamed",
            "n_HealthyInflammation",
        ]
    ],
    on=PATH_COLS,
    how="left",
    validate="one_to_one",
)


# ============================================================
# 7. Donor dominance and HHI
#
# max_donor_fraction:
#   largest fraction of cells contributed by a single donor
#
# donor_HHI:
#   sum(p_i^2)
#   lower = more evenly distributed across donors
# ============================================================

donor_counts = (
    df
    .groupby(
        PATH_COLS + ["Patient"],
        observed=True,
        dropna=False,
    )
    .size()
    .rename("n_cells_donor")
    .reset_index()
)

donor_counts[
    "path_total"
] = (
    donor_counts
    .groupby(
        PATH_COLS,
        observed=True,
        dropna=False,
    )["n_cells_donor"]
    .transform("sum")
)

donor_counts[
    "donor_fraction"
] = (
    donor_counts["n_cells_donor"]
    /
    donor_counts["path_total"]
)

donor_counts[
    "donor_fraction_sq"
] = (
    donor_counts["donor_fraction"] ** 2
)

donor_balance = (
    donor_counts
    .groupby(
        PATH_COLS,
        observed=True,
        dropna=False,
    )
    .agg(
        max_donor_cells=(
            "n_cells_donor",
            "max",
        ),
        max_donor_fraction=(
            "donor_fraction",
            "max",
        ),
        donor_HHI=(
            "donor_fraction_sq",
            "sum",
        ),
    )
    .reset_index()
)

path_stats = path_stats.merge(
    donor_balance,
    on=PATH_COLS,
    how="left",
    validate="one_to_one",
)

path_stats[
    "max_donor_fraction"
] = (
    path_stats[
        "max_donor_fraction"
    ]
    .round(4)
)

path_stats[
    "donor_HHI"
] = (
    path_stats[
        "donor_HHI"
    ]
    .round(4)
)


# ============================================================
# 8. Site composition
# ============================================================

site_counts = (
    df
    .groupby(
        PATH_COLS + ["Site"],
        observed=True,
        dropna=False,
    )
    .size()
    .unstack(
        fill_value=0
    )
    .reset_index()
)

site_counts.columns.name = None

sites = [
    "Terminal_Ileum",
    "Ascending_Colon",
    "Descending_Colon",
    "Sigmoid",
    "Rectum",
]

for site in sites:

    if site not in site_counts.columns:
        site_counts[site] = 0

site_counts = site_counts.rename(
    columns={
        site: f"n_{site}"
        for site in sites
    }
)

path_stats = path_stats.merge(
    site_counts[
        PATH_COLS
        + [
            f"n_{site}"
            for site in sites
        ]
    ],
    on=PATH_COLS,
    how="left",
    validate="one_to_one",
)


# ============================================================
# 9. Stable human-readable path key
# ============================================================

path_stats[
    "path_key"
] = (
    path_stats[
        PATH_COLS
    ]
    .astype(str)
    .agg(
        " > ".join,
        axis=1,
    )
)


# ============================================================
# 10. Detect final_analysis labels occurring in >1 path
# ============================================================

final_path_count = (
    path_stats
    .groupby(
        "final_analysis",
        observed=True,
    )
    .size()
    .rename(
        "n_hierarchy_paths"
    )
    .reset_index()
)

reused_final = (
    final_path_count.loc[
        final_path_count[
            "n_hierarchy_paths"
        ] > 1
    ]
    .copy()
)

reused_labels = set(
    reused_final[
        "final_analysis"
    ].astype(str)
)

reused_paths = (
    path_stats.loc[
        path_stats[
            "final_analysis"
        ]
        .astype(str)
        .isin(reused_labels)
    ]
    .copy()
)

reused_paths = reused_paths.merge(
    reused_final,
    on="final_analysis",
    how="left",
    validate="many_to_one",
)

reused_paths.to_csv(
    OUTPUT_DIR
    / "final_analysis_reused_paths.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 11. final_analysis-level summary
# ============================================================

final_summary = (
    path_stats
    .groupby(
        "final_analysis",
        observed=True,
    )
    .agg(
        n_hierarchy_paths=(
            "path_key",
            "nunique",
        ),
        n_cells=(
            "n_cells",
            "sum",
        ),
        max_path_donors=(
            "n_donors",
            "max",
        ),
        max_path_samples=(
            "n_samples",
            "max",
        ),
    )
    .reset_index()
)

final_summary[
    "reused_across_paths"
] = (
    final_summary[
        "n_hierarchy_paths"
    ] > 1
)

final_summary = (
    final_summary
    .sort_values(
        [
            "reused_across_paths",
            "n_cells",
        ],
        ascending=[
            False,
            False,
        ],
    )
)

final_summary.to_csv(
    OUTPUT_DIR
    / "final_analysis_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 12. Write complete path audit
# ============================================================

path_stats = (
    path_stats
    .sort_values(
        [
            "bucket",
            "sub_bucket",
            "major",
            "minor",
            "final_analysis",
        ]
    )
    .reset_index(
        drop=True
    )
)

path_stats.insert(
    0,
    "path_id",
    [
        f"GSE282122_PATH_{i:03d}"
        for i in range(
            1,
            len(path_stats) + 1,
        )
    ],
)

path_stats.to_csv(
    OUTPUT_DIR
    / "annotation_path_audit.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 13. Create manual taxonomy mapping template
#
# IMPORTANT:
# No target identity is assigned automatically.
# ============================================================

mapping_template = (
    path_stats.copy()
)

mapping_template[
    "target_identity"
] = ""

mapping_template[
    "mapping_status"
] = "UNREVIEWED"

mapping_template[
    "mapping_rationale"
] = ""

mapping_template[
    "exclude_reason"
] = ""

mapping_template[
    "review_notes"
] = ""

mapping_cols = [
    "path_id",
    "bucket",
    "sub_bucket",
    "major",
    "minor",
    "final_analysis",
    "n_cells",
    "n_samples",
    "n_donors",
    "n_CD",
    "n_UC",
    "n_Healthy",
    "n_Inflamed",
    "n_NonInflamed",
    "n_HealthyInflammation",
    "max_donor_fraction",
    "donor_HHI",
    "target_identity",
    "mapping_status",
    "mapping_rationale",
    "exclude_reason",
    "review_notes",
    "path_key",
]

mapping_template[
    mapping_cols
].to_csv(
    OUTPUT_DIR
    / "taxonomy_mapping_template.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 14. Allowed taxonomy values
# ============================================================

allowed_df = pd.DataFrame(
    ALLOWED_IDENTITIES,
    columns=[
        "target_identity",
        "description",
    ],
)

allowed_df.to_csv(
    OUTPUT_DIR
    / "taxonomy_allowed_identities.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 15. Additional coarse annotation counts
# ============================================================

for level in PATH_COLS:

    level_summary = (
        df
        .groupby(
            level,
            observed=True,
            dropna=False,
        )
        .agg(
            n_cells=("cell_id", "size"),
            n_samples=("sample_id", "nunique"),
            n_donors=("Patient", "nunique"),
        )
        .reset_index()
        .sort_values(
            "n_cells",
            ascending=False,
        )
    )

    level_summary.to_csv(
        OUTPUT_DIR
        / f"annotation_level__{level}.tsv",
        sep="\t",
        index=False,
    )


# ============================================================
# 16. Summary
# ============================================================

n_paths = len(path_stats)

n_final = (
    path_stats[
        "final_analysis"
    ]
    .nunique()
)

n_reused_labels = len(
    reused_final
)

summary.append(
    "Annotation hierarchy"
)
summary.append("-" * 76)

summary.append(
    "Hierarchy: "
    "bucket -> sub_bucket -> major -> minor -> final_analysis"
)

summary.append(
    f"Unique hierarchy paths: {n_paths}"
)

summary.append(
    f"Unique final_analysis labels: {n_final}"
)

summary.append(
    "final_analysis labels reused "
    f"across >1 path: {n_reused_labels}"
)

summary.append("")

summary.append(
    "Outputs"
)
summary.append("-" * 76)

for filename in [
    "annotation_path_audit.tsv",
    "taxonomy_mapping_template.tsv",
    "taxonomy_allowed_identities.tsv",
    "final_analysis_summary.tsv",
    "final_analysis_reused_paths.tsv",
]:

    summary.append(
        f"[PASS] {filename}"
    )

summary.append("")

summary.append(
    "Final Phase-2A status"
)
summary.append("-" * 76)

summary.append(
    "[PASS] Candidate annotation hierarchy audited."
)

summary.append(
    "[PASS] Donor/disease/inflammation/site "
    "composition calculated."
)

summary.append(
    "[PASS] Donor dominance metrics calculated."
)

summary.append(
    "[PASS] Manual taxonomy mapping template generated."
)

summary.append(
    "[PASS] No author annotation was modified."
)

summary.append(
    "[PASS] No cell was excluded or remapped."
)

write_summary(
    summary
)


# ============================================================
# Console output
# ============================================================

print()
print("=" * 76)
print(
    "GSE282122 Phase 2A taxonomy preparation completed"
)
print("=" * 76)

print()
print(
    f"Candidate cells       : {n_cells:,}"
)
print(
    f"Annotation paths      : {n_paths}"
)
print(
    f"final_analysis labels : {n_final}"
)
print(
    f"Reused final labels   : {n_reused_labels}"
)

print()
print(
    f"Output directory:"
)
print(
    f"  {OUTPUT_DIR}"
)

print()
print(
    f"Summary:"
)
print(
    f"  {SUMMARY_FILE}"
)