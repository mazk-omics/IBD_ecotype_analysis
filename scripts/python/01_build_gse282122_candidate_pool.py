#!/usr/bin/env python3

from pathlib import Path
import sys
import importlib.metadata as md

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
    / "00_schema_audit"
    / "author_obs_metadata.parquet"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "01_candidate_pool"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "candidate_pool_summary.txt"
)


# ============================================================
# Expected checkpoints
#
# These are validation targets from the independent metadata
# audit. They do NOT define the filtering result.
# ============================================================

EXPECTED_TOTAL_CELLS = 505_591
EXPECTED_TOTAL_SAMPLES = 116
EXPECTED_TOTAL_DONORS = 41

EXPECTED_BY_DISEASE = {
    "CD": {
        "cells": 206_701,
        "samples": 49,
        "donors": 16,
    },
    "UC": {
        "cells": 242_197,
        "samples": 55,
        "donors": 22,
    },
    "Healthy": {
        "cells": 56_693,
        "samples": 12,
        "donors": 3,
    },
}

SAMPLE_LEVEL_FIELDS = [
    "Patient",
    "Disease",
    "Site",
    "Treatment",
    "Disease_duration",
    "Inflammation",
    "Age",
    "Gender",
    "Ethnicity",
    "Inflammation_score",
    "Ileum_vs_Colon",
    "LibraryType",
    "CellsLoaded",
    "Match",
    "Batch",
    "Remission_status",
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
    summary.append(f"[FAIL] {message}")
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
    "Phase 1: Reference candidate pool construction"
)
summary.append("=" * 72)
summary.append("")

python_exe = Path(
    sys.executable
).resolve()

summary.append("Environment")
summary.append("-" * 72)
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
    f"pyarrow: {package_version('pyarrow')}"
)
summary.append("")


# ============================================================
# 1. Input
# ============================================================

summary.append("Input")
summary.append("-" * 72)
summary.append(
    f"Metadata: {INPUT_PARQUET}"
)

if not INPUT_PARQUET.is_file():
    fail(
        f"Input Parquet not found: {INPUT_PARQUET}",
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
# 2. Required-column validation
# ============================================================

required = [
    "obs_name",
    "sample_id",
    "Patient",
    "Disease",
    "Treatment",
    "Site",
    "Inflammation",
    "Inflammation_score",
    "Remission_status",
    "bucket",
    "sub_bucket",
    "major",
    "minor",
    "final_analysis",
]

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

summary.append("Required fields")
summary.append("-" * 72)
summary.append(
    "[PASS] All required Phase-1 fields are present."
)
summary.append("")


# ============================================================
# 3. Define derived cell_id
#
# Never modify the source H5AD. obs_name becomes the stable
# cell identifier in derived audit outputs.
# ============================================================

if df["obs_name"].isna().any():
    fail(
        "obs_name contains missing values.",
        summary,
    )

if not df["obs_name"].is_unique:
    fail(
        "obs_name is not unique.",
        summary,
    )

df = df.copy()

df.insert(
    0,
    "cell_id",
    df["obs_name"].astype("string"),
)

summary.append("Cell identifier")
summary.append("-" * 72)
summary.append(
    "[PASS] cell_id derived from unique obs_name."
)
summary.append("")


# ============================================================
# 4. Explicit semantic checks
# ============================================================

summary.append("Metadata semantic checks")
summary.append("-" * 72)

observed_diseases = set(
    df["Disease"]
    .dropna()
    .astype(str)
    .unique()
)

expected_diseases = {
    "CD",
    "UC",
    "Healthy",
}

if observed_diseases != expected_diseases:
    fail(
        "Unexpected Disease values: "
        f"{sorted(observed_diseases)}",
        summary,
    )

summary.append(
    "[PASS] Disease = CD / UC / Healthy."
)

# Healthy treatment should be missing in author metadata.
healthy = (
    df["Disease"].astype(str)
    == "Healthy"
)

ibd = (
    df["Disease"].astype(str)
    .isin(["CD", "UC"])
)

if not df.loc[
    healthy,
    "Treatment"
].isna().all():
    fail(
        "At least one Healthy cell has non-missing Treatment.",
        summary,
    )

summary.append(
    "[PASS] All Healthy cells have missing Treatment."
)

if df.loc[
    ibd,
    "Treatment"
].isna().any():
    fail(
        "At least one CD/UC cell has missing Treatment.",
        summary,
    )

summary.append(
    "[PASS] No CD/UC cell has missing Treatment."
)

ibd_treatment_values = set(
    df.loc[
        ibd,
        "Treatment",
    ]
    .astype(str)
    .unique()
)

if ibd_treatment_values != {
    "Pre",
    "Post",
}:
    fail(
        "Unexpected IBD Treatment values: "
        f"{sorted(ibd_treatment_values)}",
        summary,
    )

summary.append(
    "[PASS] IBD Treatment = Pre / Post."
)

expected_inflammation = {
    "Inflamed",
    "Non_Inflamed",
    "Healthy",
}

observed_inflammation = set(
    df["Inflammation"]
    .dropna()
    .astype(str)
    .unique()
)

if observed_inflammation != expected_inflammation:
    fail(
        "Unexpected Inflammation values: "
        f"{sorted(observed_inflammation)}",
        summary,
    )

summary.append(
    "[PASS] Inflammation categories validated."
)
summary.append("")


# ============================================================
# 5. Sample-level consistency re-check
# ============================================================

summary.append(
    "Sample-level consistency"
)
summary.append("-" * 72)

available_sample_fields = [
    x
    for x in SAMPLE_LEVEL_FIELDS
    if x in df.columns
]

grouped = df.groupby(
    "sample_id",
    observed=True,
    dropna=False,
)

for field in available_sample_fields:

    nunique = grouped[field].nunique(
        dropna=False
    )

    bad = nunique[
        nunique > 1
    ]

    if len(bad) > 0:

        bad.rename(
            "n_unique_values"
        ).reset_index().to_csv(
            OUTPUT_DIR
            / f"inconsistent_samples__{field}.tsv",
            sep="\t",
            index=False,
        )

        fail(
            f"{field}: "
            f"{len(bad)} samples contain "
            "multiple values.",
            summary,
        )

summary.append(
    "[PASS] All tested sample-level fields "
    "are internally consistent."
)
summary.append("")


# ============================================================
# 6. Build all-sample manifest
# ============================================================

sample_metadata = (
    df[
        ["sample_id"]
        + available_sample_fields
    ]
    .drop_duplicates()
)

duplicate_sample_ids = (
    sample_metadata["sample_id"]
    .duplicated()
)

if duplicate_sample_ids.any():
    fail(
        "Sample manifest still contains "
        "duplicate sample_id values.",
        summary,
    )

cell_counts_all = (
    df.groupby(
        "sample_id",
        observed=True,
    )
    .size()
    .rename("n_cells_all")
    .reset_index()
)

sample_manifest = (
    sample_metadata
    .merge(
        cell_counts_all,
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
)

sample_manifest.to_csv(
    OUTPUT_DIR
    / "sample_manifest_all.csv",
    index=False,
)


# ============================================================
# 7. Define reference eligibility
#
# Primary rule:
# Healthy
# OR
# CD/UC at Pre-treatment
# ============================================================

eligible_mask = (
    healthy
    |
    (
        ibd
        &
        (
            df["Treatment"]
            .astype("string")
            == "Pre"
        )
    )
)

df["reference_eligible"] = (
    eligible_mask
)

df["eligibility_reason"] = "Excluded"

df.loc[
    healthy,
    "eligibility_reason"
] = "Healthy_control"

df.loc[
    ibd
    &
    (
        df["Treatment"]
        .astype("string")
        == "Pre"
    ),
    "eligibility_reason"
] = "IBD_pre_treatment"

df.loc[
    ibd
    &
    (
        df["Treatment"]
        .astype("string")
        == "Post"
    ),
    "eligibility_reason"
] = "IBD_post_treatment_excluded"


# ============================================================
# 8. Candidate cells
# ============================================================

candidate = (
    df.loc[
        df["reference_eligible"]
    ]
    .copy()
)

candidate.to_parquet(
    OUTPUT_DIR
    / "candidate_cells.parquet",
    index=False,
    engine="pyarrow",
)


# ============================================================
# 9. Eligible sample manifest
# ============================================================

eligible_sample_ids = set(
    candidate["sample_id"]
    .astype(str)
    .unique()
)

sample_manifest[
    "reference_eligible"
] = (
    sample_manifest[
        "sample_id"
    ]
    .astype(str)
    .isin(
        eligible_sample_ids
    )
)

sample_manifest[
    "eligibility_reason"
] = "Excluded"

healthy_sample = (
    sample_manifest[
        "Disease"
    ]
    .astype(str)
    == "Healthy"
)

pre_ibd_sample = (
    sample_manifest[
        "Disease"
    ]
    .astype(str)
    .isin(["CD", "UC"])
    &
    (
        sample_manifest[
            "Treatment"
        ]
        .astype("string")
        == "Pre"
    )
)

post_ibd_sample = (
    sample_manifest[
        "Disease"
    ]
    .astype(str)
    .isin(["CD", "UC"])
    &
    (
        sample_manifest[
            "Treatment"
        ]
        .astype("string")
        == "Post"
    )
)

sample_manifest.loc[
    healthy_sample,
    "eligibility_reason",
] = "Healthy_control"

sample_manifest.loc[
    pre_ibd_sample,
    "eligibility_reason",
] = "IBD_pre_treatment"

sample_manifest.loc[
    post_ibd_sample,
    "eligibility_reason",
] = "IBD_post_treatment_excluded"

candidate_counts = (
    candidate.groupby(
        "sample_id",
        observed=True,
    )
    .size()
    .rename("n_candidate_cells")
    .reset_index()
)

sample_manifest = (
    sample_manifest
    .merge(
        candidate_counts,
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
)

sample_manifest[
    "n_candidate_cells"
] = (
    sample_manifest[
        "n_candidate_cells"
    ]
    .fillna(0)
    .astype(int)
)

eligible_manifest = (
    sample_manifest.loc[
        sample_manifest[
            "reference_eligible"
        ]
    ]
    .copy()
)

sample_manifest.to_csv(
    OUTPUT_DIR
    / "sample_manifest_all.csv",
    index=False,
)

eligible_manifest.to_csv(
    OUTPUT_DIR
    / "sample_manifest_eligible.csv",
    index=False,
)


# ============================================================
# 10. Summary tables
# ============================================================

by_disease = (
    candidate
    .groupby(
        "Disease",
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

by_disease.to_csv(
    OUTPUT_DIR
    / "candidate_summary_by_disease.tsv",
    sep="\t",
    index=False,
)

by_site = (
    candidate
    .groupby(
        "Site",
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

by_site.to_csv(
    OUTPUT_DIR
    / "candidate_summary_by_site.tsv",
    sep="\t",
    index=False,
)

by_inflammation = (
    candidate
    .groupby(
        "Inflammation",
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

by_inflammation.to_csv(
    OUTPUT_DIR
    / "candidate_summary_by_inflammation.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 11. Hard validation checkpoints
# ============================================================

summary.append(
    "Candidate-pool checkpoints"
)
summary.append("-" * 72)

observed_total_cells = len(
    candidate
)

observed_total_samples = (
    candidate[
        "sample_id"
    ]
    .nunique()
)

observed_total_donors = (
    candidate[
        "Patient"
    ]
    .nunique()
)

summary.append(
    f"Candidate cells: {observed_total_cells:,}"
)
summary.append(
    f"Eligible samples: {observed_total_samples:,}"
)
summary.append(
    f"Eligible donors: {observed_total_donors:,}"
)

if observed_total_cells != EXPECTED_TOTAL_CELLS:
    fail(
        "Candidate cell count mismatch: "
        f"expected {EXPECTED_TOTAL_CELLS:,}, "
        f"observed {observed_total_cells:,}.",
        summary,
    )

if observed_total_samples != EXPECTED_TOTAL_SAMPLES:
    fail(
        "Eligible sample count mismatch: "
        f"expected {EXPECTED_TOTAL_SAMPLES:,}, "
        f"observed {observed_total_samples:,}.",
        summary,
    )

if observed_total_donors != EXPECTED_TOTAL_DONORS:
    fail(
        "Eligible donor count mismatch: "
        f"expected {EXPECTED_TOTAL_DONORS:,}, "
        f"observed {observed_total_donors:,}.",
        summary,
    )

summary.append(
    "[PASS] Global checkpoints validated."
)
summary.append("")


# ============================================================
# 12. Disease-specific validation
# ============================================================

summary.append(
    "Disease-specific checkpoints"
)
summary.append("-" * 72)

for disease, expected in EXPECTED_BY_DISEASE.items():

    sub = candidate.loc[
        candidate["Disease"]
        .astype(str)
        == disease
    ]

    obs_cells = len(sub)

    obs_samples = (
        sub["sample_id"]
        .nunique()
    )

    obs_donors = (
        sub["Patient"]
        .nunique()
    )

    summary.append(
        f"{disease}: "
        f"{obs_cells:,} cells, "
        f"{obs_samples} samples, "
        f"{obs_donors} donors"
    )

    if obs_cells != expected["cells"]:
        fail(
            f"{disease} cell count mismatch.",
            summary,
        )

    if obs_samples != expected["samples"]:
        fail(
            f"{disease} sample count mismatch.",
            summary,
        )

    if obs_donors != expected["donors"]:
        fail(
            f"{disease} donor count mismatch.",
            summary,
        )

summary.append(
    "[PASS] Disease-specific checkpoints validated."
)
summary.append("")


# ============================================================
# 13. Final status
# ============================================================

summary.append(
    "Final Phase-1 status"
)
summary.append("-" * 72)

summary.append(
    "[PASS] Candidate pool constructed."
)
summary.append(
    "[PASS] Healthy controls retained."
)
summary.append(
    "[PASS] Pre-treatment CD/UC retained."
)
summary.append(
    "[PASS] Post-treatment CD/UC excluded."
)
summary.append(
    "[PASS] Original author metadata preserved."
)
summary.append(
    "[PASS] cell_id derived only in downstream metadata."
)

write_summary(
    summary
)


# ============================================================
# Console output
# ============================================================

print()
print("=" * 72)
print(
    "GSE282122 Phase 1 candidate pool completed"
)
print("=" * 72)

print()
print(
    f"Candidate cells : {observed_total_cells:,}"
)
print(
    f"Eligible samples: {observed_total_samples:,}"
)
print(
    f"Eligible donors : {observed_total_donors:,}"
)

print()
print(by_disease.to_string(index=False))

print()
print(
    f"Outputs: {OUTPUT_DIR}"
)
print(
    f"Summary: {SUMMARY_FILE}"
)