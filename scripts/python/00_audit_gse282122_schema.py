#!/usr/bin/env python3

from pathlib import Path
import sys
import platform
import importlib.metadata as md

import numpy as np
import pandas as pd
import anndata as ad
import pyarrow
import h5py


# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT = Path("/home/mazekai/IBD_EcoTyper")

# Confirmed by smoke test
H5AD = Path(
    "/home/mazekai/enteric_glia/data_scRNA/GSE282122/"
    "TAURUS_raw_counts_annotated_final.h5ad"
)

EXPECTED_ENV_ROOT = Path(
    "/home/mazekai/miniconda3/envs/ibd_refaudit"
)

AUDIT_ROOT = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
)

OUTPUT_DIR = (
    AUDIT_ROOT
    / "output"
    / "00_schema_audit"
)

LOG_DIR = AUDIT_ROOT / "logs"

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

LOG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "schema_audit_summary.txt"
)


# ============================================================
# Known dataset expectations
#
# These are used only as audit references.
# They are NOT used to filter the object.
# ============================================================

EXPECTED_N_OBS = 987_743

EXPECTED_FIELDS = [
    "cell_id",
    "sample_id",
    "Patient",
    "Disease",
    "Site",
    "Inflammation",
    "Inflammation_score",
    "Treatment",
    "Remission_status",
    "final_analysis",
    "minor",
    "major",
    "sub_bucket",
    "bucket",
]

ANNOTATION_HIERARCHY = [
    "bucket",
    "sub_bucket",
    "major",
    "minor",
    "final_analysis",
]

SAMPLE_LEVEL_FIELDS = [
    "Patient",
    "Disease",
    "Site",
    "Inflammation",
    "Inflammation_score",
    "Treatment",
    "Remission_status",
]

KEY_FIELDS = [
    "Disease",
    "Site",
    "Inflammation",
    "Treatment",
    "Remission_status",
    "bucket",
    "sub_bucket",
    "major",
    "minor",
    "final_analysis",
]


# ============================================================
# Helper functions
# ============================================================

def package_version(name):
    try:
        return md.version(name)
    except md.PackageNotFoundError:
        return "NOT INSTALLED"


def safe_nunique(series):
    return int(
        series.nunique(
            dropna=True
        )
    )


def example_values(series, n=5):
    values = (
        series
        .dropna()
        .astype("string")
        .drop_duplicates()
        .head(n)
        .tolist()
    )

    return " | ".join(
        map(str, values)
    )


def write_summary(lines):
    SUMMARY_FILE.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


# ============================================================
# 0. Environment preflight
# ============================================================

summary = []

summary.append(
    "GSE282122 Reference Audit v2"
)
summary.append(
    "Phase 0: H5AD schema and metadata audit"
)
summary.append("=" * 72)
summary.append("")

summary.append(
    "Environment preflight"
)
summary.append("-" * 72)

python_exe = Path(
    sys.executable
).resolve()

summary.append(
    f"Python executable: {python_exe}"
)

summary.append(
    f"Required environment: {EXPECTED_ENV_ROOT}"
)

if EXPECTED_ENV_ROOT not in python_exe.parents:

    summary.append(
        "[FAIL] Script is not running "
        "inside ibd_refaudit."
    )

    write_summary(summary)

    raise RuntimeError(
        "\nIncorrect Python environment.\n\n"
        "Run this script with:\n"
        "/home/mazekai/miniconda3/envs/"
        "ibd_refaudit/bin/python\n"
    )

summary.append(
    "[PASS] Correct Conda environment detected."
)
summary.append("")

versions = {
    "Python": sys.version.split()[0],
    "numpy": package_version("numpy"),
    "pandas": package_version("pandas"),
    "pyarrow": package_version("pyarrow"),
    "anndata": package_version("anndata"),
    "h5py": package_version("h5py"),
    "scipy": package_version("scipy"),
}

for name, version in versions.items():
    summary.append(
        f"{name:10s}: {version}"
    )

summary.append("")
summary.append(
    f"Platform: {platform.platform()}"
)
summary.append("")


# ============================================================
# 1. Input validation
# ============================================================

summary.append(
    "Input validation"
)
summary.append("-" * 72)

summary.append(
    f"H5AD: {H5AD}"
)

if not H5AD.is_file():

    summary.append(
        "[FAIL] H5AD file does not exist."
    )

    write_summary(summary)

    raise FileNotFoundError(
        H5AD
    )

size_gib = (
    H5AD.stat().st_size
    / (1024 ** 3)
)

summary.append(
    f"File size: {size_gib:.2f} GiB"
)

summary.append(
    "[PASS] Input file exists."
)
summary.append("")


# ============================================================
# 2. Open AnnData in backed READ-ONLY mode
# ============================================================

print(
    "[INFO] Environment preflight: PASS"
)

print(
    "[INFO] Opening H5AD with backed='r'"
)

print(
    f"[INFO] {H5AD}"
)

adata = None

try:

    adata = ad.read_h5ad(
        H5AD,
        backed="r",
    )

    if not adata.isbacked:
        raise RuntimeError(
            "AnnData was not opened "
            "in backed mode."
        )

    summary.append(
        "AnnData structure"
    )
    summary.append("-" * 72)

    summary.append(
        f"n_obs: {adata.n_obs:,}"
    )

    summary.append(
        f"n_vars: {adata.n_vars:,}"
    )

    summary.append(
        f"isbacked: {adata.isbacked}"
    )

    summary.append(
        f"X backend type: "
        f"{type(adata.X).__name__}"
    )

    summary.append(
        "layers: "
        + (
            ", ".join(
                map(
                    str,
                    adata.layers.keys(),
                )
            )
            or "<none>"
        )
    )

    summary.append(
        "obsm: "
        + (
            ", ".join(
                map(
                    str,
                    adata.obsm.keys(),
                )
            )
            or "<none>"
        )
    )

    summary.append(
        "obsp: "
        + (
            ", ".join(
                map(
                    str,
                    adata.obsp.keys(),
                )
            )
            or "<none>"
        )
    )

    summary.append(
        "uns keys: "
        + (
            ", ".join(
                map(
                    str,
                    adata.uns.keys(),
                )
            )
            or "<none>"
        )
    )

    summary.append(
        f"raw present: "
        f"{adata.raw is not None}"
    )

    summary.append("")

    # --------------------------------------------------------
    # Dataset-size audit
    # --------------------------------------------------------

    summary.append(
        "Dataset-size check"
    )
    summary.append("-" * 72)

    if adata.n_obs == EXPECTED_N_OBS:

        summary.append(
            "[PASS] n_obs matches "
            f"expected value: {EXPECTED_N_OBS:,}"
        )

    else:

        summary.append(
            "[WARNING] n_obs differs "
            "from expected value."
        )

        summary.append(
            f"Expected: {EXPECTED_N_OBS:,}"
        )

        summary.append(
            f"Observed: {adata.n_obs:,}"
        )

    summary.append("")

    # ========================================================
    # 3. Read OBS metadata only
    # ========================================================

    print(
        "[INFO] Reading .obs metadata"
    )

    obs = adata.obs.copy()

    summary.append(
        "OBS overview"
    )
    summary.append("-" * 72)

    summary.append(
        f"Rows: {len(obs):,}"
    )

    summary.append(
        f"Columns: {len(obs.columns):,}"
    )

    summary.append(
        f"obs_names unique: "
        f"{obs.index.is_unique}"
    )

    summary.append(
        "obs_names missing: "
        f"{int(pd.isna(obs.index).sum())}"
    )

    summary.append("")

    if len(obs) != adata.n_obs:
        raise RuntimeError(
            "obs row count does not "
            "match adata.n_obs."
        )

    # ========================================================
    # 4. Export complete author OBS metadata
    # ========================================================

    obs_export = obs.copy()

    index_column = "obs_name"

    if index_column in obs_export.columns:
        index_column = "__obs_name__"

    obs_export.insert(
        0,
        index_column,
        obs_export.index.astype(str),
    )

    parquet_path = (
        OUTPUT_DIR
        / "author_obs_metadata.parquet"
    )

    tsv_fallback = (
        OUTPUT_DIR
        / "author_obs_metadata.tsv.gz"
    )

    summary.append(
        "Metadata export"
    )
    summary.append("-" * 72)

    try:

        obs_export.to_parquet(
            parquet_path,
            index=False,
            engine="pyarrow",
        )

        summary.append(
            f"[PASS] Parquet: {parquet_path}"
        )

    except Exception as exc:

        summary.append(
            "[WARNING] Parquet export failed."
        )

        summary.append(
            f"Reason: "
            f"{type(exc).__name__}: {exc}"
        )

        print(
            "[WARNING] Parquet export failed; "
            "writing compressed TSV fallback."
        )

        obs_export.to_csv(
            tsv_fallback,
            sep="\t",
            index=False,
            compression="gzip",
        )

        summary.append(
            f"[PASS] TSV fallback: "
            f"{tsv_fallback}"
        )

    summary.append(
        f"Rows exported: "
        f"{len(obs_export):,}"
    )

    summary.append("")

    # ========================================================
    # 5. Complete OBS column audit
    # ========================================================

    column_rows = []

    for col in obs.columns:

        s = obs[col]

        column_rows.append(
            {
                "column": col,
                "dtype": str(s.dtype),
                "n_rows": len(s),
                "n_nonmissing": int(
                    s.notna().sum()
                ),
                "n_missing": int(
                    s.isna().sum()
                ),
                "pct_missing": round(
                    float(
                        s.isna().mean()
                        * 100
                    ),
                    4,
                ),
                "n_unique_nonmissing":
                    safe_nunique(s),
                "example_values":
                    example_values(s),
            }
        )

    column_audit = pd.DataFrame(
        column_rows
    )

    column_audit.to_csv(
        OUTPUT_DIR
        / "obs_column_audit.tsv",
        sep="\t",
        index=False,
    )

    # ========================================================
    # 6. Expected-field discovery
    # ========================================================

    field_rows = []

    for field in EXPECTED_FIELDS:

        present = (
            field in obs.columns
        )

        field_rows.append(
            {
                "field": field,
                "present": present,
                "dtype": (
                    str(
                        obs[field].dtype
                    )
                    if present
                    else ""
                ),
                "n_missing": (
                    int(
                        obs[field]
                        .isna()
                        .sum()
                    )
                    if present
                    else ""
                ),
                "n_unique": (
                    safe_nunique(
                        obs[field]
                    )
                    if present
                    else ""
                ),
            }
        )

    pd.DataFrame(
        field_rows
    ).to_csv(
        OUTPUT_DIR
        / "expected_field_audit.tsv",
        sep="\t",
        index=False,
    )

    summary.append(
        "Expected-field discovery"
    )
    summary.append("-" * 72)

    for field in EXPECTED_FIELDS:

        status = (
            "FOUND"
            if field in obs.columns
            else "ABSENT"
        )

        summary.append(
            f"[{status}] {field}"
        )

    summary.append("")

    # ========================================================
    # 7. Cell-ID audit
    # ========================================================

    summary.append(
        "Cell ID consistency"
    )
    summary.append("-" * 72)

    if "cell_id" in obs.columns:

        cell_id = (
            obs["cell_id"]
            .astype("string")
        )

        index_values = pd.Series(
            obs.index.astype(str),
            index=obs.index,
            dtype="string",
        )

        mismatch = (
            cell_id.fillna("<NA>")
            !=
            index_values.fillna("<NA>")
        )

        n_mismatch = int(
            mismatch.sum()
        )

        pd.DataFrame(
            {
                "metric": [
                    "n_cells",
                    "cell_id_unique",
                    "obs_name_unique",
                    "cell_id_vs_obs_name_mismatch",
                ],
                "value": [
                    len(obs),
                    cell_id.is_unique,
                    obs.index.is_unique,
                    n_mismatch,
                ],
            }
        ).to_csv(
            OUTPUT_DIR
            / "cell_id_audit.tsv",
            sep="\t",
            index=False,
        )

        summary.append(
            f"cell_id unique: "
            f"{cell_id.is_unique}"
        )

        summary.append(
            "cell_id vs obs_name "
            f"mismatches: {n_mismatch:,}"
        )

    else:

        summary.append(
            "[INFO] No explicit "
            "cell_id column found."
        )

        summary.append(
            "obs_names will be retained "
            "as the original cell identifier."
        )

    summary.append("")

    # ========================================================
    # 8. Value counts for important metadata fields
    # ========================================================

    count_tables = []

    for field in KEY_FIELDS:

        if field not in obs.columns:
            continue

        counts = (
            obs[field]
            .astype("string")
            .fillna("<NA>")
            .value_counts(
                dropna=False
            )
            .rename_axis("value")
            .reset_index(
                name="n_cells"
            )
        )

        counts.insert(
            0,
            "field",
            field,
        )

        count_tables.append(
            counts
        )

    if count_tables:

        key_counts = pd.concat(
            count_tables,
            ignore_index=True,
        )

        key_counts.to_csv(
            OUTPUT_DIR
            / "key_field_value_counts.tsv",
            sep="\t",
            index=False,
        )

    # ========================================================
    # 9. Annotation hierarchy audit
    # ========================================================

    hierarchy_present = [
        field
        for field
        in ANNOTATION_HIERARCHY
        if field in obs.columns
    ]

    summary.append(
        "Annotation hierarchy"
    )
    summary.append("-" * 72)

    if hierarchy_present:

        summary.append(
            " -> ".join(
                hierarchy_present
            )
        )

        hierarchy_df = (
            obs[
                hierarchy_present
            ]
            .astype("string")
            .fillna("<NA>")
            .groupby(
                hierarchy_present,
                dropna=False,
                observed=True,
            )
            .size()
            .reset_index(
                name="n_cells"
            )
            .sort_values(
                "n_cells",
                ascending=False,
            )
        )

        hierarchy_df.to_csv(
            OUTPUT_DIR
            / "annotation_hierarchy_counts.tsv",
            sep="\t",
            index=False,
        )

        summary.append(
            "Unique hierarchy "
            "combinations: "
            f"{len(hierarchy_df):,}"
        )

    else:

        summary.append(
            "[INFO] None of the expected "
            "annotation hierarchy fields "
            "were found."
        )

    summary.append("")

    # ========================================================
    # 10. Sample-level metadata consistency
    # ========================================================

    summary.append(
        "Sample-level consistency"
    )
    summary.append("-" * 72)

    if "sample_id" in obs.columns:

        sample_series = (
            obs["sample_id"]
        )

        summary.append(
            "Unique non-missing "
            "sample_id: "
            f"{sample_series.nunique(dropna=True):,}"
        )

        consistency_rows = []

        grouped = obs.groupby(
            "sample_id",
            observed=True,
            dropna=False,
        )

        for field in SAMPLE_LEVEL_FIELDS:

            if field not in obs.columns:
                continue

            nunique = (
                grouped[field]
                .nunique(
                    dropna=False
                )
            )

            bad = nunique[
                nunique > 1
            ]

            consistency_rows.append(
                {
                    "field": field,
                    "n_samples_checked":
                        len(nunique),
                    "n_samples_with_gt1_value":
                        len(bad),
                    "max_values_per_sample":
                        int(nunique.max()),
                }
            )

            status = (
                "PASS"
                if len(bad) == 0
                else "WARNING"
            )

            summary.append(
                f"[{status}] "
                f"{field}: "
                f"{len(bad)} "
                "inconsistent samples"
            )

            if len(bad) > 0:

                (
                    bad
                    .rename(
                        "n_unique_values"
                    )
                    .reset_index()
                    .to_csv(
                        OUTPUT_DIR
                        / (
                            "inconsistent_samples__"
                            f"{field}.tsv"
                        ),
                        sep="\t",
                        index=False,
                    )
                )

        pd.DataFrame(
            consistency_rows
        ).to_csv(
            OUTPUT_DIR
            / "sample_field_consistency.tsv",
            sep="\t",
            index=False,
        )

    else:

        summary.append(
            "[INFO] sample_id "
            "column not found."
        )

    summary.append("")

    # ========================================================
    # 11. Donor overview
    # ========================================================

    summary.append(
        "Donor overview"
    )
    summary.append("-" * 72)

    if "Patient" in obs.columns:

        summary.append(
            "Unique non-missing "
            "Patient: "
            f"{obs['Patient'].nunique(dropna=True):,}"
        )

    else:

        summary.append(
            "[INFO] Patient "
            "column not found."
        )

    summary.append("")

    # ========================================================
    # 12. VAR metadata audit
    # ========================================================

    print(
        "[INFO] Reading .var metadata"
    )

    var = adata.var.copy()

    var_rows = []

    for col in var.columns:

        s = var[col]

        var_rows.append(
            {
                "column": col,
                "dtype": str(
                    s.dtype
                ),
                "n_missing": int(
                    s.isna().sum()
                ),
                "n_unique":
                    safe_nunique(s),
                "example_values":
                    example_values(s),
            }
        )

    pd.DataFrame(
        var_rows
    ).to_csv(
        OUTPUT_DIR
        / "var_column_audit.tsv",
        sep="\t",
        index=False,
    )

    summary.append(
        "Gene metadata"
    )
    summary.append("-" * 72)

    summary.append(
        f"n_vars: {len(var):,}"
    )

    summary.append(
        f"var_names unique: "
        f"{var.index.is_unique}"
    )

    summary.append(
        "var columns: "
        + (
            ", ".join(
                map(
                    str,
                    var.columns,
                )
            )
            if len(var.columns)
            else "<none>"
        )
    )

    summary.append("")

    # ========================================================
    # 13. Final status
    # ========================================================

    summary.append(
        "Final Phase-0 status"
    )
    summary.append("-" * 72)

    summary.append(
        "[PASS] H5AD opened successfully "
        "in read-only backed mode."
    )

    summary.append(
        "[PASS] OBS metadata inspected."
    )

    summary.append(
        "[PASS] VAR metadata inspected."
    )

    summary.append(
        "[PASS] Expression matrix was not "
        "materialized into memory."
    )

    absent = [
        field
        for field
        in EXPECTED_FIELDS
        if field not in obs.columns
    ]

    if absent:

        summary.append("")

        summary.append(
            "[INFO] Candidate fields absent: "
            + ", ".join(absent)
        )

        summary.append(
            "This is not considered a "
            "Phase-0 failure because "
            "Phase 0 performs schema discovery."
        )

    write_summary(
        summary
    )


# ============================================================
# Always close backed H5AD
# ============================================================

finally:

    if (
        adata is not None
        and adata.isbacked
    ):
        adata.file.close()


# ============================================================
# Console report
# ============================================================

print()
print("=" * 72)
print(
    "GSE282122 Phase 0 "
    "schema audit completed"
)
print("=" * 72)
print()

print("Summary:")
print(
    f"  {SUMMARY_FILE}"
)

print()
print("Outputs:")
print(
    f"  {OUTPUT_DIR}"
)

print()
print(
    "H5AD was opened in "
    "read-only backed mode."
)

print(
    "Expression matrix was "
    "NOT loaded into memory."
)