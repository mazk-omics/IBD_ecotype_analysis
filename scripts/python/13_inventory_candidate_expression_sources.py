#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13_inventory_candidate_expression_sources.py

Read-only inventory for Glial/Tuft cross-dataset expression validation.

Purpose:
    - inspect established audit outputs
    - identify expression source files/formats
    - validate candidate metadata parquet files
    - generate a machine-readable source manifest

No expression matrix is loaded.
No previous audit result is modified.
"""

from pathlib import Path
import sys
import pandas as pd


PROJECT_ROOT = Path("/home/mazekai/IBD_EcoTyper")

OUTDIR = (
    PROJECT_ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "01_candidate_expression_validation"
    / "00_source_inventory"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


FILES = {
    "SCP1884_expression_inventory": (
        PROJECT_ROOT
        / "03_reference"
        / "SCP1884"
        / "audit_v1"
        / "output"
        / "00_metadata_schema"
        / "expression_bundle_inventory.tsv"
    ),

    "SCP1884_candidate_metadata": (
        PROJECT_ROOT
        / "03_reference"
        / "SCP1884"
        / "audit_v1"
        / "output"
        / "02_candidate_taxonomy"
        / "scp1884_candidate_cells_provisional_taxonomy.parquet"
    ),

    "SCP259_expression_audit": (
        PROJECT_ROOT
        / "03_reference"
        / "SCP259"
        / "audit_v1"
        / "output"
        / "00_source_schema"
        / "expression_bundle_audit.tsv"
    ),

    "SCP259_gene_summary": (
        PROJECT_ROOT
        / "03_reference"
        / "SCP259"
        / "audit_v1"
        / "output"
        / "00_source_schema"
        / "gene_set_summary.tsv"
    ),

    "SCP259_candidate_metadata": (
        PROJECT_ROOT
        / "03_reference"
        / "SCP259"
        / "audit_v1"
        / "output"
        / "01_candidate_taxonomy_coverage"
        / "scp259_candidate_cells_provisional_taxonomy.parquet"
    ),

    "GSE282122_candidate_metadata": (
        PROJECT_ROOT
        / "03_reference"
        / "GSE282122"
        / "audit_v2"
        / "output"
        / "04_identity_coverage_v1"
        / "mapped_candidate_cells_v1.parquet"
    ),
}


GSE282122_H5AD = Path(
    "/home/mazekai/enteric_glia/data_scRNA/"
    "GSE282122/TAURUS_raw_counts_annotated_final.h5ad"
)


def fail(msg):
    print(f"\nERROR: {msg}", file=sys.stderr)
    sys.exit(1)


print("=" * 72)
print("CANDIDATE EXPRESSION SOURCE INVENTORY")
print("=" * 72)
print()


# ============================================================
# 1. Required-file existence
# ============================================================

for name, path in FILES.items():

    status = "FOUND" if path.exists() else "MISSING"

    print(
        f"{name:35s} {status:8s} {path}"
    )

    if not path.exists():
        fail(
            f"Required established audit file missing: {path}"
        )


print()

if not GSE282122_H5AD.exists():
    fail(
        "GSE282122 H5AD missing:\n"
        f"{GSE282122_H5AD}"
    )

print(
    "GSE282122 H5AD:"
)
print(
    f"  {GSE282122_H5AD}"
)
print(
    f"  size = "
    f"{GSE282122_H5AD.stat().st_size / 1024**3:.2f} GiB"
)


# ============================================================
# 2. Inspect TSV inventories
# ============================================================

tsv_keys = [
    "SCP1884_expression_inventory",
    "SCP259_expression_audit",
    "SCP259_gene_summary",
]

schema_rows = []

for key in tsv_keys:

    path = FILES[key]

    df = pd.read_csv(
        path,
        sep="\t",
    )

    print()
    print("=" * 72)
    print(key)
    print("=" * 72)

    print(
        f"Rows: {len(df):,}"
    )

    print(
        f"Columns ({len(df.columns)}):"
    )

    for col in df.columns:
        print(f"  - {col}")

    print()
    print("Preview:")
    print(
        df.head(20).to_string(
            index=False
        )
    )

    for col in df.columns:

        schema_rows.append(
            {
                "source": key,
                "column": col,
                "dtype": str(df[col].dtype),
                "n_nonmissing": int(
                    df[col].notna().sum()
                ),
                "n_unique": int(
                    df[col].nunique(
                        dropna=True
                    )
                ),
            }
        )


schema = pd.DataFrame(
    schema_rows
)

schema.to_csv(
    OUTDIR
    / "audit_inventory_schema.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 3. Find path-like values inside audit inventories
# ============================================================

path_records = []

PATH_HINTS = (
    "path",
    "file",
    "dir",
    "matrix",
    "expression",
    "bundle",
)

for key in [
    "SCP1884_expression_inventory",
    "SCP259_expression_audit",
]:

    df = pd.read_csv(
        FILES[key],
        sep="\t",
    )

    candidate_cols = [
        col
        for col in df.columns
        if any(
            hint in col.lower()
            for hint in PATH_HINTS
        )
    ]

    for col in candidate_cols:

        for value in (
            df[col]
            .dropna()
            .astype(str)
            .unique()
        ):

            value = value.strip()

            if not value:
                continue

            p = Path(value)

            looks_absolute = (
                value.startswith("/")
            )

            exists = (
                p.exists()
                if looks_absolute
                else False
            )

            path_records.append(
                {
                    "audit_source": key,
                    "column": col,
                    "value": value,
                    "is_absolute_path":
                        looks_absolute,
                    "exists":
                        exists,
                }
            )


path_table = pd.DataFrame(
    path_records
)

if len(path_table) > 0:

    path_table.to_csv(
        OUTDIR
        / "expression_source_path_candidates.tsv",
        sep="\t",
        index=False,
    )


# ============================================================
# 4. Inspect candidate metadata parquet files
# ============================================================

candidate_keys = [
    "GSE282122_candidate_metadata",
    "SCP1884_candidate_metadata",
    "SCP259_candidate_metadata",
]

candidate_summary = []

for key in candidate_keys:

    path = FILES[key]

    df = pd.read_parquet(path)

    print()
    print("=" * 72)
    print(key)
    print("=" * 72)

    print(
        f"Rows: {len(df):,}"
    )

    print(
        f"Columns ({len(df.columns)}):"
    )

    for col in df.columns:
        print(f"  - {col}")

    identity_candidates = [
        x
        for x in [
            "coarse_identity",
            "provisional_coarse_identity",
            "reference_identity",
        ]
        if x in df.columns
    ]

    if len(identity_candidates) != 1:

        fail(
            f"{key}: expected exactly one "
            "coarse identity column, found "
            f"{identity_candidates}"
        )

    identity_col = (
        identity_candidates[0]
    )

    for candidate in [
        "Glial",
        "Tuft",
    ]:

        n = int(
            (
                df[identity_col]
                .astype(str)
                == candidate
            ).sum()
        )

        candidate_summary.append(
            {
                "source": key,
                "identity_column":
                    identity_col,
                "candidate_identity":
                    candidate,
                "n_cells": n,
            }
        )

        print(
            f"{candidate}: {n:,} cells"
        )


candidate_summary = pd.DataFrame(
    candidate_summary
)

candidate_summary.to_csv(
    OUTDIR
    / "candidate_metadata_checkpoint.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 5. Expected checkpoints
# ============================================================

expected = {
    (
        "GSE282122_candidate_metadata",
        "Glial",
    ): 7971,

    (
        "GSE282122_candidate_metadata",
        "Tuft",
    ): 2216,

    (
        "SCP1884_candidate_metadata",
        "Glial",
    ): 4834,

    (
        "SCP1884_candidate_metadata",
        "Tuft",
    ): 3001,

    (
        "SCP259_candidate_metadata",
        "Glial",
    ): 1262,

    (
        "SCP259_candidate_metadata",
        "Tuft",
    ): 899,
}


for _, row in (
    candidate_summary.iterrows()
):

    key = (
        row["source"],
        row["candidate_identity"],
    )

    observed = int(
        row["n_cells"]
    )

    expected_n = expected[key]

    if observed != expected_n:

        fail(
            f"Candidate-cell checkpoint failed "
            f"for {key}: "
            f"expected {expected_n:,}, "
            f"observed {observed:,}"
        )


print()
print(
    "PASS: all Glial/Tuft candidate "
    "metadata checkpoints reproduced"
)


# ============================================================
# 6. Save simple source manifest
# ============================================================

manifest_rows = [
    {
        "dataset": "GSE282122",
        "metadata_file":
            str(
                FILES[
                    "GSE282122_candidate_metadata"
                ]
            ),
        "expression_source":
            str(GSE282122_H5AD),
        "expression_source_type":
            "H5AD_CSR_RAW_COUNTS",
    },
    {
        "dataset": "SCP1884",
        "metadata_file":
            str(
                FILES[
                    "SCP1884_candidate_metadata"
                ]
            ),
        "expression_source":
            "SEE_EXPRESSION_BUNDLE_INVENTORY",
        "expression_source_type":
            "TO_BE_RESOLVED_FROM_AUDIT",
    },
    {
        "dataset": "SCP259",
        "metadata_file":
            str(
                FILES[
                    "SCP259_candidate_metadata"
                ]
            ),
        "expression_source":
            "SEE_EXPRESSION_BUNDLE_AUDIT",
        "expression_source_type":
            "TO_BE_RESOLVED_FROM_AUDIT",
    },
]

manifest = pd.DataFrame(
    manifest_rows
)

manifest.to_csv(
    OUTDIR
    / "candidate_expression_source_manifest.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 72)
print("OUTPUT DIRECTORY")
print("=" * 72)
print(OUTDIR)

print()
print("=" * 72)
print("SUCCESS")
print("=" * 72)
print(
    "Candidate expression-source "
    "inventory completed."
)