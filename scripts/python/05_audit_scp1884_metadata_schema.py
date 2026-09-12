#!/usr/bin/env python3

from pathlib import Path
import sys
import zipfile
import importlib.metadata as md

import pandas as pd


# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT = Path("/home/mazekai/IBD_EcoTyper")

SOURCE_ROOT = Path(
    "/home/mazekai/enteric_glia/data_scRNA/SCP1884"
)

ZIP_PATH = (
    SOURCE_ROOT
    / "co_ti_cmb_metadata.zip"
)

ZIP_MEMBER = "co_ti_cmb_metadata.txt"

V2_PATH = (
    SOURCE_ROOT
    / "metadata"
    / "scp_metadata_combined.v2.txt"
)

EXPRESSION_ROOT = (
    SOURCE_ROOT
    / "expression"
)

EXPECTED_ENV_ROOT = Path(
    "/home/mazekai/miniconda3/envs/ibd_refaudit"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "03_reference"
    / "SCP1884"
    / "audit_v1"
    / "output"
    / "00_metadata_schema"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "metadata_schema_audit_summary.txt"
)

EXPECTED_TOTAL_CELLS = 720_633


# ============================================================
# Expression bundles
# ============================================================

BUNDLES = {
    "CO_STR": {
        "site_dir": "Colon",
        "compartment_dir": "Stromal",
        "expected_genes": 28_663,
        "expected_cells": 39_433,
    },
    "CO_IMM": {
        "site_dir": "Colon",
        "compartment_dir": "Immune",
        "expected_genes": 28_663,
        "expected_cells": 152_509,
    },
    "CO_EPI": {
        "site_dir": "Colon",
        "compartment_dir": "Epithelial",
        "expected_genes": 28_663,
        "expected_cells": 97_788,
    },
    "TI_STR": {
        "site_dir": "TerminalIleum",
        "compartment_dir": "Stromal",
        "expected_genes": 28_923,
        "expected_cells": 75_695,
    },
    "TI_IMM": {
        "site_dir": "TerminalIleum",
        "compartment_dir": "Immune",
        "expected_genes": 28_923,
        "expected_cells": 201_072,
    },
    "TI_EPI": {
        "site_dir": "TerminalIleum",
        "compartment_dir": "Epithelial",
        "expected_genes": 28_923,
        "expected_cells": 154_136,
    },
}


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


def normalize_string(s):
    return (
        s.astype("string")
        .str.strip()
    )


# ============================================================
# 0. Environment
# ============================================================

summary = []

summary.append(
    "SCP1884 Reference Audit"
)
summary.append(
    "Phase 0B: Metadata provenance and schema audit"
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
    f"pyarrow: {package_version('pyarrow')}"
)
summary.append("")


# ============================================================
# 1. Input existence
# ============================================================

for f in [
    ZIP_PATH,
    V2_PATH,
]:

    if not f.is_file():
        fail(
            f"Required input missing: {f}",
            summary,
        )

summary.append("Inputs")
summary.append("-" * 78)
summary.append(
    f"Author metadata ZIP: {ZIP_PATH}"
)
summary.append(
    f"Portal V2 metadata: {V2_PATH}"
)
summary.append("")


# ============================================================
# 2. Read AUTHOR metadata
#
# IMPORTANT:
# cell ID is stored as an unnamed row-name field.
# Explicit index_col=0 is therefore required.
# ============================================================

with zipfile.ZipFile(
    ZIP_PATH
) as z:

    if ZIP_MEMBER not in z.namelist():
        fail(
            f"{ZIP_MEMBER} not found in ZIP.",
            summary,
        )

    with z.open(
        ZIP_MEMBER
    ) as fh:

        author = pd.read_csv(
            fh,
            sep="\t",
            index_col=0,
            low_memory=False,
        )


author.index = (
    author.index
    .astype(str)
)

author.index.name = "cell_id"

summary.append(
    "Author metadata"
)
summary.append("-" * 78)

summary.append(
    f"Rows: {len(author):,}"
)
summary.append(
    f"Metadata columns: {len(author.columns)}"
)
summary.append(
    "Columns: "
    + ", ".join(author.columns)
)
summary.append(
    f"cell_id unique: {author.index.is_unique}"
)
summary.append(
    f"cell_id missing: {author.index.isna().sum()}"
)

if len(author) != EXPECTED_TOTAL_CELLS:
    fail(
        "Author metadata cell count mismatch: "
        f"{len(author):,}",
        summary,
    )

if len(author.columns) != 10:
    fail(
        "Expected 10 author metadata columns; "
        f"found {len(author.columns)}.",
        summary,
    )

if not author.index.is_unique:
    fail(
        "Author metadata cell IDs are not unique.",
        summary,
    )

summary.append(
    "[PASS] Author metadata = "
    "720,633 unique cells × 10 metadata columns."
)
summary.append("")


# ============================================================
# 3. Export canonical author metadata
# ============================================================

author_export = (
    author
    .reset_index()
)

author_export.to_parquet(
    OUTPUT_DIR
    / "author_metadata.parquet",
    index=False,
    engine="pyarrow",
)


# ============================================================
# 4. Author field audit
# ============================================================

field_rows = []

for col in author.columns:

    s = author[col]

    field_rows.append(
        {
            "field": col,
            "dtype": str(s.dtype),
            "n_missing": int(
                s.isna().sum()
            ),
            "n_unique": int(
                s.nunique(
                    dropna=True
                )
            ),
            "examples": " | ".join(
                s.dropna()
                .astype(str)
                .drop_duplicates()
                .head(8)
                .tolist()
            ),
        }
    )

pd.DataFrame(
    field_rows
).to_csv(
    OUTPUT_DIR
    / "author_metadata_field_audit.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 5. Read Single Cell Portal V2 metadata
#
# First data row is the schema/type row:
# NAME == "TYPE"
# ============================================================

v2 = pd.read_csv(
    V2_PATH,
    sep="\t",
    dtype=str,
    low_memory=False,
)

summary.append(
    "Portal V2 metadata"
)
summary.append("-" * 78)

summary.append(
    f"Rows including TYPE row: {len(v2):,}"
)
summary.append(
    f"Columns: {len(v2.columns)}"
)

if "NAME" not in v2.columns:
    fail(
        "V2 metadata has no NAME column.",
        summary,
    )

type_mask = (
    v2["NAME"]
    .astype(str)
    == "TYPE"
)

n_type_rows = int(
    type_mask.sum()
)

summary.append(
    f"TYPE rows: {n_type_rows}"
)

if n_type_rows != 1:
    fail(
        f"Expected exactly one TYPE row; "
        f"found {n_type_rows}.",
        summary,
    )

v2_cells = (
    v2.loc[
        ~type_mask
    ]
    .copy()
)

if len(v2_cells) != EXPECTED_TOTAL_CELLS:
    fail(
        "V2 real-cell count mismatch: "
        f"{len(v2_cells):,}",
        summary,
    )

if (
    v2_cells["NAME"]
    .isna()
    .any()
):
    fail(
        "V2 NAME contains missing cell IDs.",
        summary,
    )

if not (
    v2_cells["NAME"]
    .is_unique
):
    fail(
        "V2 NAME is not unique.",
        summary,
    )

summary.append(
    "[PASS] V2 contains one TYPE row + "
    "720,633 unique real cells."
)
summary.append("")


# ============================================================
# 6. Cell-ID equivalence
# ============================================================

author_ids = set(
    author.index.astype(str)
)

v2_ids = set(
    v2_cells[
        "NAME"
    ].astype(str)
)

author_only = (
    author_ids
    - v2_ids
)

v2_only = (
    v2_ids
    - author_ids
)

summary.append(
    "Cell-ID comparison"
)
summary.append("-" * 78)

summary.append(
    f"Author-only IDs: {len(author_only):,}"
)

summary.append(
    f"V2-only IDs: {len(v2_only):,}"
)

if author_only or v2_only:

    pd.Series(
        sorted(author_only),
        name="cell_id",
    ).to_csv(
        OUTPUT_DIR
        / "author_only_cell_ids.tsv",
        sep="\t",
        index=False,
    )

    pd.Series(
        sorted(v2_only),
        name="cell_id",
    ).to_csv(
        OUTPUT_DIR
        / "v2_only_cell_ids.tsv",
        sep="\t",
        index=False,
    )

    fail(
        "Author metadata and V2 metadata "
        "do not contain identical cell-ID sets.",
        summary,
    )

summary.append(
    "[PASS] Author and V2 metadata contain "
    "identical 720,633 cell IDs."
)
summary.append("")


# ============================================================
# 7. Compare shared semantic fields
# ============================================================

v2_compare = (
    v2_cells
    .set_index("NAME")
    .loc[
        author.index
    ]
)

comparison_spec = {
    "PubIDSample": "biosample_id",
    "n_genes": "n_genes",
    "n_counts": "n_counts",
    "Chem": "Chem",
    "Site": "Site",
    "Type": "Type",
    "PubID": "donor_id",
    "Layer": "Layer",
}

comparison_rows = []

for author_col, v2_col in comparison_spec.items():

    if (
        author_col not in author.columns
        or
        v2_col not in v2_compare.columns
    ):
        comparison_rows.append(
            {
                "author_field": author_col,
                "v2_field": v2_col,
                "n_mismatch": "FIELD_MISSING",
            }
        )
        continue

    a = normalize_string(
        author[author_col]
    )

    b = normalize_string(
        v2_compare[v2_col]
    )

    mismatch = (
        a.fillna("<NA>")
        !=
        b.fillna("<NA>")
    )

    comparison_rows.append(
        {
            "author_field": author_col,
            "v2_field": v2_col,
            "n_mismatch": int(
                mismatch.sum()
            ),
        }
    )


comparison_df = pd.DataFrame(
    comparison_rows
)

comparison_df.to_csv(
    OUTPUT_DIR
    / "author_vs_v2_field_comparison.tsv",
    sep="\t",
    index=False,
)

summary.append(
    "Author vs V2 shared fields"
)
summary.append("-" * 78)

for row in comparison_df.itertuples(
    index=False
):
    summary.append(
        f"{row.author_field} -> "
        f"{row.v2_field}: "
        f"{row.n_mismatch} mismatches"
    )

summary.append("")


# ============================================================
# 8. anno2 vs Celltype
#
# Do NOT assume they are identical.
# Record exact correspondence.
# ============================================================

anno_compare = pd.DataFrame(
    {
        "anno2":
            author["anno2"].astype(str),
        "Celltype":
            v2_compare["Celltype"].astype(str),
    }
)

anno_pairs = (
    anno_compare
    .value_counts(
        dropna=False
    )
    .rename(
        "n_cells"
    )
    .reset_index()
    .sort_values(
        "n_cells",
        ascending=False,
    )
)

anno_pairs.to_csv(
    OUTPUT_DIR
    / "anno2_vs_v2_Celltype.tsv",
    sep="\t",
    index=False,
)

summary.append(
    "Annotation comparison"
)
summary.append("-" * 78)

summary.append(
    f"Unique author anno2: "
    f"{author['anno2'].nunique(dropna=True)}"
)

summary.append(
    f"Unique V2 Celltype: "
    f"{v2_cells['Celltype'].nunique(dropna=True)}"
)

summary.append(
    f"Unique anno2–Celltype pairs: "
    f"{len(anno_pairs)}"
)

summary.append("")


# ============================================================
# 9. Author categorical value counts
# ============================================================

key_fields = [
    "PubIDSample",
    "anno_overall",
    "Chem",
    "Site",
    "Type",
    "PubID",
    "Layer",
    "anno2",
]

count_tables = []

for field in key_fields:

    if field not in author.columns:
        continue

    x = (
        author[field]
        .astype("string")
        .fillna("<NA>")
        .value_counts(
            dropna=False
        )
        .rename_axis(
            "value"
        )
        .reset_index(
            name="n_cells"
        )
    )

    x.insert(
        0,
        "field",
        field,
    )

    count_tables.append(
        x
    )

pd.concat(
    count_tables,
    ignore_index=True,
).to_csv(
    OUTPUT_DIR
    / "author_key_field_value_counts.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 10. V2 disease metadata
# ============================================================

v2_disease_fields = [
    "donor_id",
    "disease",
    "disease__ontology_label",
    "organ",
    "organ__ontology_label",
]

available_v2_disease = [
    c
    for c in v2_disease_fields
    if c in v2_cells.columns
]

v2_disease_summary = (
    v2_cells[
        available_v2_disease
    ]
    .drop_duplicates()
)

v2_disease_summary.to_csv(
    OUTPUT_DIR
    / "v2_donor_disease_metadata.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 11. Expression-bundle inventory
# ============================================================

bundle_rows = []

for bundle, info in BUNDLES.items():

    d = (
        EXPRESSION_ROOT
        / info["site_dir"]
        / info["compartment_dir"]
    )

    raw_files = list(
        d.glob(
            "raw*.mtx"
        )
    )

    norm_file = (
        d / "matrix.mtx"
    )

    feature_files = list(
        d.glob(
            "*.features.tsv"
        )
    )

    barcode_files = list(
        d.glob(
            "*.barcodes.tsv"
        )
    )

    if (
        len(raw_files) != 1
        or
        not norm_file.is_file()
        or
        len(feature_files) != 1
        or
        len(barcode_files) != 1
    ):
        fail(
            f"Incomplete expression bundle: {bundle}",
            summary,
        )

    feature_file = feature_files[0]
    barcode_file = barcode_files[0]

    n_features = sum(
        1 for _ in open(
            feature_file,
            "r",
            encoding="utf-8",
        )
    )

    n_barcodes = sum(
        1 for _ in open(
            barcode_file,
            "r",
            encoding="utf-8",
        )
    )

    bundle_rows.append(
        {
            "bundle": bundle,
            "site": info["site_dir"],
            "compartment": info["compartment_dir"],
            "expected_genes":
                info["expected_genes"],
            "observed_features":
                n_features,
            "expected_cells":
                info["expected_cells"],
            "observed_barcodes":
                n_barcodes,
        }
    )


bundle_df = pd.DataFrame(
    bundle_rows
)

bundle_df.to_csv(
    OUTPUT_DIR
    / "expression_bundle_inventory.tsv",
    sep="\t",
    index=False,
)

for row in bundle_df.itertuples(
    index=False
):

    if (
        row.expected_genes
        !=
        row.observed_features
    ):
        fail(
            f"{row.bundle}: feature-count mismatch.",
            summary,
        )

    if (
        row.expected_cells
        !=
        row.observed_barcodes
    ):
        fail(
            f"{row.bundle}: barcode-count mismatch.",
            summary,
        )

summary.append(
    "Expression bundles"
)
summary.append("-" * 78)

summary.append(
    f"Total expected cells: "
    f"{bundle_df['expected_cells'].sum():,}"
)

summary.append(
    f"Total barcode rows: "
    f"{bundle_df['observed_barcodes'].sum():,}"
)

summary.append(
    "[PASS] Six expression bundles match "
    "expected feature/barcode dimensions."
)
summary.append("")


# ============================================================
# 12. Final
# ============================================================

summary.append(
    "Final Phase-0B status"
)
summary.append("-" * 78)

summary.append(
    "[PASS] Author metadata parsed explicitly "
    "with cell IDs as row names."
)

summary.append(
    "[PASS] Author metadata contains "
    "720,633 cells."
)

summary.append(
    "[PASS] Portal V2 metadata contains "
    "720,633 real cells plus one TYPE row."
)

summary.append(
    "[PASS] Expression bundle inventory validated."
)

summary.append(
    "[INFO] Portal V2 is retained as a "
    "supplementary semantic/disease metadata source."
)

summary.append(
    "[INFO] Author co_ti_cmb_metadata is retained "
    "as the primary annotation source of truth."
)

write_summary(
    summary
)


print()
print("=" * 78)
print(
    "SCP1884 Phase 0B metadata schema audit completed"
)
print("=" * 78)

print()
print(
    f"Author cells       : {len(author):,}"
)
print(
    f"Author columns     : {len(author.columns)}"
)
print(
    f"V2 real cells      : {len(v2_cells):,}"
)
print(
    f"Author anno2 labels: "
    f"{author['anno2'].nunique()}"
)

print()
print(
    f"Summary:"
)
print(
    f"  {SUMMARY_FILE}"
)