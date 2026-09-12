#!/usr/bin/env python3

from pathlib import Path
import sys
import importlib.metadata as md

import pandas as pd


# ============================================================
# Configuration
# ============================================================

PROJECT = Path(
    "/home/mazekai/IBD_EcoTyper"
)

RAW = (
    PROJECT
    / "01_raw_data/scRNA/SCP259/raw/SCP259"
)

ENV_ROOT = Path(
    "/home/mazekai/miniconda3/envs/ibd_refaudit"
)

OUT = (
    PROJECT
    / "03_reference/SCP259/audit_v1/output/"
    "00_source_schema"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

META_FILE = (
    RAW
    / "metadata/all.meta2.txt"
)


BUNDLES = {

    "Epithelial": {

        "matrix":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc2348/"
                "gene_sorted-Epi.matrix.mtx"
            ),

        "barcodes":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc2348/"
                "Epi.barcodes2.tsv"
            ),

        "genes":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc2348/"
                "Epi.genes.tsv"
            ),

        "tsne":
            RAW
            / "cluster/Epi.tsne.txt",
    },

    "Stromal": {

        "matrix":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc2349/"
                "gene_sorted-Fib.matrix.mtx"
            ),

        "barcodes":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc2349/"
                "Fib.barcodes2.tsv"
            ),

        "genes":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc2349/"
                "Fib.genes.tsv"
            ),

        "tsne":
            RAW
            / "cluster/Fib.tsne.txt",
    },

    "Immune": {

        "matrix":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc234a/"
                "gene_sorted-Imm.matrix.mtx"
            ),

        "barcodes":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc234a/"
                "Imm.barcodes2.tsv"
            ),

        "genes":
            RAW
            / (
                "expression/"
                "5cdc540d328cee7a2efc234a/"
                "Imm.genes.tsv"
            ),

        "tsne":
            RAW
            / "cluster/Imm.tsne.txt",
    },
}


EXPECTED_CELLS = 365_492


# ============================================================
# Helpers
# ============================================================

def fail(msg):
    raise RuntimeError(msg)


def pkg(name):

    try:
        return md.version(name)

    except md.PackageNotFoundError:
        return "NOT INSTALLED"


def line_count(path):

    n = 0

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for _ in handle:
            n += 1

    return n


def matrix_market_dimensions(path):

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        banner = handle.readline().strip()

        if not banner.startswith(
            "%%MatrixMarket"
        ):
            fail(
                f"Not MatrixMarket: {path}"
            )

        for line in handle:

            line = line.strip()

            if (
                not line
                or line.startswith("%")
            ):
                continue

            parts = line.split()

            if len(parts) != 3:
                fail(
                    f"Unexpected dimension line "
                    f"in {path}: {line}"
                )

            return tuple(
                map(
                    int,
                    parts
                )
            )

    fail(
        f"No dimension line in {path}"
    )


def matrix_value_preview(
    path,
    n_values=20,
):

    values = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        dimension_seen = False

        for line in handle:

            line = line.strip()

            if (
                not line
                or line.startswith("%")
            ):
                continue

            if not dimension_seen:

                dimension_seen = True
                continue

            parts = line.split()

            if len(parts) >= 3:

                values.append(
                    parts[2]
                )

            if len(values) >= n_values:
                break

    return values


def collapse_values(series):

    return " | ".join(
        sorted(
            set(
                series
                .dropna()
                .astype(str)
            )
        )
    )


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


print(
    "=" * 78
)

print(
    "SCP259 PHASE 0A/0B SOURCE + SCHEMA AUDIT"
)

print(
    "=" * 78
)

print(
    "Python :",
    sys.version.split()[0],
)

print(
    "pandas :",
    pkg("pandas"),
)

print(
    "pyarrow:",
    pkg("pyarrow"),
)


# ============================================================
# 1. Validate files
# ============================================================

for bundle, spec in BUNDLES.items():

    for key, path in spec.items():

        if not path.is_file():

            fail(
                f"Missing {bundle}/{key}: {path}"
            )


if not META_FILE.is_file():

    fail(
        f"Missing metadata: {META_FILE}"
    )


print()
print(
    "[PASS] All expected source files exist."
)


# ============================================================
# 2. Matrix / barcode / gene closure
# ============================================================

bundle_rows = []

barcode_frames = []

gene_sets = {}


for bundle, spec in BUNDLES.items():

    n_gene_matrix, n_cell_matrix, nnz = (
        matrix_market_dimensions(
            spec["matrix"]
        )
    )

    barcodes = (
        pd.read_csv(
            spec["barcodes"],
            sep="\t",
            header=None,
            dtype=str,
        )
        .iloc[:, 0]
        .astype(str)
    )

    genes = (
        pd.read_csv(
            spec["genes"],
            sep="\t",
            header=None,
            dtype=str,
        )
        .iloc[:, 0]
        .astype(str)
    )

    n_tsne = line_count(
        spec["tsne"]
    )

    if len(genes) != n_gene_matrix:

        fail(
            f"{bundle}: "
            f"matrix genes={n_gene_matrix:,}, "
            f"gene file={len(genes):,}."
        )

    if len(barcodes) != n_cell_matrix:

        fail(
            f"{bundle}: "
            f"matrix cells={n_cell_matrix:,}, "
            f"barcode file={len(barcodes):,}."
        )

    if n_tsne != n_cell_matrix:

        fail(
            f"{bundle}: "
            f"tSNE lines={n_tsne:,}, "
            f"matrix cells={n_cell_matrix:,}."
        )

    if not barcodes.is_unique:

        duplicates = (
            barcodes[
                barcodes.duplicated(
                    keep=False
                )
            ]
        )

        pd.DataFrame(
            {
                "cell_id":
                    duplicates
            }
        ).to_csv(
            OUT
            / f"ERROR_{bundle}_duplicate_barcodes.tsv",
            sep="\t",
            index=False,
        )

        fail(
            f"{bundle}: duplicate barcodes."
        )

    gene_sets[
        bundle
    ] = set(
        genes
    )

    barcode_frames.append(
        pd.DataFrame(
            {
                "cell_id":
                    barcodes,

                "expression_bundle":
                    bundle,
            }
        )
    )

    bundle_rows.append(
        {
            "expression_bundle":
                bundle,

            "n_genes":
                n_gene_matrix,

            "n_cells":
                n_cell_matrix,

            "nnz":
                nnz,

            "genes_unique":
                genes.nunique(),

            "barcodes_unique":
                barcodes.nunique(),

            "tsne_lines":
                n_tsne,

            "matrix_value_preview":
                " ".join(
                    matrix_value_preview(
                        spec["matrix"]
                    )
                ),
        }
    )

    print()
    print(
        f"[PASS] {bundle}: "
        f"{n_gene_matrix:,} genes × "
        f"{n_cell_matrix:,} cells; "
        f"nnz={nnz:,}"
    )


bundle_audit = pd.DataFrame(
    bundle_rows
)

bundle_audit.to_csv(
    OUT
    / "expression_bundle_audit.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 3. Global barcode audit
# ============================================================

barcode_map = pd.concat(
    barcode_frames,
    ignore_index=True,
)


if len(
    barcode_map
) != EXPECTED_CELLS:

    fail(
        f"Combined barcode rows="
        f"{len(barcode_map):,}; "
        f"expected={EXPECTED_CELLS:,}."
    )


duplicate_global = (
    barcode_map[
        "cell_id"
    ]
    .duplicated(
        keep=False
    )
)


if duplicate_global.any():

    barcode_map.loc[
        duplicate_global
    ].to_csv(
        OUT
        / "ERROR_barcode_overlap_between_bundles.tsv",
        sep="\t",
        index=False,
    )

    fail(
        "Cell IDs overlap between expression bundles."
    )


print()
print(
    "[PASS] Three expression bundles contain "
    "365,492 globally unique cell IDs."
)


# ============================================================
# 4. Gene-set relationships
# ============================================================

common_genes = set.intersection(
    *gene_sets.values()
)

union_genes = set.union(
    *gene_sets.values()
)


gene_summary = pd.DataFrame(
    [
        {
            "metric":
                "Epithelial_genes",

            "n_genes":
                len(
                    gene_sets[
                        "Epithelial"
                    ]
                ),
        },

        {
            "metric":
                "Stromal_genes",

            "n_genes":
                len(
                    gene_sets[
                        "Stromal"
                    ]
                ),
        },

        {
            "metric":
                "Immune_genes",

            "n_genes":
                len(
                    gene_sets[
                        "Immune"
                    ]
                ),
        },

        {
            "metric":
                "intersection_all_3",

            "n_genes":
                len(
                    common_genes
                ),
        },

        {
            "metric":
                "union_all_3",

            "n_genes":
                len(
                    union_genes
                ),
        },
    ]
)


gene_summary.to_csv(
    OUT
    / "gene_set_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 5. Metadata schema
# ============================================================

meta_raw = pd.read_csv(
    META_FILE,
    sep="\t",
    dtype=str,
    low_memory=False,
)


print()
print(
    "Metadata columns:"
)

print(
    list(
        meta_raw.columns
    )
)


expected_columns = [
    "NAME",
    "Cluster",
    "nGene",
    "nUMI",
    "Subject",
    "Health",
    "Location",
    "Sample",
]


missing = [
    x
    for x in expected_columns
    if x not in meta_raw.columns
]


if missing:

    fail(
        f"Missing metadata columns: {missing}"
    )


type_mask = (
    meta_raw[
        "NAME"
    ]
    .astype(str)
    .eq(
        "TYPE"
    )
)


if int(
    type_mask.sum()
) != 1:

    fail(
        "Expected exactly one TYPE schema row; "
        f"found {int(type_mask.sum())}."
    )


type_row = (
    meta_raw.loc[
        type_mask
    ]
    .copy()
)


meta = (
    meta_raw.loc[
        ~type_mask
    ]
    .copy()
)


if len(meta) != EXPECTED_CELLS:

    fail(
        f"Real metadata rows="
        f"{len(meta):,}; "
        f"expected={EXPECTED_CELLS:,}."
    )


meta[
    "NAME"
] = (
    meta[
        "NAME"
    ]
    .astype(str)
)


if not meta[
    "NAME"
].is_unique:

    fail(
        "Metadata NAME is not unique."
    )


print()
print(
    "[PASS] Metadata contains "
    "1 TYPE row + 365,492 unique real cells."
)


# ============================================================
# 6. Exact barcode ↔ metadata closure
# ============================================================

metadata_ids = set(
    meta[
        "NAME"
    ]
)

barcode_ids = set(
    barcode_map[
        "cell_id"
    ]
)


metadata_only = (
    metadata_ids
    - barcode_ids
)

barcode_only = (
    barcode_ids
    - metadata_ids
)


if metadata_only:

    pd.Series(
        sorted(
            metadata_only
        ),
        name="cell_id",
    ).to_csv(
        OUT
        / "ERROR_metadata_only_cells.tsv",
        sep="\t",
        index=False,
    )


if barcode_only:

    pd.Series(
        sorted(
            barcode_only
        ),
        name="cell_id",
    ).to_csv(
        OUT
        / "ERROR_barcode_only_cells.tsv",
        sep="\t",
        index=False,
    )


if metadata_only or barcode_only:

    fail(
        "Barcode ↔ metadata closure failed: "
        f"metadata_only={len(metadata_only)}, "
        f"barcode_only={len(barcode_only)}."
    )


print(
    "[PASS] Expression barcodes exactly equal "
    "metadata NAME cell IDs."
)


# ============================================================
# 7. Add expression-bundle provenance
# ============================================================

meta = (
    meta.merge(
        barcode_map,
        left_on="NAME",
        right_on="cell_id",
        how="left",
        validate="one_to_one",
    )
)


if meta[
    "expression_bundle"
].isna().any():

    fail(
        "Some metadata cells lack expression bundle."
    )


# ============================================================
# 8. Metadata cardinalities
# ============================================================

print()
print(
    "=" * 78
)

print(
    "METADATA CARDINALITIES"
)

print(
    "=" * 78
)


for col in [
    "Cluster",
    "Subject",
    "Health",
    "Location",
    "Sample",
]:

    print()
    print(
        f"===== {col} ====="
    )

    print(
        meta[
            col
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )


# ============================================================
# 9. Sample consistency
# ============================================================

sample_consistency_rows = []


for field in [
    "Subject",
    "Health",
    "Location",
]:

    x = (
        meta.groupby(
            "Sample",
            observed=True,
        )[field]
        .nunique(
            dropna=False
        )
    )

    sample_consistency_rows.append(
        {
            "field":
                field,

            "max_unique_per_sample":
                int(
                    x.max()
                ),

            "n_inconsistent_samples":
                int(
                    (x > 1).sum()
                ),
        }
    )


sample_consistency = pd.DataFrame(
    sample_consistency_rows
)


sample_consistency.to_csv(
    OUT
    / "sample_field_consistency.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 10. Subject consistency
# ============================================================

subject_health = (
    meta.groupby(
        "Subject",
        observed=True,
    )[
        "Health"
    ]
    .nunique(
        dropna=False
    )
)


if (
    subject_health > 1
).any():

    print()
    print(
        "[NOTE] Some subjects contain more than "
        "one Health state."
    )

else:

    print()
    print(
        "[PASS] Health is subject-consistent."
    )


# ============================================================
# 11. Sample manifest
# ============================================================

sample_manifest = (
    meta.groupby(
        "Sample",
        observed=True,
    )
    .agg(

        subject=(
            "Subject",
            "first",
        ),

        health=(
            "Health",
            "first",
        ),

        location=(
            "Location",
            "first",
        ),

        n_cells=(
            "NAME",
            "size",
        ),

        n_clusters=(
            "Cluster",
            "nunique",
        ),

        expression_bundles=(
            "expression_bundle",
            collapse_values,
        ),
    )
    .reset_index()
)


sample_manifest.to_csv(
    OUT
    / "sample_manifest.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 12. Bundle × metadata relationships
# ============================================================

bundle_location = (
    meta.groupby(
        [
            "expression_bundle",
            "Location",
        ],
        observed=True,
    )
    .agg(

        n_cells=(
            "NAME",
            "size",
        ),

        n_samples=(
            "Sample",
            "nunique",
        ),

        n_subjects=(
            "Subject",
            "nunique",
        ),
    )
    .reset_index()
)


bundle_location.to_csv(
    OUT
    / "expression_bundle_by_location.tsv",
    sep="\t",
    index=False,
)


bundle_health = (
    meta.groupby(
        [
            "expression_bundle",
            "Health",
        ],
        observed=True,
    )
    .agg(

        n_cells=(
            "NAME",
            "size",
        ),

        n_samples=(
            "Sample",
            "nunique",
        ),

        n_subjects=(
            "Subject",
            "nunique",
        ),
    )
    .reset_index()
)


bundle_health.to_csv(
    OUT
    / "expression_bundle_by_health.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 13. Cluster audit
# ============================================================

cluster_summary = (
    meta.groupby(
        "Cluster",
        observed=True,
    )
    .agg(

        n_cells=(
            "NAME",
            "size",
        ),

        n_samples=(
            "Sample",
            "nunique",
        ),

        n_subjects=(
            "Subject",
            "nunique",
        ),

        health_states=(
            "Health",
            collapse_values,
        ),

        locations=(
            "Location",
            collapse_values,
        ),

        expression_bundles=(
            "expression_bundle",
            collapse_values,
        ),
    )
    .reset_index()
    .sort_values(
        "n_cells",
        ascending=False,
    )
)


cluster_summary.to_csv(
    OUT
    / "cluster_taxonomy_audit.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 14. Save cell-level metadata
# ============================================================

meta.to_parquet(
    OUT
    / "scp259_metadata_with_expression_bundle.parquet",
    index=False,
    engine="pyarrow",
)


type_row.to_csv(
    OUT
    / "metadata_TYPE_schema_row.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 15. Console report
# ============================================================

print()
print(
    "=" * 78
)

print(
    "SCP259 PHASE 0 SUMMARY"
)

print(
    "=" * 78
)

print(
    f"Cells       : {len(meta):,}"
)

print(
    f"Subjects    : {meta['Subject'].nunique()}"
)

print(
    f"Samples     : {meta['Sample'].nunique()}"
)

print(
    f"Clusters    : {meta['Cluster'].nunique()}"
)

print()

print(
    "Expression bundles:"
)

print(
    meta[
        "expression_bundle"
    ]
    .value_counts()
    .to_string()
)

print()

print(
    "Gene sets:"
)

print(
    gene_summary.to_string(
        index=False
    )
)

print()

print(
    "Sample consistency:"
)

print(
    sample_consistency.to_string(
        index=False
    )
)

print()

print(
    "Bundle × Location:"
)

print(
    bundle_location.to_string(
        index=False
    )
)

print()

print(
    "Bundle × Health:"
)

print(
    bundle_health.to_string(
        index=False
    )
)

print()

print(
    "Cluster taxonomy:"
)

print(
    cluster_summary.to_string(
        index=False
    )
)

print()

print(
    "Matrix value previews:"
)

print(
    bundle_audit[
        [
            "expression_bundle",
            "matrix_value_preview",
        ]
    ].to_string(
        index=False
    )
)

print()

print(
    "PHASE 0 STATUS: PASS"
)

print(
    f"Outputs: {OUT}"
)
