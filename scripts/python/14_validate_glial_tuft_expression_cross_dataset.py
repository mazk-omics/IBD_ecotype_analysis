#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
14_validate_glial_tuft_expression_cross_dataset.py

Expression-level cross-dataset validation of Glial and Tuft identities.

Datasets
--------
GSE282122
SCP1884
SCP259

Principles
----------
1. Reuse canonical audited cell-level metadata.
2. Use raw counts.
3. Do NOT load multi-GB Matrix Market files entirely into RAM.
4. Build donor/site-level pseudobulk profiles.
5. Compare candidate identity against biologically appropriate compartment
   backgrounds WITHIN each dataset.
6. Compare candidate-vs-background signatures across datasets.

Glial background:
    Fibroblast
    Pericyte
    Endothelial

Tuft background:
    Colonic_absorptive
    Ileal_absorptive
    Goblet
    Stem_TA_progenitor
    Paneth
    Enteroendocrine

Outputs are evidence for ontology freezing; the script does not silently
rewrite the reference taxonomy.
"""

from pathlib import Path
import gc
import sys

import numpy as np
import pandas as pd

try:
    import h5py
except ImportError:
    raise SystemExit(
        "ERROR: h5py is required in ibd_refaudit."
    )


# ============================================================
# 0. Configuration
# ============================================================

ROOT = Path("/home/mazekai/IBD_EcoTyper")

OUT = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "01_candidate_expression_validation"
    / "01_expression_validation"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ------------------------------------------------------------
# GSE282122
# ------------------------------------------------------------

GSE_H5AD = Path(
    "/home/mazekai/enteric_glia/data_scRNA/"
    "GSE282122/"
    "TAURUS_raw_counts_annotated_final.h5ad"
)

GSE_META = (
    ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "04_identity_coverage_v1"
    / "mapped_candidate_cells_v1.parquet"
)


# ------------------------------------------------------------
# SCP1884
# ------------------------------------------------------------

SCP1884_ROOT = Path(
    "/home/mazekai/enteric_glia/data_scRNA/"
    "SCP1884/expression"
)

SCP1884_META = (
    ROOT
    / "03_reference"
    / "SCP1884"
    / "audit_v1"
    / "output"
    / "02_candidate_taxonomy"
    / "scp1884_candidate_cells_provisional_taxonomy.parquet"
)


# ------------------------------------------------------------
# SCP259
# ------------------------------------------------------------

SCP259_EXPR = (
    ROOT
    / "01_raw_data"
    / "scRNA"
    / "SCP259"
    / "raw"
    / "SCP259"
    / "expression"
)

SCP259_META = (
    ROOT
    / "03_reference"
    / "SCP259"
    / "audit_v1"
    / "output"
    / "01_candidate_taxonomy_coverage"
    / "scp259_candidate_cells_provisional_taxonomy.parquet"
)


# ------------------------------------------------------------
# Candidate definitions
# ------------------------------------------------------------

CANDIDATES = {
    "Glial": {
        "background": {
            "Fibroblast",
            "Pericyte",
            "Endothelial",
        },
        "canonical_markers": [
            "S100B",
            "SLC1A3",
            "PLP1",
            "SOX10",
            "SLC1A1",
            "GFAP",
        ],
    },

    "Tuft": {
        "background": {
            "Colonic_absorptive",
            "Ileal_absorptive",
            "Goblet",
            "Stem_TA_progenitor",
            "Paneth",
            "Enteroendocrine",
        },
        "canonical_markers": [
            "POU2F3",
            "TRPM5",
            "GNG13",
            "AVIL",
            "SH2D6",
            "IL17RB",
        ],
    },
}


# ------------------------------------------------------------
# Pseudobulk eligibility
# ------------------------------------------------------------

MIN_CANDIDATE_CELLS_PER_UNIT = 3
MIN_BACKGROUND_CELLS_PER_UNIT = 20
MIN_ELIGIBLE_UNITS = 10

MTX_CHUNK_ROWS = 1_000_000
H5AD_CELL_CHUNK = 5000

TOP_N_MARKERS = 100


# ============================================================
# 1. Utilities
# ============================================================

def fail(msg):
    raise RuntimeError(msg)


def check_exists(path):
    if not Path(path).exists():
        fail(f"Missing required file: {path}")


def clean_string(x):
    return (
        x.astype("string")
        .fillna("NA")
        .str.strip()
    )


def make_unit_key(donor, site):
    return (
        clean_string(donor)
        + "||"
        + clean_string(site)
    )


def is_technical_gene(gene):
    g = str(gene).upper()

    return (
        g.startswith("MT-")
        or g.startswith("RPL")
        or g.startswith("RPS")
        or g.startswith("HBA")
        or g.startswith("HBB")
    )


def unique_match(root, pattern):
    hits = list(
        Path(root).rglob(pattern)
    )

    if len(hits) != 1:
        fail(
            f"Expected exactly one match for {pattern} "
            f"under {root}; found {len(hits)}:\n"
            + "\n".join(map(str, hits))
        )

    return hits[0]


# ============================================================
# 2. Feature / barcode readers
# ============================================================

def read_features(path):
    """
    Headerless feature reader.

    10x-like >=2-column file:
        use column 2 as gene symbol.

    1-column file:
        use column 1.
    """

    x = pd.read_csv(
        path,
        sep="\t",
        header=None,
        dtype=str,
    )

    if x.shape[1] == 1:
        genes = x.iloc[:, 0]

    else:
        genes = x.iloc[:, 1]

    genes = (
        genes
        .astype(str)
        .str.strip()
    )

    if genes.isna().any():
        fail(
            f"Missing gene names in {path}"
        )

    if (genes == "").any():
        fail(
            f"Blank gene names in {path}"
        )

    codes, unique_genes = pd.factorize(
        genes,
        sort=False,
    )

    return (
        np.asarray(unique_genes, dtype=str),
        codes.astype(np.int64),
        len(genes),
    )


def read_barcodes(path):
    x = pd.read_csv(
        path,
        sep="\t",
        header=None,
        dtype=str,
    )

    barcodes = (
        x.iloc[:, 0]
        .astype(str)
        .str.strip()
        .to_numpy()
    )

    if len(set(barcodes)) != len(barcodes):
        fail(
            f"Duplicate barcodes in {path}"
        )

    return barcodes


# ============================================================
# 3. Matrix Market header
# ============================================================

def read_mtx_header(path):
    """
    Return:
        n_rows
        n_cols
        nnz
        number of lines to skip before entries
        MatrixMarket first line
    """

    with open(
        path,
        "rt",
        buffering=1024 * 1024,
    ) as fh:

        first = fh.readline().strip()

        if not first.startswith(
            "%%MatrixMarket"
        ):
            fail(
                f"Not a Matrix Market file: {path}"
            )

        line_count = 1

        while True:

            line = fh.readline()

            if not line:
                fail(
                    f"No dimension line in {path}"
                )

            line_count += 1
            line = line.strip()

            if (
                not line
                or line.startswith("%")
            ):
                continue

            parts = line.split()

            if len(parts) != 3:
                fail(
                    f"Invalid dimension line: "
                    f"{line}"
                )

            n_rows, n_cols, nnz = (
                map(int, parts)
            )

            return {
                "n_rows": n_rows,
                "n_cols": n_cols,
                "nnz": nnz,
                "skiprows": line_count,
                "header": first,
            }


# ============================================================
# 4. Metadata standardization
# ============================================================

def prepare_metadata(
    dataset,
    candidate,
):

    if dataset == "GSE282122":

        df = pd.read_parquet(
            GSE_META
        )

        identity_col = (
            "coarse_identity"
        )

        donor = df["Patient"]
        site = df["Site"]

        df["bundle_std"] = "H5AD"

    elif dataset == "SCP1884":

        df = pd.read_parquet(
            SCP1884_META
        )

        identity_col = (
            "provisional_coarse_identity"
        )

        donor = df["v2_donor_id"]
        site = df["site_canonical"]

        df["bundle_std"] = (
            clean_string(
                df["expression_bundle"]
            )
        )

    elif dataset == "SCP259":

        df = pd.read_parquet(
            SCP259_META
        )

        identity_col = (
            "coarse_identity"
        )

        donor = df["donor_id"]
        site = df["site"]

        df["bundle_std"] = (
            clean_string(
                df["expression_bundle"]
            )
        )

    else:
        fail(
            f"Unknown dataset {dataset}"
        )


    df["identity_std"] = (
        clean_string(
            df[identity_col]
        )
    )

    df["action_std"] = (
        clean_string(
            df[
                "primary_reference_action"
            ]
        )
    )

    df["cell_id_std"] = (
        clean_string(
            df["cell_id"]
        )
    )

    df["unit"] = make_unit_key(
        donor,
        site,
    )


    background = (
        CANDIDATES[
            candidate
        ]["background"]
    )


    is_candidate = (
        (df["identity_std"] == candidate)
        &
        (
            df["action_std"]
            == "HOLD_CROSS_DATASET"
        )
    )

    is_background = (
        df["identity_std"].isin(
            background
        )
        &
        (
            df["action_std"]
            == "KEEP"
        )
    )


    sub = df.loc[
        is_candidate
        | is_background
    ].copy()


    sub["group"] = np.where(
        sub["identity_std"]
        == candidate,
        "CANDIDATE",
        "BACKGROUND",
    )


    if sub["cell_id_std"].duplicated().any():
        bad = sub.loc[
            sub[
                "cell_id_std"
            ].duplicated(
                keep=False
            )
        ]

        fail(
            f"{dataset}/{candidate}: "
            "duplicate cell_id among selected "
            f"metadata.\n{bad.head()}"
        )


    counts = (
        sub
        .groupby(
            [
                "unit",
                "group",
            ],
            observed=True,
        )
        .size()
        .unstack(
            fill_value=0
        )
    )


    if "CANDIDATE" not in counts:
        counts[
            "CANDIDATE"
        ] = 0

    if "BACKGROUND" not in counts:
        counts[
            "BACKGROUND"
        ] = 0


    eligible = counts.loc[
        (
            counts["CANDIDATE"]
            >= MIN_CANDIDATE_CELLS_PER_UNIT
        )
        &
        (
            counts["BACKGROUND"]
            >= MIN_BACKGROUND_CELLS_PER_UNIT
        )
    ].copy()


    if len(eligible) < MIN_ELIGIBLE_UNITS:
        fail(
            f"{dataset}/{candidate}: only "
            f"{len(eligible)} eligible donor/site "
            "units after pseudobulk thresholds."
        )


    sub = sub.loc[
        sub["unit"].isin(
            eligible.index
        )
    ].copy()


    units = sorted(
        eligible.index.tolist()
    )


    group_names = []

    for unit in units:
        group_names.append(
            f"{unit}||CANDIDATE"
        )
        group_names.append(
            f"{unit}||BACKGROUND"
        )


    group_index = {
        x: i
        for i, x
        in enumerate(group_names)
    }


    sub["aggregate_key"] = (
        sub["unit"]
        + "||"
        + sub["group"]
    )

    sub["aggregate_index"] = (
        sub["aggregate_key"]
        .map(group_index)
    )


    if sub[
        "aggregate_index"
    ].isna().any():
        fail(
            f"{dataset}/{candidate}: "
            "aggregate-index assignment failed."
        )


    sub[
        "aggregate_index"
    ] = (
        sub[
            "aggregate_index"
        ].astype(int)
    )


    print(
        f"{dataset}/{candidate}: "
        f"{len(units)} eligible donor/site units; "
        f"{(sub['group'] == 'CANDIDATE').sum():,} "
        "candidate cells; "
        f"{(sub['group'] == 'BACKGROUND').sum():,} "
        "background cells"
    )


    return (
        sub,
        units,
        group_names,
        eligible,
    )


# ============================================================
# 5. Streaming Matrix Market pseudobulk
# ============================================================

def stream_mtx_pseudobulk(
    matrix_path,
    feature_path,
    barcode_path,
    metadata,
    group_names,
    bundle_name,
):

    print(
        f"\nReading bundle {bundle_name}"
    )
    print(
        f"  matrix: {matrix_path}"
    )


    header = read_mtx_header(
        matrix_path
    )

    genes, raw_to_unique, n_features = (
        read_features(
            feature_path
        )
    )

    barcodes = read_barcodes(
        barcode_path
    )


    if (
        header["n_rows"]
        != n_features
    ):
        fail(
            f"{bundle_name}: matrix rows "
            f"{header['n_rows']:,} != "
            f"features {n_features:,}"
        )


    if (
        header["n_cols"]
        != len(barcodes)
    ):
        fail(
            f"{bundle_name}: matrix cols "
            f"{header['n_cols']:,} != "
            f"barcodes {len(barcodes):,}"
        )


    # Candidate/background metadata belonging
    # to this expression bundle.

    meta = metadata.loc[
        metadata["bundle_std"]
        == bundle_name
    ].copy()


    if len(meta) == 0:
        fail(
            f"No selected metadata cells "
            f"for {bundle_name}"
        )


    barcode_index = pd.Index(
        barcodes
    )

    positions = barcode_index.get_indexer(
        meta["cell_id_std"]
    )


    missing = (
        positions < 0
    )

    if missing.any():

        examples = (
            meta.loc[
                missing,
                "cell_id_std",
            ]
            .head(20)
            .tolist()
        )

        fail(
            f"{bundle_name}: "
            f"{missing.sum():,} selected metadata "
            "cells not found in barcode file. "
            f"Examples: {examples}"
        )


    # column -> aggregate pseudobulk index
    col_to_group = np.full(
        header["n_cols"],
        -1,
        dtype=np.int32,
    )

    col_to_group[
        positions
    ] = (
        meta[
            "aggregate_index"
        ]
        .to_numpy(
            dtype=np.int32
        )
    )


    n_gene_unique = len(genes)
    n_groups = len(group_names)

    result = np.zeros(
        (
            n_gene_unique,
            n_groups,
        ),
        dtype=np.float64,
    )


    # --------------------------------------------------------
    # Chunked Matrix Market read
    # --------------------------------------------------------

    reader = pd.read_csv(
        matrix_path,
        sep=r"\s+",
        header=None,
        names=[
            "gene",
            "cell",
            "value",
        ],
        skiprows=header["skiprows"],
        chunksize=MTX_CHUNK_ROWS,
        engine="c",
    )


    observed_entries = 0
    used_entries = 0


    for chunk_i, chunk in enumerate(
        reader,
        start=1,
    ):

        observed_entries += len(
            chunk
        )

        cells = (
            chunk["cell"]
            .to_numpy(
                dtype=np.int64
            )
            - 1
        )

        groups = (
            col_to_group[
                cells
            ]
        )

        keep = (
            groups >= 0
        )


        if not np.any(keep):
            continue


        raw_gene = (
            chunk.loc[
                keep,
                "gene",
            ]
            .to_numpy(
                dtype=np.int64
            )
            - 1
        )

        gene_idx = (
            raw_to_unique[
                raw_gene
            ]
        )

        group_idx = (
            groups[
                keep
            ].astype(
                np.int64
            )
        )

        values = (
            chunk.loc[
                keep,
                "value",
            ]
            .to_numpy(
                dtype=np.float64
            )
        )


        # raw-count sanity check
        if (
            not np.all(
                np.isfinite(values)
            )
            or np.any(values < 0)
        ):
            fail(
                f"{bundle_name}: invalid "
                "matrix values."
            )


        flat_index = (
            gene_idx
            * n_groups
            + group_idx
        )

        sums = np.bincount(
            flat_index,
            weights=values,
            minlength=(
                n_gene_unique
                * n_groups
            ),
        )

        result += sums.reshape(
            n_gene_unique,
            n_groups,
        )

        used_entries += int(
            keep.sum()
        )


    if (
        observed_entries
        != header["nnz"]
    ):
        fail(
            f"{bundle_name}: expected "
            f"{header['nnz']:,} entries, read "
            f"{observed_entries:,}."
        )


    print(
        f"  PASS dimensions: "
        f"{header['n_rows']:,} genes × "
        f"{header['n_cols']:,} cells"
    )

    print(
        f"  selected nonzero entries: "
        f"{used_entries:,}"
    )


    counts = pd.DataFrame(
        result,
        index=genes,
        columns=group_names,
    )


    return counts


# ============================================================
# 6. H5AD helpers
# ============================================================

def decode_array(x):
    arr = np.asarray(x)

    if arr.dtype.kind in {
        "S",
        "O",
        "U",
    }:
        return np.array([
            (
                v.decode("utf-8")
                if isinstance(
                    v,
                    (bytes, bytearray)
                )
                else str(v)
            )
            for v in arr
        ])

    return arr.astype(str)


def read_h5ad_vector(
    h5,
    path,
):

    obj = h5[path]

    if isinstance(
        obj,
        h5py.Dataset,
    ):
        return decode_array(
            obj[:]
        )


    if isinstance(
        obj,
        h5py.Group,
    ):

        if (
            "categories" in obj
            and "codes" in obj
        ):

            cats = decode_array(
                obj[
                    "categories"
                ][:]
            )

            codes = (
                obj["codes"][:]
            )

            out = np.full(
                len(codes),
                "NA",
                dtype=object,
            )

            valid = (
                codes >= 0
            )

            out[valid] = (
                cats[
                    codes[valid]
                ]
            )

            return np.asarray(
                out,
                dtype=str,
            )


    fail(
        f"Unsupported H5AD vector: {path}"
    )


# ============================================================
# 7. GSE282122 H5AD pseudobulk
# ============================================================

def gse_pseudobulk(
    metadata,
    group_names,
):

    print(
        "\nStreaming GSE282122 H5AD..."
    )


    with h5py.File(
        GSE_H5AD,
        "r",
    ) as h5:

        if "/obs/_index" in h5:
            obs_names = (
                read_h5ad_vector(
                    h5,
                    "/obs/_index",
                )
            )
        else:
            fail(
                "GSE282122 H5AD missing "
                "/obs/_index"
            )


        if "/var/gene_symbol" in h5:
            raw_genes = (
                read_h5ad_vector(
                    h5,
                    "/var/gene_symbol",
                )
            )

        elif "/var/_index" in h5:
            raw_genes = (
                read_h5ad_vector(
                    h5,
                    "/var/_index",
                )
            )

        else:
            fail(
                "GSE282122 H5AD has no "
                "gene-symbol field."
            )


        gene_codes, genes = pd.factorize(
            pd.Series(
                raw_genes
            ).astype(str),
            sort=False,
        )

        gene_codes = (
            gene_codes.astype(
                np.int64
            )
        )

        genes = np.asarray(
            genes,
            dtype=str,
        )


        n_cells = len(
            obs_names
        )

        n_groups = len(
            group_names
        )


        obs_index = pd.Index(
            obs_names
        )

        positions = (
            obs_index.get_indexer(
                metadata[
                    "cell_id_std"
                ]
            )
        )


        missing = (
            positions < 0
        )

        if missing.any():

            examples = (
                metadata.loc[
                    missing,
                    "cell_id_std",
                ]
                .head(20)
                .tolist()
            )

            fail(
                "GSE282122 metadata/H5AD "
                f"cell mismatch: {missing.sum():,}. "
                f"Examples: {examples}"
            )


        row_to_group = np.full(
            n_cells,
            -1,
            dtype=np.int32,
        )

        row_to_group[
            positions
        ] = (
            metadata[
                "aggregate_index"
            ]
            .to_numpy(
                dtype=np.int32
            )
        )


        result = np.zeros(
            (
                len(genes),
                n_groups,
            ),
            dtype=np.float64,
        )


        indptr = (
            h5[
                "/X/indptr"
            ][:]
        )

        data_ds = h5[
            "/X/data"
        ]

        index_ds = h5[
            "/X/indices"
        ]


        selected_rows = (
            row_to_group
            >= 0
        )


        for start_row in range(
            0,
            n_cells,
            H5AD_CELL_CHUNK,
        ):

            end_row = min(
                start_row
                + H5AD_CELL_CHUNK,
                n_cells,
            )


            local_selected = np.flatnonzero(
                selected_rows[
                    start_row:end_row
                ]
            )


            if len(
                local_selected
            ) == 0:
                continue


            global_rows = (
                local_selected
                + start_row
            )


            data_start = int(
                indptr[
                    start_row
                ]
            )

            data_end = int(
                indptr[
                    end_row
                ]
            )


            chunk_data = (
                data_ds[
                    data_start:data_end
                ]
            )

            chunk_indices = (
                index_ds[
                    data_start:data_end
                ]
            )


            for row in global_rows:

                a = int(
                    indptr[row]
                    - data_start
                )

                b = int(
                    indptr[row + 1]
                    - data_start
                )


                raw_idx = (
                    chunk_indices[
                        a:b
                    ]
                )

                vals = (
                    chunk_data[
                        a:b
                    ]
                )


                if (
                    np.any(vals < 0)
                    or not np.all(
                        np.isfinite(vals)
                    )
                ):
                    fail(
                        "Invalid GSE282122 "
                        "raw-count values."
                    )


                unique_idx = (
                    gene_codes[
                        raw_idx
                    ]
                )

                grp = int(
                    row_to_group[
                        row
                    ]
                )


                np.add.at(
                    result[:, grp],
                    unique_idx,
                    vals,
                )


    print(
        "  PASS: GSE282122 selected "
        "raw counts aggregated"
    )


    return pd.DataFrame(
        result,
        index=genes,
        columns=group_names,
    )


# ============================================================
# 8. Dataset-specific expression sources
# ============================================================

def scp1884_pseudobulk(
    candidate,
    metadata,
    group_names,
):

    if candidate == "Glial":

        bundle_cfg = {
            "CO_STR": (
                SCP1884_ROOT
                / "Colon"
                / "Stromal"
            ),
            "TI_STR": (
                SCP1884_ROOT
                / "TerminalIleum"
                / "Stromal"
            ),
        }

    elif candidate == "Tuft":

        bundle_cfg = {
            "CO_EPI": (
                SCP1884_ROOT
                / "Colon"
                / "Epithelial"
            ),
            "TI_EPI": (
                SCP1884_ROOT
                / "TerminalIleum"
                / "Epithelial"
            ),
        }

    else:
        fail(candidate)


    combined = None


    for bundle, directory in (
        bundle_cfg.items()
    ):

        matrix = (
            directory
            / "raw.mtx"
        )

        features = (
            directory
            / f"{bundle}.features.tsv"
        )

        barcodes = (
            directory
            / f"{bundle}.barcodes.tsv"
        )


        for p in [
            matrix,
            features,
            barcodes,
        ]:
            check_exists(p)


        x = stream_mtx_pseudobulk(
            matrix,
            features,
            barcodes,
            metadata,
            group_names,
            bundle,
        )


        if combined is None:
            combined = x

        else:
            combined = combined.add(
                x,
                fill_value=0,
            )


        del x
        gc.collect()


    return combined


def scp259_pseudobulk(
    candidate,
    metadata,
    group_names,
):

    if candidate == "Glial":

        stem = "Fib"
        bundle = "Stromal"

    elif candidate == "Tuft":

        stem = "Epi"
        bundle = "Epithelial"

    else:
        fail(candidate)


    matrix = unique_match(
        SCP259_EXPR,
        f"gene_sorted-{stem}.matrix.mtx",
    )

    features = unique_match(
        SCP259_EXPR,
        f"{stem}.genes.tsv",
    )

    barcodes = unique_match(
        SCP259_EXPR,
        f"{stem}.barcodes2.tsv",
    )


    return stream_mtx_pseudobulk(
        matrix,
        features,
        barcodes,
        metadata,
        group_names,
        bundle,
    )


# ============================================================
# 9. Convert pseudobulk counts into expression signature
# ============================================================

def build_signature(
    counts,
    units,
):

    if (
        counts.values < 0
    ).any():
        fail(
            "Negative pseudobulk counts."
        )


    candidate_cols = [
        f"{u}||CANDIDATE"
        for u in units
    ]

    background_cols = [
        f"{u}||BACKGROUND"
        for u in units
    ]


    for col in (
        candidate_cols
        + background_cols
    ):
        if col not in counts.columns:
            fail(
                f"Missing pseudobulk column: "
                f"{col}"
            )


    cand = (
        counts[
            candidate_cols
        ]
        .to_numpy(
            dtype=np.float64
        )
    )

    bg = (
        counts[
            background_cols
        ]
        .to_numpy(
            dtype=np.float64
        )
    )


    cand_lib = (
        cand.sum(
            axis=0
        )
    )

    bg_lib = (
        bg.sum(
            axis=0
        )
    )


    if (
        np.any(cand_lib <= 0)
        or np.any(bg_lib <= 0)
    ):
        fail(
            "Zero pseudobulk library size."
        )


    cand_cpm = (
        cand
        / cand_lib
        * 1e6
    )

    bg_cpm = (
        bg
        / bg_lib
        * 1e6
    )


    cand_log = np.log2(
        cand_cpm + 1
    )

    bg_log = np.log2(
        bg_cpm + 1
    )


    unit_lfc = (
        cand_log
        - bg_log
    )


    out = pd.DataFrame(
        {
            "gene":
                counts.index.astype(str),

            "median_log2FC":
                np.median(
                    unit_lfc,
                    axis=1,
                ),

            "mean_log2FC":
                np.mean(
                    unit_lfc,
                    axis=1,
                ),

            "fraction_units_positive":
                np.mean(
                    unit_lfc > 0,
                    axis=1,
                ),

            "mean_candidate_CPM":
                np.mean(
                    cand_cpm,
                    axis=1,
                ),

            "mean_background_CPM":
                np.mean(
                    bg_cpm,
                    axis=1,
                ),
        }
    )


    out[
        "mean_max_CPM"
    ] = np.maximum(
        out[
            "mean_candidate_CPM"
        ],
        out[
            "mean_background_CPM"
        ],
    )


    out[
        "technical_gene"
    ] = (
        out["gene"]
        .map(
            is_technical_gene
        )
    )


    return out


# ============================================================
# 10. Spearman without scipy
# ============================================================

def spearman(x, y):

    x = pd.Series(
        x
    ).rank(
        method="average"
    )

    y = pd.Series(
        y
    ).rank(
        method="average"
    )

    return float(
        x.corr(y)
    )


# ============================================================
# 11. Pairwise signature comparison
# ============================================================

def compare_signatures(
    candidate,
    dataset_a,
    sig_a,
    dataset_b,
    sig_b,
):

    a = sig_a.loc[
        ~sig_a[
            "technical_gene"
        ]
    ].copy()

    b = sig_b.loc[
        ~sig_b[
            "technical_gene"
        ]
    ].copy()


    merged = a.merge(
        b,
        on="gene",
        suffixes=(
            "_a",
            "_b",
        ),
    )


    informative = merged.loc[
        (
            merged[
                "mean_max_CPM_a"
            ]
            >= 1
        )
        &
        (
            merged[
                "mean_max_CPM_b"
            ]
            >= 1
        )
        &
        (
            (
                merged[
                    "median_log2FC_a"
                ].abs()
                >= 0.25
            )
            |
            (
                merged[
                    "median_log2FC_b"
                ].abs()
                >= 0.25
            )
        )
    ].copy()


    if len(
        informative
    ) < 500:
        fail(
            f"{candidate}: too few "
            f"informative genes for "
            f"{dataset_a} vs {dataset_b}: "
            f"{len(informative)}"
        )


    rho = spearman(
        informative[
            "median_log2FC_a"
        ],
        informative[
            "median_log2FC_b"
        ],
    )


    direction_set = (
        informative.loc[
            (
                informative[
                    "median_log2FC_a"
                ].abs()
                >= 0.5
            )
            &
            (
                informative[
                    "median_log2FC_b"
                ].abs()
                >= 0.5
            )
        ]
    )


    if len(direction_set) > 0:

        direction_concordance = float(
            (
                np.sign(
                    direction_set[
                        "median_log2FC_a"
                    ]
                )
                ==
                np.sign(
                    direction_set[
                        "median_log2FC_b"
                    ]
                )
            ).mean()
        )

    else:
        direction_concordance = np.nan


    top_a = set(
        a.loc[
            (
                a[
                    "median_log2FC"
                ]
                > 0
            )
            &
            (
                a[
                    "mean_candidate_CPM"
                ]
                >= 1
            )
        ]
        .nlargest(
            TOP_N_MARKERS,
            "median_log2FC",
        )["gene"]
    )


    top_b = set(
        b.loc[
            (
                b[
                    "median_log2FC"
                ]
                > 0
            )
            &
            (
                b[
                    "mean_candidate_CPM"
                ]
                >= 1
            )
        ]
        .nlargest(
            TOP_N_MARKERS,
            "median_log2FC",
        )["gene"]
    )


    intersection = (
        top_a & top_b
    )

    union = (
        top_a | top_b
    )


    jaccard = (
        len(intersection)
        / len(union)
        if len(union)
        else np.nan
    )


    return {
        "candidate_identity":
            candidate,

        "dataset_a":
            dataset_a,

        "dataset_b":
            dataset_b,

        "n_common_genes":
            len(merged),

        "n_informative_genes":
            len(informative),

        "spearman_signature":
            rho,

        "n_direction_genes":
            len(direction_set),

        "direction_concordance":
            direction_concordance,

        "top_marker_overlap_n":
            len(intersection),

        "top_marker_jaccard":
            jaccard,

        "top_marker_overlap":
            " | ".join(
                sorted(
                    intersection
                )
            ),
    }


# ============================================================
# 12. Main
# ============================================================

print("=" * 72)
print("GLIAL / TUFT CROSS-DATASET EXPRESSION VALIDATION")
print("=" * 72)


for p in [
    GSE_H5AD,
    GSE_META,
    SCP1884_META,
    SCP259_META,
]:
    check_exists(p)


all_pairwise = []
all_marker_checks = []
all_dataset_summaries = []


for candidate in [
    "Glial",
    "Tuft",
]:

    print()
    print("=" * 72)
    print(
        f"CANDIDATE: {candidate}"
    )
    print("=" * 72)


    signatures = {}


    # --------------------------------------------------------
    # Process each dataset
    # --------------------------------------------------------

    for dataset in [
        "GSE282122",
        "SCP1884",
        "SCP259",
    ]:

        print()
        print("-" * 72)
        print(
            f"{dataset} / {candidate}"
        )
        print("-" * 72)


        meta, units, group_names, eligibility = (
            prepare_metadata(
                dataset,
                candidate,
            )
        )


        eligibility.to_csv(
            OUT
            / (
                f"{candidate}__"
                f"{dataset}__eligible_units.tsv"
            ),
            sep="\t",
        )


        if dataset == "GSE282122":

            counts = gse_pseudobulk(
                meta,
                group_names,
            )

        elif dataset == "SCP1884":

            counts = scp1884_pseudobulk(
                candidate,
                meta,
                group_names,
            )

        else:

            counts = scp259_pseudobulk(
                candidate,
                meta,
                group_names,
            )


        # Remove genes with zero counts everywhere.
        counts = counts.loc[
            counts.sum(
                axis=1
            ) > 0
        ]


        signature = build_signature(
            counts,
            units,
        )


        signature = signature.sort_values(
            "median_log2FC",
            ascending=False,
        )


        signature.to_csv(
            OUT
            / (
                f"{candidate}__"
                f"{dataset}__signature.tsv"
            ),
            sep="\t",
            index=False,
        )


        signatures[
            dataset
        ] = signature


        all_dataset_summaries.append(
            {
                "candidate_identity":
                    candidate,

                "dataset":
                    dataset,

                "n_eligible_units":
                    len(units),

                "n_candidate_cells":
                    int(
                        (
                            meta["group"]
                            == "CANDIDATE"
                        ).sum()
                    ),

                "n_background_cells":
                    int(
                        (
                            meta["group"]
                            == "BACKGROUND"
                        ).sum()
                    ),

                "n_expression_genes":
                    len(signature),

                "n_positive_genes_logFC_gt_0_5":
                    int(
                        (
                            signature[
                                "median_log2FC"
                            ]
                            >= 0.5
                        ).sum()
                    ),
            }
        )


        # ----------------------------------------------------
        # Canonical marker sanity check
        # ----------------------------------------------------

        marker_set = (
            CANDIDATES[
                candidate
            ]["canonical_markers"]
        )


        marker_table = (
            signature.loc[
                signature["gene"]
                .isin(marker_set),
                [
                    "gene",
                    "median_log2FC",
                    "fraction_units_positive",
                    "mean_candidate_CPM",
                    "mean_background_CPM",
                ],
            ]
            .copy()
        )


        marker_table[
            "candidate_identity"
        ] = candidate

        marker_table[
            "dataset"
        ] = dataset


        present_markers = set(
            marker_table["gene"]
        )


        for marker in marker_set:

            if marker not in present_markers:

                all_marker_checks.append(
                    {
                        "candidate_identity":
                            candidate,
                        "dataset":
                            dataset,
                        "gene":
                            marker,
                        "gene_present":
                            False,
                        "median_log2FC":
                            np.nan,
                        "fraction_units_positive":
                            np.nan,
                        "mean_candidate_CPM":
                            np.nan,
                        "mean_background_CPM":
                            np.nan,
                    }
                )

            else:

                row = (
                    marker_table.loc[
                        marker_table["gene"]
                        == marker
                    ]
                    .iloc[0]
                )

                all_marker_checks.append(
                    {
                        "candidate_identity":
                            candidate,
                        "dataset":
                            dataset,
                        "gene":
                            marker,
                        "gene_present":
                            True,
                        "median_log2FC":
                            row[
                                "median_log2FC"
                            ],
                        "fraction_units_positive":
                            row[
                                "fraction_units_positive"
                            ],
                        "mean_candidate_CPM":
                            row[
                                "mean_candidate_CPM"
                            ],
                        "mean_background_CPM":
                            row[
                                "mean_background_CPM"
                            ],
                    }
                )


        del counts
        del meta
        gc.collect()


    # --------------------------------------------------------
    # Cross-dataset pairwise comparisons
    # --------------------------------------------------------

    pairs = [
        (
            "GSE282122",
            "SCP1884",
        ),
        (
            "GSE282122",
            "SCP259",
        ),
        (
            "SCP1884",
            "SCP259",
        ),
    ]


    for a, b in pairs:

        result = compare_signatures(
            candidate,
            a,
            signatures[a],
            b,
            signatures[b],
        )

        all_pairwise.append(
            result
        )


# ============================================================
# 13. Save summaries
# ============================================================

dataset_summary = pd.DataFrame(
    all_dataset_summaries
)

dataset_summary.to_csv(
    OUT
    / "dataset_candidate_expression_summary.tsv",
    sep="\t",
    index=False,
)


pairwise = pd.DataFrame(
    all_pairwise
)

pairwise.to_csv(
    OUT
    / "pairwise_signature_concordance.tsv",
    sep="\t",
    index=False,
)


marker_checks = pd.DataFrame(
    all_marker_checks
)

marker_checks.to_csv(
    OUT
    / "canonical_marker_sanity_check.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 14. Candidate-level evidence summary
# ============================================================

final_rows = []


for candidate in [
    "Glial",
    "Tuft",
]:

    p = pairwise.loc[
        pairwise[
            "candidate_identity"
        ]
        == candidate
    ]


    m = marker_checks.loc[
        (
            marker_checks[
                "candidate_identity"
            ]
            == candidate
        )
        &
        (
            marker_checks[
                "gene_present"
            ]
        )
    ].copy()


    # Canonical marker support:
    # positive logFC in >= 2/3 datasets.

    marker_consensus = (
        m.assign(
            positive=(
                m[
                    "median_log2FC"
                ]
                > 0
            )
        )
        .groupby(
            "gene"
        )["positive"]
        .sum()
    )


    n_marker_consensus = int(
        (
            marker_consensus
            >= 2
        ).sum()
    )


    median_rho = float(
        p[
            "spearman_signature"
        ].median()
    )

    min_rho = float(
        p[
            "spearman_signature"
        ].min()
    )

    median_direction = float(
        p[
            "direction_concordance"
        ].median()
    )

    median_overlap = float(
        p[
            "top_marker_overlap_n"
        ].median()
    )


    # These are technical evidence grades,
    # NOT automatic biological truth.

    if (
        min_rho >= 0.15
        and median_direction >= 0.65
        and median_overlap >= 10
        and n_marker_consensus >= 3
    ):

        evidence = (
            "STRONG_CROSS_DATASET_SUPPORT"
        )

    elif (
        median_rho >= 0.15
        and median_direction >= 0.60
        and n_marker_consensus >= 2
    ):

        evidence = (
            "CROSS_DATASET_SUPPORT"
        )

    else:

        evidence = (
            "REVIEW_EXPRESSION_CONSISTENCY"
        )


    final_rows.append(
        {
            "candidate_identity":
                candidate,

            "median_pairwise_spearman":
                median_rho,

            "minimum_pairwise_spearman":
                min_rho,

            "median_direction_concordance":
                median_direction,

            "median_top100_overlap_n":
                median_overlap,

            "n_canonical_markers_positive_in_2plus_datasets":
                n_marker_consensus,

            "expression_evidence_grade":
                evidence,

            "ontology_action":
                (
                    "ELIGIBLE_FOR_FINAL_KEEP_REVIEW"
                    if evidence
                    != "REVIEW_EXPRESSION_CONSISTENCY"
                    else
                    "RETAIN_HOLD_PENDING_REVIEW"
                ),
        }
    )


final = pd.DataFrame(
    final_rows
)

final.to_csv(
    OUT
    / "candidate_expression_validation_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 15. Console summary
# ============================================================

print()
print("=" * 72)
print("DATASET-LEVEL SUMMARY")
print("=" * 72)

print(
    dataset_summary.to_string(
        index=False
    )
)


print()
print("=" * 72)
print("PAIRWISE SIGNATURE CONCORDANCE")
print("=" * 72)

print(
    pairwise[
        [
            "candidate_identity",
            "dataset_a",
            "dataset_b",
            "n_informative_genes",
            "spearman_signature",
            "direction_concordance",
            "top_marker_overlap_n",
            "top_marker_jaccard",
        ]
    ].to_string(
        index=False
    )
)


print()
print("=" * 72)
print("FINAL EXPRESSION EVIDENCE")
print("=" * 72)

print(
    final.to_string(
        index=False
    )
)


print()
print("=" * 72)
print("OUTPUT")
print("=" * 72)
print(OUT)

print()
print("=" * 72)
print("SUCCESS")
print("=" * 72)

print(
    "Glial/Tuft expression-level "
    "cross-dataset validation completed."
)