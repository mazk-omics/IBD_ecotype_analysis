#!/usr/bin/env python3

from pathlib import Path
import importlib.metadata as md
import re
import sys

import pandas as pd


# =====================================================================
# SCP1884 Phase 0C/0D FINAL AUDIT
# =====================================================================

PROJECT = Path("/home/mazekai/IBD_EcoTyper")

SOURCE = Path(
    "/home/mazekai/enteric_glia/data_scRNA/SCP1884"
)

ENV_ROOT = Path(
    "/home/mazekai/miniconda3/envs/ibd_refaudit"
)

AUTHOR_FILE = (
    PROJECT
    / "03_reference/SCP1884/audit_v1/output/"
    "00_metadata_schema/author_metadata.parquet"
)

V2_FILE = (
    SOURCE
    / "metadata/scp_metadata_combined.v2.txt"
)

EXPR = (
    SOURCE
    / "expression"
)

# Fresh directory:
# do not mix final outputs with earlier failed diagnostics.
OUT = (
    PROJECT
    / "03_reference/SCP1884/audit_v1/output/"
    "01_barcode_sample_structure_final"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY = (
    OUT
    / "audit_summary.txt"
)


# =====================================================================
# Expected dataset structure
# =====================================================================

EXPECTED_CELLS = 720_633
EXPECTED_CHANNELS = 225
EXPECTED_SAMPLES = 136
EXPECTED_DONORS = 71

EXPECTED_LAYER_CHANNELS = {
    "E": 89,
    "L": 100,
    "N": 36,
}

EXPECTED_TYPE_CELLS = {
    "NonI": 418_708,
    "Heal": 181_806,
    "Infl": 120_119,
}

EXPECTED_SITE_CELLS = {
    "TI": 353_349,
    "CO": 289_730,
    "SB": 77_554,
}


# =====================================================================
# Expression bundles
# =====================================================================

BUNDLES = {

    "CO_STR": (
        EXPR / "Colon/Stromal",
        39_433,
        "Stromal cells",
        {"CO"},
    ),

    "CO_IMM": (
        EXPR / "Colon/Immune",
        152_509,
        "Immune cells",
        {"CO"},
    ),

    "CO_EPI": (
        EXPR / "Colon/Epithelial",
        97_788,
        "Epithelial cells",
        {"CO"},
    ),

    "TI_STR": (
        EXPR / "TerminalIleum/Stromal",
        75_695,
        "Stromal cells",
        {"TI", "SB"},
    ),

    "TI_IMM": (
        EXPR / "TerminalIleum/Immune",
        201_072,
        "Immune cells",
        {"TI", "SB"},
    ),

    "TI_EPI": (
        EXPR / "TerminalIleum/Epithelial",
        154_136,
        "Epithelial cells",
        {"TI", "SB"},
    ),
}


# =====================================================================
# Helpers
# =====================================================================

issues = []


def ok(msg):

    print(
        f"[PASS] {msg}"
    )


def bad(msg):

    issues.append(
        msg
    )

    print(
        f"[FAIL] {msg}"
    )


def pkg(name):

    try:
        return md.version(
            name
        )

    except md.PackageNotFoundError:
        return "NOT INSTALLED"


def count_check(
    series,
    expected,
    name,
):

    got = (
        series
        .value_counts(
            dropna=False
        )
        .to_dict()
    )

    if got == expected:

        ok(
            f"{name} counts match expected values."
        )

    else:

        bad(
            f"{name} counts mismatch: "
            f"observed={got}; "
            f"expected={expected}"
        )


def join_unique(x):

    return "+".join(
        sorted(
            set(
                map(
                    str,
                    x,
                )
            )
        )
    )


def parse_channel(channel):

    """
    Return:

        biological_sample_id
        layer_from_name
        naming_pattern

    Critical rule
    -------------
    Numeric / A-B suffixes identify BIOLOGICAL SAMPLES,
    not technical replicates.

    Therefore:

        N105598_E1
        N105598_L1
            -> same sample N105598_1

        N105598_E2
        N105598_L2
            -> same sample N105598_2

        N130064_N3
            -> unseparated sample N130064_3

        N10_Epi_A
        N10_LP_A
            -> same sample N10_A

        N105446_E
        N105446_L
            -> same sample N105446

        N104689_N
            -> sample N104689
    """

    channel = str(
        channel
    )

    # ---------------------------------------------------------
    # Legacy convention:
    #
    # N10_Epi_A
    # N10_LP_A
    # ---------------------------------------------------------

    m = re.fullmatch(
        r"(.+)_(Epi|LP)_([AB])",
        channel,
    )

    if m:

        base, tag, idx = (
            m.groups()
        )

        layer = {
            "Epi": "E",
            "LP": "L",
        }[tag]

        biological_sample_id = (
            f"{base}_{idx}"
        )

        return (
            biological_sample_id,
            layer,
            f"{tag}_{idx}",
        )

    # ---------------------------------------------------------
    # Numbered convention:
    #
    # N105598_E1
    # N105598_L1
    # N130064_N3
    # ---------------------------------------------------------

    m = re.fullmatch(
        r"(.+)_([ELN])([1-9][0-9]*)",
        channel,
    )

    if m:

        base, layer, idx = (
            m.groups()
        )

        biological_sample_id = (
            f"{base}_{idx}"
        )

        return (
            biological_sample_id,
            layer,
            f"{layer}{idx}",
        )

    # ---------------------------------------------------------
    # Simple convention:
    #
    # N105446_E
    # N105446_L
    # N104689_N
    # ---------------------------------------------------------

    m = re.fullmatch(
        r"(.+)_([ELN])",
        channel,
    )

    if m:

        base, layer = (
            m.groups()
        )

        return (
            base,
            layer,
            layer,
        )

    return (
        None,
        None,
        "UNPARSED",
    )


# =====================================================================
# 0. Environment
# =====================================================================

print(
    "=" * 78
)

print(
    "SCP1884 Phase 0C/0D FINAL AUDIT"
)

print(
    "=" * 78
)

python_exe = Path(
    sys.executable
).resolve()

if ENV_ROOT in python_exe.parents:

    ok(
        f"Correct environment: {python_exe}"
    )

else:

    bad(
        f"Wrong Python environment: {python_exe}"
    )

print(
    f"Python={sys.version.split()[0]} "
    f"pandas={pkg('pandas')} "
    f"pyarrow={pkg('pyarrow')}"
)


for f in (
    AUTHOR_FILE,
    V2_FILE,
):

    if not f.is_file():

        raise FileNotFoundError(
            f
        )


# =====================================================================
# 1. Author metadata
# =====================================================================

author = pd.read_parquet(
    AUTHOR_FILE
)

required = [

    "cell_id",
    "PubIDSample",
    "anno_overall",
    "n_genes",
    "n_counts",
    "Chem",
    "Site",
    "Type",
    "PubID",
    "Layer",
    "anno2",
]

missing = [

    x
    for x in required
    if x not in author.columns
]

if missing:

    raise RuntimeError(
        f"Missing author columns: {missing}"
    )


for c in [

    "cell_id",
    "PubIDSample",
    "PubID",
    "Site",
    "Type",
    "Layer",
    "Chem",
    "anno_overall",
    "anno2",
]:

    author[c] = (
        author[c]
        .astype(str)
    )


if len(author) == EXPECTED_CELLS:

    ok(
        "Author metadata contains "
        "720,633 cells."
    )

else:

    bad(
        f"Author metadata rows="
        f"{len(author):,}; "
        "expected=720,633."
    )


if author[
    "cell_id"
].is_unique:

    ok(
        "Author cell_id is unique."
    )

else:

    bad(
        "Author cell_id is not unique."
    )


count_check(
    author["Type"],
    EXPECTED_TYPE_CELLS,
    "Type",
)

count_check(
    author["Site"],
    EXPECTED_SITE_CELLS,
    "Site",
)


author[
    "expression_channel_id"
] = (
    author[
        "cell_id"
    ]
    .str.rsplit(
        "-",
        n=1,
    )
    .str[0]
)


# =====================================================================
# 2. expression_channel_id <-> PubIDSample
# =====================================================================

pairs = (
    author[
        [
            "expression_channel_id",
            "PubIDSample",
        ]
    ]
    .drop_duplicates()
)


n_expr = (
    author[
        "expression_channel_id"
    ]
    .nunique()
)

n_meta = (
    author[
        "PubIDSample"
    ]
    .nunique()
)


left_max = (
    pairs.groupby(
        "expression_channel_id"
    )[
        "PubIDSample"
    ]
    .nunique()
    .max()
)


right_max = (
    pairs.groupby(
        "PubIDSample"
    )[
        "expression_channel_id"
    ]
    .nunique()
    .max()
)


if (
    n_expr,
    n_meta,
    len(pairs),
    left_max,
    right_max,
) == (
    225,
    225,
    225,
    1,
    1,
):

    ok(
        "225 expression channels <-> "
        "225 PubIDSample is strict 1:1."
    )

else:

    bad(
        "Channel mapping unexpected: "
        f"expression={n_expr}, "
        f"metadata={n_meta}, "
        f"pairs={len(pairs)}, "
        f"left_max={left_max}, "
        f"right_max={right_max}"
    )


aliases = (
    pairs.loc[
        pairs[
            "expression_channel_id"
        ]
        !=
        pairs[
            "PubIDSample"
        ]
    ]
    .copy()
)


aliases.to_csv(
    OUT
    / "channel_name_aliases.tsv",
    sep="\t",
    index=False,
)


observed_aliases = set(
    map(
        tuple,
        aliases[
            [
                "expression_channel_id",
                "PubIDSample",
            ]
        ].to_numpy(),
    )
)


if observed_aliases == {
    (
        "N104689_N2",
        "N104689_N",
    )
}:

    ok(
        "Only known channel alias exists: "
        "N104689_N2 -> N104689_N."
    )

else:

    bad(
        "Unexpected channel aliases: "
        f"{sorted(observed_aliases)}"
    )


alias_cells = (
    author.loc[
        author[
            "expression_channel_id"
        ].eq(
            "N104689_N2"
        )
        &
        author[
            "PubIDSample"
        ].eq(
            "N104689_N"
        )
    ]
    .copy()
)


if (
    len(alias_cells) == 24_434
    and alias_cells[
        "PubID"
    ].eq(
        "104689"
    ).all()
    and alias_cells[
        "Site"
    ].eq(
        "CO"
    ).all()
    and alias_cells[
        "Type"
    ].eq(
        "NonI"
    ).all()
    and alias_cells[
        "Layer"
    ].eq(
        "N"
    ).all()
):

    ok(
        "Known channel alias contains "
        "exactly 24,434 internally consistent cells."
    )

else:

    bad(
        "Known channel alias validation failed; "
        f"cells={len(alias_cells):,}."
    )


# =====================================================================
# 3. Channel manifest and biological-sample reconstruction
# =====================================================================

channel = (
    author.groupby(
        "PubIDSample",
        observed=True,
    )
    .agg(

        expression_channel_id=(
            "expression_channel_id",
            "first",
        ),

        donor_id=(
            "PubID",
            "first",
        ),

        site_original=(
            "Site",
            "first",
        ),

        type_original=(
            "Type",
            "first",
        ),

        layer_metadata=(
            "Layer",
            "first",
        ),

        chemistry=(
            "Chem",
            "first",
        ),

        n_cells=(
            "cell_id",
            "size",
        ),

        n_expression_channel=(
            "expression_channel_id",
            "nunique",
        ),

        n_donor=(
            "PubID",
            "nunique",
        ),

        n_site=(
            "Site",
            "nunique",
        ),

        n_type=(
            "Type",
            "nunique",
        ),

        n_layer=(
            "Layer",
            "nunique",
        ),
    )
    .reset_index()
    .rename(
        columns={
            "PubIDSample":
                "metadata_channel_id"
        }
    )
)


internal_cols = [

    "n_expression_channel",
    "n_donor",
    "n_site",
    "n_type",
    "n_layer",
]


for col in internal_cols:

    x = (
        channel.loc[
            channel[col] != 1
        ]
    )

    if not x.empty:

        x.to_csv(
            OUT
            / f"ERROR_channel_internal_{col}.tsv",
            sep="\t",
            index=False,
        )

        bad(
            f"{len(x)} channels are internally "
            f"inconsistent for {col}."
        )


parsed = (
    channel[
        "metadata_channel_id"
    ]
    .map(
        parse_channel
    )
)


channel[
    "biological_sample_id"
] = [
    x[0]
    for x in parsed
]


channel[
    "layer_from_name"
] = [
    x[1]
    for x in parsed
]


channel[
    "channel_pattern"
] = [
    x[2]
    for x in parsed
]


unparsed = (
    channel.loc[
        channel[
            "channel_pattern"
        ].eq(
            "UNPARSED"
        )
    ]
)


if unparsed.empty:

    ok(
        "All 225 channel names match a "
        "recognized author naming convention."
    )

else:

    unparsed.to_csv(
        OUT
        / "ERROR_unparsed_channels.tsv",
        sep="\t",
        index=False,
    )

    bad(
        f"{len(unparsed)} channels "
        "could not be parsed."
    )


layer_bad = (
    channel.loc[
        channel[
            "layer_from_name"
        ].notna()
        &
        channel[
            "layer_from_name"
        ].ne(
            channel[
                "layer_metadata"
            ]
        )
    ]
)


if layer_bad.empty:

    ok(
        "Name-derived layer agrees with "
        "metadata Layer for all 225 channels."
    )

else:

    layer_bad.to_csv(
        OUT
        / "ERROR_channel_name_vs_layer.tsv",
        sep="\t",
        index=False,
    )

    bad(
        f"{len(layer_bad)} channels "
        "disagree on layer."
    )


count_check(
    channel[
        "layer_metadata"
    ],
    EXPECTED_LAYER_CHANNELS,
    "Channel Layer",
)


pattern_counts = (
    channel.groupby(
        [
            "layer_metadata",
            "channel_pattern",
        ],
        observed=True,
    )
    .size()
    .rename(
        "n_channels"
    )
    .reset_index()
    .sort_values(
        [
            "layer_metadata",
            "channel_pattern",
        ]
    )
)


pattern_counts.to_csv(
    OUT
    / "channel_naming_patterns.tsv",
    sep="\t",
    index=False,
)


# A biological sample may have E + L,
# but never >1 channel for the same layer.

dup = (
    channel
    .duplicated(
        [
            "biological_sample_id",
            "layer_metadata",
        ],
        keep=False,
    )
)


if dup.any():

    channel.loc[
        dup
    ].to_csv(
        OUT
        / "ERROR_duplicate_sample_layer.tsv",
        sep="\t",
        index=False,
    )

    bad(
        f"{int(dup.sum())} channel rows duplicate "
        "biological_sample_id + layer."
    )

else:

    ok(
        "At most one channel per "
        "biological sample/layer."
    )


# Biological metadata must agree between
# E and L channels belonging to one sample.

for field in [

    "donor_id",
    "site_original",
    "type_original",
]:

    nuniq = (
        channel.groupby(
            "biological_sample_id",
            observed=True,
        )[field]
        .nunique(
            dropna=False
        )
    )

    bad_ids = (
        nuniq.loc[
            nuniq > 1
        ]
        .index
    )

    if len(
        bad_ids
    ):

        channel.loc[
            channel[
                "biological_sample_id"
            ].isin(
                bad_ids
            )
        ].to_csv(
            OUT
            / f"ERROR_sample_inconsistent_{field}.tsv",
            sep="\t",
            index=False,
        )

        bad(
            f"{len(bad_ids)} biological samples "
            f"contain multiple {field} values."
        )


sample_manifest = (
    channel.groupby(
        "biological_sample_id",
        observed=True,
    )
    .agg(

        donor_id=(
            "donor_id",
            "first",
        ),

        site_original=(
            "site_original",
            "first",
        ),

        type_original=(
            "type_original",
            "first",
        ),

        n_channels=(
            "metadata_channel_id",
            "nunique",
        ),

        layers=(
            "layer_metadata",
            join_unique,
        ),

        channel_patterns=(
            "channel_pattern",
            join_unique,
        ),

        n_cells=(
            "n_cells",
            "sum",
        ),
    )
    .reset_index()
)


if len(
    sample_manifest
) == EXPECTED_SAMPLES:

    ok(
        "Exactly 136 biological samples reconstructed."
    )

else:

    bad(
        f"Biological samples="
        f"{len(sample_manifest)}; "
        "expected=136."
    )


if (
    sample_manifest[
        "n_cells"
    ]
    .sum()
    == EXPECTED_CELLS
):

    ok(
        "Biological sample manifest accounts "
        "for all 720,633 cells."
    )

else:

    bad(
        "Biological sample manifest does not "
        "account for all cells."
    )


# ---------------------------------------------------------
# Publication-level sample structure:
#
# 136 biological samples =
#   36 unseparated N samples
# + 100 separated samples
#
# separated samples:
#   89 have E
#  100 have L
# ---------------------------------------------------------

invalid_mix = (
    sample_manifest.loc[
        sample_manifest[
            "layers"
        ].str.contains(
            "N",
            na=False,
        )
        &
        sample_manifest[
            "layers"
        ].str.contains(
            r"E|L",
            regex=True,
            na=False,
        )
    ]
)


if not invalid_mix.empty:

    invalid_mix.to_csv(
        OUT
        / "ERROR_samples_mix_N_with_EL.tsv",
        sep="\t",
        index=False,
    )

    bad(
        f"{len(invalid_mix)} samples "
        "mix N with E/L."
    )


n_n = int(
    sample_manifest[
        "layers"
    ]
    .eq(
        "N"
    )
    .sum()
)


n_sep = int(
    len(
        sample_manifest
    )
    - n_n
)


n_e = int(
    channel.loc[
        channel[
            "layer_metadata"
        ].eq(
            "E"
        ),
        "biological_sample_id",
    ]
    .nunique()
)


n_l = int(
    channel.loc[
        channel[
            "layer_metadata"
        ].eq(
            "L"
        ),
        "biological_sample_id",
    ]
    .nunique()
)


if (
    n_n,
    n_sep,
    n_e,
    n_l,
) == (
    36,
    100,
    89,
    100,
):

    ok(
        "Sample structure closes: "
        "36 unseparated + 100 separated; "
        "89 E + 100 L."
    )

else:

    bad(
        "Sample structure mismatch: "
        f"N={n_n}, "
        f"separated={n_sep}, "
        f"E={n_e}, "
        f"L={n_l}."
    )


# =====================================================================
# 4. Portal V2
# =====================================================================

v2 = pd.read_csv(
    V2_FILE,
    sep="\t",
    dtype=str,
    low_memory=False,
)


need_v2 = [

    "NAME",
    "donor_id",
    "disease",
    "disease__ontology_label",
    "Celltype",
]


missing = [

    x
    for x in need_v2
    if x not in v2.columns
]


if missing:

    raise RuntimeError(
        f"Missing V2 columns: {missing}"
    )


type_row = (
    v2[
        "NAME"
    ]
    .astype(str)
    .eq(
        "TYPE"
    )
)


if int(
    type_row.sum()
) != 1:

    bad(
        "Portal V2 TYPE rows="
        f"{int(type_row.sum())}; "
        "expected=1."
    )


v2 = (
    v2.loc[
        ~type_row,
        need_v2,
    ]
    .copy()
    .rename(
        columns={

            "NAME":
                "cell_id",

            "donor_id":
                "v2_donor_id",

            "disease":
                "disease_ontology_id",

            "disease__ontology_label":
                "disease_label",

            "Celltype":
                "v2_celltype",
        }
    )
)


if (
    len(v2) == EXPECTED_CELLS
    and v2[
        "cell_id"
    ].is_unique
):

    ok(
        "Portal V2 contains "
        "720,633 unique real cells."
    )

else:

    bad(
        f"Portal V2 rows={len(v2):,}; "
        "unique IDs="
        f"{v2['cell_id'].nunique():,}."
    )


# =====================================================================
# 5. Six expression bundles + barcode closure
# =====================================================================

author_idx = (
    author
    .set_index(
        "cell_id",
        drop=False,
    )
)


barcode_frames = []
bundle_rows = []


for (
    bundle,
    (
        directory,
        expected_n,
        expected_compartment,
        allowed_sites,
    ),
) in BUNDLES.items():

    files = list(
        directory.glob(
            "*.barcodes.tsv"
        )
    )

    if len(
        files
    ) != 1:

        bad(
            f"{bundle}: expected 1 barcode file; "
            f"found {len(files)}."
        )

        continue


    barcodes = (
        pd.read_csv(
            files[0],
            sep="\t",
            header=None,
            dtype=str,
        )
        .iloc[:, 0]
        .astype(str)
    )


    if len(
        barcodes
    ) != expected_n:

        bad(
            f"{bundle}: rows="
            f"{len(barcodes):,}; "
            f"expected={expected_n:,}."
        )


    if not barcodes.is_unique:

        bad(
            f"{bundle}: duplicate barcode IDs."
        )


    # IMPORTANT:
    # register each bundle exactly once,
    # before any annotation-exception logic.

    barcode_frames.append(
        pd.DataFrame(
            {
                "cell_id":
                    barcodes.to_numpy(),

                "expression_bundle":
                    bundle,
            }
        )
    )


    missing_ids = (
        ~barcodes.isin(
            author_idx.index
        )
    )


    if missing_ids.any():

        pd.DataFrame(
            {
                "cell_id":
                    barcodes.loc[
                        missing_ids
                    ]
            }
        ).to_csv(
            OUT
            / f"ERROR_{bundle}_not_in_metadata.tsv",
            sep="\t",
            index=False,
        )

        bad(
            f"{bundle}: "
            f"{int(missing_ids.sum()):,} "
            "IDs absent from metadata."
        )

        continue


    sub = (
        author_idx.loc[
            barcodes.tolist()
        ]
        .copy()
    )


    sites = set(
        sub[
            "Site"
        ]
        .unique()
    )


    if not sites.issubset(
        allowed_sites
    ):

        bad(
            f"{bundle}: "
            f"unexpected sites="
            f"{sorted(sites)}."
        )


    mismatch = (
        sub.loc[
            sub[
                "anno_overall"
            ].ne(
                expected_compartment
            )
        ]
        .copy()
    )


    known_exception = 0


    if bundle == "TI_IMM":

        valid = (
            len(
                mismatch
            ) == 376
            and mismatch[
                "anno_overall"
            ].eq(
                "Epithelial cells"
            ).all()
            and mismatch[
                "anno2"
            ].eq(
                "Epithelial cells HBB+ HBA+"
            ).all()
        )


        if valid:

            known_exception = 376

            ok(
                "TI_IMM contains exactly the "
                "known 376-cell author "
                "annotation exception."
            )

        else:

            mismatch.to_csv(
                OUT
                / "ERROR_TI_IMM_author_exception.tsv",
                sep="\t",
                index=False,
            )

            bad(
                "TI_IMM unexpected compartment "
                f"mismatches={len(mismatch):,}."
            )


    elif not mismatch.empty:

        mismatch.to_csv(
            OUT
            / f"ERROR_{bundle}_compartment.tsv",
            sep="\t",
            index=False,
        )

        bad(
            f"{bundle}: unexpected compartment "
            f"mismatches={len(mismatch):,}."
        )


    bundle_rows.append(
        {

            "bundle":
                bundle,

            "barcode_file":
                str(
                    files[0]
                ),

            "n_cells":
                len(
                    barcodes
                ),

            "n_expected_compartment":
                int(
                    sub[
                        "anno_overall"
                    ]
                    .eq(
                        expected_compartment
                    )
                    .sum()
                ),

            "n_known_exception":
                known_exception,

            "observed_sites":
                "|".join(
                    sorted(
                        sites
                    )
                ),
        }
    )


if (
    len(
        barcode_frames
    ) == 6
    and sum(
        map(
            len,
            barcode_frames,
        )
    )
    == EXPECTED_CELLS
):

    ok(
        "Six expression barcode bundles "
        "registered exactly once; "
        "total=720,633."
    )

else:

    bad(
        "Barcode registration: "
        f"bundles={len(barcode_frames)}, "
        "rows="
        f"{sum(map(len, barcode_frames)):,}."
    )


barcode_map = pd.concat(
    barcode_frames,
    ignore_index=True,
)


if barcode_map[
    "cell_id"
].duplicated().any():

    barcode_map.loc[
        barcode_map[
            "cell_id"
        ].duplicated(
            keep=False
        )
    ].to_csv(
        OUT
        / "ERROR_duplicate_barcodes_across_bundles.tsv",
        sep="\t",
        index=False,
    )

    bad(
        "Barcode IDs overlap across "
        "expression bundles."
    )


meta_ids = set(
    author[
        "cell_id"
    ]
)

expr_ids = set(
    barcode_map[
        "cell_id"
    ]
)

metadata_only = (
    meta_ids
    - expr_ids
)

expression_only = (
    expr_ids
    - meta_ids
)


if (
    not metadata_only
    and not expression_only
    and len(
        barcode_map
    ) == EXPECTED_CELLS
):

    ok(
        "720,633 expression barcodes exactly equal "
        "720,633 metadata cell IDs."
    )

else:

    bad(
        "Cell-ID closure failed: "
        f"metadata_only={len(metadata_only)}, "
        f"expression_only={len(expression_only)}, "
        f"rows={len(barcode_map):,}."
    )


pd.DataFrame(
    bundle_rows
).to_csv(
    OUT
    / "expression_bundle_audit.tsv",
    sep="\t",
    index=False,
)


# =====================================================================
# 6. Cell-level provenance + known 376-cell V2 correction
# =====================================================================

mapped = (
    author
    .merge(
        barcode_map,
        on="cell_id",
        how="left",
        validate="one_to_one",
    )
    .merge(
        v2,
        on="cell_id",
        how="left",
        validate="one_to_one",
    )
)


channel_lookup = (
    channel[
        [
            "metadata_channel_id",
            "biological_sample_id",
            "layer_from_name",
            "channel_pattern",
        ]
    ]
    .rename(
        columns={
            "metadata_channel_id":
                "PubIDSample"
        }
    )
)


mapped = mapped.merge(
    channel_lookup,
    on="PubIDSample",
    how="left",
    validate="many_to_one",
)


if mapped[
    "biological_sample_id"
].notna().all():

    ok(
        "Every cell maps to a reconstructed "
        "biological sample."
    )

else:

    bad(
        "Some cells lack biological_sample_id."
    )


if (
    mapped[
        "PubID"
    ]
    .astype(str)
    .eq(
        mapped[
            "v2_donor_id"
        ]
        .astype(str)
    )
    .all()
):

    ok(
        "Author PubID and Portal V2 donor_id "
        "agree for all cells."
    )

else:

    bad(
        "Author PubID and Portal V2 donor_id "
        "disagree for some cells."
    )


type_disease = (
    mapped.groupby(
        [
            "Type",
            "disease_label",
        ],
        observed=True,
    )
    .size()
    .rename(
        "n_cells"
    )
    .reset_index()
)


type_disease.to_csv(
    OUT
    / "type_vs_disease.tsv",
    sep="\t",
    index=False,
)


expected_pairs = {

    (
        "Heal",
        "normal",
    ),

    (
        "NonI",
        "Crohn's disease",
    ),

    (
        "Infl",
        "Crohn's disease",
    ),
}


observed_pairs = set(
    zip(
        type_disease[
            "Type"
        ],
        type_disease[
            "disease_label"
        ],
    )
)


if observed_pairs == expected_pairs:

    ok(
        "Type <-> Portal V2 disease "
        "semantics agree exactly."
    )

else:

    bad(
        "Unexpected Type/disease pairs: "
        f"{sorted(observed_pairs)}"
    )


mapped[
    "disease"
] = mapped[
    "Type"
].map(
    {
        "Heal": "Healthy",
        "NonI": "CD",
        "Infl": "CD",
    }
)


mapped[
    "inflammation"
] = mapped[
    "Type"
].map(
    {
        "Heal": "Healthy",
        "NonI": "Non_Inflamed",
        "Infl": "Inflamed",
    }
)


mapped[
    "site_canonical"
] = mapped[
    "Site"
].map(
    {
        "CO": "Colon",
        "TI": "Terminal_Ileum",
        "SB": "Small_Bowel",
    }
)


mapped[
    "author_analysis_region"
] = mapped[
    "Site"
].map(
    {
        "CO": "Colon",
        "TI": "TI_plus_SB",
        "SB": "TI_plus_SB",
    }
)


ti_exception = (
    mapped.loc[
        mapped[
            "expression_bundle"
        ].eq(
            "TI_IMM"
        )
        &
        mapped[
            "anno_overall"
        ].eq(
            "Epithelial cells"
        )
    ]
    .copy()
)


if (
    len(
        ti_exception
    ) == 376
    and ti_exception[
        "anno2"
    ].eq(
        "Epithelial cells HBB+ HBA+"
    ).all()
    and ti_exception[
        "v2_celltype"
    ].eq(
        "Monocytes HBB"
    ).all()
):

    ok(
        "All 376 TI_IMM exception cells "
        "are Portal V2 'Monocytes HBB'."
    )

else:

    bad(
        "TI_IMM V2 exception validation failed; "
        f"cells={len(ti_exception):,}."
    )


ti_exception.to_csv(
    OUT
    / "TI_IMM_376_cell_annotation_exception.tsv",
    sep="\t",
    index=False,
)


# =====================================================================
# 7. Final channel / sample / donor manifests
# =====================================================================

semantics = (
    mapped.groupby(
        "PubIDSample",
        observed=True,
    )
    .agg(

        disease=(
            "disease",
            "first",
        ),

        inflammation=(
            "inflammation",
            "first",
        ),

        site_canonical=(
            "site_canonical",
            "first",
        ),

        author_analysis_region=(
            "author_analysis_region",
            "first",
        ),
    )
    .reset_index()
    .rename(
        columns={
            "PubIDSample":
                "metadata_channel_id"
        }
    )
)


channel = channel.merge(
    semantics,
    on="metadata_channel_id",
    how="left",
    validate="one_to_one",
)


channel[
    "channel_name_status"
] = "MATCH"


channel.loc[
    channel[
        "expression_channel_id"
    ].eq(
        "N104689_N2"
    )
    &
    channel[
        "metadata_channel_id"
    ].eq(
        "N104689_N"
    ),
    "channel_name_status",
] = (
    "KNOWN_N104689_N2_ALIAS"
)


# Recheck biological consistency after
# adding canonical fields.

for field in [

    "donor_id",
    "disease",
    "inflammation",
    "site_original",
    "site_canonical",
    "author_analysis_region",
]:

    nuniq = (
        channel.groupby(
            "biological_sample_id",
            observed=True,
        )[field]
        .nunique(
            dropna=False
        )
    )

    bad_ids = (
        nuniq.loc[
            nuniq > 1
        ]
        .index
    )

    if len(
        bad_ids
    ):

        channel.loc[
            channel[
                "biological_sample_id"
            ].isin(
                bad_ids
            )
        ].to_csv(
            OUT
            / f"ERROR_final_sample_inconsistent_{field}.tsv",
            sep="\t",
            index=False,
        )

        bad(
            f"{len(bad_ids)} biological samples "
            f"inconsistent for {field}."
        )


sample_manifest = (
    channel.groupby(
        "biological_sample_id",
        observed=True,
    )
    .agg(

        donor_id=(
            "donor_id",
            "first",
        ),

        disease=(
            "disease",
            "first",
        ),

        inflammation=(
            "inflammation",
            "first",
        ),

        site_original=(
            "site_original",
            "first",
        ),

        site_canonical=(
            "site_canonical",
            "first",
        ),

        author_analysis_region=(
            "author_analysis_region",
            "first",
        ),

        n_channels=(
            "metadata_channel_id",
            "nunique",
        ),

        layers=(
            "layer_metadata",
            join_unique,
        ),

        channel_patterns=(
            "channel_pattern",
            join_unique,
        ),

        n_cells=(
            "n_cells",
            "sum",
        ),
    )
    .reset_index()
)


donor_manifest = (
    sample_manifest.groupby(
        "donor_id",
        observed=True,
    )
    .agg(

        disease=(
            "disease",
            lambda x:
            x.iloc[0]
            if x.nunique() == 1
            else "INCONSISTENT",
        ),

        n_samples=(
            "biological_sample_id",
            "nunique",
        ),

        n_cells=(
            "n_cells",
            "sum",
        ),
    )
    .reset_index()
)


if len(
    sample_manifest
) == 136:

    ok(
        "Final biological-sample manifest "
        "contains 136 samples."
    )

else:

    bad(
        "Final biological-sample manifest "
        f"contains {len(sample_manifest)} samples."
    )


n_cd = int(
    donor_manifest[
        "disease"
    ]
    .eq(
        "CD"
    )
    .sum()
)


n_h = int(
    donor_manifest[
        "disease"
    ]
    .eq(
        "Healthy"
    )
    .sum()
)


if (
    len(
        donor_manifest
    ) == 71
    and n_cd == 46
    and n_h == 25
    and not donor_manifest[
        "disease"
    ].eq(
        "INCONSISTENT"
    ).any()
):

    ok(
        "Donors close exactly: "
        "71 total = 46 CD + 25 Healthy."
    )

else:

    bad(
        "Donor structure: "
        f"total={len(donor_manifest)}, "
        f"CD={n_cd}, "
        f"Healthy={n_h}."
    )


sampling = (
    sample_manifest.groupby(
        [
            "disease",
            "inflammation",
            "site_original",
        ],
        observed=True,
    )
    .agg(

        n_samples=(
            "biological_sample_id",
            "nunique",
        ),

        n_donors=(
            "donor_id",
            "nunique",
        ),

        n_cells=(
            "n_cells",
            "sum",
        ),
    )
    .reset_index()
)


# =====================================================================
# 8. Save outputs
# =====================================================================

channel.to_csv(
    OUT
    / "channel_manifest.tsv",
    sep="\t",
    index=False,
)


sample_manifest.to_csv(
    OUT
    / "biological_sample_manifest.tsv",
    sep="\t",
    index=False,
)


donor_manifest.to_csv(
    OUT
    / "donor_manifest.tsv",
    sep="\t",
    index=False,
)


sampling.to_csv(
    OUT
    / "sampling_structure.tsv",
    sep="\t",
    index=False,
)


mapped.to_parquet(
    OUT
    / "scp1884_cells_with_final_source_structure.parquet",
    index=False,
    engine="pyarrow",
)


# =====================================================================
# 9. Final report
#
# Only one final exit point.
# =====================================================================

report = [

    "SCP1884 Phase 0C/0D FINAL AUDIT",

    "=" * 78,

    f"Cells: {len(mapped):,}",

    f"Expression channels: {n_expr}",

    f"Metadata channels: {n_meta}",

    (
        "Biological samples: "
        f"{len(sample_manifest)}"
    ),

    f"Donors: {len(donor_manifest)}",

    (
        "CD / Healthy donors: "
        f"{n_cd} / {n_h}"
    ),

    (
        "Known channel alias cells: "
        f"{len(alias_cells):,}"
    ),

    (
        "TI_IMM annotation exception cells: "
        f"{len(ti_exception):,}"
    ),

    "",

    "Channel naming patterns:",

    pattern_counts.to_string(
        index=False
    ),

    "",

    "Sampling structure:",

    sampling.to_string(
        index=False
    ),

    "",
]


if issues:

    report.append(
        f"FINAL STATUS: FAIL "
        f"({len(issues)} issue(s))"
    )

    report += [
        f"  - {x}"
        for x in issues
    ]


else:

    report += [

        "FINAL STATUS: PASS",

        (
            "  - 720,633 cell IDs "
            "closed exactly"
        ),

        (
            "  - 225 expression channels "
            "<-> 225 PubIDSample strict 1:1"
        ),

        (
            "  - all channel naming "
            "conventions parsed"
        ),

        (
            "  - 136 biological samples "
            "reconstructed"
        ),

        (
            "  - 36 unseparated + "
            "100 separated samples confirmed"
        ),

        (
            "  - 89 E + 100 L + "
            "36 N channels confirmed"
        ),

        (
            "  - 71 donors = "
            "46 CD + 25 Healthy"
        ),

        (
            "  - both known exceptions validated"
        ),
    ]


SUMMARY.write_text(
    "\n".join(
        report
    )
    + "\n",
    encoding="utf-8",
)


print()

print(
    "=" * 78
)

print(
    "FINAL RESULTS"
)

print(
    "=" * 78
)

print(
    f"Cells               : "
    f"{len(mapped):,}"
)

print(
    f"Expression channels : "
    f"{n_expr}"
)

print(
    f"Metadata channels   : "
    f"{n_meta}"
)

print(
    f"Biological samples  : "
    f"{len(sample_manifest)}"
)

print(
    f"Donors              : "
    f"{len(donor_manifest)}"
)

print(
    f"CD / Healthy donors : "
    f"{n_cd} / {n_h}"
)

print(
    f"Known alias cells   : "
    f"{len(alias_cells):,}"
)

print(
    f"TI_IMM exception    : "
    f"{len(ti_exception):,}"
)

print()

print(
    "Channel naming patterns:"
)

print(
    pattern_counts.to_string(
        index=False
    )
)

print()

print(
    "Sampling structure:"
)

print(
    sampling.to_string(
        index=False
    )
)

print()


if issues:

    print(
        f"FINAL STATUS: FAIL "
        f"({len(issues)} issue(s))"
    )

    for x in issues:

        print(
            f"  - {x}"
        )

    print()

    print(
        f"Full report: {SUMMARY}"
    )

    sys.exit(
        1
    )


print(
    "FINAL STATUS: PASS"
)

print(
    f"Full report: {SUMMARY}"
)

print(
    f"Final outputs: {OUT}"
)
