#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd


PROJECT = Path("/home/mazekai/IBD_EcoTyper")

SRC = (
    PROJECT
    / "03_reference/SCP1884/audit_v1/output/"
    "01_barcode_sample_structure_final"
)

OUT = (
    PROJECT
    / "03_reference/SCP1884/audit_v1/output/"
    "01_barcode_sample_structure_FINAL"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

CHANNEL_FILE = (
    SRC / "channel_manifest.tsv"
)

CELL_FILE = (
    SRC
    / "scp1884_cells_with_final_source_structure.parquet"
)


EXPECTED_CELLS = 720_633
EXPECTED_CHANNELS = 225
EXPECTED_SAMPLES = 136
EXPECTED_DONORS = 71

EXPECTED_LAYER_CHANNELS = {
    "E": 89,
    "L": 100,
    "N": 36,
}


def fail(msg):
    raise RuntimeError(msg)


def join_unique(x):
    return "+".join(
        sorted(
            set(
                x.astype(str)
            )
        )
    )


# ============================================================
# 1. Load finalized channel-level results from previous audit
# ============================================================

channel = pd.read_csv(
    CHANNEL_FILE,
    sep="\t",
    dtype=str,
)

channel["n_cells"] = (
    channel["n_cells"]
    .astype(int)
)

if len(channel) != EXPECTED_CHANNELS:
    fail(
        f"Expected 225 channels; found {len(channel)}."
    )

print(
    "[PASS] Loaded 225 validated channels."
)


# ============================================================
# 2. Diagnose the two unmatched E samples BEFORE correction
# ============================================================

initial = (
    channel.groupby(
        "biological_sample_id",
        observed=True,
    )
    .agg(
        layers=(
            "layer_metadata",
            join_unique,
        ),
        donor_id=(
            "donor_id",
            "first",
        ),
        site=(
            "site_original",
            "first",
        ),
        type=(
            "type_original",
            "first",
        ),
    )
    .reset_index()
)

e_only = set(
    initial.loc[
        initial["layers"] == "E",
        "biological_sample_id",
    ]
)

l_only = set(
    initial.loc[
        initial["layers"] == "L",
        "biological_sample_id",
    ]
)

expected_e_only = {
    "N127693",
    "I127693",
}

expected_partner_l = {
    "N127693_1",
    "I127693_2",
}

if e_only != expected_e_only:

    print(
        "Observed E-only sample IDs:",
        sorted(e_only),
    )

    fail(
        "Unexpected E-only sample set. "
        "Do not apply pairing correction."
    )

if not expected_partner_l.issubset(
    l_only
):
    fail(
        "Expected N127693_1 and I127693_2 "
        "L-only samples were not found."
    )

print(
    "[PASS] Exactly two unmatched E samples identified:"
)

print(
    "       N127693_E and I127693_E"
)


# ============================================================
# 3. Preserve initial reconstruction
# ============================================================

channel[
    "biological_sample_id_initial"
] = (
    channel[
        "biological_sample_id"
    ]
)

channel[
    "sample_pairing_status"
] = "STANDARD"


# ============================================================
# 4. Correct the two author naming exceptions
#
# N127693_E  + N127693_L1 -> N127693_1
# I127693_E  + I127693_L2 -> I127693_2
# ============================================================

special_map = {

    "N127693_E":
        "N127693_1",

    "N127693_L1":
        "N127693_1",

    "I127693_E":
        "I127693_2",

    "I127693_L2":
        "I127693_2",
}


for channel_id, sample_id in special_map.items():

    mask = (
        channel[
            "metadata_channel_id"
        ]
        .eq(
            channel_id
        )
    )

    if int(mask.sum()) != 1:
        fail(
            f"{channel_id}: expected one channel row; "
            f"found {int(mask.sum())}."
        )

    channel.loc[
        mask,
        "biological_sample_id",
    ] = sample_id

    channel.loc[
        mask,
        "sample_pairing_status",
    ] = (
        "KNOWN_127693_E_L_PAIRING_EXCEPTION"
    )


# ============================================================
# 5. Validate each corrected biological sample
# ============================================================

# No biological sample may contain two E or two L channels.

dup = (
    channel.duplicated(
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

    fail(
        "Some corrected samples contain "
        "multiple channels from the same layer."
    )


# donor/site/type/disease/inflammation must be identical
# within a biological sample.

consistency_fields = [
    "donor_id",
    "site_original",
    "type_original",
    "disease",
    "inflammation",
]

for field in consistency_fields:

    x = (
        channel.groupby(
            "biological_sample_id",
            observed=True,
        )[field]
        .nunique(
            dropna=False
        )
    )

    bad = x[
        x > 1
    ]

    if not bad.empty:

        channel.loc[
            channel[
                "biological_sample_id"
            ].isin(
                bad.index
            )
        ].to_csv(
            OUT
            / f"ERROR_inconsistent_{field}.tsv",
            sep="\t",
            index=False,
        )

        fail(
            f"{len(bad)} samples inconsistent "
            f"for {field}."
        )


print(
    "[PASS] Corrected E/L pairs are biologically consistent."
)


# ============================================================
# 6. Build final biological-sample manifest
# ============================================================

sample = (
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


if len(sample) != EXPECTED_SAMPLES:

    fail(
        f"Final samples={len(sample)}; "
        "expected=136."
    )


if sample["n_cells"].sum() != EXPECTED_CELLS:

    fail(
        "Sample manifest does not sum to 720,633 cells."
    )


print(
    "[PASS] Exactly 136 biological samples reconstructed."
)


# ============================================================
# 7. Validate 36 unseparated + 100 separated
# ============================================================

n_unseparated = int(
    sample[
        "layers"
    ]
    .eq(
        "N"
    )
    .sum()
)

n_separated = (
    len(sample)
    - n_unseparated
)

n_e = (
    channel.loc[
        channel[
            "layer_metadata"
        ] == "E",
        "biological_sample_id",
    ]
    .nunique()
)

n_l = (
    channel.loc[
        channel[
            "layer_metadata"
        ] == "L",
        "biological_sample_id",
    ]
    .nunique()
)

layer_counts = (
    channel[
        "layer_metadata"
    ]
    .value_counts()
    .to_dict()
)

if layer_counts != EXPECTED_LAYER_CHANNELS:

    fail(
        f"Layer channel counts mismatch: {layer_counts}"
    )

if (
    n_unseparated,
    n_separated,
    n_e,
    n_l,
) != (
    36,
    100,
    89,
    100,
):

    fail(
        "Final sample structure mismatch: "
        f"N={n_unseparated}, "
        f"separated={n_separated}, "
        f"E={n_e}, "
        f"L={n_l}."
    )


print(
    "[PASS] Sample structure closes exactly:"
)

print(
    "       36 unseparated + 100 separated"
)

print(
    "       89 E + 100 L + 36 N channels"
)


# ============================================================
# 8. Validate donor structure
# ============================================================

donor_n_disease = (
    sample.groupby(
        "donor_id",
        observed=True,
    )[
        "disease"
    ]
    .nunique()
)

if (
    donor_n_disease > 1
).any():

    fail(
        "A donor maps to multiple disease states."
    )


donor = (
    sample.groupby(
        "donor_id",
        observed=True,
    )
    .agg(

        disease=(
            "disease",
            "first",
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


n_cd = int(
    donor[
        "disease"
    ]
    .eq(
        "CD"
    )
    .sum()
)

n_healthy = int(
    donor[
        "disease"
    ]
    .eq(
        "Healthy"
    )
    .sum()
)


if (
    len(donor),
    n_cd,
    n_healthy,
) != (
    71,
    46,
    25,
):

    fail(
        "Donor structure mismatch: "
        f"total={len(donor)}, "
        f"CD={n_cd}, "
        f"Healthy={n_healthy}."
    )


print(
    "[PASS] 71 donors = 46 CD + 25 Healthy."
)


# ============================================================
# 9. Correct cell-level biological_sample_id
# ============================================================

cells = pd.read_parquet(
    CELL_FILE
)

if len(cells) != EXPECTED_CELLS:
    fail(
        f"Cell provenance rows={len(cells):,}; "
        "expected=720,633."
    )


cells[
    "biological_sample_id_initial"
] = (
    cells[
        "biological_sample_id"
    ]
)


cells[
    "sample_pairing_status"
] = "STANDARD"


for channel_id, sample_id in special_map.items():

    mask = (
        cells[
            "PubIDSample"
        ]
        .astype(str)
        .eq(
            channel_id
        )
    )

    if not mask.any():
        fail(
            f"No cells found for {channel_id}."
        )

    cells.loc[
        mask,
        "biological_sample_id",
    ] = sample_id

    cells.loc[
        mask,
        "sample_pairing_status",
    ] = (
        "KNOWN_127693_E_L_PAIRING_EXCEPTION"
    )


if (
    cells[
        "biological_sample_id"
    ]
    .nunique()
    != EXPECTED_SAMPLES
):

    fail(
        "Corrected cell-level data does not "
        "contain 136 biological samples."
    )


# ============================================================
# 10. Sampling table
# ============================================================

sampling = (
    sample.groupby(
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


# ============================================================
# 11. Save FINAL outputs
# ============================================================

channel.to_csv(
    OUT
    / "channel_manifest_FINAL.tsv",
    sep="\t",
    index=False,
)

sample.to_csv(
    OUT
    / "biological_sample_manifest_FINAL.tsv",
    sep="\t",
    index=False,
)

donor.to_csv(
    OUT
    / "donor_manifest_FINAL.tsv",
    sep="\t",
    index=False,
)

sampling.to_csv(
    OUT
    / "sampling_structure_FINAL.tsv",
    sep="\t",
    index=False,
)

cells.to_parquet(
    OUT
    / "scp1884_cells_with_source_structure_FINAL.parquet",
    index=False,
    engine="pyarrow",
)


# ============================================================
# 12. Final report
# ============================================================

report = [
    "SCP1884 PHASE 0 FINAL",
    "=" * 70,
    "",
    f"Cells: {len(cells):,}",
    f"Channels: {len(channel)}",
    f"Biological samples: {len(sample)}",
    f"Donors: {len(donor)}",
    f"CD donors: {n_cd}",
    f"Healthy donors: {n_healthy}",
    "",
    "Sample structure:",
    f"  Unseparated N samples: {n_unseparated}",
    f"  Separated samples: {n_separated}",
    f"  E channels: {layer_counts['E']}",
    f"  L channels: {layer_counts['L']}",
    f"  N channels: {layer_counts['N']}",
    "",
    "Special sample pairing:",
    "  N127693_E + N127693_L1 -> N127693_1",
    "  I127693_E + I127693_L2 -> I127693_2",
    "",
    "Previously validated exceptions:",
    "  N104689_N2 -> N104689_N channel alias",
    "  TI_IMM 376 HBB/HBA cells -> Portal V2 Monocytes HBB",
    "",
    "FINAL STATUS: PASS",
]

(
    OUT
    / "PHASE0_FINAL_SUMMARY.txt"
).write_text(
    "\n".join(report) + "\n",
    encoding="utf-8",
)


print()
print("=" * 78)
print("SCP1884 PHASE 0 FINAL RESULTS")
print("=" * 78)

print(
    f"Cells               : {len(cells):,}"
)
print(
    f"Channels            : {len(channel)}"
)
print(
    f"Biological samples  : {len(sample)}"
)
print(
    f"Donors              : {len(donor)}"
)
print(
    f"CD / Healthy donors : {n_cd} / {n_healthy}"
)

print()
print(
    "Corrected pairings:"
)
print(
    "  N127693_E + N127693_L1 -> N127693_1"
)
print(
    "  I127693_E + I127693_L2 -> I127693_2"
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
print(
    "FINAL STATUS: PASS"
)
print(
    f"Outputs: {OUT}"
)
