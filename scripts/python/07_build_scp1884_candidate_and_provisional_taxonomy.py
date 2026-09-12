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

ENV_ROOT = Path(
    "/home/mazekai/miniconda3/envs/ibd_refaudit"
)

SOURCE = (
    PROJECT
    / "03_reference/SCP1884/audit_v1/output/"
    "01_barcode_sample_structure_FINAL/"
    "scp1884_cells_with_source_structure_FINAL.parquet"
)

OUT = (
    PROJECT
    / "03_reference/SCP1884/audit_v1/output/"
    "02_candidate_taxonomy"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

EXPECTED_CELLS = 720_633
EXPECTED_SAMPLES = 136
EXPECTED_DONORS = 71
EXPECTED_ANNO2 = 66


# ============================================================
# Helpers
# ============================================================

def pkg(name):

    try:
        return md.version(name)

    except md.PackageNotFoundError:
        return "NOT INSTALLED"


def fail(msg):

    raise RuntimeError(msg)


def collapse_values(x):

    return " | ".join(
        sorted(
            set(
                x.dropna().astype(str)
            )
        )
    )


# ============================================================
# Provisional taxonomy
#
# IMPORTANT:
# This is dataset-level provisional harmonization.
# It is NOT the final cross-dataset taxonomy.
# ============================================================

def map_annotation(label):

    x = str(label)

    # --------------------------------------------------------
    # B / plasma
    # --------------------------------------------------------

    if x == "Plasma cells":
        return (
            "Plasma",
            "KEEP",
            "",
            "Plasma-cell identity."
        )

    if x.startswith("B cells"):
        return (
            "B",
            "KEEP",
            "",
            "B-cell identity."
        )

    # --------------------------------------------------------
    # T / NK / ILC
    # --------------------------------------------------------

    if (
        x.startswith("T cells CD4")
        or x == "Tregs"
        or x == "T cells Naive CD4+"
    ):
        return (
            "CD4_T",
            "KEEP",
            "",
            "CD4 T-cell lineage."
        )

    if x.startswith("T cells CD8"):
        return (
            "CD8_T",
            "KEEP",
            "",
            "CD8 T-cell lineage."
        )

    if x == "T cells OGT+":
        return (
            "T_ambiguous",
            "HOLD_CROSS_DATASET",
            "",
            "T-cell state without sufficiently specific "
            "CD4/CD8 identity from label alone."
        )

    if x == "IELs ID3+ ENTPD1+":
        return (
            "Unconventional_T",
            "HOLD_CROSS_DATASET",
            "",
            "IEL / unconventional T-cell population."
        )

    if x == "NK-like cells ID3+ ENTPD1+":
        return (
            "Unconventional_T",
            "HOLD_CROSS_DATASET",
            "",
            "NK-like / unconventional lymphocyte population."
        )

    if x == "NK cells KLRF1+ CD3G-":
        return (
            "NK",
            "KEEP",
            "",
            "Canonical NK population."
        )

    if x == "ILCs":
        return (
            "ILC",
            "KEEP",
            "",
            "Innate lymphoid cells."
        )

    # --------------------------------------------------------
    # Myeloid
    # --------------------------------------------------------

    if (
        x.startswith("Macrophages")
        or x.startswith("Monocytes")
    ):
        return (
            "Monocyte_Macrophage",
            "KEEP",
            "",
            "Monocyte/macrophage lineage."
        )

    if (
        x.startswith("DC1")
        or x.startswith("DC2")
        or x == "Mature DCs"
    ):
        return (
            "DC",
            "KEEP",
            "",
            "Dendritic-cell lineage."
        )

    if x == "Mast cells":
        return (
            "Mast",
            "KEEP",
            "",
            "Mast-cell identity."
        )

    if x.startswith("Neutrophils"):
        return (
            "Neutrophil",
            "EXCLUDE_PRIMARY",
            "",
            "Biologically valid neutrophils, but excluded "
            "from current primary-reference taxonomy."
        )

    # --------------------------------------------------------
    # Stromal
    # --------------------------------------------------------

    if (
        x.startswith("Fibroblasts")
        or x.startswith("Myofibroblasts")
        or x.startswith("Activated fibroblasts")
        or x.startswith("Inflammatory fibroblasts")
    ):
        return (
            "Fibroblast",
            "KEEP",
            "",
            "Fibroblast/myofibroblast lineage."
        )

    if x.startswith("Pericytes"):
        return (
            "Pericyte",
            "KEEP",
            "",
            "Pericyte identity."
        )

    if (
        x.startswith("Endothelial cells")
        or x == "Lymphatics"
    ):
        return (
            "Endothelial",
            "KEEP",
            "lymphatic" if x == "Lymphatics" else "",
            "Vascular/lymphatic endothelial lineage."
        )

    if x == "Glial cells":
        return (
            "Glial",
            "HOLD_CROSS_DATASET",
            "",
            "Enteric glial identity; retained for "
            "cross-dataset decision."
        )

    # --------------------------------------------------------
    # Absorptive epithelial
    #
    # Actual coarse identity is assigned per CELL using Site.
    # --------------------------------------------------------

    if x.startswith("Enterocytes"):
        return (
            "SITE_DEPENDENT_ABSORPTIVE",
            "KEEP",
            "",
            "Absorptive epithelial identity resolved using "
            "anatomical Site."
        )

    # --------------------------------------------------------
    # Secretory / progenitor epithelial
    # --------------------------------------------------------

    if x.startswith("Stem cells"):
        return (
            "Stem_TA_progenitor",
            "KEEP",
            "",
            "Stem/transit-amplifying progenitor lineage."
        )

    if x.startswith("Goblet cells"):
        return (
            "Goblet",
            "KEEP",
            "",
            "Goblet-cell identity."
        )

    if x == "Paneth cells":
        return (
            "Paneth",
            "KEEP",
            "",
            "Paneth identity; retain but review across "
            "datasets before final signature construction."
        )

    if (
        x == "L cells"
        or x == "Enterochromaffin cells"
        or x == "Enteroendocrine cells"
    ):
        return (
            "Enteroendocrine",
            "KEEP",
            "",
            "Enteroendocrine lineage."
        )

    if x == "Tuft cells":
        return (
            "Tuft",
            "HOLD_CROSS_DATASET",
            "",
            "Biologically valid tuft population; "
            "retain pending cross-dataset taxonomy."
        )

    # --------------------------------------------------------
    # Cycling populations
    # --------------------------------------------------------

    if x == "Epithelial Cycling cells":
        return (
            "Epithelial_ambiguous",
            "EXCLUDE_PRIMARY",
            "cycling",
            "Cycling epithelial state is unsuitable for "
            "primary identity signature."
        )

    if x == "Immune Cycling cells":
        return (
            "Immune_ambiguous",
            "EXCLUDE_PRIMARY",
            "cycling",
            "Cycling immune state is unsuitable for "
            "primary identity signature."
        )

    if x == "Stromal Cycling cells":
        return (
            "Stromal_ambiguous",
            "EXCLUDE_PRIMARY",
            "cycling",
            "Cycling stromal state is unsuitable for "
            "primary identity signature."
        )

    # --------------------------------------------------------
    # Ambiguous epithelial populations
    # --------------------------------------------------------

    if x == "Epithelial cells HBB+ HBA+":
        return (
            "HBB_HBA_ambiguous",
            "EXCLUDE_PRIMARY",
            "hemoglobin_program",
            "Known annotation ambiguity; Portal V2 "
            "reclassifies a subset as Monocytes HBB."
        )

    if x == "Epithelial cells METTL12+ MAFB+":
        return (
            "Epithelial_ambiguous",
            "EXCLUDE_PRIMARY",
            "",
            "Insufficiently stable canonical epithelial "
            "identity for primary reference."
        )

    # --------------------------------------------------------
    # Anything unexpected must remain visible
    # --------------------------------------------------------

    return (
        "UNMAPPED",
        "REVIEW",
        "",
        "No provisional rule defined."
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
    "Environment:",
    python_exe
)

print(
    "Python:",
    sys.version.split()[0],
)

print(
    "pandas:",
    pkg("pandas"),
)


# ============================================================
# 1. Load finalized SCP1884 metadata
# ============================================================

if not SOURCE.is_file():
    fail(
        f"Missing finalized Phase-0 file:\n{SOURCE}"
    )

df = pd.read_parquet(
    SOURCE
)

required = [
    "cell_id",
    "PubIDSample",
    "biological_sample_id",
    "PubID",
    "disease",
    "inflammation",
    "Site",
    "site_canonical",
    "author_analysis_region",
    "anno_overall",
    "anno2",
    "v2_celltype",
    "expression_bundle",
    "Layer",
    "Chem",
]

missing = [
    x
    for x in required
    if x not in df.columns
]

if missing:
    fail(
        f"Missing required columns: {missing}"
    )

if len(df) != EXPECTED_CELLS:
    fail(
        f"Cells={len(df):,}; expected={EXPECTED_CELLS:,}"
    )

if not df["cell_id"].is_unique:
    fail(
        "cell_id is not unique."
    )

if df["biological_sample_id"].nunique() != EXPECTED_SAMPLES:
    fail(
        "Expected 136 biological samples."
    )

if df["PubID"].nunique() != EXPECTED_DONORS:
    fail(
        "Expected 71 donors."
    )

print()
print("[PASS] Finalized Phase-0 structure loaded.")


# ============================================================
# 2. Phase 1 candidate eligibility
#
# SCP1884 contains only:
# Healthy
# CD Inflamed
# CD Non_Inflamed
#
# There is no treatment-phase filter comparable to GSE282122.
# Therefore all author-QC cells are eligible at the
# biological-sample selection stage.
# ============================================================

valid = (
    (
        df["disease"].eq("Healthy")
        &
        df["inflammation"].eq("Healthy")
    )
    |
    (
        df["disease"].eq("CD")
        &
        df["inflammation"].isin(
            [
                "Inflamed",
                "Non_Inflamed",
            ]
        )
    )
)

if not valid.all():

    bad = df.loc[
        ~valid
    ].copy()

    bad.to_csv(
        OUT
        / "ERROR_ineligible_metadata_states.tsv",
        sep="\t",
        index=False,
    )

    fail(
        f"{len(bad):,} cells have unexpected "
        "disease/inflammation states."
    )

candidate = df.copy()

candidate["dataset"] = "SCP1884"

candidate["candidate_status"] = (
    "ELIGIBLE"
)

candidate["candidate_reason"] = (
    candidate["disease"]
    + "_"
    + candidate["inflammation"]
)

if len(candidate) != EXPECTED_CELLS:
    fail(
        "Candidate pool does not contain all cells."
    )

print(
    "[PASS] Phase 1 candidate pool:",
    f"{len(candidate):,} cells"
)

print(
    "[PASS] Biological samples:",
    candidate["biological_sample_id"].nunique()
)

print(
    "[PASS] Donors:",
    candidate["PubID"].nunique()
)


# ============================================================
# 3. Candidate sampling summary
# ============================================================

candidate_summary = (
    candidate.groupby(
        [
            "disease",
            "inflammation",
            "Site",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),
        n_samples=(
            "biological_sample_id",
            "nunique",
        ),
        n_donors=(
            "PubID",
            "nunique",
        ),
    )
    .reset_index()
)

candidate_summary.to_csv(
    OUT
    / "candidate_sampling_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 4. Phase 2A annotation audit
# ============================================================

if candidate["anno2"].nunique() != EXPECTED_ANNO2:
    fail(
        "Expected 66 unique author anno2 labels; "
        f"found {candidate['anno2'].nunique()}."
    )

anno2_summary = (
    candidate.groupby(
        [
            "anno2",
            "anno_overall",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),
        n_samples=(
            "biological_sample_id",
            "nunique",
        ),
        n_donors=(
            "PubID",
            "nunique",
        ),
        sites=(
            "Site",
            collapse_values,
        ),
        diseases=(
            "disease",
            collapse_values,
        ),
        inflammation_states=(
            "inflammation",
            collapse_values,
        ),
        v2_celltypes=(
            "v2_celltype",
            collapse_values,
        ),
    )
    .reset_index()
)

anno2_summary[
    "fraction_of_candidate"
] = (
    anno2_summary["n_cells"]
    / len(candidate)
)

anno2_summary = (
    anno2_summary
    .sort_values(
        "n_cells",
        ascending=False,
    )
    .reset_index(
        drop=True
    )
)


# Site-specific coverage

site_counts = (
    pd.crosstab(
        candidate["anno2"],
        candidate["Site"],
    )
    .reset_index()
)

site_counts = site_counts.rename(
    columns={
        "CO": "n_CO",
        "TI": "n_TI",
        "SB": "n_SB",
    }
)

for col in [
    "n_CO",
    "n_TI",
    "n_SB",
]:

    if col not in site_counts.columns:
        site_counts[col] = 0


anno2_summary = anno2_summary.merge(
    site_counts,
    on="anno2",
    how="left",
    validate="one_to_one",
)


anno2_summary.to_csv(
    OUT
    / "anno2_taxonomy_audit.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 5. anno2 x Portal V2 audit
# ============================================================

anno2_v2 = (
    candidate.groupby(
        [
            "anno2",
            "v2_celltype",
        ],
        observed=True,
    )
    .size()
    .rename(
        "n_cells"
    )
    .reset_index()
    .sort_values(
        [
            "anno2",
            "n_cells",
        ],
        ascending=[
            True,
            False,
        ],
    )
)

anno2_v2.to_csv(
    OUT
    / "anno2_vs_portal_v2.tsv",
    sep="\t",
    index=False,
)

v2_splits = (
    anno2_v2.groupby(
        "anno2",
        observed=True,
    )[
        "v2_celltype"
    ]
    .nunique()
)

v2_split_labels = (
    v2_splits[
        v2_splits > 1
    ]
    .index
    .tolist()
)


# ============================================================
# 6. Provisional mapping table
# ============================================================

mapping = (
    anno2_summary[
        [
            "anno2",
            "anno_overall",
            "n_cells",
            "n_samples",
            "n_donors",
            "n_CO",
            "n_TI",
            "n_SB",
            "sites",
            "v2_celltypes",
        ]
    ]
    .copy()
)

mapped_values = (
    mapping["anno2"]
    .map(
        map_annotation
    )
)

mapping[
    "provisional_coarse_identity"
] = [
    x[0]
    for x in mapped_values
]

mapping[
    "primary_reference_action"
] = [
    x[1]
    for x in mapped_values
]

mapping[
    "state_flag"
] = [
    x[2]
    for x in mapped_values
]

mapping[
    "mapping_rationale"
] = [
    x[3]
    for x in mapped_values
]

mapping[
    "review_notes"
] = ""


if (
    mapping[
        "provisional_coarse_identity"
    ]
    .eq(
        "UNMAPPED"
    )
    .any()
):

    unmapped = mapping.loc[
        mapping[
            "provisional_coarse_identity"
        ]
        .eq(
            "UNMAPPED"
        )
    ]

    unmapped.to_csv(
        OUT
        / "UNMAPPED_annotations.tsv",
        sep="\t",
        index=False,
    )

    print()
    print(
        "[WARNING] Some anno2 labels are UNMAPPED:"
    )

    print(
        unmapped[
            [
                "anno2",
                "n_cells",
            ]
        ].to_string(
            index=False
        )
    )


mapping.to_csv(
    OUT
    / "taxonomy_mapping_provisional_v1.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 7. Apply provisional mapping to cells
# ============================================================

map_lookup = (
    mapping[
        [
            "anno2",
            "provisional_coarse_identity",
            "primary_reference_action",
            "state_flag",
        ]
    ]
    .set_index(
        "anno2"
    )
)

candidate = candidate.join(
    map_lookup,
    on="anno2",
)


# ------------------------------------------------------------
# Resolve absorptive identity using anatomical site.
# CO -> Colonic_absorptive
# TI/SB -> Ileal_absorptive
# ------------------------------------------------------------

absorptive = (
    candidate[
        "provisional_coarse_identity"
    ]
    .eq(
        "SITE_DEPENDENT_ABSORPTIVE"
    )
)

candidate.loc[
    absorptive
    &
    candidate[
        "Site"
    ].eq(
        "CO"
    ),
    "provisional_coarse_identity",
] = "Colonic_absorptive"

candidate.loc[
    absorptive
    &
    candidate[
        "Site"
    ].isin(
        [
            "TI",
            "SB",
        ]
    ),
    "provisional_coarse_identity",
] = "Ileal_absorptive"


# ------------------------------------------------------------
# Explicitly retain provenance of the 376-cell HBB exception.
#
# Biological interpretation:
# Portal V2 -> Monocytes HBB.
#
# Primary-reference action remains EXCLUDE_PRIMARY because the
# HBB/HBA program is undesirable for stable identity signatures.
# ------------------------------------------------------------

hbb_monocyte = (
    candidate["anno2"]
    .eq(
        "Epithelial cells HBB+ HBA+"
    )
    &
    candidate["v2_celltype"]
    .eq(
        "Monocytes HBB"
    )
)

if int(
    hbb_monocyte.sum()
) != 376:

    fail(
        "Expected 376 HBB/HBA -> Monocytes HBB cells; "
        f"found {int(hbb_monocyte.sum())}."
    )

candidate.loc[
    hbb_monocyte,
    "provisional_coarse_identity",
] = "Monocyte_Macrophage"

candidate.loc[
    hbb_monocyte,
    "state_flag",
] = "HBB_HBA_annotation_exception"


# ============================================================
# 8. Mapping summary
# ============================================================

mapping_counts = (
    candidate.groupby(
        [
            "provisional_coarse_identity",
            "primary_reference_action",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "cell_id",
            "size",
        ),
        n_samples=(
            "biological_sample_id",
            "nunique",
        ),
        n_donors=(
            "PubID",
            "nunique",
        ),
    )
    .reset_index()
    .sort_values(
        "n_cells",
        ascending=False,
    )
)

mapping_counts.to_csv(
    OUT
    / "provisional_identity_counts.tsv",
    sep="\t",
    index=False,
)


candidate.to_parquet(
    OUT
    / "scp1884_candidate_cells_provisional_taxonomy.parquet",
    index=False,
    engine="pyarrow",
)


# ============================================================
# 9. Console report
# ============================================================

print()
print("=" * 78)
print("SCP1884 PHASE 1 + PHASE 2A/2B PROVISIONAL RESULTS")
print("=" * 78)

print()
print(
    f"Candidate cells      : {len(candidate):,}"
)
print(
    f"Biological samples  : "
    f"{candidate['biological_sample_id'].nunique()}"
)
print(
    f"Donors              : "
    f"{candidate['PubID'].nunique()}"
)
print(
    f"Author anno2 labels : "
    f"{candidate['anno2'].nunique()}"
)

print()
print("Candidate sampling:")
print(
    candidate_summary.to_string(
        index=False
    )
)

print()
print("Portal V2 split labels:")
if v2_split_labels:
    for x in v2_split_labels:
        print(
            f"  - {x}"
        )
else:
    print(
        "  None"
    )

print()
print("Provisional mapping:")
print(
    mapping_counts.to_string(
        index=False
    )
)

print()
print(
    "Mapping actions:"
)

print(
    candidate[
        "primary_reference_action"
    ]
    .value_counts()
    .to_string()
)

print()
print(
    "Output directory:"
)
print(
    OUT
)
