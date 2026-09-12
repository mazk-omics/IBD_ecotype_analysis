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

INPUT_TEMPLATE = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "02_taxonomy_mapping"
    / "taxonomy_mapping_template.tsv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "03_reference"
    / "GSE282122"
    / "audit_v2"
    / "output"
    / "03_taxonomy_mapping_v1"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "taxonomy_mapping_v1_summary.txt"
)


# ============================================================
# Taxonomy definitions
# ============================================================

CORE_IDENTITIES = {
    "B",
    "Plasma",
    "CD4_T",
    "CD8_T",
    "NK",
    "Monocyte_Macrophage",
    "DC",
    "Mast",
    "Fibroblast",
    "Pericyte",
    "Endothelial",
    "Ileal_absorptive",
    "Colonic_absorptive",
    "Goblet",
    "Stem_TA_progenitor",
    "Paneth",
    "Enteroendocrine",
    "ILC",
}

PROVISIONAL_IDENTITIES = {
    "Unconventional_T",
    "Tuft",
    "M_like",
    "Glial",
    "Stromal_ambiguous",
}

ALLOWED_ACTIONS = {
    "KEEP",
    "EXCLUDE_PRIMARY",
    "HOLD_CROSS_DATASET",
}

ALLOWED_IDENTITY_STATUS = {
    "CONFIDENT",
    "PROVISIONAL",
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
    summary.append(f"[FAIL] {message}")
    write_summary(summary)
    raise RuntimeError(message)


def path_ids(start, end):
    """
    Inclusive path range:
    path_ids(13, 20)
    -> GSE282122_PATH_013 ... PATH_020
    """
    return [
        f"GSE282122_PATH_{i:03d}"
        for i in range(start, end + 1)
    ]


# ============================================================
# 0. Environment preflight
# ============================================================

summary = []

summary.append(
    "GSE282122 Reference Audit v2"
)
summary.append(
    "Phase 2B: Provisional taxonomy mapping v1"
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
summary.append("")


# ============================================================
# 1. Load Phase-2A template
# ============================================================

summary.append("Input")
summary.append("-" * 76)

summary.append(
    f"Mapping template: {INPUT_TEMPLATE}"
)

if not INPUT_TEMPLATE.is_file():
    fail(
        f"Input template not found: {INPUT_TEMPLATE}",
        summary,
    )

df = pd.read_csv(
    INPUT_TEMPLATE,
    sep="\t",
)

summary.append(
    f"Annotation paths: {len(df)}"
)
summary.append("")

if len(df) != 110:
    fail(
        f"Expected 110 annotation paths; found {len(df)}.",
        summary,
    )

if df["path_id"].duplicated().any():
    fail(
        "Duplicate path_id values detected.",
        summary,
    )


# ============================================================
# 2. New mapping fields
# ============================================================

mapping = {}

def assign(
    ids,
    coarse_identity,
    action,
    rationale,
    identity_status="CONFIDENT",
    state_flag="none",
    review_notes="",
):
    """
    Add one or more path IDs to manual mapping.
    Duplicate assignments are prohibited.
    """

    for pid in ids:

        if pid in mapping:
            raise RuntimeError(
                f"Duplicate mapping assignment: {pid}"
            )

        mapping[pid] = {
            "coarse_identity": coarse_identity,
            "identity_status": identity_status,
            "core_taxonomy_member":
                coarse_identity in CORE_IDENTITIES,
            "primary_reference_action": action,
            "state_flag": state_flag,
            "mapping_rationale": rationale,
            "review_notes": review_notes,
        }


# ============================================================
# 3. B lineage
# ============================================================

assign(
    path_ids(1, 5)
    + path_ids(7, 12),
    "B",
    "KEEP",
    (
        "Author hierarchy consistently identifies "
        "B-cell lineage; retain for coarse B identity."
    ),
)

assign(
    ["GSE282122_PATH_006"],
    "B",
    "EXCLUDE_PRIMARY",
    (
        "Author hierarchy identifies B-cell lineage, "
        "but final state is Cycling B."
    ),
    state_flag="cycling",
    review_notes=(
        "Retain B identity metadata; exclude from "
        "primary signature learning to limit "
        "cell-cycle-driven features."
    ),
)


# ============================================================
# 4. Conventional CD4 T cells
# ============================================================

assign(
    path_ids(13, 32),
    "CD4_T",
    "KEEP",
    (
        "Author hierarchy identifies conventional "
        "CD4 T-cell lineage, including helper, memory, "
        "naive and Treg states."
    ),
)


# ============================================================
# 5. Conventional CD8 T cells
# ============================================================

assign(
    path_ids(33, 40)
    + ["GSE282122_PATH_042"],
    "CD8_T",
    "KEEP",
    (
        "Author hierarchy identifies conventional "
        "CD8 T-cell lineage."
    ),
)

assign(
    ["GSE282122_PATH_041"],
    "CD8_T",
    "EXCLUDE_PRIMARY",
    (
        "Author hierarchy identifies CD8 T-cell lineage, "
        "but final state is cycling."
    ),
    state_flag="cycling",
    review_notes=(
        "Retain CD8 identity metadata; exclude from "
        "primary signature learning."
    ),
)


# ============================================================
# 6. Unconventional T cells
# ============================================================

assign(
    [
        "GSE282122_PATH_043",
        "GSE282122_PATH_045",
        "GSE282122_PATH_046",
    ],
    "Unconventional_T",
    "HOLD_CROSS_DATASET",
    (
        "MAIT and gamma-delta T cells are distinct "
        "unconventional T-cell populations and should "
        "not be forced into conventional CD4/CD8 classes."
    ),
    review_notes=(
        "Final reference inclusion depends on "
        "cross-dataset representation."
    ),
)


# ============================================================
# 7. NK
# ============================================================

assign(
    ["GSE282122_PATH_044"],
    "NK",
    "KEEP",
    (
        "Author hierarchy explicitly identifies "
        "natural killer cells."
    ),
)


# ============================================================
# 8. Ileal epithelium
# ============================================================

assign(
    [
        "GSE282122_PATH_047",
        "GSE282122_PATH_050",
        "GSE282122_PATH_052",
        "GSE282122_PATH_053",
        "GSE282122_PATH_054",
        "GSE282122_PATH_055",
    ],
    "Ileal_absorptive",
    "KEEP",
    (
        "Ileal enterocyte/BEST4-related populations "
        "mapped to ileal absorptive epithelial identity."
    ),
)

assign(
    [
        "GSE282122_PATH_048",
        "GSE282122_PATH_049",
    ],
    "Enteroendocrine",
    "KEEP",
    (
        "Author annotation identifies ileal "
        "enteroendocrine cells."
    ),
)

assign(
    [
        "GSE282122_PATH_051",
        "GSE282122_PATH_058",
    ],
    "Stem_TA_progenitor",
    "KEEP",
    (
        "Ileal TA and LGR5-positive stem populations "
        "mapped to epithelial stem/TA progenitor identity."
    ),
)

assign(
    ["GSE282122_PATH_056"],
    "Goblet",
    "KEEP",
    "Author annotation identifies ileal goblet cells.",
)

assign(
    ["GSE282122_PATH_057"],
    "Paneth",
    "KEEP",
    "Author annotation identifies ileal Paneth cells.",
)

assign(
    ["GSE282122_PATH_059"],
    "Tuft",
    "HOLD_CROSS_DATASET",
    (
        "Author annotation identifies tuft cells, "
        "which are biologically distinct but outside "
        "the current 18-class core taxonomy."
    ),
    review_notes=(
        "Assess cross-dataset abundance and annotation "
        "consistency before final taxonomy freeze."
    ),
)


# ============================================================
# 9. Non-ileal / colonic epithelium
# ============================================================

assign(
    [
        "GSE282122_PATH_060",
        "GSE282122_PATH_065",
        "GSE282122_PATH_067",
        "GSE282122_PATH_068",
        "GSE282122_PATH_069",
        "GSE282122_PATH_070",
    ],
    "Colonic_absorptive",
    "KEEP",
    (
        "Non-ileal BEST4, colonocyte and enterocyte "
        "populations mapped to colonic/non-ileal "
        "absorptive epithelial identity."
    ),
)

assign(
    [
        "GSE282122_PATH_061",
        "GSE282122_PATH_062",
        "GSE282122_PATH_063",
    ],
    "Enteroendocrine",
    "KEEP",
    (
        "Author annotation identifies non-ileal "
        "enteroendocrine populations."
    ),
)

assign(
    ["GSE282122_PATH_064"],
    "M_like",
    "HOLD_CROSS_DATASET",
    (
        "Author annotation identifies M-like cells, "
        "a distinct epithelial population outside the "
        "current core taxonomy."
    ),
    review_notes=(
        "Final inclusion requires cross-dataset review."
    ),
)

assign(
    [
        "GSE282122_PATH_066",
        "GSE282122_PATH_073",
    ],
    "Stem_TA_progenitor",
    "KEEP",
    (
        "Non-ileal TA and LGR5-positive stem populations "
        "mapped to epithelial stem/TA progenitor identity."
    ),
)

assign(
    ["GSE282122_PATH_071"],
    "Goblet",
    "KEEP",
    "Author annotation identifies non-ileal goblet cells.",
)

assign(
    ["GSE282122_PATH_072"],
    "Paneth",
    "KEEP",
    "Author annotation identifies non-ileal Paneth cells.",
)

assign(
    ["GSE282122_PATH_074"],
    "Tuft",
    "HOLD_CROSS_DATASET",
    (
        "Author annotation identifies non-ileal tuft "
        "cells outside the current core taxonomy."
    ),
    review_notes=(
        "Assess cross-dataset representation before "
        "final inclusion."
    ),
)


# ============================================================
# 10. ILC
# ============================================================

assign(
    ["GSE282122_PATH_075"],
    "ILC",
    "KEEP",
    (
        "Author hierarchy explicitly identifies "
        "innate lymphoid cells separate from NK."
    ),
)


# ============================================================
# 11. Mast
# ============================================================

assign(
    ["GSE282122_PATH_076"],
    "Mast",
    "KEEP",
    "Author hierarchy explicitly identifies mast cells.",
)


# ============================================================
# 12. Myeloid
# ============================================================

assign(
    ["GSE282122_PATH_077"],
    "Monocyte_Macrophage",
    "EXCLUDE_PRIMARY",
    (
        "Cycling MNP retains broad mononuclear phagocyte "
        "lineage identity, but cycling state may dominate "
        "signature learning."
    ),
    state_flag="cycling",
)

assign(
    path_ids(78, 82),
    "DC",
    "KEEP",
    (
        "Author hierarchy identifies conventional, "
        "LAMP3, XCR1 and plasmacytoid dendritic-cell "
        "populations."
    ),
)

assign(
    path_ids(83, 86),
    "Monocyte_Macrophage",
    "KEEP",
    (
        "Author hierarchy identifies monocyte/macrophage "
        "lineage."
    ),
)


# ============================================================
# 13. Plasma lineage
# ============================================================

assign(
    path_ids(87, 91),
    "Plasma",
    "KEEP",
    (
        "IgA plasma, IgG plasma and plasmablast "
        "populations mapped to coarse plasma identity."
    ),
)


# ============================================================
# 14. Cycling stroma
# ============================================================

assign(
    ["GSE282122_PATH_092"],
    "Stromal_ambiguous",
    "EXCLUDE_PRIMARY",
    (
        "Cycling stromal cells lack a sufficiently "
        "specific non-cycling stromal lineage assignment."
    ),
    identity_status="PROVISIONAL",
    state_flag="cycling",
    review_notes=(
        "Do not force into fibroblast, pericyte or "
        "other stromal class without additional evidence."
    ),
)


# ============================================================
# 15. Fibroblast lineage
# ============================================================

assign(
    [
        "GSE282122_PATH_093",
        "GSE282122_PATH_094",
        "GSE282122_PATH_095",
        "GSE282122_PATH_096",
        "GSE282122_PATH_097",
        "GSE282122_PATH_098",
        "GSE282122_PATH_105",
        "GSE282122_PATH_107",
    ],
    "Fibroblast",
    "KEEP",
    (
        "Author hierarchy identifies fibroblast or "
        "myofibroblast lineage; myofibroblasts are "
        "collapsed into the coarse fibroblast identity."
    ),
)


# ============================================================
# 16. Pericyte
# ============================================================

assign(
    path_ids(99, 104),
    "Pericyte",
    "KEEP",
    (
        "Author hierarchy identifies arterial/venous "
        "and NOTCH3-positive pericyte populations."
    ),
)


# ============================================================
# 17. Glial
# ============================================================

assign(
    ["GSE282122_PATH_106"],
    "Glial",
    "HOLD_CROSS_DATASET",
    (
        "Author hierarchy explicitly identifies glial "
        "cells with strong donor/sample representation."
    ),
    review_notes=(
        "Biological identity is retained. Final inclusion "
        "as an additional coarse reference identity "
        "requires SCP1884/SCP259 cross-dataset review."
    ),
)


# ============================================================
# 18. Endothelium
# ============================================================

assign(
    path_ids(108, 110),
    "Endothelial",
    "KEEP",
    (
        "Arterial, venous and lymphatic endothelial "
        "populations collapsed into coarse endothelial "
        "identity."
    ),
)


# ============================================================
# 19. Mapping completeness checks
# ============================================================

summary.append("Mapping validation")
summary.append("-" * 76)

input_ids = set(
    df["path_id"].astype(str)
)

mapped_ids = set(
    mapping.keys()
)

missing_ids = sorted(
    input_ids - mapped_ids
)

extra_ids = sorted(
    mapped_ids - input_ids
)

if missing_ids:
    fail(
        "Unmapped path IDs: "
        + ", ".join(missing_ids),
        summary,
    )

if extra_ids:
    fail(
        "Mapping contains unknown path IDs: "
        + ", ".join(extra_ids),
        summary,
    )

if len(mapping) != 110:
    fail(
        f"Expected 110 mappings; found {len(mapping)}.",
        summary,
    )

summary.append(
    "[PASS] 110/110 annotation paths mapped exactly once."
)
summary.append("")


# ============================================================
# 20. Apply path-level mapping to table
#
# Phase-2A template already contains empty placeholder
# mapping columns. Remove those placeholders before merging
# the finalized Phase-2B mapping, otherwise pandas will create
# _x / _y suffixed duplicate columns.
# ============================================================

mapping_df = (
    pd.DataFrame
    .from_dict(
        mapping,
        orient="index",
    )
    .reset_index()
    .rename(
        columns={
            "index": "path_id"
        }
    )
)

PHASE2A_PLACEHOLDER_COLUMNS = [
    "target_identity",
    "mapping_status",
    "mapping_rationale",
    "exclude_reason",
    "review_notes",
]

df_for_merge = df.drop(
    columns=PHASE2A_PLACEHOLDER_COLUMNS,
    errors="ignore",
)

out = df_for_merge.merge(
    mapping_df,
    on="path_id",
    how="left",
    validate="one_to_one",
)


# ============================================================
# 21. Vocabulary validation
# ============================================================

allowed_identities = (
    CORE_IDENTITIES
    | PROVISIONAL_IDENTITIES
)

observed_identities = set(
    out["coarse_identity"]
)

unexpected_identities = (
    observed_identities
    - allowed_identities
)

if unexpected_identities:
    fail(
        "Unexpected coarse identities: "
        + ", ".join(
            sorted(unexpected_identities)
        ),
        summary,
    )

unexpected_actions = (
    set(out["primary_reference_action"])
    - ALLOWED_ACTIONS
)

if unexpected_actions:
    fail(
        "Unexpected primary actions: "
        + ", ".join(
            sorted(unexpected_actions)
        ),
        summary,
    )

unexpected_status = (
    set(out["identity_status"])
    - ALLOWED_IDENTITY_STATUS
)

if unexpected_status:
    fail(
        "Unexpected identity statuses: "
        + ", ".join(
            sorted(unexpected_status)
        ),
        summary,
    )

summary.append("Vocabulary validation")
summary.append("-" * 76)

summary.append(
    "[PASS] Identity vocabulary validated."
)
summary.append(
    "[PASS] Reference-action vocabulary validated."
)
summary.append(
    "[PASS] Identity-status vocabulary validated."
)
summary.append("")


# ============================================================
# 22. Output mapping
# ============================================================

keep_front = [
    "path_id",
    "bucket",
    "sub_bucket",
    "major",
    "minor",
    "final_analysis",
    "n_cells",
    "n_samples",
    "n_donors",
    "max_donor_fraction",
    "donor_HHI",
    "coarse_identity",
    "identity_status",
    "core_taxonomy_member",
    "primary_reference_action",
    "state_flag",
    "mapping_rationale",
    "review_notes",
]

remaining = [
    c
    for c in out.columns
    if c not in keep_front
]

out = out[
    keep_front + remaining
]

out.to_csv(
    OUTPUT_DIR
    / "taxonomy_mapping_v1.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 23. Review subsets
# ============================================================

out.loc[
    out["primary_reference_action"]
    == "KEEP"
].to_csv(
    OUTPUT_DIR
    / "taxonomy_mapping_v1__KEEP.tsv",
    sep="\t",
    index=False,
)

out.loc[
    out["primary_reference_action"]
    == "EXCLUDE_PRIMARY"
].to_csv(
    OUTPUT_DIR
    / "taxonomy_mapping_v1__EXCLUDE_PRIMARY.tsv",
    sep="\t",
    index=False,
)

out.loc[
    out["primary_reference_action"]
    == "HOLD_CROSS_DATASET"
].to_csv(
    OUTPUT_DIR
    / "taxonomy_mapping_v1__HOLD_CROSS_DATASET.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 24. Summary by identity/action
# ============================================================

identity_summary = (
    out
    .groupby(
        [
            "coarse_identity",
            "primary_reference_action",
        ],
        observed=True,
    )
    .agg(
        n_paths=("path_id", "size"),
        n_cells=("n_cells", "sum"),
    )
    .reset_index()
    .sort_values(
        [
            "primary_reference_action",
            "coarse_identity",
        ]
    )
)

identity_summary.to_csv(
    OUTPUT_DIR
    / "taxonomy_mapping_v1_summary_by_identity.tsv",
    sep="\t",
    index=False,
)

action_summary = (
    out
    .groupby(
        "primary_reference_action",
        observed=True,
    )
    .agg(
        n_paths=("path_id", "size"),
        n_cells=("n_cells", "sum"),
    )
    .reset_index()
)

action_summary.to_csv(
    OUTPUT_DIR
    / "taxonomy_mapping_v1_summary_by_action.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 25. Final summary
# ============================================================

summary.append("Mapping overview")
summary.append("-" * 76)

for row in action_summary.itertuples(
    index=False
):
    summary.append(
        f"{row.primary_reference_action}: "
        f"{row.n_paths} paths, "
        f"{row.n_cells:,} cells"
    )

summary.append("")

summary.append(
    f"Core identities represented: "
    f"{out.loc[out['core_taxonomy_member'], 'coarse_identity'].nunique()}"
)

summary.append(
    f"Provisional identities represented: "
    f"{out.loc[~out['core_taxonomy_member'], 'coarse_identity'].nunique()}"
)

summary.append("")

summary.append("Final Phase-2B status")
summary.append("-" * 76)

summary.append(
    "[PASS] 110/110 paths manually mapped."
)
summary.append(
    "[PASS] Biological identity separated from "
    "primary-reference inclusion decision."
)
summary.append(
    "[PASS] Cycling populations retain lineage "
    "identity but are excluded from primary "
    "signature learning."
)
summary.append(
    "[PASS] Glial, unconventional T, tuft and "
    "M-like populations retained for "
    "cross-dataset review."
)
summary.append(
    "[PASS] No candidate cell metadata was modified."
)
summary.append(
    "[PASS] Mapping remains provisional until "
    "cross-dataset harmonization."
)

write_summary(summary)


# ============================================================
# Console output
# ============================================================

print()
print("=" * 76)
print(
    "GSE282122 Phase 2B taxonomy mapping v1 completed"
)
print("=" * 76)

print()
print(
    action_summary.to_string(index=False)
)

print()
print(
    identity_summary.to_string(index=False)
)

print()
print(
    f"Mapping:"
)
print(
    f"  {OUTPUT_DIR / 'taxonomy_mapping_v1.tsv'}"
)

print()
print(
    f"Summary:"
)
print(
    f"  {SUMMARY_FILE}"
)