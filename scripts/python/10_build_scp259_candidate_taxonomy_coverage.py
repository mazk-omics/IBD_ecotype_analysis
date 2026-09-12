#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


PROJECT = Path("/home/mazekai/IBD_EcoTyper")
ENV_ROOT = Path("/home/mazekai/miniconda3/envs/ibd_refaudit")

SOURCE = (
    PROJECT
    / "03_reference/SCP259/audit_v1/output/"
    "00_source_schema/"
    "scp259_metadata_with_expression_bundle.parquet"
)

OUT = (
    PROJECT
    / "03_reference/SCP259/audit_v1/output/"
    "01_candidate_taxonomy_coverage"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

EXPECTED_CELLS = 365_492
EXPECTED_DONORS = 30
EXPECTED_SAMPLES = 133
EXPECTED_CLUSTERS = 51

EXPECTED_PRIMARY = {
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
    "Colonic_absorptive",
    "Goblet",
    "Stem_TA_progenitor",
    "Enteroendocrine",
    "ILC",
}


# ============================================================
# Provisional cluster mapping
#
# coarse_identity describes biology.
# primary_reference_action describes suitability for the
# primary CIBERSORTx identity-reference matrix.
# ============================================================

MAPPING = {

    # B / plasma
    "Plasma":
        ("Plasma", "KEEP", "", "Plasma-cell identity."),

    "Follicular":
        ("B", "KEEP", "", "Follicular B-cell identity."),

    "GC":
        ("B", "KEEP", "germinal_center", "Germinal-center B cells."),

    "Cycling B":
        ("B", "EXCLUDE_PRIMARY", "cycling",
         "B lineage retained biologically; cycling state excluded from primary signature."),

    # CD4 T
    "CD4+ Memory":
        ("CD4_T", "KEEP", "memory", "CD4 T-cell lineage."),

    "CD4+ Activated Fos-hi":
        ("CD4_T", "KEEP", "activated_Fos_hi", "Activated CD4 T-cell state."),

    "CD4+ Activated Fos-lo":
        ("CD4_T", "KEEP", "activated_Fos_lo", "Activated CD4 T-cell state."),

    "Tregs":
        ("CD4_T", "KEEP", "Treg", "Regulatory CD4 T cells."),

    "CD4+ PD1+":
        ("CD4_T", "KEEP", "PD1_positive", "PD1-positive CD4 T-cell state."),

    # CD8 T
    "CD8+ LP":
        ("CD8_T", "KEEP", "LP", "CD8 T-cell lineage."),

    "CD8+ IELs":
        ("CD8_T", "KEEP", "IEL", "Author explicitly annotates CD8-positive IELs."),

    "CD8+ IL17+":
        ("CD8_T", "KEEP", "IL17_positive", "CD8-positive IL17 state."),

    "Cycling T":
        ("T_ambiguous", "EXCLUDE_PRIMARY", "cycling",
         "Cycling T-cell state lacks stable CD4/CD8 identity for primary signature."),

    # NK / ILC
    "NKs":
        ("NK", "KEEP", "", "NK identity."),

    "ILCs":
        ("ILC", "KEEP", "", "Innate lymphoid-cell identity."),

    # Myeloid
    "Macrophages":
        ("Monocyte_Macrophage", "KEEP", "macrophage", "Macrophage lineage."),

    "Inflammatory Monocytes":
        ("Monocyte_Macrophage", "KEEP", "inflammatory_monocyte",
         "Inflammatory monocyte state retained in parent lineage."),

    "Cycling Monocytes":
        ("Monocyte_Macrophage", "EXCLUDE_PRIMARY", "cycling",
         "Monocyte lineage retained biologically; cycling state excluded."),

    "DC1":
        ("DC", "KEEP", "DC1", "Dendritic-cell lineage."),

    "DC2":
        ("DC", "KEEP", "DC2", "Dendritic-cell lineage."),

    "MT-hi":
        ("Immune_ambiguous", "EXCLUDE_PRIMARY", "metallothionein_high",
         "State-dominated immune population unsuitable for identity signature."),

    # Mast
    "CD69+ Mast":
        ("Mast", "KEEP", "CD69_positive", "Mast-cell identity."),

    "CD69- Mast":
        ("Mast", "KEEP", "CD69_negative", "Mast-cell identity."),

    # Fibroblast / stromal
    "WNT2B+ Fos-lo 1":
        ("Fibroblast", "KEEP", "WNT2B_Fos_lo_1", "Fibroblast state."),

    "WNT2B+ Fos-lo 2":
        ("Fibroblast", "KEEP", "WNT2B_Fos_lo_2", "Fibroblast state."),

    "WNT2B+ Fos-hi":
        ("Fibroblast", "KEEP", "WNT2B_Fos_hi", "Fibroblast state."),

    "WNT5B+ 1":
        ("Fibroblast", "KEEP", "WNT5B_1", "Fibroblast state."),

    "WNT5B+ 2":
        ("Fibroblast", "KEEP", "WNT5B_2", "Fibroblast state."),

    "Inflammatory Fibroblasts":
        ("Fibroblast", "KEEP", "inflammatory",
         "Inflammatory state retained in fibroblast parent lineage."),

    "Myofibroblasts":
        ("Fibroblast", "KEEP", "myofibroblast",
         "Myofibroblast mapped to fibroblast coarse identity."),

    "RSPO3+":
        ("Fibroblast", "KEEP", "RSPO3_positive", "Stromal fibroblast state."),

    # Vascular
    "Endothelial":
        ("Endothelial", "KEEP", "endothelial", "Endothelial identity."),

    "Post-capillary Venules":
        ("Endothelial", "KEEP", "post_capillary_venule",
         "Vascular endothelial subtype."),

    "Microvascular":
        ("Endothelial", "KEEP", "microvascular",
         "Microvascular endothelial subtype."),

    "Pericytes":
        ("Pericyte", "KEEP", "", "Pericyte identity."),

    # Glia
    "Glia":
        ("Glial", "HOLD_CROSS_DATASET", "",
         "Enteric glial identity; retain for final cross-dataset decision."),

    # Stem / TA / progenitor
    "Stem":
        ("Stem_TA_progenitor", "KEEP", "stem", "Stem/progenitor lineage."),

    "TA 1":
        ("Stem_TA_progenitor", "KEEP", "TA1", "Transit-amplifying progenitors."),

    "TA 2":
        ("Stem_TA_progenitor", "KEEP", "TA2", "Transit-amplifying progenitors."),

    "Enterocyte Progenitors":
        ("Stem_TA_progenitor", "KEEP", "enterocyte_progenitor",
         "Enterocyte progenitor state."),

    "Secretory TA":
        ("Stem_TA_progenitor", "KEEP", "secretory_TA",
         "Secretory transit-amplifying progenitors."),

    "Cycling TA":
        ("Stem_TA_progenitor", "EXCLUDE_PRIMARY", "cycling",
         "Biological progenitor identity retained; cycling state excluded."),

    # Absorptive epithelium
    # SCP259 is a colon-mucosa dataset.
    "Immature Enterocytes 1":
        ("Colonic_absorptive", "KEEP", "immature_1",
         "Colonic absorptive epithelial lineage."),

    "Immature Enterocytes 2":
        ("Colonic_absorptive", "KEEP", "immature_2",
         "Colonic absorptive epithelial lineage."),

    "Enterocytes":
        ("Colonic_absorptive", "KEEP", "mature",
         "Colonic absorptive epithelial lineage."),

    "Best4+ Enterocytes":
        ("Colonic_absorptive", "KEEP", "BEST4_positive",
         "BEST4-positive colonic absorptive epithelial state."),

    # Goblet
    "Immature Goblet":
        ("Goblet", "KEEP", "immature", "Goblet-cell lineage."),

    "Goblet":
        ("Goblet", "KEEP", "mature", "Goblet-cell identity."),

    # Other epithelial
    "Enteroendocrine":
        ("Enteroendocrine", "KEEP", "", "Enteroendocrine identity."),

    "Tuft":
        ("Tuft", "HOLD_CROSS_DATASET", "",
         "Canonical tuft identity; retain for final cross-dataset decision."),

    "M cells":
        ("M_like", "HOLD_CROSS_DATASET", "",
         "Microfold/M-like epithelial population; retain pending harmonization."),
}


def fail(msg):
    raise RuntimeError(msg)


# ============================================================
# 0. Environment / input
# ============================================================

python_exe = Path(sys.executable).resolve()

if ENV_ROOT not in python_exe.parents:
    fail(f"Wrong environment: {python_exe}")

df = pd.read_parquet(SOURCE)

if len(df) != EXPECTED_CELLS:
    fail(
        f"Cells={len(df):,}; expected={EXPECTED_CELLS:,}"
    )

if df["NAME"].nunique() != EXPECTED_CELLS:
    fail("NAME is not unique.")

if df["Subject"].nunique() != EXPECTED_DONORS:
    fail("Expected 30 subjects.")

if df["Sample"].nunique() != EXPECTED_SAMPLES:
    fail("Expected 133 Sample labels.")

if df["Cluster"].nunique() != EXPECTED_CLUSTERS:
    fail("Expected 51 author clusters.")

observed_clusters = set(df["Cluster"])
mapping_clusters = set(MAPPING)

missing_mapping = observed_clusters - mapping_clusters
extra_mapping = mapping_clusters - observed_clusters

if missing_mapping:
    fail(
        "Unmapped clusters: "
        + ", ".join(sorted(missing_mapping))
    )

if extra_mapping:
    fail(
        "Mapping contains absent clusters: "
        + ", ".join(sorted(extra_mapping))
    )


# ============================================================
# 1. Candidate pool
#
# All downloaded cells are author-QC cells belonging to the
# UC/healthy colon atlas and are eligible at sample-selection
# stage. Identity suitability is handled separately below.
# ============================================================

candidate = df.copy()

candidate["dataset"] = "SCP259"

candidate["donor_id"] = (
    candidate["Subject"].astype(str)
)

candidate["sample_id"] = (
    candidate["Sample"].astype(str)
)

candidate["author_annotation"] = (
    candidate["Cluster"].astype(str)
)

candidate["tissue_fraction"] = (
    candidate["Location"].map(
        {
            "Epi": "Epithelial_fraction",
            "LP": "Lamina_Propria",
        }
    )
)

candidate["site"] = "Colon"

candidate["disease"] = (
    candidate["Health"].map(
        {
            "Healthy": "Healthy",
            "Inflamed": "UC",
            "Non-inflamed": "UC",
        }
    )
)

candidate["inflammation"] = (
    candidate["Health"].map(
        {
            "Healthy": "Healthy",
            "Inflamed": "Inflamed",
            "Non-inflamed": "Non_Inflamed",
        }
    )
)

if candidate[
    [
        "tissue_fraction",
        "disease",
        "inflammation",
    ]
].isna().any().any():
    fail(
        "Canonical metadata mapping produced NA."
    )

candidate["candidate_status"] = "ELIGIBLE"


# ============================================================
# 2. Apply taxonomy mapping
# ============================================================

mapping_rows = []

for cluster in sorted(MAPPING):

    coarse, action, state, rationale = (
        MAPPING[cluster]
    )

    sub = candidate.loc[
        candidate["Cluster"].eq(cluster)
    ]

    mapping_rows.append(
        {
            "author_annotation": cluster,
            "n_cells": len(sub),
            "n_samples": sub["sample_id"].nunique(),
            "n_donors": sub["donor_id"].nunique(),
            "expression_bundles":
                " | ".join(
                    sorted(
                        sub[
                            "expression_bundle"
                        ].unique()
                    )
                ),
            "coarse_identity": coarse,
            "primary_reference_action": action,
            "state_flag": state,
            "mapping_rationale": rationale,
            "review_notes": "",
        }
    )


mapping = pd.DataFrame(
    mapping_rows
)

mapping.to_csv(
    OUT
    / "cluster_mapping_provisional_v1.tsv",
    sep="\t",
    index=False,
)


lookup = (
    mapping[
        [
            "author_annotation",
            "coarse_identity",
            "primary_reference_action",
            "state_flag",
        ]
    ]
    .set_index(
        "author_annotation"
    )
)


candidate = candidate.join(
    lookup,
    on="author_annotation",
)


if candidate[
    "coarse_identity"
].isna().any():
    fail(
        "Some candidate cells are unmapped."
    )


# ============================================================
# 3. Action counts
# ============================================================

action_counts = (
    candidate[
        "primary_reference_action"
    ]
    .value_counts()
)

expected_actions = {
    "KEEP": 341_126,
    "EXCLUDE_PRIMARY": 21_764,
    "HOLD_CROSS_DATASET": 2_602,
}

if action_counts.to_dict() != expected_actions:
    fail(
        "Unexpected action totals: "
        f"{action_counts.to_dict()}"
    )


# ============================================================
# 4. Identity coverage
# ============================================================

coverage_rows = []

for (
    identity,
    action
), sub in candidate.groupby(
    [
        "coarse_identity",
        "primary_reference_action",
    ],
    observed=True,
):

    donor_counts = (
        sub[
            "donor_id"
        ]
        .value_counts()
    )

    frac = (
        donor_counts
        / donor_counts.sum()
    )

    coverage_rows.append(
        {
            "coarse_identity": identity,
            "primary_reference_action": action,
            "n_cells": len(sub),
            "n_samples": sub["sample_id"].nunique(),
            "n_donors": sub["donor_id"].nunique(),
            "median_cells_per_donor":
                float(donor_counts.median()),
            "max_donor_fraction":
                float(frac.max()),
            "donor_HHI":
                float((frac ** 2).sum()),
        }
    )


coverage = (
    pd.DataFrame(
        coverage_rows
    )
    .sort_values(
        [
            "primary_reference_action",
            "n_cells",
        ],
        ascending=[
            True,
            False,
        ],
    )
)


coverage.to_csv(
    OUT
    / "identity_coverage.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 5. Confirm expected primary taxonomy for this colon dataset
# ============================================================

keep = candidate.loc[
    candidate[
        "primary_reference_action"
    ].eq(
        "KEEP"
    )
]

observed_primary = set(
    keep["coarse_identity"]
)

if observed_primary != EXPECTED_PRIMARY:

    fail(
        "Unexpected SCP259 primary identities.\n"
        f"Observed={sorted(observed_primary)}\n"
        f"Expected={sorted(EXPECTED_PRIMARY)}"
    )


# ============================================================
# 6. Condition coverage
# ============================================================

condition = (
    candidate.groupby(
        [
            "coarse_identity",
            "primary_reference_action",
            "disease",
            "inflammation",
        ],
        observed=True,
    )
    .agg(
        n_cells=(
            "NAME",
            "size",
        ),
        n_samples=(
            "sample_id",
            "nunique",
        ),
        n_donors=(
            "donor_id",
            "nunique",
        ),
    )
    .reset_index()
)


condition.to_csv(
    OUT
    / "identity_by_condition.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 7. Save candidate metadata
# ============================================================

candidate.to_parquet(
    OUT
    / "scp259_candidate_cells_provisional_taxonomy.parquet",
    index=False,
    engine="pyarrow",
)


# ============================================================
# 8. Console
# ============================================================

print()
print("=" * 78)
print("SCP259 PHASE 1-3 RESULTS")
print("=" * 78)

print()
print(
    f"Candidate cells      : {len(candidate):,}"
)
print(
    f"Subjects/donors      : {candidate['donor_id'].nunique()}"
)
print(
    f"Sample labels        : {candidate['sample_id'].nunique()}"
)
print(
    f"Author clusters      : {candidate['author_annotation'].nunique()}"
)
print(
    f"Primary identities   : {len(observed_primary)}"
)

print()
print("Mapping actions:")
print(
    action_counts.to_string()
)

print()
print("Identity coverage:")
print(
    coverage.to_string(
        index=False
    )
)

print()
print("Final primary identities:")
for x in sorted(observed_primary):
    print("  -", x)

print()
print("PHASE 1-3 STATUS: PASS")
print(
    f"Outputs: {OUT}"
)
