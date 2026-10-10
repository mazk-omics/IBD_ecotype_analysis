"""Plan (not extract) epithelial or whole-mucosa pilots from aggregate audit reports."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


TYPE_ORDER = ["Ileal_absorptive", "Stem_TA_progenitor", "Goblet",
              "Enteroendocrine", "Tuft", "BEST4"]
WHOLE_TYPE_ORDER = TYPE_ORDER + ["B", "Plasma", "CD4_T", "CD8_T", "T_gamma_delta",
                                "ILC", "Monocyte_Macrophage", "DC", "Mast",
                                "Fibroblast", "Pericyte", "Endothelial"]
STRATA = {"Healthy_uninflamed": 100, "CD_inflamed": 50, "CD_uninflamed": 50}
PENDING_T_LABELS = {"T cells CD4- CD8-", "T cells CD4- CD8-lo"}


def map_type(label):
    if label == "Enterocytes BEST4+":
        return "BEST4"
    if label.startswith("Stem cells ") or label.startswith("Enterocytes progenitor crypt "):
        return "Stem_TA_progenitor"
    if label.startswith("Enterocytes "):
        return "Ileal_absorptive"
    if label.startswith("Goblet cells "):
        return "Goblet"
    if label in {"Enteroendocrine cells", "Enterochromaffin cells"}:
        return "Enteroendocrine"
    if label == "Tuft cells":
        return "Tuft"
    return None


def stable_rng(type_name, stratum, seed):
    # Do not use Python hash(), which changes across processes.
    digest = hashlib.sha256(f"{seed}|{type_name}|{stratum}".encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "big"))


def map_whole_type(label):
    epithelial = map_type(label)
    if epithelial is not None:
        return epithelial
    if label.startswith("B cells"):
        return "B"
    if label.startswith("Plasma B cell "):
        return "Plasma"
    if label.startswith("T cells CD4+"):
        return "CD4_T"
    if label.startswith("T cells CD8+"):
        return "CD8_T"
    if label == "T cells gd":
        return "T_gamma_delta"
    if label.startswith("Macrophages ") or label.startswith("Monocytes "):
        return "Monocyte_Macrophage"
    return {"ILCs": "ILC", "Mast cells": "Mast", "cDC1": "DC", "pDCs/cDC2": "DC",
            "Fibroblasts": "Fibroblast", "Pericyte cells": "Pericyte",
            "Endothelial cells": "Endothelial"}.get(label)


def run(args):
    source, dest = Path(args.reports), Path(args.outdir)
    mapper = map_whole_type if args.scope == "whole_mucosa" else map_type
    type_order = WHOLE_TYPE_ORDER if args.scope == "whole_mucosa" else TYPE_ORDER
    summary = json.loads((source / "summary.json").read_text(encoding="utf-8"))
    table = pd.read_csv(source / "label_category_sample_donor_coverage.tsv", sep="\t")
    assert int(table.n_cells.sum()) == summary["n_cells"]
    complete = table.reference_metadata_complete.astype(str).str.lower().eq("true")
    assert int(table.loc[complete, "n_cells"].sum()) == summary["cells_with_complete_reference_metadata"]
    labels = table[["label", "category"]].drop_duplicates().sort_values("label")
    labels["pilot_type"] = labels.label.map(mapper)
    author_epi = labels.category.isin(["Enterocyte", "Secretory", "Stem cell"])
    if args.scope == "whole_mucosa":
        assert set(labels.loc[labels.pilot_type.isna(), "label"]) == PENDING_T_LABELS, "Unexpected unmapped author labels"
    else:
        assert labels.pilot_type.notna().equals(author_epi), "Unmapped or unexpected epithelial labels"
    labels["mapping_status"] = np.where(labels.pilot_type.notna(), "provisional_pilot", "outside_epithelial_pilot")
    labels["note"] = ""
    labels.loc[labels.pilot_type.eq("Stem_TA_progenitor"), "note"] = "Crypt progenitors and stem labels provisionally merged; verify markers after extraction"
    labels.loc[labels.pilot_type.eq("Ileal_absorptive"), "note"] = "Crypt KRT20 precursors and villus labels merged; IFI27 states not separate types"
    labels.loc[labels.pilot_type.eq("Enteroendocrine"), "note"] = "Enterochromaffin and other endocrine labels provisionally merged"
    if args.scope == "whole_mucosa":
        labels.loc[labels.label.isin(PENDING_T_LABELS), "mapping_status"] = "pending_identity_review_not_a_signature_type"
        labels.loc[labels.label.isin(PENDING_T_LABELS), "note"] = "Keep for identity review; no pooled T_other signature and no forced CD4/CD8/NK assignment"
        labels.loc[labels.pilot_type.eq("T_gamma_delta"), "note"] = "Separate candidate from author gamma-delta label; verify expression before production"
    labels.loc[labels.label.eq("B cells germinal centre/plasmablasts"), "note"] = "Follows author B-cell category for pilot; verify B/plasma identity before production"
    labels.loc[labels.label.eq("T cells CD8+ TRGC2+"), "note"] = "Follows explicit author CD8 label; inspect gamma-delta/cytotoxic identity before production"
    pool = table.loc[complete & table.biopsy_type.eq("terminal ileum")].copy()
    pool["pilot_type"] = pool.label.map(mapper)
    pending = pool.loc[pool.pilot_type.isna()].copy() if args.scope == "whole_mucosa" else pool.iloc[:0].copy()
    pending_support = pending.groupby("label", as_index=False).agg(
        n_cells=("n_cells", "sum"), n_donors=("individual_id", "nunique"),
        n_samples=("sanger_sample_id", "nunique"))
    pending_support["status"] = "pending_identity_review_not_a_signature_type"
    pool = pool.loc[pool.pilot_type.notna()].copy()
    if args.scope == "whole_mucosa":
        assert int(pool.n_cells.sum() + pending.n_cells.sum()) == summary["cells_with_complete_reference_metadata"]
    pool["stratum"] = pool.disease_status + "_" + pool.ti_ses_inflamed
    assert set(pool.stratum) == set(STRATA)
    assert pool.individual_id.notna().all()
    assert pool.groupby("individual_id").stratum.nunique().max() == 1
    donors = (pool.groupby(["pilot_type", "stratum", "individual_id"], as_index=False)
              .agg(available_cells=("n_cells", "sum"), n_samples=("sanger_sample_id", "nunique")))
    donors["cap5_capacity"] = donors.available_cells.clip(upper=5)
    total_rows, strata_rows, quotas = [], [], []
    for type_name in type_order:
        sub = pool.loc[pool.pilot_type.eq(type_name)]
        total_rows.append({"pilot_type": type_name, "available_cells": int(sub.n_cells.sum()),
                           "n_donors": sub.individual_id.nunique(), "n_samples": sub.sanger_sample_id.nunique(),
                           "n_author_labels": sub.label.nunique(), "target_cells": 200})
        for stratum, target in STRATA.items():
            ds = donors.loc[donors.pilot_type.eq(type_name) & donors.stratum.eq(stratum)].copy()
            ds = ds.sort_values("individual_id").reset_index(drop=True)
            capacity = int(ds.cap5_capacity.sum())
            if capacity < target:
                raise ValueError(f"{type_name}/{stratum}: cap5 capacity {capacity} < target {target}")
            order = stable_rng(type_name, stratum, args.seed).permutation(len(ds))
            selected = np.zeros(len(ds), dtype=int)
            remaining = target
            for _ in range(5):
                for idx in order:
                    if selected[idx] < ds.at[idx, "cap5_capacity"]:
                        selected[idx] += 1
                        remaining -= 1
                        if remaining == 0:
                            break
                if remaining == 0:
                    break
            assert remaining == 0
            strata_rows.append({"pilot_type": type_name, "stratum": stratum,
                                "available_cells": int(ds.available_cells.sum()), "n_donors": len(ds),
                                "capacity_at_max5_per_donor": capacity, "target_cells": target,
                                "planned_donors": int((selected > 0).sum()),
                                "max_planned_cells_per_donor": int(selected.max())})
            for idx in np.flatnonzero(selected):
                row = ds.loc[idx]
                samples = sub.loc[sub.individual_id.eq(row.individual_id), "sanger_sample_id"].unique()
                quotas.append({"pilot_type": type_name, "stratum": stratum,
                               "individual_id": row.individual_id, "sanger_sample_ids": "|".join(sorted(samples)),
                               "available_cells": int(row.available_cells), "planned_cells": int(selected[idx]),
                               "status": "metadata_only_draft_no_cell_ids_selected"})
    totals = pd.DataFrame(total_rows)
    coverage = pd.DataFrame(strata_rows)
    quota = pd.DataFrame(quotas)
    target_total = 200 * len(type_order)
    assert int(quota.planned_cells.sum()) == target_total
    assert quota.groupby("pilot_type").planned_cells.sum().eq(200).all()
    assert quota.groupby(["pilot_type", "individual_id"]).planned_cells.sum().max() <= 5
    gene_dups = pd.read_csv(source / "duplicate_gene_symbols.tsv", sep="\t")
    assert len(gene_dups) == 48 and gene_dups.gene_symbol.nunique() == 24
    dest.mkdir(parents=True, exist_ok=True)
    for name, frame in {"author_label_to_pilot_type.tsv": labels,
                        "pilot_type_support.tsv": totals,
                        "pilot_type_stratum_support.tsv": coverage,
                        "draft_donor_sampling_quotas.tsv": quota,
                        "pending_label_support.tsv": pending_support,
                        "candidate_genes_to_exclude_ambiguous_symbols.tsv": gene_dups}.items():
        frame.to_csv(dest / name, sep="\t", index=False, encoding="utf-8")
    result = {"scope": "planning from aggregate metadata; no expression extraction or cell-level QC",
              "source_report_directory": str(source.resolve()), "seed": args.seed,
              "pilot_scope": args.scope, "candidate_cells": int(pool.n_cells.sum()),
              "mapping_revision": "20261006_no_pooled_T_other",
              "labels_pending_identity_review": sorted(pending.label.unique()),
              "cells_pending_identity_review": int(pending.n_cells.sum()),
              "ready_for_web_export": False,
              "candidate_donors": pool.individual_id.nunique(), "pilot_types": type_order,
              "target_cells_per_type": 200, "target_total_cells": target_total,
              "stratum_targets_per_type": STRATA, "max_cells_per_donor_per_type": 5,
              "draft_plan_feasible_before_cell_level_checks": True,
              "actual_cells_selected": False, "planned_distinct_donors": quota.individual_id.nunique(),
              "max_actual_draft_quota_per_donor_per_type": int(quota.groupby(["pilot_type", "individual_id"]).planned_cells.sum().max()),
              "conflicting_sample_excluded": "OTARscRNA14592808",
              "gene_handling_proposal": "Exclude both/all entries of ambiguous symbol groups in pilot; no arbitrary first locus or suffix",
              "report_sha256": hashlib.sha256((source / "label_category_sample_donor_coverage.tsv").read_bytes()).hexdigest()}
    (dest / "plan_summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(totals.to_string(index=False))
    print(coverage.to_string(index=False))
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=20261006)
    parser.add_argument("--scope", choices=["epithelial", "whole_mucosa"], default="epithelial")
    run(parser.parse_args())
