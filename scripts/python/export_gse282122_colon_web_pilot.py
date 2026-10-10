#!/usr/bin/env python3
"""Candidate colon web input, using the project's existing GSE282122 mapping.

Read-only source; bounded selected CSR reads; no clustering or clinical analysis.
Requires adjacent ibdverse_pilot_v2 helpers on PYTHONPATH.
"""
import argparse
import gzip
import hashlib
import json
import shutil
import sys
import tarfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from inspect_ibdverse_h5ad import read_column
from plan_ibdverse_epithelial_pilot import stable_rng
from export_ibdverse_cibersortx_pilot import (
    balanced_pick, selected_csr, write_qc_summary, MARKERS, sha256, dump_json,
)

VERSION = "20261006_gse282122_colon_candidate_v1"
TYPES = ["Colonic_absorptive", "Stem_TA_progenitor", "Goblet", "Enteroendocrine",
         "B", "Plasma", "CD4_T", "CD8_T", "NK", "ILC", "Monocyte_Macrophage",
         "DC", "Mast", "Fibroblast", "Pericyte", "Endothelial"]
SITES = ["Ascending_Colon", "Descending_Colon", "Sigmoid"]
QC = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]


def run(args):
    started = time.monotonic()
    source = Path(args.input).resolve(strict=True)
    candidate = Path(args.candidate).resolve(strict=True)
    out = Path(args.outdir).resolve()
    if out.exists() or out == source.parent or out == source:
        raise ValueError("Use a new output directory; refuse overwrite")
    if args.cells_per_type < 20 or args.cells_per_type % 2:
        raise ValueError("Even --cells-per-type >=20 required")
    out.mkdir(parents=True)
    dump_json(out / "run_parameters.json", {"version": VERSION, **vars(args),
              "python": sys.executable, "pandas": pd.__version__, "numpy": np.__version__,
              "source_bytes": source.stat().st_size, "candidate_sha256": sha256(candidate)})
    c = pd.read_parquet(candidate)
    assert c.obs_name.is_unique and c.cell_id.is_unique
    assert c.reference_eligible.eq(True).all()
    assert (c.Treatment.eq("Pre") | c.Disease.eq("Healthy")).all()
    assert c.Patient.notna().all() and c.Disease.isin(["CD", "UC", "Healthy"]).all()
    donors = c[["Patient", "Disease"]].drop_duplicates().sort_values("Patient").reset_index(drop=True)
    assert donors.Patient.is_unique, "Donor disease conflict"
    donors["split"] = "training"
    for disease in ["CD", "UC", "Healthy"]:
        positions = donors.index[donors.Disease.eq(disease)].to_numpy()
        n = max(1, int(round(len(positions) * .2)))
        assert n < len(positions)
        chosen = stable_rng("GSE282122_donor_holdout", disease, args.seed).choice(positions, n, replace=False)
        donors.loc[chosen, "split"] = "holdout"
    donors.to_csv(out / "donor_split.tsv", sep="\t", index=False)
    c["split"] = c.Patient.map(donors.set_index("Patient").split)
    colon = c[c.Site.isin(SITES)].copy()
    colon.groupby(["coarse_identity", "primary_reference_action", "Disease"], observed=True).agg(
        n_cells=("sample_id", "size"), n_donors=("Patient", "nunique")).reset_index().to_csv(
        out / "candidate_coverage.tsv", sep="\t", index=False)
    keep = colon.primary_reference_action.eq("KEEP") & colon.coarse_identity.isin(TYPES)
    keep &= colon.identity_status.eq("CONFIDENT")
    for field in QC:
        v = pd.to_numeric(colon[field], errors="coerce")
        keep &= np.isfinite(v)
        if field != "pct_counts_mt":
            keep &= v > 0
    assert pd.api.types.is_bool_dtype(colon.predicted_doublets.dtype), "Unexpected doublet encoding"
    keep &= ~colon.predicted_doublets
    pool = colon[keep & colon.split.eq("training")].copy()
    pool["individual_id"] = pool.Patient
    with h5py.File(source, "r") as h:
        ids = read_column(h["obs"][h["obs"].attrs.get("_index", "_index")]).astype(str)
        assert ids.is_unique
        positions = pd.Index(ids).get_indexer(pool.obs_name.astype(str))
        assert (positions >= 0).all(), "Candidate IDs absent from source"
        pool["source_cell_row"] = positions
        # Verify cached metadata against current source, without reading expression.
        for field in ["sample_id", "Patient", "Disease", "Site", "final_analysis"]:
            now = read_column(h["obs"][field]).iloc[positions].astype(str).to_numpy()
            assert np.array_equal(now, pool[field].astype(str).to_numpy()), field + " metadata mismatch"
        healthy_caps = []
        for name in TYPES:
            q = pool[pool.coarse_identity.eq(name) & pool.Disease.eq("Healthy")]
            healthy_caps.append(int(q.groupby("Patient", observed=True).size().clip(upper=5).sum()))
        healthy_target = min(10, min(healthy_caps))
        healthy_target -= healthy_target % 2
        if healthy_target < 4:
            raise ValueError("Healthy donor capacity <4 per class; review candidate scheme")
        case_target = (args.cells_per_type - healthy_target) // 2
        quotas = {"Healthy": healthy_target, "CD": case_target, "UC": case_target}
        rows, capacity = [], []
        for name in TYPES:
            for disease, target in quotas.items():
                q = pool[pool.coarse_identity.eq(name) & pool.Disease.eq(disease)]
                chosen, status = balanced_pick(q, target, name, disease, args.seed, 5)
                rows.extend(chosen)
                capacity.append(status)
        pd.DataFrame(capacity).to_csv(out / "sampling_capacity.tsv", sep="\t", index=False)
        dump_json(out / "pilot_decisions.json", {"types": TYPES, "sites": SITES,
                  "cells_per_type": args.cells_per_type, "common_disease_quotas": quotas,
                  "healthy_capacity_per_type": dict(zip(TYPES, healthy_caps)),
                  "reason": "candidate pilot; 3 healthy donors; common quota capped by training-donor capacity",
                  "excluded": {"Ileal_absorptive": "ileum-specific", "Paneth": "colonic Paneth requires separate metaplasia/identity review",
                  "Tuft_Glial_Unconventional_T_M_like": "existing mapping HOLD_CROSS_DATASET; not promoted in this candidate pilot"},
                  "BEST4": "follows existing mapping, may be inside Colonic_absorptive; no separate class asserted",
                  "final_reference_dataset_selected": False})
        if any(s["status"] != "ok" for s in capacity):
            raise ValueError("Sampling quota shortage; inspect sampling_capacity.tsv")
        rows = np.sort(np.asarray(rows, dtype=np.int64))
        assert len(rows) == args.cells_per_type * len(TYPES) and len(np.unique(rows)) == len(rows)
        selected = pool.set_index("source_cell_row", drop=False).loc[rows].reset_index(drop=True)
        assert not selected.Patient.isin(donors.loc[donors.split.eq("holdout"), "Patient"]).any()
        selected["web_column"] = ""
        for name in TYPES:
            ix = selected.index[selected.coarse_identity.eq(name)]
            selected.loc[ix, "web_column"] = [f"{name}.{i+1}" for i in range(len(ix))]
        genes = pd.DataFrame({"source_gene_id": read_column(h["var"]["gene_id"]).astype("string"),
                              "gene_symbol": read_column(h["var"]["gene_symbol"]).astype("string")})
        assert genes.source_gene_id.is_unique
        print(f"Reading {len(rows)} colon cells; {len(TYPES)} types; disease quotas {quotas}", flush=True)
        indices, values, pointers, totals, detected = selected_csr(h["X"], rows, len(genes))
    assert np.allclose(totals, selected.total_counts.to_numpy(), rtol=0, atol=.01), "UMI totals mismatch"
    assert np.array_equal(detected, selected.n_genes_by_counts.to_numpy()), "Gene detection mismatch"
    selected["extracted_total_counts"] = totals
    selected["extracted_n_genes"] = detected
    selected.to_csv(out / "selected_cell_manifest.tsv", sep="\t", index=False)
    selected.groupby(["coarse_identity", "Disease", "Site", "Inflammation"], observed=True).size().rename(
        "n_cells").reset_index().to_csv(out / "selected_composition.tsv", sep="\t", index=False)
    write_qc_summary(selected, ["coarse_identity", "Disease"], QC, out / "selected_qc_summary.tsv")
    cell = np.repeat(np.arange(len(rows), dtype=np.int32), np.diff(pointers))
    symbols = genes.gene_symbol.fillna("")
    valid_symbol = symbols.str.len().gt(0) & ~symbols.str.contains(r"[\t\r\n]")
    ambiguous = symbols.duplicated(keep=False)
    sums = np.bincount(indices, weights=values, minlength=len(genes))
    gene_keep = (valid_symbol & ~ambiguous).to_numpy() & (sums > 0)
    genes["ambiguous_symbol"] = ambiguous
    genes["web_included"] = gene_keep
    genes.to_csv(out / "gene_export_manifest.tsv", sep="\t", index=False)
    with h5py.File(out / "selected_counts.h5", "w") as t:
        t.attrs["format"] = "GSE282122 colon pilot custom CSR; not AnnData h5ad"
        g = t.create_group("counts"); g.attrs["encoding-type"] = "csr_matrix"; g.attrs["shape"] = [len(rows), len(genes)]
        for name, data in [("data", values), ("indices", indices), ("indptr", pointers)]:
            g.create_dataset(name, data=data, compression="gzip", compression_opts=1)
        strings = h5py.string_dtype("utf-8")
        for name, data in [("gene_symbol", genes.gene_symbol), ("source_gene_id", genes.source_gene_id),
                           ("cell_id", selected.obs_name), ("analysis_group", selected.coarse_identity),
                           ("web_column", selected.web_column)]:
            t.create_dataset(name, data=np.asarray(data.fillna("").astype(str), dtype=object), dtype=strings)
        t.create_dataset("source_cell_row", data=rows)
    marker_records = []
    groups = selected.groupby("coarse_identity", observed=True).groups
    for symbol in MARKERS + ["CA1", "CA2", "AQP8", "GUCA2A", "SLC26A3"]:
        for gene in genes.index[genes.gene_symbol.eq(symbol).fillna(False)]:
            counts = np.bincount(cell[indices == gene], weights=values[indices == gene], minlength=len(rows))
            for name, ix in groups.items():
                ix = np.asarray(ix)
                marker_records.append({"coarse_identity": name, "gene_symbol": symbol,
                    "source_gene_id": str(genes.at[gene, "source_gene_id"]), "web_included": bool(gene_keep[gene]),
                    "n_cells": len(ix), "detection_fraction": float((counts[ix] > 0).mean()),
                    "mean_cpm_original_library": float((counts[ix] / totals[ix] * 1e6).mean())})
    pd.DataFrame(marker_records).to_csv(out / "marker_expression_summary.tsv", sep="\t", index=False)
    use = gene_keep[indices]
    lookup = np.full(len(genes), -1, dtype=np.int32); lookup[gene_keep] = np.arange(gene_keep.sum())
    gi, ci, vv = lookup[indices[use]], cell[use], values[use]
    assert (np.bincount(ci, weights=vv, minlength=len(rows)) > 0).all()
    order = np.argsort(gi, kind="stable"); gi, ci, vv = gi[order], ci[order], vv[order]
    ptr = np.r_[0, np.cumsum(np.bincount(gi, minlength=int(gene_keep.sum())))]
    names = genes.loc[gene_keep, "gene_symbol"].astype(str).to_numpy()
    web = out / "GSE282122_Colon_16types_scRNA_reference.txt"
    with web.open("w", encoding="utf-8", newline="") as f:
        f.write("GeneSymbol\t" + "\t".join(selected.web_column) + "\n")
        for begin in range(0, len(names), 128):
            end = min(begin + 128, len(names))
            block = np.zeros((end - begin, len(rows)), dtype=np.int32)
            left, right = int(ptr[begin]), int(ptr[end])
            np.add.at(block, (gi[left:right] - begin, ci[left:right]), vv[left:right])
            for name, row in zip(names[begin:end], block):
                f.write(name + "\t" + "\t".join(map(str, row)) + "\n")
    with web.open("rb") as r, gzip.open(str(web) + ".gz", "wb", compresslevel=1) as w:
        shutil.copyfileobj(r, w)
    result = {"version": VERSION, "status": "completed", "source_opened_read_only": True,
              "web_reference_cells": len(rows), "web_reference_genes": len(names), "types": TYPES,
              "common_disease_quotas": quotas, "source_sites": SITES, "includes_rectum": False,
              "holdout_donors_used": 0, "training_donors": int(donors.split.eq("training").sum()),
              "holdout_donors": int(donors.split.eq("holdout").sum()), "actual_web_donors": selected.Patient.nunique(),
              "web_unit": "raw nonlog integer UMI counts from X; every extracted row verified against total_counts",
              "max_web_count": int(vv.max()), "web_format_valid": bool(vv.max() >= 50),
              "max_per_donor_per_type": int(selected.groupby(["coarse_identity", "Patient"], observed=True).size().max()),
              "web_input_bytes": web.stat().st_size, "web_input_sha256": sha256(web),
              "seconds_elapsed": round(time.monotonic() - started, 2), "is_finished_signature_matrix": False,
              "biological_identity_review_complete": False, "final_reference_dataset_selected": False}
    dump_json(out / "summary.json", result)
    with tarfile.open(out / "small_audit_reports.tar.gz", "w:gz") as t:
        for p in sorted(out.iterdir()):
            if p.suffix in {".json", ".tsv"}: t.add(p, arcname=p.name, recursive=False)
    print(json.dumps(result, indent=2), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--candidate", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--cells-per-type", type=int, default=100)
    p.add_argument("--seed", type=int, default=20261006)
    run(p.parse_args())


if __name__ == "__main__":
    main()
