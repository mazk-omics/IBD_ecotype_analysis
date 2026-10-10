#!/usr/bin/env python3
"""Small donor-balanced GSE282122 joint/ileum and SCP259 reference inputs.

Uses existing audited labels and source files. No source mutation, clustering,
clinical outcome analysis or expression integration. Requires existing helpers.
"""
import argparse
import gzip
import json
import shutil
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

VERSION = "20261006_three_reference_pilots_v1"
ONTOLOGY = ["B", "Plasma", "CD4_T", "CD8_T", "NK", "ILC", "Monocyte_Macrophage",
            "DC", "Mast", "Fibroblast", "Pericyte", "Endothelial", "Ileal_absorptive",
            "Colonic_absorptive", "Goblet", "Stem_TA_progenitor", "Paneth",
            "Enteroendocrine", "Glial", "Tuft"]
GSE_MODES = {"gse_joint": ONTOLOGY, "gse_ileum": [t for t in ONTOLOGY if t != "Colonic_absorptive"]}
SCP_TYPES = [t for t in ONTOLOGY if t not in {"Ileal_absorptive", "Paneth"}]
BUNDLES = {
    "Epithelial": ("5cdc540d328cee7a2efc2348", "Epi"),
    "Stromal": ("5cdc540d328cee7a2efc2349", "Fib"),
    "Immune": ("5cdc540d328cee7a2efc234a", "Imm"),
}


def fresh(path):
    p = Path(path).resolve()
    if p.exists():
        raise ValueError("Use a new output directory; refuse overwrite: " + str(p))
    p.mkdir(parents=True)
    return p


def split_scp_donors(frame, seed):
    d = frame[["individual_id", "Disease"]].drop_duplicates().sort_values("individual_id").reset_index(drop=True)
    assert d.individual_id.is_unique, "Donor disease conflict"
    d["split"] = "training"
    for disease in sorted(d.Disease.unique()):
        ix = d.index[d.Disease.eq(disease)].to_numpy()
        assert len(ix) >= 3
        n = max(1, int(round(len(ix) * .2)))
        chosen = stable_rng("SCP259_donor_holdout", disease, seed).choice(ix, n, replace=False)
        d.loc[chosen, "split"] = "holdout"
    return d


def choose_cells(pool, types, target, seed, out):
    pool = pool.reset_index(drop=True)
    pool["source_cell_row"] = np.arange(len(pool))
    rows, capacity = [], []
    for name in types:
        q = pool[pool.pilot_type.eq(name)]
        picked, status = balanced_pick(q, target, name, "all_training_donors", seed, 5)
        rows.extend(picked); capacity.append(status)
    pd.DataFrame(capacity).to_csv(out / "sampling_capacity.tsv", sep="\t", index=False)
    if any(x["status"] != "ok" for x in capacity):
        raise ValueError("Sampling shortage; inspect sampling_capacity.tsv; no quota altered")
    selected = pool.iloc[np.sort(rows)].copy().reset_index(drop=True)
    assert len(selected) == len(types) * target and selected.cell_id.is_unique
    selected["web_column"] = ""
    for name in types:
        ix = selected.index[selected.pilot_type.eq(name)]
        selected.loc[ix, "web_column"] = [f"{name}.{j+1}" for j in range(len(ix))]
    return selected


def write_outputs(out, selected, genes, csr, types, mode, decisions, started):
    indices, values, pointers, totals, detected = csr
    n_cells, n_genes = len(selected), len(genes)
    selected["extracted_total_counts"] = totals
    selected["extracted_n_genes"] = detected
    selected["total_counts_matches_obs"] = np.isclose(totals, selected.total_counts.to_numpy(), rtol=0, atol=.01)
    selected["n_genes_matches_obs"] = detected == selected.n_genes_by_counts.to_numpy()
    if mode.startswith("gse"):
        assert selected.total_counts_matches_obs.all() and selected.n_genes_matches_obs.all(), "Extracted X/OBS counts mismatch"
    selected.to_csv(out / "selected_cell_manifest.tsv", sep="\t", index=False)
    selected.groupby(["pilot_type", "Disease", "Site", "Inflammation"], observed=True).size().rename(
        "n_cells").reset_index().to_csv(out / "selected_composition.tsv", sep="\t", index=False)
    selected.groupby(["pilot_type", "author_label"], observed=True).size().rename("n_cells").reset_index().to_csv(
        out / "selected_author_label_composition.tsv", sep="\t", index=False)
    write_qc_summary(selected, ["pilot_type"], ["total_counts", "n_genes_by_counts"], out / "selected_qc_summary.tsv")
    cell = np.repeat(np.arange(n_cells, dtype=np.int32), np.diff(pointers))
    symbol = genes.gene_symbol.fillna("").astype("string")
    valid = symbol.str.len().gt(0) & ~symbol.str.contains(r"[\t\r\n]")
    duplicate = symbol.duplicated(keep=False)
    sums = np.bincount(indices, weights=values, minlength=n_genes)
    keep = (valid & ~duplicate).to_numpy() & (sums > 0)
    genes["ambiguous_symbol"] = duplicate; genes["web_included"] = keep
    genes.to_csv(out / "gene_export_manifest.tsv", sep="\t", index=False)
    with h5py.File(out / "selected_counts.h5", "w") as h:
        h.attrs["format"] = "Custom selected-cell CSR; not AnnData h5ad"
        c = h.create_group("counts"); c.attrs["shape"] = [n_cells, n_genes]; c.attrs["encoding-type"] = "csr_matrix"
        for name, data in [("indices", indices), ("data", values), ("indptr", pointers)]:
            c.create_dataset(name, data=data, compression="gzip", compression_opts=1)
        strings = h5py.string_dtype("utf-8")
        for name, data in [("gene_symbol", symbol), ("source_gene_id", genes.source_gene_id),
                           ("cell_id", selected.cell_id), ("pilot_type", selected.pilot_type),
                           ("web_column", selected.web_column)]:
            h.create_dataset(name, data=np.asarray(data.fillna("").astype(str), dtype=object), dtype=strings)
        h.create_dataset("source_cell_row", data=selected.source_cell_row.to_numpy())
    records = []
    groups = selected.groupby("pilot_type", observed=True).groups
    for marker in MARKERS + ["CA1", "CA2", "AQP8", "SLC26A3"]:
        for gene in genes.index[symbol.eq(marker).fillna(False)]:
            counts = np.bincount(cell[indices == gene], weights=values[indices == gene], minlength=n_cells)
            for name, ix in groups.items():
                ix = np.asarray(ix)
                records.append({"pilot_type": name, "gene_symbol": marker, "source_gene_id": str(genes.at[gene, "source_gene_id"]),
                    "web_included": bool(keep[gene]), "n_cells": len(ix), "detection_fraction": float((counts[ix] > 0).mean()),
                    "mean_cpm_original_library": float((counts[ix] / totals[ix] * 1e6).mean())})
    pd.DataFrame(records).to_csv(out / "marker_expression_summary.tsv", sep="\t", index=False)
    use = keep[indices]
    lookup = np.full(n_genes, -1, dtype=np.int32); lookup[keep] = np.arange(keep.sum())
    gi, ci, vv = lookup[indices[use]], cell[use], values[use]
    assert (np.bincount(ci, weights=vv, minlength=n_cells) > 0).all()
    order = np.argsort(gi, kind="stable"); gi, ci, vv = gi[order], ci[order], vv[order]
    gp = np.r_[0, np.cumsum(np.bincount(gi, minlength=int(keep.sum())))]
    names = symbol[keep].to_numpy()
    filename = {"gse_joint": "GSE282122_Joint20_scRNA_reference.txt",
                "gse_ileum": "GSE282122_Ileum19_scRNA_reference.txt",
                "scp259": "SCP259_ColonRectum18_scRNA_reference.txt"}[mode]
    web = out / filename
    print(f"Writing {len(names)} genes x {n_cells} cells: {filename}", flush=True)
    with web.open("w", encoding="utf-8", newline="") as f:
        f.write("GeneSymbol\t" + "\t".join(selected.web_column) + "\n")
        for begin in range(0, len(names), 128):
            end = min(begin+128, len(names)); left, right = map(int, gp[[begin, end]])
            block = np.zeros((end-begin, n_cells), dtype=np.int32)
            np.add.at(block, (gi[left:right]-begin, ci[left:right]), vv[left:right])
            for name, row in zip(names[begin:end], block):
                f.write(name + "\t" + "\t".join(map(str, row)) + "\n")
    with web.open("rb") as r, gzip.open(str(web)+".gz", "wb", compresslevel=1) as w:
        shutil.copyfileobj(r, w)
    result = {"version": VERSION, "status": "completed", "mode": mode,
        "source_opened_read_only": True, "web_reference_types": types,
        "web_reference_cells": n_cells, "web_reference_genes": len(names),
        "cells_per_type": int(selected.groupby("pilot_type").size().min()),
        "actual_web_donors": int(selected.individual_id.nunique()), "holdout_donors_used": 0,
        "max_per_donor_per_type": int(selected.groupby(["pilot_type", "individual_id"]).size().max()),
        "max_web_count": int(vv.max()), "web_format_valid": bool(vv.max() >= 50),
        "web_unit": "nonlog raw integer counts", "web_input_filename": filename,
        "web_input_bytes": web.stat().st_size, "web_input_sha256": sha256(web),
        "total_counts_matching_obs": int(selected.total_counts_matches_obs.sum()),
        "n_genes_matching_obs": int(selected.n_genes_matches_obs.sum()),
        "seconds_elapsed": round(time.monotonic()-started, 2), "is_finished_signature_matrix": False,
        "biological_identity_review_complete": False, **decisions}
    dump_json(out / "summary.json", result)
    with tarfile.open(out / "small_audit_reports.tar.gz", "w:gz") as t:
        for p in sorted(out.iterdir()):
            if p.suffix in {".json", ".tsv"}: t.add(p, arcname=p.name, recursive=False)
    print(json.dumps(result, indent=2), flush=True)
    return result


def run_gse(args):
    started = time.monotonic(); out = fresh(args.outdir)
    dump_json(out / "run_parameters.json", {"version": VERSION, **vars(args)})
    c = pd.read_parquet(args.candidate)
    d = pd.read_csv(args.donor_split, sep="\t")
    assert d.Patient.is_unique and c.obs_name.is_unique and c.cell_id.is_unique
    assert c.reference_eligible.eq(True).all() and (c.Treatment.eq("Pre") | c.Disease.eq("Healthy")).all()
    assert c.Patient.isin(d.Patient).all(), "Cached donor split missing candidate donors"
    assert c[["Patient", "Disease"]].astype(str).drop_duplicates().set_index("Patient").Disease.to_dict() == (
        d[["Patient", "Disease"]].astype(str).set_index("Patient").Disease.to_dict()), "Donor split/candidate disease mismatch"
    c["split"] = c.Patient.map(d.set_index("Patient").split)
    c["pilot_type"] = c.coarse_identity; c["individual_id"] = c.Patient; c["author_label"] = c.final_analysis
    promoted = c.pilot_type.isin(["Glial", "Tuft"]) & c.primary_reference_action.eq("HOLD_CROSS_DATASET")
    action = c.primary_reference_action.eq("KEEP") | promoted
    c["included_by_this_pilot"] = action & c.identity_status.eq("CONFIDENT")
    c.groupby(["pilot_type", "primary_reference_action", "included_by_this_pilot"]).size().rename("n_cells").reset_index().to_csv(
        out / "taxonomy_actions.tsv", sep="\t", index=False)
    types = GSE_MODES[args.mode]
    valid = c.included_by_this_pilot & c.pilot_type.isin(types) & c.split.eq("training")
    assert pd.api.types.is_bool_dtype(c.predicted_doublets)
    valid &= ~c.predicted_doublets
    for field in ["n_genes_by_counts", "total_counts"]:
        valid &= np.isfinite(pd.to_numeric(c[field], errors="coerce")) & c[field].gt(0)
    if args.mode == "gse_ileum": valid &= c.Site.eq("Terminal_Ileum")
    pool = c[valid].copy().sort_values("cell_id")
    d.to_csv(out / "donor_split.tsv", sep="\t", index=False)
    selected = choose_cells(pool, types, args.cells_per_type, args.seed, out)
    assert not selected.individual_id.isin(d.loc[d.split.eq("holdout"), "Patient"]).any()
    with h5py.File(args.input, "r") as h:
        ids = read_column(h["obs"][h["obs"].attrs.get("_index", "_index")]).astype(str)
        assert ids.is_unique
        rows = pd.Index(ids).get_indexer(selected.obs_name.astype(str)); assert (rows >= 0).all()
        for field in ["sample_id", "Patient", "Disease", "Site", "final_analysis"]:
            now = read_column(h["obs"][field]).iloc[rows].astype(str).to_numpy()
            assert np.array_equal(now, selected[field].astype(str).to_numpy()), field + " mismatch"
        order = np.argsort(rows); rows = rows[order]; selected = selected.iloc[order].reset_index(drop=True)
        selected["source_cell_row"] = rows
        genes = pd.DataFrame({"source_gene_id": read_column(h["var"]["gene_id"]).astype("string"),
                              "gene_symbol": read_column(h["var"]["gene_symbol"]).astype("string")})
        assert genes.source_gene_id.is_unique
        print(f"{args.mode}: extracting {len(rows)} selected cells", flush=True)
        csr = selected_csr(h["X"], rows, len(genes))
    decisions = {"training_donors": int(d.split.eq("training").sum()), "holdout_donors": int(d.split.eq("holdout").sum()),
        "same_gse_donor_split_as_colon_pilot": True, "source_sites": sorted(selected.Site.unique().tolist()),
        "sampling_policy": "fixed cells/type, donor-round-robin cap5, no forced disease/site/subtype quotas",
        "taxonomy_policy": "original20; author-confirmed Glial/Tuft promoted for this pilot only; BEST4 in corresponding absorptive type",
        "candidate_sha256": sha256(Path(args.candidate)), "donor_split_sha256": sha256(Path(args.donor_split))}
    return write_outputs(out, selected, genes, csr, types, args.mode, decisions, started)


def stream_mtx_selected(path, genes, barcodes, selected, gene_lookup, chunksize=1_000_000):
    """Scan text once, retain only selected columns; never materialize full matrix."""
    with Path(path).open() as f:
        header = f.readline().strip()
        if header != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError("Require released raw integer MatrixMarket counts")
        skip = 1
        for line in f:
            skip += 1
            if line.startswith("%"):
                continue
            nr, nc, nnz = map(int, line.split()); break
        else:
            raise ValueError("Missing MTX dimensions")
    assert nr == len(genes) and nc == len(barcodes)
    index = pd.Index(barcodes); assert index.is_unique
    source_col = index.get_indexer(selected.cell_id.astype(str)); assert (source_col >= 0).all()
    lut = np.full(nc, -1, dtype=np.int32); lut[source_col] = selected.index.to_numpy(dtype=np.int32)
    mapped_genes = np.asarray([gene_lookup[g] for g in genes], dtype=np.int32)
    all_cell, all_gene, all_values = [], [], []; scanned = 0
    reader = pd.read_csv(path, sep=r"\s+", header=None, names=["gene", "cell", "value"],
                         skiprows=skip, dtype=np.int64, chunksize=chunksize, engine="c")
    for chunk in reader:
        arr = chunk.to_numpy(); assert arr.shape[1] == 3
        if ((arr[:, 0] < 1) | (arr[:, 0] > nr) | (arr[:, 1] < 1) | (arr[:, 1] > nc)).any():
            raise ValueError("MTX index outside declared dimensions")
        if (arr[:, 2] < 0).any() or (arr[:, 2] > np.iinfo(np.int32).max).any():
            raise ValueError("MTX counts must be nonnegative bounded integers")
        ci = lut[arr[:, 1]-1]; use = (ci >= 0) & (arr[:, 2] > 0)
        all_cell.append(ci[use]); all_gene.append(mapped_genes[arr[use, 0]-1]); all_values.append(arr[use, 2].astype(np.int32))
        scanned += len(chunk)
        if scanned % 20_000_000 == 0:
            print(f"{Path(path).name}: scanned {scanned:,}/{nnz:,} entries", flush=True)
    assert scanned == nnz, "MTX truncated/extra entries"
    return np.concatenate(all_cell), np.concatenate(all_gene), np.concatenate(all_values), nnz


def coalesce_selected(cell, gene, values, n_cells, n_genes):
    keys = cell.astype(np.int64) * n_genes + gene
    order = np.argsort(keys, kind="stable"); keys, values = keys[order], values[order]
    unique, first = np.unique(keys, return_index=True)
    values = np.add.reduceat(values.astype(np.int64), first)
    assert values.max() <= np.iinfo(np.int32).max
    ci, gi = unique // n_genes, unique % n_genes
    detected = np.bincount(ci, minlength=n_cells)
    totals = np.bincount(ci, weights=values, minlength=n_cells).astype(np.int64)
    assert (totals > 0).all()
    pointers = np.r_[0, np.cumsum(detected)]
    return gi.astype(np.int32), values.astype(np.int32), pointers, totals, detected


def run_scp(args):
    started = time.monotonic(); out = fresh(args.outdir)
    dump_json(out / "run_parameters.json", {"version": VERSION, **vars(args)})
    c = pd.read_parquet(args.candidate)
    assert c.cell_id.is_unique and c.NAME.eq(c.cell_id).all()
    c["individual_id"] = c.donor_id.astype(str); c["Disease"] = c.disease; c["Inflammation"] = c.inflammation
    c["Site"] = "ColonRectum_source_not_cell_resolved"  # Epi/LP is a fraction, not an anatomic site.
    c["pilot_type"] = c.coarse_identity; c["author_label"] = c.Cluster
    c["n_genes_by_counts"] = pd.to_numeric(c.nGene, errors="raise")
    c["total_counts"] = pd.to_numeric(c.nUMI, errors="raise")
    donors = split_scp_donors(c, args.seed)
    c["split"] = c.individual_id.map(donors.set_index("individual_id").split)
    donors.to_csv(out / "donor_split.tsv", sep="\t", index=False)
    promoted = c.pilot_type.isin(["Glial", "Tuft"]) & c.primary_reference_action.eq("HOLD_CROSS_DATASET")
    action = c.primary_reference_action.eq("KEEP") | promoted
    c.groupby(["pilot_type", "primary_reference_action"]).size().rename("n_cells").reset_index().to_csv(
        out / "taxonomy_actions.tsv", sep="\t", index=False)
    valid = action & c.pilot_type.isin(SCP_TYPES) & c.split.eq("training")
    for field in ["total_counts", "n_genes_by_counts"]: valid &= np.isfinite(c[field]) & c[field].gt(0)
    selected = choose_cells(c[valid].sort_values("cell_id"), SCP_TYPES, args.cells_per_type, args.seed, out)
    assert not selected.individual_id.isin(donors.loc[donors.split.eq("holdout"), "individual_id"]).any()
    # Compare cached annotations with released text metadata before expression scan.
    raw = Path(args.raw)
    metadata = pd.read_csv(raw / "metadata/all.meta2.txt", sep="\t", dtype=str)
    metadata = metadata[metadata.NAME.ne("TYPE")]
    assert metadata.NAME.is_unique
    now = metadata.set_index("NAME").loc[selected.cell_id]
    for field in ["Cluster", "Subject", "Health", "Location", "Sample"]:
        assert np.array_equal(now[field].to_numpy(), selected[field].astype(str).to_numpy()), field + " source mismatch"
    bundles = {}; all_genes = set()
    for name, (folder, prefix) in BUNDLES.items():
        directory = raw / "expression" / folder
        genes = pd.read_csv(directory / (prefix+".genes.tsv"), sep="\t", header=None, dtype=str)
        bc = pd.read_csv(directory / (prefix+".barcodes2.tsv"), sep="\t", header=None, dtype=str)
        assert genes.shape[1] == bc.shape[1] == 1 and genes[0].is_unique and bc[0].is_unique
        genes, bc = genes[0].tolist(), bc[0].tolist()
        bundles[name] = (directory / ("gene_sorted-"+prefix+".matrix.mtx"), genes, bc)
        all_genes.update(genes)
    universe = sorted(all_genes); lookup = {g: i for i, g in enumerate(universe)}
    cells, genes, values, scans = [], [], [], []
    for name, (path, gene_names, barcodes) in bundles.items():
        q = selected[selected.expression_bundle.eq(name)]
        assert len(q) > 0
        print(f"SCP259 {name}: scan released matrix, retain {len(q)} selected cells", flush=True)
        ci, gi, vv, nnz = stream_mtx_selected(path, gene_names, barcodes, q, lookup)
        cells.append(ci); genes.append(gi); values.append(vv)
        scans.append({"bundle": name, "source_nnz_scanned": nnz, "selected_cells": len(q), "retained_entries": len(vv),
                      "gene_count": len(gene_names), "source_bytes": path.stat().st_size})
    pd.DataFrame(scans).to_csv(out / "expression_bundle_scan_audit.tsv", sep="\t", index=False)
    csr = coalesce_selected(np.concatenate(cells), np.concatenate(genes), np.concatenate(values), len(selected), len(universe))
    genes = pd.DataFrame({"gene_symbol": universe, "source_gene_id": universe})
    decisions = {"training_donors": int(donors.split.eq("training").sum()), "holdout_donors": int(donors.split.eq("holdout").sum()),
        "sampling_policy": "fixed cells/type, donor-round-robin cap5; no forced disease/subtype quotas",
        "taxonomy_policy": "original20 supported subset18; author Glial/Tuft included; BEST4 in Colonic_absorptive; no Paneth label",
        "anatomical_site_annotation": "Public cell metadata Epi/LP describes tissue fraction; colon-vs-rectum unresolved here",
        "source_gene_union": len(universe), "matrix_read_method": "single streaming scan per bundle, selected columns only",
        "candidate_sha256": sha256(Path(args.candidate))}
    return write_outputs(out, selected, genes, csr, SCP_TYPES, args.mode, decisions, started)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["gse_joint", "gse_ileum", "scp259"], required=True)
    p.add_argument("--input")
    p.add_argument("--raw")
    p.add_argument("--candidate", required=True)
    p.add_argument("--donor-split")
    p.add_argument("--outdir", required=True)
    p.add_argument("--cells-per-type", type=int, default=50)
    p.add_argument("--seed", type=int, default=20261006)
    args = p.parse_args()
    if args.cells_per_type < 10: p.error("--cells-per-type >=10 required")
    if args.mode.startswith("gse"):
        if not args.input or not args.donor_split: p.error("GSE needs --input and --donor-split")
        run_gse(args)
    else:
        if not args.raw: p.error("SCP259 needs --raw")
        run_scp(args)


if __name__ == "__main__":
    main()
