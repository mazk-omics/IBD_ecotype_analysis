#!/usr/bin/env python3
"""Read-only, bounded IBDverse pilot extraction. Requires h5py/numpy/pandas only.

Uses the completed v4 aggregate metadata audit, not another sidecar merge.
Dependencies: inspect_ibdverse_h5ad.py and plan_ibdverse_epithelial_pilot.py.
The output is a labeled scRNA reference INPUT, not a finished signature matrix.
"""
import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
import tarfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from inspect_ibdverse_h5ad import frame_columns, read_column
from plan_ibdverse_epithelial_pilot import (
    WHOLE_TYPE_ORDER, PENDING_T_LABELS, map_whole_type, stable_rng,
)

VERSION = "20261006_pilot_export_v2_qc_compat"
QC_FIELDS = ["conf_score", "n_genes_by_counts", "total_counts",
             "pct_counts_gene_group__mito_transcript"]
MARKERS = ["EPCAM", "KRT8", "ALPI", "SI", "FABP1", "BEST4", "OTOP2",
           "LGR5", "OLFM4", "MKI67", "MUC2", "CLCA1", "CHGA", "CHGB",
           "TPH1", "POU2F3", "DCLK1", "DEFA5", "DEFA6", "ITLN2", "PLA2G2A",
           "CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "TRDC", "TRGC1",
           "TRGC2", "NKG7", "GNLY", "KLRD1", "MS4A1", "CD79A", "JCHAIN",
           "MZB1", "SDC1", "LYZ", "LST1", "S100A8", "FCER1A", "CLEC9A",
           "KIT", "TPSAB1", "COL1A1", "DCN", "PDGFRB", "RGS5", "PECAM1",
           "VWF", "SOX10", "PLP1", "S100B", "GFAP", "SOD2"]


def dump_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_donors(donors, fraction, seed):
    donors = donors.copy().sort_values("individual_id")
    if donors.individual_id.duplicated().any():
        raise ValueError("A donor has conflicting strata; cannot split or balance safely")
    donors["split"] = "training"
    for stratum in sorted(donors.stratum.unique()):
        positions = donors.index[donors.stratum.eq(stratum)].to_numpy()
        rng = stable_rng("whole_donor_holdout", stratum, seed)
        n_holdout = int(round(len(positions) * fraction))
        chosen = rng.choice(positions, n_holdout, replace=False)
        donors.loc[chosen, "split"] = "holdout"
    return donors


def write_qc_summary(frame, groups, metrics, path):
    """Flat QC table; avoid pandas 2.3 groupby.describe NA-key concatenation."""
    records = []
    for key, subset in frame.groupby(groups, observed=True, dropna=False, sort=True):
        keys = key if isinstance(key, tuple) else (key,)
        identity = {name: ("<MISSING>" if pd.isna(value) else str(value))
                    for name, value in zip(groups, keys)}
        for metric in metrics:
            values = pd.to_numeric(subset[metric], errors="coerce").to_numpy(dtype=float, na_value=np.nan)
            finite = values[np.isfinite(values)]
            quantiles = np.quantile(finite, [0, .25, .5, .75, 1]) if len(finite) else [np.nan] * 5
            records.append({**identity, "metric": metric, "n_cells": len(values),
                "n_finite": len(finite), "n_missing": int(np.isnan(values).sum()),
                "n_nonfinite": int((~np.isfinite(values)).sum()),
                "mean": float(finite.mean()) if len(finite) else np.nan,
                "std": float(finite.std(ddof=1)) if len(finite) > 1 else np.nan,
                **dict(zip(["min", "q25", "median", "q75", "max"], quantiles))})
    columns = groups + ["metric", "n_cells", "n_finite", "n_missing", "n_nonfinite",
                        "mean", "std", "min", "q25", "median", "q75", "max"]
    pd.DataFrame(records, columns=columns).to_csv(path, sep="\t", index=False)


def balanced_pick(frame, target, group_name, stratum, seed, max_per_donor):
    counts = frame.groupby("individual_id", observed=True).size().sort_index()
    capacities = counts.clip(upper=max_per_donor)
    if int(capacities.sum()) < target:
        return [], {"group": group_name, "stratum": stratum, "target": target,
                    "available_cells": len(frame), "available_donors": len(counts),
                    "capacity_with_donor_cap": int(capacities.sum()), "selected_cells": 0,
                    "status": "shortage"}
    order = stable_rng(group_name, stratum, seed).permutation(len(counts))
    quotas = np.zeros(len(counts), dtype=int)
    remaining = target
    for _ in range(max_per_donor):
        for position in order:
            if quotas[position] < capacities.iloc[position]:
                quotas[position] += 1
                remaining -= 1
                if remaining == 0:
                    break
        if remaining == 0:
            break
    selected = []
    for position in np.flatnonzero(quotas):
        donor = counts.index[position]
        candidates = frame.loc[frame.individual_id.eq(donor), "source_cell_row"].to_numpy(copy=True)
        candidates.sort()
        rng = stable_rng(group_name + "|" + donor, stratum, seed)
        selected.extend(rng.choice(candidates, int(quotas[position]), replace=False).tolist())
    return selected, {"group": group_name, "stratum": stratum, "target": target,
                      "available_cells": len(frame), "available_donors": len(counts),
                      "capacity_with_donor_cap": int(capacities.sum()),
                      "selected_cells": len(selected), "selected_donors": int((quotas > 0).sum()),
                      "max_selected_cells_per_donor": int(quotas.max()), "status": "ok"}


def selected_csr(node, rows, n_genes):
    encoding = node.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding != "csr_matrix" or tuple(node.attrs["shape"])[1] != n_genes:
        raise ValueError("layers/counts must be a cell-by-gene CSR matrix")
    chunks_i, chunks_v, pointers, totals, detected = [], [], [0], [], []
    for position, row in enumerate(rows):
        start, stop = map(int, node["indptr"][int(row):int(row) + 2])
        if not 0 <= start <= stop <= len(node["data"]) or stop - start > 10_000_000:
            raise ValueError(f"Invalid/unbounded CSR pointers at cell row {row}")
        indices = np.asarray(node["indices"][start:stop], dtype=np.int64)
        values = np.asarray(node["data"][start:stop], dtype=np.float64)
        if (len(indices) and (indices.min() < 0 or indices.max() >= n_genes)):
            raise ValueError(f"Invalid gene index at cell row {row}")
        if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.rint(values), rtol=0, atol=1e-5):
            raise ValueError(f"counts must be finite nonnegative integers; failed cell row {row}")
        if len(values) and values.max() > np.iinfo(np.int32).max:
            raise ValueError(f"Selected count exceeds int32 capacity at row {row}")
        order = np.argsort(indices, kind="stable")
        indices, values = indices[order], np.rint(values[order]).astype(np.int64)
        if len(indices):
            unique, first = np.unique(indices, return_index=True)
            values = np.add.reduceat(values, first)
            indices = unique
        positive = values > 0
        indices, values = indices[positive], values[positive]
        if values.sum() <= 0:
            raise ValueError(f"Selected cell row {row} has zero library size")
        if values.max() > np.iinfo(np.int32).max:
            raise ValueError(f"Selected count exceeds int32 capacity at row {row}")
        chunks_i.append(indices.astype(np.int32))
        chunks_v.append(values.astype(np.int32))
        pointers.append(pointers[-1] + len(values))
        totals.append(int(values.sum()))
        detected.append(len(values))
        if (position + 1) % 500 == 0:
            print(f"Extracted {position + 1}/{len(rows)} sparse cells", flush=True)
    return (np.concatenate(chunks_i), np.concatenate(chunks_v),
            np.asarray(pointers, dtype=np.int64), np.asarray(totals), np.asarray(detected))


def run(args):
    started = time.monotonic()
    source = Path(args.input).expanduser().resolve(strict=True)
    coverage_path = Path(args.coverage).expanduser().resolve(strict=True)
    out = Path(args.outdir).expanduser().resolve()
    if out == source.parent or out == source or (out.exists() and any(out.iterdir())):
        raise ValueError("Use a new, empty output directory separate from the source directory")
    if args.cells_per_type < 8 or args.cells_per_type % 4:
        raise ValueError("--cells-per-type must be >=8 and divisible by 4")
    if not 0 <= args.holdout_fraction < 0.5 or not 1 <= args.max_per_donor <= 20:
        raise ValueError("Require 0 <= holdout fraction <0.5 and 1 <= max-per-donor <=20")
    if args.min_conf_score is not None and not np.isfinite(args.min_conf_score):
        raise ValueError("--min-conf-score must be finite")
    out.mkdir(parents=True, exist_ok=True)
    dump_json(out / "run_parameters.json", {"version": VERSION, **vars(args),
              "python": sys.executable, "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
              "h5py": h5py.__version__, "numpy": np.__version__, "pandas": pd.__version__,
              "source_bytes": source.stat().st_size, "coverage_sha256": sha256(coverage_path)})
    coverage = pd.read_csv(coverage_path, sep="\t", keep_default_na=False)
    cols = ["sanger_sample_id", "individual_id", "disease_status", "ti_ses_inflamed",
            "biopsy_type", "cohort", "atlas", "reference_metadata_complete"]
    metadata = coverage[cols].drop_duplicates()
    if metadata.sanger_sample_id.duplicated().any():
        raise ValueError("Aggregate audit contains inconsistent sample metadata; inspect rather than keep first")
    complete = metadata.reference_metadata_complete.astype(str).str.lower().eq("true")
    known = metadata.loc[complete & metadata.biopsy_type.eq("terminal ileum")].copy()
    known["stratum"] = known.disease_status + "_" + known.ti_ses_inflamed
    strata_targets = {"Healthy_uninflamed": args.cells_per_type // 2,
                      "CD_inflamed": args.cells_per_type // 4,
                      "CD_uninflamed": args.cells_per_type // 4}
    if set(known.stratum) != set(strata_targets) or known.individual_id.eq("").any():
        raise ValueError("Unexpected disease/inflammation strata or empty donor IDs")
    donors = known[["individual_id", "stratum"]].drop_duplicates().sort_values("individual_id")
    donors = split_donors(donors, args.holdout_fraction, args.seed)
    donors.to_csv(out / "donor_split.tsv", sep="\t", index=False)
    known = known.merge(donors, on=["individual_id", "stratum"], validate="many_to_one")

    with h5py.File(source, "r") as h:
        oi, obs_columns = frame_columns(h["obs"])
        vi, _ = frame_columns(h["var"])
        missing = set(["sanger_sample_id", "label", "category"] + QC_FIELDS) - set(obs_columns)
        if missing:
            raise ValueError(f"Required obs fields missing: {sorted(missing)}")
        n_cells, n_genes = len(h["obs"][oi]), len(h["var"][vi])
        if int(pd.to_numeric(coverage.n_cells).sum()) != n_cells:
            raise ValueError("Audit cell count differs from source; use audit from the same H5AD")
        print(f"Reading annotation/QC columns for {n_cells:,} cells; expression remains unopened", flush=True)
        obs = pd.DataFrame({name: read_column(h["obs"][name]) for name in
                            ["sanger_sample_id", "label", "category"] + QC_FIELDS})
        # Verify the whole annotation census before using the downloaded joins.
        census = obs.groupby(["sanger_sample_id", "label", "category"], observed=True).size().rename("actual_n").reset_index()
        check = coverage[["sanger_sample_id", "label", "category", "n_cells"]].merge(
            census, on=["sanger_sample_id", "label", "category"], how="outer", validate="one_to_one")
        if check[["n_cells", "actual_n"]].isna().any().any() or not check.n_cells.astype(int).eq(check.actual_n).all():
            raise ValueError("Source annotations differ from audited census")
        obs["source_cell_row"] = np.arange(n_cells, dtype=np.int64)
        obs["pilot_type"] = obs.label.map(map_whole_type).astype("string")
        unexpected = set(obs.loc[obs.pilot_type.isna(), "label"].astype(str)) - PENDING_T_LABELS
        if unexpected:
            raise ValueError(f"Unexpected unmapped labels: {sorted(unexpected)}")
        pool = obs.merge(known, on="sanger_sample_id", how="inner", validate="many_to_one")
        for name in QC_FIELDS:
            pool[name] = pd.to_numeric(pool[name], errors="coerce").astype("float64")
        valid = np.isfinite(pool[QC_FIELDS]).all(axis=1) & pool.n_genes_by_counts.gt(0) & pool.total_counts.gt(0)
        valid &= pool.pct_counts_gene_group__mito_transcript.between(0, 100)
        if args.min_conf_score is not None:
            valid &= pool.conf_score.ge(args.min_conf_score)
        pool.loc[~valid, ["source_cell_row", "sanger_sample_id", "label", "individual_id"] + QC_FIELDS].to_csv(
            out / "excluded_by_validity_or_user_confidence.tsv", sep="\t", index=False)
        eligible = pool.loc[valid & pool.split.eq("training")].copy()
        print(f"Training donors: {(donors.split == 'training').sum()}; holdout donors: {(donors.split == 'holdout').sum()}; eligible training cells: {len(eligible):,}", flush=True)
        grouped_counts = pool.groupby(["pilot_type", "label", "stratum", "split"], dropna=False, observed=True).size()
        grouped_counts.rename("n_cells").reset_index().to_csv(out / "candidate_coverage.tsv", sep="\t", index=False)
        qc_pool = pool.assign(eligibility=np.where(valid, "eligible_before_donor_split", "excluded"))
        # Pending labels are audit groups, not signature types. Name them only in
        # this QC table, keeping pilot_type missing for all sampling/export logic.
        qc_pool["qc_identity"] = qc_pool.pilot_type.fillna("PENDING::" + qc_pool.label.astype("string"))
        write_qc_summary(qc_pool, ["eligibility", "qc_identity"], QC_FIELDS,
                         out / "candidate_qc_summary.tsv")
        picked, statuses, review_rows = [], [], []
        for type_name in WHOLE_TYPE_ORDER:
            for stratum, target in strata_targets.items():
                frame = eligible.loc[eligible.pilot_type.eq(type_name) & eligible.stratum.eq(stratum)]
                rows, status = balanced_pick(frame, target, type_name, stratum, args.seed, args.max_per_donor)
                picked.extend(rows)
                statuses.append(status)
        # Small, separate identity-review panel for the two unassigned T labels.
        for label in sorted(PENDING_T_LABELS):
            frame = eligible.loc[eligible.label.eq(label)]
            capacity = int(frame.groupby("individual_id", observed=True).size().clip(upper=args.max_per_donor).sum())
            target = min(args.review_cells_per_pending_label, capacity)
            if target:
                rows, _ = balanced_pick(frame, target, label, "review_all_strata", args.seed, args.max_per_donor)
                review_rows.extend(rows)
        pd.DataFrame(statuses).to_csv(out / "sampling_capacity.tsv", sep="\t", index=False)
        shortages = [r for r in statuses if r["status"] != "ok"]
        if shortages:
            raise ValueError(f"{len(shortages)} strata cannot meet the quota after donor split/QC; see sampling_capacity.tsv")
        all_rows = np.sort(np.asarray(picked + review_rows, dtype=np.int64))
        if len(np.unique(all_rows)) != len(all_rows):
            raise ValueError("A cell was selected twice")
        selected = pool.set_index("source_cell_row", drop=False).loc[all_rows].reset_index(drop=True)
        selected["usage"] = np.where(selected.source_cell_row.isin(picked), "web_reference", "identity_review_only")
        selected["analysis_group"] = selected.pilot_type.fillna(selected.label.astype("string"))
        index_node = h["obs"][oi]
        if isinstance(index_node, h5py.Dataset):
            selected["source_cell_id"] = index_node.asstr()[all_rows]
        else:
            selected["source_cell_id"] = read_column(index_node).iloc[all_rows].to_numpy()
        genes = pd.DataFrame({"source_gene_id": read_column(h["var"][vi]).astype("string"),
                              "gene_symbol": read_column(h["var"]["gene_symbols"]).astype("string")})
        genes["source_gene_row"] = np.arange(n_genes)
        genes["ambiguous_symbol"] = genes.gene_symbol.duplicated(keep=False)
        symbol_valid = genes.gene_symbol.notna() & genes.gene_symbol.str.strip().ne("")
        symbol_valid &= ~genes.gene_symbol.str.contains(r"[\t\r\n]", na=True)
        gene_eligible = (symbol_valid & ~genes.ambiguous_symbol).to_numpy(dtype=bool)
        print(f"Reading {len(picked)} web-reference cells + {len(review_rows)} identity-review cells in source row order", flush=True)
        indices, values, indptr, totals, detected = selected_csr(h["layers/counts"], all_rows, n_genes)
    selected["extracted_total_counts"] = totals
    selected["extracted_n_genes"] = detected
    selected["total_counts_matches_obs"] = np.isclose(totals, selected.total_counts, rtol=1e-5, atol=1)
    selected["n_genes_matches_obs"] = detected == selected.n_genes_by_counts
    selected["web_column"] = ""
    web_mask = selected.usage.eq("web_reference").to_numpy()
    for type_name in WHOLE_TYPE_ORDER:
        positions = selected.index[web_mask & selected.pilot_type.eq(type_name).fillna(False).to_numpy(dtype=bool)]
        selected.loc[positions, "web_column"] = [f"{type_name}.{i + 1}" for i in range(len(positions))]
    selected.to_csv(out / "selected_cell_manifest.tsv", sep="\t", index=False)
    selected.groupby(["usage", "analysis_group", "label", "stratum"], observed=True).size().rename("n_cells").reset_index().to_csv(
        out / "selected_original_label_composition.tsv", sep="\t", index=False)
    selected.groupby(["usage", "analysis_group", "individual_id", "stratum"], observed=True).size().rename("n_cells").reset_index().to_csv(
        out / "selected_donor_contributions.tsv", sep="\t", index=False)
    write_qc_summary(selected, ["usage", "analysis_group"],
                     QC_FIELDS + ["extracted_total_counts", "extracted_n_genes"],
                     out / "selected_qc_summary.tsv")
    nnz_cell = np.repeat(np.arange(len(selected), dtype=np.int32), np.diff(indptr))
    eligible_entry = web_mask[nnz_cell] & gene_eligible[indices]
    gene_sums = np.bincount(indices[eligible_entry], weights=values[eligible_entry], minlength=n_genes)
    gene_keep = gene_eligible & (gene_sums > 0)
    genes["web_included"] = gene_keep
    genes["web_exclusion_reason"] = np.select(
        [~symbol_valid.to_numpy(dtype=bool), genes.ambiguous_symbol.to_numpy(), ~gene_keep],
        ["invalid_symbol", "ambiguous_symbol_all_loci_excluded_in_pilot", "zero_in_web_cells"], default="")
    genes.to_csv(out / "gene_export_manifest.tsv", sep="\t", index=False)
    # Retain ALL original gene loci in a small custom HDF5, including ambiguities.
    with h5py.File(out / "selected_counts.h5", "w") as target:
        target.attrs["format"] = "IBDverse pilot custom CSR; not an AnnData h5ad"
        csr = target.create_group("counts")
        csr.attrs["encoding-type"] = "csr_matrix"
        csr.attrs["shape"] = [len(selected), n_genes]
        for name, data in [("data", values), ("indices", indices), ("indptr", indptr)]:
            csr.create_dataset(name, data=data, compression="gzip", compression_opts=1)
        strings = h5py.string_dtype("utf-8")
        for name, data in [("gene_symbol", genes.gene_symbol), ("source_gene_id", genes.source_gene_id),
                           ("cell_id", selected.source_cell_id), ("usage", selected.usage),
                           ("analysis_group", selected.analysis_group)]:
            target.create_dataset(name, data=np.asarray(data.fillna("").astype(str), dtype=object), dtype=strings)
        target.create_dataset("source_cell_row", data=all_rows)
    marker_rows = []
    marker_groups = selected.groupby(["analysis_group", "stratum"], observed=True).groups
    for symbol in MARKERS:
        gene_positions = genes.index[genes.gene_symbol.eq(symbol).fillna(False)].to_numpy()
        for gene in gene_positions:
            mask = indices == gene
            cell_counts = np.bincount(nnz_cell[mask], weights=values[mask], minlength=len(selected))
            for (group, stratum), positions in marker_groups.items():
                subset = cell_counts[np.asarray(positions)]
                cpm = subset / totals[np.asarray(positions)] * 1_000_000
                marker_rows.append({"analysis_group": str(group), "stratum": stratum, "gene_symbol": symbol,
                    "source_gene_id": str(genes.at[gene, "source_gene_id"]),
                    "ambiguous_symbol": bool(genes.at[gene, "ambiguous_symbol"]),
                    "gene_exported_to_web": bool(gene_keep[gene]), "n_cells": len(positions),
                    "detection_fraction": float((subset > 0).mean()), "mean_cpm_original_library": float(cpm.mean()),
                    "note": "diagnostic panel only; no automatic identity calls; duplicate loci reported separately"})
    pd.DataFrame(marker_rows).to_csv(out / "marker_expression_summary.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_symbol": MARKERS,
                  "present_in_source_gene_universe": [bool(genes.gene_symbol.eq(m).any()) for m in MARKERS]}).to_csv(
        out / "marker_gene_presence.tsv", sep="\t", index=False)
    # Build CSC ordering only for the selected, web-eligible entries; bounded blocks.
    use = web_mask[nnz_cell] & gene_keep[indices]
    cell_lookup = np.full(len(selected), -1, dtype=np.int32)
    cell_lookup[web_mask] = np.arange(web_mask.sum(), dtype=np.int32)
    gene_lookup = np.full(n_genes, -1, dtype=np.int32)
    gene_lookup[gene_keep] = np.arange(gene_keep.sum(), dtype=np.int32)
    gi, ci, vv = gene_lookup[indices[use]], cell_lookup[nnz_cell[use]], values[use]
    web_totals = np.bincount(ci, weights=vv, minlength=int(web_mask.sum()))
    if not (web_totals > 0).all():
        raise ValueError("A web cell has no counts after symbol filtering")
    order = np.argsort(gi, kind="stable")
    gi, ci, vv = gi[order], ci[order], vv[order]
    pointers = np.concatenate([[0], np.cumsum(np.bincount(gi, minlength=int(gene_keep.sum())))])
    symbols = genes.loc[gene_keep, "gene_symbol"].astype(str).to_numpy()
    web_name = "IBDverse_TI_18types_scRNA_reference.txt"
    print(f"Writing {len(symbols):,} genes x {web_mask.sum():,} cells in bounded blocks", flush=True)
    with (out / web_name).open("w", encoding="utf-8", newline="") as stream:
        stream.write("GeneSymbol\t" + "\t".join(selected.loc[web_mask, "web_column"]) + "\n")
        for begin in range(0, len(symbols), 128):
            end = min(begin + 128, len(symbols))
            block = np.zeros((end - begin, int(web_mask.sum())), dtype=np.int32)
            left, right = int(pointers[begin]), int(pointers[end])
            np.add.at(block, (gi[left:right] - begin, ci[left:right]), vv[left:right])
            for symbol, row in zip(symbols[begin:end], block):
                stream.write(symbol + "\t" + "\t".join(map(str, row)) + "\n")
    with (out / web_name).open("rb") as raw, gzip.open(out / (web_name + ".gz"), "wb", compresslevel=1) as compressed:
        shutil.copyfileobj(raw, compressed, length=1024 * 1024)
    result = {"version": VERSION, "status": "completed", "source_opened_read_only": True,
        "n_source_cells": n_cells, "n_source_genes": n_genes,
        "web_reference_cells": int(web_mask.sum()), "web_reference_types": WHOLE_TYPE_ORDER,
        "web_reference_genes": int(gene_keep.sum()), "identity_review_cells": int((~web_mask).sum()),
        "training_donors": int(donors.split.eq("training").sum()), "holdout_donors": int(donors.split.eq("holdout").sum()),
        "actual_web_donors": selected.loc[web_mask, "individual_id"].nunique(),
        "holdout_donors_used_for_extraction": int(selected.individual_id.isin(donors.loc[donors.split.eq("holdout"), "individual_id"]).sum()),
        "max_cells_per_donor_per_type": int(selected.loc[web_mask].groupby(["pilot_type", "individual_id"], observed=True).size().max()),
        "web_unit": "raw integer UMI counts; CIBERSORTx scRNA module handles CPM normalization",
        "max_web_count": int(vv.max()), "web_format_valid": bool(vv.max() >= 50),
        "qc_policy": "author labels plus finite positive basic metrics; no new global QC or confidence threshold unless supplied",
        "biological_identity_review_complete": False, "is_finished_signature_matrix": False,
        "total_counts_matching_obs_cells": int(selected.total_counts_matches_obs.sum()),
        "n_genes_matching_obs_cells": int(selected.n_genes_matches_obs.sum()),
        "excluded_invalid_or_user_confidence_cells": int((~valid).sum()),
        "seconds_elapsed": round(time.monotonic() - started, 2), "web_input_bytes": (out / web_name).stat().st_size,
        "web_input_sha256": sha256(out / web_name),
        "notes": ["Review selected QC, label composition and marker expression before interpreting results",
                  "No Paneth/Glial/NK class inferred; pending T cells excluded from web input",
                  "Ambiguous gene symbols excluded from web pilot but all original loci retained in selected_counts.h5",
                  "If max_web_count <50, CIBERSORTx may assume log input; do not upload without revising scale"]}
    dump_json(out / "summary.json", result)
    with tarfile.open(out / "small_audit_reports.tar.gz", "w:gz") as archive:
        for path in sorted(out.iterdir()):
            if path.suffix in {".tsv", ".json"}:
                archive.add(path, arcname=path.name, recursive=False)
    print(json.dumps(result, indent=2, ensure_ascii=False), flush=True)
    print(f"Finished: {out}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    parser.add_argument("--input", required=True)
    parser.add_argument("--coverage", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cells-per-type", type=int, default=200)
    parser.add_argument("--holdout-fraction", type=float, default=0.2)
    parser.add_argument("--max-per-donor", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261006)
    parser.add_argument("--min-conf-score", type=float)
    parser.add_argument("--review-cells-per-pending-label", type=int, default=60)
    args = parser.parse_args()
    if not 0 <= args.review_cells_per_pending_label <= 200:
        parser.error("--review-cells-per-pending-label must be between 0 and 200")
    try:
        run(args)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        # Refuse invalid data, with an actionable message instead of silently changing quotas.
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        print(f"Inspect output diagnostics if created: {args.outdir}", file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
