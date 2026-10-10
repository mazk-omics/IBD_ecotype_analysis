#!/usr/bin/env python3
"""IBDverse stage 2: metadata joins and a bounded sparse-row expression check.

Run in ibd_refaudit. Keep inspect_ibdverse_h5ad.py in the same script directory.
Reads annotations, gene names, small sidecar tables and at most --row-samples
rows per expression layer. Does not select cells, change labels, or write a matrix.
Conflicting records in all three metadata tables are quarantined: shared attributes are
retained, conflicting fields remain missing. Cell totals are never duplicated.
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from inspect_ibdverse_h5ad import read_column, frame_columns

SCRIPT_VERSION = "20261006_metadata_keys_v4_all_tables"


def read_metadata(path, key, required, report_dir=None, conflict_policy="error"):
    table = pd.read_csv(path, sep="\t", dtype="string", keep_default_na=False)
    missing = set(required) - set(table.columns)
    if missing:
        raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
    stats = {"file": path.name, "key": key, "input_rows": len(table)}
    # Normalize whitespace in join identifiers only. Leave phenotype labels intact.
    original_keys = table[key].copy()
    table[key] = table[key].str.strip()
    stats["keys_whitespace_normalized"] = int(original_keys.ne(table[key]).sum())
    missing_key = table[key].str.lower().isin(["", "na", "n/a", "nan", "none", "null", "<na>"])
    stats["missing_key_rows"] = int(missing_key.sum())
    duplicate_key = table.loc[~missing_key, key].duplicated(keep=False)
    duplicates = table.loc[~missing_key].loc[duplicate_key]
    if report_dir is not None:
        table.loc[missing_key].to_csv(report_dir / f"{path.stem}.missing_keys.tsv", sep="\t", index=False)
        duplicates.to_csv(report_dir / f"{path.stem}.duplicate_keys.tsv", sep="\t", index=False)
    # A missing join identifier cannot match any valid observed sample/donor.
    # Exclude those metadata rows, recording them, then verify all cell mappings.
    table = table.loc[~missing_key].copy()
    before = len(table)
    table = table.drop_duplicates()
    stats["exact_duplicate_rows_removed"] = before - len(table)
    conflicts = table.loc[table[key].duplicated(keep=False)].copy()
    stats["conflicting_key_groups"] = int(conflicts[key].nunique())
    stats["output_rows"] = len(table)
    if report_dir is not None:
        conflicts.to_csv(report_dir / f"{path.stem}.conflicting_rows.tsv", sep="\t", index=False)
        (report_dir / f"{path.stem}.key_audit.json").write_text(
            json.dumps(stats, indent=2, ensure_ascii=False), encoding="utf-8")
    if not conflicts.empty:
        differing = []
        for identifier, group in conflicts.groupby(key, dropna=False, sort=False):
            columns = [c for c in table.columns if c != key and group[c].nunique(dropna=False) > 1]
            differing.append({"key_value": str(identifier), "differing_columns": columns,
                              "rows": len(group)})
        if report_dir is not None:
            (report_dir / f"{path.stem}.conflict_details.json").write_text(
                json.dumps(differing, indent=2, ensure_ascii=False), encoding="utf-8")
        print("Metadata key audit: " + json.dumps(stats, ensure_ascii=False), flush=True)
        print("Conflict examples: " + json.dumps(differing[:8], ensure_ascii=False), flush=True)
        if report_dir is not None:
            print(f"Diagnostic output: {report_dir}", flush=True)
        if conflict_policy != "quarantine":
            raise ValueError(f"{path.name}: {len(differing)} keys have different metadata records after exact-row deduplication; inspect diagnostics, do not keep the first row arbitrarily")
        # Preserve one keyed record with only agreed attributes. A conflicting
        # donor ID becomes NA; never choose either donor or duplicate its cells.
        table = table.loc[~table[key].duplicated(keep=False)].copy()
        table["mapping_record_status"] = "resolved"
        quarantined = []
        for identifier, group in conflicts.groupby(key, dropna=False, sort=False):
            record = {column: (group[column].iloc[0] if group[column].nunique(dropna=False) == 1 else pd.NA)
                      for column in conflicts.columns}
            record["mapping_record_status"] = "conflicting_metadata"
            quarantined.append(record)
        # All source fields were read as strings. Explicit equal dtypes avoid
        # pandas' all-NA concatenation inference and its FutureWarning.
        table = table.astype("string")
        quarantine_frame = pd.DataFrame(quarantined, columns=table.columns).astype("string")
        table = pd.concat([table, quarantine_frame], ignore_index=True)
        stats["conflict_policy"] = "quarantine_shared_attributes_only"
        stats["output_rows"] = len(table)
        if report_dir is not None:
            table.loc[table.mapping_record_status.eq("conflicting_metadata")].to_csv(
                report_dir / f"{path.stem}.quarantined_records.tsv", sep="\t", index=False)
            (report_dir / f"{path.stem}.key_audit.json").write_text(
                json.dumps(stats, indent=2, ensure_ascii=False), encoding="utf-8")
    elif conflict_policy == "quarantine":
        table["mapping_record_status"] = "resolved"
    if table[key].duplicated().any():
        raise ValueError(f"{path.name}: unresolved non-unique key after processing")
    table.attrs["key_audit"] = stats
    print("Metadata key audit: " + json.dumps(stats, ensure_ascii=False), flush=True)
    return table


def sparse_row(node, row):
    if isinstance(node, h5py.Dataset) or node.attrs.get("encoding-type") not in ("csr_matrix", b"csr_matrix"):
        raise ValueError(f"{node.name}: expected CSR encoding; refusing a full matrix read")
    start, end = map(int, node["indptr"][row:row + 2])
    # Guard against anomalous pointers before touching expression storage.
    if start < 0 or end < start or end > node["data"].shape[0] or end - start > 10_000_000:
        raise ValueError(f"{node.name}: invalid/unbounded row {row}")
    indices = node["indices"][start:end]
    values = np.asarray(node["data"][start:end], dtype=np.float64)
    order = np.argsort(indices)
    return indices[order], values[order]


def run(args):
    source = Path(args.input).expanduser().resolve(strict=True)
    out = Path(args.outdir).expanduser().resolve()
    if out == source.parent or (out.exists() and any(out.iterdir())):
        raise ValueError("Choose a separate, empty output directory")
    if not 1 <= args.row_samples <= 128:
        raise ValueError("--row-samples must be between 1 and 128")
    out.mkdir(parents=True, exist_ok=True)
    diagnostics = out / "metadata_key_audit"
    diagnostics.mkdir()
    mapping = read_metadata(source.parent / "TI_atlas-sample_individual_mapping.tsv",
                            "sanger_sample_id", ["sanger_sample_id", "individual_id", "cohort", "atlas"],
                            diagnostics, conflict_policy="quarantine")
    individuals = read_metadata(source.parent / "individual_metadata.tsv", "individual_id",
                                ["individual_id", "disease_status", "ti_ses_inflamed"],
                                diagnostics, conflict_policy="quarantine")
    samples = read_metadata(source.parent / "sample_metadata.tsv", "sanger_sample_id",
                            ["sanger_sample_id", "biopsy_type"], diagnostics, conflict_policy="quarantine")
    individuals.rename(columns={"mapping_record_status": "individual_record_status"}, inplace=True)
    samples.rename(columns={"mapping_record_status": "sample_record_status"}, inplace=True)
    summary = {"script_version": SCRIPT_VERSION,
               "source": str(source), "python": sys.executable,
               "metadata_key_audits": [t.attrs["key_audit"] for t in [mapping, individuals, samples]],
               "scope": "metadata census plus sampled sparse rows; no full expression scan",
               "cells_selected": False, "expression_rows_sampled": args.row_samples}
    with h5py.File(source, "r") as h:
        obs = pd.DataFrame({name: read_column(h["obs"][name])
                            for name in ["sanger_sample_id", "label", "category"]})
        summary["n_cells"] = len(obs)
        base = obs.groupby(["sanger_sample_id", "label", "category"],
                           observed=True, dropna=False).size().rename("n_cells").reset_index()
        # Aggregate before joining: these tables are small; validate cardinality.
        coverage = base.merge(mapping, on="sanger_sample_id", how="left", validate="many_to_one")
        coverage = coverage.merge(individuals, on="individual_id", how="left", validate="many_to_one")
        coverage = coverage.merge(samples[["sanger_sample_id", "biopsy_type", "sample_record_status"]],
                                  on="sanger_sample_id", how="left", validate="many_to_one")
        if int(coverage.n_cells.sum()) != len(obs):
            raise ValueError("Metadata join changed total cell count")
        failures = {}
        for name in ["individual_id", "disease_status", "biopsy_type"]:
            bad = coverage[name].isna() | coverage[name].eq("")
            failures[name] = {"unmatched_samples": int(coverage.loc[bad, "sanger_sample_id"].nunique()),
                              "unmatched_cells": int(coverage.loc[bad, "n_cells"].sum())}
        summary["join_failures"] = failures
        summary["all_metadata_keys_matched"] = all(x["unmatched_cells"] == 0 for x in failures.values())
        summary["n_samples"] = int(coverage.sanger_sample_id.nunique())
        summary["n_individuals_matched"] = int(coverage.individual_id.nunique())
        # Quantify each conflict even if it is absent from the observed H5AD.
        per_sample = base.groupby("sanger_sample_id", observed=True, dropna=False).n_cells.sum()
        mapping_conflicts = set(mapping.loc[mapping.mapping_record_status.eq("conflicting_metadata"), "sanger_sample_id"])
        sample_conflicts = set(samples.loc[samples.sample_record_status.eq("conflicting_metadata"), "sanger_sample_id"])
        individual_conflicts = set(individuals.loc[individuals.individual_record_status.eq("conflicting_metadata"), "individual_id"])
        individual_conflict_samples = set(mapping.loc[mapping.individual_id.isin(individual_conflicts), "sanger_sample_id"])
        conflict_rows = pd.DataFrame({"sanger_sample_id": pd.Series(
            sorted(mapping_conflicts | sample_conflicts | individual_conflict_samples), dtype="string")})
        conflict_rows["mapping_table_conflict"] = conflict_rows.sanger_sample_id.isin(mapping_conflicts)
        conflict_rows["sample_metadata_conflict"] = conflict_rows.sanger_sample_id.isin(sample_conflicts)
        conflict_rows["individual_metadata_conflict"] = conflict_rows.sanger_sample_id.isin(individual_conflict_samples)
        conflict_rows = conflict_rows.merge(mapping, on="sanger_sample_id", how="left", validate="one_to_one")
        conflict_rows = conflict_rows.merge(samples[["sanger_sample_id", "biopsy_type"]],
                                             on="sanger_sample_id", how="left", validate="one_to_one")
        conflict_rows["n_cells_in_h5ad"] = conflict_rows.sanger_sample_id.map(per_sample).fillna(0).astype("int64")
        conflict_rows["present_in_h5ad"] = conflict_rows.n_cells_in_h5ad.gt(0)
        conflict_rows.to_csv(out / "conflicting_sample_impact.tsv", sep="\t", index=False)
        summary["conflicting_samples_in_h5ad"] = int(conflict_rows.present_in_h5ad.sum())
        summary["cells_with_conflicting_sample_metadata"] = int(conflict_rows.n_cells_in_h5ad.sum())
        summary["conflicting_samples_absent_from_h5ad"] = int((~conflict_rows.present_in_h5ad).sum())
        summary["conflicting_donor_assignment_policy"] = "remain missing; never use either donor ID without evidence"
        donor_known = coverage.individual_id.notna() & coverage.individual_id.ne("")
        coverage["donor_mapping_status"] = np.where(donor_known, "resolved", "unresolved")
        summary["cells_with_resolved_donor"] = int(coverage.loc[donor_known, "n_cells"].sum())
        summary["cells_with_unresolved_donor"] = int(coverage.loc[~donor_known, "n_cells"].sum())
        disease_known = coverage.disease_status.notna() & coverage.disease_status.ne("")
        site_known = coverage.biopsy_type.notna() & coverage.biopsy_type.ne("")
        complete = donor_known & disease_known & site_known
        coverage["reference_metadata_complete"] = complete
        summary["cells_with_complete_reference_metadata"] = int(coverage.loc[complete, "n_cells"].sum())
        summary["cells_with_incomplete_reference_metadata"] = int(coverage.loc[~complete, "n_cells"].sum())
        summary["metadata_conflicts_abort_audit"] = False
        summary["cells_removed_from_h5ad"] = 0
        coverage.loc[~donor_known].to_csv(out / "unresolved_donor_cell_coverage.tsv", sep="\t", index=False)
        coverage.loc[~complete].to_csv(out / "incomplete_reference_metadata_coverage.tsv", sep="\t", index=False)
        coverage.loc[complete].groupby(["label", "category"], observed=True, dropna=False).agg(
            n_cells=("n_cells", "sum"), n_samples=("sanger_sample_id", "nunique"),
            n_individuals=("individual_id", "nunique")
        ).reset_index().to_csv(out / "label_complete_metadata_support.tsv", sep="\t", index=False)
        coverage.to_csv(out / "label_category_sample_donor_coverage.tsv", sep="\t", index=False)
        support = coverage.groupby(["label", "category"], observed=True, dropna=False).agg(
            n_cells=("n_cells", "sum"), n_samples=("sanger_sample_id", "nunique"),
            n_individuals=("individual_id", "nunique"))
        support.reset_index().to_csv(out / "label_donor_support.tsv", sep="\t", index=False)
        coverage.loc[donor_known].groupby(["label", "category"], observed=True, dropna=False).agg(
            n_cells=("n_cells", "sum"), n_samples=("sanger_sample_id", "nunique"),
            n_individuals=("individual_id", "nunique")
        ).reset_index().to_csv(out / "label_resolved_donor_support.tsv", sep="\t", index=False)
        for filename, keys in [
            ("label_disease_inflammation_coverage.tsv", ["label", "category", "disease_status", "ti_ses_inflamed"]),
            ("category_disease_inflammation_coverage.tsv", ["category", "disease_status", "ti_ses_inflamed"]),
            ("cohort_atlas_coverage.tsv", ["cohort", "atlas"]),
        ]:
            coverage.groupby(keys, observed=True, dropna=False).agg(
                n_cells=("n_cells", "sum"), n_samples=("sanger_sample_id", "nunique"),
                n_individuals=("individual_id", "nunique")
            ).reset_index().to_csv(out / filename, sep="\t", index=False)
        summary["biopsy_types"] = coverage.biopsy_type.dropna().unique().tolist()
        index, _ = frame_columns(h["var"])
        genes = pd.DataFrame({"ensembl_id": read_column(h["var"][index]),
                              "gene_symbol": read_column(h["var"]["gene_symbols"])})
        duplicate = genes[genes.gene_symbol.duplicated(keep=False)].sort_values(["gene_symbol", "ensembl_id"])
        duplicate.to_csv(out / "duplicate_gene_symbols.tsv", sep="\t", index=False)
        summary["duplicate_symbol_groups"] = int(duplicate.gene_symbol.nunique())
        summary["rows_in_duplicate_symbol_groups"] = len(duplicate)
        markers = ["DEFA5", "DEFA6", "ITLN2", "PLA2G2A", "LYZ", "MUC2", "CLCA1", "CHGA", "CHGB", "TPH1",
                   "BEST4", "OTOP2", "LGR5", "OLFM4", "MKI67", "POU2F3", "DCLK1", "SOX10", "GFAP", "PLP1", "S100B"]
        present = set(genes.gene_symbol)
        pd.DataFrame({"gene_symbol": markers, "present_in_gene_universe": [m in present for m in markers],
                      "scope": "gene annotation only; not expression evidence"}).to_csv(out / "identity_marker_gene_presence.tsv", sep="\t", index=False)
        positions = np.sort(np.random.default_rng(20261006).choice(len(obs), min(args.row_samples, len(obs)), replace=False))
        summaries = []
        for row in positions:
            ci, cv = sparse_row(h["layers/counts"], int(row))
            xi, xv = sparse_row(h["X"], int(row))
            li, lv = sparse_row(h["layers/log1p_cp10k"], int(row))
            finite = cv[np.isfinite(cv)]
            aligned_x = np.array_equal(ci, xi)
            aligned_log = np.array_equal(ci, li)
            total = float(cv.sum())
            expected = np.log1p(cv / total * 10000) if total > 0 and np.isfinite(total) else None
            summaries.append({"cell_row": int(row), "label": str(obs.iloc[row]["label"]),
                "sample": str(obs.iloc[row]["sanger_sample_id"]), "stored_counts_entries": len(cv),
                "counts_nonfinite": int((~np.isfinite(cv)).sum()), "counts_negative": int((cv < 0).sum()),
                "counts_noninteger": int((np.abs(finite - np.rint(finite)) > 1e-5).sum()),
                "counts_row_sum": total, "x_indices_equal_counts": aligned_x,
                "x_matches_counts": bool(aligned_x and np.allclose(cv, xv, rtol=1e-6, atol=1e-6)),
                "log_indices_equal_counts": aligned_log,
                "log1p_cp10k_matches_formula": bool(aligned_log and expected is not None and np.allclose(lv, expected, rtol=1e-5, atol=1e-5))})
        pd.DataFrame(summaries).to_csv(out / "sampled_expression_rows.tsv", sep="\t", index=False)
        summary["sampled_counts_nonnegative_finite_integer"] = all(
            r["counts_nonfinite"] == r["counts_negative"] == r["counts_noninteger"] == 0 for r in summaries)
        summary["sampled_X_equals_counts"] = all(r["x_matches_counts"] for r in summaries)
        summary["sampled_log1p_cp10k_matches_formula"] = all(r["log1p_cp10k_matches_formula"] for r in summaries)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"Output: {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=SCRIPT_VERSION)
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--row-samples", type=int, default=24)
    run(parser.parse_args())
