#!/usr/bin/env python3
"""Read-only H5AD structure/annotation audit; never read expression arrays.

Project environment: ibd_refaudit (environments/conda/ibd_refaudit.yml).
Requires h5py, numpy, pandas. Does not require Scanpy or AnnData.
Usage: python inspect_ibdverse_h5ad.py --input SOURCE.h5ad --outdir NEW_DIR
Optional: --group-cols cell_type donor_id disease_status
Schema: https://anndata.readthedocs.io/en/stable/fileformat-prose.html
"""
import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROLE_PATTERN = re.compile(
    r"cell.?type|annotat|label|cluster|lineage|compartment|donor|patient|"
    r"individual|subject|sample|biopsy|disease|diagnos|condition|status|"
    r"phenotype|inflamm|ses.?cd|region|site|tissue|location|cohort|qc|pass", re.I
)


def plain(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        return [plain(x) for x in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_dataset(node, limit=None):
    selection = () if limit is None else slice(0, limit)
    if h5py.check_string_dtype(node.dtype) is not None:
        return node.asstr()[selection]
    return node[selection]


def read_column(node, limit=None):
    if isinstance(node, h5py.Dataset):
        return pd.Series(read_dataset(node, limit))
    encoding = plain(node.attrs.get("encoding-type", ""))
    if encoding == "categorical" or {"codes", "categories"} <= set(node.keys()):
        codes = np.asarray(read_dataset(node["codes"], limit), dtype=np.int64)
        categories = read_dataset(node["categories"])
        return pd.Series(pd.Categorical.from_codes(codes, categories))
    if {"values", "mask"} <= set(node.keys()):
        values = pd.Series(read_dataset(node["values"], limit))
        mask = np.asarray(read_dataset(node["mask"], limit), dtype=bool)
        return values.mask(mask)
    raise ValueError(f"Unsupported annotation encoding: {node.name} ({encoding})")


def frame_columns(group):
    index = plain(group.attrs.get("_index", "_index"))
    columns = plain(group.attrs.get("column-order", []))
    if isinstance(columns, str):
        columns = [columns]
    return index, columns or [k for k in group.keys() if k != index]


def describe(node):
    result = {"path": node.name, "encoding": plain(node.attrs.get("encoding-type", ""))}
    if isinstance(node, h5py.Dataset):
        result.update(shape=list(node.shape), dtype=str(node.dtype), compression=node.compression)
    else:
        result["children"] = list(node.keys())
        if "shape" in node.attrs:
            result["shape"] = plain(node.attrs["shape"])
        result["datasets"] = {
            k: {"shape": list(v.shape), "dtype": str(v.dtype)}
            for k, v in node.items() if isinstance(v, h5py.Dataset)
        }
    return result


def audit(args):
    source = Path(args.input).expanduser().resolve(strict=True)
    out = Path(args.outdir).expanduser().resolve()
    if out == source.parent or out == source:
        raise ValueError("Use a separate audit output directory, not the source directory.")
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory is not empty; choose a new directory.")
    out.mkdir(parents=True, exist_ok=True)
    summary = {
        "source": str(source), "source_bytes": source.stat().st_size,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": {"python_executable": sys.executable, "python_version": sys.version,
                    "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
                    "h5py": h5py.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        "expression_values_read": False,
        "expression_scale": "NOT VERIFIED: layer names/dtypes do not establish raw counts or log scale",
        "errors": [],
    }
    obs_rows, value_rows, gene_rows = [], [], []
    with h5py.File(source, "r") as handle:
        summary["root_keys"] = list(handle.keys())
        obs_index, obs_cols = frame_columns(handle["obs"])
        var_index, var_cols = frame_columns(handle["var"])
        summary["n_cells"] = int(handle["obs"][obs_index].shape[0])
        summary["n_genes"] = int(handle["var"][var_index].shape[0])
        summary["obs_columns"] = obs_cols
        summary["var_columns"] = var_cols
        structures = []
        for name in ("X", "raw/X"):
            if name in handle:
                structures.append(describe(handle[name]))
        if "layers" in handle:
            structures.extend(describe(v) for v in handle["layers"].values())
        summary["expression_structures"] = structures
        summary["uns_keys"] = list(handle["uns"].keys()) if "uns" in handle else []
        groups = {}
        for column in obs_cols:
            node = handle["obs"][column]
            # Read categorical annotations and relevant identity/sample metadata.
            # Continuous QC measurements get a preview, not full-column statistics.
            full = (isinstance(node, h5py.Group)
                    or bool(ROLE_PATTERN.search(column)) or column in args.group_cols)
            try:
                series = read_column(node, None if full else 256)
                counts = series.value_counts(dropna=False)
                counts = counts[counts > 0]
                obs_rows.append({
                    "column": column, "encoding": plain(node.attrs.get("encoding-type", "")),
                    "scope": "all_cells" if full else "first_256_cells",
                    "rows_read": len(series), "n_unique_excluding_missing": int(series.nunique()),
                    "missing": int(series.isna().sum()),
                    "example_values": json.dumps([plain(v) for v in series.dropna().unique()[:12]], ensure_ascii=False),
                })
                if full:
                    for value, count in counts.iloc[:200].items():
                        value_rows.append({"column": column, "value": "<MISSING>" if pd.isna(value) else str(value),
                                           "n_cells": int(count), "top_200_only": len(counts) > 200})
                if column in args.group_cols:
                    groups[column] = series.astype("string").fillna("<MISSING>")
            except (ValueError, TypeError, KeyError) as exc:
                summary["errors"].append(f"obs/{column}: {exc}")
        for frame_name in ("var", "raw/var"):
            if frame_name not in handle:
                continue
            index, columns = frame_columns(handle[frame_name])
            for column in dict.fromkeys([index] + columns):
                try:
                    series = read_column(handle[frame_name][column])
                    gene_rows.append({"frame": frame_name, "column": column, "n_rows": len(series),
                                      "n_unique_excluding_missing": int(series.nunique()),
                                      "missing": int(series.isna().sum()),
                                      "example_values": json.dumps([plain(v) for v in series.dropna().unique()[:12]], ensure_ascii=False)})
                except (ValueError, TypeError, KeyError) as exc:
                    summary["errors"].append(f"{frame_name}/{column}: {exc}")
        if args.group_cols:
            missing = sorted(set(args.group_cols) - set(groups))
            if missing:
                raise ValueError(f"Requested grouping columns absent/unreadable: {missing}")
            table = pd.DataFrame(groups).groupby(args.group_cols, dropna=False, observed=True).size()
            table.rename("n_cells").reset_index().to_csv(out / "selected_metadata_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(obs_rows).to_csv(out / "obs_columns.tsv", sep="\t", index=False)
    pd.DataFrame(value_rows, columns=["column", "value", "n_cells", "top_200_only"]).to_csv(out / "obs_value_counts.tsv", sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(out / "gene_columns.tsv", sep="\t", index=False)
    sidecars = {}
    for name in ("sample_metadata.tsv", "individual_metadata.tsv", "TI_atlas-sample_individual_mapping.tsv"):
        path = source.parent / name
        if path.exists():
            table = pd.read_csv(path, sep="\t", nrows=0)
            sidecars[name] = {"bytes": path.stat().st_size, "columns": table.columns.tolist()}
    summary["sidecar_metadata"] = sidecars
    (out / "h5ad_structure.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    report = [f"Source: {source}", f"Cells: {summary['n_cells']:,}; genes: {summary['n_genes']:,}",
              f"Python: {sys.executable}; Conda: {summary['runtime']['conda_env']}",
              f"File bytes: {summary['source_bytes']:,}", "Expression arrays were NOT read; source opened read-only.",
              "Expression scale remains unverified; no signature or cell sampling performed.",
              "Cell/donor/disease labels must be reviewed before choosing grouping columns.",
              f"Output: {out}", f"Annotation decoding errors: {len(summary['errors'])}",
              "OBS columns: " + ", ".join(summary["obs_columns"]),
              "Expression structures: " + ", ".join(x["path"] for x in structures)]
    (out / "REPORT.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--group-cols", nargs="*", default=[])
    audit(parser.parse_args())
