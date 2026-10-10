#!/usr/bin/env python3
"""End-to-end fixture for keyed row alignment, donor holdout and counts safety."""
import argparse
import hashlib
import json
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from export_gse282122_colon_web_pilot import TYPES, run


def main():
    with tempfile.TemporaryDirectory(prefix="gse_colon_pilot_test_") as tmp:
        root = Path(tmp)
        source, candidate = root / "source.h5ad", root / "mapped.parquet"
        records = []
        for disease, n in [("CD", 12), ("UC", 12), ("Healthy", 3)]:
            for donor in range(n):
                for identity in TYPES:
                    for rep in range(6):
                        name = f"cell{len(records)}"
                        records.append({"obs_name": name, "cell_id": name, "Patient": f"{disease}{donor}",
                            "Disease": disease, "sample_id": f"{disease}{donor}sample",
                            "Site": "Ascending_Colon", "Treatment": "ND" if disease == "Healthy" else "Pre",
                            "Inflammation": "Healthy" if disease == "Healthy" else "Inflamed",
                            "reference_eligible": True, "coarse_identity": identity,
                            "identity_status": "CONFIDENT", "primary_reference_action": "KEEP",
                            "predicted_doublets": rep == 5, "final_analysis": identity + " subtype",
                            "n_genes_by_counts": 4, "total_counts": 112 + len(records) % 7, "pct_counts_mt": 1.})
        df = pd.DataFrame(records)
        n = len(df)
        # Shuffle the metadata relative to physical source rows to exercise the join.
        df.sample(frac=1, random_state=7).to_parquet(candidate, index=False)
        vals = np.column_stack([100 + np.arange(n) % 7, np.ones(n) * 10, np.ones(n), np.ones(n)]).astype(np.int32)
        with h5py.File(source, "w") as h:
            obs = h.create_group("obs"); obs.attrs["_index"] = "_index"
            strings = h5py.string_dtype("utf-8")
            obs.create_dataset("_index", data=df.obs_name.to_numpy(dtype=object), dtype=strings)
            for field in ["sample_id", "Patient", "Disease", "Site", "final_analysis"]:
                obs.create_dataset(field, data=df[field].to_numpy(dtype=object), dtype=strings)
            var = h.create_group("var")
            var.create_dataset("gene_id", data=np.asarray(["g0", "g1", "g2", "g3", "g4"], dtype=object), dtype=strings)
            var.create_dataset("gene_symbol", data=np.asarray(["EPCAM", "COL1A1", "RGS5", "RGS5", "ZERO"], dtype=object), dtype=strings)
            x = h.create_group("X"); x.attrs["encoding-type"] = "csr_matrix"; x.attrs["shape"] = [n, 5]
            x.create_dataset("data", data=vals.ravel())
            x.create_dataset("indices", data=np.tile(np.arange(4), n))
            x.create_dataset("indptr", data=np.arange(n + 1) * 4)
        before = hashlib.sha256(source.read_bytes()).hexdigest()
        def options(name):
            return argparse.Namespace(input=str(source), candidate=str(candidate), outdir=str(root/name), cells_per_type=20, seed=20261006)
        run(options("first"))
        first = root / "first"
        result = json.loads((first / "summary.json").read_text())
        assert result["web_reference_cells"] == 320 and result["web_reference_genes"] == 2
        assert result["holdout_donors_used"] == 0 and result["web_format_valid"]
        selected = pd.read_csv(first / "selected_cell_manifest.tsv", sep="\t")
        assert selected.groupby("coarse_identity").size().eq(20).all()
        assert not selected.predicted_doublets.any()
        split = pd.read_csv(first / "donor_split.tsv", sep="\t")
        assert not selected.Patient.isin(split.loc[split.split.eq("holdout"), "Patient"]).any()
        web = pd.read_csv(first / "GSE282122_Colon_16types_scRNA_reference.txt", sep="\t", index_col=0)
        assert web.columns.tolist() == selected.web_column.tolist()
        assert np.array_equal(web.loc["EPCAM"].to_numpy(), vals[selected.source_cell_row, 0])
        assert np.array_equal(web.loc["COL1A1"].to_numpy(), vals[selected.source_cell_row, 1])
        with h5py.File(first / "selected_counts.h5") as h:
            assert h["counts"].attrs["shape"].tolist() == [320, 5]
        assert hashlib.sha256(source.read_bytes()).hexdigest() == before
        run(options("repeat"))
        assert (first / "GSE282122_Colon_16types_scRNA_reference.txt").read_bytes() == (
            root / "repeat/GSE282122_Colon_16types_scRNA_reference.txt").read_bytes()
        try:
            run(options("first"))
        except ValueError as exc:
            assert "overwrite" in str(exc)
        else:
            raise AssertionError("Existing output overwritten")
        with h5py.File(source, "r+") as h:
            h["X/data"][:] = -1
        try:
            run(options("negative"))
        except ValueError as exc:
            assert "nonnegative integers" in str(exc)
        else:
            raise AssertionError("Negative counts accepted")
    print("PASS: keyed source alignment, donor holdout, deterministic sampling, duplicate-gene exclusion, doublet exclusion, no overwrite and negative-count rejection")


if __name__ == "__main__":
    main()
