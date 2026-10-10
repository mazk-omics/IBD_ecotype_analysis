"""Small synthetic CSR integration checks; no access to the real cohort needed."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from export_ibdverse_cibersortx_pilot import run, WHOLE_TYPE_ORDER


LABELS = ["Enterocytes top villus", "Stem cells MKI67+", "Goblet cells bottom villus",
          "Enteroendocrine cells", "Tuft cells", "Enterocytes BEST4+", "B cells naive",
          "Plasma B cell XBP1/CD38+", "T cells CD4+ naive", "T cells CD8+ tissue-resident",
          "T cells gd", "ILCs", "Monocytes S100A8/9+", "cDC1", "Mast cells",
          "Fibroblasts", "Pericyte cells", "Endothelial cells", "T cells CD4- CD8-", "T cells CD4- CD8-lo"]


def string_dataset(group, name, values):
    group.create_dataset(name, data=np.asarray(values, dtype=object), dtype=h5py.string_dtype("utf-8"))


def make_fixture(root):
    records, dense, coverage = [], [], []
    for disease, inflammation in [("Healthy", "uninflamed"), ("CD", "inflamed"), ("CD", "uninflamed")]:
        for donor_number in range(20):
            sample = f"s_{disease}_{inflammation}_{donor_number}"
            donor = "individual_" + sample
            for label in LABELS:
                vals = [100 + donor_number, int(label.startswith("T cells")) * 10,
                        int(label == "Enterocytes BEST4+") * 40, 2, 3, 0,
                        int(label.startswith("B cells")) * 12]
                row = {"sanger_sample_id": sample, "label": label, "category": "test_category",
                       "conf_score": 0.95, "n_genes_by_counts": int(np.count_nonzero(vals)),
                       "total_counts": sum(vals), "pct_counts_gene_group__mito_transcript": 3.0}
                records.append(row)
                dense.append(vals)
                coverage.append({"sanger_sample_id": sample, "label": label, "category": "test_category",
                    "n_cells": 1, "individual_id": donor, "disease_status": disease,
                    "ti_ses_inflamed": inflammation, "biopsy_type": "terminal ileum", "cohort": "Discovery",
                    "atlas": "not_atlas70", "reference_metadata_complete": True})
    source = root / "fixture.h5ad"
    dense = np.asarray(dense, dtype=np.int32)
    obs = pd.DataFrame(records)
    with h5py.File(source, "w") as h:
        o = h.create_group("obs")
        o.attrs["_index"] = "_index"
        o.attrs["column-order"] = np.asarray(obs.columns, dtype=h5py.string_dtype())
        string_dataset(o, "_index", [f"cell_{i}" for i in range(len(obs))])
        for name in obs:
            if obs[name].dtype == object or name == "conf_score":
                # Exercise categorical decoding as in the real source.
                cat = pd.Categorical(obs[name].astype(str) if name == "conf_score" else obs[name])
                g = o.create_group(name)
                g.attrs["encoding-type"] = "categorical"
                g.create_dataset("codes", data=cat.codes)
                string_dataset(g, "categories", cat.categories)
            else:
                o.create_dataset(name, data=obs[name].to_numpy())
        v = h.create_group("var")
        v.attrs["_index"] = "_index"
        v.attrs["column-order"] = np.asarray(["gene_symbols"], dtype=h5py.string_dtype())
        string_dataset(v, "_index", [f"ENSG_test_{i}" for i in range(7)])
        string_dataset(v, "gene_symbols", ["EPCAM", "CD3D", "BEST4", "RGS5", "RGS5", "ZERO", "CD79A"])
        indices, values, pointers = [], [], [0]
        for row in dense:
            nonzero = np.flatnonzero(row)
            indices.extend(nonzero)
            values.extend(row[nonzero])
            pointers.append(len(values))
        csr = h.create_group("layers/counts")
        csr.attrs["encoding-type"] = "csr_matrix"
        csr.attrs["shape"] = dense.shape
        csr.create_dataset("indices", data=np.asarray(indices, dtype=np.int32))
        csr.create_dataset("data", data=np.asarray(values, dtype=np.float32))
        csr.create_dataset("indptr", data=np.asarray(pointers, dtype=np.int64))
    report = root / "coverage.tsv"
    pd.DataFrame(coverage).to_csv(report, sep="\t", index=False)
    return source, report, dense


def options(source, report, out, cells=8):
    return argparse.Namespace(input=str(source), coverage=str(report), outdir=str(out),
        cells_per_type=cells, holdout_fraction=0.2, max_per_donor=5, seed=20261006,
        min_conf_score=None, review_cells_per_pending_label=4)


def main():
    with tempfile.TemporaryDirectory(prefix="ibdverse_pilot_test_") as directory:
        root = Path(directory)
        source, report, original = make_fixture(root)
        source_digest = hashlib.sha256(source.read_bytes()).hexdigest()
        run(options(source, report, root / "first"))
        first = root / "first"
        summary = json.loads((first / "summary.json").read_text())
        assert summary["web_reference_cells"] == 144 and summary["identity_review_cells"] == 8
        assert summary["holdout_donors"] == 12 and summary["training_donors"] == 48
        assert summary["holdout_donors_used_for_extraction"] == 0
        manifest = pd.read_csv(first / "selected_cell_manifest.tsv", sep="\t")
        web = manifest.loc[manifest.usage.eq("web_reference")]
        assert set(web.pilot_type) == set(WHOLE_TYPE_ORDER)
        assert web.groupby("pilot_type").size().eq(8).all()
        assert set(web.pilot_type).isdisjoint({"T_other"})
        assert manifest.total_counts_matches_obs.all() and manifest.n_genes_matches_obs.all()
        with (first / "IBDverse_TI_18types_scRNA_reference.txt").open() as stream:
            rows = list(csv.reader(stream, delimiter="\t"))
        assert rows[0][1:] == web.web_column.tolist() and len(set(rows[0][1:])) == 144
        gene_names = [r[0] for r in rows[1:]]
        assert gene_names == ["EPCAM", "CD3D", "BEST4", "CD79A"]
        matrix = np.asarray([[int(v) for v in r[1:]] for r in rows[1:]])
        expected = original[web.source_cell_row.to_numpy()][:, [0, 1, 2, 6]].T
        np.testing.assert_array_equal(matrix, expected)
        marker = pd.read_csv(first / "marker_expression_summary.tsv", sep="\t")
        rg = marker.loc[marker.gene_symbol.eq("RGS5")]
        assert rg.source_gene_id.nunique() == 2 and (~rg.gene_exported_to_web).all()
        with h5py.File(first / "selected_counts.h5", "r") as h:
            assert tuple(h["counts"].attrs["shape"]) == (152, 7)
            reconstructed = np.zeros((152, 7), dtype=int)
            ip, ix, data = h["counts/indptr"][:], h["counts/indices"][:], h["counts/data"][:]
            for i in range(152):
                reconstructed[i, ix[ip[i]:ip[i + 1]]] = data[ip[i]:ip[i + 1]]
            np.testing.assert_array_equal(reconstructed, original[manifest.source_cell_row.to_numpy()])
        run(options(source, report, root / "repeat"))
        assert (first / "IBDverse_TI_18types_scRNA_reference.txt").read_bytes() == (
            root / "repeat/IBDverse_TI_18types_scRNA_reference.txt").read_bytes()
        assert hashlib.sha256(source.read_bytes()).hexdigest() == source_digest
        try:
            run(options(source, report, root / "shortage", cells=200))
        except ValueError as exc:
            assert "cannot meet the quota" in str(exc)
        else:
            raise AssertionError("Unattainable quotas must fail with diagnostics")
        assert not (root / "shortage/IBDverse_TI_18types_scRNA_reference.txt").exists()
        try:
            run(options(source, report, first))
        except ValueError as exc:
            assert "empty output" in str(exc)
        else:
            raise AssertionError("Existing outputs must never be overwritten")
        with h5py.File(source, "r+") as h:
            ix = h["layers/counts/indices"][:]
            data = h["layers/counts/data"][:]
            data[ix == 0] = -1
            h["layers/counts/data"][:] = data
        try:
            run(options(source, report, root / "negative"))
        except ValueError as exc:
            assert "nonnegative integers" in str(exc)
        else:
            raise AssertionError("Negative counts must be rejected")
    print("PASS: exact matrix/manifest alignment, donor holdout, deterministic sampling, duplicate loci preservation, shortages, read-only source and negative-count rejection")


if __name__ == "__main__":
    main()
