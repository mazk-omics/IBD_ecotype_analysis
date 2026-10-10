"""Small fixture check: annotation decoding and zero expression-array reads."""
import argparse
import importlib.util
import json
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("inspector", Path(__file__).with_name("inspect_ibdverse_h5ad.py"))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

with tempfile.TemporaryDirectory(prefix="ibdverse_inspector_") as folder:
    base = Path(folder)
    source = base / "Full_cohort.h5ad"
    with h5py.File(source, "w") as h:
        obs = h.create_group("obs")
        obs.attrs["_index"] = "_index"
        obs.attrs["column-order"] = ["cell_type", "donor_id", "condition", "optional", "qc_metric"]
        obs.create_dataset("_index", data=np.array([f"c{i}" for i in range(5)], dtype=h5py.string_dtype()))
        cat = obs.create_group("cell_type")
        cat.attrs["encoding-type"] = "categorical"
        cat.create_dataset("categories", data=np.array(["Goblet", "Enterocyte", "Unused"], dtype=h5py.string_dtype()))
        cat.create_dataset("codes", data=[0, 0, 1, 1, -1])
        obs.create_dataset("donor_id", data=np.array(["d1", "d1", "d2", "d2", "d3"], dtype=h5py.string_dtype()))
        obs.create_dataset("condition", data=np.array(["Healthy", "Healthy", "CD", "CD", "CD"], dtype=h5py.string_dtype()))
        nullable = obs.create_group("optional")
        nullable.attrs["encoding-type"] = "nullable-integer"
        nullable.create_dataset("values", data=[1, 2, 3, 4, 5])
        nullable.create_dataset("mask", data=[False, True, False, False, False])
        obs.create_dataset("qc_metric", data=[1.1, 2.2, 3.3, 4.4, 5.5])
        var = h.create_group("var")
        var.attrs["_index"] = "_index"
        var.attrs["column-order"] = ["gene_symbol"]
        var.create_dataset("_index", data=np.array(["ENSG1", "ENSG2", "ENSG3"], dtype=h5py.string_dtype()))
        var.create_dataset("gene_symbol", data=np.array(["MUC2", "CHGA", "MUC2"], dtype=h5py.string_dtype()))
        h.create_dataset("X", data=np.zeros((5, 3), dtype=np.float32))
        csr = h.create_group("layers/counts")
        csr.attrs["encoding-type"] = "csr_matrix"
        csr.attrs["shape"] = [5, 3]
        csr.create_dataset("data", data=[2., 3.])
        csr.create_dataset("indices", data=[0, 1])
        csr.create_dataset("indptr", data=[0, 1, 1, 2, 2, 2])
    original = h5py.Dataset.__getitem__

    def reject_expression_reads(self, selection):
        if self.name == "/X" or self.name.startswith(("/layers/", "/raw/X")):
            raise AssertionError(f"Expression data read: {self.name}")
        return original(self, selection)

    h5py.Dataset.__getitem__ = reject_expression_reads
    try:
        module.audit(argparse.Namespace(input=str(source), outdir=str(base / "audit"),
                                        group_cols=["cell_type", "donor_id", "condition"]))
    finally:
        h5py.Dataset.__getitem__ = original
    summary = json.loads((base / "audit/h5ad_structure.json").read_text(encoding="utf-8"))
    assert summary["n_cells"] == 5 and summary["n_genes"] == 3
    assert summary["errors"] == [] and not summary["expression_values_read"]
    coverage = pd.read_csv(base / "audit/selected_metadata_coverage.tsv", sep="\t")
    assert coverage.n_cells.sum() == 5 and len(coverage) == 3
    assert set(coverage.cell_type) == {"Goblet", "Enterocyte", "<MISSING>"}
    fields = pd.read_csv(base / "audit/obs_columns.tsv", sep="\t").set_index("column")
    assert fields.loc["cell_type", "missing"] == 1
    assert fields.loc["optional", "missing"] == 1
    genes = pd.read_csv(base / "audit/gene_columns.tsv", sep="\t").set_index("column")
    assert genes.loc["gene_symbol", "n_unique_excluding_missing"] == 2
    print("PASS: categorical/missing metadata, duplicate gene symbols, coverage, zero expression reads")
