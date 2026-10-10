"""Verify donor join cardinality, sparse scale checks and duplicate symbols."""
import argparse
import tempfile
import warnings
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from audit_ibdverse_reference_inputs import run, read_metadata

warnings.simplefilter("error", FutureWarning)


with tempfile.TemporaryDirectory(prefix="ibdverse_reference_inputs_") as folder:
    root = Path(folder)
    counts = np.array([[2, 0, 1, 0], [0, 3, 0, 2], [1, 1, 0, 0],
                       [0, 4, 0, 1], [2, 2, 2, 0], [1, 0, 0, 4]], dtype=float)
    with h5py.File(root / "Full_cohort.h5ad", "w") as h:
        obs = h.create_group("obs")
        for name, values in {"sanger_sample_id": ["s1", "s1", "s2", "s2", "s3", "s3"],
                             "label": ["Goblet", "Goblet", "Enterocyte", "Enterocyte", "Goblet", "Enterocyte"],
                             "category": ["Secretory", "Secretory", "Enterocyte", "Enterocyte", "Secretory", "Enterocyte"]}.items():
            obs.create_dataset(name, data=np.array(values, dtype=h5py.string_dtype()))
        var = h.create_group("var")
        var.attrs["_index"] = "_index"
        var.create_dataset("_index", data=np.array(["g1", "g2", "g3", "g4"], dtype=h5py.string_dtype()))
        var.create_dataset("gene_symbols", data=np.array(["MUC2", "CHGA", "MUC2", "DEFA5"], dtype=h5py.string_dtype()))
        for name, matrix in {"X": counts, "layers/counts": counts,
                             "layers/log1p_cp10k": np.log1p(counts / counts.sum(axis=1)[:, None] * 10000)}.items():
            group = h.create_group(name)
            group.attrs["encoding-type"] = "csr_matrix"
            group.attrs["shape"] = matrix.shape
            rows, cols = np.nonzero(matrix)
            group.create_dataset("data", data=matrix[rows, cols].astype("float32"))
            group.create_dataset("indices", data=cols.astype("int64"))
            group.create_dataset("indptr", data=np.r_[0, np.cumsum((matrix != 0).sum(axis=1))].astype("int64"))
    pd.DataFrame({"sanger_sample_id": ["s1", "s2", "s3"], "individual_id": ["d1", "d1", "d2"],
                  "cohort": ["Discovery", "Discovery", "Replication"], "atlas": ["yes", "yes", "no"]}).to_csv(
        root / "TI_atlas-sample_individual_mapping.tsv", sep="\t", index=False)
    pd.DataFrame({"individual_id": ["d1", "d2"], "disease_status": ["CD", "Healthy"],
                  "ti_ses_inflamed": ["inflamed", "NA"]}).to_csv(root / "individual_metadata.tsv", sep="\t", index=False)
    pd.DataFrame({"sanger_sample_id": ["s1", "s2", "s3"], "biopsy_type": ["terminal ileum"] * 3}).to_csv(
        root / "sample_metadata.tsv", sep="\t", index=False)
    result = run(argparse.Namespace(input=str(root / "Full_cohort.h5ad"), outdir=str(root / "audit"), row_samples=6))
    import json
    summary = json.loads((root / "audit/summary.json").read_text(encoding="utf-8"))
    assert summary["all_metadata_keys_matched"] and summary["n_individuals_matched"] == 2
    assert summary["n_samples"] == 3 and summary["n_cells"] == 6
    assert summary["sampled_counts_nonnegative_finite_integer"]
    assert summary["sampled_X_equals_counts"] and summary["sampled_log1p_cp10k_matches_formula"]
    assert summary["duplicate_symbol_groups"] == 1 and summary["rows_in_duplicate_symbol_groups"] == 2
    table = pd.read_csv(root / "TI_atlas-sample_individual_mapping.tsv", sep="\t")
    pd.concat([table, table.iloc[[0]]]).to_csv(root / "bad_mapping.tsv", sep="\t", index=False)
    deduped = read_metadata(root / "bad_mapping.tsv", "sanger_sample_id", ["sanger_sample_id"])
    assert len(deduped) == 3 and deduped.attrs["key_audit"]["exact_duplicate_rows_removed"] == 1
    blank = table.iloc[[0]].copy()
    blank["sanger_sample_id"] = "NA"
    pd.concat([table, blank]).to_csv(root / "blank_mapping.tsv", sep="\t", index=False)
    normalized = read_metadata(root / "blank_mapping.tsv", "sanger_sample_id", ["sanger_sample_id"])
    assert len(normalized) == 3 and normalized.attrs["key_audit"]["missing_key_rows"] == 1
    conflict = table.iloc[[0]].copy()
    conflict["individual_id"] = "wrong_donor"
    pd.concat([table, conflict]).to_csv(root / "conflicting_mapping.tsv", sep="\t", index=False)
    diagnostics = root / "conflict_diagnostics"
    diagnostics.mkdir()
    try:
        read_metadata(root / "conflicting_mapping.tsv", "sanger_sample_id", ["sanger_sample_id"], diagnostics)
    except ValueError:
        pass
    else:
        raise AssertionError("Conflicting donor keys must be rejected")
    details = json.loads((diagnostics / "conflicting_mapping.conflict_details.json").read_text(encoding="utf-8"))
    assert details[0]["differing_columns"] == ["individual_id"]
    assert len(pd.read_csv(diagnostics / "conflicting_mapping.conflicting_rows.tsv", sep="\t")) == 2
    quarantined = read_metadata(root / "conflicting_mapping.tsv", "sanger_sample_id", ["sanger_sample_id"],
                                conflict_policy="quarantine")
    ambiguous = quarantined.loc[quarantined.sanger_sample_id.eq("s1")].iloc[0]
    assert len(quarantined) == 3 and pd.isna(ambiguous.individual_id)
    assert ambiguous.cohort == "Discovery" and ambiguous.atlas == "yes"
    # An observed ambiguous sample must remain in cell totals with no donor
    # assigned. Other samples from its candidate donor remain available.
    pd.concat([table, conflict]).to_csv(root / "TI_atlas-sample_individual_mapping.tsv", sep="\t", index=False)
    run(argparse.Namespace(input=str(root / "Full_cohort.h5ad"), outdir=str(root / "audit_ambiguous"), row_samples=2))
    unresolved = json.loads((root / "audit_ambiguous/summary.json").read_text(encoding="utf-8"))
    assert unresolved["n_cells"] == 6 and unresolved["cells_removed_from_h5ad"] == 0
    assert unresolved["cells_with_unresolved_donor"] == 2
    assert unresolved["cells_with_resolved_donor"] == 4
    assert unresolved["conflicting_samples_in_h5ad"] == 1
    assert unresolved["n_individuals_matched"] == 2
    # A conflict outside the H5AD must not make its observed mappings fail.
    external = pd.DataFrame({"sanger_sample_id": ["outside", "outside"], "individual_id": ["d1", "d2"],
                             "cohort": ["Replication", "Replication"], "atlas": ["no", "no"]})
    pd.concat([table, external]).to_csv(root / "TI_atlas-sample_individual_mapping.tsv", sep="\t", index=False)
    run(argparse.Namespace(input=str(root / "Full_cohort.h5ad"), outdir=str(root / "audit_external"), row_samples=2))
    external_summary = json.loads((root / "audit_external/summary.json").read_text(encoding="utf-8"))
    assert external_summary["all_metadata_keys_matched"]
    assert external_summary["cells_with_unresolved_donor"] == 0
    assert external_summary["conflicting_samples_absent_from_h5ad"] == 1
    # The real failure: three sample IDs have conflicting tissue attributes.
    table.to_csv(root / "TI_atlas-sample_individual_mapping.tsv", sep="\t", index=False)
    site_base = pd.DataFrame({"sanger_sample_id": ["s1", "s2", "s3"],
                              "biopsy_type": ["terminal ileum"] * 3, "chip_lot": ["lot1"] * 3})
    site_other = site_base.copy()
    site_other["biopsy_type"] = "rectum"
    site_other.loc[0, "chip_lot"] = "lot2"
    pd.concat([site_base, site_other]).to_csv(root / "sample_metadata.tsv", sep="\t", index=False)
    run(argparse.Namespace(input=str(root / "Full_cohort.h5ad"), outdir=str(root / "audit_sites"), row_samples=2))
    site_summary = json.loads((root / "audit_sites/summary.json").read_text(encoding="utf-8"))
    assert site_summary["n_cells"] == 6 and site_summary["cells_removed_from_h5ad"] == 0
    assert site_summary["cells_with_unresolved_donor"] == 0
    assert site_summary["cells_with_complete_reference_metadata"] == 0
    assert site_summary["conflicting_samples_in_h5ad"] == 3
    assert site_summary["join_failures"]["biopsy_type"]["unmatched_cells"] == 6
    # Conflicting individual phenotypes also yield a complete audit, not a crash.
    site_base.to_csv(root / "sample_metadata.tsv", sep="\t", index=False)
    individual_base = pd.read_csv(root / "individual_metadata.tsv", sep="\t", dtype="string", keep_default_na=False)
    individual_other = individual_base.iloc[[0]].copy()
    individual_other["disease_status"] = "Healthy"
    pd.concat([individual_base, individual_other]).to_csv(root / "individual_metadata.tsv", sep="\t", index=False)
    run(argparse.Namespace(input=str(root / "Full_cohort.h5ad"), outdir=str(root / "audit_individuals"), row_samples=2))
    individual_summary = json.loads((root / "audit_individuals/summary.json").read_text(encoding="utf-8"))
    assert individual_summary["n_cells"] == 6 and individual_summary["n_individuals_matched"] == 2
    assert individual_summary["join_failures"]["disease_status"]["unmatched_cells"] == 4
    assert individual_summary["cells_with_complete_reference_metadata"] == 2
    print("PASS: all-table conflicts, dtype warnings as errors, preserved cell totals, sparse normalization")
